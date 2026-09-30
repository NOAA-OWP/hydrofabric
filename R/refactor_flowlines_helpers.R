# fast_refactor_helpers.R
# Reconcile + split helpers for the refactor engine (refactor_flowlines() lives in
# refactor_flowlines.R). data.table integer-index graphs, no repeated joins in hot
# loops, furrr-parallel split.
#
# Exposes: reconcile_collapsed_flowlines() (reconcile step), .fast_split()
# /.snap_events_to_reconciled()/.lookup_pos() (engine-internal), and the
# events_to_exclude_cats() utility.
#
# Column conventions (OWP native):
#   flowpath_id, flowpath_toid, lengthkm, hydroseq, levelpathid, [totdasqkm]

# -- reconcile_collapsed_flowlines ---------------------------------------

#' @title Reconcile collapsed flowlines
#' @description Equivalent to the reference reconcile-collapsed-flowlines step
#'   but data.table-based for speed. Returns sf or data.frame depending on whether
#'   `geom` is supplied.
#' @param flines data.frame output of the collapse step (`.reactor_collapse`).
#' @param geom sf data.frame with geometry column and id column for line union.
#' @param id character name of the ID column in geom (default "flowpath_id").
#' @return sf data.frame (if geom supplied) or data.frame with columns:
#'   ID, toID, lengthkm, member_flowpath_id (list), levelpathid, hydroseq.
#' @examples
#' \dontrun{
#' collapsed <- .reactor_collapse(flines)
#' rec <- reconcile_collapsed_flowlines(collapsed, geom = flowline_sf,
#'                                      id = "flowpath_id")
#' }
#' @export
reconcile_collapsed_flowlines <- function(flines, geom = NULL, id = "flowpath_id") {
  dt <- data.table::as.data.table(flines)
  dt[, flowpath_id        := hfutils::hf_fmt_id(flowpath_id)]
  dt[, flowpath_toid      := hfutils::hf_fmt_id(flowpath_toid)]
  if (!"joined_from_flowpath_id" %in% names(dt)) dt[, joined_from_flowpath_id := NA_character_]
  if (!"joined_to_flowpath_id"   %in% names(dt)) dt[, joined_to_flowpath_id   := NA_character_]

  # -- Step 1: compute becomes -----------------------------------------------
  # Only follow join pointers that stay within the same levelpathid.
  # Cross-levelpath joins (e.g. a short headwater tributary absorbed downstream
  # into a mainstem segment) produce non-contiguous MULTILINESTRING geometries
  # and incorrect member_flowpath_id groupings.
  .lp_by_fp <- stats::setNames(dt$levelpathid, dt$flowpath_id)
  dt[, becomes := {
    from    <- joined_from_flowpath_id
    to      <- joined_to_flowpath_id
    own_lp  <- levelpathid
    from_ok <- !is.na(from) & from != "-9999" &
      !is.na(.lp_by_fp[from]) & .lp_by_fp[from] == own_lp
    to_ok   <- !is.na(to)   & to   != "-9999" &
      !is.na(.lp_by_fp[to])   & .lp_by_fp[to]   == own_lp
    ifelse(from_ok, from, ifelse(to_ok, to, flowpath_id))
  }]

  # -- Step 2: fix "joined_from target also joins" (mirrors original logic)
  joined_from <- dt[!is.na(joined_from_flowpath_id) & joined_from_flowpath_id != "-9999" &
    joined_from_flowpath_id %in% dt$becomes]
  joined_to   <- dt[!is.na(joined_to_flowpath_id)   & joined_to_flowpath_id   != "-9999" &
    joined_to_flowpath_id   %in% dt$becomes]
  joined_tofrom <- joined_to[becomes %in% joined_from$flowpath_id]
  if (nrow(joined_tofrom) > 0L) {
    from_map <- joined_from[, .(becomes_key = flowpath_id, new_becomes = becomes)]
    upd      <- from_map[joined_tofrom, on = c("becomes_key" = "becomes"), nomatch = NA]
    upd      <- upd[!is.na(new_becomes), .(flowpath_id, new_becomes)]
    dt[upd, on = "flowpath_id", becomes := new_becomes]
  }

  # -- Step 3: per-becomes aggregate (max lengthkm, min hydroseq, min levelpathid)
  dt[, lengthkm   := max(lengthkm,   na.rm = TRUE), by = becomes]
  dt[, hydroseq   := min(hydroseq,   na.rm = TRUE), by = becomes]
  dt[, levelpathid := min(levelpathid, na.rm = TRUE), by = becomes]

  # -- Step 4: part detection (e.g. "12345.2" -> orig="12345", part="2") ------
  parts <- strsplit(dt$flowpath_id, "\\.", fixed = FALSE)
  dt[, orig_flowpath_id := vapply(parts, `[`, character(1), 1)]
  dt[, part       := vapply(parts,
    function(x) if (length(x) > 1L) x[2L] else "0", character(1))]

  # -- Step 5: representative flowpath_toid per becomes group ---------------------
  # Pick a flowpath_toid pointing OUTSIDE the becomes group.  If the first row
  # (after sort) has a toid that targets another member of the SAME group, the
  # toID would self-loop.  Filter to external toIDs first; groups where every
  # member's toid is intra-group are outlets (toID = NA).
  #
  # CRITICAL: Prioritise the becomes-representative row (flowpath_id == becomes).
  # Zero-length absorbed members retain stale toid pointers after de-alias (e.g.
  # a member absorbed via skip_past ends up pointing to an upstream group rather
  # than the true downstream). The becomes-representative was never skip_past'd and
  # always holds the authoritative downstream toid for the group.
  #
  # CRITICAL: Terminal toids ("0") are always external and take highest priority.
  # When a becomes group contains a COMID that is a true network outlet (original
  # toCOMID was 0 in the reference fabric), that group IS a terminal group, its
  # representative toid must be "0" regardless of any non-terminal external toids
  # carried by other members.  NA toids at this stage are unresolvable stale
  # references (not true outlets) and are excluded as candidates.
  member_set <- split(dt$flowpath_id, dt$becomes)
  # At this stage flowpath_toid is either:
  #   "0"  : genuine terminal (original toCOMID was 0 in reference fabric)
  #   NA   : unresolvable stale reference (NOT a true outlet)
  #   other, a valid flowpath_id pointing to the downstream segment
  # Only "0" is treated as terminal; NA stale artifacts are excluded entirely.
  is_external_toid <- vapply(seq_len(nrow(dt)), function(i) {
    tc <- dt$flowpath_toid[i]
    if (!is.na(tc) && tc == "0") return(TRUE)  # genuine terminal
    if (is.na(tc)) return(FALSE)               # stale artifact -- exclude
    !(tc %in% member_set[[dt$becomes[i]]])
  }, logical(1))
  dt_ext <- dt[is_external_toid, .(flowpath_id, becomes, flowpath_toid,
    is_rep      = (flowpath_id == becomes),
    is_terminal = !is.na(flowpath_toid) &
      flowpath_toid == "0")]
  # Sort: terminals first (highest priority), then becomes-representative
  data.table::setorder(dt_ext, becomes, -is_terminal, -is_rep)
  tc_rep <- dt_ext[, .(flowpath_toid = flowpath_toid[1L]), by = becomes]
  # Groups with no external toid candidates (all members point intra-group or
  # have only stale-artifact NAs) are genuine network outlets.
  missing_bec <- setdiff(unique(dt$becomes), tc_rep$becomes)
  if (length(missing_bec) > 0)
    tc_rep <- rbind(tc_rep,
      data.table::data.table(becomes = missing_bec, flowpath_toid = "0"))
  dt[, flowpath_toid := NULL]
  dt[tc_rep, on = "becomes", flowpath_toid := flowpath_toid]

  # -- Step 6: sequential integer ID per unique becomes ---------------------
  u_becomes <- unique(dt$becomes)
  becomes_id <- data.table::data.table(
    becomes = u_becomes,
    ID      = seq_along(u_becomes)
  )
  dt[becomes_id, on = "becomes", ID := ID]

  # -- Step 7: toID via flowpath_toid -> becomes -> ID ------------------------------
  bc_of <- dt[!is.na(becomes), .(flowpath_id, bc_of = becomes)]
  # map each row's flowpath_toid to its becomes
  dt[bc_of, on = c("flowpath_toid" = "flowpath_id"), to_becomes := bc_of]
  # CRITICAL: use i.ID to reference becomes_id$ID (the downstream group's ID),
  # NOT dt$ID (the row's own ID). data.table X[Y, ...] gives X-column priority
  # when both X and Y share a column name, causing self-loops if we write just ID.
  dt[becomes_id, on = c("to_becomes" = "becomes"), toID := i.ID]
  # toID = NA for outlets (matches the reference convention)

  # -- Step 8: collapse to one row per becomes (sorted like original) --------
  dt[, part_num := suppressWarnings(as.integer(part))]
  dt[is.na(part_num), part_num := 0L]
  data.table::setorder(dt, hydroseq, -part_num)

  out <- dt[, .(
    toID         = toID[1L],
    lengthkm     = lengthkm[1L],
    levelpathid   = levelpathid[1L],
    ID           = ID[1L],
    min_hs       = min(hydroseq),
    member_flowpath_id = list(flowpath_id)
  ), by = becomes]

  # -- Step 9: new sequential hydroseq (ID_Hydroseq in original) ------------
  data.table::setorder(out, min_hs)
  out[, ID_Hydroseq := .I]

  # Fix toID to use new ID_Hydroseq space
  id_map <- out[, .(ID_old = ID, ID_new = ID_Hydroseq)]
  out[id_map, on = c("toID" = "ID_old"), toID_new := ID_new]
  # toID_new = NA for outlets (reference convention)

  # -- Step 10: levelpathid -------------------------------------------------
  # Save the original NHD levelpathid before overwriting with new integer IDs.
  out[, orig_levelpathid := levelpathid]
  # CRITICAL: use a different column name (new_lp) in the j-expression to
  # avoid overwriting the by-key column.  data.table replaces the by-column
  # values with the j-expression values when they share a name, destroying
  # the grouping key and making the subsequent join fail (all levelpathid = NA).
  lp_outlet <- out[, .(new_lp = ID_Hydroseq[which.min(min_hs)]), by = levelpathid]
  out[lp_outlet, on = "levelpathid", levelpathid := i.new_lp]

  # -- Step 11: select output ------------------------------------------------
  out <- out[, .(
    ID           = ID_Hydroseq,
    toID         = toID_new,
    lengthkm,
    member_flowpath_id,
    levelpathid,
    hydroseq     = ID_Hydroseq,
    orig_levelpathID = orig_levelpathid
  )]

  # Assign stable unique terminal nexus IDs.
  # Each outlet (toID == NA) gets a distinct value of (1:n + 1e9) so that
  # multiple outlets sharing the same physical coordinate remain distinguishable
  # through collapse -> reconcile -> aggregate -> nextgen_from_gpkg without requiring
  # any geometry-based coordinate hashing.
  term_idx <- which(is.na(out$toID))
  if (length(term_idx) > 0L)
    out$toID[term_idx] <- seq_along(term_idx) + 1e9

  # -- Step 12: attach geometry ----------------------------------------------
  if (!is.null(geom)) {
    geom_col <- attr(geom, "sf_column")
    geom2 <- geom |>
      dplyr::mutate(member_flowpath_id = as.character(.data[[id]])) |>
      dplyr::select(member_flowpath_id, dplyr::all_of(geom_col))

    # Flatten member_flowpath_id lists
    members_flat <- data.table::data.table(
      ID           = rep(out$ID,        lengths(out$member_flowpath_id)),
      member_flowpath_id = unlist(out$member_flowpath_id)
    )

    geom_join <- dplyr::left_join(
      as.data.frame(members_flat),
      as.data.frame(geom2),
      by = "member_flowpath_id"
    ) |> sf::st_as_sf()

    non_empty <- !sf::st_is_empty(sf::st_geometry(geom_join))
    geom_for_union <- geom_join[non_empty, ]
    line_types     <- c("LINESTRING", "MULTILINESTRING")
    geom_for_union <- geom_for_union[
      as.character(sf::st_geometry_type(geom_for_union)) %in% line_types, ]

    # Fast path: split into singleton groups (no union needed) and multi-member
    # groups (need st_union + st_line_merge). For AK-scale inputs (~214k rows)
    # the vast majority are singletons; only touching the multi-member subset
    # avoids the O(n) dplyr group-summarise cost on the full set.
    gfu_dt  <- data.table::as.data.table(
      data.frame(ID = geom_for_union$ID,
        row = seq_len(nrow(geom_for_union))))
    counts  <- gfu_dt[, .N, by = ID]
    solo_ids  <- counts$ID[counts$N == 1L]
    multi_ids <- counts$ID[counts$N >  1L]

    geom_col  <- attr(sf::st_geometry(geom_for_union), "class")[1]
    g_sfc     <- hfutils::rename_geometry(geom_for_union, "geometry")

    # Singletons: grab directly
    solo_rows <- gfu_dt[ID %in% solo_ids, row]
    solo_sf   <- g_sfc[solo_rows, c("ID", "geometry")]

    # Multi-member: union then line-merge
    if (length(multi_ids) > 0L) {
      multi_rows <- gfu_dt[ID %in% multi_ids, row]
      multi_sub  <- g_sfc[multi_rows, c("ID", "geometry")]
      multi_sf   <- multi_sub |>
        dplyr::group_by(ID) |>
        dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
        sf::st_as_sf()
      is_multi <- sf::st_geometry_type(multi_sf) == "MULTILINESTRING"
      if (any(is_multi))
        sf::st_geometry(multi_sf)[is_multi] <-
          sf::st_line_merge(sf::st_geometry(multi_sf)[is_multi])
      geom_union <- rbind(solo_sf, multi_sf)
    } else {
      geom_union <- solo_sf
    }

    out_sf <- dplyr::left_join(as.data.frame(out), geom_union, by = "ID") |>
      sf::st_as_sf()
    # IDs whose every member was a degenerate type (all filtered above) have no
    # entry in geom_union -> NA geometry after the join. Replace with an empty
    # LINESTRING so sf::st_length() and downstream writes don't fail.
    na_geom <- which(is.na(sf::st_geometry(out_sf)))
    if (length(na_geom) > 0L) {
      empty_ls <- sf::st_sfc(
        rep(list(sf::st_linestring()), length(na_geom)),
        crs = sf::st_crs(out_sf))
      sf::st_geometry(out_sf)[na_geom] <- empty_ls
    }
    out_sf$lengthkm <- as.numeric(sf::st_length(out_sf)) / 1000
    return(out_sf)
  }

  as.data.frame(out)
}

# -- snap events to reconciled flowpath outlets -------------------------------
# Maps each event to its reconciled flowpath ID via member_flowpath_id, then
# snaps the event geometry to the nearest flowpath outlet within snap_dist_m.
# Returns an sf POINT with columns: poi_id, hl_link, hl_reference,
# reconciled_id, orig_flowpath_id, snap_dist_m, and geometry.

.snap_events_to_reconciled <- function(events, reconciled, snap_dist_m = 250,
                                       in_crs = NULL) {
  if (is.null(events) || nrow(events) == 0L) return(NULL)
  if (is.null(reconciled) || nrow(reconciled) == 0L) return(NULL)

  # Work in a projected CRS for distance calculations
  work_crs <- sf::st_crs(reconciled)
  if (sf::st_is_longlat(reconciled)) work_crs <- sf::st_crs(5070)
  rec_proj <- sf::st_transform(reconciled, work_crs)
  ev_proj  <- sf::st_transform(events, work_crs)

  # Exclude zero-length geometries: lwgeom::st_endpoint segfaults on them
  # (start == end, e.g. AK TDX source degenerate flowlines).
  nonzero <- as.numeric(sf::st_length(rec_proj)) > 0
  rec_proj <- rec_proj[nonzero, ]
  if (nrow(rec_proj) == 0L) return(NULL)

  # Build endpoints: all outlets first, then all inlets. When a POI sits at
  # the junction between two reconciled flowpaths, it is equidistant from the
  # upstream flowpath's outlet AND the downstream flowpath's inlet. Listing
  # outlets first makes `ties.method = "first"` prefer the outlet, which
  # matches the convention that `poi_id` on flowpath X means X's
  # downstream outlet is the POI location. Without this, aggregate_along_
  # mainstems marks the wrong flowpath as POI-bearing and collapse_headwaters
  # later merges across the intended break point.
  outlets <- lwgeom::st_endpoint(sf::st_geometry(rec_proj))
  inlets  <- lwgeom::st_startpoint(sf::st_geometry(rec_proj))
  n_rec   <- nrow(rec_proj)
  endpoints <- sf::st_sf(
    rec_id  = c(rec_proj$ID, rec_proj$ID),
    ep_type = c(rep("outlet", n_rec), rep("inlet", n_rec)),
    geometry = c(outlets, inlets)
  )

  # Member lookup (fallback only): vectorized strsplit instead of lapply.
  # Also build a stripped-base index so events keyed on the unsplit parent
  # reference COMID (e.g. "6128655") can still resolve to a reconciled group that
  # only lists its split children ("6128655.1", "6128655.2").
  rec_dt  <- data.table::as.data.table(sf::st_drop_geometry(reconciled))
  splits  <- strsplit(rec_dt$member_flowpath_id, ",", fixed = TRUE)
  member_lu <- data.table::data.table(
    member = trimws(unlist(splits)),
    rec_id = rep(rec_dt$ID, lengths(splits))
  )
  member_lu_base <- data.table::copy(member_lu)
  member_lu_base[, member := sub("\\.[0-9]+$", "", member)]
  # Multiple split children of the same parent DO land in different reconciled
  # groups, that is what a POI-driven split produces. Keeping the first-seen
  # mapping picked one arbitrarily with respect to flow direction, which is how a
  # POI ended up attributed to the reach it is the INLET of. All candidates are
  # kept here and the choice is made by outlet proximity in the resolver below.
  member_lu_base <- unique(member_lu_base, by = c("member", "rec_id"))

  ev_attrs_df <- sf::st_drop_geometry(ev_proj)
  ev_fp_id    <- as.character(ev_attrs_df[["flowpath_id"]])
  has_poi  <- "poi_id"          %in% names(ev_attrs_df)
  has_link <- "hl_link"         %in% names(ev_attrs_df)
  has_ref  <- "hl_reference"    %in% names(ev_attrs_df)
  has_orig <- "orig_flowpath_id" %in% names(ev_attrs_df)

  # Member lookup: resolve the hydrolocation's original reference comid to a reconciled
  # fp. Use orig_flowpath_id (preserved before spatial snap overwrote flowpath_id)
  # when available, it carries the NHD reach assignment from the reference fabric,
  # which is authoritative even when a different fp's endpoint is spatially closer.
  # Stale prefixed IDs (e.g. "fp-27586") won't match numeric members and fall
  # through to spatial snapping harmlessly.
  ev_orig_id <- if (has_orig) as.character(ev_attrs_df[["orig_flowpath_id"]]) else ev_fp_id
  member_resolved <- member_lu$rec_id[match(ev_orig_id, member_lu$member)]

  na_mr <- is.na(member_resolved)
  if (any(na_mr)) {
    # format() pads shorter IDs when the vector contains mixed-width values
    # (e.g. 7-digit COMIDs padded when 9-digit synthetic junction IDs are present).
    # trimws() removes the padding so the match against member_lu succeeds.
    ev_fp_num <- trimws(suppressWarnings(
      format(as.numeric(ev_fp_id[na_mr]), scientific = FALSE)))
    member_resolved[na_mr] <- member_lu$rec_id[match(ev_fp_num, member_lu$member)]
  }
  na_mr <- is.na(member_resolved)
  if (any(na_mr)) {
    member_resolved[na_mr] <- member_lu_base$rec_id[
      match(ev_orig_id[na_mr], member_lu_base$member)]
  }
  # Step 3b: orig_flowpath_id may be stale/absent while flowpath_id is a split
  # child (e.g. flowpath_id=24287444, member_lu_base has 24287444 from 24287444.2).
  na_mr <- is.na(member_resolved)
  if (any(na_mr)) {
    ev_fp_base <- sub("\\.[0-9]+$", "", ev_fp_id[na_mr])
    member_resolved[na_mr] <- member_lu_base$rec_id[
      match(ev_fp_base, member_lu_base$member)]
  }

  # -- Choose among a split parent's children by OUTLET PROXIMITY ---------------
  # A POI-driven split sends the children of one reference reach into different
  # reconciled groups, so a base-name lookup has several valid candidates and
  # `match()` returns whichever is first. Keeping every candidate in the table is
  # not enough on its own, the choice has to be made per event, or every POI on a
  # split parent inherits one arbitrary child.
  #
  # Several POIs along one parent then resolve to the same child, and the
  # shared-reach merge downstream treats them as collocated and absorbs them.
  #
  # A POI is the OUTLET of the reach it belongs to, so the candidate whose outlet
  # lies nearest the event is that reach. Applied to every resolved event whose
  # parent has more than one child, so it corrects the exact-match path as well as
  # the stripped-base path.
  {
    base_cand <- split(member_lu_base$rec_id, member_lu_base$member)
    bc_names  <- names(base_cand)
    rec_ix    <- stats::setNames(seq_len(nrow(rec_proj)), as.character(rec_proj$ID))
    out_xy    <- sf::st_coordinates(outlets)[, 1:2, drop = FALSE]
    ev_xy     <- sf::st_coordinates(sf::st_geometry(ev_proj))[, 1:2, drop = FALSE]
    keys      <- sub("\\.[0-9]+$", "", ev_orig_id)
    keys[is.na(keys) | !nzchar(keys)] <- NA_character_

    n_refined <- 0L
    for (k in which(!is.na(member_resolved) & !is.na(keys))) {
      j <- match(keys[k], bc_names)
      if (is.na(j)) next
      cand <- base_cand[[j]]
      if (length(cand) < 2L) next
      ii <- rec_ix[as.character(cand)]
      ok <- !is.na(ii)
      if (sum(ok) < 2L) next
      cand <- cand[ok]; ii <- ii[ok]
      d <- (ev_xy[k, 1] - out_xy[ii, 1])^2 + (ev_xy[k, 2] - out_xy[ii, 2])^2
      best <- cand[[which.min(d)]]
      if (!identical(as.character(best), as.character(member_resolved[k]))) {
        member_resolved[k] <- best
        n_refined <- n_refined + 1L
      }
    }
    message(sprintf(paste0(
      ".snap_events_to_reconciled: %d event(s) re-attributed among the split ",
      "children of their parent reach by outlet proximity"), n_refined))
  }

  # -- Hydrolocation placement rule, applied against the RECONCILED network ------
  # events.R applies this tree at event generation, but against the REFERENCE
  # network. Refactor and reconcile then re-cut and re-aggregate the reaches, so
  # the property has to be re-established here: mapping an event to a reconciled
  # reach and taking whatever comes out leaves POIs attributed to the reach they
  # are the INLET of.
  #
  # For each resolved event, measured ALONG the line:
  #   within snap_dist_m of the outlet  -> keep it here, move it to the outlet
  #   within snap_dist_m of the inlet   -> it is the OUTLET of an incoming reach;
  #                                        reassign there, preferring the mainstem
  #                                        inflow (same levelpath), then hydroseq
  #   otherwise                         -> it belongs mid-reach and must be a split
  #                                        event; counted and reported, not moved
  {
    .rid   <- as.character(rec_proj$ID)
    .ix_of <- stats::setNames(seq_len(nrow(rec_proj)), .rid)
    .lp    <- as.character(rec_proj$levelpathid)
    .hs    <- suppressWarnings(as.numeric(rec_proj$hydroseq))
    .to    <- as.character(rec_proj$toID)
    .ups_of <- split(seq_len(nrow(rec_proj)), .to)   # reaches draining INTO each id

    # Along-line position of a point: cumulative length to its nearest vertex.
    .along <- function(i_rec, pt) {
      m <- sf::st_coordinates(sf::st_geometry(rec_proj)[[i_rec]])
      d <- cumsum(c(0, sqrt(diff(m[, 1])^2 + diff(m[, 2])^2)))
      j <- which.min((m[, 1] - pt[1])^2 + (m[, 2] - pt[2])^2)
      c(pos = d[j], len = d[length(d)])
    }

    .n_out <- 0L; .n_up <- 0L; .n_split <- 0L
    .res_i <- which(!is.na(member_resolved))
    .ev_xy <- sf::st_coordinates(sf::st_geometry(ev_proj))
    for (k in .res_i) {
      i_rec <- .ix_of[[as.character(member_resolved[k])]]
      if (is.null(i_rec) || is.na(i_rec)) next
      pl <- .along(i_rec, .ev_xy[k, 1:2])
      d_out <- pl[["len"]] - pl[["pos"]]
      d_in  <- pl[["pos"]]

      # Decide by PROXIMITY, not by branch order.
      #
      # On a reach shorter than snap_dist_m both ends pass the tolerance test, so
      # testing the outlet first settled the question by ordering rather than by
      # where the POI actually is: a POI 0.0 m from the inlet of a 207 m reach was
      # scored "at the outlet" 200 m away and stayed attributed to the reach it is
      # the INLET of. The aggregation then honoured that attribution faithfully,
      # which is how a POI ended up interior to its own flowpath while every stage
      # reported success.
      if (d_out <= snap_dist_m && d_out <= d_in) { .n_out <- .n_out + 1L; next }

      if (d_in <= snap_dist_m) {
        cand <- .ups_of[[as.character(member_resolved[k])]]
        if (length(cand)) {
          # Mainstem inflow first: the one continuing this reach's levelpath.
          same_lp <- cand[!is.na(.lp[cand]) & .lp[cand] == .lp[i_rec]]
          if (length(same_lp)) cand <- same_lp
          # Then the most downstream by hydroseq, a stable stand-in for the
          # larger-order inflow when several share a levelpath.
          pick <- cand[which.min(ifelse(is.na(.hs[cand]), Inf, .hs[cand]))]
          member_resolved[k] <- rec_proj$ID[pick]
          .n_up <- .n_up + 1L
          next
        }
        # Nothing drains in, so no upstream outlet can own this location. Keeping
        # it here is still better than calling it mid-reach when the outlet is
        # within tolerance; a headwater POI that cannot reach an outlet is
        # withdrawn later, with the reason recorded on its hydrolocation.
        if (d_out <= snap_dist_m) { .n_out <- .n_out + 1L; next }
      }
      .n_split <- .n_split + 1L
    }
    message(sprintf(paste0(
      ".snap_events_to_reconciled: placement -- %d at outlet, %d reassigned to an ",
      "incoming reach's outlet, %d mid-reach (require a split event)"),
      .n_out, .n_up, .n_split))
  }

  rec_ids   <- rep(NA_integer_, nrow(ev_proj))
  snap_geom <- sf::st_geometry(ev_proj)
  snap_type <- rep(NA_character_,  nrow(ev_proj))
  dists_near <- rep(NA_real_, nrow(ev_proj))

  # Primary: member-resolved events: use the outlet endpoint of the resolved
  # reconciled fp as the snap geometry, honouring the NHD reach assignment.
  member_ok <- !is.na(member_resolved)
  if (any(member_ok)) {
    rec_ids[member_ok] <- member_resolved[member_ok]
    # Snap geometry to the outlet of the resolved fp (POI convention: outlet)
    outlet_idx <- match(member_resolved[member_ok], rec_proj$ID)
    valid_idx  <- !is.na(outlet_idx)
    rows_ok    <- which(member_ok)[valid_idx]
    snap_geom[rows_ok] <- outlets[outlet_idx[valid_idx]]
    snap_type[member_ok] <- "member"
    dists_near[member_ok] <- as.numeric(sf::st_distance(
      sf::st_geometry(ev_proj)[member_ok],
      snap_geom[member_ok],
      by_element = TRUE))
  }

  # Secondary: spatial snap for events with no member match (stale/missing fp_id).
  # Only compute the distance matrix for events that still need spatial resolution,
  # avoids an O(n_events x 2n_rec) allocation that is fatal at AK scale (212k fps).
  spatial_rows <- which(!member_ok)
  snapped      <- rep(FALSE, nrow(ev_proj))
  best_i       <- rep(NA_integer_, nrow(ev_proj))
  if (length(spatial_rows) > 0L) {
    dist_sub <- matrix(as.numeric(sf::st_distance(ev_proj[spatial_rows, ], endpoints)),
      nrow = length(spatial_rows))
    # Prefer an OUTLET outright rather than only on an exact tie. Outlets occupy the
    # first n_rec columns; a node is both an upstream reach's outlet and a downstream
    # reach's inlet, and the two coincident points can differ by rounding, which let
    # the inlet win and attributed the POI to the reach it is the inlet of. An inlet
    # is used only when no outlet lies within the snap distance, the legitimate case
    # being a POI at the very top of a headwater reach.
    .is_out <- seq_len(ncol(dist_sub)) <= n_rec
    .d_out  <- dist_sub; .d_out[, !.is_out] <- Inf
    bi_out  <- max.col(-.d_out, ties.method = "first")
    d_out_best <- .d_out[cbind(seq_along(bi_out), bi_out)]
    bi_sub  <- max.col(-dist_sub, ties.method = "first")
    .use_out <- is.finite(d_out_best) & d_out_best <= snap_dist_m
    bi_sub[.use_out] <- bi_out[.use_out]
    dn_sub   <- dist_sub[cbind(seq_along(bi_sub), bi_sub)]
    best_i[spatial_rows]     <- bi_sub
    dists_near[spatial_rows] <- dn_sub
    snapped[spatial_rows]    <- !is.na(dn_sub) & dn_sub <= snap_dist_m
  }

  spatial_needed <- !member_ok & snapped
  if (any(spatial_needed)) {
    rec_ids[spatial_needed]   <- endpoints$rec_id[best_i[spatial_needed]]
    snap_geom[spatial_needed] <- sf::st_geometry(endpoints)[best_i[spatial_needed]]
    snap_type[spatial_needed] <- endpoints$ep_type[best_i[spatial_needed]]
  }

  # Tertiary: beyond snap threshold and no member match: record as unsnapped
  unsnapped <- which(!member_ok & !snapped)
  if (length(unsnapped) > 0L) {
    snap_type[unsnapped] <- "unsnapped"
  }

  n_unsnapped <- sum(!member_ok & !snapped)
  if (n_unsnapped > 0L)
    .hf_log_warn(sprintf(
      ".snap_events_to_reconciled: %d POI event(s) dropped -- no member match and >%gm from any endpoint",
      n_unsnapped, snap_dist_m))

  keep <- !is.na(rec_ids)
  if (!any(keep)) return(NULL)

  n_keep <- sum(keep)
  out <- sf::st_sf(
    poi_id           = if (has_poi) ev_attrs_df$poi_id[keep]       else seq_len(nrow(ev_proj))[keep],
    hl_link          = if (has_link) ev_attrs_df$hl_link[keep]      else rep(NA_character_, n_keep),
    hl_reference     = if (has_ref) ev_attrs_df$hl_reference[keep] else rep(NA_character_, n_keep),
    reconciled_id    = rec_ids[keep],
    orig_flowpath_id = ev_fp_id[keep],
    snap_dist_m      = round(dists_near[keep], 1),
    snap_type        = snap_type[keep],
    geometry         = snap_geom[keep]
  )

  out_crs <- if (!is.null(in_crs) && !is.na(in_crs)) in_crs else work_crs
  sf::st_transform(out, out_crs)
}

# -- fast_split ---------------------------------------------------------------
# Wrapper around split_flowlines that adds furrr parallelism by splitting the
# network into independent level-path subsets and processing them in parallel.

.fast_split <- function(flines, max_length, events = NULL, avoid = NULL,
                        cores = 1L, min_split_m = 10) {
  use_para <- cores > 1L && requireNamespace("furrr", quietly = TRUE)
  if (!use_para) {
    return(split_flowlines(flines,
      max_length = max_length,
      events = events,
      avoid = avoid,
      min_split_m = min_split_m))
  }

  # Split network into independent chunks by levelpathid or connected components
  chunk_col <- if ("levelpathid" %in% names(flines)) "levelpathid" else NULL

  if (!is.null(chunk_col)) {
    chunks <- split(seq_len(nrow(flines)), flines[[chunk_col]])
  } else {
    # Fall back to sequential chunks of ~500 rows each
    idx_seq <- seq_len(nrow(flines))
    chunks  <- split(idx_seq, ceiling(idx_seq / 500))
  }

  results <- furrr::future_map(chunks, function(idx) {
    sub_sf <- flines[idx, , drop = FALSE]
    sub_ev <- if (!is.null(events)) {
      events[events$flowpath_id %in% sub_sf$flowpath_id, , drop = FALSE]
    } else NULL
    tryCatch(
      split_flowlines(sub_sf, max_length = max_length,
        events = sub_ev, avoid = avoid,
        min_split_m = min_split_m),
      error = function(e) sub_sf
    )
  }, .options = furrr::furrr_options(seed = NULL))

  do.call(rbind, results)
}

# -- Helper: null-coalescing operator -----------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -- Helper: vectorized flowpath_id -> row-index lookup ------------------------------
.lookup_pos <- function(comids, pos_map) {
  out <- pos_map[as.character(comids)]
  out[is.na(out)] <- 0L
  as.integer(out)
}

# -- events_to_exclude_cats ----------------------------------------------------

#' Build exclude_cats from event/split_event data
#'
#' Returns a character vector of segment IDs that should be protected from
#' collapse because they host (or are immediately downstream of) an event.
#' Pass the result directly as the `exclude_cats` argument to
#' [refactor_flowlines()].
#'
#' @param events sf or data.frame of events.
#'   Should have a `flowpath_id` column identifying the segment each event sits on.
#' @param split_events sf or data.frame of split events. Also expected to have
#'   `flowpath_id`.
#'   When both `events` and `split_events` are provided, IDs from both are used.
#' @param flowlines Optional sf or data.frame with flowpath topology columns
#'   (`flowpath_id`/`flowpath_id` and `flowpath_toid`/`flowpath_toid`). When supplied and
#'   `protect_downstream = TRUE`, the immediate downstream neighbor of each
#'   event-hosting segment is also included: matching the reference
#'   refactor-nhdplus behaviour which adds both the event flowpath_id and its
#'   `flowpath_toid` to `exclude_cats`.
#' @param protect_downstream Logical (default `TRUE`). Include the downstream
#'   (`toid`) neighbor of each event segment when `flowlines` is provided.
#' @param only_split_candidates Logical (default `FALSE`). When `TRUE` and
#'   `split_events` has a `split_candidate` column, restrict to rows where
#'   `split_candidate == TRUE`. Useful if only confirmed split-interior events
#'   should force protection.
#' @return Named character vector of unique segment IDs to exclude from collapse.
#' @examples
#' \dontrun{
#' excl <- events_to_exclude_cats(events    = poi_events,
#'                                flowlines = flines,
#'                                protect_downstream = TRUE)
#' refactor_flowlines(flines, split_flines_meters = 2000,
#'                    collapse_flines_meters = 1000,
#'                    out_refactored = "refactored.gpkg",
#'                    out_reconciled = "reconciled.gpkg",
#'                    exclude_cats = excl)
#' }
#' @export
events_to_exclude_cats <- function(events          = NULL,
                                   split_events     = NULL,
                                   flowlines        = NULL,
                                   protect_downstream    = TRUE,
                                   only_split_candidates = FALSE) {

  id_col <- function(df) {
    nms <- names(df)
    for (cn in c("flowpath_id", "comid", "id")) if (cn %in% nms) return(cn)
    stop("Cannot find an ID column in events/split_events data.")
  }

  collect_ids <- function(df) {
    if (is.null(df) || nrow(df) == 0L) return(character(0))
    if (inherits(df, "sf")) df <- sf::st_drop_geometry(df)
    col <- id_col(df)
    vals <- as.character(df[[col]])
    vals[!is.na(vals) & nzchar(vals)]
  }

  ev_ids  <- collect_ids(events)

  se_ids <- if (!is.null(split_events) && nrow(split_events) > 0L) {
    se <- if (inherits(split_events, "sf")) sf::st_drop_geometry(split_events) else split_events
    if (only_split_candidates && "split_candidate" %in% names(se)) {
      se <- se[isTRUE(se$split_candidate) | se$split_candidate == TRUE, , drop = FALSE]
    }
    collect_ids(se)
  } else character(0)

  protect <- unique(c(ev_ids, se_ids))
  protect <- protect[!is.na(protect) & nzchar(protect)]

  if (length(protect) == 0L) return(character(0))

  if (isTRUE(protect_downstream) && !is.null(flowlines) && nrow(flowlines) > 0L) {
    fl <- if (inherits(flowlines, "sf")) sf::st_drop_geometry(flowlines) else flowlines
    # resolve id/toid column names (support both NHD and OWP conventions)
    id_c  <- id_col(fl)
    tid_c <- {
      nms <- names(fl)
      for (cn in c("flowpath_toid", "toid")) if (cn %in% nms) {
        tid_c <- cn
        break
      }
      if (!exists("tid_c")) NA_character_ else tid_c
    }
    if (!is.na(tid_c)) {
      fl_ids  <- as.character(fl[[id_c]])
      fl_tids <- as.character(fl[[tid_c]])
      toid_map <- stats::setNames(fl_tids, fl_ids)
      downstream <- toid_map[protect]
      downstream <- downstream[!is.na(downstream) & nzchar(downstream) & downstream != "0"]
      protect <- unique(c(protect, as.character(downstream)))
    }
  }

  protect
}

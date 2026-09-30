# Flowpath-only variants of the network-aggregation internals.
# These strip all divide geometry work (union_polygons, terra::makeValid,
# st_snap, terra::relate) so network aggregation runs in seconds rather than
# minutes. Divide geometry is constructed separately in construct_agg_divides().

# -- .agg_sets_fp --------------------------------------------------------------
# Flowpath-only variant of the aggregate_sets step.
# Returns list(flowpaths = sf, lookup = data.frame, set_topo_fin = data.frame).
# No divide geometry is touched.
.agg_sets_fp <- function(network_list, index_table, rec_touch = NULL) {
  # -- Make every set internally connected (single-outlet) before contracting
  # cs_group groups flowpaths by levelpathid. A levelpath can be split by another
  # levelpath's reach, so the resulting set straddles a foreign set: two members
  # exit to two different sets and contraction yields a 2-cycle (set A->B and
  # B->A). Split each set into the weakly-connected components of its INTERNAL toid
  # subgraph: edges whose tail AND head are in the same set. Every resulting set
  # is then a single connected unit draining to one outlet, so its contraction is
  # acyclic. Components never merge across original sets (an internal edge requires
  # both endpoints in the same set), so this only ever sub-divides.
  #
  # An internal edge additionally requires the two members' DIVIDES to touch
  # spatially (`rec_touch` adjacency). A levelpath whose channel is routed
  # through a non-member reach's catchment is network-connected but its member
  # divides are spatially DISCONNECTED, so their union is a multi-part polygon
  # and construct_agg_divides drops the detached part, leaving an enclosed gap.
  # Dropping the non-adjacent edge isolates the stranded piece into its own set,
  # which the collapse then folds into the neighbour it actually touches. With
  # `rec_touch = NULL` the check is skipped (network-only, legacy behaviour).
  {
    .ids   <- as.character(index_table$flowpath_id)
    .toids <- as.character(index_table$flowpath_toid)
    .set   <- as.character(index_table$set)
    .setof <- stats::setNames(.set, .ids)
    .internal <- !is.na(.setof[.toids]) & .setof[.toids] == .set
    if (!is.null(rec_touch) && !nzchar(Sys.getenv("HF_NO_CONTIG_SPLIT"))) {
      .adj <- .internal
      .wi  <- which(.internal)
      if (length(.wi)) {
        # Sever only on POSITIVE evidence of spatial non-adjacency: both member
        # reaches carry a divide (present in rec_touch) AND the downstream divide
        # is not a neighbour. If EITHER endpoint is divide-less (absent from
        # rec_touch, a connector / waterbody reach with no catchment), the
        # touch test is meaningless, so keep the edge: severing a divide-less
        # connector strands it into an isolated divide-less leaf that Stage-3
        # demotion cannot re-home (no catchment to fold into), tripping
        # no_divideless_leaf_flowpaths.
        .adj[.wi] <- vapply(.wi, function(k) {
          nb   <- rec_touch[[.ids[k]]]
          nbdn <- rec_touch[[.toids[k]]]
          if (is.null(nb) || is.null(nbdn)) return(TRUE)   # divide-less: keep edge
          .toids[k] %in% nb
        }, logical(1L))
        .n_sever <- sum(.internal & !.adj)
        if (.n_sever > 0L)
          message(sprintf(
            "[agg_sets] severed %d spatially-disconnected internal edge(s) (divides do not touch)",
            .n_sever))
      }
      .internal <- .adj
    }
    .g <- igraph::add_vertices(igraph::make_empty_graph(directed = FALSE),
      length(.ids), name = .ids)
    if (any(.internal))
      .g <- igraph::add_edges(.g, c(rbind(.ids[.internal], .toids[.internal])))
    .comp <- igraph::components(.g, mode = "weak")$membership
    n_before <- length(unique(.set))
    new_set  <- as.integer(.comp[.ids])
    n_after  <- length(unique(new_set))
    if (n_after > n_before) {
      message(sprintf(
        "[agg_sets] split %d non-contiguous set(s) into connected components (%d -> %d sets)",
        n_after - n_before, n_before, n_after))

      # PROBE: classify the SHAPE of each split. Walk each set's exit edges forward
      # through foreign reaches (<=K hops). If the path re-enters the SAME set it's
      # a "detour" (mainstem routed through foreign reaches); if it never reconnects
      # the set is "disconnected", a levelpathid spanning physically separate
      # streams (e.g. levelpath over-assignment or VPU-boundary severing).
      .toid_of <- stats::setNames(.toids, .ids)
      .lp_of   <- stats::setNames(as.character(index_table$levelpathid), .ids)
      .ncomp   <- tapply(.comp[.ids], .set, function(z) length(unique(z)))
      .split_s <- names(.ncomp)[.ncomp > 1L]
      .K <- 200L
      .hop <- integer(0)
      .egs <- character(0)
      for (.s in .split_s) {
        .mem <- .ids[.set == .s]
        .ext <- unique(.toid_of[.mem][.toid_of[.mem] %in% .ids & !(.toid_of[.mem] %in% .mem)])
        .best <- NA_integer_
        .via <- character(0)
        for (.x in .ext) {
          .cur <- .x
          .lps <- character(0)
          for (.h in seq_len(.K)) {
            if (is.na(.cur) || .cur == "0" || !(.cur %in% .ids)) break
            if (.cur %in% .mem) {
              if (is.na(.best) || .h < .best) {
                .best <- .h
                .via <- .lps
              }
              break
            }
            .lps <- c(.lps, .lp_of[.cur])
            .cur <- .toid_of[.cur]
          }
        }
        .hop <- c(.hop, .best)
        if (!is.na(.best) && length(.egs) < 4L)
          .egs <- c(.egs, sprintf("lp %s %d-hop via {%s}", .lp_of[.mem[1]], .best,
            paste(unique(.via), collapse = ",")))
      }
      .lab <- ifelse(is.na(.hop), "disconnected",
        ifelse(.hop <= 2L, "detour(<=2hop)", "long(>2hop)"))
      .tab <- sort(table(.lab), decreasing = TRUE)
      .maxlong <- if (any(!is.na(.hop) & .hop > 2L)) max(.hop[!is.na(.hop) & .hop > 2L]) else 0L
      message(sprintf("[agg_sets] split shape (%d, K=%d): %s; max long-detour=%d hop%s",
        length(.split_s), .K,
        paste(sprintf("%s=%d", names(.tab), as.integer(.tab)), collapse = " "), .maxlong,
        if (length(.egs)) paste0("; e.g. ", paste(.egs, collapse = " | ")) else ""))
      rm(.toid_of, .lp_of, .ncomp, .split_s, .K, .hop, .egs, .lab, .tab, .maxlong)
    }
    index_table$set <- new_set
    # `n` (set size) is consumed below to route single- vs multi-fp geometry
    # unions, so recompute it for the new set ids.
    index_table <- dplyr::ungroup(dplyr::mutate(
      dplyr::group_by(index_table, set), n = dplyr::n()))
  }

  set_topo <- index_table |>
    dplyr::group_by(set) |>
    dplyr::mutate(
      member_comid = paste(member_comid, collapse = ","),
      poi_id       = paste(poi_id[!is.na(poi_id)], collapse = ","),
      poi_id       = ifelse(poi_id == "", NA_character_, poi_id)) |>
    dplyr::arrange(hydroseq) |>
    dplyr::select(set, flowpath_id, flowpath_toid, levelpathid,
      poi_id, hydroseq, member_comid) |>
    dplyr::ungroup()

  set_topo_fin <- dplyr::left_join(
    dplyr::select(set_topo, set, flowpath_id = flowpath_toid,
      hydroseq, levelpathid, poi_id, member_comid),
    dplyr::select(set_topo, toset = set, flowpath_id),
    by = "flowpath_id") |>
    dplyr::group_by(set) |>
    dplyr::mutate(toset = ifelse(is.na(toset), 0L, toset)) |>
    dplyr::filter(set != toset) |>
    dplyr::slice_min(hydroseq) |>
    dplyr::ungroup() |>
    dplyr::select(set, toset, levelpathid, poi_id, member_comid)

  # -- GUARD: the set->toset contraction must be acyclic ---------------------
  # The internal-component split above guarantees single-outlet sets, so this
  # should always hold. If a cycle ever survives, fail loudly with the offending
  # SCC's sets + member fps (toid/hydroseq) rather than letting the downstream
  # accumulate abort with an opaque "cannot accumulate".
  {
    .e <- set_topo_fin[set_topo_fin$toset != 0L, c("set", "toset")]
    if (nrow(.e) > 0L) {
      .g  <- igraph::graph_from_data_frame(
        data.frame(from = as.character(.e$set), to = as.character(.e$toset)),
        directed = TRUE)
      .sc <- igraph::components(.g, mode = "strong")
      .bad <- which(.sc$csize > 1L)
      if (length(.bad) > 0L) {
        for (.b in utils::head(.bad, 3L)) {
          .sets <- names(.sc$membership)[.sc$membership == .b]
          msg <- sprintf("[agg_sets] CYCLE among %d set(s): %s",
            length(.sets), paste(utils::head(.sets, 8), collapse = " <-> "))
          for (.s in utils::head(.sets, 4L)) {
            .mem <- index_table[as.character(index_table$set) == .s,
              c("flowpath_id", "flowpath_toid", "hydroseq", "levelpathid")]
            msg <- paste0(msg, sprintf("\n  set %s (lp %s): %s", .s,
              paste(unique(.mem$levelpathid), collapse = ","),
              paste(sprintf("%s->%s(hs%s)", .mem$flowpath_id,
                .mem$flowpath_toid, .mem$hydroseq),
              collapse = " ")))
          }
          .hf_log_warn(msg)
        }
      }
    }
  }

  # union flowpath geometries per set
  single_fp <- tryCatch(
    dplyr::filter(index_table, n == 1L) |>
      dplyr::inner_join(dplyr::select(network_list$flowpaths, flowpath_id),
        by = "flowpath_id") |>
      sf::st_as_sf() |>
      dplyr::select(set),
    error = function(e) NULL)

  multi_fp <- tryCatch(
    dplyr::filter(index_table, n > 1L) |>
      dplyr::inner_join(network_list$flowpaths, by = "flowpath_id") |>
      sf::st_as_sf() |>
      dplyr::select(set) |>
      hfutils::union_linestrings("set"),
    error = function(e) NULL)

  flowpaths_out <- dplyr::bind_rows(multi_fp, single_fp) |>
    dplyr::left_join(set_topo_fin, by = "set") |>
    dplyr::rename(flowpath_id = set, flowpath_toid = toset) |>
    dplyr::mutate(flowpath_id   = as.character(flowpath_id),
      flowpath_toid = as.character(flowpath_toid))

  # Carry incremental divide area onto each aggregate before the topology-only
  # prepare_network() pass below. The stub divides have empty geometries, so
  # area cannot be recovered there and would otherwise be replaced by zero.
  area_by_set <- index_table |>
    dplyr::select(set, flowpath_id) |>
    dplyr::mutate(set = as.character(set)) |>
    dplyr::mutate(flowpath_id = as.character(flowpath_id)) |>
    dplyr::left_join(
      sf::st_drop_geometry(network_list$flowpaths) |>
        dplyr::mutate(flowpath_id = as.character(flowpath_id)) |>
        dplyr::select(flowpath_id, areasqkm),
      by = "flowpath_id") |>
    dplyr::group_by(set) |>
    dplyr::summarise(areasqkm = sum(areasqkm, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(flowpath_id = set)
  flowpaths_out <- flowpaths_out |>
    dplyr::left_join(area_by_set, by = "flowpath_id")

  # stub divide table: divide_id + flowpath_toid only, no geometry
  div_stub <- dplyr::select(set_topo_fin, divide_id = set, flowpath_toid = toset) |>
    dplyr::mutate(divide_id     = as.character(divide_id),
      flowpath_toid = as.character(flowpath_toid))

  # call prepare_network with stub divides (geometry-free) so topology is valid
  stub_sf <- sf::st_sf(
    divide_id = div_stub$divide_id,
    geometry  = sf::st_sfc(
      lapply(seq_len(nrow(div_stub)), function(i) sf::st_geometrycollection()),
      crs = sf::st_crs(network_list$flowpaths)))
  stub_sf$flowpath_toid <- div_stub$flowpath_toid

  nl_out <- prepare_network(
    list(flowpaths = flowpaths_out, divides = stub_sf))

  lookup <- dplyr::mutate(
    dplyr::rename(
      dplyr::summarise(
        dplyr::group_by(index_table, set),
        # Membership is this set's own rows plus anything they absorbed before
        # aggregation, so a collapsed reach's divide is still found by rec_id.
        flowline_id   = paste(unique(c(as.character(flowpath_id),
          expand_member_csv(absorbed_recs, unique = TRUE))), collapse = ","),
        flowline_toid = paste(flowpath_toid, collapse = ","),
        member_comid  = paste(member_comid,  collapse = ","),
        .groups = "drop"),
      flowpath_id = set),
    flowpath_id = as.character(flowpath_id))

  nl_out$lookup       <- lookup
  nl_out$set_topo_fin <- set_topo_fin
  nl_out
}

# -- .build_hw_collapse_fp -----------------------------------------------------
# Flowpath-only variant of the build_headwater_collapse step.
# Detects headwaters using graph topology only, no terra::relate, no divide sf.
# Cross-round counters for the collapse. A data.frame attribute does not survive
# the subsetting and dplyr pipeline the candidate table goes through, so the
# guard leaves its count here for .collapse_hw_fp() to report once at the end.
.hf_collapse_stats <- new.env(parent = emptyenv())

.build_hw_collapse_fp <- function(network_list, min_area_sqkm, min_length_km,
                                  rec_touch = NULL, rec_area = NULL,
                                  ideal_size_sqkm = NULL, connector_ids = NULL,
                                  max_area_sqkm = NULL) {
  fp <- sf::st_drop_geometry(network_list$flowpaths)

  toids     <- as.character(fp$flowpath_toid)
  all_ids   <- as.character(fp$flowpath_id)

  # inflow count per fp: 0 = headwater, 1 = single-inflow mid-reach, >=2 = confluence
  inflow_n  <- tabulate(match(toids[toids %in% all_ids], all_ids),
    nbins = length(all_ids))

  # IDEAL-SIZE cross-levelpath collapse (principled mode, ideal_size + rec_area
  # given): fold both headwaters (0 inflow) and single-inflow mid-reaches (1
  # inflow) DOWNSTREAM until the accumulated unit reaches ideal_size. cs_group
  # only aggregates WITHIN a levelpath, so the cross-levelpath merge of short
  # tributaries is the collapse's job. Two requirements keep this safe:
  #   (1) area is member-summed from rec_div (rec_area), recomputed each round,
  #       so a grown unit exits at ideal_size, no unbounded cascade.
  #   (2) mid-reaches collapse STRICTLY to their toid-downstream (acyclic). Only
  #       true headwaters may take the lateral fallback (no inflow -> no edge
  #       redirected -> no cycle).
  # Confluences (>=2 inflow) are never sources, so junctions are preserved.
  # LEGACY mode (no ideal_size/rec_area): headwaters only, stale fp$areasqkm,
  # min_area/min_length floor.
  ideal_mode <- !is.null(ideal_size_sqkm) && !is.null(rec_area)
  is_cand    <- if (ideal_mode) inflow_n <= 1L else inflow_n == 0L

  toid_map  <- stats::setNames(toids, all_ids)
  hw_ids    <- all_ids[is_cand]
  touches   <- toid_map[hw_ids]

  # `is_hw` (true headwater, no inflow) was carried only to gate the lateral
  # re-home, which is gone; collapse now always targets the network downstream.
  df <- data.frame(
    flowpath_id  = hw_ids,
    becomes      = touches,
    stringsAsFactors = FALSE)

  # member-summed (true aggregate-divide) area when in ideal mode
  if (ideal_mode) {
    fl_all   <- strsplit(as.character(fp$flowline_id), ",", fixed = TRUE)
    area_all <- stats::setNames(
      vapply(fl_all, function(m) sum(rec_area[trimws(m)], na.rm = TRUE), numeric(1L)),
      all_ids)
  }

  # filter to small, unprotected candidates
  fp_idx <- match(df$flowpath_id, all_ids)
  df$poi_id    <- fp$poi_id[fp_idx]
  df <- df[!is.na(df$becomes) & df$becomes != "0" &
    df$becomes != df$flowpath_id, ]
  if (ideal_mode) {
    df$areasqkm  <- area_all[df$flowpath_id]
    df$recv_area <- area_all[df$becomes]
    # Two modes, matching the two phases the OG aggregation applies within a
    # levelpath (.assign_id, then .pinch_sides / .middle_massage):
    #
    #   floor  a unit under min_area_sqkm is offered a merge first, because
    #          .cs_group cannot clear the floor across levelpaths, its repair
    #          helpers return early on a single-group chunk, which is the
    #          reason this collapse exists at all. It gets first claim on a
    #          receiver's headroom, but not a licence to exceed it.
    #   target a unit already past the floor merges only while the result stays
    #          at or under ideal_size_sqkm. This is .assign_id's rule: the
    #          element that would carry a group past the target starts a new
    #          group rather than enlarging this one.
    #
    # Testing the donor alone, as this did, admits two 9.9 km2 units merging into
    # 19.8. That is where a median 30% above target and a fifth of all units past
    # 20 km2 come from, and why the overshoot tracks stage-1 catchment size: the
    # bigger the qualifying donor, the bigger the jump.
    # One rule, applied here and again per receiver below: a merge happens only
    # while the result stays at or under the target. Being under the floor buys
    # a donor priority for the headroom that exists, not an exemption from it.
    # Three sizes, each doing one job:
    #   min_area_sqkm    a unit below it must be absorbed if anything can take it
    #   ideal_size_sqkm  what an ordinary merge aims at and stops at
    #   max_area_sqkm    the ceiling absorbing a sub-floor unit may reach, so
    #                    clearing the floor cannot manufacture an oversized unit
    #
    # Without the ceiling the two criteria cannot both be served: an
    # unconditional floor merge lets a receiver already at the target keep
    # absorbing (10 -> 12.9 -> 15.8 -> 18.7), while a hard stop at the target
    # strands every sliver whose only neighbour is full.
    amax <- if (is.null(max_area_sqkm)) Inf else max_area_sqkm
    df$below_floor <- df$areasqkm < min_area_sqkm
    total       <- df$areasqkm + df$recv_area
    fits_target <- is.finite(df$recv_area) & total <= ideal_size_sqkm
    # A unit under min_area_sqkm MUST merge. The floor is a mandate and the
    # ceiling yields to it: a unit below the floor is not a valid computational
    # unit at all, while one above the ceiling is a valid unit of the wrong size.
    # Refusing the rescue to protect the ceiling trades a real failure for a
    # a small one, and the donor has nowhere else to go: there is one
    # downstream receiver, so a refusal is permanent.
    #
    # max_area_sqkm therefore bounds ordinary merges only. Floor merges that pass
    # it are counted and reported rather than blocked, so the cost of the mandate
    # stays visible.
    df <- df[(fits_target | df$below_floor) & is.na(df$poi_id), ]
    if (nrow(df)) {
      over <- df$below_floor & is.finite(df$recv_area) &
        (df$areasqkm + df$recv_area) > amax
      if (any(over)) message(sprintf(paste0("  [floor-mandate] %d sub-%g km2 ",
        "unit(s) merged past the %g km2 ceiling: no receiver had room"),
        sum(over), min_area_sqkm, amax))
    }
  } else {
    df$areasqkm <- fp$areasqkm[fp_idx]
    df$lengthkm <- fp$lengthkm[fp_idx]
    df <- df[(df$areasqkm < min_area_sqkm | df$lengthkm < min_length_km) &
      is.na(df$poi_id), ]
  }
  df <- df[df$becomes %in% all_ids, ]
  # CONNECTOR GUARD: never collapse across a catchment <-> connector boundary
  # (a unit is a connector iff ALL its members are reference-divide-less reaches).
  if (length(connector_ids) && nrow(df)) {
    fl_m  <- stats::setNames(strsplit(as.character(fp$flowline_id), ",", fixed = TRUE), all_ids)
    is_cn <- function(id) {
      m <- trimws(fl_m[[id]])
      !is.null(m) && length(m) && all(m %in% connector_ids)
    }
    df <- df[vapply(df$flowpath_id, is_cn, logical(1L)) ==
      vapply(df$becomes,     is_cn, logical(1L)), , drop = FALSE]
    if (nrow(df) == 0L) return(df)
  }
  # NO-CHAIN GUARD: a candidate whose target is itself collapsing this round
  # would chain (A->B while B->C), so defer it; it resolves next round.
  df <- df[!df$becomes %in% df$flowpath_id, ]
  if (nrow(df) == 0L) return(df)

  # CONTIGUITY REQUIREMENT: only collapse H into its network downstream when the
  # two CATCHMENTS touch. A divide is the area draining to its flowpath's outlet,
  # so if H drains into F then H union F is the area draining to F's outlet,
  # which for a real basin is contiguous. Merging two disjoint catchments
  # therefore manufactures a catchment shape that cannot exist, and the
  # aggregate's channel would have to leave its own divide, cross another, and
  # come back.
  #
  # When they do not touch, DECLINE the collapse: H stays its own flowpath with
  # its own contiguous divide. It sits below the size floor, which is a
  # preference, not a correctness property.
  #
  # A headwater is never folded sideways into a merely-touching neighbour. The
  # move is safe only at aggregate grain, where a headwater has no inflow and so
  # redirects no edge; at member grain its reaches still drain out of the
  # receiving flowpath, which gives that flowpath a second outlet and can close
  # a cycle (fp-A -> fp-B -> fp-A). The area argument for it does not hold
  # either: hfutils::union_polygons preserves every part, and multi-part divides
  # are legitimate in construct_agg_divides.
  #
  # Measured on VPU 18 against the delivered build, three ways:
  #                          flowpaths  2nd outlets  disjoint divides  lobe km2
  #   delivered (re-home)       26,866           53                11      17.1
  #   collapse anyway           26,876            0                30      66.0
  #   decline (this)            26,898            0                10      23.1
  # Declining matches the best topology and keeps disjoint divides below the
  # delivered build's, for +32 flowpaths (+0.12%) and +19 units under 3 km2.
  if (!is.null(rec_touch) && nrow(df)) {
    fl_split <- strsplit(as.character(fp$flowline_id), ",", fixed = TRUE)
    members  <- stats::setNames(lapply(fl_split, trimws), all_ids)
    keep <- vapply(seq_len(nrow(df)), function(i) {
      hm <- members[[df$flowpath_id[i]]]
      bm <- members[[df$becomes[i]]]
      if (is.null(hm) || is.null(bm)) return(FALSE)
      any(unique(unlist(rec_touch[hm], use.names = FALSE)) %in% bm)
    }, logical(1L))
    n_decl <- sum(!keep)
    if (n_decl > 0L) message(sprintf(
      "[collapse_hw] declined %d of %d collapse(s): catchments not contiguous",
      n_decl, length(keep)))
    df <- df[keep, , drop = FALSE]
    if (nrow(df) == 0L) return(df)
  }

  # -- LEVELPATH-CONTIGUITY GUARD ---------------------------------------------
  # Do not fold a reach DOWNSTREAM across a levelpath boundary while a surviving
  # upstream flowpath on the SAME levelpath still drains into it. Each aggregate
  # flowpath's geometry is later rebuilt as its mainstem walk
  # (extend_flowpaths_to_mainstem), which keeps only members sharing the
  # aggregate's mainstem levelpath. Folding a levelpath-outlet reach into a
  # foreign-levelpath aggregate therefore buries that reach OFF the target's
  # mainstem line, and the upstream same-levelpath flowpath is left draining onto
  # geometry that is absent from the flowpath network, a non-terminal outlet
  # that touches no flowpath. Deferring the fold keeps the levelpath a contiguous
  # chain; the reach may still collapse in a later round once its same-levelpath
  # upstream has itself merged in.
  if (nrow(df) && "levelpathid" %in% names(fp)) {
    # Vectorized levelpath lookup (match() once, not a per-row named-vector scan).
    lp_vec  <- as.character(fp$levelpathid)          # positional, aligned to all_ids
    lpX     <- lp_vec[match(df$flowpath_id, all_ids)]
    lpB     <- lp_vec[match(df$becomes,     all_ids)]
    cross   <- !is.na(lpX) & !is.na(lpB) & lpX != lpB  # only cross-lp folds can orphan
    keep_lp <- rep(TRUE, nrow(df))
    if (any(cross)) {
      ups_by <- split(all_ids, as.character(fp$flowpath_toid))  # fp -> inflow fp(s)
      into_x <- split(df$flowpath_id, df$becomes)               # becomes -> folders
      for (i in which(cross)) {
        X  <- df$flowpath_id[i]
        up <- ups_by[[X]]
        up <- up[!is.na(up) & up != X]
        if (!length(up)) next
        lp_up      <- lp_vec[match(up, all_ids)]
        same_lp_up <- up[!is.na(lp_up) & lp_up == lpX[i]]
        # a same-lp upstream folding INTO X this round merges cleanly (no orphan);
        # any other same-lp upstream survives and would be stranded -> defer.
        keep_lp[i] <- length(setdiff(same_lp_up, into_x[[X]])) == 0L
      }
    }
    n_def <- sum(!keep_lp)
    if (n_def)
      message(sprintf(
        "  [levelpath-guard] deferred %d cross-levelpath collapse(s) to preserve upstream contiguity",
        n_def))
    # The caller reports whether this ever clears. A deferral is only meant to
    # last until the blocking same-levelpath upstream merges in; when that
    # upstream carries a POI it never does, and the same candidates are deferred
    # every round forever. Reported once at the end rather than looking like
    # progress each round.
    .hf_collapse_stats$deferred_levelpath <- n_def
    df <- df[keep_lp, , drop = FALSE]
    if (nrow(df) == 0L) return(df)
  }

  # -- PER-RECEIVER BUDGET ----------------------------------------------------
  # Every accepted candidate is applied in the SAME round, and the block below
  # groups by `becomes` precisely so one receiver can take several donors at
  # once. Testing each donor on its own against the receiver's pre-merge area
  # therefore does not bound the result: at a confluence two 4 km2 donors each
  # pass against a 5 km2 receiver (9 <= 10) and leave it at 13.
  #
  # So spend the receiver's headroom once. Floor-mode donors are taken first and
  # unconditionally, because the floor outranks the target; the rest are taken
  # smallest-first while the running total still fits. A donor that does not fit
  # is not lost, only deferred, the caller iterates to a fixed point and
  # re-derives areas each round, so it is reconsidered against the receiver's
  # real size next time.
  if (ideal_mode && nrow(df) > 1L) {
    ord <- order(df$becomes,
      !df$below_floor,   # sub-floor donors get first claim on the headroom
      df$areasqkm)       # then smallest first
    df  <- df[ord, , drop = FALSE]
    run <- df$recv_area[!duplicated(df$becomes)]
    names(run) <- df$becomes[!duplicated(df$becomes)]
    keep <- logical(nrow(df))
    for (i in seq_len(nrow(df))) {
      b <- df$becomes[i]; a <- df$areasqkm[i]
      # One rule for every donor: a merge that would carry the receiver past the
      # target does not happen. Sub-floor donors are ordered first so they get
      # the first claim on the headroom, but they do not get to ignore it, an
      # unconditional floor merge lets a receiver that is already at the target
      # keep absorbing, and 10 -> 12.9 -> 15.8 -> 18.7 is how a unit ends up over
      # 20 with five members and none of them over 10.
      #
      # A donor that cannot be placed is deferred, not lost: the caller iterates
      # to a fixed point. One that is never placeable stays below the floor,
      # which is a size preference; building a 20 km2 unit to avoid it trades a
      # small miss for a large one.
      # A sub-floor donor is taken whatever the running total: the floor is a
      # mandate. Ordinary donors take only what headroom is left under the target.
      if (df$below_floor[i] ||
          (is.finite(run[[b]]) && run[[b]] + a <= ideal_size_sqkm)) {
        keep[i] <- TRUE
        run[[b]] <- run[[b]] + a
      }
    }
    n_def <- sum(!keep)
    if (n_def) message(sprintf(
      "  [budget] deferred %d collapse(s): receiver had no headroom left this round",
      n_def))
    df <- df[keep, , drop = FALSE]
    if (nrow(df) == 0L) return(df)
  }

  # build member_comid and flowline_id aggregations
  fp_mc <- stats::setNames(as.character(fp$member_comid), all_ids)
  fp_fl <- stats::setNames(as.character(fp$flowline_id),  all_ids)
  fp_hl <- stats::setNames(as.character(fp$poi_id),       all_ids)

  df$mC1 <- fp_mc[df$flowpath_id]
  df$mC2 <- fp_mc[df$becomes]
  df$fl1 <- fp_fl[df$flowpath_id]
  df$fl2 <- fp_fl[df$becomes]
  df$hl2 <- fp_hl[df$becomes]

  suppressWarnings(
    dplyr::ungroup(dplyr::mutate(
      dplyr::group_by(df, becomes),
      member_comid = paste0(mC2[1L], ",", paste(mC1, collapse = ",")),
      fl_id        = paste0(fl2[1L], ",", paste(fl1, collapse = ",")),
      poi_id       = as.numeric(paste(unique(stats::na.omit(c(hl2, poi_id))),
        collapse = ",")))) |>
      dplyr::select(flowpath_id, becomes, member_comid, poi_id, fl_id))
}

# -- .collapse_hw_fp -----------------------------------------------------------
# Flowpath-only variant of the collapse_headwaters step.
# Iteratively collapses small headwaters into their downstream fp using only
# flowpath topology, no divide geometry, no terra, no union_polygons.
.collapse_hw_fp <- function(network_list, min_area_sqkm, min_length_km,
                            rec_touch = NULL, rec_area = NULL,
                            ideal_size_sqkm = NULL, connector_ids = NULL,
                            max_area_sqkm = NULL) {
  message("--- Collapse Network Inward (fp-only) ---")

  # ensure flowline_id present (added by aggregate_along_mainstems lookup)
  if (!"flowline_id" %in% names(network_list$flowpaths)) {
    if (!is.null(network_list$lookup)) {
      network_list$flowpaths <- dplyr::left_join(
        network_list$flowpaths,
        dplyr::select(network_list$lookup, flowpath_id, flowline_id),
        by = "flowpath_id")
    } else {
      network_list$flowpaths$flowline_id <- as.character(
        network_list$flowpaths$flowpath_id)
    }
  }

  .hf_collapse_stats$deferred_levelpath <- 0L
  # The guard only runs while candidates survive to reach it, so the final round
  # of a converged loop reports nothing. What matters is the count still standing
  # when the loop stopped making progress, so keep the high-water mark.
  n_stuck <- 0L
  start <- nrow(network_list$flowpaths)
  mapping_table <- .build_hw_collapse_fp(network_list, min_area_sqkm, min_length_km, rec_touch, rec_area, ideal_size_sqkm, connector_ids, max_area_sqkm)
  count <- 0L

  while (nrow(mapping_table) > 0L) {
    count <- count + 1L
    message(sprintf("  Collapsing: %d features (round %d)",
      nrow(mapping_table), count))

    fp  <- network_list$flowpaths
    ind <- match(mapping_table$becomes, as.character(fp$flowpath_id))

    fp$member_comid[ind] <- mapping_table$member_comid
    fp$flowline_id[ind]  <- mapping_table$fl_id

    # remap toids: anything pointing to a collapsed fp now points to its becomes
    old_ids  <- as.character(mapping_table$flowpath_id)
    new_ids  <- as.character(mapping_table$becomes)
    toids    <- as.character(fp$flowpath_toid)
    remap_i  <- match(toids, old_ids)
    do_remap <- !is.na(remap_i)
    if (any(do_remap))
      fp$flowpath_toid[do_remap] <- new_ids[remap_i[do_remap]]

    # drop collapsed fps
    fp <- fp[!as.character(fp$flowpath_id) %in% old_ids, , drop = FALSE]

    # stub divides to satisfy prepare_network signature
    stub_sf <- sf::st_sf(
      divide_id = as.character(fp$flowpath_id),
      geometry  = sf::st_sfc(
        lapply(seq_len(nrow(fp)), function(i) sf::st_geometrycollection()),
        crs = sf::st_crs(fp)))

    network_list <- prepare_network(
      list(flowpaths = fp, divides = stub_sf))

    mapping_table <- .build_hw_collapse_fp(network_list, min_area_sqkm, min_length_km, rec_touch, rec_area, ideal_size_sqkm, connector_ids, max_area_sqkm)
    n_stuck <- max(n_stuck, .hf_collapse_stats$deferred_levelpath %||% 0L)
  }

  network_list <- add_network_type(network_list)
  message(sprintf("  Collapsed %d features.", start - nrow(network_list$flowpaths)))
  # Anything still deferred when the loop has converged is deferred permanently:
  # the blocking same-levelpath upstream is a POI outlet, which is never a donor,
  # so the fold can never become safe. These units stay below the floor because
  # honouring the POI's position outranks the size preference.
  if (isTRUE(n_stuck > 0L))
    message(sprintf(paste0("  %d collapse(s) permanently declined: folding them ",
      "would strand a same-levelpath upstream outlet (usually a POI)"), n_stuck))
  network_list
}


# -- .build_absorb_small_fp ----------------------------------------------------
# Mirror of headwater collapse for FIXED-OUTLET units. A small unit whose outlet
# is a POI or a terminal cannot itself collapse downstream (its outlet must be
# preserved), so instead it ABSORBS its single, non-POI, contiguous UPSTREAM
# continuation: growing the unit while keeping the same outlet. Skips
# confluences (>1 upstream -> would cross a junction) and any upstream that
# carries its own POI (that POI must remain an outlet). Area is summed from the
# member rec_div polygons (`rec_area`) so the "still small?" test stays correct
# as members accumulate across rounds (fp$areasqkm is not re-derived mid-loop).
.build_absorb_small_fp <- function(network_list, min_area_sqkm, rec_touch, rec_area, connector_ids = NULL,
                                   max_area_sqkm = NULL) {
  fp      <- sf::st_drop_geometry(network_list$flowpaths)
  all_ids <- as.character(fp$flowpath_id)
  toids   <- as.character(fp$flowpath_toid)
  poi     <- stats::setNames(!is.na(fp$poi_id), all_ids)
  lp      <- stats::setNames(as.character(fp$levelpathid), all_ids)  # mainstem id

  fl_split <- strsplit(as.character(fp$flowline_id), ",", fixed = TRUE)
  members  <- stats::setNames(lapply(fl_split, trimws), all_ids)
  area_of  <- vapply(members, function(m) sum(rec_area[m], na.rm = TRUE), numeric(1L))
  # connector unit = all members are reference-divide-less reaches
  is_conn  <- if (length(connector_ids))
    vapply(members, function(m) length(m) && all(m %in% connector_ids), logical(1L)) else
    stats::setNames(rep(FALSE, length(all_ids)), all_ids)

  # MAINSTEM-continuation upstream: of a unit's upstream contributors, the one
  # sharing its levelpathid (the river continuing up). At a confluence the other
  # upstream(s) are side tributaries on different levelpaths, they keep draining
  # in (junction preserved); only the mainstem reach is absorbed. A levelpath is a
  # single path so there is at most one same-lp upstream; the tie-break is a
  # guard that should not fire. It takes the SMALLEST, following the OG
  # .middle_massage, which merges a below-minimum group into the smaller of its
  # neighbours so a floor-driven merge overshoots the target as little as it can.
  # Lets a small POI/terminal outlet grow UP ITS OWN RIVER rather than remain a
  # sub-km2 sliver stranded at a junction.
  ups    <- split(all_ids, toids)
  ms_up  <- vapply(all_ids, function(id) {
    u <- ups[[id]]
    if (length(u) == 0L) return(NA_character_)
    same <- u[!is.na(lp[u]) & lp[u] == lp[id]]
    if (length(same) == 0L) NA_character_
    else if (length(same) == 1L) same
    else same[which.min(area_of[same])]
  }, character(1L))

  is_terminal <- is.na(toids) | toids == "0" | !(toids %in% all_ids)
  # absorb-upstream eligibility: units that cannot collapse DOWNSTREAM and so
  # would otherwise stay a sub-km2 sliver: (a) fixed outlets (POI/terminal,
  # outlet must be preserved) and (b) CONFLUENCES (>=2 inflow, excluded from the
  # downstream collapse to keep the junction). All grow up their mainstem
  # continuation while side tributaries keep draining in (junction preserved).
  inflow_n    <- tabulate(match(toids[toids %in% all_ids], all_ids), nbins = length(all_ids))
  is_conf     <- inflow_n >= 2L
  eligible    <- poi | is_terminal | is_conf

  into <- character(0)
  absorb <- character(0)
  n_over <- 0L
  for (i in seq_along(all_ids)) {
    if (!eligible[i] || area_of[i] >= min_area_sqkm) next
    S <- all_ids[i]
    U <- ms_up[[S]]
    # need a mainstem upstream that is not itself a required outlet (its POI must
    # stay an outlet -> two gages in series on the mainstem are never merged),
    # and the same connector/catchment type (never absorb a divide-less reach
    # into a catchment outlet, or vice versa)
    if (is.na(U) || isTRUE(poi[U]) || isTRUE(is_conn[[S]]) != isTRUE(is_conn[[U]])) next
    sm <- members[[S]]
    um <- members[[U]]
    if (is.null(sm) || is.null(um)) next
    if (!any(unlist(rec_touch[um], use.names = FALSE) %in% sm)) next   # not contiguous
    # Same ceiling the downstream collapse honours. This absorb is floor-driven
    # and used to be unbounded, so a sub-3 outlet could pull in a large mainstem
    # upstream and land far past the target: clearing one criterion by missing
    # the other, which is the trade the ceiling exists to refuse. A unit that
    # cannot be rescued within the ceiling stays small and is reported as such.
    # Report, never refuse. This is the only path left for a fixed outlet that
    # cannot collapse downstream, so blocking it strands the unit below the floor
    # permanently, the outcome the floor mandate exists to prevent.
    if (!is.null(max_area_sqkm) &&
        is.finite(area_of[i]) && is.finite(area_of[[U]]) &&
        (area_of[i] + area_of[[U]]) > max_area_sqkm) n_over <- n_over + 1L
    into <- c(into, S)
    absorb <- c(absorb, U)
  }
  # NO-CHAIN GUARD: a unit must not be both absorbed (removed) and an absorber
  # (winner) in the same round, that happens on a chain S<-U<-UU where U is
  # itself a small confluence absorbing UU; applying both at once removes U yet
  # remaps UU's toid to U -> dangling toid ("All toIDs not present"). Defer any
  # mapping whose absorbed unit is itself absorbing this round; the
  # .absorb_small_loop iteration resolves the chain over successive passes.
  keep   <- !(absorb %in% into)
  if (n_over) message(sprintf(paste0("  [floor-mandate] %d fixed outlet(s) grew ",
    "past the %g km2 ceiling absorbing upstream: no smaller rescue existed"),
    n_over, max_area_sqkm))
  data.frame(absorb = absorb[keep], into = into[keep], stringsAsFactors = FALSE)
}

# -- .absorb_small_loop --------------------------------------------------------
# Iterate .build_absorb_small_fp, merging each upstream unit INTO its small
# fixed-outlet downstream until none remain small/mergeable. Topology: the
# absorbed upstream's own inflows are redirected to the absorbing unit; the
# absorbing unit keeps its outlet and downstream link.
.absorb_small_loop <- function(network_list, min_area_sqkm, rec_touch, rec_area, connector_ids = NULL,
                               max_area_sqkm = NULL) {
  if (is.null(rec_touch) || is.null(rec_area)) return(network_list)
  message("--- Absorb Upstream into Small Fixed-Outlet Units (fp-only) ---")
  start <- nrow(network_list$flowpaths)
  mt    <- .build_absorb_small_fp(network_list, min_area_sqkm, rec_touch, rec_area, connector_ids, max_area_sqkm)
  count <- 0L
  while (nrow(mt) > 0L) {
    count <- count + 1L
    message(sprintf("  Absorbing: %d upstream unit(s) (round %d)", nrow(mt), count))
    fp  <- network_list$flowpaths
    ids <- as.character(fp$flowpath_id)
    win <- match(mt$into, ids)
    los <- match(mt$absorb, ids)
    ok  <- !is.na(win) & !is.na(los)
    win <- win[ok]
    los <- los[ok]
    mt <- mt[ok, , drop = FALSE]
    if (nrow(mt) == 0L) break

    # winner (into) keeps its outlet; absorb the upstream unit's members
    fp$member_comid[win] <- paste(fp$member_comid[win], fp$member_comid[los], sep = ",")
    fp$flowline_id[win]  <- paste(fp$flowline_id[win],  fp$flowline_id[los],  sep = ",")

    # redirect the absorbed unit's inflows to the absorbing unit
    old_ids <- as.character(mt$absorb)
    new_ids <- as.character(mt$into)
    tt <- as.character(fp$flowpath_toid)
    ri <- match(tt, old_ids)
    dr <- !is.na(ri)
    if (any(dr)) fp$flowpath_toid[dr] <- new_ids[ri[dr]]

    fp <- fp[!ids %in% old_ids, , drop = FALSE]

    stub_sf <- sf::st_sf(
      divide_id = as.character(fp$flowpath_id),
      geometry  = sf::st_sfc(
        lapply(seq_len(nrow(fp)), function(i) sf::st_geometrycollection()),
        crs = sf::st_crs(fp)))
    network_list <- prepare_network(list(flowpaths = fp, divides = stub_sf))
    mt <- .build_absorb_small_fp(network_list, min_area_sqkm, rec_touch, rec_area, connector_ids, max_area_sqkm)
  }
  network_list <- add_network_type(network_list)
  message(sprintf("  Absorbed %d upstream unit(s).", start - nrow(network_list$flowpaths)))
  network_list
}


#' Aggregate flowpath network topology without divide geometry
#'
#' @description
#' Flowpath-only replacement for the aggregate-to-distribution step.
#' Runs mainstem aggregation and headwater collapse using purely graph/tabular
#' operations, no `union_polygons`, `terra::makeValid`, `st_snap`, or
#' `terra::relate`. Returns aggregated flowpaths, network table, and poi layer.
#'
#' The returned `divides` element is a geometry-free stub (`data.frame`) carrying
#' only `divide_id` and `flowpath_toid`. Pass it with `rec_div` to
#' [construct_agg_divides()] to build real divide geometry.
#'
#' @param flowpath Reconciled flowpath `sf` from Stage-1 `rec_fp_agg`.
#' @param divide Reconciled divide `sf` from Stage-1 `rec_div`. Used only by
#'   [prepare_network()] for initial topology setup.
#' @param pois POI `sf` from the refactor phase.
#' @param vpu VPU identifier written to network table.
#' @param max_junction,ideal_size_sqkm,min_area_sqkm,max_area_sqkm,min_length_km Aggregation
#'   parameters.
#' @param connector_ids Optional character vector of catchment-less connector
#'   reach ids (e.g. divide-less reaches over waterbodies). Such reaches are
#'   held out of the small-catchment collapse/absorb steps so they are
#'   preserved rather than merged away. `NULL` (default) treats none as
#'   connectors.
#'
#' @return Named list: `flowpaths` (sf), `divides` (stub data.frame),
#'   `network` (data.frame), `pois` (sf or NULL).
#' @examples
#' \dontrun{
#' s1_fp   <- sf::st_read("rec_fp_agg.gpkg", "flowpaths")
#' rec_div <- sf::st_read("rec_div.gpkg", "divides")
#' pois    <- sf::st_read("pois.gpkg", "pois")
#' agg <- aggregate_network(flowpath = s1_fp, divide = rec_div,
#'                          pois = pois, vpu = "01")
#' }
#' @export
aggregate_network <- function(
    flowpath,
    divide,
    pois            = NULL,
    vpu             = NULL,
    max_junction    = 2,
    ideal_size_sqkm = 10,
    min_area_sqkm   = 3,
    max_area_sqkm   = 15,
    min_length_km   = 0.25,
    connector_ids   = NULL
    ) {
  stopifnot(
    inherits(flowpath, "sf"),
    inherits(divide,   "sf"),
    all(c("flowpath_id", "flowpath_toid", "levelpathid",
      "hydroseq", "lengthkm", "member_comid") %in% names(flowpath))
  )

  # -- prepare network ---------------------------------------------------
  nl <- add_network_type(
    prepare_network(list(flowpaths = flowpath, divides = divide)))

  # -- poi join ----------------------------------------------------------
  has_pois <- !is.null(pois) && nrow(pois) > 0L
  if (has_pois) {
    names(pois) <- tolower(names(pois))
    outflows <- dplyr::distinct(
      dplyr::select(sf::st_drop_geometry(pois), poi_id, flowpath_id),
      flowpath_id, .keep_all = TRUE)
    nl$flowpaths <- dplyr::left_join(
      dplyr::mutate(nl$flowpaths, poi_id = NULL),
      sf::st_drop_geometry(outflows), by = "flowpath_id")
  } else {
    nl$flowpaths$poi_id <- NA_character_
    pois <- NULL
  }

  # POI identity trace through aggregation. A minted poi_id is a commitment, and
  # aggregation is where one can quietly stop being carried: the id lives on a
  # reconciled reach, and every merge has to move it onto the surviving aggregate.
  # Reporting the distinct-id count after each step localises a loss to the step
  # that caused it instead of leaving it visible only in the final product.
  .poi_trace <- function(x, tag) {
    if (is.null(x) || !"poi_id" %in% names(x)) return(invisible(NULL))
    v <- as.character(x$poi_id)
    v <- v[!is.na(v) & nzchar(v)]
    ids <- unique(unlist(strsplit(v, ",", fixed = TRUE)))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    message(sprintf("  [poi] %-24s rows carrying an id=%5d  distinct ids=%5d",
      tag, length(v), length(ids)))
    invisible(length(ids))
  }
  .poi_trace(nl$flowpaths, "after poi join")

  # rec_div adjacency (computed ONCE) for the headwater-collapse contiguity
  # check, the geometry-aware constraint (divide_touches) that
  # keeps every aggregate's catchment contiguous, so the divide union never goes
  # disjoint and union_polygons never silently drops area.
  rid_col <- if ("flowpath_id" %in% names(divide)) "flowpath_id" else "divide_id"
  rec_touch <- tryCatch(
    {
      rid <- as.character(divide[[rid_col]])
      ix  <- sf::st_intersects(divide)
      stats::setNames(lapply(ix, function(j) rid[j]), rid)
    },
    error = function(e) {
      .hf_log_warn("rec_div adjacency unavailable; headwater collapse runs without ",
        "contiguity check: ", conditionMessage(e))
      NULL
    })
  # per-rec_div area (member-summed inside the absorption pass, so it stays
  # correct as members accumulate: independent of any stale fp$areasqkm)
  rec_area <- tryCatch(
    stats::setNames(as.numeric(sf::st_area(divide)) / 1e6,
      as.character(divide[[rid_col]])),
    error = function(e) NULL)

  # -- mainstem aggregation (fp-only aggregate_sets) ---------------------
  message("--- Aggregate Along Mainstem ---")
  upstream_index <- dplyr::ungroup(dplyr::slice(
    dplyr::group_by(
      dplyr::filter(
        dplyr::distinct(dplyr::select(sf::st_drop_geometry(nl$flowpaths),
          flowpath_id, hl_un = poi_id)),
        !is.na(hl_un)),
      flowpath_id), 1L))

  fline <- dplyr::distinct(dplyr::left_join(
    dplyr::mutate(sf::st_drop_geometry(nl$flowpaths), hl_dn = poi_id),
    upstream_index, by = "flowpath_id"))
  # per-reach connector flag (reference areasqkm=0 waterbody reaches) so cs_group
  # breaks the levelpath at the catchment<->connector boundary and never groups
  # the two. connector_ids are in the current flowpath_id space.
  conn <- as.character(connector_ids)
  fline$has_divide <- !(as.character(fline$flowpath_id) %in% conn)
  # Reconciled ids absorbed by a pre-aggregation collapse (Step 3c.5). They are
  # members of this flowpath even though they have no row of their own, so they must
  # reach rec_id_list; construct_agg_divides looks up member divides by these ids.
  if (!"absorbed_recs" %in% names(fline))
    fline$absorbed_recs <- NA_character_

  index_table <- dplyr::select(
    dplyr::ungroup(dplyr::mutate(
      dplyr::group_by(
        dplyr::ungroup(dplyr::mutate(
          dplyr::group_by(fline, levelpathid),
          ind = .cs_group(
            areasqkm, lengthkm, max_junction, hl_dn, hl_un,
            ideal_size_sqkm, min_area_sqkm, min_length_km,
            has_divide = has_divide))),
        levelpathid, ind),
      set = dplyr::cur_group_id(), n = dplyr::n())),
    set, flowpath_id, flowpath_toid, levelpathid, hydroseq, member_comid, poi_id,
    absorbed_recs, n)

  main_agg <- .agg_sets_fp(nl, index_table, rec_touch = rec_touch)
  main_agg <- add_network_type(main_agg)
  .poi_trace(main_agg$flowpaths, "after .agg_sets_fp")
  message(sprintf("  Merged: %d -> %d flowpaths",
    nrow(nl$flowpaths), nrow(main_agg$flowpaths)))

  # -- collapse + absorb, ITERATED TO CONVERGENCE ------------------------
  # The collapse folds small non-fixed-outlet units downstream toward
  # ideal_size; the absorb grows small fixed-outlet (POI/terminal) units up
  # their mainstem. They INTERACT: an absorb redirects an absorbed unit's
  # inflows and can leave a newly-collapsible mid-reach/headwater, and a collapse
  # can free a fixed-outlet sliver to absorb. Running each once leaves those
  # stragglers, so iterate until a full absorb+collapse pass removes no flowpath.
  # Both ops only ever MERGE (flowpath count strictly decreases), so this is
  # bounded; the seq_len cap is a backstop.
  nl <- .collapse_hw_fp(main_agg, min_area_sqkm, min_length_km, rec_touch, rec_area, ideal_size_sqkm, conn, max_area_sqkm)
  .poi_trace(nl$flowpaths, "after collapse (pre-loop)")
  for (.iter in seq_len(10L)) {
    before <- nrow(nl$flowpaths)
    nl <- .absorb_small_loop(nl, min_area_sqkm, rec_touch, rec_area, conn, max_area_sqkm)
    .poi_trace(nl$flowpaths, sprintf("iter %d after absorb", .iter))
    nl <- .collapse_hw_fp(nl, min_area_sqkm, min_length_km, rec_touch, rec_area, ideal_size_sqkm, conn, max_area_sqkm)
    .poi_trace(nl$flowpaths, sprintf("iter %d after collapse", .iter))
    if (nrow(nl$flowpaths) >= before) break
  }
  rm(main_agg)

  # -- finalize flowpaths ------------------------------------------------
  div_ids <- as.character(nl$divides$divide_id)
  nl$flowpaths <- dplyr::select(
    nl$flowpaths,
    flowpath_id, flowpath_toid,
    flowline_id,
    mainstem = levelpathid,
    member_comid, poi_id, hydroseq, lengthkm, areasqkm,
    tot_drainage_areasqkm = tot_drainage_area,
    has_divide) |>
    dplyr::mutate(
      divide_id = ifelse(flowpath_id %in% div_ids, flowpath_id, NA_character_))

  # clear self-loop toids
  sl <- which(!is.na(nl$flowpaths$flowpath_toid) &
    as.character(nl$flowpaths$flowpath_id) ==
      as.character(nl$flowpaths$flowpath_toid))
  if (length(sl) > 0L) {
    .hf_log_warn(sprintf("aggregate_network: %d self-loop toid(s) cleared", length(sl)))
    nl$flowpaths$flowpath_toid[sl] <- "0"
  }

  # -- remap toids that reference an ABSORBED rec id ---------------------
  # .agg_sets_fp / .collapse_hw_fp can leave a flowpath_toid pointing at a
  # reconciled id that was merged INTO another aggregate (it lives in that
  # aggregate's flowline_id member list, no longer a standalone flowpath_id).
  # Such a toid dangles and trips check_network_validity ("All toIDs are not
  # present in network"). Redirect it to the absorbing aggregate's flowpath_id
  # the remap the aggregate-to-distribution step performed internally.
  # Unresolvable + non-terminal toids drain out of the modeled domain -> terminal
  # "0".
  {
    fp_ids_agg <- as.character(nl$flowpaths$flowpath_id)
    toids_agg  <- as.character(nl$flowpaths$flowpath_toid)
    fl_lists   <- strsplit(as.character(nl$flowpaths$flowline_id), ",", fixed = TRUE)
    rec2agg    <- stats::setNames(rep(fp_ids_agg, lengths(fl_lists)),
      trimws(unlist(fl_lists, use.names = FALSE)))
    tnum       <- suppressWarnings(as.numeric(toids_agg))
    need       <- which(!is.na(toids_agg) & nzchar(toids_agg) & toids_agg != "0" &
      !(toids_agg %in% fp_ids_agg) &
      !(is.finite(tnum) & tnum >= 1e9))
    if (length(need) > 0L) {
      mapped <- unname(rec2agg[toids_agg[need]])
      mapped[is.na(mapped) | mapped == fp_ids_agg[need]] <- "0"  # gone or self-loop -> terminal
      nl$flowpaths$flowpath_toid[need] <- mapped
      message(sprintf("  aggregate_network: remapped %d dangling toid(s) to absorbing aggregate",
        length(need)))
    }
  }

  # -- divide stub -------------------------------------------------------
  topo <- dplyr::select(sf::st_drop_geometry(nl$flowpaths), divide_id, flowpath_toid)
  matched_div_ids <- unique(stats::na.omit(nl$flowpaths$divide_id))
  div_stub <- dplyr::mutate(
    dplyr::left_join(
      data.frame(divide_id = matched_div_ids, stringsAsFactors = FALSE),
      topo, by = "divide_id"),
    has_flowline = TRUE,
    ds_id        = NA_character_,
    type         = "network")

  # -- network table -----------------------------------------------------
  network <- dplyr::distinct(
    dplyr::left_join(
      dplyr::mutate(
        .separate_longer(
          dplyr::select(
            sf::st_drop_geometry(nl$flowpaths),
            flowpath_id, flowpath_toid,
            member = member_comid,
            flowline_id, divide_id,
            dplyr::any_of("poi_id"),
            mainstem, hydroseq, lengthkm, areasqkm, tot_drainage_areasqkm),
          col = "member", delim = ","),
        hf_id_part = sapply(strsplit(member, ".", fixed = TRUE), `[`, 2L),
        hf_id_part = ifelse(is.na(hf_id_part), 1L,
          floor(as.numeric(hf_id_part))),
        hf_id      = as.numeric(sapply(strsplit(member, ".", fixed = TRUE), `[`, 1L)),
        member     = NULL,
        hf_source  = "NHDPlusV2"),
      dplyr::select(div_stub, divide_id, type, ds_id),
      by = "divide_id"))

  if (!is.null(vpu)) network$vpu <- vpu else network$vpu <- NA_character_

  # rename flowline_id -> rec_id_list to match stage2 convention
  names(network)[names(network) == "flowline_id"] <- "rec_id_list"
  if ("flowline_id" %in% names(nl$flowpaths))
    names(nl$flowpaths)[names(nl$flowpaths) == "flowline_id"] <- "rec_id_list"


  # -- poi layer ---------------------------------------------------------
  .poi_trace(nl$flowpaths, "before poi layer build")
  pois_out <- NULL
  if (has_pois) {
    agg_poi <- unique(sf::st_drop_geometry(
      nl$flowpaths[!is.na(nl$flowpaths$poi_id),
        c("flowpath_id", "poi_id"), drop = FALSE]))
    if (nrow(agg_poi) > 0L) {
      pois_df <- sf::st_drop_geometry(pois)
      names(pois_df) <- tolower(names(pois_df))
      merged  <- merge(agg_poi, pois_df, by = "poi_id",
        all.x = TRUE, suffixes = c("", "_orig"))
      merged  <- merged[!duplicated(merged$flowpath_id), , drop = FALSE]
      gm_idx  <- match(as.character(merged$poi_id), as.character(pois$poi_id))
      if (inherits(pois, "sf") && any(!is.na(gm_idx))) {
        pois_out <- sf::st_sf(merged,
          geometry = sf::st_geometry(pois)[gm_idx],
          crs      = sf::st_crs(pois))
      } else {
        pois_out <- merged
      }
    }
  }

  list(
    flowpaths = nl$flowpaths,
    divides   = div_stub,
    network   = network,
    pois      = pois_out
  )
}


#' Construct aggregated divides from Stage-1 rec_div via rec_id_list
#'
#' @description
#' Given the network table from [aggregate_network()] (which carries
#' `flowpath_id` -> `rec_id_list`) and the Stage-1 `rec_div` polygons, unions
#' the `rec_div` members for each aggregated flowpath.
#'
#' No raster, no snapping, no mapshaper. Stage-3 `fill_small_holes` and
#' `clean_gpkg_layer` handle any residual topology.
#'
#' @param network Network `data.frame` from [aggregate_network()] with
#'   `flowpath_id` and `rec_id_list`.
#' @param rec_div Stage-1 reconciled divide `sf` with `flowpath_id` (= rec_id).
#' @param flowpaths Aggregated flowpath `sf` from [aggregate_network()].
#'
#' @return `sf` with `divide_id`, `flowpath_toid`, `areasqkm`,
#'   `has_flowline`, `ds_id`, `type`, geometry (POLYGON).
#' @examples
#' \dontrun{
#' agg     <- aggregate_network(flowpath = s1_fp, divide = rec_div, pois = pois)
#' rec_div <- sf::st_read("rec_div.gpkg", "divides")
#' divides <- construct_agg_divides(network   = agg$network,
#'                                  rec_div   = rec_div,
#'                                  flowpaths = agg$flowpaths)
#' }
#' @export
construct_agg_divides <- function(network, rec_div, flowpaths) {
  stopifnot(
    is.data.frame(network),
    all(c("flowpath_id", "rec_id_list") %in% names(network)),
    inherits(rec_div, "sf"),
    "flowpath_id" %in% names(rec_div)
  )

  rd         <- rec_div
  rd$rec_id  <- as.character(rd$flowpath_id)
  target_crs <- sf::st_crs(rd)

  # Build fp -> rec_id mapping from the network rec_id_list: every member rec_id
  # of each aggregate fp. The fp's divide is the union of those members'
  # rec_div polygons (assembled below); fp attributes are joined on at the end.
  net_unique <- unique(data.frame(
    flowpath_id = as.character(network$flowpath_id),
    rec_id_list = as.character(network$rec_id_list),
    stringsAsFactors = FALSE))

  fp_to_rec <- data.table::rbindlist(lapply(seq_len(nrow(net_unique)), function(i) {
    recs <- trimws(strsplit(net_unique$rec_id_list[i], ",", fixed = TRUE)[[1]])
    recs <- recs[nzchar(recs) & recs %in% rd$rec_id]
    if (!length(recs)) return(NULL)
    data.table::data.table(flowpath_id = net_unique$flowpath_id[i], rec_id = recs)
  }))

  if (nrow(fp_to_rec) == 0L)
    stop("construct_agg_divides: no rec_id matches found in rec_div")

  fp_ids <- unique(fp_to_rec$flowpath_id)
  message(sprintf("construct_agg_divides: unioning rec_div for %d aggregated fp(s)",
    length(fp_ids)))

  # join rec_div geometry onto the expanded mapping then grouped union,
  # single vectorised summarise instead of many serial st_union calls
  rd_sub <- rd[rd$rec_id %in% fp_to_rec$rec_id, c("rec_id"), drop = FALSE]
  rd_sub$rec_id <- as.character(rd_sub$rec_id)
  fp_to_rec_df  <- as.data.frame(fp_to_rec)

  tagged <- dplyr::left_join(rd_sub, fp_to_rec_df, by = "rec_id")
  tagged <- tagged[!is.na(tagged$flowpath_id), , drop = FALSE]
  # Backstop: never feed an empty geometry to union_polygons. terra::aggregate
  # inside union_polygons drops empty geometries but keeps their attribute row,
  # producing an sf whose nrow != geometry length ("nrow dataframe does not
  # match nrow geometry"). Empty members contribute no area, so dropping them is
  # area-safe and defends construct against any upstream step that leaks one.
  if (nrow(tagged) > 0L) {
    .te <- sf::st_is_empty(tagged)
    if (any(.te)) tagged <- tagged[!.te, , drop = FALSE]
  }

  divides <- tagged |>
    dplyr::select(divide_id = flowpath_id) |>
    hfutils::union_polygons("divide_id") |>
    sf::st_as_sf()

  if (nrow(divides) == 0L)
    stop("construct_agg_divides: no divides produced")

  # Keep one row per aggregate. Casting to POLYGON explodes a multi-part divide
  # into one row per part, all sharing a divide_id, which downstream steps then
  # have to merge back.
  #
  # A divide is legitimately multi-part when its members' catchments are not
  # contiguous, and MULTIPOLYGON divides ship on every row, so nothing downstream
  # needs single-part input.
  divides <- suppressWarnings(sf::st_cast(divides, "MULTIPOLYGON"))

  # join flowpath_toid from aggregated flowpaths
  fp_attr <- sf::st_drop_geometry(flowpaths)[,
    intersect(c("flowpath_id", "flowpath_toid", "has_divide"),
      names(flowpaths)), drop = FALSE]
  fp_attr$flowpath_id <- as.character(fp_attr$flowpath_id)
  divides$divide_id   <- as.character(divides$divide_id)

  divides <- dplyr::left_join(divides, fp_attr,
    by = c("divide_id" = "flowpath_id"))
  divides$areasqkm    <- as.numeric(sf::st_area(divides)) / 1e6
  divides$has_flowline <- TRUE
  divides$ds_id        <- NA_character_
  divides$type         <- "network"

  message(sprintf("construct_agg_divides: %d divide(s) from %d rec_div polygon(s)",
    nrow(divides), nrow(rd)))
  divides
}

#' Prepare an externally-supplied hydrolocation set for embedding
#'
#' Maps a raw hydrolocation point dataset (e.g. the gage / assimilation points
#' supplied by OWP) to the canonical `hydrolocations` schema that
#' [run_pipeline()] embeds, and optionally writes it as a GeoPackage suitable for
#' the `hydrolocations_gpkg` argument.
#'
#' The POIs are **snapped onto `flowpaths`** and assigned `flowpath_id`,
#' `poi_id` (one per flowpath, co-located POIs bundled), and linear-reference
#' `pathmeas`/`reachmeas`. This pre-assignment is required: the refactor phase's spatial
#' re-snap fallback is lossy (it drops many POIs when the input lacks
#' `flowpath_id`/`poi_id`), so raw points MUST be snapped here first. Rows
#' without `hl_reference`, or farther than `max_dist` from any flowpath, are
#' dropped with a warning.
#'
#' @param x An `sf` of hydrolocation POINTs, or a path to a point GeoPackage.
#' @param flowpaths Reference flowpaths to snap onto: an `sf` or a path to the
#'   reference GeoPackage (its `flowpaths` layer is read). Required.
#' @param mapping Named character vector mapping *canonical* column ->
#'   *source* column, e.g. `c(hl_reference = "gage_id", hl_type = "type")`.
#'   Canonical columns: `hl_reference` (REQUIRED), `hl_type`, `hl_link`,
#'   `hl_class`, `hl_source`. Anything unmapped/absent is filled `NA`.
#' @param max_dist Max snap distance (m); POIs farther than this are dropped.
#' @param link_keep_da_ratio A distrusted link is restored when the nearest
#'   flowpath drains less than this fraction of the reach the source names
#'   (default 100, i.e. under 1%). Guards against relocating a large-river
#'   location onto an adjacent minor feature.
#' @param link_keep_min_da_sqkm Only apply that rescue when the named reach
#'   drains at least this much, so the comparison is meaningful. Default 100.
#' @param max_link_m Max distance a point may be from the reach its source
#'   `flowpath_id` prescribes before that link is distrusted and the nearest
#'   flowpath used instead. Reported when it triggers.
#'   Default 1000.
#' @param out_gpkg Optional output path. When given, writes the `hydrolocations`
#'   layer there and returns the path; otherwise returns the prepared `sf`.
#' @param layer Layer to read when `x` is a GeoPackage path (default: first).
#' @return The prepared hydrolocations `sf`, or `out_gpkg` when written.
#' @examples
#' \dontrun{
#' # OWP delivers gage points with columns `siteid` and `network`:
#' prepare_hydrolocations(
#'   "owp_gages.gpkg", flowpaths = "reference.gpkg",
#'   mapping  = c(hl_reference = "siteid", hl_source = "network", hl_class = "gage"),
#'   out_gpkg = "hydrolocations.gpkg")
#' # -> run_pipeline(reference_gpkg = ref, hydrolocations_gpkg = "hydrolocations.gpkg", ...)
#' }
#' @export
prepare_hydrolocations <- function(x, flowpaths, mapping = character(0),
                                   max_dist = 1000, max_link_m = 250,
                                   link_keep_da_ratio = 100,
                                   link_keep_min_da_sqkm = 100,
                                   out_gpkg = NULL, layer = NULL) {
  hl <- if (is.character(x)) {
    lyr <- if (is.null(layer)) sf::st_layers(x)$name[1L] else layer
    sf::st_read(x, layer = lyr, quiet = TRUE)
  } else x
  if (!inherits(hl, "sf")) stop("`x` must be an sf object or a GeoPackage path", call. = FALSE)

  gt <- as.character(unique(sf::st_geometry_type(hl)))
  if (!all(gt %in% c("POINT", "MULTIPOINT")))
    stop("hydrolocations must be POINT geometry; got: ", paste(gt, collapse = ", "), call. = FALSE)

  # Apply the mapping: canonical <- source column.
  if (length(mapping)) {
    if (is.null(names(mapping)) || any(!nzchar(names(mapping))))
      stop("`mapping` must be a NAMED vector: c(<canonical> = <source>)", call. = FALSE)
    for (cn in names(mapping)) {
      src <- mapping[[cn]]
      if (!src %in% names(hl))
        stop(sprintf("mapping: source column '%s' (for '%s') not found in input", src, cn), call. = FALSE)
      hl[[cn]] <- hl[[src]]
    }
  }

  if (!"hl_reference" %in% names(hl))
    stop("`hl_reference` is required -- map it, e.g. mapping = c(hl_reference = '<col>')", call. = FALSE)

  # Drop rows with no hl_reference (Stage 3 promote_hydrolocations drops them).
  ref <- trimws(as.character(hl$hl_reference))
  bad <- is.na(ref) | !nzchar(ref)
  if (any(bad)) {
    warning(sprintf("prepare_hydrolocations: dropping %d/%d hydrolocation(s) with no hl_reference",
      sum(bad), nrow(hl)), call. = FALSE)
    hl <- hl[!bad, , drop = FALSE]
  }
  if (!nrow(hl)) stop("no hydrolocations remain after hl_reference validation", call. = FALSE)

  # -- snap onto flowpaths + assign flowpath_id / poi_id / linear-reference -----
  # The pipeline's own re-snap is lossy without these, so pre-assign them here.
  fp <- if (is.character(flowpaths)) sf::st_read(flowpaths, "flowpaths", quiet = TRUE) else flowpaths
  if (!inherits(fp, "sf")) stop("`flowpaths` must be an sf or a GeoPackage path", call. = FALSE)
  hl <- sf::st_transform(hl, sf::st_crs(fp))
  nf <- sf::st_nearest_feature(hl, fp)
  d  <- as.numeric(sf::st_distance(sf::st_geometry(hl), sf::st_geometry(fp)[nf], by_element = TRUE))
  far <- is.na(d) | d > max_dist
  if (any(far)) {
    warning(sprintf("prepare_hydrolocations: dropping %d/%d hydrolocation(s) > %g m from any flowpath",
      sum(far), nrow(hl), max_dist), call. = FALSE)
    hl <- hl[!far, , drop = FALSE]; nf <- nf[!far]
  }
  if (!nrow(hl)) stop("no hydrolocations remain after the ", max_dist, " m snap filter", call. = FALSE)
  # Prefer the reach link the source already carries.
  #
  # The link is an authoritative assignment; the nearest flowpath is a geometric
  # guess, and the two disagree for a fifth of hydrolocations: mostly coastal
  # coupling nodes, which sit near shorelines where the nearest line is often the
  # wrong reach. Nearest-feature is the fallback for rows with no usable link, and
  # the linear references are computed against whichever reach is chosen so
  # measures and assignment cannot disagree.
  fp_pos <- stats::setNames(seq_len(nrow(fp)), as.character(fp$flowpath_id))
  given  <- if ("flowpath_id" %in% names(hl)) as.character(hl$flowpath_id) else
    rep(NA_character_, nrow(hl))
  gi   <- unname(fp_pos[given])

  # Check, not a silent trust: how far is each point from the reach it is
  # prescribed to? The link is authoritative about WHICH reach, but a multi-km
  # separation means the link and the coordinate disagree in the source, and
  # relocating a gage that far is worse than declining to. Such rows are reported
  # and fall back to the nearest flowpath rather than being moved across the basin.
  #
  # This is measured on the RAW coordinate, and must stay that way. `flowpath_id`
  # is not unique across reference domains, so a link can resolve onto an unrelated
  # reach by id collision; the separation between the point and the reach its link
  # names is what distinguishes a real assignment from a collision. Measuring it
  # after any snap onto that same reach reports zero by construction and the guard
  # stops discriminating entirely.
  link_far <- rep(FALSE, nrow(hl))
  if (any(!is.na(gi))) {
    j  <- which(!is.na(gi))
    dl <- as.numeric(sf::st_distance(sf::st_geometry(hl)[j],
      sf::st_geometry(fp)[gi[j]], by_element = TRUE))
    link_far[j] <- dl > max_link_m
    if (any(link_far)) {
      rr <- if ("hl_reference" %in% names(hl)) as.character(hl$hl_reference) else
        rep(NA_character_, nrow(hl))
      ex <- utils::head(rr[link_far], 4)
      message(sprintf(paste0(
        "prepare_hydrolocations: CHECK -- %d row(s) are more than %.0fm from the ",
        "reach they are prescribed to (max %.0fm); using nearest flowpath for ",
        "these instead [%s]"),
        sum(link_far), max_link_m, max(dl[link_far[j]], na.rm = TRUE),
        paste(ex, collapse = ", ")))
    } else {
      message(sprintf(paste0(
        "prepare_hydrolocations: every source link is within %.0fm of the reach ",
        "it names (max %.0fm)"), max_link_m, max(dl, na.rm = TRUE)))
    }
  }

  # Distance alone is the wrong test when the nearest flowpath is not a plausible
  # host. A location named on a large river whose nearest feature drains a tiny
  # fraction of it is a coordinate that missed the channel, not a mis-stated link:
  # no streamgage on a 400,000 km2 river sits on a 1.4 km2 ditch. Falling back
  # there moves the location off the river it measures.
  #
  # So a distrusted link is restored when the named reach is a real river and the
  # nearest candidate drains less than 1/link_keep_da_ratio of it. The point still
  # does not move: only which reach it is attributed to.
  n_rescued <- 0L
  if (any(link_far) && "totdasqkm" %in% names(fp)) {
    .da <- suppressWarnings(as.numeric(fp$totdasqkm))
    lf  <- which(link_far)
    named_da <- .da[gi[lf]]
    near_da  <- .da[nf[lf]]
    keep <- is.finite(named_da) & is.finite(near_da) &
      named_da >= link_keep_min_da_sqkm &
      near_da  <  named_da / link_keep_da_ratio
    if (any(keep)) {
      link_far[lf[keep]] <- FALSE
      n_rescued <- sum(keep)
      message(sprintf(paste0(
        "prepare_hydrolocations: kept the source link for %d distrusted row(s) ",
        "whose nearest flowpath drains less than %.2g%% of the reach they name"),
        n_rescued, 100 / link_keep_da_ratio))
    }
  }

  # Resolve the reach ONCE, before anything is moved: the prescribed link where it
  # is trustworthy, the nearest flowpath where there is no link or the link is too
  # far to believe.
  idx <- ifelse(is.na(gi) | link_far, nf, gi)
  message(sprintf(paste0(
    "prepare_hydrolocations: reach link from source for %d of %d row(s); ",
    "nearest-flowpath fallback for %d (%d unlinked, %d distrusted); ",
    "source and nearest disagree on %d"),
    sum(!is.na(gi) & !link_far), nrow(hl), sum(is.na(gi) | link_far),
    sum(is.na(gi)), sum(link_far), sum(!is.na(gi) & gi != nf)))

  # A POI on a reach with no local catchment AND no upstream drainage can only
  # become a nexus with nothing flowing into it: no divide contributes to it and
  # nothing routes through it. Such a nexus is inert in NextGen and shows up
  # downstream as an isolated location with no contributing divides.
  #
  # Both facts are known here, so the decision is made here rather than building
  # the fabric around the POI and cleaning up afterwards.
  #
  # The test requires BOTH conditions. `has_divide == 0` alone is not evidence a
  # location is meaningless, a reach can carry no local catchment and still
  # route a large basin, so dropping on that alone would remove major gauged
  # rivers. Adding `totdasqkm == 0` removes only what can never carry water.
  if (all(c("has_divide", "totdasqkm") %in% names(fp))) {
    .hdv  <- suppressWarnings(as.integer(fp$has_divide))[idx]
    .tda  <- suppressWarnings(as.numeric(fp$totdasqkm))[idx]
    .dead <- !is.na(.hdv) & .hdv == 0L & (is.na(.tda) | .tda <= 0)
    if (any(.dead)) {
      .ty <- if ("hl_type" %in% names(hl))
        paste(sprintf("%s=%d", names(table(as.character(hl$hl_type)[.dead])),
                      as.integer(table(as.character(hl$hl_type)[.dead]))),
              collapse = ", ") else "unclassified"
      message(sprintf(paste0(
        "prepare_hydrolocations: dropping %d of %d hydrolocation(s) on reaches ",
        "with no catchment and no upstream drainage (they could only become ",
        "nexuses with no inflow): %s"), sum(.dead), nrow(hl), .ty))
      hl     <- hl[!.dead, , drop = FALSE]
      idx    <- idx[!.dead]
      nf     <- nf[!.dead]
      given  <- given[!.dead]
      if (!nrow(hl))
        stop("no hydrolocations remain after the empty-nexus filter", call. = FALSE)
    }
  }

  # Snap the hydrolocation onto the reach it has been assigned, at that reach's
  # NEAREST VERTEX.
  #
  # Honouring the assignment while leaving the point on its raw coordinate is
  # incoherent: the assignment names reach A, the coordinate may lie nearer reach B,
  # and every later measure is then taken against a reach the point is not on.
  #
  # Snapping to a vertex rather than to the nearest point on the line matters: a
  # vertex is a real node of the linestring, so the position survives the splitting
  # and reconciliation that follow, whereas an interpolated point on a segment does
  # not correspond to anything in the geometry.
  snap_i <- vapply(seq_len(nrow(hl)), function(k) {
    m <- sf::st_coordinates(sf::st_geometry(fp)[[idx[k]]])
    q <- sf::st_coordinates(sf::st_geometry(hl)[k])
    which.min((m[, 1] - q[1])^2 + (m[, 2] - q[2])^2)
  }, integer(1))
  snapped <- sf::st_sfc(lapply(seq_len(nrow(hl)), function(k) {
    m <- sf::st_coordinates(sf::st_geometry(fp)[[idx[k]]])
    sf::st_point(unname(m[snap_i[k], 1:2]))
  }), crs = sf::st_crs(fp))
  moved_m <- as.numeric(sf::st_distance(sf::st_geometry(hl), snapped,
    by_element = TRUE))
  sf::st_geometry(hl) <- snapped
  message(sprintf(paste0(
    "prepare_hydrolocations: snapped to nearest vertex of the assigned reach ",
    "(median %.1fm, p90 %.0fm, max %.0fm)"),
    stats::median(moved_m, na.rm = TRUE),
    stats::quantile(moved_m, 0.9, na.rm = TRUE), max(moved_m, na.rm = TRUE)))

  # Preserve the reach the SOURCE prescribed before overwriting it.
  #
  # `flowpath_id` is resolved here and then legitimately re-pointed downstream: a
  # hydrolocation within the split floor of a reach's inlet is reassigned to the
  # upstream reach whose outlet it is, because a POI belongs to the reach it
  # TERMINATES. Overwriting in place would leave a consumer joining on the
  # supplied id with an unexplained mismatch, so the supplied value is carried
  # through as its own column.
  hl$hl_flowpath_id_source <- given
  hl$flowpath_id <- as.character(fp$flowpath_id[idx])
  if ("vpuid" %in% names(fp)) hl$vpuid <- as.character(fp$vpuid[idx])
  pm <- 1 - sf::st_line_project(sf::st_geometry(fp)[idx], sf::st_geometry(hl), normalized = TRUE)
  hl$pathmeas <- 100 * pm
  hl$reachmeas <- if (all(c("frommeas", "tomeas") %in% names(fp))) {
    fm <- fp$frommeas[idx]; tm <- fp$tomeas[idx]; (1 - pm) * (tm - fm) + fm
  } else NA_real_
  # `idx` and `snap_i` are parallel to the rows, so they have to travel with them.
  # They did not: the sort moved the rows and left both behind, and the
  # along-reach position computed below then read another row's reach and another
  # row's vertex. That value decides which hydrolocations collapse into one POI,
  # so it mattered for every row the sort displaced.
  .o     <- order(hl$flowpath_id)
  hl     <- hl[.o, ]
  idx    <- idx[.o]
  snap_i <- snap_i[.o]
  hl$hl_id  <- seq_len(nrow(hl))

  # A POI is a LOCATION, so collocation means "the same place", not "the same
  # reach". Assigning one poi_id per flowpath merged hydrolocations that merely
  # share a reference reach: two gages a kilometre apart on one reach became one
  # POI. Grouped by position instead, at 1 m, which is coordinate precision rather
  # than a search radius.
  # Collocation is along-line proximity on the shared reach, at the same 250 m used
  # by the placement rule and by split_flowlines' min_split_m, so a POI this rule
  # keeps can never be one the splitter then discards for leaving a short piece.
  COLLAPSE_M <- 250
  along <- vapply(seq_len(nrow(hl)), function(k) {
    m <- sf::st_coordinates(sf::st_geometry(fp)[[idx[k]]])
    cd <- cumsum(c(0, sqrt(diff(m[, 1])^2 + diff(m[, 2])^2)))
    cd[snap_i[k]]
  }, numeric(1))
  grp <- paste(as.character(hl$flowpath_id),
    floor(along / COLLAPSE_M), sep = "\r")
  hl$poi_id <- match(grp, unique(grp))
  message(sprintf(
    "prepare_hydrolocations: %d row(s) -> %d POI location(s) (%d collocated)",
    nrow(hl), length(unique(hl$poi_id)), nrow(hl) - length(unique(hl$poi_id))))

  # Fill any missing canonical columns and select the canonical schema.
  canon <- c("vpuid", "flowpath_id", "hl_flowpath_id_source", "poi_id", "hl_id",
             "hl_link", "hl_type",
             "hl_reference", "hl_class", "hl_source", "hl_uri", "pathmeas", "reachmeas")
  for (cn in setdiff(canon, names(hl))) hl[[cn]] <- NA
  gcol <- attr(hl, "sf_column")
  hl <- hl[, c(canon, gcol)]

  if (is.null(out_gpkg)) return(hl)
  if (file.exists(out_gpkg)) file.remove(out_gpkg)
  sf::st_write(hl, out_gpkg, "hydrolocations", quiet = TRUE)
  message(sprintf("prepare_hydrolocations: wrote %d hydrolocation(s) -> %s",
    nrow(hl), out_gpkg))
  invisible(out_gpkg)
}

#' Pipeline invariant checks
#'
#' @param stage Character. One of `"refactored"`, `"reconciled"`,
#'   `"aggregated"`, `"ngen"`.
#' @param ... Stage-specific named arguments. For `"ngen"`: `flowpaths`,
#'   `divides`, `nexus`, and optionally `network` (the aspatial network table,
#'   used to check flowline_id global uniqueness, flowline_toid referential
#'   integrity, and flowline graph acyclicity) and `lakes` (the waterbody
#'   layer, used to check that lake_id stamps intersect their lake).
#' @param strict Logical. If `TRUE` (default) any failed check throws an error.
#' @param coverage Logical. Run expensive flowpath-in-catchment coverage check.
#' @param coverage_min Minimum fraction of a flowpath inside its catchment.
#' @param area Logical. Run the geometric non-overlap check (computes the union
#'   of all network divides: moderately expensive on large VPUs). The cheap
#'   total-area, area-conservation, and shared-divide checks always run for
#'   stages that carry a `divides` layer, regardless of this flag.
#' @param area_baseline_m2 Numeric. Total area (m2) of the upstream divide set
#'   this stage must conserve (typically the sum of the input divides' areas).
#'   When supplied, `divides_area_conserved` asserts the output total stays
#'   within `area_tol_frac` of it. `NULL` (default) reports the total as info
#'   only.
#' @param coverage_baseline_m2 Numeric. Union area (m2) of the upstream divide
#'   set, i.e. the ground it covers rather than the sum of its parts. When
#'   supplied, `divides_coverage_conserved` reports the change in footprint.
#'   This is what detects lost ground: removing an overlap and punching a hole
#'   both reduce the summed area, but only a hole reduces the union. Requires
#'   `area = TRUE`, which is when the union is computed.
#' @param area_tol_frac Numeric. Fractional tolerance for the area-conservation
#'   and non-overlap checks. Default `0.005` (0.5%).
#'
#' @return Invisibly a list with `ok`, `stage`, and `checks`.
#' @examples
#' \dontrun{
#' check_invariants("ngen",
#'                  flowpaths = sf::st_read("ngen_vpu01.gpkg", "flowpaths"),
#'                  divides   = sf::st_read("ngen_vpu01.gpkg", "divides"),
#'                  nexus     = sf::st_read("ngen_vpu01.gpkg", "nexus"),
#'                  strict = FALSE, area = TRUE)
#' }
#' @importFrom sf st_area st_length st_intersection st_sfc st_crs st_transform
#'   st_is_empty st_geometry st_geometry_type st_union
#' @importFrom igraph graph_from_data_frame topo_sort as_ids V is_dag
#' @importFrom lwgeom st_endpoint st_startpoint
#' @importFrom utils head
#' @export
check_invariants <- function(stage, ..., strict = TRUE,
                             coverage = FALSE, coverage_min = 0.90,
                             area = FALSE, area_baseline_m2 = NULL,
                             coverage_baseline_m2 = NULL,
                             area_tol_frac = 0.005) {
  stage <- match.arg(stage, c("refactored", "reconciled", "aggregated", "ngen"))
  args <- list(...)
  checks <- switch(stage,
    refactored = .hf_check_refactored(args),
    reconciled = .hf_check_reconciled(args, coverage = coverage,
      coverage_min = coverage_min),
    aggregated = .hf_check_aggregated(args, coverage = coverage,
      coverage_min = coverage_min),
    ngen       = .hf_check_ngen(args, coverage = coverage,
      coverage_min = coverage_min))

  # Divide-area integrity checks: appended for every stage that carries a
  # `divides` layer. These guard the two guarantees that the rest of the
  # pipeline can silently erode: total drainage-area conservation and the
  # dendritic 1-catchment-per-flowpath rule.
  if (stage %in% c("reconciled", "aggregated", "ngen"))
    checks <- c(checks,
      .hf_area_checks(args$divides, args$flowpaths,
        area_baseline_m2     = area_baseline_m2,
        coverage_baseline_m2 = coverage_baseline_m2,
        area_tol_frac        = area_tol_frac,
        overlap              = area))

  .hf_report_checks(checks, stage, strict, .hf_input_summary(args))
}

.hf_check_refactored <- function(args) {
  refactored <- args$refactored
  reconciled <- args$reconciled
  checks <- list()

  if (!is.null(refactored)) {
    ids <- as.character(refactored$flowpath_id)
    checks$refactored_unique <- .hf_ok(
      sum(duplicated(ids)) == 0,
      sprintf("refactored flowpath_id has %d duplicate(s)", sum(duplicated(ids))))
    toids <- as.character(refactored$flowpath_toid)
    dangling <- setdiff(toids[toids != "0" & !is.na(toids)], ids)
    checks$refactored_topology <- .hf_ok(
      length(dangling) == 0,
      sprintf("%d refactored flowpath_toid values reference nonexistent IDs",
        length(dangling)))
  }

  if (!is.null(reconciled)) {
    rec_ids <- as.character(reconciled$ID)
    checks$reconciled_unique <- .hf_ok(
      sum(duplicated(rec_ids)) == 0,
      sprintf("reconciled ID has %d duplicate(s)", sum(duplicated(rec_ids))))

    mems <- lapply(reconciled$member_flowpath_id, function(m) {
      if (is.null(m) || length(m) == 0) return(character(0))
      s <- if (is.list(m)) unlist(m) else as.character(m)
      trimws(unlist(strsplit(s, ",", fixed = TRUE)))
    })
    n_empty <- sum(lengths(mems) == 0L)
    checks$reconciled_has_members <- .hf_ok(
      n_empty == 0L,
      sprintf("%d reconciled group(s) have no member flowlines", n_empty))

    if (!is.null(refactored)) {
      all_mems <- unique(unlist(mems))
      rfp_ids  <- as.character(refactored$flowpath_id)
      base_ids <- unique(sub("\\.[0-9]+$", "", rfp_ids))
      valid    <- union(rfp_ids, base_ids)
      missing  <- setdiff(all_mems, valid)
      checks$reconciled_members_exist <- .hf_ok(
        length(missing) == 0L,
        sprintf("%d reconciled member ID(s) not in refactored or split-parents (e.g. %s)",
          length(missing), paste(utils::head(missing, 3), collapse = ", ")))

      # Inverse: every refactored flowpath_id must appear in some reconciled group
      missing_inv <- setdiff(rfp_ids, all_mems)
      checks$refactored_fully_reconciled <- .hf_ok(
        length(missing_inv) == 0L,
        sprintf("%d refactored flowpath_id(s) not in any reconciled group (e.g. %s)",
          length(missing_inv),
          paste(utils::head(missing_inv, 3), collapse = ", ")))
    }
  }
  checks
}

.hf_check_reconciled <- function(args, coverage, coverage_min) {
  reconciled <- args$reconciled
  divides    <- args$divides
  checks <- list()

  if (!is.null(reconciled) && inherits(reconciled, "sf") &&
    (any(c("toID", "rec_toid") %in% names(reconciled))) && nrow(reconciled) > 0L) {
    # Reconciled reaches are stored as single-path MULTILINESTRING; merge to
    # LINESTRING so start/end nodes are well-defined before checking flow
    # direction and disconnects. Only genuinely branching multi-part geometry
    # (not a clean flowpath) is excluded.
    merged_geom <- sf::st_line_merge(sf::st_geometry(reconciled))
    ls_rows     <- as.character(sf::st_geometry_type(merged_geom)) == "LINESTRING"
    if (any(ls_rows)) {
      rec_ls     <- reconciled[ls_rows, ]
      sf::st_geometry(rec_ls) <- merged_geom[ls_rows]
      rec_ids    <- as.character(if ("ID"     %in% names(rec_ls)) rec_ls$ID     else rec_ls$rec_id)
      rid_to_idx <- stats::setNames(seq_len(nrow(rec_ls)), rec_ids)
      rtoid      <- as.character(if ("toID"   %in% names(rec_ls)) rec_ls$toID   else rec_ls$rec_toid)
      has_ds     <- !is.na(rtoid) & rtoid != "0" & rtoid %in% names(rid_to_idx)
      if (any(has_ds)) {
        eps <- sf::st_coordinates(
          lwgeom::st_endpoint(sf::st_geometry(rec_ls)))[, 1:2]
        sps <- sf::st_coordinates(
          lwgeom::st_startpoint(sf::st_geometry(rec_ls)))[, 1:2]
        ds_idx  <- rid_to_idx[rtoid[has_ds]]
        d_end   <- sqrt(rowSums((eps[has_ds, , drop = FALSE] -
          sps[ds_idx, , drop = FALSE])^2))
        d_start <- sqrt(rowSums((sps[has_ds, , drop = FALSE] -
          sps[ds_idx, , drop = FALSE])^2))
        # A reach is REVERSED only if its START connects to the downstream's
        # start (within conn_tol) while its END does not: i.e. flipping it would
        # fix the connection. Reaches where NEITHER end is near the downstream are
        # DISCONNECTED (a geometry gap), not reversed; a bare d_start < d_end would
        # mislabel that gap as a reversal. conn_tol matches the flag threshold and
        # .fix_reconciled_flowdir's connection tolerance.
        conn_tol     <- 100
        flagged      <- which(d_end > 100)
        reversed     <- sum(d_start[flagged] <= conn_tol &
          d_start[flagged] < d_end[flagged])
        disconnected <- length(flagged) - reversed
        checks$reconciled_flow_direction <- .hf_ok(
          reversed == 0L,
          sprintf("%d flowline(s) reversed; %d disconnected (geometry gap)",
            reversed, disconnected))
      }
    }
  }

  if (is.null(divides)) return(checks)

  is_lnd   <- grepl("^lnd-", as.character(divides$divide_id)) |
    is.na(divides$flowpath_id)
  net_divs <- divides[!is_lnd, ]
  lnd_divs <- divides[is_lnd, ]

  if (!is.null(reconciled)) {
    rec_ids <- as.character(reconciled$ID)
    if ("has_divide" %in% names(reconciled))
      rec_ids <- rec_ids[is.na(reconciled$has_divide) | as.logical(reconciled$has_divide)]
    div_fp  <- as.character(net_divs$flowpath_id)
    missing <- setdiff(rec_ids, div_fp)
    checks$every_flowline_has_divide <- .hf_ok(
      length(missing) == 0L,
      sprintf("%d reconciled flowline(s) have no matching network divide (e.g. %s)",
        length(missing), paste(utils::head(missing, 5), collapse = ", ")))
  }

  if (nrow(net_divs) > 0L) {
    areas  <- suppressWarnings(as.numeric(sf::st_area(net_divs)))
    n_tiny <- sum(!is.finite(areas) | areas < 1, na.rm = TRUE)
    checks$network_divides_have_area <- .hf_ok(
      n_tiny == 0L,
      sprintf("%d network divide(s) have < 1 m2 area", n_tiny))
  }

  n_rings  <- .hf_ring_counts(sf::st_geometry(net_divs))
  max_rings <- if (length(n_rings)) max(n_rings) else 0L
  checks$no_excessive_rings <- .hf_ok(
    max_rings <= 50L,
    sprintf("max interior rings in a single divide: %d (threshold: 50)", max_rings))

  checks$landscape_preserved <- .hf_info(
    sprintf("landscape divides carried: %d", nrow(lnd_divs)))

  # POI contract at reconciled grain. Reported alongside the aggregated-grain
  # figure so a POI that ends up off its outlet can be attributed to stage 1 or
  # to aggregation, rather than only being visible in the product.
  # Absence is reported rather than skipped: a check that quietly disappears when
  # its input is missing looks identical to a passing one in the record.
  if (is.null(args$pois) || is.null(reconciled)) {
    checks$poi_at_outlet_reconciled <- .hf_info(
      "no POI layer supplied at reconciled grain -- POI position not measured")
  } else {
    rid <- if ("ID" %in% names(reconciled)) "ID" else "rec_id"
    dec <- .hf_poi_decompose(args$pois, reconciled, line_id_col = rid)
    checks$poi_at_outlet_reconciled <- if (is.null(dec))
      .hf_info("POI position at reconciled grain could not be computed")
      else .hf_poi_position_info(dec, "reconciled",
        types = if ("hl_type" %in% names(args$pois)) args$pois$hl_type else NULL)
  }

  checks
}

.hf_check_aggregated <- function(args, coverage, coverage_min) {
  flowpaths <- args$flowpaths
  divides   <- args$divides
  network   <- args$network
  checks <- list()



  # -- rec_id_partition -------------------------------------------------------
  # Every reconciled reach belongs to exactly one aggregate. This is the property
  # construct_agg_divides' exactness rests on: it dissolves rec_div by member, so a
  # reach claimed twice would double-count its catchment and a reach claimed by none
  # would drop it. Asserted here rather than enforced downstream: Step
  # 4_dedup_recid used to repair it and reported "partition OK" on all 24 VPUs,
  # so the property holds and only wanted stating.
  if (!is.null(network) && all(c("flowpath_id", "rec_id_list") %in% names(network))) {
    nu  <- unique(data.frame(
      fp = as.character(network$flowpath_id),
      rl = as.character(network$rec_id_list), stringsAsFactors = FALSE))
    recs <- unlist(lapply(nu$rl, function(x)
      trimws(strsplit(x, ",", fixed = TRUE)[[1]])), use.names = FALSE)
    recs <- recs[nzchar(recs) & !is.na(recs)]
    n_dup <- sum(duplicated(recs))
    checks$rec_id_partition <- .hf_ok(n_dup == 0L, sprintf(
      "%d reconciled reach(es) claimed by more than one aggregate (of %d member links)",
      n_dup, length(recs)))
  }

  # -- aggregated_flowpath_single_part ----------------------------------------
  # Reference flowlines are all LINESTRING, so a derived line should be single-part.
  # A multi-part aggregate line means its member reaches are not contiguous, which
  # Step 4a0 used to report as "non-contiguous members suggest a topology issue"
  # before reducing the line to one component and discarding the rest: 417 km
  # across 24 VPUs, from geometry stage-3 extend then rebuilt from members anyway.
  # Keeping the diagnostic without the discard: INFO, because it is non-zero and the
  # residual is genuine reference disconnection, not something to normalise away.
  if (!is.null(flowpaths) && inherits(flowpaths, "sf")) {
    np <- tryCatch({
      g <- sf::st_geometry(flowpaths)
      vapply(seq_along(g), function(i) {
        gi <- g[[i]]
        if (inherits(gi, "MULTILINESTRING")) length(gi) else 1L
      }, integer(1L))
    }, error = function(e) NULL)
    if (!is.null(np) && length(np))
      checks$aggregated_flowpath_single_part <- .hf_info(sprintf(
        "%d of %d aggregated flowpath(s) are multi-part (%.2f%%); max %d part(s)",
        sum(np > 1L), length(np), 100 * sum(np > 1L) / length(np), max(np)))
  }

  # -- aggregated_toid_resolves -----------------------------------------------
  # Every non-terminal flowpath_toid names a surviving flowpath, and no flowpath
  # drains to itself. Steps 4_toid and 4_selfloop used to repair these and reported
  # "toids OK" and "no self-loop toids" on all 24 VPUs.
  if (!is.null(flowpaths) &&
      all(c("flowpath_id", "flowpath_toid") %in% names(flowpaths))) {
    ids   <- as.character(flowpaths$flowpath_id)
    toids <- as.character(flowpaths$flowpath_toid)
    tnum  <- suppressWarnings(as.numeric(toids))
    term  <- is.na(toids) | toids == "0" | (!is.na(tnum) & tnum >= 1e9)
    n_dang <- sum(!term & !(toids %in% ids))
    n_self <- sum(!is.na(toids) & toids != "0" & ids == toids)
    checks$aggregated_toid_resolves <- .hf_ok(n_dang == 0L && n_self == 0L, sprintf(
      "%d dangling non-terminal flowpath_toid(s) + %d self-loop(s)", n_dang, n_self))
  }

  if (!is.null(flowpaths)) {
    ids <- as.character(flowpaths$flowpath_id)
    checks$aggregated_flowpath_unique <- .hf_ok(
      sum(duplicated(ids)) == 0L,
      sprintf("aggregated flowpath_id has %d duplicate(s)", sum(duplicated(ids))))
  }

  if (!is.null(divides)) {
    dids <- as.character(divides$divide_id)
    checks$aggregated_divide_unique <- .hf_ok(
      sum(duplicated(dids)) == 0L,
      sprintf("aggregated divide_id has %d duplicate(s)", sum(duplicated(dids))))
  }

  # -- aggregated_divide_link_resolves ----------------------------------------
  # A flowpath that claims a divide must name one that exists. When it does not,
  # the catchment is not merely misplaced: there is no polygon to place. Stage
  # 3 then clears the unresolvable id and the flowpath ships as a real river
  # owning no land, which is what a consumer sees as a NULL divide.
  #
  # Stated here because this is where it is decided. The finished product is
  # always clean: Step 6d5 nulls a dangling divide_id, so an ngen-grain check
  # would pass while the loss had already happened. VPU 04 carries 69 of these,
  # 377 km of channel, the longest 94.9 km.
  if (!is.null(flowpaths) && !is.null(divides) &&
      all(c("divide_id", "has_divide") %in% names(flowpaths)) &&
      "divide_id" %in% names(divides)) {
    fdid <- as.character(flowpaths$divide_id)
    hd   <- flowpaths$has_divide %in% c(TRUE, 1L, "1", "TRUE")
    claim <- hd & !is.na(fdid) & nzchar(fdid)
    dang  <- claim & !fdid %in% as.character(divides$divide_id)
    n_dang <- sum(dang, na.rm = TRUE)
    ex <- utils::head(as.character(flowpaths$flowpath_id)[which(dang)], 3)
    checks$aggregated_divide_link_resolves <- {
      msg <- sprintf(paste0("%d of %d flowpath(s) claiming a divide name a ",
        "divide_id with no divide row (target 0)%s"),
        n_dang, sum(claim, na.rm = TRUE),
        if (length(ex)) sprintf(" (e.g. %s)", paste(ex, collapse = ", ")) else "")
      if (n_dang == 0L) .hf_ok(TRUE, msg) else .hf_info(msg)
    }
  }

  if (coverage && !is.null(flowpaths) && !is.null(divides))
    checks$flowpath_coverage <- .hf_coverage_check(flowpaths, divides, coverage_min)

  checks
}

.hf_check_ngen <- function(args, coverage, coverage_min) {
  flowpaths <- args$flowpaths
  divides   <- args$divides
  nexus     <- args$nexus
  network   <- args$network
  checks <- list()

  # -- total_da_accumulates -------------------------------------------------
  # The delineated drainage area is built by accumulating each flowline's own
  # divide area down the flowline graph, so at every reach
  #   total = incremental + sum(total of the reaches draining into it)
  # must hold exactly. That identity is what makes the number checkable
  # without an external baseline, and it is worth asserting rather than
  # trusting: the value is written at Step 6d8 and several later steps rewrite
  # the flowlines layer.
  #
  # It replaces a comparison against the reference `flowline_total_dasqkm`,
  # which cannot serve as a baseline: that column is a reference attribute
  # carried onto RECONCILED reaches, so a split fragment inherits its parent's
  # figure and stops describing the reach it sits on. Measured on VPU 01, it
  # agrees with the delineated area within 10% on 26.5% of flowlines and reads
  # under a tenth of it on 24.5%, so a tolerance against it measures the
  # carrying rather than the delineation.
  fl_da <- args$flowlines
  if (!is.null(fl_da) && nrow(fl_da) &&
      all(c("flowline_id", "flowline_toid", "total_dasqkm",
            "incremental_areasqkm") %in% names(fl_da))) {
    fid <- as.character(fl_da$flowline_id)
    fto <- as.character(fl_da$flowline_toid)
    tot <- suppressWarnings(as.numeric(fl_da$total_dasqkm))
    inc <- suppressWarnings(as.numeric(fl_da$incremental_areasqkm))
    inc[!is.finite(inc)] <- 0
    up_sum <- rep(0, length(fid))
    recv   <- match(fto, fid)
    keep   <- !is.na(recv) & is.finite(tot)
    if (any(keep)) {
      agg <- tapply(tot[keep], recv[keep], sum)
      up_sum[as.integer(names(agg))] <- as.numeric(agg)
    }
    # Scale the tolerance with the value: these are sums of thousands of
    # floating-point areas, so a fixed epsilon fails on large basins alone.
    expect <- inc + up_sum
    ok_da  <- is.finite(tot) & is.finite(expect) &
      abs(tot - expect) <= pmax(1e-6, 1e-9 * pmax(abs(tot), abs(expect)))
    n_bad  <- sum(is.finite(tot) & !ok_da)
    worst  <- if (n_bad) max(abs(tot - expect)[is.finite(tot) & !ok_da]) else 0
    checks$total_da_accumulates <- .hf_ok(n_bad == 0L,
      sprintf(paste0("%d flowline(s) where total_dasqkm != incremental + ",
        "sum(upstream totals) (worst %.6f km2)"), n_bad, worst))
  }

  if (!is.null(flowpaths)) {
    ids   <- as.character(flowpaths$flowpath_id)
    n_bad <- sum(!is.na(ids) & !startsWith(ids, "fp-"))
    checks$ngen_fp_prefix <- .hf_ok(n_bad == 0L,
      sprintf("%d ngen flowpath_id(s) do not start with 'fp-'", n_bad))

    # -- flowpath_single_part -------------------------------------------------
    # Reference flowlines are all LINESTRING, so every derived line layer must be
    # a single-part line too. A multi-part flowpath has no continuous
    # inlet-to-outlet path, which makes its flow direction, length and outlet
    # position ambiguous, and routing needs one unbroken path per feature.
    #
    # INFO, not a gate: a small remainder is genuinely multi-part, and driving it
    # to zero is tracked work. Most of the population is single-part wearing a
    # MULTILINESTRING type.
    #
    # Never cast a multi-part line to its longest part, that discards geometry
    # and hides the defect rather than resolving it.
    if (inherits(flowpaths, "sf")) {
      # Count parts per feature explicitly. lengths() is not usable here: on an
      # sfc_LINESTRING each element is a coordinate matrix, so lengths() returns
      # the coordinate count and every single-part line would read as multi-part.
      np <- tryCatch({
        g <- sf::st_geometry(flowpaths)
        vapply(seq_along(g), function(i) {
          gi <- g[[i]]
          if (inherits(gi, "MULTILINESTRING")) length(gi) else 1L
        }, integer(1L))
      }, error = function(e) NULL)
      if (!is.null(np) && length(np)) {
        gt      <- as.character(sf::st_geometry_type(flowpaths))
        n_multi <- sum(np > 1L)
        checks$flowpath_single_part <- .hf_info(sprintf(
          "%d of %d flowpath(s) are multi-part (%.2f%%); max %d part(s); stored type(s): %s",
          n_multi, length(np), 100 * n_multi / length(np), max(np),
          paste(unique(gt), collapse = "/")))
      }
    }

    # -- flowpath_flow_direction ----------------------------------------------
    # Aggregated flowpath geometry must run upstream->downstream (outlet = last
    # vertex). extend_flowpaths_to_mainstem's line_merge can flip orientation and
    # only the reconciled flowlines are otherwise checked, so validate the
    # flowpaths layer here using the reliably-oriented flowlines to locate each
    # outlet. Step 6l used to re-orient here; it reversed 0 on all 24 VPUs and was
    # removed, so this is now the assertion rather than a report on a repair.
    fl_dir <- args$flowlines
    if (!is.null(fl_dir) && inherits(fl_dir, "sf") && nrow(fl_dir) > 0L &&
        all(c("flowline_id", "flowline_toid", "flowpath_id") %in% names(fl_dir))) {
      fl2 <- fl_dir
      fl2$flowline_id   <- as.character(fl2$flowline_id)
      fl2$flowline_toid <- as.character(fl2$flowline_toid)
      fl2$flowpath_id   <- as.character(fl2$flowpath_id)
      if (sf::st_crs(fl2) != sf::st_crs(flowpaths))
        fl2 <- sf::st_transform(fl2, sf::st_crs(flowpaths))
      fl_by <- split(seq_len(nrow(fl2)), fl2$flowpath_id)
      tda   <- suppressWarnings(as.numeric(fl2$tot_drainage_areasqkm)); tda[is.na(tda)] <- 0
      fl_end <- sf::st_coordinates(lwgeom::st_endpoint(sf::st_geometry(fl2)))[, 1:2, drop = FALSE]
      gm    <- sf::st_line_merge(
        sf::st_cast(sf::st_geometry(flowpaths), "MULTILINESTRING", warn = FALSE))
      is_ls <- as.character(sf::st_geometry_type(gm)) == "LINESTRING" &
        as.numeric(sf::st_length(gm)) > 0
      reversed <- 0L
      if (any(is_ls)) {
        sp <- sf::st_coordinates(lwgeom::st_startpoint(gm[is_ls]))[, 1:2, drop = FALSE]
        ep <- sf::st_coordinates(lwgeom::st_endpoint(gm[is_ls]))[, 1:2, drop = FALSE]
        wk <- which(is_ls); j <- 0L
        for (k in wk) {
          j <- j + 1L
          r <- fl_by[[ids[k]]]; if (is.null(r)) next
          fid <- fl2$flowline_id[r]
          o   <- r[!(fl2$flowline_toid[r] %in% fid)]
          orow <- if (length(o)) o[which.max(tda[o])] else r[which.max(tda[r])]
          ov <- fl_end[orow, ]
          if (sum((sp[j, ] - ov)^2) < sum((ep[j, ] - ov)^2)) reversed <- reversed + 1L
        }
      }
      checks$flowpath_flow_direction <- .hf_ok(reversed == 0L,
        sprintf("%d flowpath(s) reversed (outlet is NOT the last vertex)", reversed))
    }

    # -- flowpath_touches_network ---------------------------------------------
    # Every non-terminal flowpath's line should TOUCH the flowpath network (share
    # a point with another flowpath, anywhere: end node OR mid-reach). A line
    # that touches nothing is a genuine geometric disconnect. Uses a 30m tolerance
    # (the DEM/snap resolution): sub-30m near-misses are noise, only >30m gaps
    # count. INFO, not gated: the residual gaps are inherited from the reference
    # hydrography; reported per build so an increase is visible.
    if (inherits(flowpaths, "sf") && nrow(flowpaths) > 0L &&
        "flowpath_toid" %in% names(flowpaths)) {
      n_float <- tryCatch({
        fpm  <- if (isTRUE(sf::st_is_longlat(flowpaths)))
          sf::st_transform(flowpaths, 5070) else flowpaths
        geo  <- sf::st_geometry(fpm)
        toid <- as.character(flowpaths$flowpath_toid)
        term <- grepl("^tnx", toid) | is.na(toid) | toid == "0"
        # within-30m via 30m buffer intersect (indexed; avoids slow pairwise dist)
        near <- sf::st_intersects(sf::st_buffer(geo, 30), geo)
        touch <- vapply(seq_along(near), function(i) any(near[[i]] != i), logical(1))
        sum(!touch & !term)
      }, error = function(e) NA_integer_)
      if (!is.na(n_float))
        checks$flowpath_touches_network <- if (n_float == 0L)
          .hf_ok(TRUE, "every non-terminal flowpath line touches the network (<=30m)")
        else
          .hf_info(sprintf(
            "%d non-terminal flowpath(s) whose line touches no other flowpath within 30m -- geometric disconnect (reference-hydrography gap)",
            n_float))
    }

    # -- upstream_index_monotone ----------------------------------------------
    # When the nested-set upstream index is present (add_upstream_index), it is a
    # near-free cross-check on topology and drainage-area accumulation: walking
    # one edge downstream, the upstream count must strictly increase (the
    # downstream flowpath sees this one plus its whole up-tree), and total
    # drainage area must not decrease. A violation means the index, the topology,
    # or the accumulated area is inconsistent. Skipped when the index is absent.
    if ("num_upstreams" %in% names(flowpaths) &&
        "flowpath_toid" %in% names(flowpaths) && !is.null(nexus)) {
      fpid   <- as.character(flowpaths$flowpath_id)
      numu   <- suppressWarnings(as.integer(flowpaths$num_upstreams))
      nexmap <- stats::setNames(as.character(nexus$nexus_toid),
        as.character(nexus$nexus_id))
      ds     <- unname(nexmap[as.character(flowpaths$flowpath_toid)])
      di     <- match(ds, fpid)                    # downstream row (NA = terminal)
      has_dn <- !is.na(di) & !is.na(numu) & !is.na(numu[di])
      bad_nu <- sum(has_dn & !(numu[di] > numu))
      bad_da <- 0L
      if ("total_dasqkm" %in% names(flowpaths)) {
        da  <- suppressWarnings(as.numeric(flowpaths$total_dasqkm))
        tol <- 1e-6
        bad_da <- sum(has_dn & is.finite(da) & is.finite(da[di]) & (da[di] + tol < da))
      }
      checks$upstream_index_monotone <- .hf_ok(bad_nu == 0L && bad_da == 0L,
        sprintf(paste0("%d edge(s) where num_upstreams does not increase ",
          "downstream; %d edge(s) where total_dasqkm decreases downstream"),
          bad_nu, bad_da))
    }

    # -- slope_valid ----------------------------------------------------------
    # Channel routing needs a strictly positive slope on every flowpath. Flag
    # negative, zero, or null slopes. When the column is entirely NA the slope
    # step was not run for this build (info, not a failure); once any slope is
    # present, every value must be > 0.
    if ("slope" %in% names(flowpaths)) {
      sl <- suppressWarnings(as.numeric(flowpaths$slope))
      if (all(is.na(sl))) {
        checks$slope_valid <- .hf_info(
          "slope not computed (all NA) for this build")
      } else {
        n_na  <- sum(is.na(sl))
        n_neg <- sum(sl < 0, na.rm = TRUE)
        n_zero<- sum(sl == 0, na.rm = TRUE)
        checks$slope_valid <- .hf_ok(
          n_na == 0L && n_neg == 0L && n_zero == 0L,
          sprintf("%d null + %d negative + %d zero slope(s) (routing needs slope > 0)",
            n_na, n_neg, n_zero))
      }
    }

    # -- So_valid -------------------------------------------------------------
    # Channel bottom slope `So` (routing attribute) must be strictly positive on
    # every feature. Source `slope` can be <=0 on reversed or degenerate segments;
    # .calculate_flow_atts floors it, so any non-positive So in the output is a
    # defect. Check wherever So lands (flowpaths or network attribute table).
    So_src <- c(
      if ("So" %in% names(flowpaths)) suppressWarnings(as.numeric(flowpaths$So)),
      if (!is.null(args$flowlines) && "So" %in% names(args$flowlines))
        suppressWarnings(as.numeric(args$flowlines$So)),
      if (!is.null(network) && "So" %in% names(network))
        suppressWarnings(as.numeric(network$So)))
    if (length(So_src) > 0L) {
      so <- So_src
      n_na  <- sum(is.na(so))
      n_neg <- sum(so < 0,  na.rm = TRUE)
      n_zero<- sum(so == 0, na.rm = TRUE)
      checks$So_valid <- .hf_ok(
        n_na == 0L && n_neg == 0L && n_zero == 0L,
        sprintf("%d null + %d negative + %d zero So / channel-slope value(s) (routing needs So > 0)",
          n_na, n_neg, n_zero))
    }
  }

  if (!is.null(divides)) {
    dids  <- as.character(divides$divide_id)
    n_bad <- sum(!is.na(dids) & !startsWith(dids, "cat-") & !startsWith(dids, "lnd-"))
    checks$ngen_cat_prefix <- .hf_ok(n_bad == 0L,
      sprintf("%d ngen divide_id(s) are not 'cat-' or 'lnd-' prefixed", n_bad))

    # -- network_divide_has_owner ---------------------------------------------
    # Every cat- divide must name an owning flowpath. This closes a silent hole
    # rather than restating divide_flowpath_consistent: .hf_network_divides
    # classifies a divide as landscape when its flowpath_id is NA, regardless of
    # its cat- prefix, so an unowned cat- divide is quietly reclassified out of
    # every area and coverage figure. Its land leaves the books with no geometry
    # lost and no check failing.
    #
    # divide_flowpath_consistent cannot see it: that check works from the
    # flowpath side, so an unowned divide no flowpath references is invisible to
    # it.
    if ("flowpath_id" %in% names(divides)) {
      own      <- as.character(divides$flowpath_id)
      is_cat   <- !is.na(dids) & startsWith(dids, "cat-")
      unowned  <- is_cat & (is.na(own) | !nzchar(own))
      n_unown  <- sum(unowned)
      a_unown  <- if (n_unown > 0L && "areasqkm" %in% names(divides))
        sum(suppressWarnings(as.numeric(divides$areasqkm[unowned])), na.rm = TRUE)
      else 0
      checks$network_divide_has_owner <- .hf_ok(n_unown == 0L, sprintf(
        "%d cat- divide(s) have no owning flowpath (%.4f km2 excluded from coverage)%s",
        n_unown, a_unown,
        if (n_unown > 0L)
          sprintf(" -- e.g. %s", paste(utils::head(dids[unowned], 3), collapse = ", "))
        else ""))
    }
  }

  if (!is.null(flowpaths) && !is.null(divides)) {
    # The flowpath<->divide link is stored redundantly in three places, written by
    # separate steps: flowpaths.divide_id, divides.flowpath_id, and network.divide_id.
    # All three must describe the same bijection. A disagreement is contested
    # ownership (a divide claimed by more than one flowpath), which is what breaks
    # the cross-VPU merge, so this gates rather than reports: the condition is
    # asserted here instead of being silently repaired downstream.
    fpd <- ifelse(is.na(flowpaths$divide_id) | flowpaths$divide_id == "", NA,
      as.character(flowpaths$divide_id))
    names(fpd) <- as.character(flowpaths$flowpath_id)
    dvf <- divides[!is.na(divides$flowpath_id) & divides$flowpath_id != "", ]
    dv_by_fp <- stats::setNames(as.character(dvf$divide_id), as.character(dvf$flowpath_id))
    aa <- as.character(flowpaths$flowpath_id)
    xx <- fpd[aa]
    yy <- dv_by_fp[aa]
    disagree <- function(a, b) {
      xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b)
    }
    n_mis <- sum(disagree(xx, yy))

    # Third copy: the aspatial network table. Compared per flowpath rather than
    # per row, since a flowpath owns many network rows that must all agree.
    n_net <- 0L
    if (!is.null(network) && all(c("flowpath_id", "divide_id") %in% names(network))) {
      nf <- as.character(network$flowpath_id)
      nd <- ifelse(is.na(network$divide_id) | network$divide_id == "", NA,
        as.character(network$divide_id))
      keep  <- !is.na(nf) & nf %in% aa
      n_net <- sum(disagree(nd[keep], fpd[nf[keep]]))
    }

    checks$divide_flowpath_consistent <- .hf_ok(n_mis == 0L && n_net == 0L, sprintf(
      paste0("%d flowpath(s) whose divide_id disagrees with divides.flowpath_id ",
        "+ %d network row(s) disagreeing with flowpaths.divide_id ",
        "(contested ownership -- construct_agg_divides)"),
      n_mis, n_net))
  }

  if (!is.null(flowpaths) && !is.null(nexus)) {
    nex_lu  <- stats::setNames(as.character(nexus$nexus_toid),
      as.character(nexus$nexus_id))
    next_fp <- nex_lu[as.character(flowpaths$flowpath_toid)]
    edge_df <- data.frame(
      flowpath_id   = as.character(flowpaths$flowpath_id),
      flowpath_toid = ifelse(is.na(next_fp), "0", as.character(next_fp)),
      stringsAsFactors = FALSE)
    is_dag_val <- tryCatch(
      hfutils::hf_network_is_dag(edge_df, id_col = "flowpath_id", toid_col = "flowpath_toid"),
      error = function(e) NA)
    if (isTRUE(is_dag_val)) {
      checks$network_is_dag <- .hf_ok(TRUE,
        sprintf("fp->nexus->fp graph is a DAG (%d flowpaths, %d connectors)",
          nrow(flowpaths), sum(!flowpaths$has_divide, na.rm = TRUE)))
    } else {
      checks$network_is_dag <- .hf_ok(FALSE, "fp->nexus->fp graph contains cycle(s)")
    }

    # -- no_divideless_leaf_flowpaths (tripwire) ------------------------------
    # A divide-less (has_divide=FALSE) headwater leaf that is neither a POI nor a
    # lake reach should have been lifted out of the flowpath grain by stage-3
    # Step 6e2 (.demote_divideless_leaf_flowpaths). A survivor means a
    # catchment-less flowpath leaked into the model. Interior divide-less reaches
    # (with upstream inflow) are legitimate connectors and excluded here.
    hd <- if ("has_divide" %in% names(flowpaths)) as.logical(flowpaths$has_divide)
      else !is.na(flowpaths$divide_id) & nzchar(as.character(flowpaths$divide_id))
    hd[is.na(hd)] <- FALSE
    poi_f  <- "poi_id"  %in% names(flowpaths) &
      !is.na(flowpaths$poi_id)  & nzchar(as.character(flowpaths$poi_id))
    lake_f <- "lake_id" %in% names(flowpaths) &
      !is.na(flowpaths$lake_id) & nzchar(as.character(flowpaths$lake_id))
    inflow_nex <- stats::setNames(as.character(nexus$nexus_id),
      as.character(nexus$nexus_toid))
    donors <- table(as.character(flowpaths$flowpath_toid))
    fin    <- unname(inflow_nex[as.character(flowpaths$flowpath_id)])
    nup    <- as.integer(donors[fin]); nup[is.na(nup)] <- 0L
    n_leaf_conn <- sum(!hd & !poi_f & !lake_f & nup == 0L)
    checks$no_divideless_leaf_flowpaths <- .hf_ok(
      n_leaf_conn == 0L,
      sprintf(paste0("%d divide-less headwater-leaf flowpath(s) survived ",
        "(should be demoted to flowline-only)"), n_leaf_conn))

    # -- flowpath_outlet_at_nexus ---------------------------------------------
    # NON-DIAGNOSTIC measurement (INFO only): reports the distance from each
    # flowpath's downstream endpoint to the nexus its flowpath_toid points to.
    # This gap is NOT a defect and does NOT converge to zero: an aggregated
    # flowpath's line is the merged mainstem stem, whose outlet is wherever that
    # stem geometrically ends; the nexus is a topological connector point and need
    # not sit on the flowpath endpoint. Kept only for visibility. The real
    # outlet-touches-inflow check is at the FLOWLINE grain
    # (flowline_toid_geometry_consistent): flowline outlets, not flowpath
    # outlets, are the ones that should coincide with their downstream reach.
    # Do NOT flip this to a hard check.
    tol_m <- 50
    has_geom <- inherits(flowpaths, "sf") && inherits(nexus, "sf") &&
      !is.null(sf::st_geometry(flowpaths)) && !is.null(sf::st_geometry(nexus))
    if (has_geom) {
      gaps <- tryCatch({
        nex_g  <- sf::st_geometry(nexus)
        nex_id <- as.character(nexus$nexus_id)
        # match each fp to its downstream nexus point
        mi <- match(as.character(flowpaths$flowpath_toid), nex_id)
        ok <- !is.na(mi)
        if (!any(ok)) numeric(0) else {
          fp_end <- lwgeom::st_endpoint(sf::st_geometry(flowpaths)[ok])
          as.numeric(sf::st_distance(fp_end, nex_g[mi[ok]], by_element = TRUE))
        }
      }, error = function(e) NULL)
      if (is.null(gaps)) {
        checks$flowpath_outlet_at_nexus <- .hf_info(
          "flowpath_outlet_at_nexus: gap not computable (geometry error)")
      } else if (length(gaps) == 0L) {
        checks$flowpath_outlet_at_nexus <- .hf_info(
          "flowpath_outlet_at_nexus: no fp->nexus links to measure")
      } else {
        n_bad <- sum(gaps > tol_m, na.rm = TRUE)
        checks$flowpath_outlet_at_nexus <- .hf_info(sprintf(
          "flowpath-outlet-to-nexus distance (expected non-zero, non-diagnostic): %d/%d > %dm (mean %.0fm, max %.0fm)",
          n_bad, length(gaps), tol_m, mean(gaps, na.rm = TRUE), max(gaps, na.rm = TRUE)))
      }
    }
  }

  # -- terminal_has_a_receiver ------------------------------------------------
  # A terminal flowpath whose outlet lies within a DEM cell of a channel draining
  # more than it does is water dropped beside a river already carrying its
  # neighbours: routing sees the flow arrive and never leave. OWP reported this
  # as tnx- flows not routed downstream.
  #
  # Reported rather than gated because the population it measures is not only
  # ours. Most terminals are genuine outlets, and many more are reaches whose
  # downstream neighbour is absent from the reference entirely, in VPU 04,
  # 2,948 of 2,999 terminals name a downstream by dnhydroseq that the input does
  # not contain. Neither is repairable here. Step 3a-pre reconnects the ones that
  # are, so this counts what is left; promote to a gate if it reaches zero.
  if (inherits(flowpaths, "sf") && !is.null(nexus) &&
      "flowpath_toid" %in% names(flowpaths)) {
    n_left <- tryCatch({
      tfp <- grepl("^tnx-", as.character(flowpaths$flowpath_toid))
      if (!any(tfp)) 0L else {
        da <- suppressWarnings(as.numeric(flowpaths$total_dasqkm %||%
          flowpaths$tot_drainage_areasqkm))
        ends <- lwgeom::st_endpoint(sf::st_geometry(flowpaths)[tfp])
        near <- suppressWarnings(sf::st_is_within_distance(
          ends, sf::st_geometry(flowpaths), dist = 30))
        wi <- which(tfp)
        sum(vapply(seq_along(wi), function(k) {
          js <- setdiff(near[[k]], wi[k])
          length(js) && any(is.finite(da[js]) & is.finite(da[wi[k]]) &
            da[js] > da[wi[k]])
        }, logical(1L)))
      }
    }, error = function(e) NA_integer_)
    checks$terminal_has_a_receiver <- if (is.na(n_left))
      .hf_info("terminal_has_a_receiver: not computable (geometry error)") else
      .hf_info(sprintf(paste0("%d terminal flowpath(s) end within 30m of a ",
        "channel draining more than they do (target 0)"), n_left))
  }

  # -- lake_membership_upstream_of_outlet -------------------------------------
  # Membership comes from the VAA `comid -> wbareacomi` mapping. What makes a
  # stamp checkable is the network, not distance: a lake's members should sit at
  # or upstream of the outlet reach RouteLink names for that lake. A Euclidean
  # bound cannot work here, the lakes are POINTS, so there is no extent to
  # scale one against.
  #
  # A gate, because the build guarantees it: unreachable members are dropped at
  # Step 6h and recorded in `lake_id_flags`. Both this check and that drop are
  # conditional on the upstream index being present, so they cannot disagree.
  #
  # Was briefly INFO, when the plan was to report these rather than remove them.
  # They are removed because the harm is not cosmetic: t-route collapses a
  # lake's members into a single node, so a member on a branch the network cannot
  # reach from the outlet fabricates a routing connection, the same class of
  # corruption as a lake-id gap, in a shape the contiguity check cannot see,
  # because the member is not BETWEEN anything.
  if (!is.null(flowpaths) && inherits(flowpaths, c("sf", "data.frame")) &&
      all(c("flowpath_id", "lake_id", "upstream_id", "num_upstreams") %in% names(flowpaths)) &&
      !is.null(args$lakes) && "flowpath_id" %in% names(args$lakes)) {
    res <- tryCatch({
      fpc <- as.character(flowpaths$flowpath_id)
      lkv <- as.character(flowpaths$lake_id)
      u   <- suppressWarnings(as.numeric(flowpaths$upstream_id))
      nu  <- suppressWarnings(as.numeric(flowpaths$num_upstreams))
      out <- stats::setNames(as.character(args$lakes$flowpath_id),
                             as.character(args$lakes$lake_id))
      k   <- which(!is.na(lkv) & nzchar(lkv))
      oi  <- match(unname(out[lkv[k]]), fpc)
      ok  <- !is.na(oi) & is.finite(u[k]) & is.finite(u[oi]) &
             u[k] >= u[oi] & u[k] <= u[oi] + nu[oi]
      list(bad = sum(!is.na(oi) & !ok), n = length(k))
    }, error = function(e) NULL)
    if (!is.null(res))
      checks$lake_membership_upstream_of_outlet <- .hf_ok(res$bad == 0L, sprintf(
        paste0("%d of %d lake member(s) are not upstream of their lake's ",
               "RouteLink outlet -- t-route collapses a lake into one node, so a ",
               "member on an unreachable branch invents a routing connection ",
               "(see `lake_id_flags`)"),
        res$bad, res$n))
  }

  # -- gage_within_one_feature -------------------------------------------------
  # A gage is the location the fabric is reported at, so it must stay where the
  # source put it: on the flowpath carrying the reach its `hl_flowpath_id_source`
  # names, or at most one feature up- or downstream of it.
  #
  # Topological distance is the test, not drainage-area agreement. A gage placed
  # one feature upstream reports that reach's area rather than the carrier's, and
  # that difference is a consequence of the outlet convention, not a defect,
  # comparing areas conflates it with a genuine misplacement. Distance separates
  # the two and needs no external data.
  #
  # INFO: a small residue sits further out, mostly two or three hops along one
  # channel where the area is close, and at least one on a different branch
  # entirely. Promoting this to a gate wants that population at zero first.
  if (!is.null(args$hydrolocations) && !is.null(flowpaths) &&
      all(c("hl_reference", "hl_flowpath_id_source", "flowpath_id") %in%
          names(args$hydrolocations)) &&
      all(c("flowpath_id", "flowpath_toid", "member_comid") %in% names(flowpaths)) &&
      !is.null(args$nexus) && "nexus_id" %in% names(args$nexus)) {
    res <- tryCatch({
      hl  <- args$hydrolocations
      fid <- as.character(flowpaths$flowpath_id)
      n2f <- stats::setNames(as.character(args$nexus$nexus_toid),
                             as.character(args$nexus$nexus_id))
      dn  <- match(unname(n2f[as.character(flowpaths$flowpath_toid)]), fid)
      up  <- split(seq_along(fid), dn)
      memb <- strsplit(as.character(flowpaths$member_comid), ",", fixed = TRUE)
      car  <- split(rep(seq_along(fid), lengths(memb)),
                    sub("[.][0-9]+$", "", trimws(unlist(memb))))
      isg <- grepl("(^|\\|)(nwis|gages)-", as.character(hl$hl_reference))
      far <- character(0); n <- 0L
      for (i in which(isg)) {
        att <- match(as.character(hl$flowpath_id[i]), fid)
        cs  <- car[[as.character(hl$hl_flowpath_id_source[i])]]
        if (is.na(att) || is.null(cs)) next
        n <- n + 1L
        if (att %in% cs) next
        one <- FALSE
        for (c1 in cs) {
          u <- up[[as.character(c1)]]
          if ((!is.null(u) && att %in% u) || (!is.na(dn[c1]) && dn[c1] == att)) {
            one <- TRUE; break
          }
        }
        if (!one) far <- c(far, as.character(hl$hl_reference[i]))
      }
      list(far = far, n = n)
    }, error = function(e) NULL)
    if (!is.null(res) && res$n > 0L)
      checks$gage_within_one_feature <- .hf_info(sprintf(
        paste0("%d of %d gage(s) sit more than one feature from the flowpath ",
               "carrying their assigned reach (%.2f%% within one)%s"),
        length(res$far), res$n, 100 * (1 - length(res$far) / res$n),
        if (length(res$far)) sprintf(" (e.g. %s)",
          paste(utils::head(res$far, 3), collapse = ", ")) else ""))
  }

  # -- lake_members_contiguous -------------------------------------------------
  # t-route collapses every flowpath sharing a lake_id into ONE waterbody node,
  # so those flowpaths must form a connected subgraph. A non-member sitting
  # between two members makes the collapse produce
  # node_L -> gap -> node_L: a cycle, and routing with waterbodies enabled fails.
  #
  # lake_id is stamped per flowpath by a length-weighted majority of its member
  # reaches, a rule that cannot see the network: a reach mostly outside the lake
  # is rejected even when both its neighbours are members.
  #
  # Checks the immediate case, a single non-member directly between two members
  # of the same lake, which is what the gap-closing step guarantees against.
  if (!is.null(flowpaths) && inherits(flowpaths, c("sf", "data.frame")) &&
      all(c("flowpath_id", "flowpath_toid", "lake_id") %in% names(flowpaths)) &&
      !is.null(args$nexus) && "nexus_id" %in% names(args$nexus)) {
    res <- tryCatch({
      fid  <- as.character(flowpaths$flowpath_id)
      lk   <- as.character(flowpaths$lake_id)
      lk[is.na(lk) | !nzchar(lk)] <- NA_character_
      n2f  <- stats::setNames(as.character(args$nexus$nexus_toid),
                              as.character(args$nexus$nexus_id))
      nxt  <- match(unname(n2f[as.character(flowpaths$flowpath_toid)]), fid)
      ups  <- split(seq_along(fid), nxt)
      bad  <- character(0)
      for (i in which(is.na(lk))) {
        dn <- nxt[i]
        if (is.na(dn) || is.na(lk[dn])) next
        u <- ups[[as.character(i)]]
        if (is.null(u)) next
        if (any(!is.na(lk[u]) & lk[u] == lk[dn])) bad <- c(bad, fid[i])
      }
      bad
    }, error = function(e) NULL)
    if (!is.null(res))
      checks$lake_members_contiguous <- .hf_ok(length(res) == 0L, sprintf(
        paste0("%d flowpath(s) break a lake into disconnected pieces -- a ",
               "non-member between two members of one lake, which t-route ",
               "collapses into a cycle%s (target 0)"),
        length(res),
        if (length(res)) sprintf(" (e.g. %s)",
          paste(utils::head(res, 3), collapse = ", ")) else ""))
  }

  # -- flowpath member contiguity and cross-grain topology agreement -----------
  # An aggregated flowpath is a contraction of a set of reconciled reaches. For
  # the contraction to mean anything the members must form ONE connected chain,
  # and the flowpath-grain topology (flowpath_toid through the nexus layer) must
  # name the same downstream flowpath the flowline-grain topology
  # (flowline_toid) actually reaches. Neither property was asserted anywhere,
  # which is how interleaved flowpaths, A owning reaches on both sides of B,
  # so the flowline walk gives A -> B -> A: shipped while every existing gate
  # passed. The declared flowpath graph is a DAG in that state; it just is not
  # the graph the flowlines describe.
  #
  # Measured on the ROUTING SPINE (flowpath_constituent) where available: the
  # non-constituent members are tributaries inside the catchment, which
  # legitimately leave the flowpath's drawn line and must not count as exits.
  if (!is.null(args$flowlines) &&
      all(c("flowline_id", "flowline_toid", "flowpath_id") %in% names(args$flowlines))) {
    res <- tryCatch({
      fl   <- args$flowlines
      flid <- as.character(fl$flowline_id)
      fpid <- as.character(fl$flowpath_id)
      dn   <- match(as.character(fl$flowline_toid), flid)
      ds_fp <- ifelse(is.na(dn), NA_character_, fpid[dn])
      spine <- if ("flowpath_constituent" %in% names(fl))
        !is.na(fl$flowpath_constituent) & as.logical(fl$flowpath_constituent)
        else rep(TRUE, length(flid))
      # An exit is a member whose downstream belongs to a DIFFERENT flowpath. A
      # member with no downstream at all (flowline_toid 0/NA, a network
      # terminal, or a dangling reference reach inside the catchment) is a chain
      # END, not a way out: counting it as an exit reports a second outlet for
      # every flowpath that happens to contain one. On VPU 01 that alone
      # inflated the count from 0 to 116.
      term  <- is.na(fl$flowline_toid) | fl$flowline_toid == 0
      exit  <- !term & !is.na(ds_fp) & ds_fp != fpid
      keep  <- !is.na(fpid)
      ex_sp <- exit & spine & keep
      list(fpid = fpid, ds_fp = ds_fp, spine = spine, keep = keep,
           n_exit_all = tapply(exit[keep], fpid[keep], sum),
           n_exit_sp = tapply(ex_sp[keep], fpid[keep], sum),
           lens      = if ("lengthkm" %in% names(fl))
             tapply(ifelse(spine[keep], as.numeric(fl$lengthkm)[keep], 0),
                    fpid[keep], sum) else NULL,
           der = tapply(seq_len(sum(keep)), fpid[keep], function(ii) {
             d <- unique(ds_fp[keep][ii][ex_sp[keep][ii]])
             if (length(d) == 1L) d else NA_character_
           }))
    }, error = function(e) NULL)

    if (!is.null(res)) {
      # -- flowpath_members_contiguous -------------------------------------
      # Counted over ALL members, not just the spine. A non-constituent
      # tributary inside the catchment drains to the mainstem WITHIN the
      # flowpath, so it is not an exit; if it leaves, the flowpath has a second
      # outlet, which is the defect. Measuring the spine alone misses exactly
      # the interleaved case (one spine exit, pointing at the wrong flowpath).
      n_disc <- sum(res$n_exit_all > 1L, na.rm = TRUE)
      n_tot  <- length(res$n_exit_all)
      ex_ids <- names(res$n_exit_all)[which(res$n_exit_all > 1L)]
      n_sp   <- sum(res$n_exit_sp > 1L, na.rm = TRUE)
      msg <- sprintf(paste0("%d of %d flowpath(s) have more than one member ",
        "reach leaving the flowpath (second outlet); %d of those on the ",
        "routing spine%s (target 0)"),
        n_disc, n_tot, n_sp,
        if (length(ex_ids)) sprintf(" (e.g. %s)",
          paste(utils::head(ex_ids, 3), collapse = ", ")) else "")
      checks$flowpath_members_contiguous <-
        if (n_disc == 0L) .hf_ok(TRUE, msg) else .hf_info(msg)

      # -- flowpath_topology_matches_flowlines -----------------------------
      if (!is.null(flowpaths) && !is.null(nexus) &&
          all(c("flowpath_id", "flowpath_toid") %in% names(flowpaths)) &&
          all(c("nexus_id", "nexus_toid") %in% names(nexus))) {
        n2f <- stats::setNames(as.character(nexus$nexus_toid),
                               as.character(nexus$nexus_id))
        ftoid    <- as.character(flowpaths$flowpath_toid)
        declared <- unname(n2f[ftoid])
        terminal <- !is.na(ftoid) & startsWith(ftoid, "tnx")
        derived  <- unname(res$der[as.character(flowpaths$flowpath_id)])
        cmp <- !terminal & !is.na(derived) & !is.na(declared)
        bad <- cmp & derived != declared
        n_bad <- sum(bad)
        bid <- as.character(flowpaths$flowpath_id)[which(bad)]
        msg <- sprintf(paste0("%d of %d non-terminal flowpath(s) declare a ",
          "downstream that the flowline network does not reach%s (target 0)"),
          n_bad, sum(cmp),
          if (length(bid)) sprintf(" (e.g. %s)",
            paste(utils::head(bid, 3), collapse = ", ")) else "")
        checks$flowpath_topology_matches_flowlines <-
          if (n_bad == 0L) .hf_ok(TRUE, msg) else .hf_info(msg)
      }

      # -- flowpath_length_matches_members ---------------------------------
      # lengthkm is the mainstem arm extend_flowpaths_to_mainstem walked, which
      # for a contiguous flowpath is exactly its constituent chain. A gap means
      # the drawn line is not the member chain.
      if (!is.null(res$lens) && !is.null(flowpaths) &&
          "lengthkm" %in% names(flowpaths)) {
        fp_len <- suppressWarnings(as.numeric(flowpaths$lengthkm))
        mem    <- unname(res$lens[as.character(flowpaths$flowpath_id)])
        rel    <- abs(fp_len - mem) / pmax(fp_len, 1e-9)
        bad    <- is.finite(rel) & fp_len > 0 & rel > 0.01
        n_bad  <- sum(bad, na.rm = TRUE)
        msg <- sprintf(paste0("%d of %d flowpath(s) have a lengthkm more than ",
          "1%% from the sum of their constituent member reaches (target 0)"),
          n_bad, sum(is.finite(rel) & fp_len > 0))
        checks$flowpath_length_matches_members <-
          if (n_bad == 0L) .hf_ok(TRUE, msg) else .hf_info(msg)
      }
    }
  }

  # -- split_siblings_share_catchment -----------------------------------------
  # Splitting a reference reach at a POI splits its catchment too: each fragment
  # drains the strip of land beside it. So when one fragment ends up with the
  # whole parent catchment and another with none, the partition did not split at
  # all, and the starved fragment ships as a real river owning no land.
  #
  # That is what a consumer sees as a flowpath with a NULL divide. The reach is
  # routed, carries drainage from upstream, and contributes no area of its own.
  # The loss happens later, when the mispaired divide is cleared as disjoint and
  # absorbed into a neighbour, but the starvation here is the cause and the only
  # point at which it can be stated without the reference in hand.
  #
  # Stated on split fragments alone, so it needs no threshold: siblings of one
  # reference reach either all carry catchment or the split failed. A reach whose
  # fragments all sit in one flowpath is not a split at product grain and is not
  # judged.
  if (!is.null(flowpaths) && all(c("member_comid", "areasqkm") %in% names(flowpaths))) {
    fpid <- as.character(flowpaths$flowpath_id)
    area <- suppressWarnings(as.numeric(flowpaths$areasqkm))
    mem  <- strsplit(as.character(flowpaths$member_comid), ",", fixed = TRUE)
    owner <- rep(fpid, lengths(mem))
    m     <- trimws(unlist(mem, use.names = FALSE))
    frag  <- grepl("[.][0-9]+$", m)
    n_bad <- 0L; n_grp <- 0L; ex <- character(0)
    if (any(frag)) {
      base <- sub("[.][0-9]+$", "", m[frag])
      own  <- owner[frag]
      a    <- stats::setNames(area, fpid)[own]
      by   <- split(seq_along(base), base)
      for (g in by) {
        fps <- unique(own[g])
        if (length(fps) < 2L) next          # not split across flowpaths
        n_grp <- n_grp + 1L
        ag <- a[g]
        af <- vapply(fps, function(f) unname(ag[own[g] == f][1L]), numeric(1L))
        if (any(af == 0, na.rm = TRUE) && any(af > 0, na.rm = TRUE)) {
          n_bad <- n_bad + 1L
          if (length(ex) < 3L) ex <- c(ex, fps[which(af == 0)[1L]])
        }
      }
    }
    msg <- sprintf(paste0("%d of %d split reference reach(es) have a fragment ",
      "with no catchment beside one that has it (target 0)%s"),
      n_bad, n_grp,
      if (length(ex)) sprintf(" (e.g. %s)", paste(ex, collapse = ", ")) else "")
    checks$split_siblings_share_catchment <-
      if (n_bad == 0L) .hf_ok(TRUE, msg) else .hf_info(msg)
  }

  # -- divides_no_enclosed_gaps ------------------------------------------------
  # Area of REAL dropped-catchment gaps: enclosed holes in the divide union
  # (holes fully surrounded by divides) that are backed by a Stage-1 rec_div,
  # i.e. a catchment that existed in the reference but whose polygon fell out of
  # aggregation. divides_area_conserved cannot see these (a gap can sit within its
  # tolerance yet still be a real hole).
  #
  # Raw enclosed-hole area is dominated by legitimate open-water / outside-domain
  # voids (large lakes, closed basins) that were never catchments, so gate each
  # hole by rec_div coverage (>= 0.5, matching .fill_enclosed_divide_gaps): only
  # holes that were actually catchments count. When `rec_div` is absent the total
  # is reported UNGATED with a clear label rather than a bare, misleading total.
  # INFO; gated behind `coverage` (full-layer terra dissolve) like other expensive
  # geometric checks.
  if (coverage && !is.null(divides) && inherits(divides, "sf") &&
      nrow(divides) > 0L) {
    rec_div <- args$rec_div
    res <- tryCatch({
      v  <- terra::aggregate(terra::vect(sf::st_make_valid(divides)))  # dissolve
      hz <- terra::erase(terra::fillHoles(v), v)                       # holes
      if (is.null(hz) || terra::nrow(hz) == 0L) list(total = 0, real = 0, gated = FALSE)
      else {
        hp <- terra::disagg(hz)
        a  <- terra::expanse(hp, unit = "m"); keep <- a >= 1
        hp <- hp[keep]; a <- a[keep]; total <- sum(a) / 1e6
        if (length(a) == 0L) list(total = 0, real = 0, gated = FALSE)
        else if (!is.null(rec_div) && inherits(rec_div, "sf") && nrow(rec_div) > 0L) {
          hs <- sf::st_as_sf(hp)
          rg <- sf::st_geometry(sf::st_make_valid(rec_div))
          if (sf::st_crs(rg) != sf::st_crs(hs)) rg <- sf::st_transform(rg, sf::st_crs(hs))
          hits <- sf::st_intersects(hs, rg); real_a <- 0
          for (i in seq_len(nrow(hs))) {
            ix <- hits[[i]]; if (!length(ix)) next
            inr <- suppressWarnings(sf::st_intersection(
              sf::st_geometry(hs)[i], sf::st_union(rg[ix])))
            if (!length(inr) || all(sf::st_is_empty(inr))) next
            if (as.numeric(sum(sf::st_area(inr))) / a[i] >= 0.5) real_a <- real_a + a[i] / 1e6
          }
          list(total = total, real = real_a, gated = TRUE)
        } else list(total = total, real = NA_real_, gated = FALSE)
      }
    }, error = function(e) NULL)
    checks$divides_no_enclosed_gaps <- if (is.null(res))
      .hf_info("divides_no_enclosed_gaps: not computable (geometry error)")
    else if (isTRUE(res$gated))
      .hf_info(sprintf(
        "%.3f km2 real dropped-catchment gap area (rec_div-gated; %.0f km2 total enclosed incl. open-water/outside-domain voids)",
        res$real, res$total))
    else
      .hf_info(sprintf(
        "%.3f km2 total enclosed gap area (UNGATED -- includes lakes/outside-domain voids; pass rec_div to isolate real dropped catchments)",
        res$total))
  }

  # -- flowpath_has_constituent ------------------------------------------------
  # Every flowpath carrying member flowlines must have >= 1 flowline flagged
  # flowpath_constituent == TRUE (the stem reach(es) its line traces).
  # flag_flowpath_constituents (build_nextgen Step 6m) guarantees this; a
  # violation means the flag pass was skipped or mis-computed for that flowpath.
  # Skipped silently when the attribute is absent (pre-Step-6m outputs).
  fl_c <- args$flowlines
  if (!is.null(fl_c) && inherits(fl_c, "sf") && nrow(fl_c) > 0L &&
      "flowpath_constituent" %in% names(fl_c) && "flowpath_id" %in% names(fl_c)) {
    by_fp <- tapply(as.logical(fl_c$flowpath_constituent),
      as.character(fl_c$flowpath_id), function(z) isTRUE(any(z, na.rm = TRUE)))
    n_empty <- sum(!by_fp)
    checks$flowpath_has_constituent <- .hf_ok(n_empty == 0L, sprintf(
      "%d flowpath(s) with member flowlines but no flowpath_constituent == TRUE",
      n_empty))
  }

  # -- divide_flowpath_exists -------------------------------------------------
  # Every divide's flowpath_id must name a flowpath that exists in the flowpaths
  # layer. A divide whose flowpath_id points to a missing flowpath is orphaned:
  # its divide_id integer no longer matches a real flowpath integer, which breaks
  # the cat-N/fp-N convention the global merge re-IDs by (yielding a duplicate
  # divide_id). Catch it here, per VPU, before the merge.
  if (!is.null(flowpaths) && !is.null(divides) &&
      "flowpath_id" %in% names(divides) && "flowpath_id" %in% names(flowpaths)) {
    dfp   <- as.character(divides$flowpath_id)
    has   <- !is.na(dfp) & nzchar(dfp)
    fpset <- as.character(flowpaths$flowpath_id)
    orphan <- has & !(dfp %in% fpset)
    checks$divide_flowpath_exists <- .hf_ok(sum(orphan) == 0L, sprintf(
      "%d divide(s) whose flowpath_id is not a flowpath in the flowpaths layer (e.g. %s)",
      sum(orphan),
      paste(utils::head(as.character(divides$divide_id)[orphan], 3), collapse = ", ")))
  }

  # -- divide_has_member_flowline ---------------------------------------------
  # The meaningful counterpart to flowpath_coverage, which is INFO because an
  # aggregated flowpath's mainstem line need not sit inside its own catchment.
  # What IS a defect: a network divide that contains none of its owning
  # flowpath's flowlines. Flowlines are the reconciled reaches and do carry
  # enforced geometry, so if not one of them touches the catchment they are
  # supposed to drain, the pairing is wrong rather than merely coarse.
  #
  # Judged within `tol_m` (30 m) rather than by strict intersection. 30 m is this
  # fabric's universal tolerance, the DEM resolution and the snap distance, and
  # a sub-30 m separation is treated as noise everywhere else it is measured, which
  # is how the multi-part flowpaths were judged. Strict intersection failed VPU 04
  # on one divide of 37302 (cat-62242, a 450 m2 coastal connector sliver) whose
  # flowline lies 0.84 m outside it: a rasterisation artefact three orders of
  # magnitude inside the tolerance, not a mispairing.
  #
  # This is a deliberate change to a gate's definition, so the reason is recorded
  # here. It weakens nothing meaningful: a genuine mispairing puts the flowline in a
  # different catchment entirely, hundreds of metres to kilometres away, and still
  # fails. The alternative considered was accepting a coverage tolerance so the
  # partition could be changed to fix the sliver, which was rejected: area
  # conservation is the fabric's foundational accounting property and loosening its
  # gate has far wider blast radius than reusing an established distance tolerance.
  .dhmf_tol <- 30
  .fl_mem <- args$flowlines
  if (!is.null(divides) && !is.null(.fl_mem) &&
      inherits(divides, "sf") && inherits(.fl_mem, "sf") &&
      "flowpath_id" %in% names(divides) && "flowpath_id" %in% names(.fl_mem)) {
    res <- tryCatch({
      net <- .hf_network_divides(divides)
      net <- net[!is.na(net$flowpath_id) & nzchar(as.character(net$flowpath_id)), ]
      if (nrow(net) == 0L) NULL else {
        fl_t <- sf::st_transform(.fl_mem, sf::st_crs(net))
        fl_by <- split(seq_len(nrow(fl_t)), as.character(fl_t$flowpath_id))
        owners <- as.character(net$flowpath_id)
        hits <- vapply(seq_len(nrow(net)), function(i) {
          ix <- fl_by[[owners[i]]]
          if (is.null(ix) || !length(ix)) return(NA)   # no flowlines: not testable
          if (any(lengths(sf::st_intersects(sf::st_geometry(fl_t)[ix],
            sf::st_geometry(net)[i])) > 0L)) return(TRUE)
          # Not intersecting: accept if within the fabric's 30 m tolerance.
          suppressWarnings(any(as.numeric(sf::st_distance(
            sf::st_geometry(fl_t)[ix], sf::st_geometry(net)[i])) <= .dhmf_tol,
            na.rm = TRUE))
        }, logical(1L))
        list(bad = which(!is.na(hits) & !hits), n = sum(!is.na(hits)),
             ids = as.character(net$divide_id))
      }
    }, error = function(e) NULL)
    if (!is.null(res))
      checks$divide_has_member_flowline <- .hf_ok(
        length(res$bad) == 0L,
        sprintf(paste0("%d of %d network divide(s) contain none of their ",
          "flowpath's flowlines within 30m (e.g. %s)"),
          length(res$bad), res$n,
          paste(utils::head(res$ids[res$bad], 3), collapse = ", ")))
  }

  # -- divide_id_names_a_flowpath ---------------------------------------------
  # The stage-4 precondition, which is NOT what divide_flowpath_exists tests.
  # Stage 4 re-IDs a cat- divide by looking up ITS OWN integer in that VPU's
  # flowpath pool, so cat-N must name a flowpath that exists. A divide can be
  # correctly owned and still fail this: cat-78392 owned by fp-78391 satisfies
  # divide_flowpath_exists, yet the merge could not map it, kept a VPU-local id,
  # and collided with another VPU's remapped divide. Assert it per VPU rather
  # than discovering it part-way through a global merge.
  if (!is.null(flowpaths) && !is.null(divides) &&
      "divide_id" %in% names(divides) && "flowpath_id" %in% names(flowpaths)) {
    did   <- as.character(divides$divide_id)
    iscat <- !is.na(did) & startsWith(did, "cat-")
    dint  <- suppressWarnings(as.integer(sub("^cat-", "", did)))
    fint  <- suppressWarnings(as.integer(sub("^fp-", "",
      as.character(flowpaths$flowpath_id))))
    bad   <- iscat & !is.na(dint) & !(dint %in% fint)
    checks$divide_id_names_a_flowpath <- .hf_ok(sum(bad) == 0L, sprintf(
      paste0("%d cat- divide(s) whose integer names no flowpath here, so the ",
        "stage-4 re-ID cannot map them (e.g. %s)"),
      sum(bad), paste(utils::head(did[bad], 3), collapse = ", ")))
  }

    if (!is.null(network)) {
      checks <- .hf_check_ngen_network(checks, network, args)
    }

  if (coverage && !is.null(flowpaths) && !is.null(divides))
    checks$flowpath_coverage <- .hf_coverage_check(flowpaths, divides, coverage_min)

  # -- flowpath_divide_assigned ---------------------------------------------
  # Every non-connector fp should have a divide assigned.  Connectors are fps
  # that legitimately have no catchment (nexus-only routing links); they are
  # identified by has_divide == FALSE *and* their toid is a nexus-id in the
  # nexus layer.  Any fp that is NOT a connector but still has divide_id=NA
  # indicates a divide was dropped by a cleanup step (mapshaper dedup, sliver
  # prune, etc.).  The count of connector fps equals
  # sum(!has_divide) in the invariant-passing state.
  if (!is.null(flowpaths) && !is.null(divides)) {
    fp_did      <- as.character(flowpaths$divide_id)
    no_div_mask <- is.na(fp_did) | !nzchar(fp_did)
    n_no_div    <- sum(no_div_mask)

    # True connectors: fps with zero drainage area (routing-only links that
    # legitimately have no catchment polygon).  A toid-prefix test (nex-/tnx-)
    # would match every ngen fp since all fps drain to a nexus, masking real
    # missing-divide gaps, so use drainage area instead.
    if ("areasqkm" %in% names(flowpaths)) {
      areas <- suppressWarnings(as.numeric(flowpaths$areasqkm))
      n_connectors <- sum(no_div_mask & !is.na(areas) & areas == 0)
    } else {
      n_connectors <- NA_integer_
    }

    n_unexplained <- if (is.na(n_connectors)) n_no_div
    else max(0L, n_no_div - n_connectors)

    checks$flowpath_divide_assigned <- .hf_ok(
      n_unexplained == 0L,
      sprintf(
        "%d fp(s) with no divide: %d connector(s) + %d unexplained gap(s)",
        n_no_div,
        if (is.na(n_connectors)) NA_integer_ else n_connectors,
        n_unexplained))
  }

  # -- POI contract, aggregated grain ----------------------------------------
  # This is the delivered layer, so these are the figures that count. All three
  # start as INFO because all three are non-zero today; each becomes a gate when
  # its class reaches zero across the domain.
  if (is.null(args$pois)) {
    # No POI layer is only an absence of information when there were no POIs to
    # carry. If hydrolocations went in and nothing came out, the layer was lost
    # somewhere in the build, and that is a defect rather than something left
    # unmeasured, an Alaska build once shipped this way and every POI check
    # stood down to a single info line while the product was unusable.
    n_in <- max(vapply(list(args$pois_input, args$hydrolocations),
      function(x) if (is.null(x)) 0L else nrow(x), integer(1L)))
    checks$pois_delivered <- .hf_ok(n_in == 0L, sprintf(
      "no POI layer in the product, from %d POI/hydrolocation row(s) upstream",
      n_in))
    checks$poi_at_outlet <- .hf_info(
      "no POI layer in the product -- POI contract not measured")
  } else {
    checks$pois_preserved <- .hf_poi_preserved(args$pois, args$pois_input,
      args$hydrolocations)
    checks$poi_ids_minted <- .hf_poi_ids_minted(args$pois, args$hydrolocations)
    checks$poi_ids_not_reminted <- .hf_poi_ids_not_reminted(args$pois,
      args$pois_input)
    checks$poi_id_is_integer <- .hf_poi_id_is_integer(args$pois,
      args$hydrolocations, args$network)

    if (!is.null(flowpaths)) {
      dec <- .hf_poi_decompose(args$pois, flowpaths, line_id_col = "flowpath_id")
      checks$poi_at_outlet <- .hf_poi_position_info(dec, "aggregated",
        types = if ("hl_types" %in% names(args$pois)) args$pois$hl_types
                else if ("hl_type" %in% names(args$pois)) args$pois$hl_type else NULL)
      checks$poi_at_flowline_outlet <- .hf_poi_at_flowline_outlet(args$pois,
        args$flowlines)
      checks$poi_location_stable <- .hf_poi_location_stable(args$pois,
        args$pois_input, args$flowlines)
      checks$poi_outlet_has_divide <- .hf_poi_outlet_has_divide(
        args$pois, flowpaths, divides, dec, line_id_col = "flowpath_id")
    }
  }

  checks$id_types_not_double <- .hf_id_types_not_double(args$gpkg)

  checks$connector_claim_valid <- .hf_connector_claim_valid(
    flowpaths, args$reference_gpkg)

  checks
}

# Identifiers must not be stored as floating point. A double id is exact only
# below 2^53 and is formatted in scientific notation by anything that round-trips
# it through text, and this fabric has already lost ids both ways, a 32-bit
# toid overflow and sci-notation COMID truncation. TEXT and INTEGER are both
# fine (prefixed ids like fp-/cat-/nex- are legitimately TEXT); REAL is not.
#
# Read from the GeoPackage's declared column types rather than from the loaded
# layers, because GDAL surfaces a correctly-typed INTEGER64 column to R as a
# double, so the R-side type cannot distinguish a stored id from a stored float.
.hf_id_types_not_double <- function(gpkg) {
  if (is.null(gpkg) || !is.character(gpkg) || !file.exists(gpkg)) return(NULL)
  id_cols <- c("flowpath_id", "flowpath_toid", "divide_id", "flowline_id",
    "flowline_toid", "nexus_id", "nexus_toid", "poi_id", "mainstem_id",
    "hl_id", "upstream_id", "hydroseq", "flowpath_hydroseq")
  tryCatch({
    con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    layers <- tryCatch(
      DBI::dbGetQuery(con, "SELECT table_name FROM gpkg_contents")$table_name,
      error = function(e) DBI::dbListTables(con))
    bad <- character(0)
    for (lyr in intersect(layers, DBI::dbListTables(con))) {
      ti <- DBI::dbGetQuery(con, sprintf('PRAGMA table_info("%s")', lyr))
      hit <- ti$name %in% id_cols &
        grepl("REAL|DOUB|FLOA", toupper(ti$type))
      if (any(hit)) bad <- c(bad, sprintf("%s.%s", lyr, ti$name[hit]))
    }
    .hf_ok(length(bad) == 0, sprintf(
      "%d id column(s) stored as floating point%s (target 0)",
      length(bad),
      if (length(bad)) sprintf(": %s", paste(utils::head(bad, 8), collapse = ", ")) else ""))
  }, error = function(e) NULL)
}


# A flowpath may be divide-less ONLY if the reference gives its reaches no real
# catchment. On gate vs INFO, see vignette("invariants").
#
# `flowpath_divide_assigned` accepts any has_divide == FALSE flowpath with a nexus
# toid as a connector, which cannot distinguish a reach that LOST its catchment
# from one that never had one. Area gates do not catch it either: conservation
# says nothing about who owns what.
#
# Measured as AREA AT RISK rather than the reference has_divide flag, because the
# flag is also set on reaches whose reference polygon is effectively empty,
# schema noise, not a lost catchment. The threshold is the 1000 m2 that
# prune_ghost_divides uses to decide a polygon is too small to be a catchment.
#
# Resolved through `member_comid`, the reference ids each flowpath aggregates,
# with any split suffix stripped. Reported unconditionally so an unavailable
# reference is distinguishable from a clean result.
.hf_connector_claim_valid <- function(flowpaths, reference_gpkg,
                                      min_area_sqkm = 0.001) {
  if (is.null(flowpaths) || !inherits(flowpaths, "sf")) return(NULL)
  if (!all(c("divide_id", "member_comid") %in% names(flowpaths))) return(NULL)
  if (is.null(reference_gpkg) || !is.character(reference_gpkg) ||
      !file.exists(reference_gpkg))
    return(.hf_info("connector claims not checked -- no reference available"))
  tryCatch({
    ref <- sf::st_read(reference_gpkg,
      query = "SELECT flowpath_id, has_divide, areasqkm FROM flowpaths",
      quiet = TRUE)
    if (!nrow(ref)) return(.hf_info("connector claims not checked -- reference empty"))
    hd <- stats::setNames(as.integer(ref$has_divide), as.character(ref$flowpath_id))
    ar <- stats::setNames(suppressWarnings(as.numeric(ref$areasqkm)),
                          as.character(ref$flowpath_id))
    nod <- which(is.na(flowpaths$divide_id) | !nzchar(as.character(flowpaths$divide_id)))
    if (!length(nod))
      return(.hf_ok(TRUE, "0 divide-less flowpath(s); no connector claim to check"))
    mem <- strsplit(as.character(flowpaths$member_comid)[nod], ",", fixed = TRUE)
    lost <- vapply(mem, function(m) {
      m <- unique(sub("[.][0-9]+$", "", trimws(m)))
      m <- m[!is.na(hd[m]) & hd[m] == 1L]
      if (!length(m)) 0 else sum(ar[m], na.rm = TRUE)
    }, numeric(1))
    bad   <- lost > min_area_sqkm
    n_bad <- sum(bad, na.rm = TRUE)
    # A `.N` suffix marks a piece of a reference reach that a POI cut. Resolving
    # it means stripping the suffix, so the area above is the WHOLE parent
    # catchment even where the sibling holds nearly all of it and only the split
    # is wrong. Measured on Alaska: of 3046.64 km2 reported this way, 2875.49
    # (94.4%) was in fact owned by the same parent's other children. The area is
    # therefore an upper bound and the COUNT is the trustworthy number; splitting
    # the two apart here keeps the figure from being read as lost ground.
    n_split <- sum(bad & vapply(mem, function(m)
      any(grepl("[.][0-9]+$", trimws(m))), logical(1)), na.rm = TRUE)
    ex <- if (n_bad) sprintf(" (e.g. %s)", paste(utils::head(
      as.character(flowpaths$flowpath_id)[nod][bad], 3), collapse = ", ")) else ""
    # INFO, not a gate, until the Alaska attribution is fixed. It fails every
    # domain today on a defect that predates this check, and a gate that always
    # fails is a gate people learn to ignore. Promote it once AK reaches zero.
    .hf_info(sprintf(paste0(
      "%d of %d divide-less flowpath(s) are NOT connectors -- %d are split ",
      "children, and the reference gives their reaches at most %.2f km2 of ",
      "catchment (upper bound: a split child is charged its whole parent, which ",
      "its sibling may already hold)%s (target 0)"),
      n_bad, length(nod), n_split, sum(lost[bad], na.rm = TRUE), ex))
  }, error = function(e) NULL)
}

# ---- helpers -----------------------------------------------------------------

# `.hf_network_is_dag` now lives in hfutils as `hf_network_is_dag` (hfutils >= 0.3.6).

.hf_ok   <- function(ok, msg) list(ok = isTRUE(ok), msg = msg, kind = "check")
.hf_info <- function(msg) list(ok = TRUE,       msg = msg, kind = "info")

# Network (catchment) divides only: drop landscape (lnd-) polygons, which are
# a separate coverage layer and not part of the conserved drainage area.
# Network (cat-) divides only. NOTE this excludes any divide that is lnd-prefixed
# OR has a NULL flowpath_id, which has a consequence worth knowing: releasing a
# divide "to unowned" removes its area from every accounting built on this
# helper, so the land silently leaves the conservation and coverage figures.
# Where area is involved, absorb a divide into a host rather than releasing it
# (see .absorb_orphan_divides()).
.hf_network_divides <- function(divides) {
  if (is.null(divides) || !inherits(divides, "sf") || nrow(divides) == 0L)
    return(divides)
  did <- as.character(divides$divide_id)
  is_lnd <- (!is.na(did) & startsWith(did, "lnd-"))
  if ("flowpath_id" %in% names(divides))
    is_lnd <- is_lnd | is.na(divides$flowpath_id)
  divides[!is_lnd, ]
}

# Divide-area integrity. Operates on network (cat-) divides only.
#   * divides_total_area       : INFO: count + total area (km2).
#   * divides_area_conserved   : CHECK (only with a baseline): output total
#                                 stays within area_tol_frac of the upstream
#                                 divide total. Guards "no lost divide area".
#   * flowpaths_no_shared_divide: CHECK (only with flowpaths): no divide_id is
#                                 claimed by more than one flowpath. Guards the
#                                 dendritic 1-catchment-per-flowpath rule.
#   * divides_non_overlapping  : CHECK (only when overlap = TRUE): summed
#                                 polygon area equals the union area, i.e. no
#                                 double-counted (overlapping) drainage area.
.hf_area_checks <- function(divides, flowpaths, area_baseline_m2,
                            area_tol_frac, overlap = FALSE,
                            coverage_baseline_m2 = NULL) {
  checks <- list()
  net <- .hf_network_divides(divides)
  if (is.null(net) || !inherits(net, "sf") || nrow(net) == 0L)
    return(checks)

  areas    <- suppressWarnings(as.numeric(sf::st_area(net)))
  total_m2 <- sum(areas, na.rm = TRUE)

  checks$divides_total_area <- .hf_info(
    sprintf("%d network divide(s), total area %.4f km2",
      nrow(net), total_m2 / 1e6))

  if (!is.null(area_baseline_m2) && is.finite(area_baseline_m2) &&
    area_baseline_m2 > 0) {
    rel <- abs(total_m2 - area_baseline_m2) / area_baseline_m2
    checks$divides_area_conserved <- .hf_ok(
      rel <= area_tol_frac,
      sprintf("divide area %.4f km2 vs baseline %.4f km2 (%.4f%% drift, tol %.3f%%)",
        total_m2 / 1e6, area_baseline_m2 / 1e6,
        100 * rel, 100 * area_tol_frac))
  }

  if (!is.null(flowpaths) && inherits(flowpaths, c("sf", "data.frame")) &&
    "divide_id" %in% names(flowpaths)) {
    fid <- as.character(flowpaths$divide_id)
    fid <- fid[!is.na(fid) & nzchar(fid)]
    dup <- unique(fid[duplicated(fid)])
    checks$flowpaths_no_shared_divide <- .hf_ok(
      length(dup) == 0L,
      sprintf("%d divide_id(s) claimed by >1 flowpath (e.g. %s)",
        length(dup), paste(utils::head(dup, 3), collapse = ", ")))
  }

  if (isTRUE(overlap)) {
    union_m2 <- tryCatch(
      sum(suppressWarnings(as.numeric(
        sf::st_area(sf::st_union(sf::st_geometry(net))))), na.rm = TRUE),
      error = function(e) NA_real_)
    if (is.finite(union_m2) && union_m2 > 0) {
      ov_m2   <- max(0, total_m2 - union_m2)
      ov_frac <- ov_m2 / union_m2
      checks$divides_non_overlapping <- .hf_ok(
        ov_frac <= area_tol_frac,
        sprintf("overlap area %.4f km2 (%.4f%% of union, tol %.3f%%)",
          ov_m2 / 1e6, 100 * ov_frac, 100 * area_tol_frac))

      # COVERAGE, as distinct from the summed area above. A sum cannot tell a
      # hole from an overlap removal: both reduce it. The union can, because only
      # a hole reduces the footprint. Reported against a union baseline so lost
      # ground is visible, a 97 km2 hole read as "ok" on the sum check because
      # it sat inside the 0.5% tolerance.
      checks$divides_total_coverage <- .hf_info(
        sprintf("divide coverage (union) %.4f km2", union_m2 / 1e6))
      if (!is.null(coverage_baseline_m2) && is.finite(coverage_baseline_m2) &&
        coverage_baseline_m2 > 0) {
        lost_m2 <- coverage_baseline_m2 - union_m2
        rel     <- abs(lost_m2) / coverage_baseline_m2
        checks$divides_coverage_conserved <- .hf_info(
          sprintf(paste0("coverage %.4f km2 vs baseline %.4f km2 ",
            "(%+.4f km2, %.4f%% drift)"),
            union_m2 / 1e6, coverage_baseline_m2 / 1e6,
            -lost_m2 / 1e6, 100 * rel))
      }
    } else {
      checks$divides_non_overlapping <- .hf_info(
        "non-overlap check skipped (union failed)")
    }
  }

  checks
}

.hf_ring_counts <- function(geoms) {
  vapply(geoms, function(g) {
    if (length(g) == 0L || sf::st_is_empty(g)) return(0L)
    if (inherits(g, "POLYGON")) return(max(0L, length(g) - 1L))
    if (inherits(g, "MULTIPOLYGON")) return(sum(vapply(g, function(p)
      max(0L, length(p) - 1L), integer(1L))))
    0L
  }, integer(1L))
}

.hf_coverage_check <- function(flowpaths, divides, coverage_min) {
  has_fp <- !is.na(flowpaths$divide_id) & nchar(as.character(flowpaths$divide_id)) > 0L
  fp_sub <- flowpaths[has_fp, ]
  if (nrow(fp_sub) == 0L) return(.hf_info("no flowpath-divide pairs to check"))

  fp_t  <- sf::st_transform(fp_sub, sf::st_crs(divides))
  d_ix  <- match(as.character(fp_t$divide_id), as.character(divides$divide_id))
  matched <- !is.na(d_ix)
  if (!any(matched))
    return(.hf_ok(FALSE, "no flowpaths matched a divide by divide_id"))

  fp_g      <- sf::st_geometry(fp_t)
  d_g       <- sf::st_geometry(divides)
  total_len <- as.numeric(sf::st_length(fp_g))
  inside    <- numeric(length(total_len))
  for (i in which(matched)) {
    inter <- tryCatch(
      sf::st_intersection(fp_g[[i]], d_g[[d_ix[i]]]),
      error = function(e) NULL)
    if (is.null(inter) || length(inter) == 0L) next
    inside[i] <- as.numeric(sf::st_length(sf::st_sfc(inter, crs = sf::st_crs(fp_t))))
  }
  frac <- ifelse(total_len > 0, inside / total_len, NA_real_)
  bad  <- which(matched & !is.na(frac) & frac < coverage_min)
  ok   <- length(bad) == 0L
  msg  <- if (ok)
    sprintf("all %d flowpaths have >= %.0f%% coverage by their catchment",
      sum(matched), 100 * coverage_min)
  else
    sprintf("%d flowpath(s) below %.0f%% coverage (worst: %s at %.0f%%)",
      length(bad), 100 * coverage_min,
      as.character(fp_sub$flowpath_id[bad[which.min(frac[bad])]]),
      100 * min(frac[bad], na.rm = TRUE))
  # INFO, not a check. An aggregated flowpath's line is a mainstem walk through
  # its members, so it is not required to lie inside its own catchment, the
  # same reasoning that made flowpath_outlet_at_nexus a measurement rather than
  # a gate. Worse, this only ever passed because Step 6j cleared the offenders
  # and Step 6k then deleted them, so a green verdict here meant "the evidence
  # was removed", not "the fabric is sound". The gate that does carry meaning is
  # divide_has_member_flowline below: a catchment containing none of its
  # flowpath's flowlines is a real defect.
  .hf_info(msg)
}

# What the suite was actually given.
#
# A check whose input is absent produces no entry at all, in R, assigning NULL
# to a list element removes it, so invariants.json cannot distinguish "passed",
# "does not apply at this stage", and "never ran because nothing was supplied".
# Of 66 checks defined for the ngen stage, 44 were recorded for VPU 01; the rest
# left no trace. Recording the inputs makes every absence interpretable without a
# catalog of checks to drift out of date.
.hf_input_summary <- function(args) {
  keep <- setdiff(names(args), "")
  if (!length(keep)) return(list())
  out <- lapply(keep, function(k) {
    v <- args[[k]]
    if (is.null(v)) NULL
    else if (is.data.frame(v)) nrow(v)
    else if (is.atomic(v) && length(v) == 1L) v
    else length(v)
  })
  stats::setNames(out, keep)
}

.hf_report_checks <- function(checks, stage, strict, inputs = NULL) {
  ok_all <- all(vapply(checks, function(x) isTRUE(x$ok), logical(1L)))
  for (nm in names(checks)) {
    res  <- checks[[nm]]
    icon <- if (identical(res$kind, "info")) "i"
    else if (isTRUE(res$ok)) "ok"
    else                              "x"
    tag  <- sprintf("[invariants:%s]", stage)
    line <- sprintf("%s %s %-32s %s", tag, icon, nm, res$msg)
    if (isTRUE(res$ok) || identical(res$kind, "info")) message(line)
    else if (strict) stop(line, call. = FALSE)
    else message(line)
  }
  invisible(list(ok = ok_all, stage = stage, checks = checks,
                 inputs = inputs))
}


# Flowline-grain checks, run only when the aspatial `network` table is supplied.
# Split out of .hf_check_ngen: these all read `network` and nothing else reads it,
# so they form a group that can be read and extended on its own.
.hf_check_ngen_network <- function(checks, network, args) {
  net_df <- if (inherits(network, "sf")) sf::st_drop_geometry(network)
  else as.data.frame(network)

  # -- flowline_id global uniqueness ----------------------------------------
  # After stage4 re-ID, each flowline_id should appear in exactly one VPU.
  # Collisions mean the per-VPU integer offsets were not applied to the
  # network table, causing downstream subsetting graph traversal to find the wrong
  # upstream set.
  if ("vpuid" %in% names(net_df)) {
    fl_vpu <- unique(net_df[!is.na(net_df$flowline_id),
      c("flowline_id", "vpuid"),
      drop = FALSE])
    n_dup <- sum(duplicated(fl_vpu$flowline_id))
    checks$flowline_id_globally_unique <- .hf_ok(
      n_dup == 0L,
      sprintf(
        "%d flowline_id collision(s) across VPUs -- re-run stage4 re-ID",
        n_dup))
  }

  # -- flowline_toid referential integrity ----------------------------------
  # Every non-zero, non-NA flowline_toid must reference a flowline_id that
  # exists in the table.  Dangling toids indicate a mismatch between the
  # flowlines and network layers.
  fl_pairs <- unique(net_df[!is.na(net_df$flowline_id) &
    !is.na(net_df$flowline_toid),
  c("flowline_id", "flowline_toid"),
  drop = FALSE])
  # Scientific-safe id formatting: as.character(1e5) is "1e+05", and a
  # double-typed flowline_toid then never matches an integer-typed
  # flowline_id (only bites once ids cross 1e5, e.g. large VPUs). Compare via
  # format(scientific=FALSE) on both sides.
  .fmt_id <- function(x) format(as.numeric(x), scientific = FALSE, trim = TRUE)
  id_set     <- unique(.fmt_id(fl_pairs$flowline_id))
  # Exclude terminal sentinels (0 = no downstream, >= 1e9 = cross-domain outlet)
  non_term   <- fl_pairs$flowline_toid[fl_pairs$flowline_toid != 0 &
    fl_pairs$flowline_toid < 1e9]
  n_dangling <- sum(!.fmt_id(non_term) %in% id_set)
  checks$flowline_toid_valid <- .hf_ok(
    n_dangling == 0L,
    sprintf("%d flowline_toid value(s) reference a non-existent flowline_id",
      n_dangling))

  # -- flowline directed graph is a DAG -------------------------------------
  is_dag_fl <- tryCatch(
    {
      edges <- fl_pairs[fl_pairs$flowline_toid != 0, , drop = FALSE]
      if (nrow(edges) == 0L) TRUE
      else {
        g <- igraph::graph_from_data_frame(
          data.frame(from = as.character(edges$flowline_id),
            to   = as.character(edges$flowline_toid),
            stringsAsFactors = FALSE),
          directed = TRUE)
        igraph::is_dag(g)
      }
    },
    error = function(e) NA)
  checks$flowline_graph_dag <- .hf_ok(
    isTRUE(is_dag_fl),
    if (isTRUE(is_dag_fl))
      sprintf("flowline graph is a DAG (%d edges)", nrow(fl_pairs))
    else "flowline graph contains cycle(s)")

  # -- flowline_toid geometry consistency -----------------------------------
  # A flowline's DOWNSTREAM (end) vertex must touch the flowline its
  # flowline_toid names: anywhere on that line, within 30m. This is the real
  # outlet-touches-inflow check: flowlines are the fine-grain, reliably-oriented
  # reaches, and a corrupt link (toid naming a geometrically distant reach) makes
  # extend_flowpaths_to_mainstem's walk veer.
  # LINE metric (not START): a tributary flowline joins its downstream mainstem
  # reach MID-reach at a confluence, so its outlet touches the downstream LINE
  # but not that reach's START vertex. Measuring end->downstream-START over-counts
  # every legitimate confluence; measuring end->nearest-point-on-downstream-LINE
  # (the "touch anywhere, not just end nodes" convention) reports only true gaps.
  if (inherits(args$flowlines, "sf") && nrow(args$flowlines) > 0L &&
      all(c("flowline_id", "flowline_toid") %in% names(args$flowlines))) {
    n_corrupt <- tryCatch({
      flg   <- args$flowlines
      gt    <- as.character(sf::st_geometry_type(flg))
      ok_ls <- gt == "LINESTRING"
      fid   <- as.character(flg$flowline_id)
      ti    <- stats::setNames(seq_len(nrow(flg)), fid)[as.character(flg$flowline_toid)]
      has   <- ok_ls & !is.na(ti) & ok_ls[ifelse(is.na(ti), 1L, ti)]
      if (!any(has)) 0L else {
        g  <- sf::st_geometry(flg)
        ep <- lwgeom::st_endpoint(g[has])       # each flowline's outlet vertex
        dn <- g[ti[has]]                         # its downstream flowline LINE
        d  <- suppressWarnings(as.numeric(sf::st_distance(ep, dn, by_element = TRUE)))
        sum(d > 30, na.rm = TRUE)
      }
    }, error = function(e) NA_integer_)
    if (!is.na(n_corrupt))
      checks$flowline_toid_geometry_consistent <- if (n_corrupt == 0L)
        .hf_ok(TRUE, "all flowline_toid links coincide with geometry (<=30m)")
      else
        # INFO, not a gate: these are geometry gaps inherited from the input
        # reference hydrography (the toid topology is correct; the reach end is
        # drawn >30m from its downstream reach start, a reference-build snap
        # gap present in both TDX and NHD sources). Reported per build for
        # visibility, not failed, since it is not a pipeline defect.
        .hf_info(sprintf(
          "%d flowline(s) whose flowline_toid names a reach >30m from the flowline's downstream end -- reference-hydrography geometry gap (topology correct)",
          n_corrupt))
  }

  # -- routeable_network_closed ----------------------------------------------
  # OWP routes on the `routeable` flowlines (the flowpath-constituent reaches
  # plus the incremental reaches that bridge each coincident reach to the next
  # coincident reach downstream; see flag_routeable_flowlines). That subgraph
  # must be closed under flowline_toid: every non-terminal routeable flowline
  # drains to a routeable flowline, so routing on the routeable set never
  # dead-ends in a filtered-out reach. CHECK: should be 0 by construction.
  #
  # Fallback for outputs that predate the `routeable` flag: report the raw
  # coincident-set gap as INFO (coincident flowline whose flowline_toid is a
  # non-coincident reach), which flag_routeable_flowlines is what closes.
  if (inherits(args$flowlines, "sf") && nrow(args$flowlines) > 0L &&
      all(c("flowline_id", "flowline_toid") %in% names(args$flowlines))) {
    fl2  <- sf::st_drop_geometry(args$flowlines)
    fid  <- as.character(fl2$flowline_id)
    ftid <- as.character(fl2$flowline_toid)
    term <- is.na(ftid) | ftid == "" | ftid == "0" | !(ftid %in% fid)
    if ("routeable" %in% names(fl2)) {
      rt   <- as.logical(fl2$routeable); rt[is.na(rt)] <- FALSE
      rt_of <- stats::setNames(rt, fid)
      dn_rt <- unname(rt_of[ftid])
      n_open <- tryCatch(sum(rt & !term & !is.na(dn_rt) & !dn_rt), error = function(e) NA_integer_)
      if (!is.na(n_open))
        checks$routeable_network_closed <- .hf_ok(n_open == 0L, sprintf(
          "%d routeable flowline(s) whose flowline_toid is not routeable (routing set not closed)",
          n_open))
    } else if ("flowpath_constituent" %in% names(fl2)) {
      coin <- as.logical(fl2$flowpath_constituent); coin[is.na(coin)] <- FALSE
      dn_coin <- unname(stats::setNames(coin, fid)[ftid])
      n_break <- tryCatch(sum(coin & !term & !is.na(dn_coin) & !dn_coin), error = function(e) NA_integer_)
      if (!is.na(n_break))
        checks$coincident_flowline_toid_closed <- .hf_info(sprintf(
          paste0("%d coincident flowline(s) whose flowline_toid is a NON-coincident reach ",
            "-- routing-set gap (no `routeable` flag; run flag_routeable_flowlines)"),
          n_break))
    }
  }

  # -- network covers all input reference_ids --------------------------------
  # Every reference COMID that entered from stage1 refactored should appear as a
  # reference_id in the network table.  Missing entries mean a flowpath was
  # dropped somewhere in the aggregation pipeline.
  if (!is.null(args$refactored) && "reference_id" %in% names(net_df)) {
    rfp_ids  <- unique(sub("\\.[0-9]+$", "",
      as.character(args$refactored$flowpath_id)))
    net_refs <- unique(as.character(net_df$reference_id[!is.na(net_df$reference_id)]))
    missing  <- setdiff(rfp_ids, net_refs)
    checks$network_covers_all_reference_ids <- .hf_ok(
      length(missing) == 0L,
      sprintf("%d/%d refactored reference_id(s) missing from network (e.g. %s)",
        length(missing), length(rfp_ids),
        paste(utils::head(missing, 3), collapse = ", ")))
  }

  # -- network covers all input flowlines ------------------------------------
  # Every flowline_id present in the ngen flowlines layer must appear in the
  # network table.  Gaps mean a flowline was generated but never wired into
  # the lineage trace.
  if (!is.null(args$flowlines) && "flowline_id" %in% names(net_df)) {
    fl_df   <- if (inherits(args$flowlines, "sf"))
      sf::st_drop_geometry(args$flowlines)
    else as.data.frame(args$flowlines)
    fl_ids  <- unique(as.character(fl_df$flowline_id[!is.na(fl_df$flowline_id)]))
    net_fls <- unique(as.character(net_df$flowline_id[!is.na(net_df$flowline_id)]))
    missing <- setdiff(fl_ids, net_fls)
    checks$network_covers_all_flowlines <- .hf_ok(
      length(missing) == 0L,
      sprintf("%d/%d flowline_id(s) missing from network (e.g. %s)",
        length(missing), length(fl_ids),
        paste(utils::head(missing, 3), collapse = ", ")))

    # -- flowline_id uniqueness within the flowlines layer -------------------
    # The flowlines layer is per-segment; flowline_id is the segment id and
    # must be unique. Duplicates mean it was set to the reconciled grouping
    # (reconciled_id) instead of the per-segment id.
    n_dup_fl <- sum(duplicated(fl_df$flowline_id[!is.na(fl_df$flowline_id)]))
    checks$flowlines_id_unique <- .hf_ok(
      n_dup_fl == 0L,
      sprintf(paste0("%d duplicate flowline_id(s) in the flowlines layer -- ",
        "should be unique per segment (reconciled grouping is reconciled_id)"),
        n_dup_fl))

    # -- flowline reference-attribution plausibility ------------------------
    # Guards against levelpath-grain crosswalk inflation, where every fragment
    # on a mainstem inherits the WHOLE levelpath's COMID set so the same
    # reference COMIDs are replicated across many sibling flowlines (members
    # inherited, not owned). Its signature is member SHARING, not raw count: a
    # densely-segmented reference legitimately tiles a short reconciled line
    # with many tiny PRIVATE COMIDs (each owned by exactly one flowline, their
    # lengths summing to the line). Flag a sub-km flowline only when it carries
    # many reference members that are mostly shared with other flowlines, the
    # inheritance footprint. Membership must be keyed by the reconciled line
    # (reconciled_id), not the levelpath.
    if ("reference_id" %in% names(net_df) && "lengthkm" %in% names(fl_df)) {
      rid  <- as.character(net_df$reference_id)
      flid <- as.character(net_df$flowline_id)
      keep <- !is.na(rid) & !is.na(flid)
      rid <- rid[keep]; flid <- flid[keep]
      pair <- !duplicated(paste(flid, rid))      # unique (flowline, member) edges
      flid_u <- flid[pair]; rid_u <- rid[pair]
      fl_per_ref <- tapply(flid_u, rid_u, function(x) length(unique(x)))
      shared_u   <- fl_per_ref[rid_u] > 1L       # member also claimed elsewhere?
      n_ref      <- tapply(rid_u, flid_u, length)
      frac_share <- tapply(shared_u, flid_u, mean)
      len_by_fl  <- stats::setNames(as.numeric(fl_df$lengthkm),
        as.character(fl_df$flowline_id))
      lens <- len_by_fl[names(n_ref)]
      bad  <- names(n_ref)[!is.na(lens) & lens < 1 & n_ref > 10L &
        frac_share[names(n_ref)] > 0.5]
      checks$flowline_reference_plausible <- .hf_ok(
        length(bad) == 0L,
        sprintf(paste0("%d sub-km flowline(s) carry >10 reference COMIDs that ",
          "are mostly shared across sibling flowlines -- levelpath-grain ",
          "crosswalk inflation (members inherited, not owned by reconciled_id; ",
          "e.g. flowline_id %s)"),
          length(bad), paste(utils::head(bad, 3), collapse = ", ")))
    }
  }
  checks
}

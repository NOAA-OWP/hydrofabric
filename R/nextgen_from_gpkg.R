#' Build a NextGen hydrofabric from a single GeoPackage
#'
#' @description
#' A self-contained function that reads all inputs--flowpaths, divides, network,
#' pois, flowlines, and events--from a single GeoPackage produced by the
#' `add_flowlines` step of the pipeline.  No separate `flowline_gpkg` is
#' required because `add_flowlines` writes full-resolution split segments (with
#' `id`, `toid`, `flowpath_id`, `areasqkm`, `fac_value`, etc.) directly into
#' the primary GeoPackage.
#'
#' @details
#' Topology steps:
#' 1. Build nexus geometries from flowpath endpoints.
#' 2. Identify terminal and headwater nexus nodes.
#' 3. Assign prefixed IDs to flowpaths, divides, and nexuses.
#' 4. Join the sub-flowline segments (from the `flowlines` layer) to their
#'    parent aggregated flowpaths via `flowlines$flowpath_id`.
#' 5. Compute `flowline_incr` (FAC-derived incremental drainage area in km2
#'    per sub-flowline: `total_dasqkm[i] - max(total_dasqkm[upstream donors])`).
#' 6. Build `incremental_areas` from divides joined to flowlines.
#' 7. Join hydrolocations from the `events` layer into the network table.
#'
#' @param flowpath_gpkg  Path to a GeoPackage containing at minimum:
#'   `flowpaths`, `divides`, `network`.  Should also contain `pois`,
#'   `flowlines` (from `add_flowlines`), and `events`.
#' @param collapse_method Character, `"upstream"` (default) or `"downstream"`.
#'   Controls nexus ID keying at confluences.  `"upstream"` keys the nexus on
#'   the **receiving** flowpath (`nex-X` where X = toid) so all donors to fp X
#'   share one nexus: fewer nexuses, shared at confluences.  `"downstream"`
#'   keys the nexus on the **donating** flowpath (`nex-X` where X = flowpath_id)
#'   so each donor gets its own outlet nexus, one nexus per fp.  Terminal
#'   nexuses are unaffected by this parameter.
#' @param nexus_prefix Character prefix for nexus IDs.  Default `"nex-"`.
#' @param terminal_nexus_prefix Prefix for terminal nexus IDs.  Default
#'   `"tnx-"`.
#' @param divide_prefix Prefix for divide IDs.  Default `"cat-"`.
#' @param flowpath_prefix Prefix for flowpath IDs.  Default `"fp-"`.
#' @param term_add Numeric offset added to terminal `flowpath_toid` values.
#'   Default `1e9`.
#' @param outfile Optional path for output GeoPackage.  If `NULL` (default)
#'   returns the list invisibly instead of writing.
#'
#' @return If `outfile` is specified, writes the result to that GeoPackage
#'   and returns `NULL` invisibly.  Otherwise
#'   returns a named list with elements `nexus`, `flowpaths`, `divides`,
#'   `flowlines`, `incremental_areas`, `pois` (when present), and `network`.
#'
#' @examples
#' \dontrun{
#' nextgen_from_gpkg("out_vpu01.gpkg", outfile = "ngen_vpu01.gpkg")
#' ng <- nextgen_from_gpkg("out_vpu01.gpkg")
#' }
#'
#' @export
#'
#' @importFrom sf st_read st_drop_geometry st_as_sf st_coordinates st_crs
#' @importFrom dplyr select mutate filter group_by ungroup distinct arrange
#'   left_join right_join semi_join bind_rows bind_cols any_of starts_with
#'   desc n cur_group_id
#' @importFrom hfutils node_geometry get_node write_hydrofabric add_lengthkm layer_exists

nextgen_from_gpkg <- function(flowpath_gpkg,
                              collapse_method          = "upstream",
                              nexus_prefix             = "nex-",
                              terminal_nexus_prefix    = "tnx-",
                              divide_prefix            = "cat-",
                              flowpath_prefix          = "fp-",
                              term_add                 = 1e9,
                              outfile                  = NULL) {

  collapse_method <- match.arg(collapse_method, c("upstream", "downstream"))

  cli::cli_alert("\n--- Applying NextGen topology ---\n")

  # ---- 1. Read flowpaths and divides by explicit layer name -----------------
  # (avoids read_hydrofabric auto-detecting the 'flowlines' layer)
  fp_raw  <- sf::st_read(flowpath_gpkg, layer = "flowpaths", quiet = TRUE)
  div_raw <- sf::st_read(flowpath_gpkg, layer = "divides",   quiet = TRUE)

  # Harmonise CRS
  if (!identical(sf::st_crs(fp_raw), sf::st_crs(div_raw))) {
    div_raw <- sf::st_transform(div_raw, sf::st_crs(fp_raw))
  }

  fp_network_list <- list(flowpaths = fp_raw, divides = div_raw)

  ngen_network_list <- list()

  # ---- 2. Build nexus topology from attributes --------------------------------
  # Terminal flowpaths carry flowpath_toid >= term_add (stable IDs provided
  # by aggregate_to_distribution).  All nexus IDs are derived
  # purely from the attribute table.

  fp_attr <- sf::st_drop_geometry(fp_network_list$flowpaths)
  # Coerce IDs to numeric for reliable comparisons: GeoPackage may store as character.
  fp_attr$flowpath_toid <- as.numeric(fp_attr$flowpath_toid)
  fp_attr$flowpath_id   <- as.numeric(fp_attr$flowpath_id)

  # Format a numeric ID as a plain integer string: prevents scientific notation
  # (e.g. 1e+05 instead of 100000) when IDs reach five or more digits.
  # sprintf("%.0f") formats each element independently, avoiding the vector-
  # level padding that formatC("d") applies when values span different widths.
  # Coerce to numeric first to handle character IDs stored as strings.
  .fmt_int <- function(x) sprintf("%.0f", as.numeric(x))

  # Non-terminal nexus geometry: place at the END of the mainstem donor
  # (lowest-hydroseq fp whose toid == this nexus id). That is the actual
  # confluence point where the donor delivers flow. The same approach is used
  # for cross-VPU nexuses where the receiving fp is in an adjacent VPU.
  non_term <- fp_attr[fp_attr$flowpath_toid != 0 &
    fp_attr$flowpath_toid < term_add, , drop = FALSE]

  fp_sf        <- fp_network_list$flowpaths
  fp_ids_num   <- as.numeric(fp_sf$flowpath_id)

  # Unique receiving flowpath IDs for non-terminal nexuses (raw numeric,
  # matching fp_ids_num which is pre-prefix at this stage).
  nterm_toids  <- unique(non_term$flowpath_toid)  # numeric receiving fp IDs
  recv_idx     <- match(nterm_toids, fp_ids_num)

  # Intra-VPU nexuses: receiving fp is in this dataset.
  # Place the nexus at the end of the mainstem donor (lowest hydroseq), which is
  # the actual confluence point. This matches the cross-VPU fallback approach.
  intra_mask   <- !is.na(recv_idx)
  intra_toids  <- nterm_toids[intra_mask]

  intra_donors <- non_term[non_term$flowpath_toid %in% intra_toids, , drop = FALSE]
  intra_repr   <- do.call(rbind, lapply(
    split(intra_donors, intra_donors$flowpath_toid), function(g) {
      g[which.min(g$hydroseq), , drop = FALSE]
    }))
  intra_repr_ordered <- intra_repr[match(intra_toids, intra_repr$flowpath_toid), , drop = FALSE]
  intra_sf        <- fp_sf[match(intra_repr_ordered$flowpath_id, fp_ids_num), ]
  recv_geom_intra <- hfutils::get_node(intra_sf, "end")

  # Guard the empty case: a basin that aggregates to a single terminal flowpath
  # has no intra-VPU non-terminal nexus, so intra_toids is length 0 while
  # get_node() on the empty donor still yields a length-1 geometry: building an
  # st_sf from those mismatched lengths errors. Return NULL when there are none.
  nex_sf_intra <- if (length(intra_toids) > 0L) sf::st_sf(
    data.frame(
      nexus_id   = if (collapse_method == "upstream")
        paste0(nexus_prefix, .fmt_int(intra_toids))
      else
        paste0(nexus_prefix, .fmt_int(intra_repr_ordered$flowpath_id)),
      nexus_toid = paste0(flowpath_prefix, .fmt_int(intra_toids)),
      stringsAsFactors = FALSE),
    geometry = sf::st_geometry(recv_geom_intra)) else NULL

  # Cross-VPU nexuses: receiving fp is in an adjacent VPU, not in this dataset.
  # Same approach: end of the lowest-hydroseq donor.
  cross_toids <- nterm_toids[!intra_mask]
  nex_sf_cross <- if (length(cross_toids) > 0L) {
    cross_donors <- non_term[non_term$flowpath_toid %in% cross_toids, , drop = FALSE]
    cross_repr   <- do.call(rbind, lapply(
      split(cross_donors, cross_donors$flowpath_toid), function(g) {
        g[which.min(g$hydroseq), , drop = FALSE]
      }))
    cross_sf        <- fp_sf[match(cross_repr$flowpath_id, fp_ids_num), ]
    recv_geom_cross <- hfutils::get_node(cross_sf, "end")
    sf::st_sf(
      data.frame(
        nexus_id   = if (collapse_method == "upstream")
          paste0(nexus_prefix, .fmt_int(cross_repr$flowpath_toid))
        else
          paste0(nexus_prefix, .fmt_int(cross_repr$flowpath_id)),
        nexus_toid = paste0(flowpath_prefix, .fmt_int(cross_repr$flowpath_toid)),
        stringsAsFactors = FALSE),
      geometry = sf::st_geometry(recv_geom_cross))
  } else NULL
  # Combine non-terminal nexus pieces, dropping any that are empty (NULL).
  nex_sf_nterm <- do.call(rbind, Filter(Negate(is.null),
    list(nex_sf_intra, nex_sf_cross)))

  # Terminal nexuses: no receiving fp exists (outlet to ocean / domain edge).
  # The aggregation step encodes terminal toids as 0; legacy code may encode as >= term_add.
  # When toid == 0, each terminal fp is its own outlet: key by flowpath_id.
  # When toid >= term_add, group by flowpath_toid (shared synthetic outlet id).
  term_fps  <- fp_attr[fp_attr$flowpath_toid == 0 | fp_attr$flowpath_toid >= term_add, , drop = FALSE]

  zero_toid <- term_fps[term_fps$flowpath_toid == 0, , drop = FALSE]
  high_toid <- term_fps[term_fps$flowpath_toid >= term_add, , drop = FALSE]

  make_term_sf <- function(fps, key_col) {
    if (nrow(fps) == 0L) return(NULL)
    repr <- do.call(rbind, lapply(split(fps, fps[[key_col]]), function(g)
      g[which.min(g$hydroseq), , drop = FALSE]))
    sf_sub  <- fp_sf[match(repr$flowpath_id, fp_ids_num), ]
    geom    <- hfutils::get_node(sf_sub, "end")
    sf::st_sf(
      data.frame(nexus_id   = paste0(terminal_nexus_prefix, .fmt_int(repr[[key_col]])),
        nexus_toid = NA_character_,
        stringsAsFactors = FALSE),
      geometry = sf::st_geometry(geom))
  }

  nex_sf_zero <- make_term_sf(zero_toid, "flowpath_id")
  nex_sf_high <- make_term_sf(high_toid, "flowpath_toid")
  nex_sf_term <- do.call(rbind, Filter(Negate(is.null), list(nex_sf_zero, nex_sf_high)))

  ngen_network_list$nexus <- do.call(rbind, Filter(Negate(is.null),
    list(nex_sf_nterm, nex_sf_term)))

  # ---- 3. (nexus built above) -----------------------------------------------

  # ---- 4. Prefix flowpath / divide IDs -------------------------------------
  ngen_network_list$flowpaths <- dplyr::mutate(
    fp_network_list$flowpaths,
    .orig_fp_id   = flowpath_id,
    flowpath_id   = ifelse(is.na(flowpath_id), NA,
      paste0(flowpath_prefix, .fmt_int(flowpath_id))),
    # Terminal fps get tnx- prefix; non-terminal fps get nex- prefix.
    # "upstream": nexus keyed on toid (receiver): nex-<toid>
    # "downstream": nexus keyed on own flowpath_id (donor): nex-<flowpath_id>
    flowpath_toid = dplyr::case_when(
      is.na(flowpath_toid)                                  ~ NA_character_,
      as.numeric(flowpath_toid) == 0                        ~ paste0(terminal_nexus_prefix, .fmt_int(.orig_fp_id)),
      as.numeric(flowpath_toid) >= term_add                 ~ paste0(terminal_nexus_prefix, .fmt_int(flowpath_toid)),
      collapse_method == "upstream"                         ~ paste0(nexus_prefix, .fmt_int(flowpath_toid)),
      TRUE                                                  ~ paste0(nexus_prefix, .fmt_int(.orig_fp_id))),
    .orig_fp_id   = NULL,
    divide_id     = ifelse(is.na(divide_id), NA,
      paste0(divide_prefix, .fmt_int(divide_id))))

  ngen_network_list$divides <- dplyr::left_join(
    dplyr::mutate(fp_network_list$divides,
      divide_id     = ifelse(is.na(divide_id), NA,
        paste0(divide_prefix, .fmt_int(divide_id))),
      flowpath_toid = NULL),
    dplyr::select(sf::st_drop_geometry(ngen_network_list$flowpaths),
      divide_id, flowpath_id, flowpath_toid),
    by = "divide_id")

  # ---- 5. Join sub-flowline segments ----------------------------------------
  has_flowlines <- hfutils::layer_exists(flowpath_gpkg, "flowlines")
  if (has_flowlines) {
    cli::cli_alert_success("\n--- Flowline topology provided ---\n")

    fl_raw <- sf::st_read(flowpath_gpkg, layer = "flowlines", quiet = TRUE)

    # Prefix the parent flowpath IDs on the segments to match the prefixed
    # ngen_network_list$flowpaths, then join on flowpath_id.
    fl_raw <- dplyr::mutate(fl_raw,
      flowpath_id   = ifelse(is.na(flowpath_id), NA,
        paste0(flowpath_prefix, .fmt_int(flowpath_id))))
    # flowpath_toid is a flowpath-level attribute: every sub-flowline shares its
    # parent flowpath's downstream nexus. Take it straight from the already
    # terminal-aware prefixed flowpaths layer (toid 0 -> tnx-, >=term_add -> tnx-,
    # else nex-). Re-prefixing the raw numeric here instead would emit "nex-0" for
    # terminal flowpaths and disagree with the flowpaths/nexus layers, and with
    # the network table, which is built from these flowlines.
    fp_toid_lkp <- stats::setNames(ngen_network_list$flowpaths$flowpath_toid,
      ngen_network_list$flowpaths$flowpath_id)
    fl_raw$flowpath_toid <- unname(fp_toid_lkp[fl_raw$flowpath_id])

    fl_joined <- fl_raw

    # flowline_incr: FAC-derived incremental contributing area in km2. This is a
    # fallback and does not conserve basin area: the sub-flowline donor graph is
    # often incomplete, so most segments end up with incremental == total_dasqkm.
    # Step 6d8 supersedes it with a rec_div-anchored, length-weighted value that
    # conserves by construction, and this value survives only when rec_div is
    # unavailable.
    #
    # Over the cross-flowpath toid graph, donor_max[j] = max(total drainage of all
    # i with toid[i] == id[j]). Head sub-flowlines resolve correctly because their
    # donors are the outlet sub-flowlines of the upstream flowpath. A result above
    # 110% of the parent flowpath's areasqkm is set NA, the signature of a FAC gap
    # at a domain boundary.
    local({
      da    <- fl_joined$tot_drainage_areasqkm
      id    <- fl_joined$id
      toid  <- fl_joined$toid
      fp_id <- fl_joined$flowpath_id

      valid <- !is.na(toid) & !is.na(da) & toid != 0L
      donor_max <- rep(NA_real_, length(id))
      if (any(valid)) {
        donors <- data.table::data.table(
          recipient = toid[valid],
          da        = da[valid]
        )[, .(donor_max = max(da, na.rm = TRUE)), by = recipient]
        donor_max <- donors$donor_max[match(id, donors$recipient)]
      }

      raw_incr <- pmax(0, da - ifelse(is.na(donor_max), 0, donor_max),
        na.rm = FALSE)

      # Backstop: outside-CONUS or domain-edge FAC cells yield near-zero
      # tot_drainage_areasqkm, producing inflated incremental values when a
      # large-drainage donor is subtracted. Cap at 110% of parent fp areasqkm.
      fp_attr <- sf::st_drop_geometry(fp_network_list$flowpaths)
      fp_area <- stats::setNames(fp_attr$areasqkm,
        as.character(fp_attr$flowpath_id))[as.character(fp_id)]
      bad <- !is.na(raw_incr) & !is.na(fp_area) &
        raw_incr > pmax(fp_area * 1.1, 0.01, na.rm = TRUE)
      raw_incr[bad] <- NA_real_
      fl_joined$flowline_incr <<- raw_incr
    })

    fl_joined$incremental_area_id <- fl_joined$id

    # The flowlines layer is at the FDR/FAC segment grain, one distinct,
    # contiguous segment per row, so flowline_id IS the (unique) segment id.
    # The reconciled-line id is carried as reconciled_id for the
    # COMID -> reconciled -> flowpath lineage. Using rec_id as flowline_id would
    # duplicate it across every segment of a reconciled line and make the
    # flowpaths<->flowlines join many-to-many, inflating the network.
    # id / toid are seg_id-space; rec_id / rec_toid the reconciled grouping.
    fl_joined$flowline_id     <- fl_joined$id
    fl_joined$flowline_toid   <- if ("toid"     %in% names(fl_joined)) fl_joined$toid     else NA_integer_
    fl_joined$reconciled_id   <- if ("rec_id"   %in% names(fl_joined)) fl_joined$rec_id   else NA_integer_
    fl_joined$reconciled_toid <- if ("rec_toid" %in% names(fl_joined)) fl_joined$rec_toid else NA_integer_

    ngen_network_list$flowlines <- dplyr::select(
      fl_joined,
      flowline_id,
      flowline_toid,
      reconciled_id,
      reconciled_toid,
      flowline_incr,
      flowline_hydroseq     = dplyr::any_of("hydroseq"),
      flowpath_id,
      flowpath_toid,
      incremental_area_id,
      dplyr::any_of(c("fac_value", "fac_source", "tot_drainage_areasqkm",
        "lengthkm", "levelpathid", "orig_levelpathID")))

    # incremental_areas layer removed: without true per-sub FDR watershed
    # delineation every sub-flowline would inherit the same parent catchment
    # polygon, exactly duplicating the divides layer. Callers that need
    # per-sub incremental areas should join flowlines -> divides by
    # flowpath_id directly.
  }

  # ---- 6. POIs --------------------------------------------------------------
  if (hfutils::layer_exists(flowpath_gpkg, "pois")) {
    ngen_network_list$pois <- dplyr::mutate(
      dplyr::select(sf::st_read(flowpath_gpkg, "pois", quiet = TRUE),
        poi_id, flowpath_id),
      flowpath_id = ifelse(is.na(flowpath_id), NA,
        paste0(flowpath_prefix, .fmt_int(flowpath_id))))
  }

  # ---- 7. Network table -----------------------------------------------------
  net <- dplyr::mutate(
    dplyr::select(sf::st_read(flowpath_gpkg, "network", quiet = TRUE),
      flowpath_id, reference_id = hf_id),
    flowpath_id = ifelse(is.na(flowpath_id), NA,
      paste0(flowpath_prefix, .fmt_int(flowpath_id))))

  net <- dplyr::right_join(
    dplyr::mutate(
      .separate_longer(
        sf::st_drop_geometry(
          dplyr::select(ngen_network_list$flowpaths,
            flowpath_id, flowpath_toid, hydroseq, divide_id,
            reconciled_id = rec_id_list, poi_id,
            flowpath_lengthkm = lengthkm,
            areasqkm, tot_drainage_areasqkm)),
        col   = "reconciled_id",
        delim = ","),
      reconciled_id = as.integer(reconciled_id)),
    net,
    by           = "flowpath_id",
    relationship = "many-to-many",
    na_matches   = "never")

  if (has_flowlines) {
    fl_net <- dplyr::mutate(
      ngen_network_list$flowlines,
      flowline_lengthkm = hfutils::add_lengthkm(ngen_network_list$flowlines))

    net <- dplyr::right_join(
      sf::st_drop_geometry(
        dplyr::select(fl_net,
          flowline_id, flowline_toid, reconciled_id, flowline_hydroseq,
          flowline_lengthkm, incremental_area_id)),
      net,
      by           = "reconciled_id",
      relationship = "many-to-many",
      na_matches   = "never")

    # Join hydrolocations from events layer
    if (hfutils::layer_exists(flowpath_gpkg, "events")) {
      hl <- dplyr::distinct(
        dplyr::select(sf::st_read(flowpath_gpkg, "events", quiet = TRUE),
          dplyr::any_of(c("poi_id", "hl_reference",
            "hl_class", "hl_source", "hl_uri"))))
      hl$poi_id <- as.character(hl$poi_id)
      net <- dplyr::right_join(hl, net, by = "poi_id",
        relationship = "many-to-many", na_matches = "never")
    }
  }

  ngen_network_list$network <- dplyr::distinct(
    dplyr::select(net,
      dplyr::starts_with("flowpath"),
      "hydroseq", "divide_id",
      dplyr::starts_with("flowline"),
      "incremental_area_id",
      "poi_id",
      dplyr::starts_with("hl"),
      "reference_id",
      dplyr::everything()))

  # ---- 8. Write or return --------------------------------------------------
  if (!hfutils::hf_network_is_dag(ngen_network_list$flowpaths,
    "flowpath_id", "flowpath_toid")) {
    stop("nextgen_from_gpkg: fp->nex->fp graph is not a DAG -- fix upstream in refactor/reconcile")
  }

  if (!is.null(outfile)) {
    hfutils::write_hydrofabric(ngen_network_list, outfile = outfile,
      enforce_dm = FALSE)
    invisible(NULL)
  } else {
    invisible(ngen_network_list)
  }
}

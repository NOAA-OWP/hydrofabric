#' Annotate Flowpaths and Divides with Type Info
#'
#' Adds boolean flags `has_divide` to flowpaths and `has_flowpath` to divides.
#' @param network_list List of `flowpaths` and `divides` (sf objects).
#' @return Updated `network_list` with type info columns.
#' @examples
#' \dontrun{
#' nl <- list(flowpaths = sf::st_read("out_vpu01.gpkg", "flowpaths"),
#'            divides   = sf::st_read("out_vpu01.gpkg", "divides"))
#' nl <- add_network_type(nl)
#' }
#' @importFrom dplyr mutate filter
#' @export

add_network_type <- function(network_list) {

  network_list$flowpaths <- network_list$flowpaths |>
    dplyr::mutate(has_divide = flowpath_id %in% network_list$divides$divide_id) |>
    dplyr::filter(!duplicated(flowpath_id))

  network_list$divides <- network_list$divides |>
    dplyr::mutate(has_flowpath = divide_id %in% network_list$flowpaths$flowpath_id) |>
    dplyr::filter(!duplicated(divide_id))

  return(network_list)
}

#' Prepare Hydrologic Network
#'
#' @description
#' Adds hydrosequence (`hydroseq`), geometric measures (area/length),
#' and total contributing drainage area (`tot_drainage_area`) to the
#' `flowpaths` element of a network list. Also standardizes geometry column
#' names, normalizes names to lowercase, and drops duplicates.
#'
#' @details
#' - `tot_drainage_area` is computed by accumulating per-flowpath `areasqkm`
#'   downstream; `NA` areas are treated as zero to enable accumulation.
#' - The function expects a list with `flowpaths` and `divides` `sf` objects.
#' - If duplicate rows or duplicate IDs are found, they are removed with a
#'   warning.
#'
#' @param network_list A list containing `sf` objects:
#'   \itemize{
#'     \item `flowpaths`: flowpath/flowline features with at least
#'       `flowpath_id`, `flowpath_toid`.
#'     \item `divides`: catchment/divide features (ideally with `divide_id`).
#'   }
#'
#' @return The input list with updated `flowpaths` and `divides` `sf` objects.
#'
#' @examples
#' \dontrun{
#' nl <- list(flowpaths = sf::st_read("rec_fp_agg.gpkg", "flowpaths"),
#'            divides   = sf::st_read("rec_div.gpkg", "divides"))
#' nl <- prepare_network(nl)
#' }
#'
#' @export
#'
#' @importFrom sf st_as_sf st_drop_geometry st_geometry<- st_geometry_type st_crs
#' @importFrom dplyr select mutate filter right_join distinct any_of across
#' @importFrom cli cli_alert_info cli_alert_warning cli_alert_success cli_abort
#' @importFrom hfutils get_hydroseq accumulate_downstream read_hydrofabric add_measures rename_geometry

prepare_network <- function(network_list) {
  # ---- validate --------------------------------------------------------------
  if (!is.list(network_list) ||
    !all(c("flowpaths", "divides") %in% names(network_list))) {
    cli::cli_abort("`network_list` must be a list with elements `flowpaths` and `divides`.")
  }
  if (!inherits(network_list$flowpaths, "sf")) {
    cli::cli_abort("`network_list$flowpaths` must be an {.`sf`} object.")
  }
  if (!inherits(network_list$divides, "sf")) {
    cli::cli_abort("`network_list$divides` must be an {.`sf`} object.")
  }

  # ---- helpers ---------------------------------------------------------------
  # .rename_geometry <- function(x, nm = "geom") {
  #   # rename the active geometry column to `nm` (no-op if already named)
  #   gcol <- attr(x, "sf_column")
  #   if (!identical(gcol, nm)) {
  #     sf::st_geometry(x) <- sf::st_geometry(x)
  #     attr(x, "sf_column") <- nm
  #     names(x)[names(x) == gcol] <- nm
  #   }
  #   x
  # }

  .require_cols <- function(x, cols, where) {
    miss <- setdiff(cols, names(x))
    if (length(miss)) {
      cli::cli_abort(c(
        "Missing {qty(length(miss))} required column{?s} in `{where}`:",
        "i" = "{.val {miss}}"
      ))
    }
  }

  .dedup_by <- function(x, id_col = NULL, what = "rows") {
    if (!is.null(id_col) && id_col %in% names(x)) {
      dup_idx <- duplicated(x[[id_col]])
      n <- sum(dup_idx, na.rm = TRUE)
      if (n) {
        ids <- unique(x[[id_col]][dup_idx])
        cli::cli_alert_warning(
          "Dropping {n} duplicate {what} by `{id_col}`. Example IDs: {paste(utils::head(ids, 5), collapse = ', ')}"
        )
        x <- dplyr::filter(x, !duplicated(.data[[id_col]]))
      }
      return(x)
    }
    # Fallback: whole-row duplicate removal (ignoring geometry to be robust)
    non_geom <- unique(c(setdiff(names(x), attr(x, "sf_column"))))
    dup_idx <- duplicated(sf::st_drop_geometry(x)[, non_geom, drop = FALSE])
    n <- sum(dup_idx, na.rm = TRUE)
    if (n) {
      cli::cli_alert_warning("Dropping {n} duplicate {what} (row-wise, ignoring geometry).")
      x <- x[!dup_idx, ]
    }
    x
  }

  # ---- standardize names & geometry -----------------------------------------
  network_list$flowpaths <- hfutils::rename_geometry(network_list$flowpaths, "geom")
  network_list$divides   <- hfutils::rename_geometry(network_list$divides,   "geom")

  names(network_list$flowpaths) <- tolower(names(network_list$flowpaths))
  names(network_list$divides)   <- tolower(names(network_list$divides))

  # Accept ID/toID as aliases for flowpath_id/flowpath_toid
  fp_nms <- names(network_list$flowpaths)
  if (!"flowpath_id" %in% fp_nms && "id" %in% fp_nms)
    names(network_list$flowpaths)[fp_nms == "id"] <- "flowpath_id"
  if (!"flowpath_toid" %in% fp_nms && "toid" %in% fp_nms)
    names(network_list$flowpaths)[names(network_list$flowpaths) == "toid"] <- "flowpath_toid"

  network_list$flowpaths <- .hf_as_character_cols(
    network_list$flowpaths,
    cols = c("flowpath_id", "flowpath_toid", "levelpathid")
  )
  network_list$divides <- .hf_as_character_cols(
    network_list$divides,
    cols = c("divide_id", "flowpath_id")
  )

  # ---- drop duplicates -------------------------------------------------------
  network_list$divides   <- .dedup_by(network_list$divides,   id_col = "divide_id",   what = "catchments")
  network_list$flowpaths <- .dedup_by(network_list$flowpaths, id_col = "flowpath_id", what = "flowpaths")

  # Additional duplicate protection on flowpath_id, if present
  if ("flowpath_id" %in% names(network_list$flowpaths)) {
    n <- sum(duplicated(network_list$flowpaths$flowpath_id), na.rm = TRUE)
    if (n) {
      cli::cli_alert_warning("Dropping {n} duplicate flowpaths by `flowpath_id`.")
      network_list$flowpaths <- dplyr::filter(
        network_list$flowpaths,
        !duplicated(.data$flowpath_id)
      )
    }
  }

  # ---- detect and break cycles (IBTs, divergences) --------------------------
  dag_chk <- hfutils::hf_assert_network_dag(
    network_list$flowpaths, id_col = "flowpath_id", toid_col = "flowpath_toid")
  if (!isTRUE(dag_chk$is_dag)) {
    cli::cli_alert_warning(paste0("prepare_network: ", dag_chk$message,
      "; breaking cycles by severing one edge per SCC"))
    network_list$flowpaths <- hfutils::hf_break_cycles(
      network_list$flowpaths, id_col = "flowpath_id", toid_col = "flowpath_toid")
  }

  # ---- add hydrosequence (topological order) --------------------------------
  .require_cols(x = network_list$flowpaths, cols = c("flowpath_id", "flowpath_toid"), where = "flowpaths")

  # Remove any stale topo columns, compute topo sort, and join back
  network_list$flowpaths <- network_list$flowpaths |>
    dplyr::select(-dplyr::any_of(c("topo_sort", "hydroseq"))) |>
    (\(fp) {
      hs_try <- tryCatch(hfutils::get_hydroseq(fp), warning = function(w) w, error = function(e) e)
      hs <- NULL

      if (inherits(hs_try, "warning") || inherits(hs_try, "error")) {
        cli::cli_alert_warning("hfutils::get_hydroseq fallback triggered: using internal safe hydroseq")
      } else {
        hs <- suppressWarnings(as.numeric(hs_try))
      }

      bad <- is.null(hs) || length(hs) != nrow(fp) || anyNA(hs) || any(!is.finite(hs))
      if (isTRUE(bad)) {
        hs <- hfutils::hf_recompute_hydroseq(fp, id_col = "flowpath_id", toid_col = "flowpath_toid")$hydroseq
      }

      dplyr::mutate(fp, hydroseq = hs)
    })()

  # ---- add measures (area/length), normalize area NAs to 0 ------------------
  network_list$flowpaths$divide_id <- NULL

  # Aggregation rebuilds topology several times with geometry-free divide
  # stubs. In that case add_measures() cannot recover area from geometry and
  # would replace a valid, already-carried flowpath area with zero. Preserve
  # the carried area for topology-only passes; use the divide geometry on the
  # initial/prepared pass.
  has_divide_geometry <- any(!sf::st_is_empty(sf::st_geometry(network_list$divides)))
  carried_area <- network_list$flowpaths$areasqkm
  if (has_divide_geometry) {
    network_list <- hfutils::add_measures(
      flowpaths = network_list$flowpaths,
      divides   = network_list$divides
    )
  } else {
    network_list$flowpaths$lengthkm <- hfutils::add_lengthkm(network_list$flowpaths)
    if (!is.null(carried_area)) network_list$flowpaths$areasqkm <- carried_area
  }

  # Ensure areasqkm exists for accumulation; treat NA as 0
  if (!"areasqkm" %in% names(network_list$flowpaths)) {
    cli::cli_alert_info("`areasqkm` not found on flowpaths; creating as 0 for accumulation.")
    network_list$flowpaths$areasqkm <- 0
  } else {
    network_list$flowpaths <- dplyr::mutate(
      network_list$flowpaths,
      areasqkm = ifelse(is.na(.data$areasqkm), 0, .data$areasqkm)
    )
  }

  # ---- total contributing drainage area -------------------------------------
  network_list$flowpaths$tot_drainage_area <- tryCatch(
    hfutils::accumulate_downstream(network_list$flowpaths, attr = "areasqkm"),
    error = function(e) {
      if (!grepl("cycle", conditionMessage(e), ignore.case = TRUE)) stop(e)
      cli::cli_alert_warning(paste0(
        "prepare_network: accumulate_downstream detected cycles; ",
        "breaking cycles and retrying"))
      network_list$flowpaths <<- hfutils::hf_break_cycles(
        network_list$flowpaths, id_col = "flowpath_id", toid_col = "flowpath_toid")
      hfutils::accumulate_downstream(network_list$flowpaths, attr = "areasqkm")
    }
  )

  # ---- final validation ------------------------------------------------------
  check_network_validity(network_list)

  cli::cli_alert_success("Hydrologic network prepared: added `hydroseq`, measures, and `tot_drainage_area`.")

  network_list
}

#' Check Network Validity
#' Validates a flowpath and catchment network
#' @param network_list list containing flowline and catchment `sf` objects
#' @param term_cut cutoff integer to define terminal IDs
#' @param check logical; if FALSE, skips checks and returns input
#' @return a list containing flowline and catchment `sf` objects
#' @noRd
#' @importFrom dplyr mutate select left_join
#' @importFrom sf st_drop_geometry

check_network_validity <- function(network_list,
                                   term_cut = 1e9,
                                   check = TRUE) {

  flowpaths <- network_list$flowpaths
  div <- network_list$divides

  names(flowpaths) <- tolower(names(flowpaths))
  names(div) <- tolower(names(div))

  flowpaths$flowpath_toid <- ifelse(is.na(flowpaths$flowpath_toid), 0, flowpaths$flowpath_toid)

  if (!"ibt" %in% names(flowpaths)) {
    flowpaths$ibt <- FALSE
  }

  if (!check) {
    return(list(flowpaths = flowpaths, divides = div))
  }

  DAG <- hfutils::hf_network_is_dag(flowpaths, id_col = "flowpath_id", toid_col = "flowpath_toid")
  toid_chr <- as.character(flowpaths$flowpath_toid)
  toid_num <- suppressWarnings(as.numeric(toid_chr))
  CONNECTION <- all(
    toid_chr %in% as.character(flowpaths$flowpath_id) |
      (!is.na(toid_num) & toid_num > term_cut) |
      toid_chr == "0" |
      toid_chr == 0 |
      flowpaths$ibt
  )

  if (all(DAG, CONNECTION)) {
    return(list(flowpaths = flowpaths, divides = div))
  } else {
    if (!DAG) {
      stop("Network is not a graph.")
    }
    if (!CONNECTION) {
      stop("All toIDs are not present in network.")
    }
  }
}

# Validating a hydrolocation set before a build reads it.
#
# The reference gate catches damaged topology and the grid gate catches a grid
# meant for elsewhere; neither looks at the hydrolocations, which is where this
# delivery's input defects have actually lived. Four kinds turned up in one set:
# a layer named something the pipeline does not read, identifiers indexed against
# a different network, vpuid values naming units no build runs, and rows
# duplicated by a merge. Each is cheap to detect here and expensive to discover
# an hour into a build, or after it.

#' Validate a hydrolocation set before building from it
#'
#' @description
#' Checks a hydrolocation GeoPackage against what the pipeline requires of it.
#' Gated checks are properties a build depends on; informational ones report
#' shape without judging it.
#'
#' Gated: the layer is named so the pipeline can read it; geometry is POINT and
#' non-empty; the CRS is defined; `hl_reference` and `flowpath_id` are populated;
#' no row is an exact duplicate of another; every `vpuid` names a unit some
#' supplied reference carries; and every `flowpath_id` resolves to a reach in its
#' domain's reference.
#'
#' Reported: how many rows share an `hl_reference` (collocation is by design,
#' since several hydrolocations may collapse onto one POI), how many arrive with
#' a `poi_id`, and the mix of `hl_type`.
#'
#' @param hl_gpkg Path to the hydrolocation GeoPackage.
#' @param layer Layer to read. Default `"hydrolocations"`, which is the name the
#'   pipeline reads; a set whose layer is called anything else is rejected unless
#'   this is set explicitly, because the build cannot find it.
#' @param references Named list of references to resolve against, each a list of
#'   `ids` (character `flowpath_id`) and `vpuids`. Built by
#'   [hl_reference_index()] when not supplied directly.
#' @param check_distance If `TRUE` (default), also measure how far each point
#'   lies from the reach it names. Resolution is not placement: an identifier can
#'   exist and still name a reach on the other side of the domain, which is what
#'   a statewide point set indexed against a clipped reference produces.
#' @param max_dist_m Reporting threshold for that distance, in metres. Default
#'   `200`. The result is informational, not a gate: across the NWM v4 inputs the
#'   population is nowhere zero, and a threshold tight enough to catch genuine
#'   mis-attribution also flags several thousand forecast points that sit a few
#'   hundred metres off the reference line and snap cleanly. Read the spread the
#'   check reports, a median near zero with a long tail is snapping noise, a
#'   median in the kilometres is a point set indexed against the wrong domain.
#' @param strict If `TRUE`, a failed gate raises an error. Default `FALSE`
#'   reports every check and returns.
#' @returns Invisibly, a list with `ok`, `stage` and one entry per check.
#' @examples
#' \dontrun{
#' refs <- list(conus = hl_reference_index("conus_reference.gpkg"))
#' check_hydrolocation_input("merged_hydrolocations.gpkg", references = refs)
#' }
#' @export
check_hydrolocation_input <- function(hl_gpkg, layer = "hydrolocations",
                                      references = list(), check_distance = TRUE,
                                      max_dist_m = 200, strict = FALSE) {
  stopifnot(is.character(hl_gpkg), length(hl_gpkg) == 1L, file.exists(hl_gpkg))
  checks <- list()

  # The layer name is a gate, not a convenience. The per-VPU slice reads a layer
  # called `hydrolocations`; a set that names it otherwise fails mid-build with
  # an st_write error on a NULL, which is a long way from the cause.
  present <- tryCatch(sf::st_layers(hl_gpkg)$name, error = function(e) character(0))
  checks$hl_layer_readable <- .hf_ok(layer %in% present, sprintf(
    "layer '%s' %s (present: %s)", layer,
    if (layer %in% present) "found" else "NOT FOUND -- the build reads this name",
    paste(present, collapse = ", ")))
  if (!(layer %in% present)) return(.hf_report_hl(checks, strict))

  hl   <- sf::read_sf(hl_gpkg, layer, quiet = TRUE)
  geom <- hl
  n    <- nrow(hl)

  gt <- unique(as.character(sf::st_geometry_type(hl)))
  checks$hl_geometry_is_point <- .hf_ok(all(gt %in% c("POINT", "MULTIPOINT")),
    sprintf("geometry: %s", paste(gt, collapse = ",")))
  ne <- sum(sf::st_is_empty(hl))
  checks$hl_no_empty_geometry <- .hf_ok(ne == 0L, sprintf("%d empty", ne))
  checks$hl_crs_defined <- .hf_ok(!is.na(sf::st_crs(hl)$Name),
    sprintf("CRS: %s", sf::st_crs(hl)$Name))

  hl <- sf::st_drop_geometry(hl)
  hr <- trimws(as.character(hl$hl_reference))
  nm <- sum(is.na(hr) | hr %in% c("", "NA"))
  checks$hl_reference_present <- .hf_ok(nm == 0L,
    sprintf("%d of %d row(s) have no hl_reference (dropped downstream)", nm, n))

  nfp <- sum(is.na(hl$flowpath_id))
  checks$hl_flowpath_id_assigned <- .hf_ok(nfp == 0L, sprintf(
    "%d of %d row(s) have no flowpath_id (stage 1 expects it pre-assigned)", nfp, n))

  # Exact duplication is not collocation: same reference, reach, type and class.
  key <- paste(hl$vpuid, hl$flowpath_id, hr, hl$hl_type, hl$hl_class, sep = "\r")
  nd  <- sum(duplicated(key))
  checks$hl_no_duplicate_rows <- .hf_ok(nd == 0L, sprintf(
    "%d exact duplicate row(s) over %d distinct key(s)", nd,
    length(unique(key[duplicated(key)]))))

  known <- unique(unlist(lapply(references, `[[`, "vpuids")))
  if (length(known)) {
    bad  <- setdiff(unique(as.character(hl$vpuid)), known)
    nbad <- sum(as.character(hl$vpuid) %in% bad)
    checks$hl_vpuid_known <- .hf_ok(length(bad) == 0L, sprintf(
      "%d row(s) carry a vpuid no reference has%s", nbad,
      if (length(bad)) paste0(": ", paste(sort(bad), collapse = ", ")) else ""))
  }

  for (d in names(references)) {
    ref <- references[[d]]
    sel <- as.character(hl$vpuid) %in% ref$vpuids
    if (!any(sel)) {
      checks[[paste0("hl_indexed_", d)]] <- .hf_info(sprintf("no rows tagged for %s", d))
      next
    }
    h   <- as.character(hl$flowpath_id[sel])
    hit <- sum(h %in% ref$ids)
    checks[[paste0("hl_indexed_", d)]] <- .hf_ok(hit == length(h), sprintf(
      "%d of %d %s row(s) resolve to a reference reach (%.1f%%)",
      hit, length(h), d, 100 * hit / length(h)))

    if (isTRUE(check_distance)) {
      dd <- .hl_named_distance(geom[sel, ], h, ref)
      if (!is.null(dd)) {
        far <- sum(dd > max_dist_m, na.rm = TRUE)
        q <- stats::quantile(dd, c(0.5, 0.95, 1), na.rm = TRUE)
        checks[[paste0("hl_near_named_reach_", d)]] <- .hf_info(sprintf(
          "%d of %d %s row(s) sit more than %g m from the reach they name (median %.0f m, p95 %.0f m, max %.0f m)",
          far, sum(!is.na(dd)), d, max_dist_m, q[1], q[2], q[3]))
      }
    }
  }

  dup_ref <- sum(duplicated(hr))
  checks$hl_reference_shared <- .hf_info(sprintf(
    "%d row(s) share an hl_reference, over %d distinct reference(s)",
    dup_ref, length(unique(hr[duplicated(hr)]))))
  np <- sum(!is.na(hl$poi_id))
  checks$hl_poi_id_prepopulated <- .hf_info(sprintf("%d of %d row(s) arrive with a poi_id%s",
    np, n, if (np > 0L && np < n) sprintf(" (vpuid: %s)",
      paste(sort(unique(hl$vpuid[!is.na(hl$poi_id)])), collapse = ",")) else ""))
  tt <- sort(table(hl$hl_type), decreasing = TRUE)
  checks$hl_type_mix <- .hf_info(paste(sprintf("%s=%d", names(tt), tt), collapse = " "))

  .hf_report_hl(checks, strict)
}

#' Index a reference for hydrolocation checking
#'
#' @description
#' Reads the `flowpath_id` and `vpuid` values a reference carries, in the shape
#' [check_hydrolocation_input()] expects. Attribute-only, so it does not read
#' geometry.
#'
#' @param reference_gpkg Path to a reference fabric.
#' @param layer Flowpaths layer name. Default `"flowpaths"`.
#' @returns A list with `ids` and `vpuids`.
#' @examples
#' \dontrun{
#' hl_reference_index("conus_reference.gpkg")
#' }
#' @export
hl_reference_index <- function(reference_gpkg, layer = "flowpaths") {
  stopifnot(file.exists(reference_gpkg))
  con <- DBI::dbConnect(RSQLite::SQLite(), reference_gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  cols <- names(DBI::dbGetQuery(con, sprintf('SELECT * FROM "%s" LIMIT 0', layer)))
  has_v <- "vpuid" %in% cols
  d <- DBI::dbGetQuery(con, sprintf('SELECT flowpath_id AS id%s FROM "%s"',
    if (has_v) ", vpuid" else "", layer))
  list(ids = as.character(d$id),
       vpuids = if (has_v) sort(unique(as.character(d$vpuid))) else character(0),
       path = reference_gpkg, layer = layer)
}

# How far is each point from the reach it names?
#
# Resolution is not placement. An identifier can exist in the reference and name
# a reach on the other side of the domain, which is what a statewide point set
# indexed against a clipped reference produces: valid ids, a median 129 km away.
# Only the named reach is read, not the whole layer, so the cost is one filtered
# query rather than a nearest-neighbour search over millions of reaches.
.hl_named_distance <- function(pts, ids, ref) {
  if (is.null(ref$path) || !file.exists(ref$path)) return(NULL)
  want <- unique(ids[!is.na(ids)])
  if (!length(want)) return(NULL)
  chunk <- split(want, ceiling(seq_along(want) / 900))
  geo <- NULL
  for (k in chunk) {
    q <- sprintf('SELECT flowpath_id, geom FROM "%s" WHERE flowpath_id IN (%s)',
      ref$layer, paste(sprintf("'%s'", k), collapse = ","))
    g <- tryCatch(sf::read_sf(ref$path, query = q, quiet = TRUE), error = function(e) NULL)
    if (is.null(g) || !nrow(g)) next
    geo <- if (is.null(geo)) g else rbind(geo, g)
  }
  if (is.null(geo) || !nrow(geo)) return(NULL)
  pts <- sf::st_transform(pts, sf::st_crs(geo))
  i <- match(as.character(ids), as.character(geo$flowpath_id))
  d <- rep(NA_real_, length(ids))
  ok <- !is.na(i)
  if (any(ok)) d[ok] <- as.numeric(sf::st_distance(
    sf::st_geometry(pts)[ok], sf::st_geometry(geo)[i[ok]], by_element = TRUE))
  d
}

.hf_report_hl <- function(checks, strict) {
  ok_all <- all(vapply(checks, function(x) isTRUE(x$ok), logical(1L)))
  for (nm in names(checks)) {
    res  <- checks[[nm]]
    icon <- if (identical(res$kind, "info")) "i" else if (isTRUE(res$ok)) "ok" else "x"
    line <- sprintf("[hydrolocations] %s %-26s %s", icon, nm, res$msg)
    if (isTRUE(res$ok) || identical(res$kind, "info")) message(line)
    else if (strict) stop(line, call. = FALSE) else message(line)
  }
  invisible(list(ok = ok_all, stage = "hydrolocations", checks = checks))
}

#' Run the end-to-end hydrofabric pipeline for a single VPU
#'
#' @description
#' Local convenience wrapper that runs the refactor, aggregate and nextgen phases sequentially on the
#' same filesystem. Each stage is also exported and can be called directly.
#'
#' @param vpuid VPU identifier (character, e.g. `"01"`, `"03N"`, `"21"`).
#' @param reference_gpkg Path to the reference fabric (multi-VPU or
#'   single-VPU; extraction is auto-detected).
#' @param fdr_path,fac_path Paths to the FDR / FAC rasters.
#' @param fallback_fdr_path,fallback_fac_path Optional paths to a fallback
#'   FDR/FAC pair.  `NULL` disables fallback (default).
#' @param exclude_hl_types Hydrolocation `hl_type` values to drop before the
#'   build sees them. The default removes waterbody INLETS (`wbin`), 33,926 of
#'   the 69,036 supplied locations. An inlet names no computational boundary a
#'   consumer asks for, and each one costs twice: it pins its own unit against
#'   merging and breaks the levelpath, stranding the unit below it too, which is
#'   why they dominate the sub-minimum population.
#'
#'   Waterbody OUTLETS (`wbout`) are kept. An outlet is a real routing boundary
#'  , it is where a lake discharges, so it earns its break in a way an inlet
#'   does not. Gages (`nwis`, `rfc`) and terminal/coastal outlets (the `owp_*`
#'   types) are kept for the same reason. Lake handling is independent of this
#'   either way: lakes come from `lakes_path`, a separate input.
#'
#'   `NULL` or `character(0)` keeps every type. Filtering is by row, so a
#'   location carrying a wanted type as well survives on that type.
#' @param hydrolocations_gpkg Optional path to a separate GeoPackage holding the
#'   `hydrolocations` (POIs) + `pois` layers. `NULL` (default) reads them from
#'   `reference_gpkg` (all-in-one). When supplied, the reference need only carry
#'   `flowpaths` + `divides`, and POIs are sourced from this file.
#' @param out_dir Directory for outputs. A `vpu_{vpuid}/` subfolder is
#'   created; the stage outputs all land there so this call is
#'   backwards-compatible with the monolithic pre-split runs.
#' @param n_cores Worker count for parallel steps. Default
#'   `detectCores() - 2`.
#' @param from,through Stage range to run, by name (`"refactor"`, `"aggregate"`,
#'   `"nextgen"`, in order). Defaults run the full pipeline. Use
#'   `through = "aggregate"` to stop before the ngen output, or
#'   `from = "aggregate"` to resume from an existing refactored fabric.
#' @param stages Optional integer stage numbers (e.g. `c(2, 3)`) as an alias for
#'   `from`/`through`: when supplied, `from`/`through` are set to the min/max
#'   stage. `NULL` (default) uses `from`/`through` directly.
#' @param skip_if_exists If `TRUE`, skip a VPU whose target output (the ngen gpkg
#'   when `through = "nextgen"`, otherwise the aggregate output gpkg) already exists.
#'   Default `FALSE` (always rebuild). This is a plain already-built check, not
#'   content-hash reuse.
#' @param verbose Pass through to individual steps for extra logging.
#' @param hydraulics_path Optional path to a RiverML channel table parquet,
#'   either its v2 output (`tw_bf_m`, `y_bf_m`, `r_bf`, roughness, `slope_m_m`)
#'   or the earlier channel-geometry table. When supplied, routing attributes (`flowpath-attributes`,
#'   `flowline-attributes`) are computed in Step 6g, and a `slope_m_m` column
#'   also supplies flowpath slope when `slope_parquet` is `NULL`. `NULL` skips
#'   Step 6g (default).
#' @param lakes_path Optional path to the NWM lakes GeoPackage; forwarded to
#'   [nextgen_fabric()] (Step 6h). `NULL` skips lakes (default).
#' @param lake_mask_path Optional path to the lake-mask parquet (reach ->
#'   waterbody membership); forwarded to [nextgen_fabric()]. Required when
#'   `lakes_path` is supplied.
#' @param slope_parquet Optional path to a per-reference-id slope parquet;
#'   forwarded to [nextgen_fabric()]. `NULL` (default) takes slope from
#'   `hydraulics_path` when that table carries a `slope_m_m` column, and
#'   otherwise leaves slope `NA`.
#' @param use_checkpoints If `TRUE`, allow stage checkpoint reuse where
#'   supported (currently the refactor phase). Default `FALSE` recomputes from scratch.
#' @param check_input If `TRUE` (default), validate the inputs before building:
#'   [check_reference_input()] aborts on damaged identifiers or unresolvable
#'   topology, and [check_grid_coverage()] aborts when a primary grid does not
#'   cover the domain, which is what a configuration inheriting another domain's
#'   grid looks like.
#' @param ... Stage tuning parameters, forwarded to whichever stage declares
#'   each (e.g. `min_area_sqkm = 5`, `coverage_min = 0.95`). See
#'   [refactor_fabric()], [aggregate_fabric()], [nextgen_fabric()] for the full
#'   set and defaults.
#'
#' @return Invisibly, a list with `vpuid`, `status`, `ngen_gpkg`,
#'   `out_gpkg`, `timings`, `invariants`, `skipped`.
#'
#' @examples
#' \dontrun{
#' run_pipeline(vpuid          = "01",
#'              reference_gpkg = "reference_01.gpkg",
#'              fdr_path       = "fdr.tif",
#'              fac_path       = "fac.tif",
#'              out_dir        = "out",
#'              through        = "nextgen")
#' }
#'
#' @export
run_pipeline <- function(vpuid,
                         reference_gpkg,
                         fdr_path,
                         fac_path,
                         fallback_fdr_path = NULL,
                         fallback_fac_path = NULL,
                         hydrolocations_gpkg = NULL,
                         exclude_hl_types = "wbin",
                         out_dir        = "out",
                         from           = "refactor",
                         through        = "nextgen",
                         stages         = NULL,
                         n_cores        = NULL,
                         skip_if_exists = FALSE,
                         verbose        = FALSE,
                         hydraulics_path = NULL,
                         lakes_path      = NULL,
                         lake_mask_path  = NULL,
                         slope_parquet   = NULL,
                         use_checkpoints = FALSE,
                         check_input    = TRUE,
                         ...) {
  stopifnot(
    is.character(vpuid), length(vpuid) == 1L, nzchar(vpuid),
    file.exists(reference_gpkg),
    .raster_accessible(fdr_path),
    .raster_accessible(fac_path))

  # Validate the reference before spending a build on it. Truncated identifiers
  # and unresolvable topology are silent at read time but produce a fabric that
  # passes its own invariants while being unaggregated, so the cost of finding
  # them afterwards is a whole run.
  if (isTRUE(check_input)) {
    check_reference_input(reference_gpkg, strict = TRUE)
    # A grid meant for another domain reads fine and covers nothing here, which
    # the partition absorbs rather than reports. Compare extents up front.
    check_grid_coverage(reference_gpkg, fdr_path, fac_path,
      fallback_fdr_path, fallback_fac_path, vpuid = vpuid, strict = TRUE)
  }

  if (is.null(n_cores))
    n_cores <- max(1L, parallel::detectCores() - 2L)

  vpu_dir   <- file.path(out_dir, paste0("vpu_", vpuid))
  dir.create(vpu_dir, recursive = TRUE, showWarnings = FALSE)
  work_dir  <- file.path(vpu_dir, "work")
  out_gpkg  <- file.path(work_dir, sprintf("out_vpu%s.gpkg",  vpuid))
  ngen_gpkg <- file.path(vpu_dir, sprintf("ngen_vpu%s.gpkg", vpuid))

  # Resolve the stage range -> integer `stages` (feeds the N %in% stages gating below).
  .STAGES <- c("refactor", "aggregate", "nextgen")
  # `stages` (e.g. c(2, 3)) is an alias that sets from/through by stage number,
  # so callers can request a stage subset without naming the endpoints.
  if (!is.null(stages)) {
    .si <- suppressWarnings(as.integer(stages))
    .si <- .si[!is.na(.si) & .si >= 1L & .si <= length(.STAGES)]
    if (length(.si)) { from <- .STAGES[min(.si)]; through <- .STAGES[max(.si)] }
  }
  from    <- match.arg(from,    .STAGES)
  through <- match.arg(through, .STAGES)
  i1 <- match(from, .STAGES); i2 <- match(through, .STAGES)
  if (i1 > i2) stop("`from` must not come after `through`", call. = FALSE)
  stages <- seq.int(i1, i2)

  # Tuning params flow via ... to whichever stage(s) declare them.
  dots    <- list(...)
  known   <- unique(unlist(lapply(
    list(refactor_fabric, aggregate_fabric, nextgen_fabric), .stage_param_names)))
  unknown <- setdiff(names(dots), known)
  if (length(unknown))
    .hf_log_warn("run_pipeline: ignoring unrecognized parameter(s): ",
      paste(unknown, collapse = ", "))
  p1 <- .stage_tuning(refactor_fabric,  dots)
  p2 <- .stage_tuning(aggregate_fabric, dots)
  p3 <- .stage_tuning(nextgen_fabric,   dots)

  # -- Version-guard checkpoint reuse ----------------------------------------
  # Honour use_checkpoints only when the existing checkpoints were produced by
  # the SAME build + params. A reinstall (code fix) or param change flips the
  # fingerprint and forces a clean recompute, so a stale intermediate is never
  # silently reused (an OOM rerun on unchanged code resumes; a bug-fix rerun does
  # not). Then (re)stamp so the next same-build rerun can resume.
  .fp <- .pipeline_fingerprint(list(p1 = p1, p2 = p2, p3 = p3))
  if (use_checkpoints && !.checkpoints_valid(work_dir, .fp)) {
    if (dir.exists(work_dir) &&
        length(list.files(work_dir, pattern = "\\.gpkg$")) > 0L)
      message(sprintf(
        "[%s] pipeline fingerprint changed -> ignoring stale checkpoints (clean recompute)",
        vpuid))
    use_checkpoints <- FALSE
  }
  .stamp_fingerprint(work_dir, .fp)

  # hydrofabric *builds*. Provenance-gated reuse (content-hash skip) and output
  # stamping live in a separate layer; here skip_if_exists is a plain
  # already-built check on the target output.
  skip_target <- if (3L %in% stages) ngen_gpkg else out_gpkg
  if (isTRUE(skip_if_exists) && file.exists(skip_target)) {
    message(sprintf("[%s] skipping -- %s output already exists: %s",
      vpuid, if (3L %in% stages) "ngen" else "stage2", skip_target))
    return(invisible(list(vpuid = vpuid, status = "ok",
      ngen_gpkg = if (3L %in% stages) ngen_gpkg else NA_character_,
      out_gpkg = out_gpkg,
      timings = list(), invariants = NULL, skipped = TRUE)))
  }

  # All stages share one vpu_dir (fabric_dir == out_dir/vpu_<id>); run_pipeline()
  # keeps the historic single-folder layout (ref / out / ngen gpkgs together).
  t_total <- proc.time()[["elapsed"]]

  # Returns TRUE if any check_invariants result in the list has ok = FALSE.
  .inv_failed <- function(invs) {
    any(vapply(invs, function(x) is.list(x) && isFALSE(x$ok), logical(1)))
  }

  r1 <- if (1L %in% stages)
    do.call(refactor_fabric, c(list(
      vpuid = vpuid,
      reference_gpkg    = reference_gpkg,
      fdr_path          = fdr_path,
      fac_path          = fac_path,
      fallback_fdr_path = fallback_fdr_path,
      fallback_fac_path = fallback_fac_path,
      hydrolocations_gpkg = hydrolocations_gpkg,
      exclude_hl_types = exclude_hl_types,
      out_dir  = out_dir,
      n_cores  = n_cores,
      use_checkpoints = use_checkpoints), p1))
  else list(timings = list(), invariants = list())

  if (1L %in% stages && .inv_failed(r1$invariants)) {
    message(sprintf("[%s] Stage 1 invariants failed -- halting before Stage 2/3", vpuid))
    all_invariants <- r1$invariants
    jsonlite::write_json(all_invariants,
      file.path(work_dir, "invariants.json"),
      auto_unbox = TRUE, pretty = TRUE, force = TRUE)
    timings <- c(r1$timings, TOTAL = proc.time()[["elapsed"]] - t_total)
    saveRDS(timings, file.path(work_dir, sprintf("timings_vpu%s.rds", vpuid)))
    return(invisible(list(vpuid = vpuid, status = "failed",
      ngen_gpkg = NA_character_, out_gpkg = out_gpkg,
      timings = timings, invariants = all_invariants,
      skipped = FALSE)))
  }

  r2 <- if (2L %in% stages)
    do.call(aggregate_fabric, c(list(vpuid = vpuid,
      fabric_dir = vpu_dir, out_dir = out_dir,
      verbose = verbose), p2))
  else list(timings = list(), invariants = NULL)

  all_invariants <- c(r1$invariants,
    if (!is.null(r2$invariants))
      list(aggregated = r2$invariants))

  if (2L %in% stages && .inv_failed(list(r2$invariants))) {
    message(sprintf("[%s] Stage 2 invariants failed -- halting before Stage 3", vpuid))
    jsonlite::write_json(all_invariants,
      file.path(work_dir, "invariants.json"),
      auto_unbox = TRUE, pretty = TRUE, force = TRUE)
    timings <- c(r1$timings, r2$timings, TOTAL = proc.time()[["elapsed"]] - t_total)
    saveRDS(timings, file.path(work_dir, sprintf("timings_vpu%s.rds", vpuid)))
    return(invisible(list(vpuid = vpuid, status = "failed",
      ngen_gpkg = NA_character_, out_gpkg = out_gpkg,
      timings = timings, invariants = all_invariants,
      skipped = FALSE)))
  }

  if (3L %in% stages) {
    r3 <- do.call(nextgen_fabric, c(list(vpuid = vpuid,
      fabric_dir = vpu_dir,
      fdr_path          = fdr_path,
      fac_path          = fac_path,
      fallback_fdr_path = fallback_fdr_path,
      fallback_fac_path = fallback_fac_path,
      out_dir  = out_dir,
      n_cores  = n_cores,
      hydraulics_path = hydraulics_path,
      lakes_path      = lakes_path,
      lake_mask_path  = lake_mask_path,
      slope_parquet   = slope_parquet,
      reference_gpkg  = reference_gpkg,
      verbose = verbose), p3))
    all_invariants <- c(all_invariants, list(ngen = r3$invariants))
    ngen_out <- ngen_gpkg
    timings  <- c(r1$timings, r2$timings, r3$timings)
  } else {
    ngen_out <- NA_character_
    timings  <- c(r1$timings, r2$timings)
    # Stage 3 normally writes invariants.json; write it here instead so the
    # file always reflects all stages that actually ran.
    jsonlite::write_json(all_invariants,
      file.path(work_dir, "invariants.json"),
      auto_unbox = TRUE, pretty = TRUE, force = TRUE)
  }

  timings[["TOTAL"]] <- proc.time()[["elapsed"]] - t_total
  saveRDS(timings, file.path(work_dir, sprintf("timings_vpu%s.rds", vpuid)))

  # Stamp every GeoPackage this run wrote, not only the two headline outputs.
  # The work directory holds derived data under the same license, and an
  # intermediate that leaves the machine (handed to a colleague, attached to a
  # ticket, copied into a subset) carries no notice unless it was stamped where
  # it was written. Listing the directory rather than naming files means a new
  # intermediate is covered the day it appears.
  .hf_stamp_outputs(
    unique(c(ngen_out, out_gpkg,
             list.files(work_dir, pattern = "[.]gpkg$", full.names = TRUE))),
    vpuid = vpuid,
    grid_paths = c(fdr_path, fac_path, fallback_fdr_path, fallback_fac_path),
    lakes_path = lakes_path)

  invisible(list(
    vpuid      = vpuid,
    status     = "ok",
    ngen_gpkg  = ngen_out,
    out_gpkg   = out_gpkg,
    timings    = timings,
    invariants = all_invariants,
    skipped    = FALSE))
}

# --- internals --------------------------------------------------------

# Forward run_pipeline(...) tuning params to the stage function(s) that declare
# them. A stage's "tuning" params are its formals minus the fixed control args
# (paths, dirs, workers). Defaults live in each stage function's signature.
.stage_control_args <- c(
  "vpuid", "reference_gpkg", "fdr_path", "fac_path", "fallback_fdr_path",
  "fallback_fac_path", "hydrolocations_gpkg", "exclude_hl_types", "fabric_dir",
  "out_dir", "n_cores", "use_checkpoints", "verbose", "hydraulics_path",
  "lakes_path", "lake_mask_path", "slope_parquet", "...")

.stage_param_names <- function(fn)
  setdiff(names(formals(fn)), .stage_control_args)

.stage_tuning <- function(fn, dots) {
  keys     <- .stage_param_names(fn)
  defaults <- lapply(formals(fn)[keys], eval, envir = environment(fn))
  utils::modifyList(defaults, dots[intersect(names(dots), keys)])
}

# Data-license / attribution policy (see THIRD_PARTY_NOTICES.md).
#
# License is UNIFORM: every deliverable is CC BY-SA 4.0. This is forced by the
# most restrictive input (TDX-Hydro is CC BY-SA 4.0, ShareAlike) and adopted
# everywhere so the obligation is always satisfiable; Lynker controls the
# reference-fabric tier and elects the same license for consistency.
#
# Attribution is PER-SOURCE, and the exceptions are the domains that touch
# neither upstream source. Hawaii (20) and PRVI (21) are built from NHDPlus HR
# throughout, network and grids alike, which is a US Government work in the
# public domain and requires no attribution, so they carry the Lynker Corporation
# credit alone. Every other domain touches TDX-Hydro and HydroSHEDS v2 in some
# capacity, so both credits are owed. This governs credit only, the license is
# uniform across the product.
#
# Listing the exceptions rather than the users is deliberate. A new domain
# arrives credited and has to be argued down, which fails safe; the reverse
# ships uncredited until someone notices.
.HF_PUBLIC_DOMAIN_VPUS <- c("20", "21")  # Hawaii, PRVI

.HF_ATTR_LS  <- paste0(
  "Derived from the Lynker Spatial Reference Fabric \u00a9 Lynker Corporation, ",
  "licensed under CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/). ",
  "Aggregated to the NextGen hydrofabric data model by Lynker Spatial.")
# NGA's BY-SA notice: credit + license link + an accurate (not minimizing)
# statement of changes. The fuller account of Lynker's work lives in the
# provenance block / dataset description, not in this upstream-credit line.
.HF_ATTR_TDX <- paste0(
  "Contains TDX-Hydro \u00a9 2023 National Geospatial-Intelligence Agency (NGA), ",
  "licensed under CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/), ",
  "used as a base network. Modified and substantially transformed by Lynker ",
  "Spatial. Derived data \u00a9 Lynker Corporation, distributed under CC BY-SA 4.0.")

# Elevation-grid credits, keyed on the raster actually used.
#
# The FDR/FAC grids carry their own licenses and both current sources are CC BY 4.0,
# which requires credit. Deriving the notice from the raster PATHS rather than from a
# flag applies the same principle as the NGA credit above: it cannot be forgotten on
# a build that uses the grid, and cannot leak onto one that does not.
#
# The GEDTM30 credit is retained although no current domain reads that grid: the
# notice is keyed on the configured paths, so it re-arms by itself if a build
# names GEDTM30 again.
.HF_ATTR_GEDTM <- paste0(
  "Contains GEDTM30 terrain derivatives \u00a9 OpenLandMap / Hengl et al., ",
  "licensed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).")
.HF_ATTR_HYDROSHEDS <- paste0(
  "Contains HydroSHEDS v2 flow-direction and flow-accumulation data \u00a9 World ",
  "Wildlife Fund, Inc., licensed under CC BY 4.0 ",
  "(https://creativecommons.org/licenses/by/4.0/); derived from TanDEM-X ",
  "\u00a9 DLR.")

# Lake-polygon credit, keyed on the lake set actually read. The lake layer
# records its origin in a `source` column, so the credit follows the data rather
# than a domain rule: Alaska's lakes are HydroLAKES today, and a domain that
# switches to HydroLAKES is credited without anyone remembering to add it.
.HF_ATTR_HYDROLAKES <- paste0(
  "Contains HydroLAKES v1.0 lake polygons \u00a9 World Wildlife Fund, Inc. ",
  "(Messager et al., 2016), licensed under CC BY 4.0 ",
  "(https://creativecommons.org/licenses/by/4.0/).")

.hf_lakes_attribution <- function(lakes_path) {
  if (is.null(lakes_path) || is.na(lakes_path) || !file.exists(lakes_path))
    return(character(0))
  if (grepl("hydrolakes", tolower(basename(lakes_path)), fixed = TRUE))
    return(.HF_ATTR_HYDROLAKES)
  src <- tryCatch({
    con <- DBI::dbConnect(RSQLite::SQLite(), lakes_path)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    lyr <- DBI::dbGetQuery(con,
      "SELECT table_name FROM gpkg_contents WHERE data_type = 'features'")$table_name
    unlist(lapply(lyr, function(l)
      if ("source" %in% DBI::dbListFields(con, l))
        DBI::dbGetQuery(con, sprintf('SELECT DISTINCT source FROM "%s"', l))$source))
  }, error = function(e) character(0))
  if (any(grepl("hydrolakes", tolower(src), fixed = TRUE))) .HF_ATTR_HYDROLAKES
  else character(0)
}

# Match on the distributed file/directory names of each product. Kept deliberately
# loose (the VRT wrapper may be renamed locally, as gedtm_fac.vrt is) so a rename
# cannot silently drop a required credit.
# Grid and lake credits carried forward from a merged product's constituent VPUs.
#
# Stage 4 reads no rasters or lake sets, so it cannot derive those credits from paths
# the way a per-VPU build does. Reading them back out of the inputs it merges
# keeps the obligation with the data: a domain product is a redistribution of
# whatever its VPUs contain, so a CC BY 4.0 grid credited on a VPU must appear on
# the domain too. Matching the exact constants, rather than parsing prose, means
# a credit can only be carried forward if this package wrote it.
.hf_grid_attribution_inherited <- function(gpkgs) {
  known <- c(.HF_ATTR_GEDTM, .HF_ATTR_HYDROSHEDS, .HF_ATTR_HYDROLAKES)
  found <- character(0)
  for (g in gpkgs[file.exists(gpkgs)]) {
    # Each connection is closed in its own iteration. on.exit() would register
    # against this function's exit, so every registration would fire on the last
    # connection and leak the rest.
    con <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), g), error = function(e) NULL)
    if (is.null(con)) next
    a <- tryCatch({
      q <- DBI::dbGetQuery(con, paste0(
        "SELECT metadata FROM gpkg_metadata ",
        "WHERE mime_type = 'application/json' AND metadata LIKE '%attribution%'"))
      if (!nrow(q)) "" else
        paste(vapply(q$metadata, function(m)
          tryCatch(jsonlite::fromJSON(m)$attribution %||% "",
            error = function(e) ""), character(1L)), collapse = " ")
    }, error = function(e) "")
    try(DBI::dbDisconnect(con), silent = TRUE)
    found <- c(found, known[vapply(known, grepl, logical(1L), x = a, fixed = TRUE)])
  }
  unique(found)
}

.hf_grid_attribution <- function(paths) {
  p <- tolower(paths[!is.na(paths) & nzchar(paths)])
  if (!length(p)) return(character(0))
  out <- character(0)
  if (any(grepl("gedtm", p, fixed = TRUE)))        out <- c(out, .HF_ATTR_GEDTM)
  if (any(grepl("hydrosheds", p, fixed = TRUE)) ||
      any(grepl("_acc_1s_|_dir_1s_", p)))          out <- c(out, .HF_ATTR_HYDROSHEDS)
  out
}

# Lynker's value-add, recorded in the provenance block of every deliverable
# (not just TDX builds). This is where the processing story lives, keeping the
# upstream attribution line factual and avoiding a sprawl of sidecar files.
.HF_DESCRIPTION <- paste0(
  "NextGen hydrofabric produced by Lynker Spatial: flowline refactoring and ID ",
  "reconciliation, network topology reconstruction, catchment (divide) ",
  "derivation and area-conserving repair, POI/nexus placement, and level-path ",
  "aggregation. The sources each product draws on are named in its attribution.")

# Resolve the attribution string for a build. An explicit HF_DATA_ATTRIBUTION
# always wins (operator override); otherwise it is derived from the source:
# TDX-sourced VPUs credit NGA, everything else is Lynker-Spatial-only. Deriving
# it here (rather than relying on a driver to set the env var) makes the NGA
# credit impossible to forget on an AK build and impossible to leak onto a
# non-TDX one.
# One data license across the national product.
#
# CC BY-SA is owed where TDX-Hydro is in the lineage, since it is ShareAlike and
# an adaptation must carry it. Hawaii and PRVI inherit no such obligation, they
# derive from NHDPlus HR, a US Government work in the public domain, so their
# term is a choice rather than a requirement, and the choice is to match the rest
# of the product.
#
# Licensing a derivative of public-domain material is routine and lawful; what
# would not be is claiming the source itself, or implying it is anything other
# than freely available. NOTICE states plainly that NHDPlus HR is public domain
# and remains so, which is the line that matters. CC0 was considered and rejected:
# it is an irrevocable waiver of whatever rights subsist in the selection and
# aggregation, which is a larger step than simply not asserting them.
.hf_data_license <- function(vpuid = NA_character_) {
  env <- Sys.getenv("HF_DATA_LICENSE", "")
  if (nzchar(env)) return(env)
  "CC-BY-SA-4.0"
}

.hf_attribution <- function(vpuid = NA_character_) {
  env <- Sys.getenv("HF_DATA_ATTRIBUTION", "")
  if (nzchar(env)) return(env)
  # Match on the numeric VPU prefix so sub-VPU tiles (e.g. AK 19a/19b/19c)
  # resolve to their parent domain's attribution.
  v <- sub("[A-Za-z].*", "", toupper(trimws(as.character(vpuid))))
  # An unknown vpuid credits both sources, since the credit required depends
  # on the domain.
  if (!is.na(vpuid) && nzchar(v) && v %in% .HF_PUBLIC_DOMAIN_VPUS)
    return(.HF_ATTR_LS)
  paste(.HF_ATTR_TDX, .HF_ATTR_HYDROSHEDS)
}

# Stamp the hydrofabric release version + build provenance into output
# GeoPackages via hfutils. Best-effort: a metadata failure must never fail a
# build, so each stamp is wrapped in tryCatch. The release version is taken from
# HF_BUILD_VERSION when the driver supplies it, else the installed package
# version; HF_GIT_SHA is recorded as provenance when present.
.hf_stamp_outputs <- function(gpkgs, vpuid = NA_character_,
                              description = NULL, attribution = NULL,
                              grid_paths = character(0), lakes_path = NULL,
                              license = NULL) {
  gpkgs <- gpkgs[!is.na(gpkgs) & nzchar(gpkgs)]
  gpkgs <- gpkgs[file.exists(gpkgs)]
  if (!length(gpkgs)) return(invisible())

  ver <- Sys.getenv("HF_BUILD_VERSION", "")
  if (!nzchar(ver)) ver <- as.character(utils::packageVersion("hydrofabric"))
  sha <- Sys.getenv("HF_GIT_SHA", "")
  # Uniform license; per-source attribution (policy documented above). A merged
  # multi-VPU product (stage 4) passes `attribution` explicitly, it spans many
  # VPUs, so a single vpuid can't derive the credit, and leaves vpuid unset.
  # A merged product spans many VPUs, so it passes its license explicitly rather
  # than having one derived from a vpuid it does not have.
  lic <- if (!is.null(license)) license else .hf_data_license(vpuid)
  if (is.null(description)) description <- .HF_DESCRIPTION
  if (is.null(attribution)) attribution <- .hf_attribution(vpuid)
  # Append the grid credits for whichever elevation products this build read.
  # The domain policy already credits HydroSHEDS everywhere but Hawaii and PRVI,
  # so only add a notice the string does not already carry: a credit repeated
  # twice in one attribution reads as a defect in the stamp.
  ga <- c(.hf_grid_attribution(grid_paths), .hf_lakes_attribution(lakes_path))
  ga <- ga[!vapply(ga, grepl, logical(1L), x = attribution, fixed = TRUE)]
  if (length(ga)) attribution <- paste(c(attribution, ga), collapse = " ")

  prov <- list(
    software    = paste("hydrofabric", as.character(utils::packageVersion("hydrofabric"))),
    git_sha     = if (nzchar(sha)) sha else NULL,
    build_date  = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    description = description,
    attribution = attribution,
    vpuid       = if (!is.na(vpuid)) vpuid else NULL)

  for (g in gpkgs)
    tryCatch(hfutils::gpkg_set_version(g, ver, provenance = prov, license = lic),
      error = function(e)
        message(sprintf("[%s] version stamp skipped for %s: %s",
          vpuid, basename(g), conditionMessage(e))))
  invisible()
}

# A poi_id is a single integer by contract, and the same value is the key that
# joins a POI event to its reconciled reach. The AK OWP hydrolocation set ships
# a text handle ("wbin-2-409089") instead, on the hydrolocations and pois layers
# alike.
#
# Both layers are renumbered together, here at ingestion, so every table derived
# downstream carries the same integer. Rewriting one table further in strands its
# POIs on an NA reach and the layer is lost outright, the handle is not just an
# identifier, it is the join key. Numbering starts above the highest integer
# already in use so a renumbered id cannot collide with a real one.
.normalize_poi_id <- function(lyrs) {
  has  <- function(x) !is.null(x) && "poi_id" %in% names(x)
  keys <- c("hydrolocations", "pois")
  vals <- unique(unlist(lapply(lyrs[keys],
    function(x) if (has(x)) as.character(x$poi_id))))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (!length(vals)) return(lyrs)
  ints <- suppressWarnings(as.integer(vals))
  bad  <- vals[is.na(ints)]
  if (!length(bad)) return(lyrs)
  map <- stats::setNames(
    as.character(max(c(0L, ints[!is.na(ints)])) + seq_along(bad)), bad)
  message(sprintf(
    "Step 0: renumbered %d non-integer poi_id value(s) (e.g. %s -> %s)",
    length(bad), bad[1L], map[[1L]]))
  for (l in keys) {
    if (!has(lyrs[[l]])) next
    v <- as.character(lyrs[[l]]$poi_id)
    i <- match(v, names(map))
    v[!is.na(i)] <- unname(map[i[!is.na(i)]])
    lyrs[[l]]$poi_id <- v
  }
  lyrs
}

# Assemble the single-VPU input slice Stage 1 reads. Sources flowpaths/divides
# from `reference_gpkg` and hydrolocations/pois from `hydrolocations_gpkg` (which
# may be the SAME file, the default, or a separate hydrolocations GeoPackage).
# Each layer is sliced by `vpuid` when that layer carries the column, else read
# whole. Fast path: a single-file reference that is already one VPU (or has no
# vpuid) is returned as-is, exactly as before.
# Drop hydrolocations of unwanted `hl_type` before anything downstream sees
# them. Filtering by ROW rather than by location: a location recorded under two
# types keeps the one that was asked for, and only locations left with no rows
# disappear. Applied to `pois` too when that layer is present and typed, so the
# two stay in step: `.normalize_poi_id` renumbers them together immediately
# after, and a poi_id surviving in one layer but not the other would desync it.
.drop_hl_types <- function(poi_lyrs, exclude_hl_types, vpuid = NA_character_) {
  hl <- poi_lyrs[["hydrolocations"]]
  if (is.null(hl) || !nrow(hl) || !"hl_type" %in% names(hl)) return(poi_lyrs)
  drop <- as.character(hl$hl_type) %in% exclude_hl_types
  if (!any(drop)) return(poi_lyrs)
  poi_lyrs[["hydrolocations"]] <- hl[!drop, , drop = FALSE]
  # poi_id is not populated on the way in: .normalize_poi_id assigns it
  # immediately after this, so it identifies a location only when the input
  # already carried one. Matching on an all-NA column would take the pois layer
  # with it, so the companion filter runs only when there is something to match.
  kept_poi <- unique(stats::na.omit(as.character(hl$poi_id)[!drop]))
  po <- poi_lyrs[["pois"]]
  if (length(kept_poi) && !is.null(po) && nrow(po) && "poi_id" %in% names(po))
    poi_lyrs[["pois"]] <- po[as.character(po$poi_id) %in% kept_poi, , drop = FALSE]
  message(sprintf("[stage1][%s] dropped %d hydrolocation(s) of type %s; %d kept",
    vpuid, sum(drop), paste(exclude_hl_types, collapse = "/"), sum(!drop)))
  # An allow list, not a deny list: the build consumes the hl_types it is
  # specified to consume, so a type absent from the specification is excluded
  # by default rather than passed through.
  keep <- poi_lyrs[["hydrolocations"]]
  if ("hl_class" %in% names(keep)) {
    tc <- sort(table(as.character(keep$hl_class)), decreasing = TRUE)
    message(sprintf("[stage1][%s] hydrolocation classes kept: %s", vpuid,
      paste(sprintf("%s=%d", names(tc), as.integer(tc)), collapse = " ")))
  }
  poi_lyrs
}

.extract_vpu_slice_if_needed <- function(reference_gpkg, out_slice, vpuid,
                                         .t_fun, hydrolocations_gpkg = NULL,
                                         exclude_hl_types = NULL) {
  if (is.null(hydrolocations_gpkg)) hydrolocations_gpkg <- reference_gpkg
  exclude_hl_types <- exclude_hl_types[!is.na(exclude_hl_types)]
  separate_hl <- !identical(
    normalizePath(hydrolocations_gpkg, mustWork = FALSE),
    normalizePath(reference_gpkg,      mustWork = FALSE))

  con  <- DBI::dbConnect(RSQLite::SQLite(), reference_gpkg)
  cols <- DBI::dbListFields(con, "flowpaths")
  has_vpuid <- "vpuid" %in% cols
  distinct_vpus <- if (has_vpuid)
    DBI::dbGetQuery(con, "SELECT DISTINCT vpuid FROM flowpaths")$vpuid else character(0)
  DBI::dbDisconnect(con)
  ref_single <- !has_vpuid ||
    (length(distinct_vpus) == 1L && distinct_vpus[1L] == vpuid)
  # Only the all-in-one, already-single-VPU case can be used verbatim, and
  # only when nothing has to be filtered out of it on the way through.
  if (!separate_hl && ref_single && !length(exclude_hl_types))
    return(reference_gpkg)

  if (file.exists(out_slice)) return(out_slice)

  # Read a layer from `g`, sliced to `vpuid` when that layer carries vpuid.
  # Returns NULL when the layer is absent (e.g. a reference with no pois).
  read_slice <- function(g, lyr) {
    if (!hfutils::layer_exists(g, lyr)) return(NULL)
    lyr_cols <- names(sf::st_read(g,
      query = sprintf('SELECT * FROM "%s" LIMIT 0', lyr), quiet = TRUE))
    if ("vpuid" %in% lyr_cols)
      sf::st_read(g, query = sprintf('SELECT * FROM "%s" WHERE vpuid = \'%s\'',
        lyr, vpuid), quiet = TRUE)
    else
      sf::st_read(g, lyr, quiet = TRUE)
  }

  .t_fun(sprintf("Step 0: extract VPU slice (%s)", vpuid), {
    srcs <- list(flowpaths = reference_gpkg, divides = reference_gpkg,
      hydrolocations = hydrolocations_gpkg, pois = hydrolocations_gpkg)
    # flowpaths and divides are streamed one at a time, they are the large
    # layers, and holding them alongside each other buys nothing.
    for (lyr in c("flowpaths", "divides")) {
      x <- read_slice(srcs[[lyr]], lyr)
      if (is.null(x)) next
      # Same Stage 0 normalization as presplit_reference(): keep the two
      # entry points identical so Stage 1 sees the same geometry either way.
      x <- if (lyr == "flowpaths") .normalize_flowpaths(x) else .normalize_divides(x)
      sf::st_write(x, out_slice, lyr, quiet = TRUE, delete_layer = TRUE)
      rm(x)
    }
    # The two POI layers share one key, so they are renumbered together and are
    # the only pair that has to be in hand at the same time.
    poi_lyrs <- stats::setNames(
      lapply(c("hydrolocations", "pois"), function(l) read_slice(srcs[[l]], l)),
      c("hydrolocations", "pois"))
    if (length(exclude_hl_types)) poi_lyrs <- .drop_hl_types(poi_lyrs,
      exclude_hl_types, vpuid)
    poi_lyrs <- .normalize_poi_id(poi_lyrs)
    for (lyr in names(poi_lyrs)) {
      if (is.null(poi_lyrs[[lyr]])) next
      sf::st_write(poi_lyrs[[lyr]], out_slice, lyr, quiet = TRUE,
        delete_layer = TRUE)
    }
  })
  out_slice
}

#' First/last vertex of every LINESTRING in an sfc
#' @return list(sx, sy, ex, ey) numeric vectors, one entry per feature.
#' @keywords internal
#' @noRd
.line_ends <- function(g) {
  # Align endpoints to feature row: st_coordinates() OMITS empty geometries, so a
  # naive co[fr, ] would shift every row after an empty feature and shorten the
  # result. Index back through the feature column and leave empties as NA (the
  # caller already guards NA endpoints).
  n  <- length(g)
  co <- sf::st_coordinates(g)
  sx <- sy <- ex <- ey <- rep(NA_real_, n)
  if (nrow(co)) {
    fi <- co[, ncol(co)]                            # last index col = feature row
    fr <- !duplicated(fi)
    lr <- !duplicated(fi, fromLast = TRUE)
    rows_f <- as.integer(fi[fr]); rows_l <- as.integer(fi[lr])
    sx[rows_f] <- co[fr, "X"]; sy[rows_f] <- co[fr, "Y"]
    ex[rows_l] <- co[lr, "X"]; ey[rows_l] <- co[lr, "Y"]
  }
  list(sx = sx, sy = sy, ex = ex, ey = ey)
}

#' Normalize reference flowpaths: single LINESTRING + topology-driven outlet-last
#'
#' @description
#' Stage 0 geometry normalization for the `flowpaths` layer:
#'   1. **Single-part**: `st_line_merge()` joins contiguous parts of any
#'      MULTILINESTRING; a genuinely-disjoint remainder is reduced to its
#'      longest part (a flowpath is one line). Coordinates are otherwise
#'      untouched.
#'   2. **Orientation**, every line is stored mouth-last. Two oracles:
#'      * *fast path*: where `outlet_X/Y` coincides (within 1 m) with the
#'        start vertex but not the end, `st_reverse()`.
#'      * *topology fallback*, for lines whose `outlet_X/Y` sits on **neither**
#'        endpoint (e.g. `ls_targets` flowpaths whose outlet records a
#'        downstream junction, not a vertex), orient by connectivity: the
#'        endpoint that meets the downstream (`flowpath_toid`) flowpath within
#'        `connect_tol` metres is the mouth and must be the END. Lines where
#'        neither endpoint connects are left as-is and reported (no silent
#'        skip).
#'
#' @param fp An sf flowpaths object.
#' @param connect_tol Metres within which an endpoint is considered to meet the
#'   downstream flowpath. Default 50.
#' @return `fp` (single LINESTRING, mouth-last), row count preserved 1:1, with an
#'   `"orientation_report"` attribute summarising the topology fallback.
#' @keywords internal
#' @noRd
.normalize_flowpaths <- function(fp, connect_tol = 50) {
  if (nrow(fp) == 0L) return(fp)
  crs <- sf::st_crs(fp)

  # 1. Single LINESTRING per flowpath. Operate on a detached sfc and rebuild via
  # st_sfc() so the sfc class reflects the actual (single-part) types.
  geo   <- sf::st_geometry(fp)
  multi <- which(as.character(sf::st_geometry_type(geo)) == "MULTILINESTRING")
  if (length(multi)) {
    geo[multi] <- sf::st_line_merge(geo[multi])
    still <- which(as.character(sf::st_geometry_type(geo)) == "MULTILINESTRING")
    for (i in still) {
      parts    <- sf::st_cast(geo[i], "LINESTRING", warn = FALSE)
      geo[[i]] <- parts[[which.max(as.numeric(sf::st_length(parts)))]]
    }
    if (length(still))
      message(sprintf("  [geom-fix] reduced %d disjoint multilinestring(s) to longest part",
        length(still)))
    sf::st_geometry(fp) <- sf::st_sfc(geo, crs = crs)
  }
  fp <- sf::st_cast(fp, "LINESTRING", warn = FALSE)

  if (!inherits(sf::st_geometry(fp), "sfc_LINESTRING")) return(fp)

  # 2a. Fast orientation via outlet_X/Y (the on-line case: ~all NHD flowpaths).
  has_outlet <- all(c("outlet_X", "outlet_Y") %in% names(fp))
  ep  <- .line_ends(sf::st_geometry(fp))
  tol <- 1
  on_start <- on_end <- rep(FALSE, nrow(fp))
  if (has_outlet) {
    on_start <- abs(ep$sx - fp$outlet_X) < tol & abs(ep$sy - fp$outlet_Y) < tol
    on_end   <- abs(ep$ex - fp$outlet_X) < tol & abs(ep$ey - fp$outlet_Y) < tol
    on_start[is.na(on_start)] <- FALSE
    on_end[is.na(on_end)] <- FALSE
    rev_outlet <- which(on_start & !on_end)
    if (length(rev_outlet)) {
      g <- sf::st_geometry(fp)
      g[rev_outlet] <- sf::st_reverse(g[rev_outlet])
      sf::st_geometry(fp) <- g
      ep <- .line_ends(sf::st_geometry(fp))          # refresh after reversal
      message(sprintf("  [geom-fix] reversed %d flowpath(s) to outlet-last (outlet_X/Y on start vertex)",
        length(rev_outlet)))
    }
  }

  # 2b. Topology fallback for lines whose outlet sits on neither endpoint.
  report <- list(reversed = integer(0), confirmed = integer(0), unresolved = integer(0))
  residual <- which(!on_start & !on_end)
  if (length(residual) && all(c("flowpath_id", "flowpath_toid") %in% names(fp))) {
    id2row <- stats::setNames(seq_len(nrow(fp)), as.character(fp$flowpath_id))
    geo    <- sf::st_geometry(fp)
    rev_topo <- integer(0)
    for (r in residual) {
      # empty/degenerate geometry (e.g. zero-length reach collapsed by line_merge)
      # has NA endpoints: cannot orient by connectivity, leave unresolved.
      if (is.na(ep$ex[r]) || is.na(ep$sx[r])) {
        report$unresolved <- c(report$unresolved, r)
        next
      }
      dn <- id2row[as.character(fp$flowpath_toid[r])]
      if (is.na(dn)) {
        report$unresolved <- c(report$unresolved, r)
        next
      }
      pe <- sf::st_sfc(sf::st_point(c(ep$ex[r], ep$ey[r])), crs = crs)
      ps <- sf::st_sfc(sf::st_point(c(ep$sx[r], ep$sy[r])), crs = crs)
      d_end   <- as.numeric(sf::st_distance(pe, geo[dn]))
      d_start <- as.numeric(sf::st_distance(ps, geo[dn]))
      if (d_start <= connect_tol && d_end > connect_tol) {
        rev_topo <- c(rev_topo, r)                   # start is the mouth -> reverse
      } else if (d_end <= connect_tol) {
        report$confirmed <- c(report$confirmed, r)   # already mouth-last
      } else {
        report$unresolved <- c(report$unresolved, r) # neither endpoint connects
      }
    }
    if (length(rev_topo)) {
      geo[rev_topo] <- sf::st_reverse(geo[rev_topo])
      sf::st_geometry(fp) <- geo
      report$reversed <- rev_topo
      message(sprintf("  [geom-fix] reversed %d flowpath(s) by downstream connectivity (outlet_X/Y off-line)",
        length(rev_topo)))
    }
    if (length(report$unresolved))
      message(sprintf(paste0(
        "  [geom-fix] %d flowpath(s) could not be oriented (outlet off-line, ",
        "no downstream connection within %dm): %s"),
        length(report$unresolved), connect_tol,
        paste(utils::head(fp$flowpath_id[report$unresolved], 20L), collapse = ", ")))
  }

  attr(fp, "orientation_report") <- lapply(report, function(ix) as.character(fp$flowpath_id[ix]))
  fp
}

#' Normalize reference divides: single POLYGON per divide, slivers merged in
#'
#' @description
#' Stage 0 geometry normalization for the `divides` layer. For each multipart
#' divide it keeps the **largest part** as the divide's polygon and re-attaches
#' the remaining parts ("slivers") rather than dropping them:
#'   * a sliver is dissolved into the neighbouring polygon with which it shares
#'     the **longest boundary** (area is conserved, it moves to that neighbour);
#'   * a sliver that touches nothing along a positive-length edge (a genuinely
#'     disjoint catchment, e.g. an island) is kept attached to its **own**
#'     divide, which therefore remains MULTIPOLYGON.
#'
#' Total divide area is conserved exactly (a positional
#' `full[multi, ] <- st_cast(...)` write would mis-align the geometry column on
#' row explosion and silently lose area).
#'
#' @param dv An sf divides object.
#' @param id_col Unique per-divide identifier column. Default `"divide_id"`.
#' @return `dv` with one (largest) polygon per divide, row count preserved 1:1.
#' @keywords internal
#' @noRd
.normalize_divides <- function(dv, id_col = "divide_id") {
  if (nrow(dv) == 0L) return(dv)
  crs <- sf::st_crs(dv)

  # Part counts (cheap: only inspect MULTIPOLYGON rows).
  gt     <- as.character(sf::st_geometry_type(dv))
  mp     <- which(gt == "MULTIPOLYGON")
  nparts <- rep(1L, nrow(dv))
  if (length(mp)) nparts[mp] <- vapply(sf::st_geometry(dv)[mp], length, integer(1L))

  multi_idx <- which(nparts > 1L)
  if (length(multi_idx) == 0L) {
    # All single-part: lossless 1:1 retype (single MULTIPOLYGON -> POLYGON).
    return(sf::st_cast(dv, "POLYGON", warn = FALSE))
  }
  if (anyDuplicated(dv[[id_col]]))
    stop(sprintf(".normalize_divides: '%s' must be unique per divide", id_col))

  orig_geom <- sf::st_geometry(dv)

  # Reduce each divide to its largest part; collect non-largest parts as slivers.
  prim      <- sf::st_geometry(dv)               # list of sfg, edited in place
  sliv_list <- list()
  sliv_src  <- integer(0)                        # donor ROW index for each sliver
  for (i in mp) {
    g <- prim[[i]]
    # unwrap single-part MULTIPOLYGON -> POLYGON sfg (g[[1]] is a bare ring list)
    if (nparts[i] == 1L) {
      prim[[i]] <- sf::st_polygon(g[[1L]])
      next
    }
    parts <- sf::st_cast(sf::st_sfc(g, crs = crs), "POLYGON", warn = FALSE)
    k     <- which.max(as.numeric(sf::st_area(parts)))
    prim[[i]] <- parts[[k]]
    for (j in seq_along(parts)[-k]) {
      sliv_list[[length(sliv_list) + 1L]] <- parts[[j]]
      sliv_src <- c(sliv_src, i)
    }
  }
  prim_sfc <- sf::st_sfc(prim, crs = crs)

  # Assign each sliver to the primary with the longest shared boundary.
  n_merge <- 0L
  n_iso <- 0L
  dest <- integer(0)
  if (length(sliv_list)) {
    sliv_sfc  <- sf::st_sfc(sliv_list, crs = crs)
    hits      <- suppressWarnings(sf::st_intersects(sliv_sfc, prim_sfc))
    prim_bnd  <- sf::st_boundary(prim_sfc)
    sliv_bnd  <- sf::st_boundary(sliv_sfc)
    target    <- integer(length(sliv_list))      # primary row index, 0 = isolated
    for (s in seq_along(sliv_list)) {
      cand <- hits[[s]]
      cand <- cand[cand != sliv_src[s]]           # ignore own primary (disjoint -> 0 edge)
      if (!length(cand)) {
        target[s] <- 0L
        next
      }
      shared <- suppressWarnings(sf::st_intersection(prim_bnd[cand], sliv_bnd[s]))
      ln     <- suppressWarnings(as.numeric(sf::st_length(shared)))
      ln[is.na(ln)] <- 0
      if (length(ln) == 0L || max(ln) <= 0) {
        target[s] <- 0L
        next
      }
      target[s] <- cand[which.max(ln)]
    }
    # Merge: edge-sharing slivers -> neighbour; isolated slivers -> own divide.
    dest    <- ifelse(target > 0L, target, sliv_src)
    n_merge <- sum(target > 0L)
    n_iso <- sum(target == 0L)
    for (t in unique(dest)) {
      sset      <- which(dest == t)
      merged    <- sf::st_union(sf::st_sfc(c(list(prim[[t]]), sliv_list[sset]), crs = crs))
      prim[[t]] <- sf::st_geometry(merged)[[1L]]
    }
  }

  sf::st_geometry(dv) <- sf::st_sfc(prim, crs = crs)

  # Area conservation over every touched row (donors + merge targets): nothing
  # is dropped, so the delta must be ~0 to floating-point.
  touched  <- union(multi_idx, unique(dest))
  a_before <- sum(as.numeric(sf::st_area(orig_geom[touched])))
  a_after  <- sum(as.numeric(sf::st_area(sf::st_geometry(dv)[touched])))
  d_km2    <- (a_after - a_before) / 1e6
  if (abs(d_km2) > 1e-3)
    .hf_log_warn(sprintf(".normalize_divides: area delta %.6f km2 over %d touched divide(s)",
      d_km2, length(touched)))
  message(sprintf(
    "  [geom-fix] %d multipart divide(s): merged %d sliver(s) into longest-edge neighbour, %d kept (isolated); area delta %.6f km2",
    length(multi_idx), n_merge, n_iso, d_km2))
  dv
}

#' Pre-split a multi-VPU reference fabric into per-VPU slices
#'
#' @description
#' Reads each layer of `reference_gpkg` **once**, groups rows by `vpuid`, and
#' writes one `ref_vpu{id}.gpkg` per VPU into `<out_dir>/vpu_{id}/work/`.
#' When these files exist, refactor Step 0 is an instant skip: eliminating the
#' per-VPU random read against the large source file.
#'
#' @param reference_gpkg Path to the multi-VPU reference GeoPackage.
#' @param out_dir Parent output directory (same as `run_pipeline`).
#' @param vpuids Character vector of VPU IDs to split out. `NULL` (default)
#'   splits all VPUs found in the `flowpaths` layer.
#' @param layers Layers to extract. Default matches refactor Step 0.
#' @param overwrite If `TRUE`, overwrite existing per-VPU slices. Default `FALSE`.
#' @return Invisibly, a named character vector `vpuid -> path` for every file written.
#' @examples
#' \dontrun{
#' presplit_reference(reference_gpkg = "reference_conus.gpkg",
#'                    out_dir        = "out",
#'                    vpuids         = c("01", "02"))
#' }
#' @export
presplit_reference <- function(reference_gpkg,
                               out_dir,
                               vpuids   = NULL,
                               layers   = c("flowpaths", "divides",
                                 "hydrolocations", "pois"),
                               overwrite = FALSE) {
  stopifnot(file.exists(reference_gpkg), is.character(out_dir), nzchar(out_dir))

  # Discover which VPUs are present
  con <- DBI::dbConnect(RSQLite::SQLite(), reference_gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  cols <- DBI::dbListFields(con, "flowpaths")
  if (!"vpuid" %in% cols) {
    message("presplit_reference: no 'vpuid' column in flowpaths -- nothing to split")
    return(invisible(character(0L)))
  }
  all_vpus <- DBI::dbGetQuery(con, "SELECT DISTINCT vpuid FROM flowpaths")$vpuid
  DBI::dbDisconnect(con)
  on.exit()

  target_vpus <- if (!is.null(vpuids)) intersect(vpuids, all_vpus) else all_vpus
  if (length(target_vpus) == 0L) {
    message("presplit_reference: no matching VPUs found")
    return(invisible(character(0L)))
  }

  # Determine which per-VPU files still need writing
  out_paths <- stats::setNames(
    file.path(out_dir, paste0("vpu_", target_vpus), "work",
      sprintf("ref_vpu%s.gpkg", target_vpus)),
    target_vpus)

  needs_write <- if (overwrite) target_vpus else
    target_vpus[!file.exists(out_paths[target_vpus])]

  if (length(needs_write) == 0L) {
    message("presplit_reference: all per-VPU slices already exist -- nothing to do")
    return(invisible(out_paths))
  }

  message(sprintf("presplit_reference: splitting %d VPU(s) from %s",
    length(needs_write), basename(reference_gpkg)))

  for (vpu in needs_write) {
    dir.create(dirname(out_paths[vpu]), recursive = TRUE, showWarnings = FALSE)
    if (overwrite && file.exists(out_paths[vpu])) file.remove(out_paths[vpu])
  }

  # Read each layer once, split by vpuid, write per-VPU files
  for (lyr in layers) {
    t0 <- proc.time()[["elapsed"]]
    full <- tryCatch(
      sf::st_read(reference_gpkg, lyr, quiet = TRUE,
        options = "GEOMETRY_COLUMN=geom"),
      error = function(e) {
        # Retry without geometry hint (handles layers whose geometry col IS named "geometry")
        tryCatch(sf::st_read(reference_gpkg, lyr, quiet = TRUE),
          error = function(e2) {
            message(sprintf("  layer '%s': %s -- skipped", lyr, conditionMessage(e2)))
            NULL
          })
      })
    if (is.null(full)) next

    if (!"vpuid" %in% names(full)) {
      message(sprintf("  layer '%s': no vpuid column -- skipped", lyr))
      next
    }

    full$vpuid <- as.character(full$vpuid)

    # Stage 0 geometry normalization, one routine, shared with the in-Stage-1
    # fallback (.extract_vpu_slice_if_needed). Flowpaths: outlet-last + single
    # LINESTRING. Divides: single POLYGON, slivers merged into longest-edge
    # neighbour (area-conserving). Both preserve row count 1:1.
    if (lyr == "flowpaths") full <- .normalize_flowpaths(full)
    if (lyr == "divides") full <- .normalize_divides(full)

    for (vpu in needs_write) {
      slice <- full[full$vpuid == vpu, ]
      sf::st_write(slice, out_paths[vpu], lyr, quiet = TRUE, delete_layer = TRUE)
    }
    elapsed <- proc.time()[["elapsed"]] - t0
    message(sprintf("  layer '%-20s': %6d rows -> %d files  [%.1fs]",
      lyr, nrow(full), length(needs_write), elapsed))
  }

  message(sprintf("presplit_reference: done -- %d VPU slices written to %s",
    length(needs_write), out_dir))
  invisible(out_paths[needs_write])
}

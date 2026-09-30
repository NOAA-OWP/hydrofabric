#!/usr/bin/env Rscript
# Build NWM v4 hydrofabric VPUs from the staged input set.
#
#   Rscript run_vpu.R 21                # one VPU (PRVI, the quickest check)
#   Rscript run_vpu.R 01 02 03N         # several, in series
#   Rscript run_vpu.R all               # every VPU, then the domain merge
#
# Run from a local copy of the staged inputs, or point HF_INPUT_ROOT at it.
# Output goes to HF_OUT_DIR (default ./nwmv4).
#
# Every path comes from input_manifest.json, and every input is checked against
# its recorded MD5 before anything is built: a truncated or half-synced file
# would otherwise surface hours later as a geometry error.

suppressPackageStartupMessages({
  library(hydrofabric)
  library(jsonlite)
})

ROOT <- normalizePath(Sys.getenv("HF_INPUT_ROOT", getwd()))
OUT  <- normalizePath(Sys.getenv("HF_OUT_DIR", file.path(getwd(), "nwmv4")),
                      mustWork = FALSE)

# Record the commit that built the products, taken from the installed package
# when the environment does not set it.
if (!nzchar(Sys.getenv("HF_GIT_SHA"))) {
  sha <- utils::packageDescription("hydrofabric")$RemoteSha
  if (!is.null(sha)) Sys.setenv(HF_GIT_SHA = substr(sha, 1, 7))
}

MANIFEST <- file.path(ROOT, "input_manifest.json")
if (!file.exists(MANIFEST))
  stop("manifest not found: ", MANIFEST,
       "\n  run from the staged inputs, or set HF_INPUT_ROOT", call. = FALSE)

# Manifest paths are relative to the directory holding it.
setwd(ROOT)
hf_verify_manifest(MANIFEST, strict = TRUE)
man <- fromJSON(MANIFEST)$inputs

# A VRT is only an index; the tiles it names must be present beside it.
tiles <- unlist(lapply(man$requires, function(r) if (length(r)) r$path))
missing <- tiles[!file.exists(tiles)]
if (length(missing))
  stop(sprintf("%d grid tile(s) missing, first: %s", length(missing),
               missing[[1]]), call. = FALSE)
message(sprintf("[inputs] %d role(s) and %d tile(s) verified", nrow(man),
                length(tiles)))

path <- function(role) {
  p <- man$path[man$role == role]
  if (!length(p)) stop("manifest has no role: ", role, call. = FALSE)
  normalizePath(p)
}

CONUS <- c("01","02","03N","03S","03W","04","05","06","07","08","09",
           "10L","10U","11","12","13","14","15","16","17","18")
ALL   <- c(CONUS, "19", "20", "21")
args  <- commandArgs(trailingOnly = TRUE)
vpus  <- if (!length(args) || identical(args, "all")) ALL else args
unknown <- setdiff(vpus, ALL)
if (length(unknown)) stop("unknown VPU(s): ", paste(unknown, collapse = ", "),
                          call. = FALSE)

# Inputs per domain. Flowpath slope is read from the channel table, which
# carries it. Only CONUS has a fallback grid.
domain_inputs <- function(v) {
  d <- if (v == "19") "ak" else if (v == "20") "hi" else if (v == "21") "prvi" else "conus"
  sfx <- if (d == "conus") "" else paste0("_", d)
  list(
    reference_gpkg      = path(paste0("reference_", d)),
    hydrolocations_gpkg = path("hydrolocations"),
    fdr_path            = path(paste0("fdr_", d)),
    fac_path            = path(paste0("fac_", d)),
    fallback_fdr_path   = if (d == "conus") path("fdr_conus_fallback"),
    fallback_fac_path   = if (d == "conus") path("fac_conus_fallback"),
    hydraulics_path     = path(paste0("hydraulics", sfx)),
    lakes_path          = path(paste0("lakes", sfx)),
    lake_mask_path      = path(paste0("lake_mask", sfx)))
}

message(sprintf("[run] %d VPU(s): %s", length(vpus), paste(vpus, collapse = ", ")))
message(sprintf("[run] output to %s", OUT))

res <- do.call(rbind, lapply(vpus, function(v) {
  message(sprintf("\n[run] VPU %s", v))
  t0 <- proc.time()[["elapsed"]]
  ok <- tryCatch({
    do.call(run_pipeline, c(list(vpuid = v, out_dir = OUT), domain_inputs(v)))
    TRUE
  }, error = function(e) {
    message("[run] VPU ", v, " FAILED: ", conditionMessage(e)); FALSE
  })
  data.frame(vpuid = v, ok = ok,
             minutes = round((proc.time()[["elapsed"]] - t0) / 60, 1))
}))
print(res, row.names = FALSE)
if (!all(res$ok)) quit(status = 1L)

# The merge assigns one id pool across every VPU it is given, so it runs only
# when all 24 were built in this call.
if (setequal(vpus, ALL)) {
  gpkg <- file.path(OUT, paste0("vpu_", vpus), paste0("ngen_vpu", vpus, ".gpkg"))
  names(gpkg) <- vpus
  merge_vpus(gpkg_map = as.list(gpkg), out_dir = file.path(OUT, "global"),
             overwrite = TRUE, strict = TRUE)
}

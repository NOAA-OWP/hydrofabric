# =====================================================================
# report_layer_effect
# =====================================================================
# Effect reporting for steps that mutate a divides layer in place.
#
# Several repair steps: notably the two hfutils::clean_gpkg_layer calls at
# stage-2 Step 4b and stage-3 Step 6d: logged only their run time, so their
# effect on the fabric could not be measured from a run. A step whose effect is
# unknown can neither be trusted nor removed, which is the whole reason these
# helpers exist.
#
# Three numbers are reported, and the third is the one that matters. A row count
# says whether polygons disappeared. A summed area cannot distinguish a hole from
# an overlap being removed, because both reduce it. Only the union (the footprint)
# separates them: overlap removal leaves the union unchanged, a hole shrinks it.

# Summarise a divides layer. `union` is optional because the union of tens of
# thousands of polygons is the expensive part; skip it where only the cheap
# signal is wanted.
.divide_metrics <- function(x, union = TRUE) {
  out <- list(n = NA_integer_, sum_m2 = NA_real_, union_m2 = NA_real_)
  if (is.null(x) || !inherits(x, "sf") || nrow(x) == 0L) return(out)
  net <- tryCatch(.hf_network_divides(x), error = function(e) NULL)
  if (is.null(net) || nrow(net) == 0L) return(out)
  out$n      <- nrow(net)
  out$sum_m2 <- sum(suppressWarnings(as.numeric(sf::st_area(net))), na.rm = TRUE)
  if (isTRUE(union)) {
    # Via geos rather than sf::st_union. On the largest VPU (36,626 divides) the
    # sf path is heavy enough that adding it at two call sites killed the run
    # outright, the process died with no R error after 3.6 hours, the signature
    # of an out-of-memory kill. geos keeps the whole thing in one C-level
    # collection and returns the same number.
    out$union_m2 <- tryCatch({
      g <- geos::as_geos_geometry(sf::st_geometry(net))
      bad <- !geos::geos_is_valid(g); bad[is.na(bad)] <- TRUE
      if (any(bad)) g[bad] <- geos::geos_make_valid(g[bad])
      u <- geos::geos_unary_union(geos::geos_make_collection(g))
      as.numeric(geos::geos_area(u))
    }, error = function(e) NA_real_)
  }
  out
}

# Emit the before/after delta for one step. Silent when nothing measurable
# changed, so a no-op step stays quiet in the log rather than adding noise,
# but it says so explicitly, because "no line" and "no effect" must not look the
# same to whoever reads the log later.
.report_divide_effect <- function(label, before, after, vpuid = NA_character_) {
  tag <- if (is.na(vpuid)) "" else sprintf("[%s] ", vpuid)
  if (is.na(before$n) || is.na(after$n)) {
    message(sprintf("%s%s: effect not measurable", tag, label))
    return(invisible(NULL))
  }
  d_n   <- after$n - before$n
  d_sum <- (after$sum_m2 - before$sum_m2) / 1e6
  d_uni <- if (is.finite(before$union_m2) && is.finite(after$union_m2))
    (after$union_m2 - before$union_m2) / 1e6 else NA_real_

  # Overlap is sum - union; a drop in it is the intended work of a clean step.
  ov_before <- if (is.finite(before$union_m2))
    max(0, before$sum_m2 - before$union_m2) / 1e6 else NA_real_
  ov_after <- if (is.finite(after$union_m2))
    max(0, after$sum_m2 - after$union_m2) / 1e6 else NA_real_

  nothing <- d_n == 0L && abs(d_sum) < 1e-6 &&
    (is.na(d_uni) || abs(d_uni) < 1e-6)
  if (nothing) {
    message(sprintf("%s%s: no effect (%d divides unchanged)", tag, label, after$n))
    return(invisible(NULL))
  }
  message(sprintf(
    "%s%s: divides %+d (%d -> %d) | area %+.4f km2 | coverage %s | overlap %s",
    tag, label, d_n, before$n, after$n, d_sum,
    if (is.na(d_uni)) "n/a" else sprintf("%+.4f km2", d_uni),
    if (is.na(ov_before) || is.na(ov_after)) "n/a" else
      sprintf("%.4f -> %.4f km2", ov_before, ov_after)))
  invisible(NULL)
}

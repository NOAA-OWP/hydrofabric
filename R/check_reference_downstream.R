#' Find reference reaches whose downstream neighbour is missing
#'
#' @description
#' A reach with no downstream is an outlet: it reaches the sea, or a closed
#' basin, and the network legitimately stops. A reach that *names* a downstream
#' the file does not contain is something else. The network was truncated, and
#' the reach now reads as an outlet to anything that consumes it.
#'
#' Routing is where it shows. Every reach with no downstream becomes a terminal
#' nexus, and flow arriving there is not routed on, so a truncated network
#' underreports flow at each cut. The pipeline cannot repair it: there is no
#' downstream reach to route into.
#'
#' The reference carries the evidence itself, in the value-added attributes it
#' inherits. `flowpath_toid` says where a reach drains and `dnhydroseq` names the
#' same neighbour by its `hydroseq`. A true outlet has
#' neither. A truncated reach still carries `dnhydroseq`, pointing at a reach
#' that was removed, and that dangling pointer is what this reports.
#'
#' Observed on the CONUS reference at the Canadian border: of VPU 04's 2,999
#' terminal nexuses, 2,948 trace to a reach naming a downstream by `dnhydroseq`
#' that exists nowhere in the file. VPU 01, whose terminals are mostly the New
#' England coastline, has 9. The `terminalfl` column cannot be used to tell the
#' two apart, because it is unpopulated for `NWM_v3_0_CONUS_Prelim_ls_mod`: 14
#' reaches of that source's 109,850 carry the flag.
#'
#' Report-only, and deliberately not wired into `run_pipeline()`. The coverage
#' belongs to the reference and is not ours to extend, and this asks a question
#' about a reference release rather than about a VPU.
#'
#' @param reference_gpkg Path to the reference GeoPackage.
#' @param layer Flowline layer name.
#' @param id,toid,dn,seq Column names for the identifier, the downstream
#'   identifier, the downstream hydrologic sequence, and the hydrologic sequence.
#' @param vpuid Optional VPU to restrict the *reported* reaches to. The search
#'   for the named downstream always spans the whole file, since a reach may
#'   legitimately drain into the next VPU.
#' @param strict If `TRUE`, warn when any are found.
#' @returns Invisibly, a data frame of `id`, `dangling` (which pointer dangles)
#'   and `names_id`, the neighbour that is missing.
#' @export
check_reference_downstream_resolves <- function(reference_gpkg,
                                                layer = "flowpaths",
                                                id    = "flowpath_id",
                                                toid  = "flowpath_toid",
                                                dn    = "dnhydroseq",
                                                seq   = "hydroseq",
                                                vpuid = NA_character_,
                                                strict = FALSE) {
  stopifnot(is.character(reference_gpkg), length(reference_gpkg) == 1L,
    file.exists(reference_gpkg))

  cols <- tryCatch(names(sf::st_read(reference_gpkg,
    query = sprintf('SELECT * FROM "%s" LIMIT 0', layer), quiet = TRUE)),
    error = function(e) NULL)
  if (is.null(cols)) {
    message("check_reference_downstream_resolves: layer not readable; skipped")
    return(invisible(data.frame()))
  }
  want <- c(id, toid, dn, seq)
  have <- want[want %in% cols]
  if (!all(c(id, toid) %in% have)) {
    message("check_reference_downstream_resolves: no id/toid columns; skipped")
    return(invisible(data.frame()))
  }
  sel <- paste(sprintf('"%s"', have), collapse = ", ")
  vp  <- if (!is.na(vpuid) && "vpuid" %in% cols) ', "vpuid"' else ""

  # The whole file, because a reach may legitimately drain into the next VPU and
  # restricting the search would report every such reach as truncated.
  x <- tryCatch(sf::st_read(reference_gpkg,
    query = sprintf('SELECT %s%s FROM "%s"', sel, vp, layer), quiet = TRUE),
    error = function(e) NULL)
  if (is.null(x) || !nrow(x)) {
    message("check_reference_downstream_resolves: no features read; skipped")
    return(invisible(data.frame()))
  }
  x <- sf::st_drop_geometry(x)

  ids  <- as.character(x[[id]])
  tos  <- as.character(x[[toid]])
  keep <- if (!is.na(vpuid) && "vpuid" %in% names(x))
    as.character(x$vpuid) == vpuid else rep(TRUE, nrow(x))
  keep[is.na(keep)] <- FALSE

  # (1) toid names an id that is not in the file.
  live_to <- !is.na(tos) & nzchar(tos) & tos != "0"
  bad_to  <- which(keep & live_to & !(tos %in% ids))

  # (2) toid says outlet, but dnhydroseq still names a neighbour that is gone.
  #
  # Only when dnhydroseq can hold the values it has to point at. On the CONUS
  # reference it cannot: hydroseq reaches 13,226,027,766 and 2,609,771 of
  # 2,725,867 reaches exceed the 32-bit limit, while dnhydroseq is stored as a
  # 32-bit integer, so 2,149,436 of them have wrapped to a negative number. Every
  # wrapped pointer then "fails to resolve" for a reason that has nothing to do
  # with coverage, and counting them would report an overflow as a truncated
  # network. A negative value is the tell, since no hydroseq is negative.
  bad_dn <- integer(0); dn_usable <- FALSE
  if (all(c(dn, seq) %in% names(x))) {
    dnv <- suppressWarnings(as.numeric(x[[dn]]))
    sqv <- suppressWarnings(as.numeric(x[[seq]]))
    n_neg <- sum(dnv < 0, na.rm = TRUE)
    n_wide <- sum(sqv > 2147483647, na.rm = TRUE)
    dn_usable <- n_neg == 0L && n_wide == 0L
    if (dn_usable) {
      bad_dn <- which(keep & !live_to & is.finite(dnv) & dnv != 0 & !(dnv %in% sqv))
    } else {
      message(sprintf(paste0("check_reference_downstream_resolves: %s is unusable ",
        "on this input -- %d value(s) negative and %d %s value(s) beyond the ",
        "32-bit limit it is stored in, so the pointers have overflowed; ",
        "reporting %s only"), dn, n_neg, n_wide, seq, toid))
    }
  }

  # Each block is built only when it has rows: a scalar `dangling` cannot be
  # recycled against zero rows.
  out <- rbind(
    if (length(bad_to)) data.frame(id = ids[bad_to], dangling = toid,
      names_id = tos[bad_to], stringsAsFactors = FALSE) else NULL,
    if (length(bad_dn)) data.frame(id = ids[bad_dn], dangling = dn,
      names_id = as.character(suppressWarnings(as.numeric(x[[dn]])[bad_dn])),
      stringsAsFactors = FALSE) else NULL)
  if (is.null(out)) out <- data.frame(id = character(0), dangling = character(0),
    names_id = character(0), stringsAsFactors = FALSE)

  n_tested <- sum(keep)
  message(sprintf(paste0("check_reference_downstream_resolves: %d of %d reach(es) ",
    "name a downstream the file does not contain (%d via %s%s)"),
    nrow(out), n_tested, length(bad_to), toid,
    if (dn_usable) sprintf(", %d via %s", length(bad_dn), dn) else
      sprintf(", %s not usable", dn)))
  if (nrow(out))
    message(sprintf("  examples: %s",
      paste(sprintf("%s -> %s (%s)", utils::head(out$id, 3),
        utils::head(out$names_id, 3), utils::head(out$dangling, 3)),
        collapse = ", ")))
  if (strict && nrow(out))
    warning("reference reaches name a downstream that is absent", call. = FALSE)

  invisible(out)
}

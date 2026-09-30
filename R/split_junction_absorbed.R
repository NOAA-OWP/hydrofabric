# =====================================================================
# split_junction_absorbed
# =====================================================================
# Keeps the routing network closed under toid: a reach receiving flow from another
# flowpath must lie on its own flowpath's mainstem, or it carries no line and
# water arrives with nowhere to go. See "Where the routing network lives" in
# the package documentation.
#
# The condition is the RECEIVING reach being off-mainstem, not the inflow being
# interior. Interior inflow onto an on-mainstem reach is ordinary aggregation and
# is left alone.
#
# Runs in stage 2 after membership is final and BEFORE construct_agg_divides, so
# divides are dissolved from the corrected member partition exactly once and area
# is conserved by construction.
#
# Each junction is an endpoint join at the receiving reach's inlet, so the split
# falls on a boundary that already exists between two member reaches. No reach
# geometry is cut, and a promoted chain is a simple path, so its merged line is
# single-part.

# Find chains of absorbed reaches that receive routed flow from another aggregate.
# `mem` is a list of member rec_ids per aggregate id; `lp`/`to` are reach-level
# levelpath and downstream maps. Returns a list of list(agg, recs) most-upstream
# first.
.junction_absorbed_chains <- function(mem, agg_ms, lp, to) {
  agg_of <- unlist(lapply(names(mem), function(a)
    stats::setNames(rep(a, length(mem[[a]])), mem[[a]])), use.names = TRUE)
  recs   <- names(agg_of)
  on_ms  <- !is.na(agg_ms[agg_of[recs]]) & lp[recs] == unname(agg_ms[agg_of[recs]])
  names(on_ms) <- recs

  dn      <- unname(to[recs]); names(dn) <- recs
  has_dn  <- !is.na(dn) & dn %in% recs
  ext <- recs[has_dn &
    agg_of[dn[recs]] != agg_of[recs] &
    !on_ms[dn[recs]] &
    on_ms[recs]]
  ext <- ext[!is.na(ext)]
  if (!length(ext)) return(list())

  seen <- character(0); out <- list()
  for (u in ext) {
    b <- dn[[u]]; a <- agg_of[[b]]
    ch <- character(0); cur <- b
    repeat {
      if (is.na(cur) || !(cur %in% recs)) break
      if (on_ms[[cur]] || agg_of[[cur]] != a) break
      ch <- c(ch, cur)
      nx <- dn[[cur]]
      if (is.na(nx)) break
      cur <- nx
      if (length(ch) > 500L) break
    }
    if (!length(ch)) next
    k <- paste(sort(ch), collapse = ",")
    if (k %in% seen) next
    seen <- c(seen, k)
    out[[length(out) + 1L]] <- list(agg = a, recs = ch)
  }
  out
}

#' Split absorbed chains that receive routed flow into their own flowpaths
#'
#' @param gpkg Stage-2 GeoPackage holding `flowpaths` and `network`, before
#'   `construct_agg_divides` has run.
#' @param s1_fp Reconciled reaches (`sf`) carrying `flowpath_id`, `flowpath_toid`,
#'   `levelpathid` and geometry.
#' @param vpuid VPU label for messages.
#' @returns Invisibly, the number of chains split out.
#' @keywords internal
.split_junction_absorbed <- function(gpkg, s1_fp, vpuid = NA_character_) {
  fp  <- sf::st_read(gpkg, "flowpaths", quiet = TRUE)
  net <- as.data.frame(sf::st_read(gpkg, "network", quiet = TRUE))
  if (!all(c("flowpath_id", "rec_id_list") %in% names(net))) {
    message(sprintf("promote[%s]: no rec_id_list; skipped", vpuid))
    return(invisible(0L))
  }
  fp$flowpath_id   <- as.character(fp$flowpath_id)
  fp$flowpath_toid <- as.character(fp$flowpath_toid)

  nu  <- unique(data.frame(
    fp  = as.character(net$flowpath_id),
    rl  = as.character(net$rec_id_list), stringsAsFactors = FALSE))
  mem <- stats::setNames(lapply(nu$rl, function(x) expand_member_csv(x, unique = TRUE)), nu$fp)
  mem <- mem[lengths(mem) > 0L]
  if (!length(mem)) return(invisible(0L))

  r_id <- as.character(s1_fp$flowpath_id)
  lp   <- stats::setNames(as.character(s1_fp$levelpathid), r_id)
  to   <- stats::setNames(as.character(s1_fp$flowpath_toid), r_id)
  ms   <- stats::setNames(as.character(fp$mainstem), fp$flowpath_id)

  chains <- .junction_absorbed_chains(mem, ms, lp, to)
  if (!length(chains)) {
    message(sprintf("promote[%s]: no absorbed chains receive routed flow", vpuid))
    return(invisible(0L))
  }

  next_id <- max(suppressWarnings(as.integer(fp$flowpath_id)), na.rm = TRUE)
  r_ix    <- stats::setNames(seq_len(nrow(s1_fp)), r_id)
  new_fp  <- vector("list", length(chains))
  net$rec_id_list <- as.character(net$rec_id_list)
  net$flowpath_id <- as.character(net$flowpath_id)

  for (k in seq_along(chains)) {
    ch      <- chains[[k]]
    orig    <- ch$agg
    next_id <- next_id + 1L
    nid     <- as.character(next_id)

    # Membership moves. construct_agg_divides dissolves rec_div by these lists, so
    # the divide split follows exactly: nothing is cut, created or destroyed.
    keep <- setdiff(mem[[orig]], ch$recs)
    net$rec_id_list[net$flowpath_id == orig] <- compose_member_csv(keep)
    mem[[orig]] <- keep

    row <- fp[match(orig, fp$flowpath_id), , drop = FALSE]
    row$flowpath_id   <- nid
    row$flowpath_toid <- orig                 # the chain drains into its parent
    row$mainstem      <- unname(lp[ch$recs[length(ch$recs)]])
    if ("divide_id"    %in% names(row)) row$divide_id <- nid
    if ("poi_id"       %in% names(row)) row$poi_id <- NA
    if ("member_comid" %in% names(row)) row$member_comid <- NA_character_
    if ("has_divide"   %in% names(row)) row$has_divide <- TRUE
    gi <- unname(r_ix[ch$recs]); gi <- gi[!is.na(gi)]
    g  <- tryCatch(sf::st_line_merge(sf::st_union(sf::st_geometry(s1_fp)[gi])),
      error = function(e) sf::st_geometry(s1_fp)[gi][1])
    sf::st_geometry(row) <- sf::st_sfc(if (inherits(g, "sfc")) g[[1]] else g,
      crs = sf::st_crs(fp))
    new_fp[[k]] <- row

    nr <- net[net$flowpath_id == orig, , drop = FALSE][1, , drop = FALSE]
    nr$flowpath_id   <- nid
    nr$rec_id_list   <- compose_member_csv(ch$recs)
    if ("flowpath_toid" %in% names(nr)) nr$flowpath_toid <- orig
    net <- rbind(net, nr)

    # Whatever drained into the chain now drains into the new flowpath.
    up <- names(which(vapply(mem, function(m)
      any(!is.na(to[m]) & to[m] %in% ch$recs), logical(1L))))
    up <- setdiff(up, c(orig, nid))
    if (length(up)) {
      fp$flowpath_toid[fp$flowpath_id %in% up] <- nid
      if ("flowpath_toid" %in% names(net))
        net$flowpath_toid[net$flowpath_id %in% up] <- nid
    }
    mem[[nid]] <- ch$recs
  }

  fp <- rbind(fp, do.call(rbind, new_fp))
  if ("lengthkm" %in% names(fp))
    fp$lengthkm <- as.numeric(sf::st_length(sf::st_transform(fp, 5070))) / 1000

  sf::st_write(.drop_list_cols(fp), gpkg, "flowpaths", delete_layer = TRUE, quiet = TRUE)
  con <- DBI::dbConnect(RSQLite::SQLite(), gpkg)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbWriteTable(con, "network", net, overwrite = TRUE)

  message(sprintf(paste0(
    "promote[%s]: split %d absorbed chain(s) (%d reach(es)) into their own ",
    "flowpaths; flowpaths %d -> %d"),
    vpuid, length(chains), sum(vapply(chains, function(c) length(c$recs), integer(1))),
    nrow(fp) - length(chains), nrow(fp)))
  invisible(length(chains))
}

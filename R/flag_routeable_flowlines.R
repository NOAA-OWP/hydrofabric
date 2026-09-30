#' Flag the routeable flowlines (the routing network)
#'
#' @description
#' Adds a logical `routeable` attribute to the `flowlines` layer marking the
#' reaches in the routing network: the flowpath-constituent (coincident)
#' flowlines plus the incremental reaches that bridge a coincident flowline to
#' the next coincident flowline downstream. Routing follows the native
#' `flowline_toid`; the `routeable == TRUE` subgraph is closed under it (every
#' non-terminal routeable flowline drains to a routeable flowline), so a
#' consumer can route on the routeable set without a coincident reach draining
#' into a non-coincident connector that is filtered out ("flying water").
#'
#' A bridge reach is a real channel carrying flow between mainstems; it is
#' non-coincident only because aggregation traced a different mainstem. When a
#' coincident flowline never rejoins a coincident flowline downstream, the
#' bridge is extended to the terminal instead, so the reach's flow exits the
#' domain rather than draining into a filtered-out connector.
#'
#' Requires `flowpath_constituent` on the flowlines layer (see
#' [flag_flowpath_constituents()], which must run first).
#'
#' @param input_gpkg Path to a GeoPackage containing `flowlines_layer`.
#' @param flowlines_layer Layer name. Default `"flowlines"`.
#'
#' @return Invisibly, the number of routeable flowlines. Adds `routeable`
#'   (logical) to `flowlines_layer`, overwriting it in place.
#' @examples
#' \dontrun{
#' flag_routeable_flowlines("ngen_vpu01.gpkg")
#' }
#' @export
flag_routeable_flowlines <- function(input_gpkg, flowlines_layer = "flowlines") {
  fl <- sf::st_read(input_gpkg, flowlines_layer, quiet = TRUE)
  if (!nrow(fl)) return(invisible(0L))
  if (!("flowpath_constituent" %in% names(fl)))
    stop("flag_routeable_flowlines requires flowpath_constituent on the ",
      "flowlines layer; run flag_flowpath_constituents first", call. = FALSE)

  id   <- as.character(fl$flowline_id)
  toid <- as.character(fl$flowline_toid)
  coin <- as.logical(fl$flowpath_constituent); coin[is.na(coin)] <- FALSE
  idx  <- stats::setNames(seq_along(id), id)          # id -> row

  routeable <- coin                                   # coincident reaches route
  di   <- unname(idx[toid])                           # downstream row (NA = terminal)
  term <- is.na(di)
  dn_coin <- rep(FALSE, length(id)); dn_coin[!term] <- coin[di[!term]]

  # A break is a coincident reach whose immediate downstream is a non-coincident
  # reach. Bridge it: walk flowline_toid, marking intervening non-coincident
  # reaches routeable. The walk stops at the next coincident reach or at a
  # terminal; either way the intervening reaches are marked so the coincident
  # reach's flow reaches its downstream coincident reach or exits the domain.
  # Marking to the terminal too keeps the routeable set closed under toid (a
  # coincident outlet never drains into a filtered-out connector).
  brk <- which(coin & !term & !dn_coin)
  n_bridge <- 0L
  for (k in brk) {
    path <- integer(0); cur <- di[k]; seen <- integer(0)
    while (!is.na(cur) && !(cur %in% seen)) {
      if (coin[cur]) break            # reached the next coincident reach
      path <- c(path, cur); seen <- c(seen, cur); cur <- di[cur]
    }
    if (length(path)) { routeable[path] <- TRUE; n_bridge <- n_bridge + length(path) }
  }

  fl$routeable <- routeable
  sf::st_write(fl, input_gpkg, flowlines_layer, delete_layer = TRUE, quiet = TRUE)
  message(sprintf(
    "flag_routeable_flowlines: %d of %d flowline(s) routeable (%d bridge reach(es))",
    sum(routeable), length(id), n_bridge))
  invisible(sum(routeable))
}

# =====================================================================
# The POI contract
# =====================================================================
# Hydrolocations impose a five-step chain: they arrive with hl_reference as
# identity, snap to a location that mints one poi_id, split their reference reach
# where the location is interior, stay fixed and stay an OUTLET at every derived
# grain, and carry a catchment unless the reach is a reference connector.
#
# The chain, the outlet convention and what happens when a location cannot be
# split: vignette("poi-contract").
#
# Checks map onto the chain: pois_preserved, poi_ids_minted, poi_ids_not_reminted
# and poi_id_is_integer cover steps 1-2; poi_location_stable covers 4 against the
# location stage 1 established; poi_at_outlet and poi_at_flowline_outlet cover
# 4's outlet half at aggregated and reconciled grain; poi_outlet_has_divide
# covers 5.
#
# Measurement rules the checks depend on:
#
#   Distance is point-to-point, against the reach the POI is attributed to and
#   joined by id: never the nearest reach of any kind. The exception is
#   poi_at_flowline_outlet, where nearest IS the question.
#
#   `tol_m` is 30 m throughout: the DEM resolution, and this fabric's snap
#   tolerance.

# Position of each POI along the reach it is attributed to. Returns NULL when the
# inputs are absent so callers can stay unconditional.
.hf_poi_decompose <- function(pois, lines, line_id_col, poi_id_col = "flowpath_id",
                              tol_m = 30) {
  if (is.null(pois) || is.null(lines)) return(NULL)
  if (!inherits(pois, "sf") || !inherits(lines, "sf")) return(NULL)
  if (!nrow(pois) || !nrow(lines)) return(NULL)
  if (!poi_id_col %in% names(pois) || !line_id_col %in% names(lines)) return(NULL)

  tryCatch({
    if (sf::st_crs(pois) != sf::st_crs(lines))
      pois <- sf::st_transform(pois, sf::st_crs(lines))

    ix <- match(as.character(pois[[poi_id_col]]), as.character(lines[[line_id_col]]))
    lg <- sf::st_geometry(lines)
    pg <- sf::st_geometry(pois)

    n <- nrow(pois)
    d_out <- rep(NA_real_, n)
    d_in  <- rep(NA_real_, n)
    frac  <- rep(NA_real_, n)
    for (i in seq_len(n)) {
      if (is.na(ix[i])) next
      m <- sf::st_coordinates(lg[[ix[i]]])
      if (is.null(dim(m)) || nrow(m) < 2L) next
      # Cumulative distance along the stored vertex order. Reconciled and
      # aggregated lines carry enforced flow direction, so the last vertex is the
      # outlet and this is a real along-line measure rather than a straight line.
      cd <- cumsum(c(0, sqrt(diff(m[, 1])^2 + diff(m[, 2])^2)))
      tot <- cd[length(cd)]
      if (!is.finite(tot) || tot <= 0) next
      q <- sf::st_coordinates(pg[i])
      j <- which.min((m[, 1] - q[1])^2 + (m[, 2] - q[2])^2)
      # Point-to-point to the outlet vertex, not along-line from the projected
      # position. Projecting first discards any off-line displacement, so a POI
      # sitting 70 m to the side of its reach's end was scored as "at the outlet".
      # The assertion is that the POI COINCIDES with the outlet, which is a
      # distance between two points.
      d_out[i] <- sqrt((m[nrow(m), 1] - q[1])^2 + (m[nrow(m), 2] - q[2])^2)
      d_in[i]  <- sqrt((m[1, 1] - q[1])^2 + (m[1, 2] - q[2])^2)
      frac[i]  <- cd[j] / tot
    }

    matched <- !is.na(d_out)
    list(
      n         = n,
      unmatched = sum(!matched),
      at_outlet = sum(matched & d_out <= tol_m),
      at_inlet  = sum(matched & d_out > tol_m & d_in <= tol_m),
      interior  = sum(matched & d_out > tol_m & d_in > tol_m),
      max_dist  = if (any(matched)) max(d_out[matched]) else NA_real_,
      d_out     = d_out,
      d_in      = d_in,
      frac      = frac,
      idx       = ix
    )
  }, error = function(e) NULL)
}

# Requirement 2, reported per grain.
.hf_poi_position_info <- function(dec, grain, types = NULL, tol_m = 30) {
  if (is.null(dec) || dec$n == 0L) return(NULL)
  off <- dec$at_inlet + dec$interior + dec$unmatched

  # Name the types behind the residual. An at-inlet terminal and an interior gage
  # are different problems, and the type is what says which.
  tstr <- ""
  if (off > 0L && !is.null(types)) {
    bad <- which(is.na(dec$d_out) | (dec$d_out > tol_m))
    if (length(bad)) {
      tt <- sort(table(as.character(types[bad])), decreasing = TRUE)
      tt <- utils::head(tt, 4L)
      tstr <- sprintf(" [%s]", paste(sprintf("%s=%d", names(tt), as.integer(tt)),
        collapse = ", "))
    }
  }

  .hf_info(sprintf(paste0(
    "%s grain: %d of %d POI(s) at their reach's outlet (%.1f%%); ",
    "off: %d at an inlet, %d interior, %d unmatched; max %.0fm%s (target 100%% at outlet)"),
    grain, dec$at_outlet, dec$n, 100 * dec$at_outlet / dec$n,
    dec$at_inlet, dec$interior, dec$unmatched,
    if (is.finite(dec$max_dist)) dec$max_dist else 0, tstr))
}

# Requirement 1. Nothing may drop a POI between stage 1 and the product.
#
# Checked on hl_reference, not poi_id: a collocated pair legitimately delivers one
# poi_id and retires the other, so the retired id can never be found and every
# collapse would score as a loss. See "Collapse is not loss" in
# vignette("poi-contract").
#
# Two numbers, because they mean different things: references absent from the
# product entirely, and poi_ids retired by collocation whose references survive.
.hf_poi_preserved <- function(pois, pois_input, hydrolocations = NULL) {
  if (is.null(pois) || is.null(pois_input)) return(NULL)
  if (!"poi_id" %in% names(pois) || !"poi_id" %in% names(pois_input)) return(NULL)

  ini   <- as.character(pois_input$poi_id)
  ok_in <- !is.na(ini) & nzchar(ini)
  ini_u <- unique(ini[ok_in])
  if (!length(ini_u)) return(NULL)
  out  <- unique(as.character(pois$poi_id))
  gone <- setdiff(ini_u, out)

  .refs <- function(x) {
    if (is.null(x) || !"hl_reference" %in% names(x)) return(character(0))
    r <- unique(trimws(unlist(strsplit(as.character(x$hl_reference), "|",
      fixed = TRUE))))
    r[!is.na(r) & nzchar(r)]
  }
  prod_refs <- unique(c(.refs(pois), .refs(hydrolocations)))

  # A retired poi_id is explained when every hl_reference it carried is still in
  # the product. Anything else is a real loss.
  lost <- gone
  if (length(gone) && "hl_reference" %in% names(pois_input) && length(prod_refs)) {
    explained <- vapply(gone, function(p) {
      r <- .refs(pois_input[ok_in & ini == p, , drop = FALSE])
      length(r) > 0L && all(r %in% prod_refs)
    }, logical(1))
    lost <- gone[!explained]
  }

  # Two different reasons a poi_id can be explained-but-absent, and they are not
  # interchangeable: a POI withdrawn because the network cannot place it at an
  # outlet is a documented limitation, while a collocation collapse is the contract
  # working. The hydrolocation records the former in `poi_drop_reason`.
  withdrawn <- character(0)
  if (length(gone) && !is.null(hydrolocations) &&
      all(c("poi_id", "poi_drop_reason") %in% names(hydrolocations))) {
    wid <- as.character(hydrolocations$poi_id)[
      !is.na(hydrolocations$poi_drop_reason)]
    withdrawn <- intersect(setdiff(gone, lost), unique(wid))
  }
  collapsed <- length(gone) - length(lost) - length(withdrawn)

  types <- if (length(lost) && "hl_type" %in% names(pois_input)) {
    sel <- ok_in & ini %in% lost
    tt  <- sort(table(as.character(pois_input$hl_type)[sel]), decreasing = TRUE)
    sprintf(" [%s]", paste(sprintf("%s=%d", names(tt), as.integer(tt)),
      collapse = ", "))
  } else ""

  .hf_info(sprintf(paste0(
    "%d of %d stage-1 POI(s) lost -- hl_reference absent from the product%s; ",
    "%d withdrawn as unplaceable and %d retired by collocation, all with their ",
    "hydrolocations and references intact (target 0 lost)"),
    length(lost), length(ini_u), types, length(withdrawn), collapsed))
}

# Requirement 3. Having landed on an outlet, that reach must carry a divide,
# unless the reach never had one.
#
# A POI indexed to a connector outlet legitimately has no catchment. Connectors
# are the reference's zero-area waterbody and artificial paths: areasqkm = 0,
# has_divide = 0, no row in the reference divides layer, so no catchment exists
# for them without inventing area. Hydrolocation types index onto them by their
# nature: a waterbody outlet and a coastal coupling node sit on exactly those
# features.
#
# The assertion is therefore over catchment-bearing reaches, with connector POIs
# counted separately rather than excluded silently, a growing connector count is
# worth noticing, just not a failure.
.hf_poi_outlet_has_divide <- function(pois, lines, divides, dec, line_id_col,
                                     tol_m = 30) {
  if (is.null(dec) || is.null(lines) || is.null(divides)) return(NULL)
  if (!"flowpath_id" %in% names(divides)) return(NULL)
  at <- which(!is.na(dec$d_out) & dec$d_out <= tol_m)
  if (!length(at)) return(NULL)

  own_i <- dec$idx[at]
  owner <- as.character(lines[[line_id_col]])[own_i]
  dv    <- unique(as.character(divides$flowpath_id))
  no_div <- is.na(owner) | !nzchar(owner) | !(owner %in% dv)

  # has_divide = 0 marks a reach that carries no catchment by construction.
  hd <- if ("has_divide" %in% names(lines))
    suppressWarnings(as.integer(lines$has_divide))[own_i] else
    rep(NA_integer_, length(owner))
  connector <- !is.na(hd) & hd == 0L

  unexplained <- no_div & !connector
  .hf_ok(sum(unexplained) == 0L, sprintf(paste0(
    "%d of %d POI(s) at an outlet on a catchment-bearing reach have no divide; ",
    "%d more sit on a connector (no divide by construction, accepted)"),
    sum(unexplained), length(at), sum(no_div & connector)))
}

# A minted POI id is a commitment: it must be present, unique, and permanent.
# These two enter as GATES rather than measurements because both are already at
# zero: gating them now stops a regression rather than recording one, and the
# aggregation work that follows is exactly what would introduce a NA'd or
# re-derived id.
.hf_poi_ids_minted <- function(pois, hydrolocations = NULL) {
  bad_na <- 0L; bad_dup <- 0L; n <- 0L; seen <- FALSE
  for (d in list(pois, hydrolocations)) {
    if (is.null(d) || !"poi_id" %in% names(d)) next
    seen <- TRUE
    x <- as.character(d$poi_id)
    n       <- n + length(x)
    bad_na  <- bad_na  + sum(is.na(x) | !nzchar(x))
    keep    <- x[!is.na(x) & nzchar(x)]
    bad_dup <- bad_dup + sum(duplicated(keep))
  }
  if (!seen) return(NULL)
  .hf_ok(bad_na == 0L && bad_dup == 0L, sprintf(
    "%d of %d delivered POI row(s) carry no minted poi_id; %d duplicate id(s)",
    bad_na, n, bad_dup))
}

# An id in the product that was never minted upstream means something downstream
# invented an identity, which breaks the lineage the id exists to carry.
.hf_poi_ids_not_reminted <- function(pois, pois_input) {
  if (is.null(pois) || is.null(pois_input)) return(NULL)
  if (!"poi_id" %in% names(pois) || !"poi_id" %in% names(pois_input)) return(NULL)
  mint <- as.character(pois_input$poi_id)
  mint <- unique(mint[!is.na(mint) & nzchar(mint)])
  if (!length(mint)) return(NULL)
  out <- as.character(pois$poi_id)
  out <- unique(out[!is.na(out) & nzchar(out)])
  new <- setdiff(out, mint)
  .hf_ok(length(new) == 0L, sprintf(
    "%d of %d delivered poi_id(s) were never minted upstream%s",
    length(new), length(out),
    if (length(new)) sprintf(" (e.g. %s)", paste(utils::head(new, 3), collapse = ", ")) else ""))
}

# A POI id is a single integer. Co-located hydrolocations merge into ONE poi_id
# and the lineage keeps every hl_id, so a comma-joined or non-integer poi_id is a
# defect, not a representation of several POIs. Two aggregation sites build
# `poi_id` with paste(collapse = ","), and one coerces the result with
# as.numeric(), which yields NA for a multi-id string. This is what would break
# first if a merge ever put two distinct POIs on one aggregate.
.hf_poi_id_is_integer <- function(pois, hydrolocations = NULL, network = NULL) {
  bad_comma <- 0L; bad_int <- 0L; n <- 0L; seen <- FALSE; ex <- character(0)
  for (d in list(pois, hydrolocations, network)) {
    if (is.null(d) || !"poi_id" %in% names(d)) next
    seen <- TRUE
    v <- as.character(d$poi_id)
    v <- v[!is.na(v) & nzchar(v)]
    if (!length(v)) next
    n <- n + length(v)
    cj <- grepl(",", v, fixed = TRUE)
    bad_comma <- bad_comma + sum(cj)
    bad_int <- bad_int + sum(is.na(suppressWarnings(as.integer(v))))
    if (any(cj) && !length(ex)) ex <- utils::head(v[cj], 2L)
  }
  if (!seen) return(NULL)
  .hf_ok(bad_comma == 0L && bad_int == 0L, sprintf(
    "%d of %d poi_id value(s) are comma-joined and %d are non-integer%s",
    bad_comma, n, bad_int,
    if (length(ex)) sprintf(" (e.g. %s)", paste(ex, collapse = ", ")) else ""))
}

# A POI must coincide with the outlet of a post-split derivative (contract step 4).
#
# Measured against the flowlines layer independently of the aggregate join, because
# a POI marks a flowline boundary and the flowline on its upstream side can belong
# to a different aggregate. Point-to-point, for the same reason as above.
.hf_poi_at_flowline_outlet <- function(pois, flowlines, tol_m = 30) {
  if (is.null(pois) || is.null(flowlines)) return(NULL)
  if (!inherits(pois, "sf") || !inherits(flowlines, "sf")) return(NULL)
  if (!nrow(pois) || !nrow(flowlines)) return(NULL)
  tryCatch({
    p <- if (sf::st_crs(pois) != sf::st_crs(flowlines))
      sf::st_transform(pois, sf::st_crs(flowlines)) else pois
    outs <- sf::st_sfc(lapply(sf::st_geometry(flowlines), function(g) {
      m <- sf::st_coordinates(g); sf::st_point(unname(m[nrow(m), 1:2]))
    }), crs = sf::st_crs(flowlines))
    nf <- sf::st_nearest_feature(sf::st_geometry(p), outs)
    d  <- as.numeric(sf::st_distance(sf::st_geometry(p), outs[nf],
      by_element = TRUE))
    off <- d > tol_m
    types <- if (any(off)) {
      tc <- if ("hl_types" %in% names(p)) p$hl_types else
            if ("hl_type" %in% names(p)) p$hl_type else NULL
      if (is.null(tc)) "" else {
        tt <- sort(table(as.character(tc)[off]), decreasing = TRUE)
        sprintf(" [%s]", paste(sprintf("%s=%d", names(tt), as.integer(tt)),
          collapse = ", "))
      }
    } else ""
    .hf_info(sprintf(paste0(
      "%d of %d POI(s) coincide with a reconciled flowline outlet (%.1f%% exact, ",
      "%.1f%% within %.0fm); %d beyond%s (target 100%%)"),
      sum(d <= 1), length(d), 100 * mean(d <= 1), 100 * mean(!off), tol_m,
      sum(off), types))
  }, error = function(e) NULL)
}

# Contract step 4: a delivered POI must still represent the location stage 1
# established. Every other property can look right while the point itself has
# moved, which relocates what the fabric is reported at.
#
# Two measurement rules make the number mean anything.
#
# A poi_id can cover several stage-1 events, bundled by collocation and kilometres
# apart along their shared reach, so distance is taken to the NEAREST event
# carrying that poi_id. Otherwise this measures the spread between hydrolocations
# rather than movement of the POI.
#
# Movement is not instability. A POI placed at the outlet of the reach it sits on
# has been placed, not drifted, so those are separated out. The residual: moved
# and landed nowhere in particular: is the number worth reading.
.hf_poi_location_stable <- function(pois, pois_input, flowlines = NULL,
                                    tol_m = 30) {
  if (is.null(pois) || is.null(pois_input)) return(NULL)
  if (!inherits(pois, "sf") || !inherits(pois_input, "sf")) return(NULL)
  if (!"poi_id" %in% names(pois) || !"poi_id" %in% names(pois_input)) return(NULL)
  tryCatch({
    a <- pois_input[!is.na(pois_input$poi_id) &
      nzchar(as.character(pois_input$poi_id)), ]
    if (!nrow(a)) return(NULL)
    b <- if (sf::st_crs(pois) != sf::st_crs(a))
      sf::st_transform(pois, sf::st_crs(a)) else pois
    aid  <- as.character(a$poi_id)
    bid  <- as.character(b$poi_id)
    keep <- which(bid %in% aid)
    if (!length(keep)) return(NULL)

    # Nearest stage-1 event sharing the poi_id, per delivered POI.
    grp <- split(seq_len(nrow(a)), aid)
    d <- vapply(keep, function(i) {
      g <- grp[[bid[i]]]
      min(as.numeric(sf::st_distance(sf::st_geometry(b)[i],
        sf::st_geometry(a)[g])), na.rm = TRUE)
    }, numeric(1))
    moved <- d > tol_m

    # Placement accounts for a move when the POI landed on a flowline outlet.
    at_outlet <- rep(NA, length(keep))
    if (!is.null(flowlines) && inherits(flowlines, "sf") && nrow(flowlines)) {
      outs <- sf::st_sfc(lapply(sf::st_geometry(flowlines), function(g) {
        m <- sf::st_coordinates(g); sf::st_point(unname(m[nrow(m), 1:2]))
      }), crs = sf::st_crs(flowlines))
      if (sf::st_crs(outs) != sf::st_crs(b)) outs <- sf::st_transform(outs, sf::st_crs(b))
      nf <- sf::st_nearest_feature(sf::st_geometry(b)[keep], outs)
      at_outlet <- as.numeric(sf::st_distance(sf::st_geometry(b)[keep],
        outs[nf], by_element = TRUE)) <= 1
    }
    placed      <- sum(moved & .hf_is_true(at_outlet), na.rm = TRUE)
    unexplained <- sum(moved & !.hf_is_true(at_outlet), na.rm = TRUE)

    .hf_info(sprintf(paste0(
      "%d of %d POI(s) moved more than %.0fm from the location stage 1 ",
      "established (max %.0fm; %d moved at all); %d of those sit on a flowline ",
      "outlet (placed by design), leaving %d unexplained (target 0)"),
      sum(moved, na.rm = TRUE), length(d), tol_m,
      if (any(is.finite(d))) max(d, na.rm = TRUE) else 0,
      sum(d > 1, na.rm = TRUE), placed, unexplained))
  }, error = function(e) NULL)
}

# TRUE only where x is TRUE; NA (outlet position unknown) counts as not-placed so
# the unexplained residual is never understated.
.hf_is_true <- function(x) !is.na(x) & x

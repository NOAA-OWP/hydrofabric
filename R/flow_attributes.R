# =============================================================================
# Channel / routing flow attributes
# -----------------------------------------------------------------------------
# `.calculate_flow_atts` and `.apply_geomorphic_proxy` implement the RiverML
# flow-attribute model with a geomorphic-proxy fallback. Self-contained (no
# runtime dependency): `where()` is qualified as `dplyr::where()`. Used by
# Stage 3 Step 6g (routing attributes).
# =============================================================================

# Dingman's generalised cross-section, as an alternative hydraulics source.
#
# The OWP bankfull tables describe a channel by top width, depth and a shape
# exponent r rather than by area: r = 1 is a triangle, large r approaches a
# rectangle. Bankfull area follows,
#
#     A = W * D * r / (r + 1)
#
# which is WD/2 at r = 1 and tends to WD as r grows. That is exactly the bf_area
# the trapezoid inversion below needs, so a table in this form is renamed to the
# canonical columns and given its area; the rest of the derivation is unchanged.
# Two sources, one set of formulas.
# Column names for a Dingman-form channel table. The geometry columns appear
# with and without an `owp_` prefix, so each is resolved from its candidates.
.HF_DINGMAN_COLS <- list(
  tw    = c("owp_tw_bf_m", "tw_bf_m"),
  y     = c("owp_y_bf_m",  "y_bf_m"),
  r     = c("owp_r_bf",    "r_bf"),
  # `n` takes the in-channel roughness and `nCC` the overbank. A table with a
  # single roughness column supplies both.
  n     = c("n_in_channel", "n_single"),
  n_cc  = c("n_out_channel"),
  slope = c("slope_m_m"))

.hf_pick_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit)) hit[[1]] else NA_character_
}

.is_dingman_hydraulics <- function(x) {
  nms <- names(x)
  !anyNA(vapply(.HF_DINGMAN_COLS[c("tw", "y", "r")],
                function(cands) .hf_pick_col(nms, cands), character(1)))
}

.as_canonical_hydraulics <- function(x, lengths = NULL) {
  if (!.is_dingman_hydraulics(x)) return(x)
  d   <- if (inherits(x, "ArrowObject")) dplyr::collect(x) else x
  key <- if ("feature_id" %in% names(d)) "feature_id" else "flowpath_id"
  col <- lapply(.HF_DINGMAN_COLS, function(cands) .hf_pick_col(names(d), cands))
  tw  <- d[[col$tw]]
  y   <- d[[col$y]]
  r   <- d[[col$r]]
  out <- dplyr::tibble(
    feature_id          = hfutils::hf_fmt_id(d[[key]]),
    owp_tw_bf           = tw,
    owp_y_bf            = y,
    bf_area             = tw * y * r / (r + 1),
    owp_roughness_bathy = if (!is.na(col$n)) d[[col$n]] else NA_real_,
    slope               = if (!is.na(col$slope)) d[[col$slope]] else NA_real_)
  if (!is.na(col$n_cc)) out$owp_roughness_overbank <- d[[col$n_cc]]
  if (!is.null(lengths)) {
    lengths <- dplyr::rename_with(lengths, ~ "feature_id",
                                  .cols = dplyr::any_of(c("reference_id", "flowpath_id")))
    lengths$feature_id <- hfutils::hf_fmt_id(lengths$feature_id)
    out <- dplyr::left_join(out, lengths[, c("feature_id", "lengthkm")],
                            by = "feature_id")
  }
  # Weights the RiverML table carries inline and this one does not. Equal
  # weighting is stated rather than assumed silently.
  if (!"lengthkm" %in% names(out)) {
    message("calculate_flow_atts: no reach lengths supplied; ",
            "weighting the reaches behind each feature equally")
    out$lengthkm <- 1
  }
  out
}

#' Compute channel routing attributes (RiverML + fallback chain)
#'
#' For each hydrofabric feature, derives trapezoidal-channel routing parameters
#' (BtmWdth, TopWdth, ChSlp, n, nCC, So, TopWdthCC, Muskingum MusX/MusK) from a
#' RiverML hydraulics table, falling back to a routelink table and finally a
#' geomorphic proxy (`.apply_geomorphic_proxy`) for features with no estimate.
#'
#' @param network Network table/sf with `hf_id`, length, and (for the proxy)
#'   `total_dasqkm` / `slope` columns.
#' @param hydraulics RiverML hydraulics table (data.frame or Arrow object),
#'   keyed by `feature_id`.
#' @param hf_id Join id column: `"flowline_id"` or `"flowpath_id"`.
#' @param attribute_cols RiverML columns aggregated per `hf_id`.
#' @param weight_col Length-weight column for the weighted means.
#' @param flowpaths Optional flowpaths sf supplying `total_dasqkm`/`slope` to the
#'   proxy fallback when the network lacks them.
#' @param routelink Optional routelink table (tier-2 fallback).
#' @param proxy_params Optional named list overriding `.apply_geomorphic_proxy`
#'   coefficients.
#' @param reach_lengths Optional data frame of `reference_id`/`flowpath_id` and
#'   `lengthkm`, weighting the reaches behind each feature when the hydraulics
#'   table carries no length of its own. Only consulted for a Dingman-form
#'   table; the RiverML table supplies its own weights.
#' @return A data frame of routing attributes joined to `network` length, one
#'   row per `hf_id`.
#' @keywords internal
#' @noRd
.calculate_flow_atts <- function(network, hydraulics, hf_id = "flowline_id", attribute_cols = c("owp_y_bf",
      "owp_tw_bf", "bf_area", "owp_roughness_bathy", "slope"),
    weight_col = "lengthkm", flowpaths = NULL, routelink = NULL,
    proxy_params = list(), reach_lengths = NULL) {
  # A Dingman-form table is renamed and given its bf_area here, so everything
  # below sees one schema regardless of which source supplied it.
  hydraulics <- .as_canonical_hydraulics(hydraulics, reach_lengths)
  net_ref_col <- if ("reference_id" %in% names(network))
    "reference_id"
  else "hf_id"
  # Reference ids are matched as formatted strings. NHD-HR ids exceed 2^31,
  # which as.integer() returns as NA rather than erroring on.
  net_ids <- hfutils::hf_fmt_id(unique(network[[net_ref_col]]))
  # `attribute_cols` names the RiverML schema, so an overbank column is
  # selected separately when the source carries one.
  has_overbank <- "owp_roughness_overbank" %in% names(hydraulics)
  sel_cols <- c(attribute_cols, if (has_overbank) "owp_roughness_overbank")
  map <- dplyr::filter(dplyr::mutate(dplyr::collect(dplyr::select(hydraulics,
    reference_id = feature_id, dplyr::all_of(sel_cols),
    small_w = dplyr::all_of(weight_col))),
  reference_id = hfutils::hf_fmt_id(reference_id)), reference_id %in%
    net_ids)
  net_for_idx <- network
  if (!"hl_class" %in% names(net_for_idx)) {
    net_for_idx$hl_class <- dplyr::if_else(!is.na(net_for_idx$hl_reference) &
      nzchar(net_for_idx$hl_reference), "gage", "other")
  }
  net_for_idx$reference_id <- hfutils::hf_fmt_id(net_for_idx[[net_ref_col]])
  agg_cols <- sel_cols

  idx <- net_for_idx |>
    dplyr::select(reference_id, hl_class, hl_reference, dplyr::all_of(hf_id)) |>
    dplyr::mutate(hl_reference = ifelse(hl_class == "gage", hl_reference, NA)) |>
    dplyr::select(-hl_class) |>
    dplyr::distinct() |>
    dplyr::left_join(map, by = "reference_id") |>
    dplyr::group_by(!!dplyr::sym(hf_id)) |>
    dplyr::mutate(
      hl_reference = paste(stats::na.omit(hl_reference), collapse = ",")) |>
    dplyr::summarise(
      dplyr::across(.cols = dplyr::all_of(agg_cols),
        .fns = ~ weighted.mean(.x, w = small_w, na.rm = TRUE),
        .names = "{.col}"),
      hl_reference = hl_reference[1], .groups = "drop") |>
    dplyr::select(dplyr::all_of(hf_id), n = owp_roughness_bathy,
      Tw = owp_tw_bf, bf_area, Y = owp_y_bf, So = slope, hl_reference,
      dplyr::any_of("owp_roughness_overbank")) |>
    dplyr::mutate(
      dplyr::across(dplyr::where(is.logical), as.double),
      Bw        = pmax((2 * bf_area / Y) - Tw, 0.1),
      BtmWdth   = pmin(Bw, Tw),
      TopWdth   = pmax(Bw, Tw),
      nCC       = if (has_overbank)
                    dplyr::coalesce(owp_roughness_overbank, 2 * n) else 2 * n,
      TopWdthCC = 3 * TopWdth,
      ChSlp     = pmax(((TopWdth - BtmWdth) / 2) / Y, 1e-05),
      BtmWdth   = dplyr::if_else(ChSlp == 1e-05,
                    TopWdth - (2 * ChSlp * Y), BtmWdth),
      MusX      = 0.2,
      MusK      = 3600,
      attr_source = dplyr::if_else(!is.na(BtmWdth), "riverml", NA_character_)) |>
    dplyr::select(dplyr::all_of(hf_id), BtmWdth, TopWdth, ChSlp, n, nCC, So,
      TopWdthCC, hl_reference, MusX, MusK, attr_source)
  missing_mask <- is.na(idx$BtmWdth)
  n_missing <- sum(missing_mask)
  if (n_missing > 0) {
    message(sprintf("calculate_flow_atts: %d %s(s) have no RiverML estimate; applying fallback chain",
      n_missing, hf_id))
    missing_ids <- idx[[hf_id]][missing_mask]
    tier2 <- if (!is.null(routelink)) {
      rl <- if (inherits(routelink, "ArrowObject")) {
        dplyr::collect(routelink)
      } else {
        routelink
      }
      dplyr::mutate(dplyr::filter(rl, .data[[hf_id]] %in%
        missing_ids), attr_source = "routelink")
    } else {
      dplyr::tibble()
    }
    still_missing <- setdiff(missing_ids, if (nrow(tier2) >
      0)
      tier2[[hf_id]]
    else character(0))
    tier3 <- if (length(still_missing) > 0) {
      proxy_net <- dplyr::ungroup(dplyr::slice_max(dplyr::group_by(dplyr::select(dplyr::filter(sf::st_drop_geometry(network),
        .data[[hf_id]] %in% still_missing), dplyr::all_of(hf_id),
      dplyr::any_of(c("flowpath_id", "total_dasqkm",
        "slope"))), .data[[hf_id]]), dplyr::pick(dplyr::any_of("total_dasqkm")),
      n = 1, with_ties = FALSE))
      needed_from_fp <- setdiff(c("total_dasqkm", "slope"),
        names(proxy_net))
      if (!is.null(flowpaths) && length(needed_from_fp) >
        0 && "flowpath_id" %in% names(proxy_net)) {
        fp_cols <- dplyr::distinct(dplyr::select(sf::st_drop_geometry(flowpaths),
          flowpath_id, dplyr::any_of(needed_from_fp)))
        proxy_net <- dplyr::left_join(proxy_net, fp_cols,
          by = "flowpath_id")
      }
      if ("slope" %in% names(proxy_net)) {
        proxy_net <- dplyr::rename(proxy_net, So = slope)
      } else if (!"So" %in% names(proxy_net)) {
        proxy_net$So <- NA_real_
      }
      proxy_args <- c(list(df = proxy_net), proxy_params)
      result3 <- do.call(.apply_geomorphic_proxy, proxy_args)
      if (!hf_id %in% names(result3)) {
        result3[[hf_id]] <- proxy_net[[hf_id]]
      }
      dplyr::select(result3, dplyr::all_of(hf_id), BtmWdth,
        TopWdth, ChSlp, n, nCC, So, TopWdthCC, hl_reference,
        MusX, MusK, attr_source)
    } else {
      dplyr::tibble()
    }
    patch <- dplyr::bind_rows(tier2, tier3)
    if (nrow(patch) > 0) {
      idx <- dplyr::rows_patch(idx, patch, by = hf_id,
        unmatched = "ignore")
    }
  }
  # Floor channel slope (So) to a positive default for ALL tiers. Source slope can
  # be zero or negative on reversed and degenerate-length segments, and a
  # non-positive So is invalid for routing. One floor here covers every tier.
  idx$So <- dplyr::if_else(!is.na(idx$So) & idx$So > 0, idx$So, 1e-04)
  length_column_pref <- gsub("id", "lengthkm", hf_id)
  length_column <- if (length_column_pref %in% names(network)) {
    length_column_pref
  } else if ("lengthkm" %in% names(network)) {
    "lengthkm"
  } else {
    stop("cannot find a length column for hf_id = '", hf_id,
      "'; tried '", length_column_pref, "' and 'lengthkm'")
  }
  dplyr::right_join(dplyr::mutate(dplyr::distinct(dplyr::select(network,
    dplyr::all_of(hf_id), lengthkm = dplyr::all_of(length_column))),
  lengthm = lengthkm * 1000, lengthkm = NULL), idx, by = hf_id)
}

#' Geomorphic-proxy channel attributes (tier-3 fallback)
#'
#' Hydraulic-geometry estimate of trapezoidal-channel parameters from drainage
#' area when no RiverML/routelink attributes exist.
#'
#' @param df Data frame with `total_dasqkm` (and optionally `So`).
#' @param a_w,b_w,a_y,b_y Hydraulic-geometry coefficients (width/depth vs area).
#' @param n_default,slope_default,da_min Defaults for roughness, slope, min area.
#' @return `df` with proxy routing attributes and `attr_source = "geomorphic_proxy"`.
#' @keywords internal
#' @noRd
.apply_geomorphic_proxy <- function(df, a_w = 2.56, b_w = 0.45, a_y = 0.27, b_y = 0.294,
    n_default = 0.035, slope_default = 1e-04, da_min = 0.01) {
  dplyr::select(dplyr::mutate(df, total_dasqkm = dplyr::if_else(is.na(total_dasqkm) |
    total_dasqkm <= 0, da_min, total_dasqkm), Tw = a_w *
    total_dasqkm^b_w, Y = a_y * total_dasqkm^b_y, Bw = pmax(0.1,
    0.5 * Tw), So = dplyr::if_else(!is.na(So) & So > 0, So,
    slope_default), n = n_default, BtmWdth = pmin(Bw, Tw),
  TopWdth = pmax(Bw, Tw), nCC = 2 * n, TopWdthCC = 3 *
    TopWdth, ChSlp = pmax(((TopWdth - BtmWdth) / 2) / Y,
    1e-05), BtmWdth = dplyr::if_else(ChSlp == 1e-05,
    TopWdth - (2 * ChSlp * Y), BtmWdth), MusX = 0.2,
  MusK = 3600, hl_reference = NA_character_, attr_source = "geomorphic_proxy"),
  dplyr::any_of(c("flowline_id", "flowpath_id")), BtmWdth,
  TopWdth, ChSlp, n, nCC, So, TopWdthCC, hl_reference,
  MusX, MusK, attr_source)
}

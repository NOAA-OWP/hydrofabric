# =============================================================================
# Catchment-size grouping (vendored from an internal aggregation utility)
# -----------------------------------------------------------------------------
# .cs_group and its helpers (.splitAt, .assign_id, .pinch_sides,
# .middle_massage) were vendored verbatim to drop the external aggregation
# runtime dependency. Logic is unchanged; only the names are made
# package-private and internal calls rewired to the dotted names. .cs_group
# assigns contiguous flowlines along a levelpath to size-balanced aggregate
# groups (target ~ ideal_size_sqkm), honouring exclusion masks and a maximum
# junction count. Used by aggregate_network() (Step 4).
# =============================================================================

#' Split a vector at the given positions (after each pos index)
#' @keywords internal
#' @noRd
.splitAt <- function(x, pos) {
  unname(split(x, cumsum(seq_along(x) %in% (pos + 1))))
}

#' Greedy cumulative-area group assignment (resets when athres exceeded)
#' @keywords internal
#' @noRd
.assign_id <- function(x, athres) {
  cumsum <- 0
  group <- 1
  result <- numeric()
  for (i in seq_along(x)) {
    cumsum <- cumsum + x[i]
    if (cumsum > athres) {
      group <- group + 1
      cumsum <- x[i]
    }
    result[i] <- group
  }
  return(result)
}

#' Merge below-threshold end groups into their inward neighbour
#' @keywords internal
#' @noRd
.pinch_sides <- function(x, ind, thres) {
  group_sums <- tapply(x, ind, sum)
  group_ids <- as.numeric(names(group_sums))
  if (length(group_sums) == 1) {
    return(ind)
  }
  if (group_sums[1] < thres)
    names(group_sums)[1] <- names(group_sums)[2]
  if (group_sums[length(group_sums)] < thres)
    names(group_sums)[length(group_sums)] <- names(group_sums)[length(group_sums) -
      1]
  as.numeric(names(group_sums))[match(ind, group_ids)]
}

#' Merge below-threshold interior groups into their smaller neighbour
#' @keywords internal
#' @noRd
.middle_massage <- function(x, index_values, threshold) {
  group_sums <- tapply(x, index_values, sum)
  group_ids <- as.numeric(names(group_sums))
  if (length(group_sums) == 1) {
    return(index_values)
  }
  small_groups <- which(group_sums < threshold)
  for (grp in small_groups) {
    if (grp == 1 || grp == length(group_sums))
      next
    neighbors <- c(grp - 1, grp + 1)
    target <- neighbors[which.min(group_sums[neighbors])]
    names(group_sums)[grp] <- names(group_sums)[target]
  }
  as.numeric(names(group_sums))[match(index_values, group_ids)]
}

#' Assign contiguous flowlines to size-balanced aggregate groups
#'
#' @param areas,lengths Per-flowline incremental area (sqkm) and length (km),
#'   in downstream order along a levelpath.
#' @param max_junction Max features per junction group (NULL/NA/Inf = unbounded).
#' @param exclude_dn,exclude_un Downstream/upstream exclusion masks (non-NA marks
#'   a forced break).
#' @param ideal_size_sqkm Target aggregate catchment area.
#' @param amin,lmin Minimum area / length below which a group is merged.
#' @return Integer group id per input flowline, offset by 1e9 * chunk index.
#' @keywords internal
#' @noRd
.cs_group <- function(areas, lengths, max_junction, exclude_dn, exclude_un,
                      ideal_size_sqkm, amin, lmin, has_divide = NULL) {
  areas[is.na(areas)] <- 0
  lengths[is.na(lengths)] <- 0
  if (length(areas) == 1)
    return(1)
  breaks <- sort(c(which(!is.na(exclude_dn)), which(!is.na(exclude_un)) -
    1))
  # Break the levelpath at every has_divide TRANSITION so catchment-bearing
  # reaches are never grouped with divide-less connector reaches (reference
  # areasqkm=0 waterbody / artificial paths). Consecutive connectors still group
  # together; consecutive catchment reaches group normally.
  if (!is.null(has_divide)) {
    hd <- as.integer(has_divide)
    hd[is.na(hd)] <- 1L
    breaks <- sort(unique(c(breaks, which(diff(hd) != 0L))))
  }
  area_chunks <- .splitAt(areas, breaks)
  length_chunks <- .splitAt(lengths, breaks)
  stopifnot(length(area_chunks) == length(length_chunks))
  initial <- lapply(area_chunks, .assign_id, athres = ideal_size_sqkm)
  pinched_area <- mapply(.pinch_sides, area_chunks, initial,
    MoreArgs = list(thres = amin), SIMPLIFY = FALSE)
  pinched_len <- mapply(.pinch_sides, length_chunks, pinched_area,
    MoreArgs = list(thres = lmin), SIMPLIFY = FALSE)
  massaged_area <- mapply(.middle_massage, area_chunks, pinched_len,
    MoreArgs = list(threshold = amin), SIMPLIFY = FALSE)
  massaged_len <- mapply(.middle_massage, length_chunks, massaged_area,
    MoreArgs = list(threshold = lmin), SIMPLIFY = FALSE)
  enforce_max_junction <- function(ids, max_junction) {
    if (is.null(max_junction) || is.na(max_junction) || is.infinite(max_junction)) {
      return(ids)
    }
    max_features <- max(1, max_junction + 1)
    r <- rle(ids)
    new_ids <- integer(length(ids))
    idx <- 1L
    local_group_counter <- 0
    for (k in seq_along(r$lengths)) {
      L <- r$lengths[k]
      nseg <- ceiling(L / max_features)
      for (s in seq_len(nseg)) {
        len_s <- if (s < nseg)
          max_features
        else L - max_features * (nseg - 1)
        local_group_counter <- local_group_counter +
          1
        new_ids[idx:(idx + len_s - 1)] <- local_group_counter
        idx <- idx + len_s
      }
    }
    new_ids
  }
  massaged_len <- lapply(massaged_len, enforce_max_junction,
    max_junction = max_junction)
  for (i in seq_along(massaged_len)) {
    massaged_len[[i]] <- massaged_len[[i]] + 1e+09 * i
  }
  unlist(massaged_len)
}

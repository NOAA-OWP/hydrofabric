ls_env <- function(env) { ls(pos = env) }

#' Conflicts between the hydrofabric and other packages
#'
#' This function lists all the conflicts between packages in the hydrofabric
#' and other packages that you have loaded.
#'
#' There are four conflicts that are deliberately ignored: \code{intersect},
#' \code{union}, \code{setequal}, and \code{setdiff} from dplyr. These functions
#' make the base equivalents generic, so shouldn't negatively affect any
#' existing code.
#'
#' @export
#' @examples
#' hydrofabric_conflicts()

hydrofabric_conflicts <- function() {
  
  envs <- grep("^package:", base::search(), value = TRUE)
  envs <- purrr::set_names(envs)
  objs <- invert(lapply(envs, ls_env))
  
  conflicts <- purrr::keep(objs, ~ length(.x) > 1)
  
  tidy_names <- paste0("package:", hydrofabric_packages())
  conflicts <- purrr::keep(conflicts, ~ any(.x %in% tidy_names))
  
  conflict_funs <- purrr::imap(conflicts, confirm_conflict)
  conflict_funs <- purrr::compact(conflict_funs)
  
  structure(conflict_funs, class = "hydrofabric_conflicts")
}

hydrofabric_conflict_message <- function(x) {
  if (length(x) == 0) return("")
  
  header <- cli::rule(
    left = crayon::bold("Conflicts"),
    right = "hydrofabric_conflicts()"
  )
  
  pkgs <- x %>% purrr::map(~ gsub("^package:", "", .))
  others <- pkgs %>% purrr::map(`[`, -1)
  other_calls <- purrr::map2_chr(
    others, names(others),
    ~ paste0(crayon::blue(.x), "::", .y, "()", collapse = ", ")
  )
  
  winner <- pkgs %>% purrr::map_chr(1)
  funs <- format(paste0(crayon::blue(winner), "::", crayon::green(paste0(names(x), "()"))))
  bullets <- paste0(
    crayon::red(cli::symbol$cross), " ", funs,
    " masks ", other_calls,
    collapse = "\n"
  )
  
  paste0(header, "\n", bullets)
}

#' @export
print.hydrofabric_conflicts <- function(x, ..., startup = FALSE) {
  cli::cat_line(hydrofabric_conflict_message(x))
}

#' @importFrom magrittr %>%
#' @importFrom purrr map keep

confirm_conflict <- function(packages, name) {
  # Only look at functions
  objs <- packages %>%
    purrr::map(~ get(name, pos = .)) %>%
    purrr::keep(is.function)
  
  if (length(objs) <= 1)
    return()
  
  # Remove identical functions
  objs <- objs[!duplicated(objs)]
  packages <- packages[!duplicated(packages)]
  if (length(objs) == 1)
    return()
  
  packages
}


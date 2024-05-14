ls_env <- function(env) {
  x <- ls(pos = env)
  
  # intersect, setdiff, setequal, union come from generics
  if (env %in% c("package:dplyr", "package:lubridate")) {
    x <- setdiff(x, c("intersect", "setdiff", "setequal", "union"))
  }
  
  if (env == "package:lubridate") {
    x <- setdiff(x, c(
      "as.difftime", # lubridate makes into an S4 generic
      "date"         # matches base behaviour
    ))
  }
  
  x
}

#' Conflicts between the hydrofabric and other packages
#'
#' This function lists all the conflicts between packages in the hydrofabric
#' and other packages that you have loaded.
#'
#' @export
#' @examples
#' hydrofabric_conflicts()

hydrofabric_conflicts <- function(only = NULL) {
  
  envs <- grep("^package:", base::search(), value = TRUE)
  envs <- purrr::set_names(envs)
  
  if (!is.null(only)) {
    only <- union(only, core)
    envs <- envs[names(envs) %in% paste0("package:", only)]
  }
  
  objs <- invert(lapply(envs, ls_env))
  
  conflicts <- purrr::keep(objs, ~ length(.x) > 1)
  
  tidy_names <- paste0("package:", hydrofabric_packages(include_self = FALSE))
  #tidy_names = tidy_names[tidy_names != "package:hydrofabric"]
  conflicts <- purrr::keep(conflicts, ~ any(.x %in% tidy_names))

  
  conflict_funs <- purrr::imap(conflicts, confirm_conflict)
  conflict_funs <- purrr::compact(conflict_funs)
  
  c = lapply(1:length(conflict_funs), function(x){
    conflict_funs[[x]] = conflict_funs[[x]][!conflict_funs[[x]] %in% 
                                              c("package:hydrofabric", 
                                                "package:base",
                                                'package:stats',
                                                'package:graphics',
                                                'package:utils',
                                                'package:grDevices',
                                                'package:testthat')]
  })
  
  names(c) = names(conflict_funs)
  c = c[lengths(c) > 1]
  
  structure(c, class = "hydrofabric_conflicts")
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
  invisible(x)
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


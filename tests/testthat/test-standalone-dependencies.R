## Standalone dependencies
##
## `hydrofabric` must build and run from publicly available packages alone. It is
## the core library; orchestration, publishing and subsetting live in separate
## packages that depend on it, and the dependency must not run the other way,
## a core that reaches back into a downstream layer cannot be installed on its
## own.
##
## This is stated as an allowlist rather than a list of things to avoid, so a
## dependency added later has to be added here deliberately. Undeclared `pkg::`
## calls are already caught by R CMD check; these tests add the
## declared-dependency and source-usage gates on top.

PUBLIC <- c(
  # base / recommended
  "stats", "utils", "parallel", "methods", "tools", "grDevices", "graphics",
  "grid", "splines", "compiler", "datasets",
  # CRAN
  "cli", "data.table", "DBI", "dplyr", "furrr", "future", "geos", "igraph",
  "jsonlite", "logger", "lwgeom", "RSQLite", "sf", "terra",
  "arrow", "knitr", "rmarkdown", "testthat", "withr",
  # Apache-2.0, public: github.com/lynker-spatial/hfutils
  "hfutils"
)

.hf_pkg_root <- function() {
  root <- tryCatch(find.package("hydrofabric"), error = function(e) NA_character_)
  if (is.na(root) || !file.exists(file.path(root, "DESCRIPTION")))
    root <- testthat::test_path("..", "..")
  root
}

.declared <- function() {
  d <- read.dcf(file.path(.hf_pkg_root(), "DESCRIPTION"))
  fields <- intersect(colnames(d),
                      c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances"))
  if (!length(fields)) return(character(0))
  raw <- unlist(strsplit(paste(d[, fields], collapse = ","), ","))
  raw <- sub("\\(.*", "", raw)                       # drop version constraints
  raw <- gsub("[[:space:]]", "", raw)
  setdiff(unique(raw[nzchar(raw)]), "R")             # R itself is not a package
}

test_that("every declared dependency is a public package", {
  extra <- setdiff(.declared(), PUBLIC)
  expect_identical(extra, character(0),
    label = paste0("dependencies not on the public allowlist: ",
                   paste(extra, collapse = ", ")))
})

test_that("NAMESPACE imports only public packages", {
  ns <- file.path(.hf_pkg_root(), "NAMESPACE")
  skip_if_not(file.exists(ns))
  txt <- readLines(ns, warn = FALSE)
  imp <- regmatches(txt, regexpr("^import(From)?\\(([^,)]+)", txt))
  imp <- sub("^import(From)?\\(", "", imp)
  imp <- gsub("[\"']", "", imp)
  extra <- setdiff(unique(imp[nzchar(imp)]), PUBLIC)
  expect_identical(extra, character(0),
    label = paste0("NAMESPACE imports outside the allowlist: ",
                   paste(extra, collapse = ", ")))
})

test_that("no R/ source reaches for a package outside the allowlist", {
  r_dir <- file.path(.hf_pkg_root(), "R")
  skip_if_not(dir.exists(r_dir), "R/ source not present (installed-package check)")
  used <- character(0)
  for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
    lines <- readLines(f, warn = FALSE)
    code  <- lines[!grepl("^[[:space:]]*#", lines)]   # drop comments, incl. roxygen
    # pkg::fn and pkg:::fn
    m <- regmatches(code, gregexpr("[a-zA-Z][a-zA-Z0-9._]*(?=:::?)", code, perl = TRUE))
    used <- c(used, unlist(m))
    # library(pkg) / requireNamespace("pkg") and friends
    m2 <- regmatches(code, gregexpr(
      "(?<=library\\(|require\\(|requireNamespace\\(|loadNamespace\\(|attachNamespace\\()['\"]?[a-zA-Z][a-zA-Z0-9._]*",
      code, perl = TRUE))
    used <- c(used, gsub("['\"]", "", unlist(m2)))
  }
  extra <- setdiff(unique(used), c(PUBLIC, "hydrofabric"))
  expect_identical(extra, character(0),
    label = paste0("R/ source uses packages outside the allowlist: ",
                   paste(extra, collapse = ", ")))
})

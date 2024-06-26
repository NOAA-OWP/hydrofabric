% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sqlite.R
\name{as_sqlite}
\alias{as_sqlite}
\title{Read gpkg as SQLite}
\usage{
as_sqlite(gpkg, lyr = NULL, ignore = "gpkg_|rtree_|sqlite_")
}
\arguments{
\item{gpkg}{path}

\item{lyr}{layer of GPKG. If NULL (default the database tables are printed)}

\item{ignore}{pattern for layers to be ignored description}
}
\value{
an S4 object that inherits from DBIConnection. This object is used to communicate with the database engine.
}
\description{
Read gpkg as SQLite
}

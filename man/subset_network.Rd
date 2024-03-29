% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset_network.R
\name{subset_network}
\alias{subset_network}
\title{Subset Hydrofabric Network}
\usage{
subset_network(
  id = NULL,
  comid = NULL,
  hl_uri = NULL,
  nldi_feature = NULL,
  xy = NULL,
  bbox = NULL,
  base_s3 = "s3://lynker-spatial/v20/",
  base_dir = NULL,
  lyrs = c("divides", "nexus", "flowpaths", "network", "hydrolocations",
    "reference_flowline", "reference_catchment", "refactored_flowpaths",
    "refactored_divides"),
  areasqkm = NULL,
  pathlengthkm = NULL,
  ms_pathlengthkm = NULL,
  outfile = NULL,
  cache_dir = NULL,
  qml_dir = system.file("qml", package = "hydrofabric"),
  cache_overwrite = FALSE
)
}
\arguments{
\item{id}{hydrofabric id (relevant only to nextgen fabrics)}

\item{comid}{NHDPlusV2 COMID}

\item{hl_uri}{hydrolocation URI (relevant only to nextgen fabrics)}

\item{nldi_feature}{list with names 'featureSource' and 'featureID' where 'featureSource' is derived from the "source" column of the response of dataRetrieval::get_nldi_sources() and the 'featureID' is a known identifier from the specified 'featureSource'.}

\item{xy}{Location given as vector of XY in CRS 4326 (long, lat)}

\item{bbox}{a numeric vector of length four, with xmin, ymin, xmax and ymax values}

\item{base_s3}{the base hydrofabric directory to access in Lynker's s3}

\item{base_dir}{the base hydrofabric directory}

\item{lyrs}{layers to extract. Default is all possible in the hydrofabric GPKG data model}

\item{areasqkm}{the maximum upstream area to subset}

\item{pathlengthkm}{the maximum upstream path length to subset}

\item{ms_pathlengthkm}{the maximum upstream mainstem path length to subset}

\item{outfile}{file path to write to. Must have ".gpkg" extension}

\item{cache_dir}{should data be cached to a local directory? Will speed up multiple subsets in the same region}

\item{cache_overwrite}{description. Should a cached file be overwritten}
}
\value{
file path (outfile) or list of features
}
\description{
Subset Hydrofabric Network
}

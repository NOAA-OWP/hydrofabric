subset <- function(
    id      = NULL,
    comid   = NULL,
    hl_id   = NULL,
    layers  = c("divides", "nexus", "flowpaths", "network", "hydrolocations"),
    version = c("pre-release", "v1.0")
) {
    layers  <- match.arg(layers, several.ok = TRUE)
    version <- match.arg(version)


    pattern <- paste(
        "https://nextgen-hydrofabric.s3.amazonaws.com",
        version,
        "nextgen_{vpu}.gpkg",
        sep = "/"
    )

    logger::log_info(glue::glue("[subset] Using pattern: {pattern}"))

    hf_tmp <- tempfile(fileext = ".gpkg")
    on.exit(unlink(hf_tmp))
    hydrofabric::subset_network(
        id = id,
        comid = comid,
        hl_id = hl_id,
        network = "/hydrofabric/data/conus_net.parquet",
        pattern = pattern,
        lyrs = layers,
        export_gpkg = hf_tmp
    )
    hf_data <- readr::read_file_raw(hf_tmp)

    stringr::str_remove_all(jsonlite::base64_enc(hf_data), "\n")
}

lambdr::start_lambda()

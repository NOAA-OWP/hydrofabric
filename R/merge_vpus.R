#' Merge per-VPU fabrics into domain products with globally unique IDs
#'
#' Assigns globally unique IDs across every VPU in `gpkg_map`, then writes one
#' GeoPackage per domain:
#' - `superconus_ngen.gpkg`: VPUs 01-18
#' - `ak_ngen.gpkg`: VPU 19
#' - `hi_ngen.gpkg`: VPU 20
#' - `prvi_ngen.gpkg`: VPU 21
#' - `guam_ngen.gpkg`: VPU 22, only when it is in `gpkg_map` (not part of the
#'   NWM v4 delivery)
#'
#' IDs are assigned in contiguous blocks in VPU order, so a flowpath id such as
#' `fp-450003` belongs to exactly one VPU in any domain.
#'
#' @param gpkg_map Named list mapping vpuid -> absolute path to `ngen_vpu*.gpkg`.
#'   All VPUs that should participate in the global ID pool must be present.
#'   Example: `list("01" = "/path/ngen_vpu01.gpkg", "19" = "/path/ngen_vpu19.gpkg")`
#' @param out_dir Directory where output gpkgs are written.
#' @param mods_csv Path to the cross-VPU COMID modifications CSV (CONUS only).
#'   Defaults to the vendored copy in `inst/extdata/rpu_vpu_out.csv`.
#' @param overwrite If `TRUE`, delete existing output files before writing.
#' @param strict If `TRUE` (default), the merge is fail-fast: unresolved
#'   cross-VPU topology and invariant failures abort the run.
#' @param verbose Emit per-layer progress messages.
#' @param grouping How to partition VPUs into output domains: `"vpu"` (default)
#'   groups by the leading integer of `vpuid` (country/VPU domain); `"tdx_continent"`
#'   groups by the leading digit of the TDX-Hydro HydroBASINS level-2 id (continent).
#' @param target_crs Output coordinate reference system for the merged domain
#'   products (default `"EPSG:4326"`).
#' @return Invisibly, a list with `outfiles`, `fp_map`, `poi_map`,
#'   `cross_vpu_fix`, `invariants` (one per domain product).
#' @examples
#' \dontrun{
#' gpkg_map <- list("01" = "out/vpu_01/ngen_vpu01.gpkg",
#'                  "19" = "out/vpu_19/ngen_vpu19.gpkg")
#' merge_vpus(gpkg_map, out_dir = "out/global", overwrite = TRUE)
#' }
#' @export
merge_vpus <- function(
    gpkg_map,
    out_dir,
    mods_csv  = system.file("extdata", "rpu_vpu_out.csv", package = "hydrofabric"),
    overwrite = FALSE,
    strict    = TRUE,
    verbose   = FALSE,
    grouping  = c("vpu", "tdx_continent"),
    target_crs = "EPSG:4326"
    ) {
  grouping <- match.arg(grouping)
  stopifnot(is.list(gpkg_map), length(gpkg_map) >= 1L, !is.null(names(gpkg_map)))

  term_add <- 1e9

  # Canonical VPU order for block assignment: determines which VPU's IDs
  # come first globally. Sort numerically by the leading integer part.
  all_vpus   <- names(gpkg_map)
  vpu_order  <- all_vpus[order(as.numeric(sub("[A-Za-z].*", "", all_vpus)),
    all_vpus)]

  # Domain membership: two PARALLEL, non-conflated schemes (one merge run each):
  #  - "vpu" (default): leading-INTEGER of the vpuid -> country/VPU domain
  #    (CONUS 1-18, AK 19, HI 20, PRVI 21, Guam 22); sub-VPU tiles (19a/19b) and
  #    lettered CONUS VPUs (03N) regroup via their numeric prefix. This is the
  #    NHDPlus-VPU + oCONUS-country product set.
  #  - "tdx_continent": leading-DIGIT of the TDX-Hydro HydroBASINS level-2 id ->
  #    continent (1 Africa .. 9 Greenland). Built below at Step 4. The 10-digit
  #    basin ids overflow 32-bit integers, so that path uses string ops only.
  # The two never mix: a continent run passes grouping="tdx_continent" with a
  # TDX gpkg_map; a country run uses the default with a VPU gpkg_map.
  if (grouping == "vpu") {
    vpu_num    <- as.integer(sub("[A-Za-z].*", "", vpu_order))
    conus_vpus <- vpu_order[!is.na(vpu_num) & vpu_num >= 1L & vpu_num <= 18L]
    ak_vpus    <- vpu_order[!is.na(vpu_num) & vpu_num == 19L]
    hi_vpus    <- vpu_order[!is.na(vpu_num) & vpu_num == 20L]
    prvi_vpus  <- vpu_order[!is.na(vpu_num) & vpu_num == 21L]
    guam_vpus  <- vpu_order[!is.na(vpu_num) & vpu_num == 22L]   # Guam/Marianas (HUC 22)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # -- Step 1: Validate --------------------------------------------------------
  req_layers <- c("flowpaths", "divides", "nexus", "network")
  # pois / hydrolocations are absent for gauge-less fabrics (most of the global,
  # non-US world: TDX basins with no NWIS gages). Optional; re-ID is guarded.
  opt_layers <- c("pois", "hydrolocations", "flowlines",
    "flowpath-attributes", "flowline-attributes", "lakes")

  for (vpu in vpu_order) {
    g     <- gpkg_map[[vpu]]
    if (!file.exists(g))
      stop(sprintf("VPU %s: file not found: %s", vpu, g))
    avail <- sf::st_layers(g)$name
    miss  <- setdiff(req_layers, avail)
    if (length(miss) > 0L)
      stop(sprintf("VPU %s: required layer(s) missing: %s",
        vpu, paste(miss, collapse = ", ")))
  }

  # -- Step 2: Build global integer mappings across ALL VPUs ------------------
  message(sprintf("[stage4] Building global ID pool across %d VPUs ...",
    length(vpu_order)))

  fp_pool_list  <- vector("list", length(vpu_order))
  poi_pool_list <- vector("list", length(vpu_order))
  lnd_pool_list <- vector("list", length(vpu_order))

  for (i in seq_along(vpu_order)) {
    vpu <- vpu_order[[i]]
    g   <- gpkg_map[[vpu]]
    fps <- hfutils::as_ogr(g, "flowpaths") |>
      dplyr::select(flowpath_id, flowpath_toid, poi_id) |>
      dplyr::collect()

    fp_int  <- as.integer(sub("^fp-",  "", fps$flowpath_id))
    nex_int <- suppressWarnings(as.integer(sub("^nex-", "", fps$flowpath_toid)))
    # Large-int terminal nexuses use "tnx-" prefix; include their integers too.
    tnx_int <- suppressWarnings(as.integer(sub("^tnx-", "", fps$flowpath_toid)))
    pool    <- sort(unique(stats::na.omit(c(
      fp_int,
      nex_int[!is.na(nex_int) & nex_int >= term_add],
      tnx_int[!is.na(tnx_int) & tnx_int >= term_add]
    ))))

    fp_pool_list[[i]]  <- data.frame(vpuid = vpu, old_int = pool,
      stringsAsFactors = FALSE)

    # Pool POI ids from the hydrolocations layer as well as the flowpaths layer.
    #
    # A POI the network cannot place at an outlet is withdrawn from `pois` and its
    # poi_id cleared from `flowpaths`, while the hydrolocation is kept with a
    # `poi_drop_reason`. Pooling from flowpaths alone leaves those ids out of
    # poi_map, and the 4f re-ID then maps them to NA: silently stripping the
    # poi_id off exactly the hydrolocations whose lineage the withdrawal exists to
    # preserve. Both sources are unioned so every id that survives anywhere in the
    # product gets a global integer.
    poi_hl <- tryCatch(
      suppressWarnings(as.integer(hfutils::as_ogr(g, "hydrolocations") |>
        dplyr::select(poi_id) |> dplyr::collect() |> (\(x) x$poi_id)())),
      error = function(e) integer(0))
    poi_int <- sort(unique(stats::na.omit(c(
      suppressWarnings(as.integer(fps$poi_id)), poi_hl))))
    # gauge-less fabrics have no pois -> empty; skip (recycling vpu against 0 rows
    # errors), same as the lnd pool below.
    poi_pool_list[[i]] <- if (length(poi_int))
      data.frame(vpuid = vpu, old_int = poi_int, stringsAsFactors = FALSE) else NULL

    # Landscape divides (lnd-) are appended post-aggregation and are NOT keyed by
    # a flowpath_id, so they miss the fp pool. Pool them separately so they too
    # get globally-unique ids (lnd- and cat- share an int space only within their
    # own prefix, so an independent lnd pool is sufficient).
    lnd_int <- tryCatch(
      sort(unique(stats::na.omit(suppressWarnings(as.integer(sub("^lnd-", "",
        grep("^lnd-", hfutils::as_ogr(g, "divides") |>
          dplyr::select(divide_id) |> dplyr::collect() |>
          (\(x) x$divide_id)(), value = TRUE))))))),
      error = function(e) integer(0))
    lnd_pool_list[[i]] <- if (length(lnd_int))
      data.frame(vpuid = vpu, old_int = lnd_int, stringsAsFactors = FALSE) else NULL
  }

  # Assign contiguous global integer blocks in VPU order
  fp_all  <- do.call(rbind, fp_pool_list)
  fp_map  <- fp_all |>
    dplyr::mutate(new_int = seq_len(dplyr::n()))

  poi_all <- do.call(rbind, poi_pool_list)
  poi_map <- if (!is.null(poi_all) && nrow(poi_all) > 0L)
    dplyr::mutate(poi_all, new_int = seq_len(dplyr::n())) else
    data.frame(vpuid = character(0), old_int = integer(0), new_int = integer(0))

  lnd_all <- do.call(rbind, lnd_pool_list)
  lnd_map <- if (!is.null(lnd_all) && nrow(lnd_all) > 0L)
    dplyr::mutate(lnd_all, new_int = seq_len(dplyr::n())) else
    data.frame(vpuid = character(0), old_int = integer(0), new_int = integer(0))

  # Lookup helper: vectorised (vpuid, old_int) -> new_int
  .map_ints <- function(old_ints, vpuids, map) {
    key <- paste(vpuids, old_ints, sep = "|||")
    ref <- paste(map$vpuid, map$old_int, sep = "|||")
    map$new_int[match(key, ref)]
  }

  # Safe paste: preserves NA inputs as NA instead of "prefix-NA"
  .pfx <- function(prefix, new_ints) {
    dplyr::if_else(is.na(new_ints), NA_character_, paste0(prefix, new_ints))
  }

  .div_pfx <- function(x) {
    out <- rep(NA_character_, length(x))
    out[!is.na(x) & startsWith(x, "cat-")] <- "cat-"
    out[!is.na(x) & startsWith(x, "lnd-")] <- "lnd-"
    out
  }

  message(sprintf("[stage4] Global pool: %d fp/nex integers, %d poi integers",
    nrow(fp_map), nrow(poi_map)))

  # -- Step 3: Cross-VPU topology fix ----------------------------------------
  # Re-stitch flow across boundaries severed during tiling, driven by a seam
  # table (VPUID,toVPUID,COMID,toCOMID) matched to global flowpaths via
  # member_comid. Grouping-agnostic: the NHDPlus VPU product feeds NHD's
  # inter-VPU tocomid table; the TDX continent product feeds the mega-basin
  # confluence cuts emitted by rf.subvpu.tile(split_large=TRUE). The seam VPUs
  # are read from the table itself (intersected with what we're merging), so no
  # domain is hardcoded.
  message("[stage4] Building cross-VPU topology fix ...")
  cross_vpu_fix <- data.frame(from_fp_int = integer(0), to_fp_int = integer(0))

  if (nzchar(mods_csv) && file.exists(mods_csv)) {
    mods <- utils::read.csv(mods_csv, stringsAsFactors = FALSE) |>
      dplyr::mutate(VPUID = as.character(VPUID), toVPUID = as.character(toVPUID)) |>
      dplyr::filter(!is.na(VPUID), !is.na(toVPUID), VPUID != toVPUID,
        VPUID %in% vpu_order, toVPUID %in% vpu_order) |>
      dplyr::mutate(from_comid = as.character(COMID),
        to_comid   = as.character(toCOMID))
    seam_vpus <- intersect(union(mods$VPUID, mods$toVPUID), vpu_order)

    if (nrow(mods) > 0L && length(seam_vpus) > 1L) {
      # Build member_comid -> global fp int across all seam-involved VPUs
      mc_rows <- lapply(seam_vpus, function(vpu) {
        g   <- gpkg_map[[vpu]]
        fps <- hfutils::as_ogr(g, "flowpaths") |>
          dplyr::select(flowpath_id, member_comid) |>
          dplyr::collect()
        fps$vpuid  <- vpu
        fps$fp_int <- as.integer(sub("^fp-", "", fps$flowpath_id))
        fps
      })
      mc_df <- do.call(rbind, mc_rows) |>
        dplyr::left_join(fp_map |> dplyr::select(vpuid, old_int, new_int),
          by = c("vpuid", fp_int = "old_int")) |>
        dplyr::rename(new_fp_int = new_int) |>
        .separate_longer("member_comid", ",") |>
        dplyr::mutate(member_comid = trimws(member_comid)) |>
        dplyr::filter(nzchar(member_comid))

      # Match each seam endpoint WITHIN ITS OWN tile (vpuid), not just by COMID.
      # Reference ids are only basin-unique (every basin is offset by 1e9, so
      # 1000xxxxxx values collide ACROSS basins); a COMID-only join lets a seam
      # false-match a flowpath in another basin, minting spurious cross-VPU
      # fixes. Scoping by vpuid makes each match 1:1.
      from_new <- dplyr::left_join(
        data.frame(member_comid = mods$from_comid, vpuid = mods$VPUID, stringsAsFactors = FALSE),
        mc_df |> dplyr::select(member_comid, vpuid, new_fp_int),
        by = c("member_comid", "vpuid")
      )$new_fp_int

      to_new <- dplyr::left_join(
        data.frame(member_comid = mods$to_comid, vpuid = mods$toVPUID, stringsAsFactors = FALSE),
        mc_df |> dplyr::select(member_comid, vpuid, new_fp_int),
        by = c("member_comid", "vpuid")
      )$new_fp_int

      cross_vpu_fix <- data.frame(from_fp_int = from_new, to_fp_int = to_new,
        stringsAsFactors = FALSE) |>
        dplyr::filter(!is.na(from_fp_int), !is.na(to_fp_int)) |>
        dplyr::distinct()

      message(sprintf("[stage4] %d cross-VPU topology fix(es)", nrow(cross_vpu_fix)))
    }
  }

  # -- Step 4: Re-ID and write each domain product ---------------------------
  if (grouping == "tdx_continent") {
    # STANDARD 7-CONTINENT products. HydroBASINS level-1 regions (the leading
    # digit of the 10-digit basin id) fold to standard continents:
    #   Asia          = Siberia(3) + Asia(4)
    #   North America = NorthAmerica(7) + Arctic/Alaska(8) + Greenland(9)
    #   Africa(1) Europe(2) Oceania(5) South America(6) map 1:1
    # (Antarctica has no TDX-Hydro network.) Recycled TDX oCONUS builds carry
    # short VPU-style ids, their fabric is already TDX (the AK recipe = TDX +
    # GEDTM30), so they slot straight in: AK(19)/PRVI(21) -> North America,
    # HI(20)/Guam(22) -> Oceania. String ops only (10-digit ids overflow int32).
    .std_continent <- function(vpuid) {
      base <- sub("[A-Za-z].*$", "", vpuid)
      if (nchar(base) <= 2L)                       # oCONUS VPU-style id (recycled TDX)
        return(switch(base, "19" = "northamerica", "20" = "oceania",
          "21" = "northamerica", "22" = "oceania", NA_character_))
      c("1"="africa", "2"="europe", "3"="asia", "4"="asia", "5"="oceania",
        "6"="southamerica", "7"="northamerica", "8"="northamerica",
        "9"="northamerica")[[substr(base, 1L, 1L)]]
    }
    continent <- vapply(vpu_order, .std_continent, character(1L))
    na_ix     <- is.na(continent)
    if (any(na_ix))
      continent[na_ix] <- paste0("region_", substr(sub("[A-Za-z].*", "", vpu_order[na_ix]), 1L, 1L))
    domain_specs <- list()
    for (nm in sort(unique(continent)))
      domain_specs[[nm]] <- list(vpus = vpu_order[continent == nm],
        outfile = file.path(out_dir, paste0(nm, "_ngen.gpkg")))
  } else {
    domain_specs <- list(
      superconus = list(vpus = conus_vpus, outfile = file.path(out_dir, "superconus_ngen.gpkg")),
      ak         = list(vpus = ak_vpus,    outfile = file.path(out_dir, "ak_ngen.gpkg")),
      prvi       = list(vpus = prvi_vpus,  outfile = file.path(out_dir, "prvi_ngen.gpkg"))
    )
    # HI is a single-VPU product, its IDs are globally unique after re-ID but
    # it is written in-place (the existing ngen_vpu20.gpkg is updated).
    if (length(hi_vpus) > 0L) {
      domain_specs$hi <- list(vpus = hi_vpus,
        outfile = file.path(out_dir, "hi_ngen.gpkg"))
    }
    if (length(guam_vpus) > 0L) {
      domain_specs$guam <- list(vpus = guam_vpus,
        outfile = file.path(out_dir, "guam_ngen.gpkg"))
    }
  }

  # Drop empty domains (e.g. if AK or PRVI not supplied)
  domain_specs <- domain_specs[vapply(domain_specs, function(d) length(d$vpus) > 0L,
    logical(1L))]

  # Pre-cache which (gpkg, layer) pairs have a vpuid column.
  # Avoids reading full tables inside the per-VPU merge loop.
  # Uses ogrinfo via system2 (fast metadata-only query).
  .has_vpuid <- local({
    cache <- new.env(hash = TRUE, parent = emptyenv())
    function(gpkg_path, lyr) {
      key <- paste0(gpkg_path, "|||", lyr)
      if (!exists(key, envir = cache)) {
        cols <- tryCatch(
          {
            out <- system2("ogrinfo", c("-al", "-so",
              shQuote(gpkg_path), lyr),
            stdout = TRUE, stderr = FALSE)
            # Field lines: "fieldname: Type (width)": start at column 0 with
            # a lowercase letter or underscore, followed by ": ".
            flds <- grep("^[a-z_][A-Za-z0-9_]*:", out, value = TRUE)
            sub(":.*", "", flds)
          },
          error = function(e) character(0L))
        assign(key, "vpuid" %in% cols, envir = cache)
      }
      get(key, envir = cache)
    }
  })

  all_results <- list()
  global_fl_base <- 0L   # cumulative flowline count across prior domains (global flowline_id offset)

  for (domain in names(domain_specs)) {
    ngen_net_inv <- NULL
    spec  <- domain_specs[[domain]]
    vpus  <- spec$vpus
    ofile <- spec$outfile

    if (file.exists(ofile)) {
      if (!isTRUE(overwrite))
        stop("outfile already exists (use overwrite = TRUE): ", ofile)
      file.remove(ofile)
    }

    message(sprintf("\n[stage4][%s] Writing %d VPU(s) to %s",
      domain, length(vpus), basename(ofile)))

    # -- 4a: Merge raw layers via ogr2ogr (fast GDAL-native append) -----------
    # ogr2ogr uses the first file's schema as template; columns absent in
    # subsequent files are filled with NULL. VPU 01 goes first for CONUS so its
    # fuller schema (hf_id / total_dasqkm) supplies the template columns.
    #
    # Some source VPUs lack a vpuid column in their spatial layers. We inject it
    # via -sql so the re-ID lookup (which keys on vpuid) works after the merge.
    all_lyr <- unique(c(req_layers, opt_layers))

    # Source gpkg layer names are static across the merge; query each once
    # instead of re-opening every source file inside the per-layer loop.
    src_layers <- lapply(gpkg_map, function(g) sf::st_layers(g)$name)

    for (lyr in all_lyr) {
      has_lyr <- Filter(function(vpu) lyr %in% src_layers[[vpu]], vpus)
      if (length(has_lyr) == 0L) next

      # First VPU of first layer: plain create (no -update, establishes the file).
      # All other writes: -update -append [-addfields] to add rows/columns without
      # destroying existing layers. Track whether ofile has been created yet.
      for (i in seq_along(has_lyr)) {
        vpu_i  <- has_lyr[[i]]
        src    <- gpkg_map[[vpu_i]]

        file_exists <- file.exists(ofile)
        # Layer exists in ofile already (earlier layer wrote it, or prior VPU)?
        lyr_exists  <- file_exists &&
          (lyr %in% sf::st_layers(ofile)$name)

        # Detect whether source layer has vpuid; if not, inject it via -sql.
        needs_vpuid_inject <- !.has_vpuid(src, lyr)

        # When vpuid must be injected, use -sql (which replaces layer selection).
        # When not needed, pass the layer name positionally after dst/src.
        # ogr2ogr syntax: ogr2ogr [opts] <dst> <src> [layer]
        # -sql and -dialect must precede dst/src.
        if (needs_vpuid_inject) {
          sql_clause  <- sprintf("SELECT *, '%s' AS vpuid FROM \"%s\"", vpu_i, lyr)
          # Also rename the result layer to the canonical layer name via -nln
          inject_opts <- c("-sql", shQuote(sql_clause),
            "-dialect", "OGRSQL",
            "-nln", lyr)
          lyr_arg     <- character(0L)  # no positional layer when using -sql
        } else {
          inject_opts <- character(0L)
          lyr_arg     <- lyr
        }

        # Reproject spatial layers to the common target CRS on append. Per-tile
        # builds live in local UTM zones (14N/15N/18N...); without -t_srs,
        # ogr2ogr -append concatenates raw coordinates under the first tile's CRS
        # label, displacing every off-zone tile. Aspatial layers carry no
        # geometry -> skip.
        # pois and hydrolocations belong on this list. Left off it they kept the
        # per-VPU projected CRS while the network layers were reprojected, so every
        # domain product shipped two coordinate systems, and those two are exactly
        # the layers a consumer joins to the network. The check that should have
        # caught it was reading only three of the seven spatial layers.
        srs_opt <- if (lyr %in% c("flowpaths", "divides", "nexus", "flowlines",
                                  "lakes", "pois", "hydrolocations"))
          c("-t_srs", target_crs) else character(0L)

        args <- if (!file_exists) {
          c(inject_opts, srs_opt, "-f", "GPKG", shQuote(ofile), shQuote(src), lyr_arg)
        } else if (!lyr_exists) {
          c(inject_opts, srs_opt, "-f", "GPKG", "-update",
            shQuote(ofile), shQuote(src), lyr_arg)
        } else {
          c(inject_opts, srs_opt, "-f", "GPKG", "-update", "-append", "-addfields",
            shQuote(ofile), shQuote(src), lyr_arg)
        }

        ret    <- system2("ogr2ogr", args, stdout = TRUE, stderr = TRUE)
        status <- attr(ret, "status")
        if (!is.null(status) && status != 0L)
          stop(sprintf("[stage4][%s] ogr2ogr failed for %s VPU %s:\n%s",
            domain, lyr, vpu_i, paste(ret, collapse = "\n")))
      }
      if (verbose)
        message(sprintf("[stage4][%s]   merged: %s (%d VPUs)", domain, lyr,
          length(has_lyr)))
    }

    # -- 4b: Re-ID flowpaths -------------------------------------------------
    fps <- sf::st_read(ofile, "flowpaths", quiet = TRUE)
    fp_int      <- as.integer(sub("^fp-",  "", fps$flowpath_id))
    # flowpath_toid can be "nex-N" (regular) or "tnx-N" (terminal nexus)
    toid_pfx    <- dplyr::if_else(startsWith(fps$flowpath_toid, "tnx-"), "tnx-", "nex-")
    nex_int     <- as.integer(sub("^(nex|tnx)-", "", fps$flowpath_toid, perl = TRUE))
    div_int <- as.integer(sub("^(cat-|lnd-)", "",
      dplyr::if_else(is.na(fps$divide_id),
        NA_character_, fps$divide_id),
      perl = TRUE))
    div_pfx <- .div_pfx(fps$divide_id)
    poi_int <- suppressWarnings(as.integer(fps$poi_id))

    new_fp  <- .map_ints(fp_int,  fps$vpuid, fp_map)
    new_nex <- .map_ints(nex_int, fps$vpuid, fp_map)
    new_div <- .map_ints(div_int, fps$vpuid, fp_map)
    new_poi <- .map_ints(poi_int, fps$vpuid, poi_map)
    # mainstem_id is a reference levelpath id (a global, cross-VPU id space),
    # NOT a minted per-VPU flowpath id, so it is carried through unchanged.
    # Remapping it via fp_map would collide it onto unrelated flowpath ids and
    # fragment every mainstem.

    if (nrow(cross_vpu_fix) > 0L) {
      # Global fp ints are unique across all VPUs, so a match only fires in the
      # domain that actually owns the seam's upstream flowpath (superconus for
      # the NHDPlus product, the relevant continent for the TDX product).
      fix_ix    <- match(new_fp, cross_vpu_fix$from_fp_int)
      valid_fix <- !is.na(fix_ix)
      new_nex[valid_fix] <- cross_vpu_fix$to_fp_int[fix_ix[valid_fix]]
      # A cross-VPU-fixed flowpath now drains into an interior downstream nexus
      # across the VPU boundary, so it is no longer a within-VPU terminal: its
      # toid prefix must be "nex-" (the minted/interior nexus), not the "tnx-"
      # it carried as a terminal. Without this the toid points at a non-existent
      # tnx-<id> and fails the unresolved-toid integrity check.
      toid_pfx[valid_fix] <- "nex-"
      fps$ibt <- valid_fix
    }

    fps$flowpath_id   <- paste0("fp-",  new_fp)
    fps$flowpath_toid <- paste0(toid_pfx, new_nex)
    fps$divide_id     <- dplyr::if_else(!is.na(new_div),
      paste0(div_pfx, new_div), NA_character_)
    fps$poi_id        <- new_poi
    fps$mainstem_id   <- as.numeric(fps$mainstem_id)
    sf::st_write(fps, ofile, "flowpaths", delete_layer = TRUE, quiet = TRUE)

    # -- 4c: Re-ID nexus -----------------------------------------------------
    nex <- sf::st_read(ofile, "nexus", quiet = TRUE)
    # nexus_id can be "nex-N" (regular) or "tnx-N" (terminal nexus: preserve prefix)
    nex_id_pfx  <- dplyr::if_else(startsWith(nex$nexus_id, "tnx-"), "tnx-", "nex-")
    nex_id_int  <- as.integer(sub("^(nex|tnx)-", "", nex$nexus_id, perl = TRUE))
    nex$nexus_id <- paste0(nex_id_pfx, .map_ints(nex_id_int, nex$vpuid, fp_map))
    # nexus_toid is "fp-N" or NA (terminal); preserve NA instead of producing "fp-NA"
    nex_toid_int <- .map_ints(
      as.integer(sub("^fp-", "", nex$nexus_toid)), nex$vpuid, fp_map)
    nex$nexus_toid <- dplyr::if_else(!is.na(nex_toid_int),
      paste0("fp-", nex_toid_int), NA_character_)

    if (domain == "superconus" && nrow(cross_vpu_fix) > 0L) {
      needed_nex  <- paste0("nex-", cross_vpu_fix$to_fp_int)
      missing_nex <- setdiff(needed_nex, nex$nexus_id)
      if (length(missing_nex) > 0L) {
        # A cross-VPU fix points an upstream boundary flowpath at nex-<to>, where
        # fp-<to> is the downstream flowpath. When fp-<to> is a HEADWATER in its
        # own VPU (its upstream area lives across the VPU boundary) it never had
        # an inlet nexus, so nex-<to> is absent. Mint it: a nexus that drains
        # INTO fp-<to>, located at fp-<to>'s inlet (upstream endpoint), completing
        # the chain fp-<from> -> nex-<to> -> fp-<to>. If fp-<to> itself is absent
        # the reference is truly dangling -> strict stop (a real integrity error).
        to_fids <- paste0("fp-", as.integer(sub("^nex-", "", missing_nex)))
        fp_ix   <- match(to_fids, fps$flowpath_id)
        truly_missing <- missing_nex[is.na(fp_ix)]
        if (length(truly_missing) > 0L) {
          msg <- sprintf(
            "[stage4][%s] %d cross-VPU nexus target flowpath(s) absent (e.g. %s)",
            domain, length(truly_missing),
            paste(utils::head(truly_missing, 10), collapse = ", "))
          if (isTRUE(strict)) stop(msg) else .hf_log_warn(msg)
        }
        ok <- !is.na(fp_ix)
        if (any(ok)) {
          gcol  <- attr(nex, "sf_column")
          inlet <- lwgeom::st_startpoint(
            sf::st_line_merge(sf::st_geometry(fps)[fp_ix[ok]]))
          minted <- sf::st_sf(
            nexus_id   = missing_nex[ok],
            nexus_toid = to_fids[ok],
            vpuid      = as.character(fps$vpuid[fp_ix[ok]]),
            geometry   = sf::st_sfc(inlet, crs = sf::st_crs(nex)),
            stringsAsFactors = FALSE)
          names(minted)[names(minted) == "geometry"] <- gcol
          sf::st_geometry(minted) <- gcol
          nex <- rbind(nex[, names(minted)], minted)
          message(sprintf(
            "[stage4][%s] minted %d cross-VPU inlet nexus (e.g. %s -> %s)",
            domain, sum(ok), missing_nex[ok][1], to_fids[ok][1]))
        }
      }
    }
    sf::st_write(nex, ofile, "nexus", delete_layer = TRUE, quiet = TRUE)

    # -- 4d: Re-ID divides ---------------------------------------------------
    div <- sf::st_read(ofile, "divides", quiet = TRUE)
    ddiv_int <- as.integer(sub("^(cat-|lnd-)", "", div$divide_id, perl = TRUE))
    ddiv_pfx <- .div_pfx(div$divide_id)
    is_lnd   <- !is.na(div$divide_id) & startsWith(div$divide_id, "lnd-")
    # cat- divides re-ID via the global fp pool; lnd- landscape divides via their
    # own global lnd pool (so they are globally unique across domains too).
    new_div_int <- .map_ints(ddiv_int, div$vpuid, fp_map)
    new_lnd_int <- .map_ints(ddiv_int, div$vpuid, lnd_map)
    new_div_int[is_lnd] <- new_lnd_int[is_lnd]
    # A `cat-` divide whose own integer names no flowpath in its VPU cannot be
    # mapped from the pool. Falling through with the LOCAL id silently collides
    # with another VPU's remapped divide, so fall back to the integer of the
    # flowpath that owns it, the authoritative link, and say so when it
    # happens, since it means an upstream step broke the cat-N/fp-N convention.
    own_new <- .map_ints(
      suppressWarnings(as.integer(sub("^fp-", "", div$flowpath_id))),
      div$vpuid, fp_map)
    need <- !is_lnd & is.na(new_div_int) & !is.na(own_new)
    if (any(need)) {
      .hf_log_warn(sprintf(paste0("[stage4][%s] %d divide_id(s) do not name a flowpath ",
        "in their VPU; re-IDing via their owner instead (e.g. %s)"),
        domain, sum(need), paste(utils::head(div$divide_id[need], 3),
          collapse = ", ")))
      new_div_int[need] <- own_new[need]
    }
    mapped_ok   <- !is.na(new_div_int)
    if (any(!mapped_ok & !is.na(div$divide_id)))
      .hf_log_warn(sprintf(paste0("[stage4][%s] %d divide_id(s) resolved by no pool and ",
        "keep a VPU-local id; they may collide across VPUs"),
        domain, sum(!mapped_ok & !is.na(div$divide_id))))
    div$divide_id <- dplyr::if_else(
      mapped_ok,
      paste0(ddiv_pfx, new_div_int),
      div$divide_id)
    div$flowpath_id <- .pfx("fp-", .map_ints(
      as.integer(sub("^fp-", "", div$flowpath_id)), div$vpuid, fp_map))
    sf::st_write(div, ofile, "divides", delete_layer = TRUE, quiet = TRUE)

    avail <- sf::st_layers(ofile)$name

    # -- 4e: Re-ID pois (absent for gauge-less fabrics) ----------------------
    if ("pois" %in% avail) {
      poi <- sf::st_read(ofile, "pois", quiet = TRUE)
      poi$poi_id      <- .map_ints(as.integer(poi$poi_id), poi$vpuid, poi_map)
      poi$flowpath_id <- .pfx("fp-", .map_ints(
        as.integer(sub("^fp-", "", poi$flowpath_id)), poi$vpuid, fp_map))
      sf::st_write(poi, ofile, "pois", delete_layer = TRUE, quiet = TRUE)
    }

    fl_offsets <- NULL  # populated in 4g; reused in 4h and 4j
    dom_fl_span <- 0L   # this domain's total flowline count (advances global_fl_base)

    # -- 4f: Re-ID hydrolocations --------------------------------------------
    if ("hydrolocations" %in% avail) {
      hl <- sf::st_read(ofile, "hydrolocations", quiet = TRUE)
      if (nrow(hl) > 0L) {
        hl$poi_id      <- .map_ints(as.integer(hl$poi_id), hl$vpuid, poi_map)
        hl$flowpath_id <- .pfx("fp-", .map_ints(
          as.integer(sub("^fp-", "", hl$flowpath_id)), hl$vpuid, fp_map))
        sf::st_write(hl, ofile, "hydrolocations", delete_layer = TRUE, quiet = TRUE)
      }
    }

    # -- 4g: Re-ID flowlines -------------------------------------------------
    if ("flowlines" %in% avail) {
      fl <- sf::st_read(ofile, "flowlines", quiet = TRUE)
      if (nrow(fl) > 0L) {
        # Ensure integer arithmetic (some domains store flowline_id as character)
        fl$flowline_id   <- as.integer(fl$flowline_id)
        fl$flowline_toid <- as.integer(fl$flowline_toid)
        fl_vpus    <- sort(unique(fl$vpuid))
        fl_counts  <- tapply(fl$flowline_id, fl$vpuid,
          function(x) max(x, na.rm = TRUE))
        fl_offsets <- c(0L, cumsum(as.integer(fl_counts[fl_vpus[-length(fl_vpus)]])))
        names(fl_offsets) <- fl_vpus   # fl_offsets also used in 4h and 4j
        # Cross-domain global offset: shift this domain's flowline ids past all
        # prior domains so flowline_id is globally unique across domain products.
        fl_offsets  <- fl_offsets + global_fl_base
        dom_fl_span <- sum(as.integer(fl_counts), na.rm = TRUE)
        fl$flowline_id   <- fl$flowline_id + fl_offsets[fl$vpuid]
        term_mask        <- !is.na(fl$flowline_toid) & fl$flowline_toid == 0L
        fl$flowline_toid[!term_mask] <- fl$flowline_toid[!term_mask] +
          fl_offsets[fl$vpuid[!term_mask]]
        fl$flowpath_id   <- .pfx("fp-", .map_ints(
          as.integer(sub("^fp-",  "", fl$flowpath_id)),  fl$vpuid, fp_map))
        fl_toid_pfx  <- dplyr::if_else(startsWith(fl$flowpath_toid, "tnx-"), "tnx-", "nex-")
        fl_toid_int  <- as.integer(sub("^(nex|tnx)-", "", fl$flowpath_toid, perl = TRUE))
        fl$flowpath_toid <- dplyr::if_else(
          is.na(fl$flowpath_toid), NA_character_,
          .pfx(fl_toid_pfx, .map_ints(fl_toid_int, fl$vpuid, fp_map)))
        sf::st_write(fl, ofile, "flowlines", delete_layer = TRUE, quiet = TRUE)
      }
    }

    # -- 4h: Re-ID network ---------------------------------------------------
    net <- sf::st_read(ofile, "network", quiet = TRUE)
    net_div_int <- as.integer(sub("^(cat-|lnd-)", "", net$divide_id, perl = TRUE))
    net_div_pfx <- .div_pfx(net$divide_id)
    net$flowpath_id   <- .pfx("fp-", .map_ints(
      as.integer(sub("^fp-",  "", net$flowpath_id)),  net$vpuid, fp_map))
    # Re-ID toid from per-VPU data first (preserves tnx-/nex- prefix), then
    # sync with the already-written flowpaths layer as the authority.  This
    # reconciles cross-VPU connections and any other per-VPU inconsistencies
    # where the network table had different toid values than the flowpaths layer.
    net_toid_pfx  <- dplyr::if_else(startsWith(net$flowpath_toid, "tnx-"), "tnx-", "nex-")
    net_toid_int  <- as.integer(sub("^(nex|tnx)-", "", net$flowpath_toid, perl = TRUE))
    net$flowpath_toid <- dplyr::if_else(
      is.na(net$flowpath_toid), NA_character_,
      .pfx(net_toid_pfx, .map_ints(net_toid_int, net$vpuid, fp_map)))
    # Sync network flowpath_toid from the authoritative flowpaths layer
    fps_lookup <- local({
      # attribute-only read: sf::st_read would materialize every geometry
      # (~40GB on the largest continents) only to drop it immediately.
      .c <- DBI::dbConnect(RSQLite::SQLite(), ofile)
      on.exit(DBI::dbDisconnect(.c))
      DBI::dbGetQuery(.c, "SELECT flowpath_id, flowpath_toid FROM flowpaths")
    })
    fp_toid_vec <- stats::setNames(fps_lookup$flowpath_toid, fps_lookup$flowpath_id)
    has_fp <- !is.na(net$flowpath_id)
    net$flowpath_toid[has_fp] <- fp_toid_vec[net$flowpath_id[has_fp]]
    net_is_lnd     <- !is.na(net$divide_id) & startsWith(net$divide_id, "lnd-")
    net_new_div    <- .map_ints(net_div_int, net$vpuid, fp_map)
    net_new_lnd    <- .map_ints(net_div_int, net$vpuid, lnd_map)
    net_new_div[net_is_lnd] <- net_new_lnd[net_is_lnd]
    net$divide_id     <- dplyr::if_else(
      is.na(net$divide_id), NA_character_,
      .pfx(net_div_pfx, net_new_div))
    net$poi_id        <- .map_ints(suppressWarnings(as.integer(net$poi_id)),
      net$vpuid, poi_map)
    # Re-ID flowline_id / flowline_toid with the same per-VPU integer offsets used
    # for the flowlines layer (section 4g).  Without this, the network parquet
    # retains VPU-local integers that collide across VPUs and cause downstream subsetting's
    # graph traversal to find incorrect upstream sets.
    if (!is.null(fl_offsets) && "vpuid" %in% names(net)) {
      net$flowline_id   <- as.integer(net$flowline_id) + fl_offsets[net$vpuid]
      term_mask_net     <- !is.na(net$flowline_toid) & net$flowline_toid == 0
      net$flowline_toid[!term_mask_net] <- net$flowline_toid[!term_mask_net] +
        fl_offsets[net$vpuid[!term_mask_net]]
    }
    net <- dplyr::select(net, -dplyr::any_of(c("hydroseq", "total_dasqkm")))

    # Keep the result. This ran and was logged but never recorded, so the
    # invariants.json beside a merged product described only the merge suite,
    # the flowline id/toid/DAG checks and the POI contract were absent from the
    # audit trail of a build that had in fact run them.
    ngen_net_inv <- check_invariants("ngen", network = net, strict = FALSE)

    sf::st_write(net, ofile, "network", delete_layer = TRUE, quiet = TRUE)

    # -- 4i: Re-ID flowpath-attributes ---------------------------------------
    if ("flowpath-attributes" %in% avail) {
      fa <- sf::st_read(ofile, "flowpath-attributes", quiet = TRUE)
      if (nrow(fa) > 0L) {
        fa$flowpath_id <- .pfx("fp-", .map_ints(
          as.integer(sub("^fp-", "", fa$flowpath_id)), fa$vpuid, fp_map))
        sf::st_write(fa, ofile, "flowpath-attributes",
          delete_layer = TRUE, quiet = TRUE)
      }
    }

    # -- 4j: Re-ID flowline-attributes ---------------------------------------
    if ("flowline-attributes" %in% avail) {
      fla <- sf::st_read(ofile, "flowline-attributes", quiet = TRUE)
      if (nrow(fla) > 0L) {
        # flowline_id is a sub-flowline integer (one row per segment, 86k+ per VPU).
        # Use the same per-VPU offsets as the flowlines layer, NOT fp_map.
        fla_vpuid <- if ("vpuid" %in% names(fla)) fla$vpuid else rep(vpus[[1]], nrow(fla))
        fla_fl_int <- as.integer(sub("^fp-|^fl-", "", fla$flowline_id))
        if (!is.null(fl_offsets)) {
          fla$flowline_id <- fla_fl_int + fl_offsets[fla_vpuid]
        } else {
          fla$flowline_id <- fla_fl_int
        }
        sf::st_write(fla, ofile, "flowline-attributes",
          delete_layer = TRUE, quiet = TRUE)
      }
    }

    # -- 4k: Re-ID lakes -----------------------------------------------------
    if ("lakes" %in% avail) {
      lk <- sf::st_read(ofile, "lakes", quiet = TRUE)
      if (nrow(lk) > 0L) {
        lk_vpuid <- if ("vpuid" %in% names(lk)) lk$vpuid else rep(vpus[[1]], nrow(lk))
        lk$flowpath_id <- paste0("fp-", .map_ints(
          as.integer(sub("^fp-", "", lk$flowpath_id)), lk_vpuid, fp_map))
        sf::st_write(lk, ofile, "lakes", delete_layer = TRUE, quiet = TRUE)
      }
    }

    # -- 4l: Domain-scale hydroseq + total_dasqkm ----------------------------
    message(sprintf("[stage4][%s] Recomputing hydroseq + total_dasqkm ...", domain))
    # Read flowpaths ATTRIBUTES ONLY (no geometry). sf::st_read materializes every
    # multilinestring in memory (~40GB for Africa's 1.74M fps) and OOM-kills the
    # largest continents; the recompute only needs id/toid/areasqkm. Results are
    # written back via SQL UPDATE below, leaving geometry untouched.
    .con <- DBI::dbConnect(RSQLite::SQLite(), ofile)
    fps <- DBI::dbGetQuery(.con,
      "SELECT flowpath_id, flowpath_toid, areasqkm, mainstem_id, lengthkm FROM flowpaths")

    # flowpath_toid points at a NEXUS (nex-/tnx-), not the downstream flowpath.
    # Resolve it to the downstream flowpath id via the nexus layer
    # (nexus_id -> nexus_toid), exactly as stage-3 finalize_data_model does.
    # Without this, hydroseq/accumulate_downstream see fp -> nexus edges where
    # the nexus is a sink: hydroseq disconnects and total_dasqkm collapses to
    # each flowpath's own areasqkm (no cross-reach accumulation).
    nx_l    <- DBI::dbGetQuery(.con, "SELECT nexus_id, nexus_toid FROM nexus")
    nex2fp  <- stats::setNames(as.character(nx_l$nexus_toid),
      as.character(nx_l$nexus_id))
    fp_down <- unname(nex2fp[as.character(fps$flowpath_toid)])  # downstream fp-id; NA at terminals

    fps$hydroseq <- local({
      id_chr <- as.character(fps$flowpath_id)
      el  <- data.frame(from = id_chr, to = fp_down, stringsAsFactors = FALSE)
      el  <- el[!is.na(el$to) & el$to %in% id_chr, ]
      # Include ALL flowpaths as vertices, not just those appearing in an edge.
      # Otherwise an isolated single-reach basin (terminal AND headwater, common
      # on islands like HI/PRVI) is absent from the edge-built graph, omitted by
      # topo_sort, and left with NA hydroseq.
      g   <- igraph::graph_from_data_frame(el, directed = TRUE,
        vertices = data.frame(name = id_chr, stringsAsFactors = FALSE))
      topo <- igraph::as_ids(igraph::topo_sort(g, mode = "out"))
      rank <- seq_along(topo)
      names(rank) <- topo
      rank[id_chr]
    })

    # Accumulate over the resolved fp -> downstream-fp topology.
    fps_acc <- fps
    fps_acc$flowpath_toid <- fp_down
    fps$total_dasqkm <- hfutils::accumulate_downstream(fps_acc, attr = "areasqkm")
    # lengthkm is already correct from the per-tile UTM build (carried attribute);
    # NOT recomputed: add_lengthkm() needs geometry (defeats the attr-only read)
    # and would compute geographic lengths on the EPSG:4326 merge.

    # Write hydroseq + total_dasqkm back via SQL UPDATE, leaving the flowpaths
    # geometry untouched. A full st_write(delete_layer) would drop and rewrite
    # every geometry and can OOM on the largest continents.
    .have <- DBI::dbListFields(.con, "flowpaths")
    if (!"hydroseq"     %in% .have) DBI::dbExecute(.con, "ALTER TABLE flowpaths ADD COLUMN hydroseq REAL")
    if (!"total_dasqkm" %in% .have) DBI::dbExecute(.con, "ALTER TABLE flowpaths ADD COLUMN total_dasqkm REAL")
    # GeoPackage rtree triggers on flowpaths call ST_IsEmpty (a Spatialite fn that
    # plain RSQLite lacks): SQLite compiles a table's triggers when preparing any
    # UPDATE on it, so the UPDATE fails ("no such function: ST_IsEmpty") even though
    # we change no geometry. Drop the triggers around the UPDATE and restore them
    # verbatim from sqlite_master; the rtree index stays valid (geometry untouched).
    .trg <- DBI::dbGetQuery(.con,
      "SELECT name, sql FROM sqlite_master WHERE type='trigger' AND tbl_name='flowpaths'")
    for (.nm in .trg$name) DBI::dbExecute(.con, sprintf('DROP TRIGGER IF EXISTS "%s"', .nm))
    DBI::dbWriteTable(.con, "_fp_upd", data.frame(
      fpid = as.character(fps$flowpath_id),
      hs   = as.numeric(fps$hydroseq),
      da   = as.numeric(fps$total_dasqkm),
      stringsAsFactors = FALSE), overwrite = TRUE)
    DBI::dbExecute(.con, "CREATE INDEX _fp_upd_ix ON _fp_upd(fpid)")
    DBI::dbExecute(.con, paste(
      "UPDATE flowpaths SET",
      "hydroseq     = (SELECT hs FROM _fp_upd WHERE _fp_upd.fpid = flowpaths.flowpath_id),",
      "total_dasqkm = (SELECT da FROM _fp_upd WHERE _fp_upd.fpid = flowpaths.flowpath_id)"))
    DBI::dbExecute(.con, "DROP TABLE _fp_upd")
    for (.s in .trg$sql) if (!is.na(.s) && nzchar(.s))
      tryCatch(DBI::dbExecute(.con, .s), error = function(e)
        .hf_log_warn(sprintf("[stage4][%s] could not restore trigger: %s", domain, conditionMessage(e))))
    DBI::dbDisconnect(.con)

    net <- sf::st_read(ofile, "network", quiet = TRUE)
    # network is flowline-keyed (multiple rows per flowpath_id); use a named
    # vector lookup instead of join to avoid cartesian fanout.
    fp_hydroseq    <- stats::setNames(fps$hydroseq,    fps$flowpath_id)
    fp_total_dasqkm <- stats::setNames(fps$total_dasqkm, fps$flowpath_id)
    net$hydroseq    <- fp_hydroseq[net$flowpath_id]
    net$total_dasqkm <- fp_total_dasqkm[net$flowpath_id]
    if ("mainstem_id" %in% names(net)) {
      fp_ms <- stats::setNames(fps$mainstem_id, fps$flowpath_id)
      net$mainstem_id <- fp_ms[net$flowpath_id]
    }
    sf::st_write(net, ofile, "network", delete_layer = TRUE, quiet = TRUE)

    # -- Global nested-set upstream index ------------------------------------
    # Recompute the index on the MERGED network so upstream_id / num_upstreams
    # are globally unique and contiguous across VPUs. The per-VPU values carried
    # in from finalize collide across VPUs and cannot represent cross-VPU basins;
    # numbering the stitched national network once fixes both. Attribute-only.
    message(sprintf("[stage4][%s] Adding global upstream index ...", domain))
    tryCatch(add_upstream_index(ofile),
      error = function(e) .hf_log_warn(sprintf(
        "[stage4][%s] upstream index skipped: %s", domain, conditionMessage(e))))

    # -- 4m: Invariants ------------------------------------------------------
    message(sprintf("[stage4][%s] Checking invariants ...", domain))

    # attribute-only reads (see 4l): the invariants below are all id/area checks;
    # loading geometry here would re-introduce the ~40GB OOM on large continents.
    # hf_check_merge_invariants uses areasqkm/total_dasqkm attributes only, and its
    # CRS check soft-passes for data.frames (CRS already enforced via -t_srs).
    .con <- DBI::dbConnect(RSQLite::SQLite(), ofile)
    fps <- DBI::dbGetQuery(.con, "SELECT flowpath_id, flowpath_toid, divide_id, areasqkm, total_dasqkm, lengthkm, vpuid, hydroseq, mainstem_id FROM flowpaths")
    nex <- DBI::dbGetQuery(.con, "SELECT nexus_id, nexus_toid FROM nexus")
    div <- DBI::dbGetQuery(.con, "SELECT divide_id, areasqkm FROM divides")
    DBI::dbDisconnect(.con)

    checks <- list()

    # -- ID uniqueness --------------------------------------------------------
    n_dup_fp <- sum(duplicated(fps$flowpath_id))
    checks$no_dup_flowpath_id <- list(kind = "check", ok = n_dup_fp == 0L, n = n_dup_fp)
    if (n_dup_fp > 0L)
      .hf_log_warn(sprintf("[stage4][%s] %d duplicate flowpath_id(s)", domain, n_dup_fp))

    n_dup_div <- sum(duplicated(stats::na.omit(div$divide_id)))
    checks$no_dup_divide_id <- list(kind = "check", ok = n_dup_div == 0L, n = n_dup_div)
    if (n_dup_div > 0L)
      .hf_log_warn(sprintf("[stage4][%s] %d duplicate divide_id(s)", domain, n_dup_div))

    n_dup_nex <- sum(duplicated(nex$nexus_id))
    checks$no_dup_nexus_id <- list(kind = "check", ok = n_dup_nex == 0L, n = n_dup_nex)
    if (n_dup_nex > 0L)
      .hf_log_warn(sprintf("[stage4][%s] %d duplicate nexus_id(s)", domain, n_dup_nex))

    # -- Prefix integrity -----------------------------------------------------
    n_bad_fp_pfx  <- sum(!startsWith(fps$flowpath_id, "fp-"), na.rm = TRUE)
    checks$fp_prefix <- list(kind = "check", ok = n_bad_fp_pfx == 0L, n = n_bad_fp_pfx)
    if (n_bad_fp_pfx > 0L)
      .hf_log_warn(sprintf("[stage4][%s] %d flowpath_id(s) missing 'fp-' prefix",
        domain, n_bad_fp_pfx))

    n_bad_div_pfx <- sum(
      !is.na(div$divide_id) &
        !startsWith(div$divide_id, "cat-") &
        !startsWith(div$divide_id, "lnd-"),
      na.rm = TRUE)
    checks$div_prefix <- list(kind = "check", ok = n_bad_div_pfx == 0L, n = n_bad_div_pfx)
    if (n_bad_div_pfx > 0L)
      .hf_log_warn(sprintf("[stage4][%s] %d divide_id(s) missing 'cat-'/'lnd-' prefix",
        domain, n_bad_div_pfx))

    # Terminal nexuses use "tnx-" prefix; both "nex-" and "tnx-" are valid.
    n_bad_nex_pfx <- sum(!startsWith(nex$nexus_id, "nex-") &
      !startsWith(nex$nexus_id, "tnx-"), na.rm = TRUE)
    checks$nex_prefix <- list(kind = "check", ok = n_bad_nex_pfx == 0L, n = n_bad_nex_pfx)
    if (n_bad_nex_pfx > 0L)
      .hf_log_warn(sprintf("[stage4][%s] %d nexus_id(s) missing 'nex-'/'tnx-' prefix",
        domain, n_bad_nex_pfx))

    # -- Topology resolution --------------------------------------------------
    orphan_toid <- setdiff(fps$flowpath_toid, c(nex$nexus_id, NA_character_))
    checks$toids_resolve <- list(kind = "check", ok = length(orphan_toid) == 0L, n = length(orphan_toid))
    if (length(orphan_toid) > 0L)
      .hf_log_warn(sprintf("[stage4][%s] %d unresolved flowpath_toid(s)", domain,
        length(orphan_toid)))

    orphan_nex_toid <- setdiff(nex$nexus_toid, c(fps$flowpath_id, NA_character_))
    checks$nexus_toids_resolve <- list(kind = "check", ok = length(orphan_nex_toid) == 0L,
      n  = length(orphan_nex_toid))
    if (length(orphan_nex_toid) > 0L)
      .hf_log_warn(sprintf("[stage4][%s] %d unresolved nexus_toid(s)", domain,
        length(orphan_nex_toid)))

    # -- Foreign key: fps$divide_id -> divides ---------------------------------
    fp_div_ids  <- stats::na.omit(fps$divide_id)
    orphan_fp_div <- sum(!fp_div_ids %in% div$divide_id)
    checks$fp_divide_fk <- list(kind = "check", ok = orphan_fp_div == 0L, n = orphan_fp_div)
    if (orphan_fp_div > 0L)
      .hf_log_warn(sprintf("[stage4][%s] %d flowpath divide_id(s) not in divides layer",
        domain, orphan_fp_div))

    # -- No two flowpaths sharing one divide -----------------------------------
    fp_div_nonNA <- fps$divide_id[!is.na(fps$divide_id)]
    n_shared_div <- sum(duplicated(fp_div_nonNA))
    checks$no_shared_divide <- list(kind = "check", ok = n_shared_div == 0L, n = n_shared_div)
    if (n_shared_div > 0L)
      .hf_log_warn(sprintf("[stage4][%s] %d divide_id(s) claimed by >1 flowpath",
        domain, n_shared_div))

    # -- DAG ------------------------------------------------------------------
    is_dag <- tryCatch(
      withCallingHandlers(
        hfutils::hf_network_is_dag(fps, id_col = "flowpath_id", toid_col = "flowpath_toid"),
        warning = function(w) {
          if (grepl("NAs introduced by coercion", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
      ),
      error = function(e) FALSE
    )
    checks$is_dag <- list(kind = "check", ok = is_dag)
    if (!is_dag)
      .hf_log_warn(sprintf("[stage4][%s] network is NOT a DAG", domain))

    # Conservation + completeness invariants (hfutils): catch SILENT data loss
    # the structural checks above can't see: a dropped VPU/feature, lost
    # catchment area, or zeroed drainage area. Input expectations are summed over
    # this domain's per-VPU source gpkgs. strict=FALSE here (warn + record); the
    # domain's own strict gate below decides whether to abort.
    expected  <- tryCatch(.merge_input_expected(gpkg_map, vpus),
      error = function(e) NULL)
    merge_inv <- tryCatch(
      hfutils::hf_check_merge_invariants(
        list(flowpaths = fps, divides = div, nexus = nex),
        expected = expected, stage = sprintf("merge:%s", domain), strict = FALSE,
        gpkg = ofile),
      error = function(e) {
        .hf_log_warn(conditionMessage(e))
        list(ok = TRUE, checks = list())
      })
    checks <- c(checks, merge_inv$checks)
    # Fold in the ngen suite that ran over the re-IDed network earlier in this
    # domain: flowline id uniqueness, toid resolution, the flowline DAG and the
    # POI contract. It was logged but never recorded, so invariants.json described
    # only part of what the build actually asserted.
    if (!is.null(ngen_net_inv) && !is.null(ngen_net_inv$checks))
      checks <- c(checks, ngen_net_inv$checks[
        setdiff(names(ngen_net_inv$checks), names(checks))])

    all_ok <- all(vapply(checks, function(x) isTRUE(x$ok), logical(1L)))

    # Total network length (miles) for flowpaths (aggregated reaches) and
    # flowlines (sub-flowline segments), summed per domain network.
    .km_to_mi <- 0.621371192
    fp_len_mi <- if ("lengthkm" %in% names(fps))
      round(sum(fps$lengthkm, na.rm = TRUE) * .km_to_mi, 1) else NA_real_
    fl_len_mi <- NA_real_
    if (hfutils::layer_exists(ofile, "flowlines")) {
      # attribute-only length sum (flowlines is the largest layer: geometry load
      # would OOM); SUM in SQL avoids pulling the column into R at all.
      .con <- DBI::dbConnect(RSQLite::SQLite(), ofile)
      fll_len <- tryCatch(
        DBI::dbGetQuery(.con, "SELECT SUM(lengthkm) AS s FROM flowlines")$s,
        error = function(e) NA_real_)
      DBI::dbDisconnect(.con)
      if (length(fll_len) && !is.na(fll_len))
        fl_len_mi <- round(fll_len * .km_to_mi, 1)
    }

    inv <- list(
      timestamp       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      domain          = domain,
      ok              = all_ok,
      checks          = checks,
      n_flowpaths     = nrow(fps),
      n_divides       = nrow(div),
      n_nexus         = nrow(nex),
      flowpath_length_mi = fp_len_mi,
      flowline_length_mi = fl_len_mi,
      n_cross_vpu_fix = if (domain == "superconus") nrow(cross_vpu_fix) else 0L,
      vpus            = vpus
    )

    inv_path <- file.path(out_dir,
      sprintf("invariants_%s.json", domain))
    writeLines(jsonlite::toJSON(inv, auto_unbox = TRUE, pretty = TRUE), inv_path)

    message(sprintf(
      "[stage4][%s] Done. %d flowpaths (%s mi) | %d divides | %d nexus | %s flowline mi. %s",
      domain, nrow(fps), format(fp_len_mi, big.mark = ","), nrow(div), nrow(nex),
      format(fl_len_mi, big.mark = ","),
      if (all_ok) "All invariants OK." else "INVARIANT FAILURES -- see warnings."
    ))

    if (isTRUE(strict) && !all_ok) {
      stop(sprintf("[stage4][%s] strict mode: invariants failed; refusing to write non-canonical output", domain))
    }

    # Stamp the merged domain product with the same license/attribution/provenance
    # block the per-VPU build applies: stage 4 is a build output too, so it is
    # licensed by construction. Credit NGA when the domain contains a TDX-Hydro
    # VPU (Alaska); otherwise Lynker-Spatial-only. vpuid is left unset because the
    # product spans many VPUs.
    # Stage 4 reads no rasters, so the elevation credits are carried forward from
    # the VPU products being merged rather than derived from paths. Without this
    # a CC BY 4.0 grid credited on every constituent VPU is absent from the
    # domain product, which is the artifact that actually ships.
    # A domain credits what any of its VPUs owes. Hawaii and PRVI are the only
    # public-domain-throughout domains, so a domain is Lynker-only when EVERY
    # constituent VPU is one of those.
    dom_prefix <- sub("[A-Za-z].*", "", as.character(vpus))
    dom_attr   <- if (all(dom_prefix %in% .HF_PUBLIC_DOMAIN_VPUS))
      .HF_ATTR_LS else .hf_attribution(vpus[!dom_prefix %in% .HF_PUBLIC_DOMAIN_VPUS][1L])
    inherited  <- .hf_grid_attribution_inherited(unlist(gpkg_map[as.character(vpus)]))
    inherited  <- inherited[!vapply(inherited, grepl, logical(1L),
      x = dom_attr, fixed = TRUE)]
    if (length(inherited))
      dom_attr <- paste(c(dom_attr, inherited), collapse = " ")
    .hf_stamp_outputs(ofile, attribution = dom_attr,
      license = .hf_data_license(NA_character_))

    # Type integer-valued id columns as INTEGER + index the toid join keys. Pure
    # SQL (no geometry read), so it is safe on the largest domains; prefixed
    # string ids (fp-/cat-/nex-) stay TEXT via the value guard.
    enforce_integer_ids(ofile)

    all_results[[domain]] <- list(outfile = ofile, inv_path = inv_path,
      invariants = inv)
    global_fl_base <- global_fl_base + dom_fl_span   # advance global flowline offset
  }

  invisible(list(
    outfiles      = lapply(all_results, `[[`, "outfile"),
    fp_map        = fp_map,
    poi_map       = poi_map,
    cross_vpu_fix = cross_vpu_fix,
    invariants    = lapply(all_results, `[[`, "invariants")
  ))
}

# Sum input expectations across a domain's per-VPU source gpkgs, for the
# conservation invariants (hf_check_merge_invariants). Lightweight: count rows +
# sum divide area via DBI (no geometry read). Returns NULL if no inputs resolve.
.merge_input_expected <- function(gpkg_map, vpus) {
  paths <- unlist(gpkg_map[as.character(vpus)], use.names = FALSE)
  paths <- paths[!is.na(paths) & nzchar(paths) & file.exists(paths)]
  if (!length(paths)) return(NULL)
  nfp <- 0L
  ndv <- 0L
  area <- 0
  for (p in paths) {
    con  <- DBI::dbConnect(RSQLite::SQLite(), p)
    tabs <- DBI::dbListTables(con)
    if ("flowpaths" %in% tabs)
      nfp <- nfp + DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM flowpaths")$n
    if ("divides" %in% tabs) {
      ndv <- ndv + DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM divides")$n
      a <- tryCatch(DBI::dbGetQuery(con, "SELECT SUM(areasqkm) a FROM divides")$a,
        error = function(e) NA_real_)
      if (!is.na(a)) area <- area + a
    }
    DBI::dbDisconnect(con)
  }
  list(vpus = as.character(vpus), n_flowpaths = nfp, n_divides = ndv, area_sqkm = area)
}

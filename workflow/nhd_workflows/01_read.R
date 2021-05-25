source('workflow/utils.R')

#############################################################################
out_path <- "workflow/nhd_workflows/cache/ngen_01a-0.gpkg"
RPU = "01a"
dir = "/Users/mjohnson/nhd_rds"
#############################################################################
# Catchments only carry around geometry and ID
in_fl <- readRDS(file.path(dir, 'nhdplus_flowline_update.rds')) %>%
  filter(RPUID == RPU) %>%
  select(COMID, Hydroseq, geometry = Shape, order = StreamOrde,
         LevelPathI, toCOMID, TerminalFl) %>%
  mutate(lengthkm = as.numeric(set_units(st_length(.), "km"))) %>%
  setNames(tolower(names(.)))

in_cat <- readRDS(file.path(dir, 'nhdplus_catchment.rds')) %>%
  select(comid = FEATUREID) %>%
  filter(comid %in% in_fl$comid) %>%
  rename(geometry = Shape) %>%
  mutate(areasqkm = as.numeric(set_units(st_area(.), "km2"))) %>%
  setNames(tolower(names(.)))

# Output ------------------------------------------------------------------

make_plot(in_fl, in_cat, "Raw NHDPlus HF") %>%
  ggsave(filename = "workflow/nhd_workflows/cache/img/00-raw-gf-v20.png",
         units = "in", height = 4, width = 8)

write_sf(in_fl, out_path,  "flowline")
write_sf(in_cat, out_path, "catchment")

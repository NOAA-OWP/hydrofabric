source('workflow/utils.R')

#############################################################################
out_path <- "workflow/nhd_workflows/cache/ngen_01a-0.gpkg"
RPU      <- "01a"
dir      <- "/Users/mjohnson/nhd_rds"
##### <- ####################################################################
# Catchments only carry around geometry and ID
in_fl <- readRDS(file.path(dir, 'nhdplus_flowline_update.rds')) %>%
  filter(RPUID == RPU) %>% 
  select(ToNode, 
         FromNode, 
         COMID, 
         #StartFlag, TerminalFl, 
         Hydroseq,  levelpath = LevelPathI, 
         order    = StreamOrde,
         geometry = Shape) %>%
  st_transform(5070) %>% 
  mutate(lengthkm = as.numeric(units::set_units(st_length(.), "km"))) %>%
  janitor::clean_names()

in_fl$from_comid = in_fl$comid[match(in_fl$from_node, in_fl$to_node)]
in_fl$from_comid = ifelse(is.na(in_fl$from_comid), -1, in_fl$from_comid )

in_fl$to_comid   = in_fl$comid[match(in_fl$to_node,  in_fl$from_node)]
in_fl$to_comid   = ifelse(is.na(in_fl$to_comid), -1,  in_fl$to_comid )

in_cat <- readRDS(file.path(dir, 'nhdplus_catchment.rds')) %>%
  select(comid = FEATUREID) %>%
  filter(comid %in% in_fl$comid) %>%
  rename(geometry = Shape) %>%
  st_transform(5070) %>% 
  mutate(areasqkm = as.numeric(units::set_units(sf::st_area(.), "km2"))) %>%
  janitor::clean_names()

 # Output ------------------------------------------------------------------

make_plot(in_fl, in_cat, "Raw NHDPlus HF") %>%
  ggsave(filename = "workflow/nhd_workflows/cache/img/00-raw-nhd.png",
         units = "in", height = 4, width = 8)

write_sf(in_fl, out_path,  "flowline")
write_sf(in_cat, out_path, "catchment")

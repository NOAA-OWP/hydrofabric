#############################################################################
library('sfarrow')
out_path <- "workflow/graph_wf/basedata/"
RPU      <- "01a"
dir      <- "/Users/mjohnson/nhd_rds"
#############################################################################
# Catchments only carry around geometry and ID
in_fl <- readRDS(file.path(dir, 'nhdplus_flowline_update.rds')) %>%
  filter(RPUID == RPU) %>% 
  select(ID = COMID, 
         hydroseq = Hydroseq,  
         levelpath = LevelPathI, 
         geometry = Shape) %>%
  mutate(member_comid = ID) %>% 
  st_transform(5070) %>% 
  flowpaths_to_linestrings()

outfile = paste0(out_path, RPU,'-base-nhd-flowpaths.parquet')
sfarrow::st_write_parquet(obj=in_fl, dsn=outfile)


system.time({
  arrow::open_dataset("data/test.parquet") %>% 
    filter(levelpath == 2238605) %>% 
    read_sf_dataset()
})

in_cat <- readRDS(file.path(dir, 'nhdplus_catchment.rds')) %>%
  select(ID = FEATUREID, geometry = Shape) %>%
  inner_join(st_drop_geometry(in_fl)) %>% 
  st_transform(5070)

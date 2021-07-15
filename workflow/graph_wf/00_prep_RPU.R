#############################################################################
library('sfarrow')
out_path <- "workflow/graph_wf/basedata/"
RPU      <- "01a"
dir      <- "/Users/mjohnson/nhd_rds"

master_set_fl = c("ID", "levelpath", 'hydroseq', "comid", 'lengthkm')
master_set_fl = c("ID", 'areasqkm')
#############################################################################
# Catchments only carry around geometry and ID
in_fl <- readRDS(file.path(dir, 'nhdplus_flowline_update.rds')) %>%
  filter(RPUID == RPU) %>% 
  select(ID = COMID, 
         hydroseq = Hydroseq,  
         levelpath = LevelPathI, 
         geometry = Shape) %>%
  mutate(comid = ID) %>% 
  st_transform(5070) %>% 
  flowpaths_to_linestrings()

outfile = paste0(out_path, RPU,'-base-nhd-flowpaths.parquet')
sfarrow::st_write_parquet(obj=in_fl, dsn=outfile)

in_cat <- readRDS(file.path(dir, 'nhdplus_catchment.rds')) %>%
  select(ID = FEATUREID, geometry = Shape) %>%
  inner_join(st_drop_geometry(in_fl)) %>% 
  st_transform(5070)

outfile_cat = paste0(out_path, RPU,'-base-nhd-cats.parquet')
sfarrow::st_write_parquet(obj=in_cat, dsn=outfile_cat)


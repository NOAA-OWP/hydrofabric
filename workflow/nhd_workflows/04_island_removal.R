source('workflow/utils.R')
#############################################################################
path       <- "workflow/nhd_workflows/cache/ngen_01a-2.gpkg"
out_path   <- "workflow/nhd_workflows/cache/ngen_01a-3.gpkg"
#############################################################################

# Preprocessing .... ------------------------------------------------------

cat = read_sf(path, layer = "catchments")
fl  = read_sf(path, layer = "flowpaths") %>%
  select(comid, tocomid, levelpathi, hydroseq, order, member_COMID)

# Interior Ring Removal and Fill ------------------------------------------
# example = 15192 (outside), 15193 (interior)

inters = st_intersects(cat)

interior_rings = data.frame(ID = cat$comid,
                            n = lengths(inters)) %>%
  filter(n == 2)

int_list = inters[which(cat$comid %in% interior_rings$ID)]

all = cat$comid[unlist(int_list)]

exteriors = all[!all %in% interior_rings$ID]

message("Dissolving ", nrow(interior_rings), " islands...")

outs = list()
for(i in 1:nrow(interior_rings)){
  outs[[i]] = cat[int_list[[i]],] %>%
    mutate(area = st_area(.), id = 1) %>%
    arrange(-area) %>%
    group_by(id) %>%
    summarize(comid = comid[1])
}

new_exteriors = bind_rows(outs) %>%
  select(-id) %>%
  group_by(comid) %>%
  summarise()

new_cat = cat %>%
  filter(!comid %in% all) %>%
  bind_rows(new_exteriors) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6))

new_fl = fl %>%
  filter(comid %in% new_cat$comid) %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3))

#--------------------------------------------------
#############################################################################

make_plot(new_fl, new_cat, "Island Removal HF") %>%
  ggsave(filename = "workflow/nhd_workflows/cache/img/03-ngen-island-remove.png",
         units = "in", height = 4, width = 8)

#if(all(validate_network(new_fl, new_cat))){
unlink(out_path)
write_sf(new_cat, out_path, "catchments")
write_sf(new_fl, out_path, "flowpaths")
#} else {
#  message("Error in: ", basename(rstudioapi::getSourceEditorContext()$path))
#}





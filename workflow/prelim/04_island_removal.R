source('workflow/utils.R')

#############################################################################
path       <- "workflow/cache/ngen_01a-3.gpkg"
out_path   <- "workflow/cache/ngen_01a-4.gpkg"
#############################################################################

# Preprocessing .... ------------------------------------------------------

cat = read_sf(path, layer = "catchments")
fl  = read_sf(path, layer = "flowpaths")

# Interior Ring Removal and Fill ------------------------------------------
# example = 15192 (outside), 15193 (interior)

inters = st_intersects(cat)

interior_rings = data.frame(ID = cat$ID,
                            n = lengths(inters)) %>%
  filter(n == 2)

int_list = inters[which(cat$ID %in% interior_rings$ID)]

all = cat$ID[unlist(int_list)]

exteriors = all[!all %in% interior_rings$ID]

message("Dissolving ", nrow(interior_rings), " islands...")

outs = list()
for(i in 1:nrow(interior_rings)){
  outs[[i]] = cat[int_list[[i]],] %>%
    mutate(area = st_area(.), id = 1) %>%
    arrange(-area) %>%
    group_by(id) %>%
    summarize(ID = ID[1])
}

new_exteriors = bind_rows(outs) %>%
  select(-id) %>%
  group_by(ID) %>%
  summarise()

new_cat = cat %>%
  filter(!ID %in% all) %>%
  bind_rows(new_exteriors) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6))

new_fl = fl %>%
  filter(ID %in% new_cat$ID) %>%
  inner_join(st_drop_geometry(new_cat), by = 'ID')

#--------------------------------------------------
#############################################################################

make_plot(new_fl, new_cat, "Island Removal 2 HF") %>%
  ggsave(filename = "workflow/cache/img/04-ngen-island-remove.png",
         units = "in", height = 4, width = 8)


if(all(validate_network(new_fl, new_cat))){
  unlink(out_path)
  write_sf(new_cat, out_path, "catchments")
  write_sf(new_fl, out_path, "flowpaths")
} else {
  message("Error in: ", basename(rstudioapi::getSourceEditorContext()$path))
}


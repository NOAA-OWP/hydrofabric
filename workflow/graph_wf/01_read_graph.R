source('workflow/utils.R')

#############################################################################
out_path <- "workflow/graph_wf/cache/ngen_01a-0.gpkg"
RPU      <- "01a"
dir      <- "workflow/graph_wf/basedata/"
##### <- ####################################################################

in_fl  <- sfarrow::st_read_parquet(paste0(dir, RPU,'-base-nhd-flowpaths.parquet'))
in_cat <- sfarrow::st_read_parquet(paste0(dir, RPU,'-base-nhd-cats.parquet'))


### NEED TO DROP and MERGE FLOWLINES W/O CATS
# BUILD GRAPH NET
network = as_sfnetwork(in_fl) %>% 
  activate('edges') %>% 
  st_as_sf()

meta = select(st_drop_geometry(network), ID, hydroseq, toLP=levelpath, comid)

id_map = select(st_drop_geometry(network), from, to, ID, levelpath) %>% 
  # FIND FL without CATS
  filter(!ID %in% in_cat$ID) %>% 
  # Identify fromID from fromNODE
  left_join(select(st_drop_geometry(network), fromID  = ID, tmpto = to),   
            by = c("from" = "tmpto")) %>% 
  # Identify toID from toNODE
  left_join(select(st_drop_geometry(network), toID    = ID, tmpfrom = from),
            by = c("to" = "tmpfrom")) %>%  
  # If to/from ID are not in cats, then set to NA
  mutate(toID = ifelse(toID %in% in_cat$ID, toID, NA),
         fromID = ifelse(fromID %in% in_cat$ID, fromID, NA)) %>% 
  # Group By ID to determine unique inflows and outflow
  group_by(ID) %>% 
  mutate(fromCount = sum(!is.na(fromID)),
         toCount   = sum(!is.na(toID))) %>% 
  # Define first pass merge into
  mutate(send_to = ifelse(toCount == 1, toID, fromID)) %>% 
  # drop
  select(-from, -to) %>% 
  #ungroup
  ungroup()

id_map2 = id_map %>%
  group_by(ID) %>% 
  inner_join(meta, by = c("send_to" = "ID")) %>% 
  mutate(n = n()) %>% 
  #mutate(member_comid = paste(ID, send_to, sep = ",")) %>% 
  ungroup() 

purgeable = filter(id_map, !ID %in% id_map2$ID)

# Here we start to find cases where multiple IDs are being 
# sent to the same toID
tmp = id_map2 %>% 
  group_by(send_to) %>% 
  mutate(n = n())   %>% 
  ungroup() %>% 
  filter(n > 1) 

# TRY reversing to_send to selection
rev = tmp %>% 
  mutate(send_to = ifelse(fromCount == 1, fromID, toID)) %>% 
  group_by(send_to) %>% 
  mutate(n = n()) 

#IF FROMs are NA and going to the same toID, drop the shorter one
stubborn = filter(rev, n > 1) %>% 
  filter(sum(is.na(fromID)) != n) %>% 
  ungroup()

rm_stubborn = filter(rev, n > 1) %>% 
  filter(sum(is.na(fromID)) == n) %>% 
  ungroup()

rev2 = filter(rev, n == 1) 

to_merge = filter(id_map2, !ID %in% tmp$ID) %>% 
  bind_rows(rev2, stubborn) %>% 
  mutate(geom  = NA,
         comid = NA)

for(i in 1:nrow(to_merge)){
  to_merge$geom[i] = build_flow_line(filter(network, ID == to_merge$ID[i])$geometry,
                                     filter(network, ID == to_merge$send_to[i])$geometry)
  to_merge$comid[i] = paste(to_merge$ID[i], to_merge$send_to[i], sep = ",")
}

to_merge2 = to_merge %>% 
  group_by(send_to) %>% 
  mutate(n = n()) %>% 
  st_as_sf(crs = 5070)

singles = filter(to_merge2, n == 1) 

combs   = filter(to_merge2, n  > 1) %>% 
  group_by(send_to) %>% 
  summarize(comid     = paste(comid, collapse = ","),
            levelpath = levelpath[1],
            hydroseq  = hydroseq[1]) %>% 
  flowpaths_to_linestrings()

combs$comid = sapply(strsplit(combs$comid, ","), 
                             function(x) paste(unique(x), collapse = ","))

cleaned = bind_rows(singles, combs) %>% 
  mutate(ID = NULL) %>% 
  rename(ID = send_to, geometry = geom)

og_net = network %>% 
  filter(!ID %in% to_merge$ID) %>% 
  filter(!ID %in% to_merge$send_to) %>% 
  filter(!ID %in% purgeable$ID) %>% 
  filter(!ID %in% rm_stubborn$ID) %>% 
  mutate(comid = as.character(comid))

ooo = bind_rows(og_net, cleaned)

xx = prep_ngen(fl = ooo, cat = in_cat, ID_col = "ID")

 # Output ------------------------------------------------------------------
# make_plot(xx$fl, xx$cat, "Raw NHDPlus HF") %>%
#   ggsave(filename = "workflow/nhd_workflows/cache/img/00-raw-nhd.png",
#          units = "in", height = 4, width = 8)

unlink(out_path)
write_sf(xx$nex, out_path,  "nex")
write_sf(xx$fl, out_path, "flowline")
write_sf(xx$cat, out_path, "catchment")

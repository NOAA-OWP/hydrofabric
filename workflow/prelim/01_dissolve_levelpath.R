source('workflow/utils.R')

#############################################################################
path         <- "workflow/cache/ngen_01a-0.gpkg"
out_path     <- "workflow/cache/ngen_01a-1.gpkg"
ideal_size   <-  10   # in km2
#############################################################################

# Preprocessing .... ------------------------------------------------------
cat = read_sf(path, "catchments")

fl  = read_sf(path, "flowpaths") %>%
  select(-areasqkm) %>%
  left_join(st_drop_geometry(cat), by = "ID") %>%
  #needed to handle the flowlines with missing cats
  mutate(areasqkm =  ifelse(is.na(areasqkm), 0, areasqkm))

singles = fl %>%
  group_by(LevelPathID) %>%
  mutate(n = n()) %>%
  filter(n == 1) %>%
  mutate(member_ID = ID) %>%
  select(
    c(
      "LevelPathID",
      "Hydroseq",
      "member_ID",
      "ID" ,
      "stream_order",
      "member_COMID",
      "toID",
      "member_ID"
    )
  ) %>%
  ungroup()

multis = fl %>%
  st_drop_geometry() %>%
  group_by(LevelPathID) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  arrange(-Hydroseq) %>%
  mutate(ind = cs_group(.data$areasqkm, ideal_size)) %>%
  arrange(Hydroseq) %>%
  ungroup()

mapping = multis %>%
  group_by(LevelPathID, ind) %>%
  summarize(
    Hydroseq = first(Hydroseq),
    member_ID    = paste(ID, collapse = ","),
    ID = first(ID),
    stream_order = first(stream_order),
    member_COMID = paste(member_COMID, collapse = ",")
  ) %>%
  mutate(toID = lead(ID)) %>%
  ungroup()

tab1 =  multis %>%
  left_join(select(fl, ID, geom), by = "ID") %>%
  data.table()

fl_xx = tab1[, .(geom = st_union(geom)), by = list(LevelPathID, ind)]

out_fl = st_as_sf(fl_xx) %>%
  left_join(mapping, by = c("LevelPathID", 'ind')) %>%
  bind_rows(mutate(singles, member_ID = as.character(member_ID))) %>%
  flowpaths_to_linestrings()

out_fl$toID = ifelse(is.na(out_fl$toID),
                     fl$toID[match(out_fl$ID, fl$ID)],
                     out_fl$toID)

to_fix = filter(out_fl,!toID %in% ID) %>%
  filter(!is.na(toID))

good = filter(out_fl,!ID %in% to_fix$ID)

ofl_drop = st_drop_geometry(out_fl)

ind = lapply(1:nrow(to_fix), function(x) {
  grep(
    paste0("(?<![^,])", to_fix$toID[x], "(?![^,])"),
    ofl_drop$member_ID,
    value = FALSE,
    perl = TRUE
  )
})

to_fix$toID = out_fl$ID[unlist(ind)]

rr = bind_rows(good, to_fix) %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3))

######################

single_cat = singles %>%
  st_drop_geometry() %>%
  inner_join(cat, by = "ID") %>%
  st_as_sf()

multi_cat = multis %>%
  inner_join(select(cat, ID, geom), by = "ID") %>%
  data.table()

xx = multi_cat[, .(geom = st_union(geom)), by = list(LevelPathID, ind)]

out_cat = st_as_sf(xx) %>%
  left_join(mapping, by = c("LevelPathID", 'ind')) %>%
  select(ID) %>%
  bind_rows(select(single_cat, ID)) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6))

#--------------------------------------------------
#############################################################################

make_plot(rr, out_cat, "LevelPath Dissolve HF") %>%
  ggsave(
    filename = "workflow/cache/img/02-ngen-levelpath-dissolve.png",
    units = "in",
    height = 4,
    width = 8
  )

message("Dropped: ", nrow(fl) - nrow(rr), " features")

#if (all(validate_network(rr, out_cat)) {
unlink(out_path)
write_sf(out_cat, out_path, "catchments")
write_sf(rr, out_path, "flowpaths")
# } else {
#   message("Error in: ",
#           basename(rstudioapi::getSourceEditorContext()$path))
# }

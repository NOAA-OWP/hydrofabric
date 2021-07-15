source('workflow/utils.R')

#############################################################################
path         <- "workflow/nhd_workflows/cache/ngen_01a-1.gpkg"
out_path     <- "workflow/nhd_workflows/cache/ngen_01a-2.gpkg"
ideal_size   <-  10   # in km2
#############################################################################

# Preprocessing .... ------------------------------------------------------
cat = read_sf(path, "catchments")
fl  = read_sf(path, "flowpaths") %>%
  mutate(member_COMID = comid) %>% 
  filter(order > 0) %>%
  left_join(st_drop_geometry(cat), by = "comid") %>%
  mutate(areasqkm = ifelse(is.na(areasqkm), 0, areasqkm))

singles = fl %>%
  group_by(levelpathi) %>%
  mutate(n = n()) %>%
  filter(n == 1)  %>%
  ungroup()

multis = fl %>%
  st_drop_geometry() %>%
  group_by(levelpathi) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  arrange(-hydroseq) %>%
  mutate(ind = cs_group(.data$areasqkm, ideal_size)) %>%
  arrange(hydroseq) %>%
  ungroup()

mapping = multis %>%
  group_by(levelpathi, ind) %>%
  summarize(
    hydroseq = first(hydroseq),
    member_COMID = paste(comid, collapse = ","),
    comid = first(comid),
    order = first(order),
    terminalfl = max(terminalfl),
    n = n()
  ) %>%
  mutate(tocomid = lead(comid)) %>%
  ungroup()


tab1 =  multis %>%
  left_join(select(fl, comid, geom), by = "comid") %>%
  data.table()

fl_xx = tab1[, .(geom = st_union(geom)), by = list(levelpathi, ind)]

out_fl = st_as_sf(fl_xx) %>%
  left_join(mapping, by = c("levelpathi", 'ind')) %>%
  bind_rows(mutate(singles, member_COMID = as.character(member_COMID))) %>%
  flowpaths_to_linestrings()

out_fl$tocomid = ifelse(is.na(out_fl$tocomid),
                     fl$tocomid[match(out_fl$comid, fl$comid)],
                     out_fl$tocomid)

to_fix = filter(out_fl,!tocomid %in% comid) %>%
  filter(!is.na(tocomid))

good = filter(out_fl,!comid %in% to_fix$comid)

ofl_drop = st_drop_geometry(out_fl)

ind = lapply(1:nrow(to_fix), function(x) {
  tmp = grep(
    paste0("(?<![^,])", to_fix$tocomid[x], "(?![^,])"),
    ofl_drop$member_COMID,
    value = FALSE,
    perl = TRUE
  )

  ifelse(length(tmp) == 0, x, tmp)
})


to_fix$tocomid = out_fl$comid[unlist(ind)]

rr = bind_rows(good, to_fix) %>%
  mutate(lengthkm = as.numeric(st_length(.)/1e3))

######################

single_cat = singles %>%
  st_drop_geometry() %>%
  inner_join(cat, by = "comid") %>%
  st_as_sf()

multi_cat = multis %>%
  inner_join(select(cat, comid, geom), by = "comid") %>%
  data.table()

xx = multi_cat[, .(geom = st_union(geom)), by = list(levelpathi, ind)]

out_cat = st_as_sf(xx) %>%
  left_join(mapping, by = c("levelpathi", 'ind')) %>%
  select(comid) %>%
  bind_rows(select(single_cat, comid)) %>%
  mutate(areasqkm = as.numeric(st_area(.)/1e6))

#--------------------------------------------------
#############################################################################

make_plot(rr, out_cat, "LevelPath Dissolve HF") %>%
  ggsave(
    filename = "workflow/nhd_workflows/cache/img/02-ngen-levelpath-dissolve.png",
    units = "in",
    height = 4,
    width = 8
  )

message("Dropped: ", nrow(fl) - nrow(rr), " features")

unlink(out_path)
write_sf(out_cat, out_path, "catchments")
write_sf(rr, out_path, "flowpaths")


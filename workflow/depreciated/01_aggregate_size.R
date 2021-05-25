library(dplyr)
library(sf)

# Notes: odd projections....
#        some duplicated featrues
#        overly complex spatial topologies
#        more flowpaths then catchments

gpkg = "/Users/mjohnson/github/gfv2.0/workspace/cache/01a.gpkg"

cat = read_sf(gpkg, layer = "divides") %>%
  st_set_crs(5070) %>%
  filter(!duplicated(.)) %>%
  rmapshaper::ms_simplify(.9) %>%
  mutate(area = as.numeric(st_area(.)/1e6))

which(duplicated(cat$ID))

fl  = read_sf(gpkg, layer = "reconciled") %>%
  st_set_crs(4326) %>%
  st_transform(5070) %>%
  filter(!duplicated(.)) %>%
  rmapshaper::ms_simplify(.9)

fl  = filter(fl, ID %in% cat$ID)
cat = filter(cat, ID %in% fl$ID)

#############################################################################
max_size   = 15    # in km2
max_length = 20    # in km
#############################################################################

look_up_table <- st_drop_geometry(fl) %>%
  group_by(LevelPathID) %>%
  tally() %>%
  filter(n > 1)

single_lp_fl  <- filter(fl2, !LevelPathID %in% look_up_table$LevelPathID) %>%
  as_tibble() %>%
  arrange(ID)

single_lp_cat  <- filter(cat2, !LevelPathID %in% look_up_table$LevelPathID) %>%
  as_tibble() %>%
  select(ID, geometry) %>%
  arrange(ID)

multi_lp_fl <- filter(fl2,  LevelPathID %in% look_up_table$LevelPathID) %>%
  as_tibble() %>%
  arrange(ID)

multi_lp_cat <- filter(cat2,  LevelPathID %in% look_up_table$LevelPathID) %>%
  as_tibble() %>%
  select(ID, geometry) %>%
  arrange(ID)

lps = unique(multi_lp_fl$LevelPathID)
out  <- list()

system.time({

for (i in 1:length(lps)) {

  this.lp <- filter(multi_lp_fl, LevelPathID == lps[i]) %>%
    arrange(-Hydroseq) %>%
    select(-geometry)

  areas   <- this.lp$area
  lengths <- this.lp$LENGTHKM
  indexes <- c(); count <- 1

  while (length(areas) > 0) {
    index <- max(1, min(sum(cumsum(areas) <= max_size),
                sum(cumsum(lengths)      <= max_length)))
    inds    <- 1:index
    areas   <- areas[-inds]
    lengths <- lengths[-inds]
    indexes <- c(indexes, rep(count, length(inds)))
    count   <- count + 1
  }

  out[[i]] = mutate(this.lp, inds = as.integer(indexes))
}


})

agg_mappings = do.call(rbind, out) %>%
  arrange(ID) %>%
  select(ID2 = ID, inds, LevelPathID)

o = cbind(select(multi_lp_fl, -LevelPathID), agg_mappings) %>%
  st_as_sf()
sum(o$ID == o$ID2) == nrow(multi_lp_fl)

o2 = cbind(multi_lp_cat, agg_mappings) %>%
  st_as_sf()
sum(o2$ID == o2$ID2) == nrow(multi_lp_cat)

system.time({
  tmp2 = group_by(o, LevelPathID, inds) %>%
    summarize(outlet = which.min(Hydroseq),
              ID = min(ID),
              toID = toID[outlet],
              Hydroseq = Hydroseq[outlet],
              member_COMID = paste(member_COMID, collapse = ",")) %>%
    ungroup() %>%
    select(-inds, -outlet)
})

bool = (st_geometry_type(st_geometry(tmp2)) == "MULTILINESTRING")
multis = tmp2[bool, ]
multis$geometry = st_line_merge(multis$geometry)
singles = tmp2[!bool, ]

rr = bind_rows(multis, singles, single_lp_fl) %>%
  mutate(lengthkm = as.numeric(st_length(.))/1e3)

system.time({
  tmp3 = group_by(st_as_sf(o2), LevelPathID, inds) %>%
    summarize(ID = min(ID)) %>%
    ungroup() %>%
    select(-inds)
})

rr2 = bind_rows(tmp3, single_lp_cat) %>%
  mutate(area = as.numeric(st_area(.))/1e6)

#
# make_plot  = function(fl, cat, title, min, ideal, max){
#   library(patchwork)
#   caty=.25
#   fly = .4
#   ggplot(cat, aes(x = area)) +
#     geom_density() +
#     geom_vline(xintercept = min, col = "gray20") +
#     geom_vline(xintercept = ideal, col = "darkred") +
#     geom_vline(xintercept = max, col = "gray20") +
#     xlim(0, max+5) +
#     ylim(0, caty) +
#     labs(title = title,
#          subtitle = paste(nrow(cat), "basins/flowlines")) +
#     theme_bw() +
#   ggplot(fl, aes(x = lengthkm)) +
#     geom_density() +
#     geom_vline(xintercept = .6) +
#     geom_vline(xintercept = 5) +
#     xlim(0, 15) +
#     ylim(0, fly) +
#     theme_bw()
# }
#
#
# make_plot(select(fl, lengthkm = LENGTHKM), cat, "GF2.0", 5,10,15) /
# make_plot(rr, rr2, "Agg pass", 5,10,15)

path = "data/rf01.gpkg"
unlink(path)

st_write(rr2, path, "catch")
st_write(select(rr,  -LENGTHKM), path, "flowpath", append = TRUE)



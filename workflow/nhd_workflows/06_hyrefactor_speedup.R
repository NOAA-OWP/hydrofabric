out_path   <- "workflow/nhd_workflows/cache/ngen_01a-4.gpkg"

outlets = read_sf(out_path, "POIS") %>% 
  st_drop_geometry()

rec = read_sf("workflow/cache/gfv20/01a.gpkg", 'reconciled')
div = read_sf("workflow/cache/gfv20/01a.gpkg", 'divides')

rec_no_geom = st_drop_geometry(rec) 

ind = lapply(1:nrow(outlets), function(x) {
  tmp = grep(
    paste0("(?<![^,])", outlets$comid[x], "(?![^,])"),
    rec_no_geom$member_COMID,
    value = FALSE,
    perl = TRUE
  )
  ifelse(length(tmp) == 0, x, tmp)
})

outlets$ID = rec$ID[unlist(ind)]
outlets$toID = rec$toID[unlist(ind)]

o = select(outlets, toID, ID) %>%
  mutate(type = case_when(!is.na(toID)~ "outlet",
                          is.na(toID) ~ "terminal")) %>% 
  select(-toID)

rec_outlets = filter(rec, is.na(toID)) %>% 
  select(ID) %>%
  st_drop_geometry() %>% 
  mutate(type = "terminal") %>% 
  bind_rows(o) %>% 
  filter(!duplicated(.))


hopeful = hyRefactor::aggregate_catchments(
  flowpath = rec ,
  divide   = div , 
  outlets  = rec_outlets
)

######################### 
flowpath = rec
divide = div
outlets = rec_outlets
zero_order = NULL
coastal_cats = NULL
da_thresh = NA
only_larger = FALSE
######################### 

if (!is.null(zero_order)) {
    if (is.null(coastal_cats)) 
      stop("must supply coastal_cats with zero order")
    if (st_crs(coastal_cats) != st_crs(divide)) 
      st_transform(coastal_cats, st_crs(divide))
    zero_flowpath <- filter(flowpath, member_COMID %in% do.call(c, 
                                                                zero_order))
    flowpath <- filter(flowpath, !ID %in% zero_flowpath$ID)
    coastal <- lapply(zero_order, function(x, coastal_cats, 
                                           divide) {
      zero_cats <- filter(coastal_cats, FEATUREID %in% 
                            x)
      zero_div <- filter(divide, member_COMID %in% as.character(x))
      st_union(c(st_geometry(zero_cats), st_geometry(zero_div)))[[1]]
    }, coastal_cats = coastal_cats, divide = divide)
    coastal <- st_sfc(coastal, crs = st_crs(divide))
    coastal <- st_sf(ID = names(coastal), geom = coastal)
} else {
    coastal <- NULL
}

if (any(!outlets$ID %in% flowpath$ID))  { stop("Outlet IDs must all be in flowpaths.") }
if (any(!is.na(flowpath$toID[which(flowpath$ID %in% outlets[outlets$type ==   "terminal", ]$ID)]))) {
    stop("Terminal paths must have an NA toID")
}

lps <- hyRefactor:::get_lps(flowpath)
outlets <- hyRefactor:::make_outlets_valid(outlets, flowpath, lps, da_thresh = da_thresh, 
                              only_larger = only_larger) %>% 
  distinct()


if (any(remove_head_div <- !flowpath$ID %in% flowpath$toID & !flowpath$ID %in% divide$ID)) {
    remove_fpaths <- filter(flowpath, remove_head_div)
    message(paste("removing", nrow(remove_fpaths), "headwater/diversion flowlines without catchments."))
    flowpath <- filter(flowpath, !ID %in% remove_fpaths$ID)
  }

lps <- hyRefactor:::get_lps(flowpath)
outlets <- mutate(outlets, ID = paste0("cat-", ID))
divide <- mutate(divide, ID = paste0("cat-", ID))
catchment <- flowpath %>% mutate(toID = ifelse(is.na(toID),  -ID, toID))
nexus <- left_join(select(st_set_geometry(catchment, NULL), 
                            toID = ID), select(st_set_geometry(catchment, NULL), fromID = ID, toID), by = "toID") %>% 
  mutate(nexID = paste0("nex-",  toID), fromID = paste0("cat-", fromID), toID = paste0("cat-",  toID)) %>% 
  select(nexID, fromID, toID)

catchment <- mutate(st_set_geometry(catchment, NULL), 
                    fromID = paste0("nex-",  ID), 
                    toID = paste0("nex-", toID), 
                    ID = paste0("cat-", ID)) %>% 
  select(fromID, toID, cat_ID = ID)

cat_graph <- graph_from_data_frame(d = catchment, directed = TRUE)
outlets <- outlets %>% left_join(select(catchment, 
                                        nexID_stem = fromID, 
                                        ID = cat_ID), by = "ID") %>% 
  left_join(select(catchment,  nexID_terminal = toID, ID = cat_ID), by = "ID") %>% 
  mutate(nexID = case_when(type ==  "outlet" ~ nexID_stem, 
                            type == "terminal" ~ nexID_terminal)) %>% 
  select(ID, type, nexID) %>%
  distinct()

  cat_graph_sort_verts <- topo_sort(cat_graph)
  outlet_verts <- cat_graph_sort_verts[names(cat_graph_sort_verts) %in% 
                                         outlets$nexID]
  outlets <- outlets[match(names(outlet_verts), outlets$nexID), 
  ]
  cat_sets <- data.frame(ID = outlets$ID, nexID = outlets$nexID, 
                         set = I(rep(list(list()), nrow(outlets))), geom = I(rep(list(list()), 
                                                                                 nrow(outlets))), stringsAsFactors = FALSE)
  fline_sets <- data.frame(ID = outlets$nexID, set = I(rep(list(list()), 
                                                           nrow(outlets))), geom = I(rep(list(list()), nrow(outlets))), 
                           stringsAsFactors = FALSE)
  verts <- V(cat_graph)
  us_verts <- c()
  for (cat in seq_len(nrow(cat_sets))) {
    if (cat%%10 == 0) 
      message(paste(cat, "of", nrow(cat_sets)))
    outlet <- filter(outlets, ID == cat_sets$ID[cat])
    ut <- bfs(graph = cat_graph, root = cat_sets$nexID[cat], 
              neimode = "in", order = TRUE, unreachable = FALSE, 
              restricted = verts)
    outlet_id <- as.integer(gsub("^cat-", "", outlets$ID[cat]))
    head_id <- lps$head_ID[which(lps$ID == outlet_id)]
    head_id <- head_id[head_id != outlet_id]
    if (length(head_id) == 0) 
      head_id <- outlet_id
    abort_code <- tryCatch({
      um <- find_um(us_verts, cat_graph, cat_id = cat_sets$nexID[cat], 
                    head_id = head_id, outlet_type = filter(outlets, 
                                                            nexID == cat_sets$nexID[cat])$type)
      abort_code <- ""
    }, error = function(e) e)
    if (abort_code != "") {
      if (!is.na(post_mortem_file)) {
        save(list = ls(), file = post_mortem_file)
        stop(paste("Upstream Main error, post mortem file:", 
                   post_mortem_file))
      }
      else {
        stop(paste("Upstream Main error:", abort_code))
      }
    }
    if (length(us_verts) > 0) {
      vert <- us_verts[which(um %in% us_verts)]
      if (length(vert) == 1) {
        us_verts <- us_verts[!us_verts %in% vert]
      }
    }
    um <- as.numeric(gsub("^nex-", "", um))
    if (0 %in% um) 
      um[um == 0] <- as.numeric(gsub("^cat-", "", outlets[outlets$nexID == 
                                                            "nex-0", ]$ID))
    fline_sets$geom[[cat]] <- filter(flowpath, ID %in% um) %>% 
      st_geometry() %>% st_cast("LINESTRING") %>% st_union()
    abort_code <- tryCatch({
      if (length(cat) > 0 && length(st_geometry_type(fline_sets$geom[[cat]])) > 
          0 && st_geometry_type(fline_sets$geom[[cat]]) == 
          "MULTILINESTRING") {
        fline_sets$geom[[cat]] <- st_line_merge(fline_sets$geom[[cat]])
      }
      fline_sets$geom[[cat]] <- fline_sets$geom[[cat]][[1]]
      fline_sets$set[[cat]] <- um
      ut_verts <- ut$order[!is.na(ut$order) & names(ut$order) != 
                             cat_sets$nexID[cat]]
      cat_sets$set[[cat]] <- filter(catchment, fromID %in% 
                                      names(ut_verts))$cat_ID
      remove <- head_of(cat_graph, unlist(incident_edges(cat_graph, 
                                                         ut_verts, "in")))
      verts <- verts[!verts %in% remove]
      if (outlet$type == "outlet") {
        cat_sets$set[[cat]] <- c(cat_sets$set[[cat]], 
                                 outlet$ID)
        verts <- verts[!names(verts) == outlet$nexID]
      }
      geom <- try(st_union(st_geometry(filter(divide, ID %in% 
                                                unlist(cat_sets$set[cat])))), silent = TRUE)
      if (inherits(geom, "try-error")) {
        geom <- try(st_union(st_make_valid(st_geometry(filter(divide, 
                                                              ID %in% unlist(cat_sets$set[cat]))))))
      }
      if (length(geom) > 0) {
        cat_sets$geom[[cat]] <- geom[[1]]
      }
      else {
        cat_sets$geom[[cat]] <- st_multipolygon()
      }
      abort_code <- ""
    }, error = function(e) e)
    if (abort_code != "") {
      if (!is.na(post_mortem_file)) {
        save(list = ls(), file = post_mortem_file)
        stop(paste("error getting geometry type or with union, post mortem file:", 
                   post_mortem_file))
      }
      else {
        stop(paste("error getting geometry type or with union for line merge:", 
                   abort_code))
      }
    }
    cat_sets$set[[cat]] <- as.numeric(gsub("^cat-", "", cat_sets$set[[cat]]))
    us_verts <- c(us_verts, cat_sets$nexID[cat])
  }
  cat_sets$geom <- st_cast(st_sfc(cat_sets$geom, crs = st_crs(divide)), 
                           "MULTIPOLYGON")
  cat_sets <- select(cat_sets, -nexID)
  cat_sets <- st_sf(cat_sets)
  cat_sets[["ID"]] <- as.numeric(gsub("^cat-", "", outlets$ID))
  fline_sets$geom <- st_sfc(fline_sets$geom, crs = st_crs(flowpath))
  fline_sets <- st_sf(fline_sets)
  fline_sets[["ID"]] <- as.numeric(gsub("^cat-", "", outlets$ID))
  sets <- tidyr::unnest(st_set_geometry(fline_sets, NULL), 
                        cols = c(set))
  next_id <- sets %>% left_join(select(st_set_geometry(flowpath, 
                                                       NULL), ID, toID), by = c(set = "ID")) %>% group_by(ID) %>% 
    filter(!toID %in% set) %>% select(ID, toID) %>% ungroup() %>% 
    distinct() %>% left_join(select(sets, set_toID = ID, 
                                    set), by = c(toID = "set")) %>% select(ID, toID = set_toID)
  fline_sets <- left_join(fline_sets, next_id, by = "ID")
  cat_sets <- left_join(cat_sets, next_id, by = "ID")
  return(list(cat_sets = cat_sets, fline_sets = fline_sets, 
              coastal_sets = coastal))
}
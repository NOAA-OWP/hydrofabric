source("runners/config.R")

modifications = read.csv(gs_file) %>% 
  filter(VPUID != toVPUID) %>% 
  rename(from = COMID, to = toCOMID)

meta = assign_global_identifiers(gpkgs = pipeline$uniform, 
                                 outfiles = pipeline$uniform_global,
                                 modifications = modifications)

for(i in 1:nrow(pipeline)){
  try(append_style(pipeline$uniform[i], layer_names = c("flowpaths", "divides", "hydrolocations")), silent = TRUE)
  try(append_style(pipeline$uniform_global[i], layer_names = c("flowpaths", "divides", "hydrolocations")), silent = TRUE)
}


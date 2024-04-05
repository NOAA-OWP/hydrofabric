source("runners/config.R")

ow = TRUE

for(i in 1:length(vpus)){
  
  # You need Science base access for this to run. If you have a login, 
  # set up a secret/sb_tools.R file with the following three lines:
  # 
  # sb_username = "..."
  # sb_password = "..."
  # sbtools::authenticate_sb(sb_username, sb_password)

  #source('runners/secret/sb_tools.R')
  
  get_hydrofabric(VPU = vpus[i], type = "reference", dir =  base_reference, overwrite = ow)
  
  get_hydrofabric(VPU = vpus[i], type = "refactor",  dir = base_refactored, overwrite = ow)
  
  message(vpus[i])
}


if(!file.exists(huc12_cw)){
  meta = get_bucket_df("lynker-spatial", prefix = "tabular-resources", region = "us-west-2") |>
    filter(grepl(basename(huc12_cw), Key))
  
  save_object(object = meta$Key, bucket = meta$Bucket, file = huc12_cw, region = "us-west-2")
}



# Gridmet downloads 

var = c("pet", "pr", "tmmn", "tmmx")
year = c(2000:2020)
outdir = '/Volumes/Transcend/ngen/climate/gridmet'
grid = expand.grid(var, year)
urls = paste0('https://www.northwestknowledge.net/metdata/data/',
              grid$Var1,
              '_',
              grid$Var2,
              '.nc')

out = file.path(outdir, grid$Var1, paste0(grid$Var1,'_', grid$Var2,'.nc'))

fs::dir_create(file.path(outdir, var))

lapply(1:length(urls), function(i){
  if(!file.exists(out[i])){ httr::GET(urls[i], write_disk(out[i]), progress()) }
})


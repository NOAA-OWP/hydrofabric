
sf::sf_use_s2(FALSE)

tmp = tempfile(fileext = ".geojson")
httr::GET("https://earth-info.nga.mil/php/download.php?file=hydrobasins_level2", httr::write_disk(tmp))

xx = sf::read_sf('/Users/mjohnson/Downloads/hydrobasins_level2.geojson')

xx2 = st_make_valid(ms_simplify(xx))

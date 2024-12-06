
library(hexSticker)
library(showtext)

imgurl <- "man/figures/imgfile.png"
font_add_google("Rammetto One", "rammetto")
showtext_auto()

hexSticker::sticker(imgurl, 
                    package="hydrofabric",
                    p_family = "rammetto",
                    #p_fontface = "static",
                    y = 1.45,
                    p_size=15, 
                    p_color = "gold",
                    s_x=1, 
                    s_y=.8, 
                    s_width=.75,
                    h_color = "brown",
                    h_size = 2,
                    h_fill = "dodgerblue",
                    filename="man/figures/logo.png")

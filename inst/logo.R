
library(hexSticker)
library(showtext)

imgurl <- "man/figures/imgfile.png"

hexSticker::sticker(imgurl, 
                    package="hydrofabric",
                    y = 1.5,
                    p_size=17, 
                    p_color = "brown",
                    s_x=1, 
                    s_y=.8, 
                    s_width=.75,
                    h_color = "dodgerblue",
                    h_fill = "gray90",
                    filename="man/figures/logo.png")

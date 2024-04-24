## code to prepare `sticker`
library(hexSticker)
s <- sticker(subplot  = "vignettes/heart-white.png",
             package  ="heRmes",
             h_color  = "#f3b906",
             h_fill   = "#002855",
             p_color  = "white",
             p_size   = 20,
             p_y      = 1.5,
             s_x      = 1,
             s_y      = 0.75,
             s_width  = 0.45,
             filename ="man/figures/hex.png")

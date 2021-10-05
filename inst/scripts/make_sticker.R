library(hexSticker)

imgurl <- "monaLisa_logo_v1_transparent.png"
colr_pname <- "#5e6769"
colr_border <- "#df663e"
colr_fill <- "#fce2d9"

for (frmt in c("png", "pdf")) {
    stk <- sticker(subplot = imgurl, s_x = 1, s_y = 0.9, ## subplot name, location
                   s_width = 0.85, asp = 320/400, ## subplot width, aspect ratio
                   package = "monaLisa", p_x = 1, p_y = 1.55, p_color = colr_pname, ## package name, location, color
                   p_family = "Aller_Lt", p_size = 7.2, ## package font, size
                   h_fill = colr_fill, h_color = colr_border, h_size = 1.8, ## color to fill hexagon, hexagon border color, border width
                   spotlight = TRUE, l_x = 1, l_y = 1.4, ## spotlight position
                   filename = paste0("../www/monaLisa.", frmt),
                   url = "www.bioconductor.org")
    stk <- stk + geom_url("Image credit: http://vectorish.com/lisa-simpson.html",
                          x = 1, y = 0.45, color = "grey", 
                          family = "Aller_Rg", size = 0.75, angle = 0, hjust = 0.8)
    save_sticker(paste0("../www/monaLisa.", frmt), stk, dpi = 300)
}

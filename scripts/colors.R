#' Color palette to be used in the analyses

library(RColorBrewer)
colors_main <- c("#A5D7E5", "#FEC567")
colors_multi <- c("#a5d7e5",
                  "#fec567",
                  "#db7f9b",
                  "#b1babe",
                  "#afdd91",
                  "#c48ec0",
                  "#be8f87",
                  "#b2ecd4",
                  "#ffefba",
                  "#e3c6ea",
                  "#acb791",
                  "#f4b9c1",
                  "#a9c7cd",
                  "#77a1bf",
                  "#e0c79b",
                  "#fb7272",
                  "grey50")
colors_insulin <- c("#badce6", "#54a2b8")
colors_mval <- c("5.71" = brewer.pal(5, "YlGn")[2],
                 "7" = brewer.pal(5, "YlGn")[3],
                 "9.29" = brewer.pal(5, "YlGn")[4],
                 "12.4" = brewer.pal(5, "YlGn")[5]) 


# color_func_blue <- colorRampPalette(colors = c("#ccebf0", "#3083a7"))  # "#3097a7"
color_func_blue <- colorRampPalette(colors = c("#C6DBEF", "#075a84"))  # "#3097a7"
color_func_orange <- colorRampPalette(colors = c("#ffe2b4", "#e58f01"))
color_func_pink <- colorRampPalette(colors = c("#edbbcb", "#ca436c"))
color_func_green <- colorRampPalette(colors = c("#daefcb", "#639343"))  # "#77c544"
color_func_purp <- colorRampPalette(colors = c("#ecd9f1", "#967bce"))

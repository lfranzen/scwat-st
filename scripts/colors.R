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
colors_mval <- c("5.71" = "#a6dbbb",
                 "7" = "#3cb67b",
                 "9.29" = "#359566",
                 "12.4" = "#2e7651")
colors_adi <- c("Adipocyte 1" = "#BFABC3",
                "Adipocyte 2" = "#907198",
                "Adipocyte 3" = "#62376E")

color_low <- "#441153"
color_low2 <- "#62376e"
color_low3 <- "#772A82"
color_high <- "#3CB67B"
color_high2 <- "#1D793D"
colors_ramp <- colorRampPalette(colors = c(color_low, color_high))

#' Update witth 'viridis' colors
color_func_blue <- colorRampPalette(colors = c("#C6DBEF", "#075a84"))
color_func_orange <- colorRampPalette(colors = c("#F3E55C", "#E8602D"))
color_func_pink <- colorRampPalette(colors = c("#bfabc3", "#62376e"))
color_func_green <- colorRampPalette(colors = c("#a6dbbb", "#359566"))
color_func_purp <- colorRampPalette(colors = c("#edbbcb", "#c480a3"))

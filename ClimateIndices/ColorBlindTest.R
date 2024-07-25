#install.packages("PNWColors")
library(PNWColors)
#install.packages("dichromat")
library(dichromat)
#install.packages("recolorize")
library(recolorize)
library(nord)
library(LaCroixColoR)
install.packages("devtools")
devtools::install_github("johannesbjork/LaCroixColoR")

red_green_colors <- lacroix_palette("Berry",n=6, type = "continuous")


# convert to the three dichromacy approximations
trichromacy<- dichromat(red_green_colors, type = "trichromacy")
protan <- dichromat(red_green_colors, type = "protan")
deutan <- dichromat(red_green_colors, type = "deutan")
tritan <- dichromat(red_green_colors, type = "tritan")

# plot for comparison
layout(matrix(1:4, nrow = 4)); par(mar = rep(1, 4))
plotColorPalette(red_green_colors)
red_green_colors
recolorize::plotColorPalette(protan, main = "Protanopia")
recolorize::plotColorPalette(deutan, main = "Deutanopia")
recolorize::plotColorPalette(tritan, main = "Tritanopia")

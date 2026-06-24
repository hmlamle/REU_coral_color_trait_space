# HML test 
# blank file

# Kellen test. can push to repo?

library(colordistance)
library(ggplot2)
library(tidyverse)

setwd("C:\\Users\\kddea\\Desktop\\EPASeagrantCoralImages\\Trial1\\T8\\ExportedT8\\Trimmedandmeasured")
plotPixels('cave-ssid15.jpg', 
           lower = c(0, 0, 0), upper = c(0.19, 0.19, 0.19),
           color.space = 'rgb', ref.white = 'D65')

# image_histogram <- getLabHist('cave-ssid15.jpg',
#                               lower = c(0, 0, 0), upper = c(0.1, 0.1, 0.1),
#                               a.bounds = c(-128, 127), b.bounds = c(-128, 127),
#                               bins = c(4, 4, 4), plotting = TRUE, ref.white = 'D65')
# print(image_histogram)

plotPixels('DB-sid18.jpg', 
           lower = c(0, 0, 0), upper = c(0.19, 0.19, 0.19),
           color.space = 'rgb', ref.white = 'D65')

m = loadImage('emerald-sid02.jpg',
              lower = c(0, 0, 0), upper = c(0.1, 0.1, 0.1))

emerald2_rgbs = m$filtered.rgb.2d   |>  
  as.data.frame() |> 
  mutate(vari = (r - g)/(r + g - b))

ggplot(emerald2_rgbs, aes(r, g)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(emerald2_rgbs, aes(g, b)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(emerald2_rgbs, aes(r, vari)) +
  geom_point() +
  geom_smooth(method = "lm")

# EPA/SEAGRANT Coral Spectral Project
# Script to extract RGB pixel data from coral frag ROIs 
# Hannah-Marie Lamle and Kellen Dean
# 6/24/26



library(colordistance)
library(ggplot2)
library(tidyverse)

# dont need to set working directory because frag ROIs are in folder in project
# setwd("C:\\Users\\kddea\\Desktop\\EPASeagrantCoralImages\\Trial1\\T8\\ExportedT8\\Trimmedandmeasured")

# -------------- Playing with the package functions ------------------------
## --------------- Plotting pixels in color space: -------------------------

plotPixels('cave-ssid15.jpg', 
           lower = c(0, 0, 0), upper = c(0.19, 0.19, 0.19),
           color.space = 'rgb', ref.white = 'D65')
plotPixels('DB-sid18.jpg', 
           lower = c(0, 0, 0), upper = c(0.19, 0.19, 0.19),
           color.space = 'rgb', ref.white = 'D65')

## ---------------------- Color histograms: --------------------------------
cave15_hist <- getLabHist('cave-ssid15.jpg',
                              lower = c(0, 0, 0), upper = c(0.1, 0.1, 0.1),
                              a.bounds = c(-128, 127), b.bounds = c(-128, 127),
                              bins = c(4, 4, 4), plotting = TRUE, ref.white = 'D65')
print(cave15_hist)


# ------------------------ Extract RGB values ----------------------------


# create an object for the loaded image
frag = loadImage('fragment_ROIs/emerald-sid03.jpg',
              lower = c(0, 0, 0), upper = c(0.1, 0.1, 0.1))

# now extract the RGB values per pixel:
frag = frag$filtered.rgb.2d   |>  # using RGB values 
  as.data.frame() |>  # converting this list into a dataframe
  mutate(vari = (r - g)/(r + g - b)) # adding a column to the DF for the VARI index

# how correlated is the data? 

ggplot(frag, aes(r, g)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(frag, aes(g, b)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(frag, aes(r, vari)) +
  geom_point() +
  geom_smooth(method = "lm")

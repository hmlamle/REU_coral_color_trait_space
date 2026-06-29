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


# ------------- Trial for batch extracting RGB pixels (Ferrara 2026) --------
# This code is copied and modified from the GitHub repo for the paper 
# Red, green, blue color indices as proxy for Symbiodiniaceae cell 
# density and chlorophyll content during coral bleaching" 
# https://github.com/ErikFerrara/Coral_Color_Index/blob/main/Scripts/Extract_RGB_values/RGB-values_from_jpg.R


# Title: "Extracting RBG values from images in JPG format"
library(png)
library(here)
library(tidyverse)
library(imager)
library(grid)


# get the list of PNG files in the directory
# NOTE: Other useres will have to change to their local file path in line 74
jpg_files <- list.files(path = "C:/Users/hanna/OneDrive - Florida International University/1.Research/2.Bleaching_spectral_traits/3.analysis/REU_coral_color_trait_space/fragment_ROIs", # Insert the file path
                        recursive = TRUE,
                        pattern = ".jpg", 
                        full.names = T)

# Function to extract the RGB values
RGB_values <- function(file_path) {
  img <- load.image (file_path)
  
  #RGB channel values converted into vectors  
  img_df <- data.frame(
    Red = as.vector(img[, , 1]),
    Green = as.vector(img[, , 2]),
    Blue = as.vector(img[, , 3])
  )
  
  img_df_scaled <- img_df %>%
    mutate(across(c(Red, Green, Blue), ~ round(. * 255, 2))) %>% # Scale the data in a 0-255 system
    filter(Red != 0, Blue != 0, Green != 0) # Remove all pure white pixels since they are from the background. 
  # It might be necessary to change the 255 to 0 if the background is black
  # ^^^ For HML and KD images, needed to change background to black (change values to 0)
  
  summary_stats <- img_df_scaled %>%
    summarise(across(c(Red, Green, Blue), list(mean = ~ round(mean(.), 2), SD = ~ round(sd(.), 2)))) %>% 
    rename(
      R = Red_mean,
      G = Green_mean,
      B = Blue_mean,
      R_SD = Red_SD,
      G_SD = Green_SD,
      B_SD = Blue_SD
    )
  
  # Remove the .png from the file name
  file_name <- gsub("\\.jpg", "", basename(file_path))
  
  return(data.frame(ID = file_name, summary_stats))
  
}

#Apply the function to all files and pull the data in a dataset

#change the name of the forder results_df.short according to your needs

Results.df <- lapply(jpg_files, function(file_path) {
  stats <- RGB_values(file_path)
}) %>% 
  bind_rows()


#write out the file for use elsewhere
write.csv(Results.df, file = "frag_RGBs.csv")



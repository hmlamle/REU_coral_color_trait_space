# EPA/SEAGRANT Coral Spectral Project
# Analysis of coral RGB indices and photosynthetic data
# Hannah-Marie Lamle and Kellen Dean
# 6/29/26



# -------------------- Prepare data  --------------------------
## -------------- Load packages and data ----------------------

library(tidyverse)
library(ggplot2)
library(readxl)
library(purrr)
library(broom)

# load joined PAM and RGB dataframe (we did this by hand due to having 
# multiple frags and because of the naming difference between the 
# fragment_ID column in the two dataframes):

CHROMA <- read_excel("CHROMA.xlsx")

## ---------------- Calculate indices -----------------------

# based on the results of Ferrera et al 2026, they found that
# the best performing indices for coral to predict symbiont 
# density and chlorophyll was the red index, grayscale and 
# red green sum. We're going to use all the indices, but 
# may filter down to the best performing and see if this is 
# consistent with what they reported on for 1. different corals and 
# 2. for a different coral health metric 

# This code is copied and modified from the GitHub repo for the paper 
# Red, green, blue color indices as proxy for Symbiodiniaceae cell 
# density and chlorophyll content during coral bleaching" 
# https://github.com/ErikFerrara/Coral_Color_Index/blob/main/Scripts/Extract_RGB_values/RGB-values_from_jpg.R



CHROMA <- CHROMA %>% 
  mutate(across(c(2:4), as.factor)) %>%    # changed the columns for our df
  mutate(across(c(5:10), as.numeric)) %>%  # added this for the RGB data that came in as characters 
  dplyr::mutate(
    Grayscale = round(0.299 * R + 0.587 * G + 0.114 * B, 2),
    r = R / (R + G + B),
    g = G / (R + G + B),
    b = B / (R + G + B),
    `R+G` = R + G,
    `R+B` = R + B,
    `R-G` = R - G,
    `R-B` = R - B,
    `R/G` = ifelse(G > 0, R / G, NA_real_),
    `G/R` = ifelse(R > 0, G / R, NA_real_),
    `(R-G)/(R+G)` = (R - G) / (R + G),
    `(R-B)/(R+B)` = (R - B) / (R + B),
    `(G-R)/(G+R)` = (G - R) / (G + R),
    `(G-B)/(G+B)` = (G - B) / (G + B),
    `R/(G+B)` = ifelse((G + B) > 0, R / (G + B), NA_real_),
    `(G-R)/(R+G-B)` = (G-R)/(R+G-B)
  )

# transform the data to long format: 

CHROMA_long <- CHROMA %>%
  ungroup() %>% 
  pivot_longer(cols = -c(fragment_id, tank, treatment, site, 'F', M, 'F:', Y, R_SD, G_SD, B_SD),
               names_to = "Index", values_to = "Index_values") %>% 
  mutate(Index = factor(Index, levels = c("R", 
                                          "G", 
                                          "B", 
                                          "r", 
                                          "g", 
                                          "b", 
                                          "Grayscale", 
                                          "R+G", 
                                          "R+B",  
                                          "R-G", 
                                          "R-B", 
                                          "R/G",  
                                          "G/R", 
                                          "(R-G)/(R+G)", 
                                          "(R-B)/(R+B)", 
                                          "(G-R)/(G+R)", 
                                          "(G-B)/(G+B)", 
                                          "R/(G+B)",
                                          "(G-R)/(R+G-B)"
  ) ))

str(CHROMA_long)


# ------------------- Run Correlations: ------------------------------
## ------------------- Linear Model stats: ---------------------------

regressions <- CHROMA_long %>%
  group_by(Index) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Y ~ Index_values, data = .x)),
    coefficients = map(model, tidy),
    model_stats = map(model, glance)
  )


slopes <- regressions %>%
  select(Index, coefficients) %>%
  unnest(coefficients) %>%
  filter(term == "Index_values")


r2 <- regressions %>%
  select(Index, model_stats) %>%
  unnest(model_stats) %>%
  select(Index, r.squared, adj.r.squared)

results <- left_join(slopes, r2, by = "Index") %>%
  arrange(desc(adj.r.squared))


##  Plotting: 

red_channel <- ggplot(CHROMA, aes(R, Y)) +
  geom_point() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of red channel with photosynthetic efficiency",
       y = "Fv/Fm")
red_channel
ggsave("figs/redchannel_PAM.png", red_channel, width = 6, height = 5, units = "in")

grayscale <- ggplot(CHROMA, aes(Grayscale, Y)) +
  geom_point() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of Grayscale with photosynthetic efficiency",
       y = "Fv/Fm")
grayscale
ggsave("figs/grayscale_PAM.png", grayscale, width = 6, height = 5, units = "in")

red_green <- ggplot(CHROMA, aes(R+G, Y)) +
  geom_point() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of Red + Green Index with photosynthetic efficiency",
       y = "Fv/Fm")
red_green
ggsave("figs/redgreen_PAM.png", red_green, width = 6, height = 5, units = "in")

red_norm <- ggplot(CHROMA, aes(r, Y)) +
  geom_point() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of normalized red with photosynthetic efficiency",
       y = "Fv/Fm")
red_norm

VARI <- ggplot(CHROMA, aes((G-R)/(R+G-B), Y)) +
  geom_point() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of VARI Index with photosynthetic efficiency",
       y = "Fv/Fm")
VARI
ggsave("figs/VARI_PAM.png", VARI, width = 6, height = 5, units = "in")

red_minus_green <- ggplot(CHROMA, aes(`R-G`, Y, colour = site)) +
  geom_point() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of Red minus green index with photosynthetic efficiency",
       y = "Fv/Fm") +
  facet_wrap(~site)
red_minus_green
ggsave("figs/VARI_PAM.png", VARI, width = 6, height = 5, units = "in")

red_minus_blue <- ggplot(CHROMA, aes(`R-B`, Y, colour = site)) +
  geom_point() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of Red minus blue index with photosynthetic efficiency",
       y = "Fv/Fm") +
  facet_wrap(~site)
red_minus_blue


CHROMA %>%
  count(site)

## ------------------------ Site-wise comparisons: ---------------------







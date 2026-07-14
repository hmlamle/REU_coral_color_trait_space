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
library(ggpubr)
library(scales)
library(dplyr)
library(multcompView)
library(gt)
library(webshot2)
library(viridis)

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
    `(G-R)/(R+G-B)` = (G-R)/(R+G-B),
    Y = `Y`/1000
  )
CHROMA <- CHROMA %>%    #changes JR to RJ
  mutate(site = case_match(site, "JR's Ledge" ~ "RJ's Ledge", .default = site))

site_order <- c(
  "Caves", "BAR", "Hollywood", "Dania", "Graceland",
  "Reverse Reef", "South Canyon", "Dodge Island",
  "RJ's Ledge", "Emerald Reef"
)

CHROMA <- CHROMA %>%
  mutate(site = factor(site, levels = site_order))

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

correlationdata <- results %>%
  select(-term, -estimate, -std.error, -statistic, -adj.r.squared)

correlationdata <- correlationdata %>%
  mutate(`Index` = as.character(`Index`)) 

correlationdata %>%
  ungroup() %>%   # remove tibble grouping
  gt() %>%
  tab_style(
    style = cell_text(color = "red"),
    locations = cells_body(rows = 3)) %>%
  tab_style(
    style = cell_text(color = "blue"),
    locations = cells_body(rows = 2)) %>%
  tab_style(
    style = cell_text(color = "green"),
    locations = cells_body(rows = 1)) %>%
  cols_label(
    Index = "RGB Index",
    p.value = "P-value",
    r.squared = html("R&sup2;") 
    ) %>%
  fmt_number(
    columns = p.value,
    decimals = 3
  ) %>%
  fmt_number(
    columns = r.squared,
    decimals = 2
  ) %>%
  cols_width(
    Index ~ px(90),
    p.value ~ px(90),
    r.squared ~ px(90)
  ) %>%
  tab_header(
    title = "Correlation of PAM Measurements with RGB Indices"
  ) %>%
  tab_options(
    table.border.top.style = "solid",
    table.border.bottom.style = "solid",
    table.margin.left = px(5),
    table.margin.right = px(5),
    table.font.size = px(14)
  ) %>%
  opt_table_outline() %>%
gtsave("PAM_RGB_Correlation_Table.png", vwidth = 1650, vheight = 4800, path = "C:/Users/kddea/Desktop/EPASeagrantCoralImages/REU_coral_color_trait_space/figs")


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

## ---------------- Combine CHROMA DF with lat and long ----------------------

CHROMA <- read_excel("CHROMA.xlsx")
location <- read_excel("FINAL_SITES.xlsx")
TempandDepth <- read_excel("TEMPANDDEPTH.xlsx")

CL <- merge(CHROMA, location, by = "site")

combined_df <- merge(CL, TempandDepth, by = "fragment_id") #not sure why but this introduces more observations
#----------------------- Plot PAM with Latitude -------------------------------

pamlat <- ggplot(combined_df, aes(longitude, Y, color = site)) +
  geom_point() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of Latitude With Photosynthetic Efficiency",
       y = "Fv/Fm",
       x = "Latitude")
pamlat

pamlatbox <- ggplot(combined_df, aes(longitude, Y, color = site)) +
  geom_boxplot() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of Latitude With Photosynthetic Efficiency",
       y = "Fv/Fm",
       x = "Latitude")
pamlatbox
#---------------------- Plot PAM with initial depth ---------------------------


cleandepth <- combined_df %>% #remove na values
  filter(depth_m != "na") %>%
  filter(depth_m != "NA")

cleandepth$depth_m <- as.numeric(cleandepth$depth_m) #convert to numeric

pamdepth <- ggplot(cleandepth, aes(depth_m, Y, color = site)) +
  geom_point() +
  geom_smooth(aes(group = 1), method = "lm", color = "black") +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +
  labs(title = "Correlation of Initial Depth With Photosynthetic Efficiency",
       y = "Fv/Fm",
       x = "Initial Depth (m)") +
  stat_cor(aes(group = 1, label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), #adds R^2 and p value
           label.x = 2.5, label.y = 300,)
pamdepth

ggsave("figs/DEPTH_PAM.png", pamdepth, width = 10, height = 8, units = "in")

#--------------------- Plot PAM with initial temp ------------------------------

pamtemp <- ggplot(combined_df, aes(temp_C, Y, color = site)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(title = "Correlation of Initial Temperature With Photosynthetic Efficiency",
       y = "Fv/Fm",
       x = "Initial Temperature (C)") +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
           label.x = 24, label.y = 2)
pamtemp


#--------------------PAM ANOVA and Tukey----------------------------------------

  anova_pam <- aov(Y ~ site, data = CHROMA)
  summary(anova_pam)

pamtukey <- TukeyHSD(anova_pam)

pamcld <- multcompLetters4(anova_pam, pamtukey)

pamletter_data <- CHROMA %>%
  group_by(site) %>%
  summarise(
    third_q = quantile(Y, 0.75),
    .groups = 'drop'
  )
pamletters_df <- data.frame(
  site = names (pamcld$site$Letters), 
  Letters = pamcld$site$Letters
                            )

pamletter_data <- left_join(pamletter_data, pamletters_df, by = "site")


pam <- ggplot(CHROMA, aes(x = Y, y = site)) +
  annotate("rect",
           xmin = 0.5, xmax = Inf,
           ymin = -Inf, ymax = Inf,
           fill = "forestgreen", alpha = 0.15) +
  
  annotate("rect",
           xmin = -Inf, xmax = 0.5,
           ymin = -Inf, ymax = Inf,
           fill = "tomato", alpha = 0.15)+
 geom_boxplot(aes(fill = site), alpha = 1) +
  theme_minimal() +
  geom_text(data = pamletter_data, 
            aes(x = third_q + 0.007, y = site, label = Letters), # Adjust y offset as needed
            size = 5, fontface = "bold") + 
  labs(title = "Comparing Photosynthetic Efficiency by Site",
       y = 'site', x = "Fv/Fm") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 25),      # Changes main title size
    axis.title = element_text(size = 22),      # Changes both X and Y axis labels size
    axis.text = element_text(size = 18),     # Changes both X and Y axis tick marks size
  ) +
  scale_fill_viridis_d(option = "turbo") +
  scale_y_discrete(limits = rev)
pam

ggsave("figs/pambysite.png", pam, width = 12, height = 8, units = "in")
#-------------------Red ANOVA and Tukey-----------------------------------------
cleanred <- colordata %>% #remove na values for red
  filter(R != "na") %>%
  filter(R != "NA")

anova_red <- aov(R ~ site, data = cleanred)
summary(anova_red)

TukeyHSD(anova_red)

#------------------Green ANOVA and Tukey----------------------------------------
cleangreen <- colordata %>% #remove na values for green
  filter(G != "na") %>%
  filter(G != "NA")

anova_green <- aov(B ~ site, data = cleangreen)
summary(anova_green)

TukeyHSD(anova_green)

#------------------Blue ANOVA and Tukey-----------------------------------------
cleanblue <- colordata %>% #remove na values for blue
  filter(B != "na") %>%
  filter(B != "NA")

anova_blue <- aov(B ~ site, data = cleanblue)
summary(anova_blue)

TukeyHSD(anova_blue)



#------------------R-G ANOVA and Tukey------------------------------------------
cleanrminusg <- CHROMA %>%
  filter(R-G != "na") %>%
  filter(R-G != "NA")

anova_redminusgreen <- aov(R-G ~ site, data = cleanrminusg)
summary(anova_redminusgreen)

rminusgtukey <- TukeyHSD(anova_redminusgreen)

rminusgcld <- multcompLetters4(anova_redminusgreen, rminusgtukey)

rminusgletter_data <- cleanrminusg %>%
  group_by(site) %>%
  summarise(
   # max_val = max(`R-G`),
    third_q = quantile(`R-G`, 0.75),
    .groups = 'drop'
  )
rminusgletters_df <- data.frame(site = names (rminusgcld$site$Letters),
  Letters = rminusgcld$site$Letters)
#rminusgletters_df$site <- rownames(rminusgletters_df)

rminusgletter_data <- left_join(rminusgletter_data, rminusgletters_df, by = "site")

rminusg <- ggplot(CHROMA, aes(x = `R-G`, y = site)) +
  geom_boxplot(aes(fill = site), alpha = 1) +
  theme_minimal() +
  geom_text(data = rminusgletter_data, 
            aes(x = third_q + 0.8, y = site, label = Letters), # Adjust y offset as needed
            size = 5, fontface = "bold") + 
  labs(title = "Red Minus Green Index by Site",
       y = "Site", x = "R-G") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 25),      # Changes main title size
    axis.title = element_text(size = 22),      # Changes both X and Y axis labels size
    axis.text = element_text(size = 18)      # Changes both X and Y axis tick marks size
  ) +
  scale_fill_viridis_d(option = "turbo") +
  scale_y_discrete(limits = rev)
rminusg

ggsave("figs/rminusgbysite.png", rminusg, width = 12, height = 8, units = "in")

#------------------------Indices PAM Plots-------------------------------------------
red_minus_green <- ggplot(cleanrminusg, aes(`R-G`, Y)) + #red minus green with pam r^2 = 0.11`R-G`
  geom_point() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of Red Minus Green Index With Photosynthetic Efficiency",
       y = "Fv/Fm") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 25, color = "green"),      # Changes main title size
    axis.title = element_text(size = 22),      # Changes both X and Y axis labels size
    axis.text = element_text(size = 18)        # Changes both X and Y axis tick marks size
  ) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), #Adds R^2 and P value
           label.x = 35, label.y = 0.3, size = 6)
red_minus_green


red_minus_green <- red_minus_green + #adds color ranges
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 0.5, ymax = Inf,
           fill = "forestgreen", alpha = 0.15) +

  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = 0.50,
           fill = "tomato", alpha = 0.15)

bands <- data.frame(
  ymin = c(-Inf, 0.5),
  ymax = c(0.50, Inf),
  Health = factor(
    c("Stressed", "Healthy"),
    levels = c("Healthy", "Stressed")
  )
)
red_minus_green <- red_minus_green + #adds legend
  geom_rect(
    data = bands,
    aes(xmin = -Inf, xmax = Inf,
        ymin = ymin, ymax = ymax,
        fill = Health),
    inherit.aes = FALSE,
    alpha = 0.15
  ) +
  scale_fill_manual(
    name = "Photosynthetic Health",
    values = c(
      "Healthy" = "forestgreen",
      "Stressed" = "tomato"
    )
  )
red_minus_green
ggsave("figs/R-G_PAM.png", red_minus_green, width = 12, height = 8, units = "in")


red_minus_blue <- ggplot(cleanrminusg, aes(`R-B`, Y)) + #red minus blue with pam r^2 
  geom_point() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of Red Minus Blue Index With Photosynthetic Efficiency",
       y = "Fv/Fm") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 25, color = "blue"),      # Changes main title size
    axis.title = element_text(size = 22),      # Changes both X and Y axis labels size
    axis.text = element_text(size = 18)       # Changes both X and Y axis tick marks size
  ) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), #Adds R^2 and P value
           label.x = 60, label.y = 0.3, size = 6)
red_minus_blue


red_minus_blue <- red_minus_blue + #adds color ranges
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 0.5, ymax = Inf,
           fill = "forestgreen", alpha = 0.15) +
  
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = 0.50,
           fill = "tomato", alpha = 0.15)

bands <- data.frame(
  ymin = c(-Inf, 0.5),
  ymax = c(0.50, Inf),
  Health = factor(
    c("Stressed", "Healthy"),
    levels = c("Healthy", "Stressed")
  )
)
red_minus_blue <- red_minus_blue + #adds legend
  geom_rect(
    data = bands,
    aes(xmin = -Inf, xmax = Inf,
        ymin = ymin, ymax = ymax,
        fill = Health),
    inherit.aes = FALSE,
    alpha = 0.15
  ) +
  scale_fill_manual(
    name = "Photosynthetic Health",
    values = c(
      "Healthy" = "forestgreen",
      "Stressed" = "tomato"
    )
  )
red_minus_blue
ggsave("figs/R-B_PAM.png", red_minus_blue, width = 12, height = 10, units = "in")

red <- ggplot(cleanrminusg, aes(`R`, Y)) + #red 
  geom_point() +
  geom_smooth(method = lm) +
  labs(title = "Correlation of Red Index With Photosynthetic Efficiency",
       y = "Fv/Fm") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 25, color = "red"),      # Changes main title size
    axis.title = element_text(size = 22),      # Changes both X and Y axis labels size
    axis.text = element_text(size = 18)        # Changes both X and Y axis tick marks size
  ) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), #Adds R^2 and P value
           label.x = 150, label.y = 0.3, size = 6)
red


red <- red + #adds color ranges
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 0.5, ymax = Inf,
           fill = "forestgreen", alpha = 0.15) +
  
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = 0.50,
           fill = "tomato", alpha = 0.15)

bands <- data.frame(
  ymin = c(-Inf, 0.5),
  ymax = c(0.50, Inf),
  Health = factor(
    c("Stressed", "Healthy"),
    levels = c("Healthy", "Stressed")
  )
)
red <- red + #adds legend
  geom_rect(
    data = bands,
    aes(xmin = -Inf, xmax = Inf,
        ymin = ymin, ymax = ymax,
        fill = Health),
    inherit.aes = FALSE,
    alpha = 0.15
  ) +
  scale_fill_manual(
    name = "Photosynthetic Health",
    values = c(
      "Healthy" = "forestgreen",
      "Stressed" = "tomato"
    )
  )
red
ggsave("figs/REDPAM.png", red, width = 12, height = 10, units = "in")

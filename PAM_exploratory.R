 # EPA/SEAGRANT Coral Spectral Project
# Time Point 1 fragment photosynthetic efficiency exploratory analysis
# Hannah-Marie Lamle 
# 6/25/26

library(tidyverse)
library(ggplot2)
library(readxl)

master <- read.csv("final_ssid_assignments_06062026.csv")
pam <- read_excel("PAM_TP1.xlsx") 


# ---------------- Data Distributions: -------------------------------

summary(pam)
str(pam) # tank is numerical but we need it to be a character for plotting purposes
pam$tank <- as.character(pam$tank)

# PAM distributions by tank: 
pam_dist_tank <- ggplot(pam, aes(tank, Y, fill=tank)) +
  geom_violin()

pam_dist_tank
ggsave("figs/PAM_TP1_tank.png", pam_dist_tank, width = 8, height = 5, units = "in")

# PAM distributions by site: 

pam2 <- left_join(pam, master, by = "fragment_ID")
# have to add in some site NA's (oops)
pam2$site[30] <- "JR's Ledge"
pam2$site[31] <- "JR's Ledge"
pam2$site[127] <- "Dodge Island"
pam2$site[253] <- "JR's Ledge"
pam2$site[271] <- "Emerald Reef"
pam2$site[285] <- "South Canyon"

pam_dist_site <- ggplot(pam2, aes(site, Y, fill=site)) +
  geom_violin() + 
  labs(title = "Photosynthetic Efficiency at TP1", 
       y = "Fv/Fm",
       x = "site") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

pam_dist_site
ggsave("figs/PAM_TP1.png", pam_dist_site, width = 8, height = 5, units = "in")



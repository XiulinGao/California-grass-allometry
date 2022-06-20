#grass-allometry.R

## fit allometry model for 12 California grasses

library(tidyverse)
library(ggplot2)
library(poweRlaw)

mass <- read.csv("../data/biomass.csv", stringsAsFactors = FALSE)
size <- read.csv("../data/plant-size.csv", stringsAsFactors = FALSE)

mass <- mass %>%  mutate(id = paste(spcode, rep, sep="-"))

mass_wide <- mass %>% dplyr::select(-extra.stem) %>% 
  pivot_wider(names_from = organ, values_from = biomass)
mass_wide <- filter(mass_wide, id!="AVBA-30D-2")
#mass_wide <- filter(mass_wide, id!="MURI2-45D-3")

#extrastem <- filter(mass, !is.na(extra.stem)) 
#extrastem <- extrastem %>% select(-c(organ, biomass))

#mass_wide1 <- mass_wide %>% left_join(extrastem, by = c("spcode", "rep","id"))

mass_wide <- mass_wide %>% mutate(across(c("stem","root","leaf","extra_stem","repro"), as.character)) %>% 
  mutate(across(c("stem","root","leaf","extra_stem","repro"), as.numeric))


#mass_wide$root <- as.numeric(as.character(mass_wide1$root))
#mass_wide$leaf <- as.numeric(as.character(mass_wide1$leaf))
#mass_wide$stem <- as.numeric(as.character(mass_wide1$stem))
#mass_wide$repro <- as.numeric(as.character(mass_wide1$repro))


mass_wide <- mass_wide %>% group_by(spcode, rep,id) %>% 
  mutate(stem = sum(stem, extra_stem, na.rm=TRUE),
         agb = sum(stem, leaf, repro, na.rm=TRUE)) %>% ungroup() %>% 
  dplyr::select(-extra_stem)

size <- size %>% mutate(id = paste(spcode, rep, sep="-"))
size <- size %>% dplyr::select(-date) %>% group_by(spcode, rep, id) %>% 
  mutate(dia_base = mean(dia_base1, dia_base2),
         dia_can = mean(dia_can1, dia_can2)) %>% ungroup()
g_2_kg = 0.001
cm2_2_m2 = 0.0001
cm_2_m = 0.01
all <- mass_wide %>% left_join(size, by = c("spcode", "rep","id"))
all <- all %>% mutate(carea = pi*(dia_can*dia_can/4),
                      carea2 = (dia_can1*dia_can2*pi)/400)
all <- all %>% mutate(height = height*cm_2_m,
                      stem = stem*g_2_kg,
                      stem_c = stem*2,
                      root = root*g_2_kg,
                      leaf = leaf*g_2_kg,
                      repro = repro*g_2_kg,
                      agb = agb*g_2_kg,
                      carea = carea*cm2_2_m2)


all_log <- all %>% mutate(across(c("stem","stem_c","root","leaf","repro","agb",
                           "dia_base","dia_can","length","height","carea","carea2"),
                           ~log(.x)))
brdi <- all_log %>% filter(spcode =="BRDI2")
brdi_raw <- all %>% filter(spcode =="BRDI2")

sepu <- all_log %>% filter(spcode=="SEPU8")
vumy <- all_log %>% filter(spcode=="VUMY")
vumy_raw <- all %>% filter(spcode=="VUMY") %>% 
  mutate(l2fr = root/leaf,
         l2fr_mean = mean(l2fr))

ggplot(all_log, aes(leaf, root)) + geom_point() +
  facet_wrap(.~spcode, ncol = 4) + 
  labs(x = "Basal diameter",
       y = "Canopy area")
ggsave("../ba-agb-sp-unlog.jpeg", width = 14, height=12, unit = "cm",
       dpi=300)                        

ggplot(all_log, aes(dia_base, root, color = spcode)) + geom_point() +
  labs(x = "Diameter",
       y = "Root")
ggsave("../ba-agb-all-unlog.jpeg", width = 10, height=8, unit = "cm",
       dpi=300)     


###allometry models
leafm_mod <- lm(leaf ~ dia_base, data = vumy)
summary(leafm_mod) #p1: 0.0001764859; p2:1.64345
brdi <- all_log %>% filter(spcode=="BRDI2")

lf_mod <- lm(leaf ~ dia_base, data = brdi)
summary(lf_mod)
brdi <- brdi %>% mutate(
  lf_predic = dia_base*1.7288-8.5782)
brdi <- brdi %>% mutate(err = (leaf - lf_predic)^2,
                                err_sum = sum(err,na.rm=TRUE),
                                mse = err_sum/23)

h_mod <- lm(height ~ dia_base, data = brdi)
summary(h_mod) 
brdi <- brdi %>% mutate(
  hgt_predic = dia_base*0.49816-2.08603 )
brdi <- brdi %>% mutate(err_hgt = (height - hgt_predic)^2,
                                errhgt_sum = sum(err_hgt,na.rm=TRUE),
                                mse_hgt = errhgt_sum/23)


lfAll_mod <- lm(leaf ~ dia_base, data=all_log)#for height it is very species specific
summary(lfAll_mod)

agbw_mod  <- lm(stem ~ height, data = brdi)
summary(agbw_mod)
brdi <- brdi %>% mutate(
  agb_predic = 3.7512*height - 1.0498)
brdi <- brdi %>% mutate(err_agb = (stem - agb_predic)^2,
                                erragb_sum = sum(err_agb,na.rm=TRUE),
                                mse_agb = erragb_sum/23)


agbAll_mod <- lm(stem ~ dia_base, data=all_log)
summary(agbAll_mod)

ca_mod <- lm(carea ~ dia_base, data = brdi)
summary(ca_mod)

brdi <- brdi %>% mutate(
  ca_predic = dia_base*1.3496 - 4.3136)
brdi <- brdi %>% mutate(err_ca = (carea - ca_predic)^2,
                        errca_sum = sum(err_ca,na.rm=TRUE),
                        mse_ca = errca_sum/23)

ca2_mod <- lm(carea2 ~ dia_base, data=vumy)
summary(ca2_mod)

brdi_raw <- brdi_raw %>% mutate(l2fr = root/leaf,
                                l2dr_mean = mean(l2fr, na.rm=TRUE))
brdi_raw <- brdi_raw %>% mutate(ca_p1 = dia_base^1.7288)

ca_mod2 <- lm(carea ~ ca_p1, data=brdi_raw)
summary(ca_mod2)
all <- all %>% mutate(l2fr = root/leaf,
                      l2fr_mean = mean(l2fr,na.rm=TRUE))

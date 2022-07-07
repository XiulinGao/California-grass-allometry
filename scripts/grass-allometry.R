#grass-allometry.R

## fit allometry model for 12 California grasses

library(tidyverse)
library(ggplot2)
library(poweRlaw)

mass <- read.csv("../data/biomass.csv", stringsAsFactors = FALSE)
size <- read.csv("../data/plant-size.csv", stringsAsFactors = FALSE)
sp <- read.csv("../data/species.csv",stringsAsFactors = FALSE)
  
  
mass <- mass %>%  mutate(id = paste(spcode, rep, sep="-"))

mass_wide <- mass %>% 
  pivot_wider(names_from = organ, values_from = biomass)


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
                      carea = carea*cm2_2_m2,
                      cradius = dia_can/2)
all <- all %>% left_join(sp, by="spcode")


all_log <- all %>% mutate(across(c("stem","stem_c","root","leaf","repro","agb",
                           "dia_base","dia_can","length","height","carea","carea2",
                           "cradius"),
                           ~log(.x)))


## some overall plots

ggplot(all_log, aes(leaf, repro)) + geom_point() +
  facet_wrap(.~spcode, ncol = 4) + 
  labs(x = "Log leaf biomass",
       y = "Log reproductive biomass")

ggplot(all_log, aes(leaf, repro,colour=spcode)) + geom_point() +
  labs(x = "Log leaf biomass",
       y = "Log reproductive biomass")


ggsave("../results/lm-repro-all-log.jpeg", width = 14, height=12, unit = "cm",
       dpi=300)                        

ggplot(all_log, aes(dia_base, height, color = species)) + geom_point() +
  labs(x = "Log basal diameter",
       y = "Log height")

ggplot(all_log, aes(dia_base, height)) + geom_point()+
  facet_wrap(.~spcode,ncol=4)+
  labs(x = "Log basal diameter",
       y = "Log height")

ggsave("../results/ba-hgt-all-log.jpeg", width = 12, height=8, unit = "cm",
       dpi=300) 


ggplot(all_log, aes(leaf, repro, color = photo)) + geom_point() +
  labs(x = "Log leaf biomass",
       y = "Log reproductive biomass")

ggsave("../results/lm-repro-photo-log.jpeg", width = 12, height=8, unit = "cm",
       dpi=300) 

### subset data to annuals

anuals <- all_log %>% filter(growth=="annual")

## define allometry funtion for each organ

allom_lf <- function(df){
  lm(leaf ~ dia_base, data=df)
}

allom_stem <- function(df){
  lm(stem ~ leaf, data=df)
}

allom_rot <- function(df){
  lm(root ~ leaf, data=df)
}

allom_ca <- function(df){
  lm(carea ~ dia_base, data=df)
}

allom_hgt <- function(df){
  lm(height ~ dia_base, data=df)
}

## currently, species-specific models are fit, what about a overall fit?

mod_df <- anuals %>% group_by(spcode) %>% nest()
mod_df <- mod_df %>% mutate(lf_mod = map(data, allom_lf),
                            hgt_mod = map(data,allom_hgt),
                            sm_mod = map(data,allom_stem),
                            rot_mod = map(data,allom_rot),
                            ca_mod = map(data,allom_ca)) 

mod_df <- mod_df %>% mutate(lf = map(lf_mod,broom::glance),
                            hgt = map(hgt_mod,broom::glance),
                            sm = map(sm_mod, broom::glance),
                            rot = map(rot_mod, broom::glance),
                            ca = map(ca_mod,broom::glance)) %>% 
  unnest(c(lf,hgt,sm,rot,ca),
         names_sep=".")
  
mod_df <- mod_df %>% mutate(lf_tidy = map(lf_mod, broom::tidy),
                            hgt_tidy = map(hgt_mod, broom::tidy),
                            sm_tidy = map(sm_mod, broom::tidy),
                            rot_tidy = map(rot_mod, broom::tidy),
                            ca_tidy = map(ca_mod, broom::tidy))

## flat model parameter tables
 mod_params <- mod_df %>% select(spcode, data, lf_tidy,hgt_tidy,sm_tidy,rot_tidy,ca_tidy) %>% 
   unnest(c(lf_tidy,hgt_tidy,sm_tidy,rot_tidy,ca_tidy),names_sep=".")
 
## let's first not correct for intercepts, just use the raw params to predict
## biomass and then do a comparison to observation
 
 test <- mod_params %>% unnest(data, keep_empty = TRUE)
 intcpt <- test %>% filter(lf_tidy.term=="(Intercept)")
 slope <- test %>% filter(lf_tidy.term=="dia_base")

slope <- slope %>% select(spcode,species,photo,growth,id,stem,root,leaf,repro,agb,dia_base,
                          carea,height,lf_tidy.estimate,hgt_tidy.estimate,
                          sm_tidy.estimate,rot_tidy.estimate,ca_tidy.estimate) %>% 
  rename(lf_slope = lf_tidy.estimate,
         hgt_slope = hgt_tidy.estimate,
         sm_slope = sm_tidy.estimate,
         rot_slope = rot_tidy.estimate,
         ca_slope=ca_tidy.estimate)

intcpt <- intcpt %>% select(spcode,species,photo,growth,id,lf_tidy.estimate,hgt_tidy.estimate,
                            sm_tidy.estimate,rot_tidy.estimate,ca_tidy.estimate) %>% 
  rename(lf_incpt = lf_tidy.estimate,
         hgt_incpt = hgt_tidy.estimate,
         sm_incpt = sm_tidy.estimate,
         rot_incpt = rot_tidy.estimate,
         ca_incpt = ca_tidy.estimate)

obs_predic <- slope %>% left_join(intcpt, by=c("spcode","species","photo",
                                               "growth","id")) %>% ungroup()

## correction to intercepts can be made here, we'll skip that for now

obs_predic <- obs_predic %>% mutate(across(c("lf_incpt","hgt_incpt",
                                             "sm_incpt","rot_incpt",
                                             "ca_incpt","leaf","stem",
                                             "root","carea","dia_base",
                                             "height"),exp))

## all in 2-parameter power law
obs_predic <- obs_predic %>% group_by(spcode) %>% 
  mutate(lf_pred = lf_incpt*(dia_base^lf_slope),
         hgt_pred = hgt_incpt*(dia_base^hgt_slope),
         sm_pred = sm_incpt*(leaf^sm_slope),
         rot_pred = rot_incpt*(leaf^rot_slope),
         ca_pred = ca_incpt*(dia_base^ca_slope))

## convert biomass to g
kg_2_g = 1000
obs_predic <- obs_predic %>% mutate(leaf = leaf*kg_2_g,
                                    stem = stem*kg_2_g,
                                    root = root*kg_2_g,
                                    lf_pred = lf_pred*kg_2_g,
                                    sm_pred = sm_pred*kg_2_g,
                                    rot_pred = rot_pred*kg_2_g)

## do a overall plot to compare observation and prediction

## leaf
ggplot(obs_predic,aes(dia_base,leaf)) + geom_point(color="black")+
  geom_point(aes(dia_base,lf_pred),color="red") +
  facet_wrap(.~species, ncol=2) +
  labs(x="Basal diameter(cm)",
       y="Leaf biomass (kg)")
ggsave("../results/model-comp-leaf-sp.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)
## height
ggplot(obs_predic,aes(dia_base,height)) + geom_point(color="black")+
  geom_line(aes(dia_base,hgt_pred),color="red") +
  facet_wrap(.~species, ncol=2) +
  labs(x="Basal diameter(cm)",
       y="Plant height (m)")
ggsave("../results/model-comp-hgt-sp.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

## stem

ggplot(obs_predic,aes(leaf,stem)) + geom_point(color="black")+
  geom_line(aes(leaf,sm_pred),color="red") +
  facet_wrap(.~species, ncol=2) +
  labs(x="Leaf biomass (kg)",
       y="Stem biomass (kg)")
ggsave("../results/model-comp-stem-sp.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

#root

ggplot(obs_predic,aes(leaf,root)) + geom_point(color="black")+
  geom_line(aes(leaf,rot_pred),color="red") +
  facet_wrap(.~species, ncol=2) +
  labs(x="Leaf biomass (kg)",
       y="Root biomass (kg)")
ggsave("../results/model-comp-root-sp.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

## canopy area

ggplot(obs_predic,aes(dia_base,carea)) + geom_point(color="black")+
  geom_line(aes(dia_base,ca_pred),color="red") +
  facet_wrap(.~species, ncol=2) +
  labs(x="Basal daimeter (cm)",
       y="Canopy area (m2)")
ggsave("../results/model-comp-carea-sp.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

## what about an over-all fit?

lf_all <- lm(leaf ~ dia_base, data=anuals) #1.8023,0.0002477724
sm_all <- lm(stem ~ leaf, data=anuals) #1.17965, 2.846699
rot_all <- lm(root ~ leaf, data=anuals) #0.94227, 0.7273264
hgt_all <- lm(height ~ dia_base, data=anuals) #0.43231, 0.2118197
ca_all <- lm(carea ~ dia_base, data=anuals) #1.05456, 0.03446055

all_anls <- all %>% filter(growth=="annual") %>% select(-c(dia_base1,
                                                           dia_base2,
                                                           dia_can1,
                                                           dia_can2))
all_anls <- all_anls %>% mutate(lf_pred = 0.0002477724*(dia_base^1.8023),
                                sm_pred = 2.846699*(leaf^1.17965),
                                rot_pred=0.7273264*(leaf^0.94227),
                                hgt_pred = 0.2118197*(dia_base^0.43231),
                                ca_pred = 0.03446055*(dia_base^1.05456))
#plots
ggplot(all_anls,aes(dia_base,carea)) + geom_point(color="black")+
  geom_line(aes(dia_base,ca_pred),color="red") +
  labs(x="Basal daimeter (cm)",
       y="Canopy area (m2)")
ggsave("../results/model-comp-carea-all.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

ggplot(all_anls,aes(leaf,root)) + geom_point(color="black")+
  geom_line(aes(leaf,rot_pred),color="red") +
  labs(x="Leaf biomass (kg)",
       y="Root biomass (kg)")

ggsave("../results/model-comp-root-all.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

ggplot(all_anls,aes(leaf,stem)) + geom_point(color="black")+
  geom_line(aes(leaf,sm_pred),color="red") +
  labs(x="Leaf biomass (kg)",
       y="Stem biomass (kg)")
ggsave("../results/model-comp-stem-all.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

ggplot(all_anls,aes(dia_base,leaf)) + geom_point(color="black")+
  geom_line(aes(dia_base,lf_pred),color="red") +
  labs(x="Basal diameter (cm)",
       y="Leaf biomass (kg)")
ggsave("../results/model-comp-leaf-all.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

ggplot(all_anls,aes(dia_base,height)) + geom_point(color="black")+
  geom_line(aes(dia_base,hgt_pred),color="red") +
  labs(x="Basal diameter (cm)",
       y="Plant height (m)")
ggsave("../results/model-comp-hgt-all.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)


brdi <- all_log %>% filter(spcode =="BRDI2")
brdi_raw <- all %>% filter(spcode =="BRDI2")

sepu <- all_log %>% filter(spcode=="SEPU8")
vumy <- all_log %>% filter(spcode=="VUMY")
vumy_raw <- all %>% filter(spcode=="VUMY") %>% 
  mutate(l2fr = root/leaf,
         l2fr_mean = mean(l2fr))
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

#grass-allometry.R

## fit allometry model for 12 California grasses

library(tidyverse)
library(RColorBrewer)
require(wesanderson)
library(lmtest)
#library(ggplot2)
#library(poweRlaw)

source("./ggplot-theme.R")

mass <- read.csv("../data/biomass.csv", stringsAsFactors = FALSE)
size <- read.csv("../data/plant-size.csv", stringsAsFactors = FALSE)
sp <- read.csv("../data/species.csv",stringsAsFactors = FALSE)
  
  
mass <- mass %>%  mutate(id = paste(spcode, rep, sep="-"))

mass_wide <- mass %>% 
  pivot_wider(names_from = organ, values_from = biomass)


mass_wide <- mass_wide %>% group_by(spcode, rep,id) %>% 
  mutate(stem  = sum(stem, extra_stem, na.rm=TRUE),
         agb   = sum(stem, leaf, repro, na.rm=TRUE),
         tmass = sum(agb,root, na.rm=TRUE))%>% ungroup() %>% 
  dplyr::select(-extra_stem)

size <- size %>% mutate(id = paste(spcode, rep, sep="-"))
size <- size %>% dplyr::select(-date) %>% group_by(spcode, rep, id) %>% 
  mutate(dia_base = mean(dia_base1, dia_base2),
         dia_can = mean(dia_can1, dia_can2)) %>% ungroup()
g_2_kg = 0.001
cm2_2_m2 = 0.0001
cm_2_m = 0.01
mass_2_C = 2
kg_2_g = 1000

all <- mass_wide %>% left_join(size, by = c("spcode", "rep","id"))
all <- all %>% mutate(carea   = pi*(dia_can*dia_can/4),
                      carea2  = (dia_can1*dia_can2*pi)/400)
all <- all %>% mutate(height     = height*cm_2_m,
                      stem       = stem*g_2_kg,
                      stem_c     = stem*mass_2_C,
                      root       = root*g_2_kg,
                      root_c     = root*mass_2_C,
                      leaf       = leaf*g_2_kg,
                      leaf_c     = leaf*mass_2_C,
                      repro      = repro*g_2_kg,
                      repro_c    = repro*mass_2_C,
                      agb        = agb*g_2_kg,
                      agb_c      = agb*mass_2_C,
                      tmass      = tmass*g_2_kg,
                      tmass_c    = tmass*mass_2_C,
                      carea      = carea*cm2_2_m2,
                      cradius    = dia_can/2,
                      lf2rt      = root_c/leaf_c,
                      seed_alloc = repro_c/tmass_c)
all <- all %>% left_join(sp, by="spcode")


all_log <- all %>% mutate(across(c("stem","stem_c","root","root_c",
                                   "leaf","leaf_c","repro","repro_c",
                                   "agb","agb_c","dia_base","dia_can",
                                   "length","height","carea","carea2"),
                           ~log(.x)))
## some overall plots

ggplot(all_log, aes(leaf, stem)) + geom_point() +
  facet_wrap(.~spcode, ncol = 4) + 
  labs(x = "Log leaf biomass",
       y = "Log stem biomass") +
  prestheme +
  theme(axis.text.x=element_text(angle=90, vjust=0.5))

ggsave("../results/lf-sm-sp-log.jpeg", width = 16, height=14, unit = "cm",
       dpi=300)  

ggplot(all_log, aes(leaf, stem,colour=species)) + geom_point() +
  labs(x = "Log leaf biomass",
       y = "Log stem biomass")+
  scale_color_brewer(palette="Paired") + 
  prestheme                            +
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 14, face="italic"),
        axis.text.x=element_text(angle=90, vjust=0.5))  

ggsave("../results/lf-sm-all-log.jpeg", width = 14, height=10, unit = "cm",
       dpi=300)                        

ggplot(all_log, aes(dia_base, height, color = species)) + geom_point() +
  labs(x = "Log basal diameter",
       y = "Log height") +
  scale_color_brewer(palette="Paired") + 
  prestheme                            +
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 14, face="italic"),
        axis.text.x=element_text(angle=90, vjust=0.5))

ggsave("../results/ba-hgt-all-log.jpeg", width = 14, height=10, unit = "cm",
       dpi=300)

ggplot(all_log, aes(dia_base, height)) + geom_point()+
  facet_wrap(.~spcode,ncol=4)+
  labs(x = "Log basal diameter",
       y = "Log height") +
  prestheme +
  theme(axis.text.x=element_text(angle=90, vjust=0.5))

ggsave("../results/ba-hgt-sp-log.jpeg", width = 16, height=14, unit = "cm",
       dpi=300) 


ggplot(all_log, aes(dia_base, leaf, color = photo)) + geom_point() +
  labs(x = "Log basal diameter",
       y = "Log leaf biomass") + 
  scale_color_manual(values=schwilkcolors[1:2]) +
  prestheme +
   theme(legend.title=element_blank(),
        legend.position="bottom")
  

ggsave("../results/ba-lf-photo-log.jpeg", width = 10, height=8, unit = "cm",
       dpi=300) 

### subset data to annuals

anuals <- all_log %>% filter(growth=="annual")
c3anual <- anuals %>% filter(photo=="C3")

## determine species effects first

null_lf    <- lm(leaf_c ~ dia_base, data=c3anual)
full_lf     <- lm(leaf_c ~ dia_base*spcode, data=c3anual)
#sp_lf   <- lm(leaf_c ~ dia_base + spcode, data=c3anual)
lrtest(null_lf, full_lf) #basal diameter effect on leafC vary between species

ggplot(c3anual, aes(dia_base, leaf_c, color = short)) + geom_point() +
  labs(x = "Log basal diameter",
       y = "Log leaf carbon") +
  geom_smooth(method="lm", se=FALSE)+
  geom_abline(slope=1.91512,intercept=-7.92442, color="black", size=1.5, linetype="dashed")+
  scale_color_manual(values=schwilkcolors[1:4]) + 
  prestheme                            +
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 14, face="italic"))

ggsave("../results/lfmod-c3annual.jpeg", width = 14, height=10, unit = "cm",
       dpi=300)

### leaf-root
null_rt    <- lm(root_c ~ leaf_c, data=c3anual)
full_rt     <- lm(root_c ~ leaf_c*spcode, data=c3anual)
#sp_rt   <- lm(root_c ~ leaf_c + spcode, data=c3anual)
lrtest(null_rt, full_rt) #basal diameter effect on leafC vary between species
anova(full_rt)
ggplot(c3anual, aes(leaf_c, root_c, color = short)) + geom_point() +
  labs(x = "Log leaf carbon",
       y = "Log root carbon") +
  geom_smooth(method="lm", se=FALSE)+
  geom_abline(slope=0.90775,intercept=-0.52793, color="black", size=1.5, linetype="dashed")+
  scale_color_manual(values=schwilkcolors[1:4]) + 
  prestheme                            +
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 14, face="italic"))

ggsave("../results/rtmod-c3annual.jpeg", width = 14, height=10, unit = "cm",
       dpi=300)


## leaf-stem biomass
null_sm    <- lm(stem_c ~ leaf_c, data=c3anual)
full_sm     <- lm(stem_c ~ leaf_c*spcode, data=c3anual)

lrtest(null_sm, full_sm) #basal diameter effect on leafC vary between species
anova(full_sm)
ggplot(c3anual, aes(leaf_c, stem_c, color = short)) + geom_point() +
  labs(x = "Log leaf carbon",
       y = "Log stem carbon") +
  geom_smooth(method="lm", se=FALSE)+
  geom_abline(slope=1.22508,intercept=0.92215, color="black", size=1.5, linetype="dashed")+
  scale_color_manual(values=schwilkcolors[1:4]) + 
  prestheme                            +
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 14, face="italic"))

ggsave("../results/smmod-c3annual.jpeg", width = 14, height=10, unit = "cm",
       dpi=300)

## basal diameter-canopy area
null_ca    <- lm(carea ~ dia_base, data=c3anual)
full_ca     <- lm(carea ~ dia_base*spcode, data=c3anual)
lrtest(null_ca, full_ca)

ggplot(c3anual, aes(dia_base, carea, color = short)) + geom_point() +
  labs(x = "Log basal diameter",
       y = "Log canopy area") +
  geom_smooth(method="lm", se=FALSE)+
  geom_abline(slope=1.2435,intercept=-3.9323, color="black", size=1.5, linetype="dashed")+
  scale_color_manual(values=schwilkcolors[1:4]) + 
  prestheme                            +
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 14, face="italic"))

ggsave("../results/camod-c3annual.jpeg", width = 14, height=10, unit = "cm",
       dpi=300)

## basal diameter-height
null_hgt    <- lm(height ~ dia_base, data=c3anual)
full_hgt     <- lm(height ~ dia_base*spcode, data=c3anual)
lrtest(null_hgt, full_hgt)

ggplot(c3anual, aes(dia_base, height, color = short)) + geom_point() +
  labs(x = "Log basal diameter",
       y = "Log height") +
  geom_smooth(method="lm", se=FALSE)+
  geom_abline(slope=0.69951,intercept=-2.09103, color="black", size=1.5, linetype="dashed")+
  scale_color_manual(values=schwilkcolors[1:4]) + 
  prestheme                            +
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size = 14, face="italic"))

ggsave("../results/hgtmod-c3annual.jpeg", width = 14, height=10, unit = "cm",
       dpi=300)



## define allometry funtions

allom_lf <- function(df){
  lm(leaf_c ~ dia_base, data=df) 
}

allom_stem <- function(df){
  lm(stem_c ~ leaf_c, data=df)
}

allom_rot <- function(df){
  lm(root_c ~ leaf_c, data=df)
}

allom_ca <- function(df){
  lm(carea ~ dia_base, data=df)
}

allom_hgt <- function(df){
  lm(height ~ dia_base, data=df)
}

 
## species-specific models 

mod_df <- anuals %>% group_by(spcode) %>% nest()
mod_df <- mod_df %>% mutate(lf_mod = map(data, allom_lf),
                            hgt_mod = map(data,allom_hgt),
                            sm_mod = map(data,allom_stem),
                            rot_mod = map(data,allom_rot),
                            ca_mod = map(data,allom_ca),
                            lf_pred = map2(.x=lf_mod, .y=data,~predict(object=.x,newdata=.y)),
                            hgt_pred = map2(.x=hgt_mod, .y=data,~predict(object=.x,newdata=.y)),
                            sm_pred = map2(.x=sm_mod, .y=data,~predict(object=.x,newdata=.y)),
                            rot_pred = map2(.x=rot_mod, .y=data,~predict(object=.x,newdata=.y)),
                            ca_pred = map2(.x=ca_mod, .y=data,~predict(object=.x,newdata=.y))) 

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
 mod_params <- mod_df %>% select(spcode, data, lf_tidy,hgt_tidy,sm_tidy,rot_tidy,ca_tidy,
                                 lf_pred,hgt_pred,sm_pred,rot_pred,ca_pred) %>% 
   unnest(c(lf_tidy,hgt_tidy,sm_tidy,rot_tidy,ca_tidy),names_sep=".")
 
 
 
 params_df <- mod_params %>% unnest(c(data,lf_pred,hgt_pred,sm_pred,rot_pred,ca_pred),keep_empty = TRUE)
 
 intcpt <- params_df %>% filter(lf_tidy.term=="(Intercept)")
 slope <- params_df %>% filter(lf_tidy.term=="dia_base")

 slope <- slope %>% select(spcode,short,photo,growth,id,stem_c,root_c,leaf_c,repro_c,agb_c,
                           sm_pred,rot_pred,lf_pred,hgt_pred,ca_pred,dia_base,
                          carea,height,lf_tidy.estimate,hgt_tidy.estimate,
                          sm_tidy.estimate,rot_tidy.estimate,ca_tidy.estimate) %>% 
  rename(lf_slope = lf_tidy.estimate,
         hgt_slope = hgt_tidy.estimate,
         sm_slope = sm_tidy.estimate,
         rot_slope = rot_tidy.estimate,
         ca_slope=ca_tidy.estimate)

intcpt <- intcpt %>% select(spcode,short,photo,growth,id,lf_tidy.estimate,hgt_tidy.estimate,
                            sm_tidy.estimate,rot_tidy.estimate,ca_tidy.estimate) %>% 
  rename(lf_incpt = lf_tidy.estimate,
         hgt_incpt = hgt_tidy.estimate,
         sm_incpt = sm_tidy.estimate,
         rot_incpt = rot_tidy.estimate,
         ca_incpt = ca_tidy.estimate)

obs_predic <- slope %>% left_join(intcpt, by=c("spcode","short","photo",
                                               "growth","id")) %>% ungroup()

## let's correct for intercept, see ref: 
## https://www.jstor.org/stable/1937343?seq=1&cid=pdf-reference#references_tab_contents

obs_predic <- obs_predic %>% group_by(spcode) %>% 
                             mutate(lf_n = length(!is.na(leaf_c)),
                                    sm_n = length(!is.na(stem_c)),
                                    rt_n = length(!is.na(root_c)),
                                    hgt_n = length(!is.na(height)),
                                    ca_n  = length(!is.na(carea)),
                                    lf_sme = (leaf_c-lf_pred)*(leaf_c-lf_pred),
                                    sm_sme = (stem_c-sm_pred)*(stem_c-sm_pred),
                                    rt_sme = (root_c-rot_pred)*(root_c-rot_pred),
                                    hgt_sme = (height-hgt_pred)*(height-hgt_pred),
                                    ca_sme = (carea-ca_pred)*(carea-ca_pred))
                                    
 obs_predic <- obs_predic %>% mutate(lf_cf = sum(lf_sme, na.rm=TRUE)/(lf_n-2)/2,
                                     sm_cf = sum(sm_sme,na.rm=TRUE)/(sm_n-2)/2,
                                     rt_cf = sum(rt_sme, na.rm=TRUE)/(rt_n-2)/2,
                                     hgt_cf = sum(hgt_sme, na.rm=TRUE)/(hgt_n-2)/2,
                                     ca_cf = sum(ca_sme, na.rm=TRUE)/(ca_n-2)/2) %>% 
                                     ungroup()
 
  
obs_predic <- obs_predic %>% mutate(across(c("lf_incpt","hgt_incpt",
                                             "sm_incpt","rot_incpt",
                                             "ca_incpt","leaf_c","stem_c",
                                             "root_c","carea","dia_base",
                                             "height","lf_cf","sm_cf","rt_cf",
                                             "hgt_cf","ca_cf"),exp))

## update intercept
obs_predic <- obs_predic %>% mutate(lf_incpt = lf_incpt*lf_cf,
                                    hgt_incpt = hgt_incpt*hgt_cf,
                                    sm_incpt = sm_incpt*sm_cf,
                                    rot_incpt = rot_incpt*rt_cf,
                                    ca_incpt = ca_incpt*ca_cf)

## all in 2-parameter power law
obs_predic <- obs_predic %>% group_by(spcode) %>% 
  mutate(lf_pred = lf_incpt*(dia_base^lf_slope),
         hgt_pred = hgt_incpt*(dia_base^hgt_slope),
         sm_pred = sm_incpt*(leaf_c^sm_slope),
         rt_pred = rot_incpt*(leaf_c^rot_slope),
         ca_pred = ca_incpt*(dia_base^ca_slope))

## convert carbon to g

obs_predic <- obs_predic %>% mutate(leaf_c = leaf_c*kg_2_g,
                                    stem_c = stem_c*kg_2_g,
                                    root_c = root_c*kg_2_g,
                                    lf_pred = lf_pred*kg_2_g,
                                    sm_pred = sm_pred*kg_2_g,
                                    rot_pred = rot_pred*kg_2_g)

## do a overall plot to compare observation and prediction

## leaf
ggplot(obs_predic,aes(dia_base,leaf_c)) + geom_point(color="black")+
  geom_line(aes(dia_base,lf_pred),color="red") +
  facet_wrap(.~short, ncol=2) +
  labs(x="Basal diameter(cm)",
       y="Leaf carbon (g)") + prestheme
ggsave("../results/model-comp-leaf-sp.jpeg",width = 16, height=14, unit = "cm",
       dpi=300)
## height
ggplot(obs_predic,aes(dia_base,height)) + geom_point(color="black")+
  geom_line(aes(dia_base,hgt_pred),color="red") +
  facet_wrap(.~short, ncol=2) +
  labs(x="Basal diameter(cm)",
       y="Plant height (m)")+ prestheme
ggsave("../results/model-comp-hgt-sp.jpeg",width = 16, height=14, unit = "cm",
       dpi=300)

## stem

ggplot(obs_predic,aes(leaf_c,stem_c)) + geom_point(color="black")+
  geom_line(aes(leaf_c,sm_pred),color="red") +
  facet_wrap(.~short, ncol=2) +
  labs(x="Leaf carbon (g)",
       y="Stem carbon (g)") + prestheme
ggsave("../results/model-comp-stem-sp.jpeg",width = 16, height=14, unit = "cm",
       dpi=300)

#root

ggplot(obs_predic,aes(leaf_c,root_c)) + geom_point(color="black")+
  geom_line(aes(leaf_c,rot_pred),color="red") +
  facet_wrap(.~short, ncol=2) +
  labs(x="Leaf carbon (g)",
       y="Root carbon (g)") + prestheme
ggsave("../results/model-comp-root-sp.jpeg",width = 16, height=14, unit = "cm",
       dpi=300)

## canopy area

ggplot(obs_predic,aes(dia_base,carea)) + geom_point(color="black")+
  geom_line(aes(dia_base,ca_pred),color="red") +
  facet_wrap(.~short, ncol=2) +
  labs(x="Basal daimeter (cm)",
       y="Canopy area (m2)") + prestheme
ggsave("../results/model-comp-carea-sp.jpeg",width = 16, height=14, unit = "cm",
       dpi=300)

## what about an over-all fit?

lf_all <- lm(leaf_c ~ dia_base, data=anuals) #1.8023,0.0002477724
sm_all <- lm(stem_c ~ leaf_c, data=anuals) #1.17965, 2.846699
rot_all <- lm(root_c ~ leaf_c, data=anuals) #0.94227, 0.7273264
hgt_all <- lm(height ~ dia_base, data=anuals) #0.43231, 0.2118197
ca_all <- lm(carea ~ dia_base, data=anuals) #1.05456, 0.03446055

## model params
lf_slp <- lf_all$coefficients[2] 
lf_cpt <- lf_all$coefficients[1]
sm_slp <- sm_all$coefficients[2]
sm_cpt <- sm_all$coefficients[1]
rot_slp <- rot_all$coefficients[2]
rot_cpt <- rot_all$coefficients[1]
hgt_slp <- hgt_all$coefficients[2]
hgt_cpt <- hgt_all$coefficients[1]
ca_slp <- ca_all$coefficients[2]
ca_cpt <- ca_all$coefficients[1]

anuals$lf_pred <- predict(lf_all, newdata=anuals)
anuals$sm_pred <- predict(sm_all, newdata=anuals)
anuals$rot_pred <- predict(rot_all, newdata=anuals)
anuals$hgt_pred <- predict(hgt_all, newdata=anuals)
anuals$ca_pred <- predict(ca_all, newdata=anuals)

anuals <- anuals %>% mutate(lf_n = length(!is.na(leaf_c)),
                            sm_n = length(!is.na(stem_c)),
                            rot_n = length(!is.na(root_c)),
                            hgt_n = length(!is.na(height)),
                            ca_n = length(!is.na(carea)),
                            lf_sme = (leaf_c-lf_pred)*(leaf_c-lf_pred),
                            sm_sme = (stem_c-sm_pred)*(stem_c-sm_pred),
                            rt_sme = (root_c-rot_pred)*(root_c-rot_pred),
                            hgt_sme = (height-hgt_pred)*(height-hgt_pred),
                            ca_sme = (carea-ca_pred)*(carea-ca_pred))
anuals <- anuals %>% mutate (lf_cf = sum(lf_sme, na.rm=TRUE)/(lf_n-2)/2,
                            sm_cf = sum(sm_sme,na.rm=TRUE)/(sm_n-2)/2,
                            rt_cf = sum(rt_sme, na.rm=TRUE)/(rot_n-2)/2,
                            hgt_cf = sum(hgt_sme, na.rm=TRUE)/(hgt_n-2)/2,
                            ca_cf = sum(ca_sme, na.rm=TRUE)/(ca_n-2)/2) %>% 
                            ungroup()

anuals <- anuals %>% mutate(lf_slope = lf_slp,
                            lf_incpt = lf_cpt,
                            sm_slope = sm_slp,
                            sm_incpt = sm_cpt,
                            rot_slope = rot_slp,
                            rot_incpt = rot_cpt,
                            hgt_slope = hgt_slp,
                            hgt_incpt = hgt_cpt,
                            ca_slope = ca_slp,
                            ca_incpt = ca_cpt
                            )

anuals <- anuals %>% mutate(across(c("lf_incpt","hgt_incpt",
                                             "sm_incpt","rot_incpt",
                                             "ca_incpt","leaf_c","stem_c",
                                             "root_c","carea","dia_base",
                                             "height","lf_cf","sm_cf","rt_cf",
                                             "hgt_cf","ca_cf"),exp))

## update intercept
anuals <- anuals %>% mutate(lf_incpt = lf_incpt*lf_cf,
                                    hgt_incpt = hgt_incpt*hgt_cf,
                                    sm_incpt = sm_incpt*sm_cf,
                                    rot_incpt = rot_incpt*rt_cf,
                                    ca_incpt = ca_incpt*ca_cf)


anuals <- anuals %>% group_by(spcode) %>% 
  mutate(lf_pred = lf_incpt*(dia_base^lf_slope),
         hgt_pred = hgt_incpt*(dia_base^hgt_slope),
         sm_pred = sm_incpt*(leaf_c^sm_slope),
         rot_pred = rot_incpt*(leaf_c^rot_slope),
         ca_pred = ca_incpt*(dia_base^ca_slope))

## convert carbon to g

anuals <- anuals %>% mutate(leaf_c = leaf_c*kg_2_g,
                                    stem_c = stem_c*kg_2_g,
                                    root_c = root_c*kg_2_g,
                                    lf_pred = lf_pred*kg_2_g,
                                    sm_pred = sm_pred*kg_2_g,
                                    rot_pred = rot_pred*kg_2_g)

#plots
ggplot(anuals,aes(dia_base,carea)) + geom_point(color="black")+
  geom_line(aes(dia_base,ca_pred),color="red") +
  labs(x="Basal daimeter (cm)",
       y="Canopy area (m2)") + prestheme
ggsave("../results/model-comp-carea-all.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

ggplot(anuals,aes(leaf_c,root_c)) + geom_point(color="black")+
  geom_line(aes(leaf_c,rot_pred),color="red") +
  labs(x="Leaf carbon (g)",
       y="Root carbon (g)") + prestheme

ggsave("../results/model-comp-root-all.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

ggplot(anuals,aes(leaf_c,stem_c)) + geom_point(color="black")+
  geom_line(aes(leaf_c,sm_pred),color="red") +
  labs(x="Leaf carbon (g)",
       y="Stem carbon (g)") + prestheme
ggsave("../results/model-comp-stem-all.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

ggplot(anuals,aes(dia_base,leaf_c)) + geom_point(color="black")+
  geom_line(aes(dia_base,lf_pred),color="red") +
  labs(x="Basal diameter (cm)",
       y="Leaf carbon (g)") + prestheme
ggsave("../results/model-comp-leaf-all.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

ggplot(anuals,aes(dia_base,height)) + geom_point(color="black")+
  geom_line(aes(dia_base,hgt_pred),color="red") +
  labs(x="Basal diameter (cm)",
       y="Plant height (m)")+prestheme
ggsave("../results/model-comp-hgt-all.jpeg",width = 12, height=10, unit = "cm",
       dpi=300)

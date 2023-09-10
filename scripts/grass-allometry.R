## grass-allometry.R

## fit allometry model for 12 California grasses
## 1. model selection to find the best allometry model for each organ
## 2. test some allometry hypotheses 

setwd("~/California-grass-allometry/scripts")

library(tidyverse)
library(lmtest)
library(lme4)
library(car)
library(emmeans)



g_2_kg = 0.001
cm2_2_m2 = 0.0001
cm_2_m = 0.01
mass_2_C = 2
kg_2_g = 1000
wooden = 0.001
day2sec = 86400


mass       <- read.csv("../data/biomass.csv"    ,stringsAsFactors = FALSE)
size       <- read.csv("../data/plant-size.csv" ,stringsAsFactors = FALSE)
sp         <- read.csv("../data/species.csv"    ,stringsAsFactors = FALSE)
sla_mass   <- read.csv("../data/SLA-biomass.csv",stringsAsFactors = FALSE)
sla_area   <- read.csv("../data/leaf-area.csv"  ,stringsAsFactors = FALSE) 


sla        <- sla_area                                              %>%  
              left_join(sla_mass, by=c("spcode","rep"))             %>% 
              mutate(sla = area/biomass)                            %>% 
              group_by(spcode)                                      %>% 
              summarise(sla = mean(sla,na.rm=TRUE))                 %>% 
              ungroup()                                             %>% 
              mutate(sla = sla*cm2_2_m2)

  
mass <- mass                                                        %>%  
        mutate(id = paste(spcode, rep, sep="-"))

mass_wide <- mass                                                   %>% 
             pivot_wider(names_from = organ, values_from = biomass)


mass_wide <- mass_wide                                              %>% 
             group_by(spcode, rep,id)                               %>% 
             mutate(stem  = sum(stem, extra_stem, na.rm=TRUE)
                   ,agb   = sum(stem, leaf, repro, na.rm=TRUE)
                   ,tmass = sum(agb,root, na.rm=TRUE))              %>% 
             ungroup()                                              %>% 
             dplyr::select(-extra_stem)

size <- size                                        %>% 
        mutate(id = paste(spcode, rep, sep="-"))    %>% 
        group_by(spcode,rep,id)                     %>% 
        mutate(dia_base = mean(dia_base1,dia_base2)
              ,dia_can  = mean(dia_can1,dia_can2  ) 
              ,date     = lubridate::mdy(date)    ) %>% 
        ungroup()


all <- mass_wide                                            %>% 
       left_join(size, by = c("spcode", "rep","id"))        %>% 
       left_join(sla,  by =   "spcode")                     %>% 
       left_join (sp,  by =   "spcode")                     %>% 
       filter(!spcode %in% c("POSE","FECA"))
        
all <- all                                                  %>% 
       mutate(carea      = pi*(dia_can*dia_can/4  )
             ,barea      = pi*(dia_base*dia_base/4)      
             ,height     = height*cm_2_m
             ,stem       = stem*g_2_kg
             ,stem_c     = stem*mass_2_C
             ,root       = root*g_2_kg
             ,root_c     = root*mass_2_C
             ,leaf       = leaf*g_2_kg
             ,leaf_c     = leaf*mass_2_C
             ,repro      = repro*g_2_kg
             ,repro_c    = repro*mass_2_C
             ,agb        = agb*g_2_kg
             ,agb_c      = agb*mass_2_C
             ,tmass      = tmass*g_2_kg
             ,tmass_c    = tmass*mass_2_C
             ,barea      = barea*cm2_2_m2
             ,carea      = carea*cm2_2_m2
             ,cradius    = dia_can/2
             ,lf2rt      = root_c/leaf_c
             ,agb2rt     = root_c/agb_c
             ,seed_alloc = repro_c/agb_c
             ,seed_alloct= repro_c/tmass_c
             ,lf2sm      = leaf_c/stem_c
             ,rt2lf      = root_c/leaf_c
             ,ca2ba      = barea/carea
             ,tlarea     = leaf*sla)


all_log <- all                                                  %>% 
           mutate(across(c("stem","stem_c","root","root_c"    , 
                           "leaf","leaf_c","repro","repro_c"  ,
                           "agb","agb_c","tmass","tmass_c"    ,
                           "dia_base","dia_can","length"      ,
                           "height","barea","carea","lf2sm"   ,
                           "rt2lf","agb2rt","ca2ba","sla"     ,
                           "tlarea") , ~log(.x)))


####### 1. Allometry model selection and species-specific variation #######

#### 1a.height allometry model ####

hgt_null   <- lm(height ~ dia_base, all_log)
hgt_full1  <- lm(height ~ dia_base*spcode,all_log)
hgt_full2  <- lm(height ~ dia_base*(photo+growth),all_log)

AIC(hgt_null,hgt_full1, hgt_full2)
hgt_poshoc  <- emmeans(hgt_full1, pairwise ~ dia_base*spcode)


hgt_allom <- hgt_full1
summary(hgt_allom)
anova(hgt_allom)
anova(hgt_full2)

#### 1b.canopy area allometry ####

#ca_null1  <- lm(carea ~ dia_base*height,all_log)
ca_null2  <- lm(carea ~ dia_base + height, all_log)
ca_null3  <- lm(carea ~ dia_base, all_log)
ca_null4  <- lm(carea ~ height, all_log)
#ca_full1  <- lm(carea ~ dia_base*height*spcode,all_log)
ca_full2  <- lm(carea ~ (dia_base + height)*spcode, all_log)
ca_full3  <- lm(carea ~ dia_base*spcode, all_log)
ca_full4  <- lm(carea ~ height*spcode, all_log)
ca_full5  <- lm(carea ~ (dia_base + height)*(photo + growth),all_log)
ca_full6  <- lm(carea ~ dia_base*(photo + growth),all_log)
ca_full7  <- lm(carea ~ height*(photo + growth),all_log)
AIC(ca_null2,ca_null3,ca_null4,ca_full2,ca_full3,ca_full4,ca_full5,ca_full6,ca_full7)
ca_poshoc <- emmeans(ca_full2, pairwise ~ dia_base:height:spcode)


ca_allom <- ca_full2
summary(ca_allom)
anova(ca_allom)
anova(ca_full5)

####1c. leaf allometry ####
lf_null1  <- lm(leaf ~ carea, all_log)
lf_null2  <- lm(leaf ~ dia_base, all_log)
lf_null3  <- lm(leaf ~ height, all_log)
lf_null4  <- lm(leaf ~ dia_base + height, all_log)
lf_null5  <- lm(leaf ~ dia_base + carea, all_log)
lf_null6  <- lm(leaf ~ carea + height, all_log )
lf_null7  <- lm(leaf ~ carea + dia_base + height,all_log)


#lf_full1  <- lm(leaf ~ dia_base*height*spcode, all_log)
lf_full1  <- lm(leaf ~ carea*spcode, all_log)
lf_full2  <- lm(leaf ~ carea*(photo + growth), all_log)
lf_full3  <- lm(leaf ~ dia_base*spcode,all_log)
lf_full4  <- lm(leaf ~ dia_base*(growth + photo), all_log)
lf_full5  <- lm(leaf ~ height*spcode, all_log)
lf_full6  <- lm(leaf ~ height*(growth + photo),all_log)
lf_full7  <- lm(leaf ~ (dia_base + height)*spcode,all_log)
lf_full8  <- lm(leaf ~ (dia_base + height)*(growth + photo),all_log)
lf_full9  <- lm(leaf ~ (dia_base + carea)*spcode,all_log)
lf_full10  <- lm(leaf ~ (dia_base + carea)*(growth + photo),all_log)
lf_full11  <- lm(leaf ~ (carea + height)*spcode,all_log)
lf_full12  <- lm(leaf ~ (carea + height)*(growth + photo),all_log)
lf_full13  <- lm(leaf ~ (dia_base + height + carea)*spcode,all_log)
lf_full14  <- lm(leaf ~ (dia_base + height + carea)*(growth + photo),all_log)


AIC(lf_null1,lf_null2,lf_null3,lf_null4,lf_full1,lf_full2,lf_full3,lf_full4,lf_full5,lf_full6,
    lf_full7,lf_full8,lf_full9,lf_full10,lf_full11,lf_full12,lf_full13,lf_full14)

lf_allom  <- lf_full7
anova(lf_allom)
summary(lf_allom)
lf_poshoc <- emmeans(lf_allom, pairwise ~ dia_base:height:spcode)
lf_poshoc
plot(lf_poshoc)
anova(lf_full8)
## best model is lf_full13, but it does not increase R2 more than 1 comparing to lf_full7
## including canopy area is not that useful

#### 1d.stem allometry ####
sm_null1 <- lm(stem ~ dia_base, all_log)
sm_null2 <- lm(stem ~ height, all_log)
sm_null3 <- lm(stem ~ carea, all_log)
sm_null4 <- lm(stem ~ dia_base + height, all_log)
sm_null5 <- lm(stem ~ dia_base + carea, all_log)
sm_null6 <- lm(stem ~ height + carea, all_log)
sm_null7 <- lm(stem ~ height + dia_base + carea, all_log)
#sm_null4 <- lm(stem ~ dia_base*height,all_log)
sm_full1 <- lm(stem ~ dia_base*spcode, all_log)
sm_full2 <- lm(stem ~ height*spcode, all_log)
sm_full3 <- lm(stem ~ carea*spcode, all_log)
sm_full4 <- lm(stem ~ (dia_base + height)*spcode, all_log)
sm_full5 <- lm(stem ~ (dia_base + carea)*spcode, all_log)
sm_full6 <- lm(stem ~ (height + carea)*spcode, all_log)
sm_full7 <- lm(stem ~ (height + dia_base + carea)*spcode,all_log)
#sm_full4 <- lm(stem ~ dia_base*height*spcode,all_log)
sm_full8 <- lm(stem ~ dia_base*(photo + growth),all_log)
sm_full9 <- lm(stem ~ height*(photo + growth),all_log)
sm_full10 <- lm(stem ~ carea*(photo + growth),all_log)
sm_full11 <- lm(stem ~ (dia_base + height)*(photo+growth),all_log)
sm_full12 <- lm(stem ~ (dia_base + carea)*(photo+growth),all_log)
sm_full13 <- lm(stem ~ (height+ carea)*(photo+growth),all_log)
sm_full14 <- lm(stem ~ (dia_base + height + carea)*(photo+growth),all_log)

### sm_full7 is the best but did not improve R2 more than 1 comapring to sm_full4

AIC(sm_null1,sm_null2,sm_null3,sm_null4,sm_null5,sm_null6,sm_null7,sm_full1,sm_full2,sm_full3,
    sm_full4,sm_full5,sm_full6,sm_full7,sm_full8,sm_full9,sm_full10,sm_full11,sm_full12,
    sm_full13,sm_full14)

sm_allom <- sm_full4
summary(sm_allom)
anova(sm_allom)
sm_poshoc <- emmeans(sm_allom,pairwise ~ dia_base:height:spcode)
sm_poshoc
plot(sm_poshoc)
anova(sm_full11)

#### 1e. root allometry ####
rt_null1 <- lm(root ~ leaf,all_log)
rt_null2 <- lm(root ~ carea, all_log)
rt_null3 <- lm(root ~ dia_base, all_log)
rt_null4 <- lm(root ~ height, all_log)
rt_null5 <- lm(root ~ dia_base + height, all_log)
rt_null6 <- lm(root ~ dia_base + carea, all_log)
rt_null7 <- lm(root ~ height + carea, all_log)
rt_null8 <- lm(root ~ dia_base + height + carea, all_log)

rt_full1 <- lm(root ~ leaf*spcode,all_log)
rt_full2 <- lm(root ~ carea*spcode,all_log)
rt_full3 <- lm(root ~ dia_base*spcode,all_log)
rt_full4 <- lm(root ~ height*spcode, all_log)
rt_full5 <- lm(root ~ (dia_base+height)*spcode,all_log)
rt_full6 <- lm(root ~ (dia_base + carea)*spcode,all_log)
rt_full7 <- lm(root ~ (height + carea)*spcode,all_log)
rt_full8 <- lm(root ~ (dia_base + height + carea)*spcode,all_log)
rt_full9 <- lm(root ~ leaf*(photo + growth),all_log)
rt_full10 <- lm(root ~ carea*(photo + growth),all_log)
rt_full11 <- lm(root ~ dia_base*(photo + growth),all_log)
rt_full12 <- lm(root ~ height*(growth + photo),all_log)
rt_full13 <- lm(root ~ (dia_base+height)*(growth+photo),all_log)
rt_full14 <- lm(root ~ (dia_base+carea)*(growth+photo),all_log)
rt_full15 <- lm(root ~ (carea+height)*(growth+photo),all_log)
rt_full16 <- lm(root ~ (dia_base+height+carea)*(growth+photo),all_log)

AIC(rt_null1,rt_null2,rt_null3,rt_null4,rt_null5,rt_null6,rt_null7,rt_null8,
    rt_full1,rt_full2,rt_full3,rt_full4,rt_full5,rt_full6,rt_full7,rt_full8,
    rt_full9,rt_full10,rt_full11,rt_full12,rt_full13,rt_full14,rt_full15,
    rt_full16)

rt_allom <- rt_full1
summary(rt_allom)
anova(rt_allom)
rt_poshoc <- emmeans(rt_allom,pairwise ~ leaf:spcode)
rt_poshoc
plot(rt_poshoc)
anova(rt_full9)


#### 1f. AGB allometry ####

agb_null1 <- lm(agb ~ dia_base, all_log)
agb_null2 <- lm(agb ~ carea, all_log)
agb_null3 <- lm(agb ~ height, all_log)
agb_null4 <- lm(agb ~ dia_base + height, all_log)
agb_null5 <- lm(agb ~ dia_base + carea, all_log)
agb_null6 <- lm(agb ~ height + carea, all_log)
agb_null7 <- lm(agb ~ height + dia_base + carea, all_log)

agb_full1 <- lm(agb ~ dia_base*spcode, all_log)
agb_full2 <- lm(agb ~ carea*spcode, all_log)
agb_full3 <- lm(agb ~ height*spcode,all_log)
agb_full4 <- lm(agb ~ (dia_base + height)*spcode, all_log)
agb_full5 <- lm(agb ~ (dia_base + carea)*spcode, all_log)
agb_full6 <- lm(agb ~ (height + carea)*spcode, all_log)
agb_full7 <- lm(agb ~ (height + dia_base + carea)*spcode,all_log)

agb_full8 <- lm(agb ~ (height + dia_base + carea)*(photo + growth),all_log)
agb_full9 <- lm(agb ~ (height + dia_base)*(photo + growth),all_log)

AIC(agb_null1,agb_null2,agb_null3,agb_null4,agb_null5,agb_null6,agb_null7,agb_full1,agb_full2,agb_full3,agb_full4,agb_full5,agb_full6,agb_full7,
    agb_full8,agb_full9)
summary(agb_full7)
summary(agb_full4)



## still, adding canopy area won't improve R2 by more than 1

#### 1g. total biomass allometry ####

tm_null1 <- lm(tmass ~ dia_base, all_log)
tm_null2 <- lm(tmass ~ carea, all_log)
tm_null3 <- lm(tmass ~ height, all_log)
tm_null4 <- lm(tmass ~ dia_base + height, all_log)
tm_null5 <- lm(tmass ~ dia_base + carea, all_log)
tm_null6 <- lm(tmass ~ height + carea, all_log)
tm_null7 <- lm(tmass ~ height + dia_base + carea, all_log)

tm_full1 <- lm(tmass ~ dia_base*spcode, all_log)
tm_full2 <- lm(tmass ~ carea*spcode, all_log)
tm_full3 <- lm(tmass ~ height*spcode,all_log)
tm_full4 <- lm(tmass ~ (dia_base + height)*spcode, all_log)
tm_full5 <- lm(tmass ~ (dia_base + carea)*spcode, all_log)
tm_full6 <- lm(tmass ~ (height + carea)*spcode, all_log)
tm_full7 <- lm(tmass ~ (height + dia_base + carea)*spcode,all_log)

tm_full8 <- lm(tmass ~ (height + dia_base + carea)*(photo + growth),all_log)
tm_full9 <- lm(tmass ~ (height + dia_base)*(photo + growth),all_log)

AIC(tm_null1,tm_null2,tm_null3,tm_null4,tm_null5,tm_null6,tm_null7,tm_full1,tm_full2,tm_full3,tm_full4,tm_full5,tm_full6,tm_full7,
    tm_full8,tm_full9)

summary(tm_full7)
summary(tm_full4)

### adding canopy area does not improve R2 by more than 1

####### 2.Hypothesis testing #######

#### 2a. root allocation is greater in perennial, especially C4 perennials #### 

#mature_df      <- all_log %>% filter(!is.na(repro))

frt_hmod    <- lmer(rt2lf ~ growth + photo + (1|spcode)
                   , data = all_log
                   , REML = TRUE)

summary(frt_hmod)
Anova(frt_hmod,type="2",test.statistic = "F") 

#when testing interaction effect, use this setting for Anova





####### H2: fast growlers invest more in reproduction #######

## growth rate

growth <- function(df){
  lm(tmass ~ lagday, data=df) 
}

grow_df   <- all_log                                                %>% 
             group_by(spcode)                                       %>% 
             mutate(lagsec = as.numeric(difftime(date,min(date)))
                   ,lagday = lagsec/day2sec
                   ,lagday = ifelse(lagday==0,1,lagday-1))          %>% 
             select(-lagsec)                                        %>% 
             ungroup()
  
grow_df   <- grow_df                                     %>% 
             group_by(spcode)                            %>% 
             nest()                                      %>% 
             mutate(grow_mod = map(data,growth))

grow_df   <- grow_df                                      %>% 
             mutate(stat = map(grow_mod,broom::glance))   %>% 
             unnest(stat)

grow_df   <- grow_df                                      %>% 
             mutate(mod_tidy = map(grow_mod,broom::tidy)) %>% 
             unnest(mod_tidy,names_sep="_")               

grow_coef <- grow_df                          %>% 
             filter(mod_tidy_term =="lagday") %>% 
             select(spcode
                   ,mod_tidy_estimate)        %>% 
             rename(grate = mod_tidy_estimate)

repro_df   <- all                                                    %>% 
              filter(!is.na(seed_alloct))                            %>%
              left_join(grow_coef,by="spcode")                       %>%
              group_by(spcode)                                       %>% 
              mutate(repro_sd  = sd(seed_alloc)
                    ,repro_se  = repro_sd/sqrt(length(seed_alloc))
                    ,repro_ave = mean(seed_alloc))                   %>% 
              ungroup()                       

                

repro_mod  <- lmer(seed_alloct ~ grate + (1|spcode)
                   ,data=repro_df
                   ,REML=TRUE)
summary(repro_mod)
Anova(repro_mod,type="2",test.statistic = "F") 
#when no interaction effect is test, use type 2

# reproduction allocation (as %above-ground biomass)
# is positively related to growth rate of above-ground biomass


######## H3: leaf/stem ratio is influenced by      #######
####### height (-) and specific leaf area (-)      #######  
#zscore <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm = TRUE) 

h3df         <- all_log                                    %>% 
                select(spcode,height,carea,sla,lf2sm)     
                
lsratio_mod  <- lmer(lf2sm ~ sla*carea + (1|spcode)
                    , data = h3df
                    , REML = TRUE)

summary(lsratio_mod)
Anova(lsratio_mod,type="3",test.statistic = "F")



####### Fit species-specific allometry with intercept correction #######

## define allometry funtions

allom_lf <- function(df){
  lm(leaf ~ dia_base + height, data=df) 
}

allom_stem <- function(df){
  lm(stem ~ dia_base + height, data=df)
}

allom_rot <- function(df){
  lm(root ~ leaf, data=df)
}

allom_ca <- function(df){
  lm(carea ~ dia_base + height, data=df)
}

allom_hgt <- function(df){
  lm(height ~ dia_base, data=df)
}

ca_sm     <- function(df){
  lm(stem ~ carea + height, data=df)
}

ca_agb   <- function(df){
  lm(leaf ~ carea + height, data=df)
}

ca_tm   <- function(df){
  lm(tmass ~ carea + height, data=df)
}

## species-specific models 

mod_df <- all_log %>% group_by(spcode) %>% nest()
mod_df <- mod_df %>% mutate(lf_mod = map(data, allom_lf),
                            hgt_mod = map(data,allom_hgt),
                            sm_mod = map(data,allom_stem),
                            rot_mod = map(data,allom_rot),
                            ca_mod = map(data,allom_ca),
                            sm_ca  = map(data,ca_sm),
                            agb_ca = map(data,ca_agb),
                            tm_ca  = map(data,ca_tm),
                            lf_pred = map2(.x=lf_mod, .y=data,~predict(object=.x,newdata=.y)),
                            hgt_pred = map2(.x=hgt_mod, .y=data,~predict(object=.x,newdata=.y)),
                            sm_pred = map2(.x=sm_mod, .y=data,~predict(object=.x,newdata=.y)),
                            rot_pred = map2(.x=rot_mod, .y=data,~predict(object=.x,newdata=.y)),
                            ca_pred = map2(.x=ca_mod, .y=data,~predict(object=.x,newdata=.y)),
                            smca_pred = map2(.x=sm_ca,.y=data,~predict(object=.x,newdata=.y)),
                            agbca_pred = map2(.x=agb_ca,.y=data,~predict(object=.x,newdata=.y)),
                            tmca_pred = map2(.x=tm_ca,.y=data,~predict(object=.x,newdata=.y))) 

##glance for model statistic
mod_df <- mod_df %>% mutate(lf = map(lf_mod,broom::glance),
                            hgt = map(hgt_mod,broom::glance),
                            sm = map(sm_mod, broom::glance),
                            rot = map(rot_mod, broom::glance),
                            ca = map(ca_mod,broom::glance),
                            smca = map(sm_ca,broom::glance),
                            agbca = map(agb_ca,broom::glance),
                            tmca = map(tm_ca,broom::glance)) %>% 
  unnest(c(lf,hgt,sm,rot,ca,smca,agbca,tmca),
         names_sep=".")

mod_df <- mod_df %>% mutate(lf_tidy = map(lf_mod, broom::tidy),
                            hgt_tidy = map(hgt_mod, broom::tidy),
                            sm_tidy = map(sm_mod, broom::tidy),
                            rot_tidy = map(rot_mod, broom::tidy),
                            ca_tidy = map(ca_mod, broom::tidy),
                            smca_tidy = map(sm_ca,broom::tidy),
                            agbca_tidy = map(agb_ca,broom::tidy),
                            tmca_tidy = map(tm_ca,broom::tidy))

## flat model parameter tables
rthgt_params <- mod_df %>% dplyr::select(spcode,data,hgt_tidy,
                                       rot_tidy,hgt_pred,rot_pred)          %>% 
  unnest(c(hgt_tidy,rot_tidy),names_sep=".", keep_empty=TRUE)

lsc_params <- mod_df %>% dplyr::select(spcode, data, sm_tidy,lf_tidy,
                                       ca_tidy,sm_pred,lf_pred,ca_pred)     %>% 
  unnest(c(sm_tidy,lf_tidy,ca_tidy),names_sep=".")

ca_params <- mod_df %>% dplyr::select(spcode,data,smca_tidy,agbca_tidy,
                                      tmca_tidy,smca_pred,agbca_pred,tmca_pred) %>% 
  unnest(c(smca_tidy,agbca_tidy,tmca_tidy),names_sep=".")


rthgt_df <- rthgt_params                               %>% 
            unnest(c(data,hgt_pred,
                     rot_pred),
                  keep_empty = TRUE)

lsc_df  <- lsc_params                                  %>% 
           unnest(c(data,sm_pred,
                    lf_pred,ca_pred),
                  keep_empty = TRUE)

ca_df   <- ca_params                                   %>% 
           unnest(c(data,smca_pred,agbca_pred,tmca_pred),
                  keep_empty = TRUE)


rthgt_intcpt <- rthgt_df %>% filter(hgt_tidy.term=="(Intercept)")
rthgt_slope  <- rthgt_df %>% filter(hgt_tidy.term=="dia_base")

lsc_intcpt <- lsc_df  %>% filter(sm_tidy.term=="(Intercept)")
lsc_slope <-  lsc_df  %>% filter(sm_tidy.term %in% c("dia_base","height"))

ca_intcpt  <- ca_df   %>% filter(smca_tidy.term=="(Intercept)")
ca_slope   <- ca_df   %>%  filter(smca_tidy.term %in% c("carea","height"))

rthgt_slope <- rthgt_slope                             %>% 
  dplyr::select(spcode,short,photo,growth,id,dia_base,
                height,root,repro,agb,
                rot_pred,hgt_pred,
                hgt_tidy.estimate,rot_tidy.estimate)   %>% 
         rename(hgt_slope = hgt_tidy.estimate,
                rot_slope = rot_tidy.estimate)

rthgt_intcpt <- rthgt_intcpt                            %>% 
                dplyr::select(spcode,short,photo,
                              growth,id,hgt_tidy.estimate,
                              rot_tidy.estimate)        %>% 
                rename(hgt_incpt = hgt_tidy.estimate,
                       rot_incpt = rot_tidy.estimate)

lsc_slope  <- lsc_slope                                   %>% 
             dplyr::select(spcode,short,photo,
                           growth,id,stem,leaf,
                           carea,sm_pred,lf_pred,
                           ca_pred,sm_tidy.estimate,
                           sm_tidy.term,lf_tidy.estimate,
                           lf_tidy.term,ca_tidy.estimate,
                           ca_tidy.term)                 %>% 
  pivot_wider(names_from  = sm_tidy.term,
              values_from = c(sm_tidy.estimate,
                              lf_tidy.estimate,
                              ca_tidy.estimate),
              names_sep   = ".")                         %>%
 rename(dbh_sm = sm_tidy.estimate.dia_base,
        hgt_sm = sm_tidy.estimate.height,
        dbh_lf = lf_tidy.estimate.dia_base,
        hgt_lf = lf_tidy.estimate.height,
        dbh_ca = ca_tidy.estimate.dia_base,
        hgt_ca = ca_tidy.estimate.height)               %>% 
  group_by(spcode)                                      %>% 
  mutate(dbh_sm = unique(dbh_sm[!is.na(dbh_sm)]),
         hgt_sm = unique(hgt_sm[!is.na(hgt_sm)]),
         dbh_lf = unique(dbh_lf[!is.na(dbh_lf)]),
         hgt_lf = unique(hgt_lf[!is.na(hgt_lf)]),
         dbh_ca = unique(dbh_ca[!is.na(dbh_ca)]),
         hgt_ca = unique(hgt_ca[!is.na(hgt_ca)]))        %>% 
  ungroup()                                              %>% 
  select(-c(lf_tidy.term,ca_tidy.term))                  %>% 
  distinct()



  

lsc_intcpt <- lsc_intcpt                                 %>% 
              dplyr::select(spcode,short,photo,growth,id,
                            sm_tidy.estimate,
                            lf_tidy.estimate, 
                            ca_tidy.estimate)            %>% 
  rename(sm_incpt =  sm_tidy.estimate,
         lf_incpt = lf_tidy.estimate,
         ca_incpt = ca_tidy.estimate)

lsc_all  = lsc_intcpt                                    %>%  
           left_join(lsc_slope,by=c("spcode","short",
                                    "photo","growth","id"))

obs_predic <- rthgt_slope                               %>% 
              left_join(rthgt_intcpt, 
                        by=c("spcode","short","photo",
                             "growth","id"))            %>% 
              ungroup()

obs_predic <- obs_predic                                %>% 
              left_join(lsc_all,by=c("spcode","short","photo",
                                     "growth","id"))    %>% 
              ungroup()


ca_slope  <- ca_slope                                   %>% 
  dplyr::select(spcode,short,photo,
                growth,id,stem,agb,tmass,
                carea,height,smca_pred,agbca_pred,
                tmca_pred,smca_tidy.estimate,
                smca_tidy.term,agbca_tidy.estimate,
                agbca_tidy.term,tmca_tidy.estimate,
                tmca_tidy.term)                 %>% 
  pivot_wider(names_from  = smca_tidy.term,
              values_from = c(smca_tidy.estimate,
                              agbca_tidy.estimate,
                              tmca_tidy.estimate),
              names_sep   = ".")                         %>%
  rename(ca_sm = smca_tidy.estimate.carea,
         hgt_sm = smca_tidy.estimate.height,
         ca_agb = agbca_tidy.estimate.carea,
         hgt_agb = agbca_tidy.estimate.height,
         ca_tm = tmca_tidy.estimate.carea,
         hgt_tm = tmca_tidy.estimate.height)               %>% 
  group_by(spcode)                                      %>% 
  mutate(ca_sm = unique(ca_sm[!is.na(ca_sm)]),
         hgt_sm = unique(hgt_sm[!is.na(hgt_sm)]),
         ca_agb = unique(ca_agb[!is.na(ca_agb)]),
         hgt_agb = unique(hgt_agb[!is.na(hgt_agb)]),
         ca_tm = unique(ca_tm[!is.na(ca_tm)]),
         hgt_tm = unique(hgt_tm[!is.na(hgt_tm)]))        %>% 
  ungroup()                                              %>% 
  select(-c(agbca_tidy.term,tmca_tidy.term))                  %>% 
  distinct()




ca_intcpt <- ca_intcpt                                 %>% 
  dplyr::select(spcode,short,photo,growth,id,
                smca_tidy.estimate,
                agbca_tidy.estimate, 
                tmca_tidy.estimate)            %>% 
  rename(sm_incpt =  smca_tidy.estimate,
         agb_incpt = agbca_tidy.estimate,
         tm_incpt = tmca_tidy.estimate)

casp_predic <- ca_slope %>% left_join(ca_intcpt,by=c("spcode","short","photo",
                                                       "growth","id")) %>% ungroup()

## let's correct for intercept, see ref: 
## https://www.jstor.org/stable/1937343?seq=1&cid=pdf-reference#references_tab_contents

obs_predic <- obs_predic %>% group_by(spcode) %>% 
  mutate(lf_n = length(!is.na(leaf)),
         sm_n = length(!is.na(stem)),
         rt_n = length(!is.na(root)),
         hgt_n = length(!is.na(height)),
         ca_n  = length(!is.na(carea)),
         lf_sme = (leaf-lf_pred)*(leaf-lf_pred),
         sm_sme = (stem-sm_pred)*(stem-sm_pred),
         rt_sme = (root-rot_pred)*(root-rot_pred),
         hgt_sme = (height-hgt_pred)*(height-hgt_pred),
         ca_sme = (carea-ca_pred)*(carea-ca_pred))

obs_predic <- obs_predic %>% mutate(lf_cf = sum(lf_sme, na.rm=TRUE)/(lf_n-2)/2,
                                    sm_cf = sum(sm_sme,na.rm=TRUE)/(sm_n-3)/2,
                                    rt_cf = sum(rt_sme, na.rm=TRUE)/(rt_n-2)/2,
                                    hgt_cf = sum(hgt_sme, na.rm=TRUE)/(hgt_n-2)/2,
                                    ca_cf = sum(ca_sme, na.rm=TRUE)/(ca_n-2)/2) %>% 
  ungroup()


obs_predic <- obs_predic %>% mutate(across(c("lf_incpt","hgt_incpt",
                                             "sm_incpt","rot_incpt",
                                             "ca_incpt","leaf","stem",
                                             "root","carea","dia_base",
                                             "height","lf_cf","sm_cf","rt_cf",
                                             "hgt_cf","ca_cf"),exp))

## update intercept
obs_predic <- obs_predic %>% mutate(lf_incpt = lf_incpt*lf_cf,
                                    hgt_incpt = hgt_incpt*hgt_cf,
                                    sm_incpt = sm_incpt*sm_cf,
                                    rot_incpt = rot_incpt*rt_cf,
                                    ca_incpt = ca_incpt*ca_cf)

## update prediction in un-logged form
obs_predic <- obs_predic %>% group_by(spcode)                  %>% 
  mutate(lf_pred = lf_incpt*(dia_base^dbh_lf)*(height^hgt_lf),
         hgt_pred = hgt_incpt*(dia_base^hgt_slope),
         sm_pred = sm_incpt*(dia_base^dbh_sm)*(height^hgt_sm),
         rt_pred = rot_incpt*(leaf^rot_slope),
         ca_pred = ca_incpt*(dia_base^dbh_ca)*(height^hgt_ca)) %>% 
         ungroup()

## convert carbon to g

obs_predic <- obs_predic %>% mutate(leaf = leaf*kg_2_g,
                                    stem = stem*kg_2_g,
                                    root = root*kg_2_g,
                                    lf_pred = lf_pred*kg_2_g,
                                    sm_pred = sm_pred*kg_2_g,
                                    rot_pred = rot_pred*kg_2_g)

allom_params = obs_predic             %>% 
               select(spcode
                     ,short
                     ,photo
                     ,growth
                     ,hgt_slope
                     ,rot_slope
                     ,hgt_incpt
                     ,rot_incpt
                     ,sm_incpt
                     ,lf_incpt
                     ,ca_incpt
                     ,dbh_sm
                     ,hgt_sm
                     ,dbh_lf
                     ,hgt_lf
                     ,dbh_ca
                     ,hgt_ca)        %>% 
  distinct()

### canopy area based models ###
casp_predic <- casp_predic %>% group_by(spcode) %>% 
  mutate(
         sm_n = length(!is.na(stem)),
         agb_n = length(!is.na(agb)),
         tm_n  = length(!is.na(tmass)),
         sm_sme = (stem-smca_pred)*(stem-smca_pred),
         agb_sme = (agb-agbca_pred)*(agb-agbca_pred),
         tm_sme = (tmass-tmca_pred)*(tmass-tmca_pred))

casp_predic <- casp_predic %>% mutate(
                                    sm_cf = sum(sm_sme,na.rm=TRUE)/(sm_n-3)/2,
                                    agb_cf = sum(agb_sme, na.rm=TRUE)/(agb_n-2)/2,
                                    tm_cf = sum(tm_sme, na.rm=TRUE)/(tm_n-2)/2) %>% 
  ungroup()


casp_predic <- casp_predic %>% mutate(across(c("sm_incpt","agb_incpt",
                                             "tm_incpt","stem","agb",
                                             "tmass","carea",
                                             "height","sm_cf","agb_cf",
                                             "tm_cf"),exp))

## update intercept
casp_predic <- casp_predic %>% mutate(
                                    sm_incpt = sm_incpt*sm_cf,
                                    agb_incpt = agb_incpt*agb_cf,
                                    tm_incpt = tm_incpt*tm_cf)

## update prediction in un-logged form
casp_predic <- casp_predic %>% group_by(spcode)                  %>% 
  mutate(smca_pred = sm_incpt*(carea^ca_sm)*(height^hgt_sm),
         agbca_pred = agb_incpt*(carea^ca_agb)*(height^hgt_agb),
         tmca_pred = tm_incpt*(carea^ca_tm)*(height^hgt_tm)) %>% 
  ungroup()

## convert carbon to g

casp_predic <- casp_predic %>% mutate(
                                    stem = stem*kg_2_g,
                                    agb  = agb*kg_2_g,
                                    tmass = tmass*kg_2_g,
                                    smca_pred = smca_pred*kg_2_g,
                                    agbca_pred = agbca_pred*kg_2_g,
                                    tmca_pred = tmca_pred*kg_2_g)

casp_params = casp_predic             %>% 
  select(spcode
         ,short
         ,photo
         ,growth
         ,ca_sm
         ,hgt_sm
         ,sm_incpt
         ,ca_agb
         ,hgt_agb
         ,agb_incpt
         ,ca_tm
         ,hgt_tm
         ,tm_incpt)        %>% 
  distinct()


#### fit allometric equations for function types ####

pf_df <- all_log %>% group_by(growth,photo) %>% nest()
pf_df <- pf_df %>% mutate(  lf_mod  = map(data, allom_lf),
                            hgt_mod = map(data,allom_hgt),
                            sm_mod  = map(data,allom_stem),
                            rot_mod = map(data,allom_rot),
                            ca_mod  = map(data,allom_ca),
                            sm_ca   = map(data,ca_sm),
                            agb_ca  = map(data,ca_agb),
                            tm_ca   = map(data,ca_tm),
                            lf_pred = map2(.x=lf_mod, .y=data,~predict(object=.x,newdata=.y)),
                            hgt_pred = map2(.x=hgt_mod, .y=data,~predict(object=.x,newdata=.y)),
                            sm_pred = map2(.x=sm_mod, .y=data,~predict(object=.x,newdata=.y)),
                            rot_pred = map2(.x=rot_mod, .y=data,~predict(object=.x,newdata=.y)),
                            ca_pred = map2(.x=ca_mod, .y=data,~predict(object=.x,newdata=.y)),
                            smca_pred = map2(.x=sm_ca, .y=data,~predict(object=.x,newdata=.y)),
                            agbca_pred = map2(.x=agb_ca,.y=data,~predict(object=.x,newdata=.y)),
                            tmca_pred  = map2(.x=tm_ca,.y=data,~predict(object=.x,newdata=.y))) 

##glance for model statistic
pf_df <- pf_df %>% mutate(  lf = map(lf_mod,broom::glance),
                            hgt = map(hgt_mod,broom::glance),
                            sm = map(sm_mod, broom::glance),
                            rot = map(rot_mod, broom::glance),
                            ca = map(ca_mod,broom::glance),
                            smca = map(sm_ca,broom::glance),
                            agbca = map(agb_ca,broom::glance),
                            tmca  = map(tm_ca,broom::glance)) %>% 
  unnest(c(lf,hgt,sm,rot,ca,smca,agbca,tmca),
         names_sep=".")

pf_df <- pf_df %>% mutate(  lf_tidy = map(lf_mod, broom::tidy),
                            hgt_tidy = map(hgt_mod, broom::tidy),
                            sm_tidy = map(sm_mod, broom::tidy),
                            rot_tidy = map(rot_mod, broom::tidy),
                            ca_tidy = map(ca_mod, broom::tidy),
                            smca_tidy = map(sm_ca, broom::tidy),
                            agbca_tidy = map(agb_ca,broom::tidy),
                            tmca_tidy  = map(tm_ca, broom::tidy))

## flat model parameter tables
rthgt_pfparams <- pf_df %>% dplyr::select(growth,photo,data,hgt_tidy,
                                         rot_tidy,hgt_pred,rot_pred)          %>% 
  unnest(c(hgt_tidy,rot_tidy),names_sep=".", keep_empty=TRUE)

lsc_pfparams <- pf_df %>% dplyr::select(growth, photo, data, sm_tidy,lf_tidy,
                                       ca_tidy,sm_pred,lf_pred,ca_pred)     %>% 
  unnest(c(sm_tidy,lf_tidy,ca_tidy),names_sep=".")

camass_pfparams  <- pf_df %>% dplyr::select(growth, photo, data, smca_tidy, 
                                            agbca_tidy,tmca_tidy, smca_pred,
                                            agbca_pred,tmca_pred)           %>% 
                    unnest(c(smca_tidy,agbca_tidy,tmca_tidy),names_sep=".")


rthgt_pfdf <- rthgt_pfparams                            %>% 
  unnest(c(data,hgt_pred,
           rot_pred),
         keep_empty = TRUE)

lsc_pfdf  <- lsc_pfparams                               %>% 
  unnest(c(data,sm_pred,
           lf_pred,ca_pred),
         keep_empty = TRUE)

camass_pfdf <- camass_pfparams                          %>% 
  unnest(c(data,smca_pred,
           agbca_pred, tmca_pred),
         keep_empty=TRUE)
  


rthgt_pfintcpt <- rthgt_pfdf %>% filter(hgt_tidy.term=="(Intercept)")
rthgt_pfslope  <- rthgt_pfdf %>% filter(hgt_tidy.term=="dia_base")

lsc_pfintcpt <- lsc_pfdf  %>% filter(sm_tidy.term=="(Intercept)")
lsc_pfslope <-  lsc_pfdf  %>% filter(sm_tidy.term %in% c("dia_base","height"))

camass_pfintcpt <- camass_pfdf %>% filter(smca_tidy.term=="(Intercept)")
camass_pfslope  <- camass_pfdf %>% filter(smca_tidy.term %in% c("carea","height"))

rthgt_pfslope <- rthgt_pfslope                          %>% 
  dplyr::select(spcode,short,photo,growth,id,dia_base,
                height,root,repro,agb,
                rot_pred,hgt_pred,
                hgt_tidy.estimate,rot_tidy.estimate)   %>% 
  rename(hgt_slope = hgt_tidy.estimate,
         rot_slope = rot_tidy.estimate)

rthgt_pfintcpt <- rthgt_pfintcpt                       %>% 
  dplyr::select(spcode,short,photo,
                growth,id,hgt_tidy.estimate,
                rot_tidy.estimate)        %>% 
  rename(hgt_incpt = hgt_tidy.estimate,
         rot_incpt = rot_tidy.estimate)

lsc_pfslope  <- lsc_pfslope                             %>% 
  dplyr::select(spcode,short,photo,
                growth,id,stem,leaf,
                carea,sm_pred,lf_pred,
                ca_pred,sm_tidy.estimate,
                sm_tidy.term,lf_tidy.estimate,
                lf_tidy.term,ca_tidy.estimate,
                ca_tidy.term)                 %>% 
  pivot_wider(names_from  = sm_tidy.term,
              values_from = c(sm_tidy.estimate,
                              lf_tidy.estimate,
                              ca_tidy.estimate),
              names_sep   = ".")                         %>%
  rename(dbh_sm = sm_tidy.estimate.dia_base,
         hgt_sm = sm_tidy.estimate.height,
         dbh_lf = lf_tidy.estimate.dia_base,
         hgt_lf = lf_tidy.estimate.height,
         dbh_ca = ca_tidy.estimate.dia_base,
         hgt_ca = ca_tidy.estimate.height)               %>% 
  group_by(growth,photo)                                %>% 
  mutate(dbh_sm = unique(dbh_sm[!is.na(dbh_sm)]),
         hgt_sm = unique(hgt_sm[!is.na(hgt_sm)]),
         dbh_lf = unique(dbh_lf[!is.na(dbh_lf)]),
         hgt_lf = unique(hgt_lf[!is.na(hgt_lf)]),
         dbh_ca = unique(dbh_ca[!is.na(dbh_ca)]),
         hgt_ca = unique(hgt_ca[!is.na(hgt_ca)]))        %>% 
  ungroup()                                              %>% 
  select(-c(lf_tidy.term,ca_tidy.term))                  %>% 
  distinct()



lsc_pfintcpt <- lsc_pfintcpt                              %>% 
  dplyr::select(spcode,short,photo,growth,id,
                sm_tidy.estimate,
                lf_tidy.estimate, 
                ca_tidy.estimate)            %>% 
  rename(sm_incpt =  sm_tidy.estimate,
         lf_incpt = lf_tidy.estimate,
         ca_incpt = ca_tidy.estimate)

lsc_pfall  = lsc_pfintcpt                               %>%  
  left_join(lsc_pfslope,by=c("spcode","short",
                           "photo","growth","id"))

camass_pfslope  <- camass_pfslope                       %>% 
  dplyr::select(spcode,short,photo,
                growth,id,stem,agb,tmass,
                carea,height,smca_pred,
                agbca_pred,tmca_pred, 
                smca_tidy.estimate,
                smca_tidy.term,agbca_tidy.estimate,
                agbca_tidy.term,
                tmca_tidy.term,tmca_tidy.estimate)      %>% 
  pivot_wider(names_from  = smca_tidy.term,
              values_from = c(smca_tidy.estimate,
                              agbca_tidy.estimate,
                              tmca_tidy.estimate),
              names_sep   = ".")                         %>%
  rename(ca_sm  = smca_tidy.estimate.carea,
         hgt_sm = smca_tidy.estimate.height,
         ca_agb = agbca_tidy.estimate.carea,
         hgt_agb = agbca_tidy.estimate.height,
         ca_tm   = tmca_tidy.estimate.carea,
         hgt_tm  = tmca_tidy.estimate.height)             %>% 
  group_by(growth,photo)                                  %>% 
  mutate(ca_sm  = unique(ca_sm[!is.na(ca_sm)]),
         hgt_sm  = unique(hgt_sm[!is.na(hgt_sm)]),
         ca_agb  = unique(ca_agb[!is.na(ca_agb)]),
         hgt_agb = unique(hgt_agb[!is.na(hgt_agb)]),
         ca_tm   = unique(ca_tm[!is.na(ca_tm)]),
         hgt_tm  = unique(hgt_tm[!is.na(hgt_tm)]))        %>% 
  ungroup()  %>% 
  select(-c(agbca_tidy.term,tmca_tidy.term))              %>% 
  distinct()   

camass_pfintcpt <- camass_pfintcpt                        %>% 
  dplyr::select(spcode,short,photo,growth,id,
                smca_tidy.estimate,
                agbca_tidy.estimate, 
                tmca_tidy.estimate)            %>% 
  rename(sm_incpt =  smca_tidy.estimate,
         agb_incpt = agbca_tidy.estimate,
         tm_incpt =  tmca_tidy.estimate)



obs_pfpredic <- rthgt_pfslope                               %>% 
  left_join(rthgt_pfintcpt, 
            by=c("spcode","short","photo",
                 "growth","id"))            %>% 
  ungroup()

obs_pfpredic <- obs_pfpredic                                %>% 
  left_join(lsc_pfall,by=c("spcode","short","photo",
                         "growth","id"))    %>% 
  ungroup()

camass_pfpredic <- camass_pfslope               %>% 
  left_join(camass_pfintcpt,
            by = c("spcode","short",
                   "photo","growth","id"))

## let's correct for intercept, see ref: 


### DBH and plant height based models ###
obs_pfpredic <- obs_pfpredic %>% group_by(growth,photo) %>% 
  mutate(lf_n = length(!is.na(leaf)),
         sm_n = length(!is.na(stem)),
         rt_n = length(!is.na(root)),
         hgt_n = length(!is.na(height)),
         ca_n  = length(!is.na(carea)),
         lf_sme = (leaf-lf_pred)*(leaf-lf_pred),
         sm_sme = (stem-sm_pred)*(stem-sm_pred),
         rt_sme = (root-rot_pred)*(root-rot_pred),
         hgt_sme = (height-hgt_pred)*(height-hgt_pred),
         ca_sme = (carea-ca_pred)*(carea-ca_pred))

obs_pfpredic <- obs_pfpredic %>% mutate(lf_cf = sum(lf_sme, na.rm=TRUE)/(lf_n-2)/2,
                                    sm_cf = sum(sm_sme,na.rm=TRUE)/(sm_n-3)/2,
                                    rt_cf = sum(rt_sme, na.rm=TRUE)/(rt_n-2)/2,
                                    hgt_cf = sum(hgt_sme, na.rm=TRUE)/(hgt_n-2)/2,
                                    ca_cf = sum(ca_sme, na.rm=TRUE)/(ca_n-2)/2) %>% 
  ungroup()


obs_pfpredic <- obs_pfpredic %>% mutate(across(c("lf_incpt","hgt_incpt",
                                             "sm_incpt","rot_incpt",
                                             "ca_incpt","leaf","stem",
                                             "root","carea","dia_base",
                                             "height","lf_cf","sm_cf","rt_cf",
                                             "hgt_cf","ca_cf"),exp))

## update intercept
obs_pfpredic <- obs_pfpredic %>% mutate(lf_incpt = lf_incpt*lf_cf,
                                    hgt_incpt = hgt_incpt*hgt_cf,
                                    sm_incpt = sm_incpt*sm_cf,
                                    rot_incpt = rot_incpt*rt_cf,
                                    ca_incpt = ca_incpt*ca_cf)




## update prediction in un-logged form
obs_pfpredic <- obs_pfpredic %>% group_by(growth,photo)             %>% 
  mutate(lf_pred = lf_incpt*(dia_base^dbh_lf)*(height^hgt_lf),
         hgt_pred = hgt_incpt*(dia_base^hgt_slope),
         sm_pred = sm_incpt*(dia_base^dbh_sm)*(height^hgt_sm),
         rt_pred = rot_incpt*(leaf^rot_slope),
         ca_pred = ca_incpt*(dia_base^dbh_ca)*(height^hgt_ca)) %>% 
  ungroup()

## convert carbon to g

obs_pfpredic <- obs_pfpredic %>% mutate(leaf = leaf*kg_2_g,
                                    stem = stem*kg_2_g,
                                    root = root*kg_2_g,
                                    lf_pred = lf_pred*kg_2_g,
                                    sm_pred = sm_pred*kg_2_g,
                                    rot_pred = rot_pred*kg_2_g)

allom_pfparams = obs_pfpredic             %>% 
                 select(spcode
                ,short
                ,photo
                ,growth
                ,hgt_slope
                ,rot_slope
                ,hgt_incpt
                ,rot_incpt
                ,sm_incpt
                ,lf_incpt
                ,ca_incpt
                ,dbh_sm
                ,hgt_sm
                ,dbh_lf
                ,hgt_lf
                ,dbh_ca
                ,hgt_ca)        %>% 
                distinct()

### canopy area based models ###

camass_pfpredic <- camass_pfpredic %>% group_by(growth,photo) %>% 
  mutate(sm_n = length(!is.na(stem)),
         agb_n = length(!is.na(agb)),
         tm_n = length(!is.na(tmass)),
         sm_sme  = (stem-smca_pred)*(stem-smca_pred),
         agb_sme = (agb-agbca_pred)*(agb-agbca_pred),
         tm_sme  = (tmass-tmca_pred)*(tmass-tmca_pred))

camass_pfpredic <- camass_pfpredic %>% mutate(
                                        sm_cf = sum(sm_sme,na.rm=TRUE)/(sm_n-3)/2,
                                        agb_cf = sum(agb_sme, na.rm=TRUE)/(agb_n-2)/2,
                                        tm_cf  = sum(tm_sme, na.rm=TRUE)/(tm_n-2)/2) %>% 
  ungroup()


camass_pfpredic <- camass_pfpredic %>% mutate(across(c(
                                                 "sm_incpt","agb_incpt",
                                                 "tm_incpt","stem", "agb",
                                                 "tmass","carea","height",
                                                 "sm_cf","agb_cf","tm_cf"),exp))

## update intercept
camass_pfpredic <- camass_pfpredic %>% mutate(sm_incpt = sm_incpt*sm_cf,
                                        agb_incpt = agb_incpt*agb_cf,
                                        tm_incpt = tm_incpt*tm_cf)




## update prediction in un-logged form
camass_pfpredic <- camass_pfpredic %>% group_by(growth,photo)             %>% 
  mutate(smca_pred = sm_incpt*(carea^ca_sm)*(height^hgt_sm),
         agbca_pred = agb_incpt*(carea^ca_agb)*(height^hgt_agb),
         tmca_pred = tm_incpt*(carea^ca_tm)*(height^hgt_tm)) %>% 
  ungroup()

## convert carbon to g

camass_pfpredic <- camass_pfpredic %>% mutate(
                                        stem = stem*kg_2_g,
                                        agb  = agb*kg_2_g,
                                        tmass = tmass*kg_2_g,
                                        smca_pred = smca_pred*kg_2_g,
                                        agbca_pred = agbca_pred*kg_2_g,
                                        tmca_pred = tmca_pred*kg_2_g)

camass_pfparams = camass_pfpredic             %>% 
  select(spcode
         ,short
         ,photo
         ,growth
         ,ca_sm
         ,hgt_sm
         ,sm_incpt
         ,ca_agb
         ,hgt_agb
         ,agb_incpt
         ,ca_tm
         ,hgt_tm
         ,tm_incpt)        %>% 
  distinct()



#### fit allometric equations for all grasses ####

grass_df <- all_log %>% mutate(pft = "grass") %>% group_by(pft) %>% nest()
grass_df <- grass_df %>% mutate(lf_mod = map(data, allom_lf),
                                hgt_mod = map(data,allom_hgt),
                                sm_mod = map(data,allom_stem),
                                rot_mod = map(data,allom_rot),
                                ca_mod = map(data,allom_ca),
                                sm_ca  = map(data,ca_sm),
                                agb_ca = map(data,ca_agb),
                                tm_ca  = map(data,ca_tm),
                                lf_pred = map2(.x=lf_mod, .y=data,~predict(object=.x,newdata=.y)),
                                hgt_pred = map2(.x=hgt_mod, .y=data,~predict(object=.x,newdata=.y)),
                                sm_pred = map2(.x=sm_mod, .y=data,~predict(object=.x,newdata=.y)),
                                rot_pred = map2(.x=rot_mod, .y=data,~predict(object=.x,newdata=.y)),
                                ca_pred = map2(.x=ca_mod, .y=data,~predict(object=.x,newdata=.y)),
                                smca_pred = map2(.x=sm_ca, .y=data, ~predict(object=.x,newdata=.y)),
                                agbca_pred = map2(.x=agb_ca, .y=data, ~predict(object=.x,newdata=.y)),
                                tmca_pred  = map2(.x=tm_ca, .y=data, ~predict(object=.x,newdata=.y))) 

##glance for model statistic
grass_df <- grass_df %>% mutate(lf = map(lf_mod,broom::glance),
                                hgt = map(hgt_mod,broom::glance),
                                sm = map(sm_mod, broom::glance),
                                rot = map(rot_mod, broom::glance),
                                ca = map(ca_mod,broom::glance),
                                smca = map(sm_ca,broom::glance),
                                agbca = map(agb_ca,broom::glance),
                                tmca  = map(tm_ca,broom::glance)) %>% 
                         unnest(c(lf,hgt,sm,rot,ca,smca,agbca,tmca), names_sep=".")

grass_df <- grass_df %>% mutate(lf_tidy = map(lf_mod, broom::tidy),
                                hgt_tidy = map(hgt_mod, broom::tidy),
                                sm_tidy = map(sm_mod, broom::tidy),
                                rot_tidy = map(rot_mod, broom::tidy),
                                ca_tidy = map(ca_mod, broom::tidy),
                                smca_tidy = map(sm_ca,broom::tidy),
                                agbca_tidy = map(agb_ca,broom::tidy),
                                tmca_tidy  = map(tm_ca,broom::tidy))

## flat model parameter tables
rthgt_grassparams <- grass_df %>% dplyr::select(pft,data,hgt_tidy,
                                                rot_tidy,hgt_pred,rot_pred) %>% 
                                  unnest(c(hgt_tidy,rot_tidy),names_sep=".", 
                                           keep_empty=TRUE)

lsc_grassparams <- grass_df %>% dplyr::select(pft, data, sm_tidy,lf_tidy,
                                           ca_tidy,sm_pred,lf_pred,ca_pred) %>% 
                             unnest(c(sm_tidy,lf_tidy,ca_tidy),names_sep=".")

ca_grassparams  <- grass_df %>% dplyr::select(pft,data,smca_tidy,agbca_tidy,
                                              tmca_tidy,smca_pred,agbca_pred,tmca_pred) %>% 
                                unnest(c(smca_tidy,agbca_tidy,tmca_tidy),names_sep=".")


rthgt_grassdf <- rthgt_grassparams                      %>% 
                 unnest(c(data,hgt_pred,rot_pred),
                        keep_empty = TRUE)

lsc_grassdf  <- lsc_grassparams                         %>% 
                unnest(c(data,sm_pred,lf_pred,ca_pred),
                       keep_empty = TRUE)
ca_grassdf   <- ca_grassparams                          %>% 
                unnest(c(data,smca_pred,agbca_pred,tmca_pred),
                       keep_empty=TRUE)


rthgt_grassintcpt <- rthgt_grassdf %>% filter(hgt_tidy.term=="(Intercept)")
rthgt_grassslope  <- rthgt_grassdf %>% filter(hgt_tidy.term=="dia_base")

lsc_grassintcpt <- lsc_grassdf  %>% filter(sm_tidy.term=="(Intercept)")
lsc_grassslope <-  lsc_grassdf  %>% filter(sm_tidy.term %in% c("dia_base","height"))

ca_grassintcpt  <- ca_grassdf   %>% filter(smca_tidy.term=="(Intercept)")
ca_grassslope   <- ca_grassdf   %>% filter(smca_tidy.term %in% c("carea","height"))

rthgt_grassslope <- rthgt_grassslope                   %>% 
  dplyr::select(pft,id,dia_base,
                height,root,repro,agb,
                rot_pred,hgt_pred,
                hgt_tidy.estimate,rot_tidy.estimate)   %>% 
  rename(hgt_slope = hgt_tidy.estimate,
         rot_slope = rot_tidy.estimate)

rthgt_grassintcpt <- rthgt_grassintcpt                 %>% 
  dplyr::select(pft,id,hgt_tidy.estimate,
                rot_tidy.estimate)                     %>% 
  rename(hgt_incpt = hgt_tidy.estimate,
         rot_incpt = rot_tidy.estimate)

lsc_grassslope  <- lsc_grassslope                      %>% 
  dplyr::select(pft,id,stem,leaf,carea,sm_pred,lf_pred,
                ca_pred,sm_tidy.estimate,
                sm_tidy.term,lf_tidy.estimate,
                lf_tidy.term,ca_tidy.estimate,
                ca_tidy.term)                          %>% 
  pivot_wider(names_from  = sm_tidy.term,
              values_from = c(sm_tidy.estimate,
                              lf_tidy.estimate,
                              ca_tidy.estimate),
              names_sep   = ".")                         %>%
  rename(dbh_sm = sm_tidy.estimate.dia_base,
         hgt_sm = sm_tidy.estimate.height,
         dbh_lf = lf_tidy.estimate.dia_base,
         hgt_lf = lf_tidy.estimate.height,
         dbh_ca = ca_tidy.estimate.dia_base,
         hgt_ca = ca_tidy.estimate.height)               %>% 
  group_by(pft)                                          %>% 
  mutate(dbh_sm = unique(dbh_sm[!is.na(dbh_sm)]),
         hgt_sm = unique(hgt_sm[!is.na(hgt_sm)]),
         dbh_lf = unique(dbh_lf[!is.na(dbh_lf)]),
         hgt_lf = unique(hgt_lf[!is.na(hgt_lf)]),
         dbh_ca = unique(dbh_ca[!is.na(dbh_ca)]),
         hgt_ca = unique(hgt_ca[!is.na(hgt_ca)]))        %>% 
  ungroup()                                              %>% 
  select(-c(lf_tidy.term,ca_tidy.term))                  %>% 
  distinct()





lsc_grassintcpt <- lsc_grassintcpt                      %>% 
  dplyr::select(pft,id,sm_tidy.estimate,
                lf_tidy.estimate, 
                ca_tidy.estimate)                       %>% 
  rename(sm_incpt =  sm_tidy.estimate,
         lf_incpt = lf_tidy.estimate,
         ca_incpt = ca_tidy.estimate)

lsc_grassall  = lsc_grassintcpt                             %>%  
  left_join(lsc_grassslope,by=c("pft","id"))

obs_grasspredic <- rthgt_grassslope                         %>% 
                   left_join(rthgt_grassintcpt, 
                             by=c("pft","id"))              %>% 
                   ungroup()

obs_grasspredic <- obs_grasspredic                          %>% 
  left_join(lsc_grassall,by=c("pft","id"))                  %>% 
  ungroup()


ca_grassslope  <- ca_grassslope                      %>% 
  dplyr::select(pft,id,stem,agb,tmass,carea,height,
                smca_pred,agbca_pred,tmca_pred,
                smca_tidy.estimate,
                smca_tidy.term,agbca_tidy.estimate,
                agbca_tidy.term,tmca_tidy.estimate,
                tmca_tidy.term)                       %>% 
  pivot_wider(names_from  = smca_tidy.term,
              values_from = c(smca_tidy.estimate,
                              agbca_tidy.estimate,
                              tmca_tidy.estimate),
              names_sep   = ".")                         %>%
  rename(ca_sm = smca_tidy.estimate.carea,
         hgt_sm = smca_tidy.estimate.height,
         ca_agb = agbca_tidy.estimate.carea,
         hgt_agb = agbca_tidy.estimate.height,
         ca_tm = tmca_tidy.estimate.carea,
         hgt_tm = tmca_tidy.estimate.height)             %>% 
  group_by(pft)                                          %>% 
  mutate(ca_sm = unique(ca_sm[!is.na(ca_sm)]),
         hgt_sm = unique(hgt_sm[!is.na(hgt_sm)]),
         ca_agb = unique(ca_agb[!is.na(ca_agb)]),
         hgt_agb = unique(hgt_agb[!is.na(hgt_agb)]),
         ca_tm = unique(ca_tm[!is.na(ca_tm)]),
         hgt_tm = unique(hgt_tm[!is.na(hgt_tm)]))        %>% 
  ungroup()                                              %>% 
  select(-c(agbca_tidy.term,tmca_tidy.term))                  %>% 
  distinct()




ca_grassintcpt <- ca_grassintcpt                      %>% 
  dplyr::select(pft,id,smca_tidy.estimate,
                agbca_tidy.estimate, 
                tmca_tidy.estimate)                       %>% 
  rename(sm_incpt =  smca_tidy.estimate,
         agb_incpt = agbca_tidy.estimate,
         tm_incpt = tmca_tidy.estimate)

ca_grasspredic <- ca_grassslope                         %>% 
  left_join(ca_grassintcpt, 
            by=c("pft","id"))              %>% 
  ungroup()



## correct for intercept: 

obs_grasspredic <- obs_grasspredic %>% group_by(pft) %>% 
  mutate(lf_n = length(!is.na(leaf)),
         sm_n = length(!is.na(stem)),
         rt_n = length(!is.na(root)),
         hgt_n = length(!is.na(height)),
         ca_n  = length(!is.na(carea)),
         lf_sme = (leaf-lf_pred)*(leaf-lf_pred),
         sm_sme = (stem-sm_pred)*(stem-sm_pred),
         rt_sme = (root-rot_pred)*(root-rot_pred),
         hgt_sme = (height-hgt_pred)*(height-hgt_pred),
         ca_sme = (carea-ca_pred)*(carea-ca_pred))

obs_grasspredic <- obs_grasspredic %>% mutate(lf_cf = sum(lf_sme, na.rm=TRUE)/(lf_n-2)/2,
                                              sm_cf = sum(sm_sme,na.rm=TRUE)/(sm_n-3)/2,
                                              rt_cf = sum(rt_sme, na.rm=TRUE)/(rt_n-2)/2,
                                              hgt_cf = sum(hgt_sme, na.rm=TRUE)/(hgt_n-2)/2,
                                              ca_cf = sum(ca_sme, na.rm=TRUE)/(ca_n-2)/2) %>% 
                                      ungroup()


obs_grasspredic <- obs_grasspredic %>% mutate(across(c("lf_incpt","hgt_incpt",
                                                       "sm_incpt","rot_incpt",
                                                       "ca_incpt","leaf","stem",
                                                       "root","carea","dia_base",
                                                       "height","lf_cf","sm_cf","rt_cf",
                                                       "hgt_cf","ca_cf"),exp))

## update intercept
obs_grasspredic <- obs_grasspredic %>% mutate(lf_incpt = lf_incpt*lf_cf,
                                              hgt_incpt = hgt_incpt*hgt_cf,
                                              sm_incpt = sm_incpt*sm_cf,
                                              rot_incpt = rot_incpt*rt_cf,
                                              ca_incpt = ca_incpt*ca_cf)

## update prediction in un-logged form
obs_grasspredic <- obs_grasspredic                                               %>%   
                   mutate(lf_pred = lf_incpt*(dia_base^dbh_lf)*(height^hgt_lf),
                          hgt_pred = hgt_incpt*(dia_base^hgt_slope),
                          sm_pred = sm_incpt*(dia_base^dbh_sm)*(height^hgt_sm),
                          rt_pred = rot_incpt*(leaf^rot_slope),
                          ca_pred = ca_incpt*(dia_base^dbh_ca)*(height^hgt_ca)) 
 
## convert carbon to g

obs_grasspredic <- obs_grasspredic %>% mutate(leaf = leaf*kg_2_g,
                                              stem = stem*kg_2_g,
                                              root = root*kg_2_g,
                                              lf_pred = lf_pred*kg_2_g,
                                              sm_pred = sm_pred*kg_2_g,
                                              rot_pred = rot_pred*kg_2_g)

allom_grassparams = obs_grasspredic             %>% 
                    select(pft
                          ,hgt_slope
                          ,rot_slope
                          ,hgt_incpt
                          ,rot_incpt
                          ,sm_incpt
                          ,lf_incpt
                          ,ca_incpt
                          ,dbh_sm
                          ,hgt_sm
                          ,dbh_lf
                          ,hgt_lf
                          ,dbh_ca
                          ,hgt_ca)        %>% 
                   distinct()


#### fit allometric models for all grasses using canopy area and height ####

ca_grasspredic     <- ca_grasspredic %>% group_by(pft) %>% 
  mutate(sm_n = length(!is.na(stem)),
         agb_n = length(!is.na(agb)),
         tm_n = length(!is.na(tmass)),
         sm_sme  = (stem-smca_pred)*(stem-smca_pred),
         agb_sme = (agb-agbca_pred)*(agb-agbca_pred),
         tm_sme  = (tmass-tmca_pred)*(tmass-tmca_pred))

ca_grasspredic <- ca_grasspredic %>% mutate(
  sm_cf = sum(sm_sme,na.rm=TRUE)/(sm_n-3)/2,
  agb_cf = sum(agb_sme, na.rm=TRUE)/(agb_n-2)/2,
  tm_cf  = sum(tm_sme, na.rm=TRUE)/(tm_n-2)/2) %>% 
  ungroup()


ca_grasspredic <- ca_grasspredic %>% mutate(across(c(
  "sm_incpt","agb_incpt",
  "tm_incpt","stem", "agb",
  "tmass","carea","height",
  "sm_cf","agb_cf","tm_cf"),exp))

## update intercept
ca_grasspredic <- ca_grasspredic %>% mutate(sm_incpt = sm_incpt*sm_cf,
                                              agb_incpt = agb_incpt*agb_cf,
                                              tm_incpt = tm_incpt*tm_cf)




## update prediction in un-logged form
ca_grasspredic <- ca_grasspredic %>% group_by(pft)             %>% 
  mutate(smca_pred = sm_incpt*(carea^ca_sm)*(height^hgt_sm),
         agbca_pred = agb_incpt*(carea^ca_agb)*(height^hgt_agb),
         tmca_pred = tm_incpt*(carea^ca_tm)*(height^hgt_tm)) %>% 
  ungroup()

## convert carbon to g

ca_grasspredic <- ca_grasspredic %>% mutate(
  stem = stem*kg_2_g,
  agb  = agb*kg_2_g,
  tmass = tmass*kg_2_g,
  smca_pred = smca_pred*kg_2_g,
  agbca_pred = agbca_pred*kg_2_g,
  tmca_pred = tmca_pred*kg_2_g)

ca_grassparams = ca_grasspredic             %>% 
  select( pft
         ,ca_sm
         ,hgt_sm
         ,sm_incpt
         ,ca_agb
         ,hgt_agb
         ,agb_incpt
         ,ca_tm
         ,hgt_tm
         ,tm_incpt)        %>% 
  distinct()

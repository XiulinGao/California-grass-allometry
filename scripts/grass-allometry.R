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
       filter(spcode!="POSE")
        
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
             ,lf2sm      = leaf_c/stem_c
             ,ca2ba      = barea/carea
             ,tlarea     = leaf*sla)


all_log <- all                                                  %>% 
           mutate(across(c("stem","stem_c","root","root_c"    , 
                           "leaf","leaf_c","repro","repro_c"  ,
                           "agb","agb_c","tmass","tmass_c"    ,
                           "dia_base","dia_can","length"      ,
                           "height","barea","carea","lf2sm"   ,
                           "ca2ba","sla","tlarea")            ,
                           ~log(.x)))


####### 1. Allometry model selection and species-specific variation #######

#### 1a.height allometry model ####

hgt_null  <- lm(height ~ dia_base, all_log)
hgt_full  <- lm(height ~ dia_base*spcode,all_log)
AIC(hgt_null,hgt_full)
hgt_poshoc  <- emmeans(hgt_full, pairwise ~ dia_base*spcode)


hgt_allom <- hgt_full
summary(hgt_allom)
anova(hgt_allom)

#### 1b.canopy area allometry ####

ca_null1  <- lm(carea ~ dia_base + height, all_log)
ca_null2  <- lm(carea ~ dia_base, all_log)
ca_null3  <- lm(carea ~ height, all_log)
ca_full1  <- lm(carea ~ (dia_base + height)*spcode, all_log)
ca_full2  <- lm(carea ~ dia_base*spcode, all_log)
ca_full3  <- lm(carea ~ height*spcode, all_log)
AIC(ca_null1,ca_null2,ca_null3,ca_full1,ca_full2,ca_full3)
ca_poshoc <- emmeans(ca_full1, pairwise ~ dia_base:height:spcode)
ca_poshoc

ca_allom <- ca_full1
summary(ca_allom)
anova(ca_allom)

####1c. leaf allometry ####
lf_null1  <- lm(leaf ~ dia_base, all_log)
lf_null2  <- lm(leaf ~ dia_base + height, all_log)
lf_null3  <- lm(leaf ~ height, all_log)
lf_full1  <- lm(leaf ~ dia_base*spcode, all_log)
lf_full2  <- lm(leaf ~ (dia_base + height)*spcode,all_log)
lf_full3  <- lm(leaf ~ height*spcode, all_log)
AIC(lf_null1,lf_null2,lf_null3,lf_full1,lf_full2,lf_full3)

lf_allom  <- lf_full2
anova(lf_allom)
summary(lf_allom)
lf_poshoc <- emmeans(lf_allom, pairwise ~ dia_base:height:spcode)
lf_poshoc
plot(lf_poshoc)

#### 1d.stem allometry ####
sm_null1 <- lm(stem ~ dia_base, all_log)
sm_null2 <- lm(stem ~ height, all_log)
sm_null3 <- lm(stem ~ dia_base + height, all_log)
sm_full1 <- lm(stem ~ dia_base*spcode, all_log)
sm_full2 <- lm(stem ~ height*spcode, all_log)
sm_full3 <- lm(stem ~ (dia_base + height)*spcode, all_log)
AIC(sm_null1,sm_null2,sm_null3,sm_full1,sm_full2,sm_full3)

sm_allom <- sm_full3
summary(sm_allom)
anova(sm_allom)
sm_poshoc <- emmeans(sm_allom,pairwise ~ dia_base:height:spcode)
sm_poshoc
plot(sm_poshoc)

#### 1e. root allometry ####
rt_null1 <- lm(root ~ leaf,all_log)
rt_null2 <- lm(root ~ agb, all_log)
rt_null3 <- lm(root ~ dia_base, all_log)
rt_full1 <- lm(root ~ leaf*spcode,all_log)
rt_full2 <- lm(root ~ agb*spcode,all_log)
rt_full3 <- lm(root ~ dia_base*spcode,all_log)
AIC(rt_null1,rt_null2,rt_null3,rt_full1,rt_full2,rt_full3)

rt_allom <- rt_full1
summary(rt_allom)
anova(rt_allom)
rt_poshoc <- emmeans(rt_allom,pairwise ~ leaf:spcode)
rt_poshoc
plot(rt_poshoc)

####### 2.Hypothesis testing #######

#### 2a. root allocation is greater in perennial #### 

#mature_df      <- all_log %>% filter(!is.na(repro))

frt_hmod    <- lmer(root ~ leaf*growth*photo + (1|spcode)
                 , data = all_log
                 , REML = TRUE)

summary(frt_hmod)
Anova(frt_hmod,type="3",test.statistic = "F") 
#when testing interaction effect, use this setting for Anova




## perennial plants have smaller roots and allocate less to below-ground than annual plants at same leaf biomass level
## C4 plants allocate more to below-ground than C3 plants at same leaf biomass level, but overall there's no difference
## between C3 and C4 plants regarding root biomass.But these statements only apply to mature plants. 
## when all data are used (including seedlings and smaller plants), there's ontogenetic effects that we only
## see higher root allocation at later stage in perennial plants. 



####### H2: fast growlers invest more in reproduction #######

## growth rate

growth <- function(df){
  lm(agb ~ lagday, data=df) 
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
              filter(!is.na(seed_alloc))                             %>%
              left_join(grow_coef,by="spcode")                       %>%
              group_by(spcode)                                       %>% 
              mutate(repro_sd  = sd(seed_alloc)
                    ,repro_se  = repro_sd/sqrt(length(seed_alloc))
                    ,repro_ave = mean(seed_alloc))                   %>% 
              ungroup()                       

                

repro_mod  <- lmer(seed_alloc ~ grate + (1|spcode)
                   ,data=repro_df
                   ,REML=TRUE)
summary(repro_mod)
Anova(repro_mod,type="3",test.statistic = "F") 


# reproduction allocation (as %above-ground biomass)
# is positively related to growth rate of above-ground biomass


######## H3: leaf/stem ratio is influenced by      #######
####### height (-) and specific leaf area (-)      #######  
#zscore <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm = TRUE) 

h3df         <- all_log                               %>% 
                select(spcode,height,carea,sla,lf2sm)        
                
lsratio_mod  <- lmer(lf2sm ~ carea*sla + (1|spcode)
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


## species-specific models 

mod_df <- all_log %>% group_by(spcode) %>% nest()
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

##glance for model statistic
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
rthgt_params <- mod_df %>% dplyr::select(spcode,data,hgt_tidy,
                                       rot_tidy,hgt_pred,rot_pred)          %>% 
  unnest(c(hgt_tidy,rot_tidy),names_sep=".", keep_empty=TRUE)

lsc_params <- mod_df %>% dplyr::select(spcode, data, sm_tidy,lf_tidy,
                                       ca_tidy,sm_pred,lf_pred,ca_pred)     %>% 
  unnest(c(sm_tidy,lf_tidy,ca_tidy),names_sep=".")


rthgt_df <- rthgt_params                               %>% 
            unnest(c(data,hgt_pred,
                     rot_pred),
                  keep_empty = TRUE)

lsc_df  <- lsc_params                                  %>% 
           unnest(c(data,sm_pred,
                    lf_pred,ca_pred),
                  keep_empty = TRUE)


rthgt_intcpt <- rthgt_df %>% filter(hgt_tidy.term=="(Intercept)")
rthgt_slope  <- rthgt_df %>% filter(hgt_tidy.term=="dia_base")

lsc_intcpt <- lsc_df  %>% filter(sm_tidy.term=="(Intercept)")
lsc_slope <-  lsc_df  %>% filter(sm_tidy.term %in% c("dia_base","height"))

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



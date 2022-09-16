## fig-tbl.R

## script to make figures and tables for grass
## allometry

library(ggplot2)
library(RColorBrewer)
library(patchwork)
#require(wesanderson)


source("./ggplot-theme.R")
source("./grass-allometry.R")

RESULTS     = "../results/"

hgt_labels  =  c("Short","Medium","Tall")

all_log     = all_log                                           %>% 
              mutate(hgroup = cut( height
                                   ,breaks = c(-Inf, -1.8, -0.9, Inf)
                                   ,labels = hgt_labels)
                    ,slagrp = cut(sla
                                 ,breaks   = c(-Inf, -3.91, Inf)
                                 ,labels   = c("SLA <= 0.02~m^2~g^{-1}", 
                                               "SLA > 0.02~m^2~g^{-1}")))

allom_fit1     = all_log                                        %>% 
                 select(spcode,short,dia_base,height,hgroup)    %>% 
                 mutate(carea = predict(ca_allom,newdata=.)
                       ,leaf  = predict(lf_allom,newdata=.)
                       ,stem  = predict(sm_allom,newdata=.))
allom_fit2     = all_log                                        %>% 
                 select(spcode,short,dia_base,leaf)             %>% 
                 mutate(height = predict(hgt_allom,newdata=.)
                       ,root   = predict(rt_allom,newdata=.))
                 

hypo_rootfit   = all_log                                        %>% 
                 select(spcode,short,growth,photo,
                       leaf)                                    %>%
                 mutate(root = predict(frt_hmod,newdata=.))

repro_fit     = repro_df                                          %>% 
                select(-seed_alloc)                               %>% 
                mutate(seed_alloc = predict(repro_mod,newdata=.))
         
lfsm_fit      = all_log                                           %>% 
                select(spcode,short,sla,height,carea,slagrp)      %>% 
                mutate(lf2sm = predict(lsratio_mod,newdata=.))


###### Fig1. allometric relationships showing species effects ######

design = "
123
456
" 

fig1a = ggplot(data = all_log, aes(dia_base,carea,colour=short))  +
        geom_point(aes(shape=hgroup),alpha=0.5)                   +
        geom_smooth(data=allom_fit1,method="lm",se=FALSE)         +
        scale_color_brewer(palette="Paired")                      +
        pubtheme.nogridlines                                      +
        theme(legend.position = "none")                           +
        labs(x="",y="ln(Canopy area)")

        
fig1b = ggplot(data = all_log, aes(dia_base,leaf,colour=short))    +
        geom_point(aes(shape=hgroup),alpha=0.5)                    +
        geom_smooth(data=allom_fit1,method="lm",se=FALSE)          +
        scale_color_brewer(palette="Paired")                       +
        pubtheme.nogridlines                                       +
        theme(legend.position = "none")                            +
        labs(x="",y="ln(Leaf biomass)")


fig1c = ggplot(data = all_log, aes(dia_base,stem,colour=short))   +
        geom_point(aes(shape=hgroup),alpha=0.5)                   +
        geom_smooth(data=allom_fit1,method="lm",se=FALSE)         +
        scale_color_brewer(palette="Paired")                      +
        pubtheme.nogridlines                                      +
        theme(legend.position = "bottom")                         +
        labs(x="ln(Basal diameter)",y="ln(Stem biomass)")         +
        guides(colour = guide_legend( ncol=2
              ,label.theme = element_text(face="italic")),
               shape  = guide_legend("Height group"))

fig1d =  ggplot(data = all_log, aes(dia_base,height,colour=short))  +
         geom_point(alpha = 0.5,shape=16)                           +
         geom_smooth(data=allom_fit2,method="lm",se=FALSE)          +
         scale_color_brewer(palette="Paired")                       +
         pubtheme.nogridlines                                       +
         theme(legend.position = "none")                            +
         labs(x="ln(Basal diameter)",y="ln(Plant height)")

fig1e = ggplot(data = all_log, aes(dia_base,root,colour=short))    +
        geom_point(alpha = 0.5,shape=16)                           +
        geom_smooth(data=allom_fit2,method="lm",se=FALSE)          +
        scale_color_brewer(palette="Paired")                       +
        pubtheme.nogridlines                                       +
        theme(legend.position = "none")                            +
        labs(x="ln(Basal diameter)",y="ln(Root biomass)")






fig1 = fig1a + fig1b + fig1c + fig1d + fig1e                      +
       guide_area()                                               +
       plot_layout(design = design, guides = "collect")           +
       plot_annotation(tag_levels = 'a') &
       theme(plot.margin = unit(c(4, 4, 4, 4), "pt"),
             plot.tag.position = c(0, 1),
             plot.tag = element_text(hjust = -2, vjust = 0.4, size = smsize+1),
             legend.title = element_blank())
fig1

ggsave(fig1, file = file.path(RESULTS, "fig1.pdf"), width = col2*1.75 , height= col2, 
       units="px", dpi = 350)

####### Fig2. root allocation between annual perennial C3 C4 #######


fig2  = ggplot(data=all_log,aes(leaf,root,colour=growth))   +
        geom_point(size=1,shape=16,alpha=0.6)               +
        geom_smooth(data=hypo_rootfit,method="lm",se=FALSE) +
        scale_color_manual(values=schwilkcolors[c(1,3)])    +
        pubtheme.nogridlines                                +
        theme(legend.position = "right"
              ,legend.title = element_blank())              +
        labs(x="ln(Leaf biomass)",y="ln(Root biomass)") 
fig2

ggsave(fig2, file = file.path(RESULTS, "fig2.pdf"), width = col1*1.5 , height= 0.85*col1, 
       units="px", dpi = 300) 


####### Fig3. fast-growing plants allocate more to reproduction #######

fig3   = ggplot(data=repro_df,aes(grate,repro_ave))                        +
         geom_errorbar(aes(ymin = repro_ave-repro_se
                          ,ymax = repro_ave+repro_se))                     +
         geom_point(shape=16,size=1)                                       +
         geom_smooth(data=repro_fit,method="lm",se=FALSE,color="black")    +
         pubtheme.nogridlines                                              +
         theme(legend.position = "right"
              ,legend.title = element_blank())                             +
         labs(x="Aboveground biomass growth rate"
             ,y="Reproduction allocation") 
fig3

ggsave(fig3, file = file.path(RESULTS, "fig3.pdf"), width = col1*1.5 , height= 0.85*col1, 
       units="px", dpi = 300) 

####### Fig4. Biomass partition between stem and foliage #######
fig4   = ggplot(data=all_log,aes(height,lf2sm,colour=slagrp))             +
         geom_point(shape=16,size=1,alpha=0.6)                            +
         scale_color_manual(values=schwilkcolors[c(1,3)]
                           ,labels=parse_format())                        +
         geom_smooth(data=lfsm_fit,method="lm",se=FALSE)                  +
         pubtheme.nogridlines                                             +
         theme(legend.position = "right"
              ,legend.title = element_blank())                            +
         labs(x="ln(Plant height)"
             ,y="ln(Leaf to stem ratio)") 
fig4

ggsave(fig4, file = file.path(RESULTS, "fig4.pdf"), width = col1*1.5 , height= 0.85*col1, 
       units="px", dpi = 300)       


####### S1. root allocation in mature plants only #######
maturt_fit   = mature_df                                   %>% 
               select(spcode,short,growth,photo,leaf)      %>% 
               mutate(root = predict(maturt_mod,newdata=.))

s1           = ggplot(data=mature_df,aes(leaf,root,colour=photo))         +
               geom_point(shape=16,size=1,alpha=0.6)                      +
               #facet_wrap(.~growth)                                       +
               scale_colour_manual(values=schwilkcolors[c(1,3)])          +
               geom_smooth(data=maturt_fit,method="lm",se=FALSE)          +
               pubtheme.nogridlines                                       +
               theme(legend.position = "right"
                    ,legend.title = element_blank())                      +
               labs(x="ln(Leaf biomass)"
                   ,y="ln(Root biomass)") 
s1               

ggsave(s1, file = file.path(RESULTS, "s1.pdf"), width = col1*1.5 , height= 0.85*col1, 
       units="px", dpi = 300)    


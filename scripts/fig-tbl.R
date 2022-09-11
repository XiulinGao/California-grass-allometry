## fig-tbl.R

## script to make figures and tables for grass
## allometry

library(ggplot2)
library(RColorBrewer)
library(patchwork)
#require(wesanderson)


source("./ggplot-theme.R")
source("../grass-allometry.R")

RESULTS     = "../results/"

hgt_labels  =  c("Short","Medium","Tall")
all_log     = all_log                                 %>% 
              mutate(hgroup = cut( height
                                   ,breaks = c(-Inf, -1.8, -0.9, Inf)
                                   ,labels = hgt_labels))
allom_fit   = all_log                                          %>% 
              select(spcode,short,date,growth,photo,
                     dia_base,height,leaf)                     %>% 
              mutate(root = predict(frt_hmod,newdata=.))       
              


###### Fig1. allometric relationships showing species effects ######

design = "
123
456
" 
fig1a =  ggplot(data = all_log, aes(dia_base,height,colour=short))  +
         geom_point(alpha = 0.5,shape=16)                           +
         geom_smooth(method="lm",se=FALSE)                          +
         scale_color_brewer(palette="Paired")                       +
         pubtheme.nogridlines                                       +
         theme(legend.position = "none")                            +
         labs(x="",y="Plant height")


        
fig1b = ggplot(data = all_log, aes(dia_base,leaf,colour=short))    +
        geom_point(alpha = 0.5,shape=16)                           +
        geom_smooth(method="lm",se=FALSE)                          +
        scale_color_brewer(palette="Paired")                       +
        pubtheme.nogridlines                                       +
        theme(legend.position = "none")                            +
        labs(x="",y="Leaf biomass")

fig1c = ggplot(data = all_log, aes(dia_base,root,colour=short))    +
        geom_point(alpha = 0.5,shape=16)                           +
        geom_smooth(method="lm",se=FALSE)                          +
        scale_color_brewer(palette="Paired")                       +
        pubtheme.nogridlines                                       +
        theme(legend.position = "none")                            +
        labs(x="Basal diameter",y="Root biomass")


fig1d = ggplot(data = all_log, aes(dia_base,carea,colour=short))  +
        geom_point(aes(shape=hgroup),alpha=0.5)                   +
        geom_smooth(method="lm",se=FALSE)                         +
        scale_color_brewer(palette="Paired")                      +
        pubtheme.nogridlines                                      +
        theme(legend.position = "none")                           +
        labs(x="Basal diameter",y="Canopy area")

fig1e = ggplot(data = all_log, aes(dia_base,stem,colour=short))   +
        geom_point(aes(shape=hgroup),alpha=0.5)                   +
        geom_smooth(method="lm",se=FALSE)                         +
        scale_color_brewer(palette="Paired")                      +
        pubtheme.nogridlines                                      +
        theme(legend.position = "bottom")                         +
        labs(x="Basal diameter",y="Stem biomass")                 +
        guides(colour = guide_legend( ncol=2
                       ,label.theme = element_text(face="italic")),
               shape  = guide_legend("Height group"))

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

####### Fig2. root allocation between annual perennial #######


fig2  = ggplot(data=all_log,aes(leaf,root,colour=growth))   +
        geom_point(shape=16,size=1)                         +
        #facet_wrap(.~photo)                                 +
        geom_smooth(method="lm",se=FALSE)                   +
        scale_color_manual(values=schwilkcolors)            +
        pubtheme.nogridlines                                +
        theme(legend.position = "right"
              ,legend.title = element_blank())              +
        labs(x="Leaf biomass",y="Root biomass") 
fig2

ggsave(fig2, file = file.path(RESULTS, "fig2.pdf"), width = col1*1.5 , height= 0.85*col1, 
       units="px", dpi = 300) 

####### Fig3. fast-growing plants allocate more to reproduction #######

fig3   = ggplot(data=repro_df,aes(grate,repro_ave))          +
         geom_errorbar(aes(ymin = repro_ave-repro_se
                          ,ymax = repro_ave+repro_se))       +
         geom_point(shape=16,size=1)                         +
         geom_smooth(method="lm",se=FALSE,color="black")     +
         pubtheme.nogridlines                                +
         theme(legend.position = "right"
              ,legend.title = element_blank())               +
         labs(x="Aboveground biomass growth rate"
             ,y="Reproduction allocation") 
fig3

ggsave(fig3, file = file.path(RESULTS, "fig3.pdf"), width = col1*1.5 , height= 0.85*col1, 
       units="px", dpi = 300) 

####### Fig4. Biomass partition between stem and foliage #######






## fig-tbl.R

## script to make figures and tables for grass
## allometry study

library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(corrplot)
library(scales)
#library(sjPlot)
#require(wesanderson)

setwd("~/California-grass-allometry/scripts")
source("./ggplot-theme.R")
source("./grass-allometry.R")

RESULTS     = "../results/"

hgt_labels  =  c("Short","Medium","Tall")

all         =   all                                            %>% 
                mutate(hgroup = cut( height
                                   ,breaks = c(-Inf, exp(-1.8), exp(-0.9), Inf)
                                   #, breaks = c(-Inf, exp(-1.92),exp(-1.14), exp(-0.36), Inf)
                                   ,labels = hgt_labels)
                      ,slagrp = cut(sla
                                   ,breaks   = c(-Inf, exp(-3.91), Inf)
                                   ,labels   = c("SLA <= 0.02~m^2~g^{-1}", 
                                               "SLA > 0.02~m^2~g^{-1}")))

all_log     =   all_log                                        %>% 
                mutate(hgroup = cut( height
                      ,breaks = c(-Inf, -1.8, -0.9, Inf)
                      #,breaks = c(-Inf, -1.92, -1.14, -0.36, Inf)
                      ,labels = hgt_labels)
         ,slagrp = cut(sla
                      ,breaks   = c(-Inf, -3.91, Inf)
                      ,labels   = c("SLA <= 0.02~m^2~g^{-1}", 
                                    "SLA > 0.02~m^2~g^{-1}")))

allom_fit1     = all_log                                        %>% 
                 select(spcode,short,dia_base,height,hgroup)    %>% 
                 mutate(carea = predict(ca_allom,newdata=.)
                       ,leaf  = predict(lf_allom,newdata=.)
                       ,stem  = predict(sm_allom,newdata=.))    %>% 
                 mutate(across(c("dia_base","height","carea"
                                ,"leaf","stem"),exp))       
allom_fit2     = all_log                                        %>% 
                 select(spcode,short,dia_base,leaf)             %>% 
                 mutate(height = predict(hgt_allom,newdata=.)
                       ,root   = predict(rt_allom,newdata=.))   %>% 
                 mutate(across(c("dia_base","leaf","height",
                                 "root"),exp))
                 

hypo_rootfit   = all_log                                        %>% 
                 select(spcode,short,growth,photo)              %>%
                 mutate(rt2lf = predict(frt_hmod,newdata=.))    %>% 
                 mutate(across(c("rt2lf"),exp))

repro_fit     = repro_df                                          %>% 
                select(-seed_alloc)                               %>% 
                mutate(seed_alloc = predict(repro_mod,newdata=.)) 
         
lfsm_fit      = all_log                                             %>% 
                select(spcode,short,sla,height,carea,slagrp,hgroup) %>% 
                mutate(lf2sm = predict(lsratio_mod,newdata=.))      %>% 
                mutate(across(c("sla","height","carea",
                                "lf2sm"),exp))


###### Fig1. allometric relationships showing species effects ######

design = "
123
456
" 
mycolor = brewer.pal(n=12, "Paired")[c(1:10,12)]
        
fig1a = ggplot(data = all, aes(dia_base,leaf,colour=short))        +
        geom_point(aes(shape=hgroup),alpha=0.5)                    +
        geom_smooth(data=allom_fit1,method="lm",se=FALSE)          +
        #scale_color_brewer(palette="Paired")                       +
        scale_color_manual(values=mycolor)                         +
        scale_x_continuous(trans="log2"
                          ,breaks=breaks_log(n=6)
                          ,labels=number_format(accuracy=0.1))     +
        scale_y_continuous(trans="log2"
                          ,breaks=breaks_log(n=6)
                          ,labels=label_scientific(digits=1))      +
        pubtheme.nogridlines                                       +
        theme(legend.position = "none")                            +
        labs(x="",y="Leaf biomass (kg)")


fig1b = ggplot(data = all, aes(dia_base,stem,colour=short))       +
        geom_point(aes(shape=hgroup),alpha=0.5)                   +
        geom_smooth(data=allom_fit1,method="lm",se=FALSE)         +
        #scale_color_brewer(palette="Paired")                      +
        scale_color_manual(values=mycolor)                        +
        scale_x_continuous(trans="log2"
                          ,breaks=breaks_log(n=6)
                          ,labels=number_format(accuracy=0.1))    +
        scale_y_continuous(trans="log2"
                          ,breaks=breaks_log(n=6)
                          ,labels=label_scientific(digits=1))     +
        pubtheme.nogridlines                                      +
        theme(legend.position = "bottom")                         +
        labs(x="",y="Stem biomass (kg)")         +
        guides(colour = guide_legend( ncol=2
              ,label.theme = element_text(face="italic",
                                          size=smsize)),
               shape  = guide_legend("Height group"))

fig1c = ggplot(data = all, aes(leaf,root,colour=short))      +
  geom_point(alpha = 0.5,shape=16)                           +
  geom_smooth(data=allom_fit2,method="lm",se=FALSE)          +
  #scale_color_brewer(palette="Paired")                       +
  scale_color_manual(values=mycolor)                         +
  scale_x_continuous(trans="log2"
                     ,breaks=breaks_log(n=6)
                     ,labels=label_scientific(digits=1))     +
  scale_y_continuous(trans="log2"
                     ,breaks=breaks_log(n=6)
                     ,labels=label_scientific(digits=1))     +
  pubtheme.nogridlines                                       +
  theme(legend.position = "none")                            +
  labs(x="Leaf Biomass (kg)",y="Root biomass (kg)")


fig1d = ggplot(data = all, aes(dia_base,carea,colour=short))+
  geom_point(aes(shape=hgroup),alpha=0.5)                   +
  geom_smooth(data=allom_fit1,method="lm",se=FALSE)         +
  #scale_color_brewer(palette="Paired")                      +
  scale_color_manual(values=mycolor)                        +
  scale_x_continuous(trans="log2"
                     ,breaks=breaks_log(n=6)
                     ,labels=number_format(accuracy=0.1))   +
  scale_y_continuous(trans="log2"
                     ,breaks=breaks_log(n=6)
                     ,labels=label_scientific(digits=1))    +
  pubtheme.nogridlines                                      +
  theme(legend.position = "none")                           +
  labs(x="Basal diameter (cm)",y=expression("Canopy area"~"("*m^2*")"))



fig1e =  ggplot(data = all, aes(dia_base,height,colour=short))      +
         geom_point(alpha = 0.5,shape=16)                           +
         geom_smooth(data=allom_fit2,method="lm",se=FALSE)          +
         #scale_color_brewer(palette="Paired")                       +
         scale_color_manual(values=mycolor)                         +
         scale_x_continuous(trans="log2"
                           ,breaks=breaks_log(n=6)
                           ,labels=number_format(accuracy=0.1))     +
         scale_y_continuous(trans="log2"
                           ,breaks=breaks_log(n=5)
                           ,labels=number_format(accuracy=0.1))     +
         pubtheme.nogridlines                                       +
         theme(legend.position = "none")                            +
         labs(x="Basal diameter (cm)",y="Plant height (m)")







fig1 = fig1a + fig1b + fig1c + fig1d + fig1e                      +
       guide_area()                                               +
       plot_layout(design = design, guides = "collect")           +
       plot_annotation(tag_levels = 'a') &
       theme(plot.margin = unit(c(4, 4, 4, 4), "pt"),
             plot.tag.position = c(0, 1),
             plot.tag = element_text(hjust = -2, vjust = 0.4, size = Ltextsize),
             legend.title = element_blank())
fig1

#ggsave(fig1, file = file.path(RESULTS, "fig1_unlog.pdf"), width = col2*1.8 , height= col2, 
#       units="cm", dpi = 350)
ggsave(fig1, file = file.path(RESULTS, "fig1_unlog.tiff"), width = col2 , height= 0.7*col2, 
       units="cm", dpi = 600)



####### Fig2. traits correlation matrix #######

cor_df = obs_predic                %>% 
         select( spcode
               , hgt_slope
               , rot_slope
               , dbh_sm
               , hgt_sm
               , dbh_lf
               , hgt_lf
               , dbh_ca
               , hgt_ca
               , hgt_incpt
               , rot_incpt
               , sm_incpt
               , lf_incpt
               , ca_incpt)         %>% 
        distinct()                 %>% 
        left_join(sla,by="spcode") %>% 
        select(-spcode)

plot_nams = c( "$ H[dbh-exp]", "$ Frt[lf-exp]", "$ Stem[dbh-exp]"
             , "$ Stem[hgt-exp]", "$ Leaf[dbh-exp]", "$ Leaf[hgt-exp]"
             , "$ CanA[dbh-exp]", "$ CanA[hgt-exp]", "$ H[coef]"
             , "$ Frt[coef]", "$ Stem[coef]","$ Leaf[coef]"
             , "$ CanA[coef]", "SLA")       
corrs     = cor(cor_df)
sigtest   = cor.mtest(cor_df, conf.level=0.95, method="spearman")
sig_p     = sigtest$p
colnames(corrs) <- plot_nams
rownames(corrs) <- plot_nams
colnames(sig_p) <- plot_nams
rownames(sig_p) <- plot_nams


f_output  = "fig2.pdf"
COL2(diverging = c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "RdYlBu"), n = 200)
colors <- scales::alpha(colorRampPalette(c("#E57A77","white", "#7ca1cc"))(10), alpha = 1)

#pdf(file = file.path(RESULTS,f_output),width = col1*1.5/300 , height= col1*0.85/300)
pdf(file = file.path(RESULTS,f_output),width = col1 , height= col1*0.7)
par(family=fontfamily)
fig2 = corrplot(corrs
        ,p.mat       = sig_p
        ,sig.level   = 0.05
        ,insig       = "label_sig"
        #,insig       = "p-value"
        ,pch.cex     = 1.5
        ,pch.col     = "black"
        ,method      = "ellipse"
        ,type        = "lower"
        ,diag        =FALSE
        ,col = colors
        #,addCoef.col = "grey50"
        ,mar         = c(0,0,1,0)
        ,order       = "original"
        #,number.cex  = 0.25
        ,tl.cex      = 1.3
        ,tl.col      = "black"
        ,cl.cex      = 1)
dev.off()
 
####### Fig3. root allocation between annual perennial C3 C4 #######


#fig3  = ggplot(data=all,aes(leaf,root))                       +
#        geom_point(size=1,shape=16)                           +
#        geom_smooth(data=hypo_rootfit,method="lm",
#                    color="black",se=FALSE) +
#        geom_abline(slope = 1, intercept=0,color="red",
#                    linetype = "dashed",size=0.8,alpha=0.6)   +
  
        #scale_color_manual(values=schwilkcolors[c(1,4)])    +
#        scale_x_continuous(trans="log2"
#                          ,breaks=breaks_log(n=6)
#                          ,labels=label_scientific(digits=1)) +
#        scale_y_continuous(trans="log2"
#                          ,breaks=breaks_log(n=6)
#                          ,labels=label_scientific(digits=1)) +
#        pubtheme.nogridlines                                  +
#        theme(legend.position = "right"
#              ,legend.title = element_blank())                +
#        labs(x="Leaf biomass (kg)",y="Root biomass (kg)") 

fig3 = ggplot(data=all, aes(photo,rt2lf))                     +
       geom_boxplot(aes(fill=growth))                         +
       geom_point(position=position_dodge(width=0.75),
                  aes(group=growth),size=0.8)                 +
       scale_y_continuous(trans="log2"
                     ,breaks=breaks_log(n=6)
                     ,labels=number_format(accuracy=0.01))    +
       scale_fill_manual(values=c("#D53E4F","#ABD9E9"))       +
       pubtheme.nogridlines                                   +
       theme(legend.position = "right"
            ,legend.title    = element_blank())               +
       labs(x="",y="Root to leaf ratio")
       
fig3

ggsave(fig3, file = file.path(RESULTS, "fig3.tiff"), width = col1 , height= 0.7*col1, 
       units="cm", dpi = 600) 


####### Fig4. fast-growing plants allocate more to reproduction #######

fig4   = ggplot(data=repro_df,aes(grate,repro_ave, color=growth))          +
         geom_errorbar(aes(ymin = repro_ave-repro_se
                          ,ymax = repro_ave+repro_se))                     +
         geom_point(shape=16,size=0.8)                                     +
         geom_smooth(data=repro_fit,method="lm",se=FALSE,color="black")    +
         #scale_color_brewer(palette="Paired")                              +
         scale_color_manual(values = c("#D53E4F","#ABD9E9"))               +
         pubtheme.nogridlines                                              +
         theme(legend.position = "right"
              ,legend.title = element_blank())                             +
         labs(x=expression("Total biomass growth rate"~"("~Ln(kg)~day^-1~")")
             ,y="Reproductive allocation")                                 +
         guides(colour = guide_legend( ncol=1))
                                #,label.theme = element_text(face="italic")))
fig4

ggsave(fig4, file = file.path(RESULTS, "fig4.tiff"), width = col1 , height= 0.7*col1, 
       units="cm", dpi =600) 

####### Fig5. Biomass partition between stem and foliage #######
fig5   = ggplot(data=all,aes(sla,lf2sm,colour=hgroup))                    +
         geom_point(shape=16,size=0.8,alpha=0.6)                          +
         scale_color_manual(values=c("#D53E4F","#A8B6CC","#ABD9E9")
                           ,labels=parse_format())                        +
         geom_smooth(data=lfsm_fit,method="lm",se=FALSE)                  +
         scale_x_continuous(trans="log2"
                     ,breaks=breaks_log(n=4)
                     ,labels=number_format(accuracy=0.001))               +
                      scale_y_continuous(trans="log2"
                     ,breaks=breaks_log(n=6)
                     ,labels=number_format(accuracy=0.1))                 +

         pubtheme.nogridlines                                             +
         theme(legend.position = "bottom"
              ,legend.title = element_blank())                            +
         labs(x= expression("Specific leaf area"~ (m^2~{g^-1}))
             ,y="Leaf to stem ratio") 

fig5

ggsave(fig5, file = file.path(RESULTS, "fig5.tiff"), width = col1, height= 0.7*col1, 
       units="cm", dpi = 600)       


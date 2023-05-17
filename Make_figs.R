rm(list=ls())
source("./ABC_drylands_function.R")

# ---------------------------- Main figures ------------------------------
## >> PCA levels of aggregation data ----


d_biocom=read.table("../Data_new/biocom_data.csv",sep=";")
d_sim=read.table("../Data_new/All_new_sim2.csv",sep=";")%>%
  mutate(., Pooling=recode_factor(Pooling,"1"="Model, no change",
                                  "2" = "Model, x2","3" = "Model, x3","4" = "Model, x4","5" = "Model, x5"))

set.seed(123)
d=rbind(stat_sim[,-c(1:2,15)]%>%dplyr::sample_n(., 400000),
        d_biocom[,c(14:ncol(d_biocom))]%>%
          add_column(.,Pooling='Data'))
sumstat_name=colnames(d)[1:11]
res.comp=imputePCA(d[,which(colnames(d) %in% sumstat_name)],ncp=3,scale = T) 

if ("completeObs" %in% names(res.comp)){
  res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
}else {
  res.pca=PCA(res.comp, ncp = 3,  graph=F)
}

axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

for (i in 1:3){
  assign(paste0("p",i),
         d%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_point(aes(x = PC1, y = PC2, color = Pooling,fill=Pooling,size=Pooling),alpha=.5)+
           scale_size_manual(values=c(rep(.5,5),1.5))+
           scale_color_manual(values=c(my_pal(5),"black"))+
           scale_fill_manual(values=c(my_pal(5),"black"))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="Change in resolution",fill="")+
           ggtitle("")+guides(shape="none")+
           theme_classic()+theme(legend.position = "bottom")+
           guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")
  )
}

p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                      p2+theme(legend.position = "none"),
                      p3+theme(legend.position = "none"),
                      ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.1))

ggsave(paste0("../Figures/Final_figs/PCA_spatial_resolution_model_and_data.pdf"),p,width = 11,height = 5)







## >> Prediction distance tipping point ----


d=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")
keep_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
#example with site 3

site=279
pred=read.table(paste0("../Data_new/Prediction/Dist_tipping_",site,".csv"),sep=",")%>%
  filter(., V1 != 0)
colnames(pred)=c("p","q","cover")

index=0;pred$ID_sim=NA
for (x in 1:nrow(pred)){
  if (pred$p[x]==0.005){
    index=index+1
  }
  pred$ID_sim[x]=index
}

pinfered=pred%>%
  group_by(., ID_sim,q)%>%
  dplyr::summarise(.,.groups = "keep",pinfered=max(p),cover=max(cover))

#example of bifu for a site

p1=ggplot(NULL)+
  geom_line(data=pred%>%filter(., abs(q-median(q)) < quantile(abs(q-median(q)),.6),p>.45)%>%
              add_column(., median_or_not=sapply(1:nrow(.),function(x){
                if (.$q[x]==median(.$q)){return("50")
                } else if (.$q[x]==quantile(.$q,.25)) {
                  return("25")
                } else if (.$q[x]==quantile(.$q,.75)) {
                  return("75")
                } else {
                  return("No")
                }
              })),
            aes(p,y=cover,group=ID_sim,color=median_or_not,size=median_or_not))+
  geom_point(data=pinfered%>%filter(., abs(q-median(q)) < quantile(abs(q-median(q)),.5)),
             aes(x=pinfered,y=cover),color=alpha("red",.3))+
  scale_color_manual(values=(c("No"="gray",'50'="#A564C3","25"="#EABBFF","75"="#EABBFF")))+
  scale_size_manual(values=c("No"=.1,"50"=1.3,"25"=1,"75"=1))+
  xlim(0.46,.62)+
  geom_segment(aes(x = .6, y = .45, xend = .5, yend = .45),arrow = arrow(length = unit(0.3, "cm")),color="black")+
  geom_text(aes(x = .55, y = .47,label="Increasing stress"),color="black")+
  the_theme+
  labs(y="Vegetation cover",x="Reproduction parameter (p)",color="Median of q ? ")+
  guides(color = guide_legend(override.aes = list(size = 1.5)),size="none")

landscape_vege=ggplot(Get_empirical_site(site)%>%melt(.))+
  geom_tile(aes(x=Var1,y=Var2,fill=as.factor(value)))+
  scale_fill_manual(values=rev(c("black","white")))+
  theme_transparent()+
  theme(legend.position = "none")

p1=p1 +
  annotation_custom(grob=ggplotGrob(landscape_vege),
                    ymin = 0, ymax=0.27, xmin=.55, xmax=.64)

# all sites
d_summarized=d%>%
  group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",abs_dis50=quantile(pinfer-pcrit,na.rm = T,.5),
                   abs_dis25=quantile(pinfer-pcrit,na.rm = T,.25),
                   abs_dis75=quantile(pinfer-pcrit,na.rm = T,.75),
                   relativ_dis50=quantile((pinfer-pcrit)/pinfer,na.rm = T,.5),
                   relativ_dis25=quantile((pinfer-pcrit)/pinfer,na.rm = T,.25),
                   relativ_dis75=quantile((pinfer-pcrit)/pinfer,na.rm = T,.75),
                   Size_tipping50=quantile(Size_tipping,na.rm = T,.5),
                   Size_tipping25=quantile(Size_tipping,na.rm = T,.25),
                   Size_tipping75=quantile(Size_tipping,na.rm = T,.75))%>%
  filter(., Site %in% keep_sites)

p2=ggplot(d_summarized%>%melt(., measure.vars=c("Sand","MF","aridity")))+
  geom_pointrange(aes(value,ymin=relativ_dis25,ymax=relativ_dis75,y=relativ_dis50),shape=21,color="black",fill="gray")+
  geom_smooth(aes(value,y=relativ_dis50),fill="#D38DEF",color="#782898",alpha=.3)+
  facet_wrap(.~variable,scales = "free")+
  the_theme+
  labs(y="Relative dist. to desert state",x="Driver value")

p3=ggplot(d_summarized%>%melt(., measure.vars=c("Sand","MF","aridity")))+
  geom_pointrange(aes(value,ymin=Size_tipping25,ymax=Size_tipping75,y=Size_tipping50),shape=21,color="black",fill="gray")+
  facet_wrap(.~variable,scales = "free")+
  the_theme+
  labs(y="Size shift",x="Driver value")


ggsave("../Figures/Final_figs/Dist_tipping.pdf",
       ggarrange(ggarrange(ggplot()+theme_void(),
                           p1,
                           ggplot()+theme_void(),ncol=3,widths = c(.4,1,.4),labels = c("",LETTERS[1],"")),
                 p2,nrow=2,labels=c("",LETTERS[2]),heights = c(1.2,1)),width = 7,height = 7)



## >> Climatic projections ----

keep_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
d=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")
list_clim=list.files("../Data_new/Climatic_data",pattern = ".csv")#[-grep("mean",list.files("../Data_new/Climatic_data",pattern = ".csv"))]

proj_clim=as.data.frame(matrix(NA,length(keep_sites),length(list_clim)))
for (x in 1:ncol(proj_clim)){
  proj_clim[,x]=read.table(paste0("../Data_new/Climatic_data/",list_clim[x]),sep=";")[keep_sites,1]
}
colnames(proj_clim)=c("Aridity RCP 4.5","Aridity RCP 8.5","Temperature RCP 4.5","Temperature RCP 8.5")

# summarizing information in each site
d_summarized=d%>%
  group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",abs_dis50=quantile(pinfer-pcrit,na.rm = T,.5),
                   abs_dis25=quantile(pinfer-pcrit,na.rm = T,.25),
                   abs_dis75=quantile(pinfer-pcrit,na.rm = T,.75),
                   relativ_dis50=quantile((pinfer-pcrit)/pinfer,na.rm = T,.5),
                   relativ_dis25=quantile((pinfer-pcrit)/pinfer,na.rm = T,.25),
                   relativ_dis75=quantile((pinfer-pcrit)/pinfer,na.rm = T,.75),
                   Size_tipping50=quantile(Size_tipping,na.rm = T,.5),
                   Size_tipping25=quantile(Size_tipping,na.rm = T,.25),
                   Size_tipping75=quantile(Size_tipping,na.rm = T,.75))%>%
  filter(., Site %in% keep_sites)%>%
  add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[keep_sites])

d_summarized=cbind(d_summarized,proj_clim)


#first pair correlation between size tipping and distance to desert state


p1=ggplot(d_summarized)+
  geom_pointrange(aes(x=Size_tipping50,ymin=relativ_dis25,ymax=relativ_dis75,y=relativ_dis50,
                      fill=Cover,color=Cover),
                  shape=21)+
  geom_pointrange(aes(x=Size_tipping50,xmin=Size_tipping25,xmax=Size_tipping75,y=relativ_dis50,
                      fill=Cover,color=Cover),
                  shape=21)+
  the_theme+
  scale_fill_gradientn(colours = colorRampPalette(c("#D3EFD3","#90C390","#47A747","#126312"))(100))+
  scale_color_gradientn(colours = colorRampPalette(c("#D3EFD3","#90C390","#47A747","#126312"))(100))+
  labs(y="Distance to desertification",x="Height of the tipping point",fill="Vegetation cover")+
  guides(color="none")+
  theme(legend.position = c(.8, .5),legend.key.size = unit(.5, 'cm'))




#Relative distance to desert state

#defining the colors
mat_tmp=Get_color_classif()
d_summarized$sort_relativ_dist=unlist(sapply(d_summarized$relativ_dis50,function(x){
  return(which(sort(d_summarized$relativ_dis50)==x)[1])
}))
d_summarized$sort_aridity_85=unlist(sapply(d_summarized$`Aridity RCP 8.5`,function(x){
  return(which(sort(d_summarized$`Aridity RCP 8.5`)==x)[1])
}))
d_summarized$color=unlist(sapply(1:nrow(d_summarized),function(x){
  return(mat_tmp[d_summarized$sort_aridity_85[x],d_summarized$sort_relativ_dist[x]])
}))

melting_d=d_summarized%>%
  melt(., measure.vars=colnames(d_summarized)[15:(ncol(d_summarized)-3)])
cols=melting_d$color[1:nrow(d_summarized)];names(cols)=1:nrow(d_summarized)


p2=ggplot(melting_d%>%
            filter(., variable %in% c("Aridity RCP 8.5"))%>%arrange(., relativ_dis50))+
  geom_pointrange(aes(value,ymin=relativ_dis25,ymax=relativ_dis75,y=relativ_dis50,
                      fill=as.character(ID),color=as.character(ID)),
                  shape=21)+
  the_theme+
  theme(legend.position = "none")+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(y="Distance to desertification",x="Vulnerability to climate change")


mat_tmp=Get_color_classif(50)

cols=as.character(mat_tmp);names(cols)=paste0(1:length(cols))

p0_1=ggplot(melt(mat_tmp)%>%add_column(., Val_id=1:nrow(.)))+
  geom_tile(aes(Var1,Var2,fill=as.factor(Val_id)))+
  scale_fill_manual(values = cols)+
  theme_transparent()+
  theme(legend.position = "none",text = element_text(face = "plain"))

p2=p2+annotation_custom(grob=ggplotGrob(p0_1),
                        xmin = -0.0002, xmax=0.00035, ymin=.65, ymax=1.05)+
  geom_richtext(aes(x=0.000075,y=.66,label="Vulnerability"),
                label.size = NA,size=2.5,family = "NewCenturySchoolbook")+
  geom_richtext(aes(x=0.000075,y=.60,label="to climate change"),
                label.size = NA,size=2.5,family = "NewCenturySchoolbook")+
  geom_richtext(aes(x=-0.00018,y=.83,label="Dist. to desert state (b)",angle=90),
                label.size = NA,size=2.5,family = "NewCenturySchoolbook")+
  geom_richtext(aes(x=-0.00028,y=.85,label="Height tip. point (c)",angle=90),
                label.size = NA,size=2.5,family = "NewCenturySchoolbook")

# 
# p0_2=ggplot(d_summarized)+
#   geom_histogram(aes(relativ_dis50),color="black",fill="gray",alpha=.5)+
#   the_theme+theme(legend.position = "none")+
#   labs(x="Distance to desertification",y="Count")
# 
# p1=p1+annotation_custom(grob=ggplotGrob(p0_2),
#                         xmin = 0.00107, xmax=0.00177, ymin=.65, ymax=1.05)+
#   geom_richtext(aes(x=0.0011,y=1,label="b"),
#                 label.size = NA,size=5,fontface="bold")
# 


#Size of the tipping point
#defining the colors
mat_tmp=Get_color_classif()
d_summarized$sort_size_tipping=unlist(sapply(d_summarized$Size_tipping50,function(x){
  return(which(sort(d_summarized$Size_tipping50)==x)[1])
}))
d_summarized$sort_aridity_85=unlist(sapply(d_summarized$`Aridity RCP 8.5`,function(x){
  return(which(sort(d_summarized$`Aridity RCP 8.5`)==x)[1])
}))
d_summarized$color=unlist(sapply(1:nrow(d_summarized),function(x){
  return(mat_tmp[d_summarized$sort_aridity_85[x],d_summarized$sort_size_tipping[x]])
}))
melting_d=d_summarized%>%
  melt(., measure.vars=colnames(d_summarized)[15:(ncol(d_summarized)-3)])
cols=melting_d$color[1:nrow(d_summarized)];names(cols)=1:nrow(d_summarized)


p3=ggplot(melting_d%>%
            filter(., variable %in% c("Aridity RCP 8.5"))%>%arrange(., relativ_dis50))+
  geom_pointrange(aes(value,ymin=Size_tipping25,ymax=Size_tipping75,y=Size_tipping50,
                      fill=as.character(ID),color=as.character(ID)),
                  shape=21)+
  the_theme+
  theme(legend.position = "none")+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(y="Height of the tipping point",x="Vulnerability to climate change")


mat_tmp=Get_color_classif(50)

cols=as.character(mat_tmp);names(cols)=paste0(1:length(cols))

p0_1=ggplot(melt(mat_tmp)%>%add_column(., Val_id=1:nrow(.)))+
  geom_tile(aes(Var1,Var2,fill=as.factor(Val_id)))+
  scale_fill_manual(values = cols)+
  theme_transparent()+
  theme(legend.position = "none",text = element_text(face = "plain"))

# p3=p3+annotation_custom(grob=ggplotGrob(p0_1),
#                         xmin = 0, xmax=0.00055, ymin=.65, ymax=1.05)+
#   geom_richtext(aes(x=0.000275,y=.65,label="Vulnerability to \n climate change"),
#                 label.size = NA,size=.65,family = "NewCenturySchoolbook")+
#   geom_richtext(aes(x=-0.00001,y=.85,label="Height tipping point",angle=90),
#                 label.size = NA,size=3,family = "NewCenturySchoolbook")
# 
# p0_2=ggplot(d_summarized)+
#   geom_histogram(aes(Size_tipping50),color="black",fill="gray",alpha=.5)+
#   the_theme+theme(legend.position = "none")+
#   labs(x="Height of the tipping point",y="Count")
# 
# p2=p2+annotation_custom(grob=ggplotGrob(p0_2),
#                         xmin = 0.00107, xmax=0.00177, ymin=.3, ymax=.6)+
#   geom_richtext(aes(x=0.0011,y=1,label="b"),
#                 label.size = NA,size=5,fontface="bold")



p_tot=ggarrange(p1,p2,p3,nrow=3,labels=letters[1:3])
ggsave("../Figures/Final_figs/Predicting_stability.pdf",p_tot,width = 5,height = 9)



## >> Comparing predicted resilience metrics with parameters and site cover ----


keep_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
d_biocom=read.table("../Data_new/biocom_data.csv",sep=";")

d=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")
d2=read.table("../Data_new/posterior_param.csv",sep=";",header=T)
d2=tibble(Site=1:345,mean_p=apply(d2[,1:345],2,mean),sd_p=apply(d2[,1:345],2,sd),
          mean_q=apply(d2[,346:690],2,mean),sd_q=apply(d2[,346:690],2,sd),
          Plot_n=d_biocom$Plot_n)%>%
  filter(., Site %in% keep_sites)

d_summarized=d%>%
  group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",abs_dis50=quantile(pinfer-pcrit,na.rm = T,.5),
                   abs_dis25=quantile(pinfer-pcrit,na.rm = T,.25),
                   abs_dis75=quantile(pinfer-pcrit,na.rm = T,.75),
                   relativ_dis50=quantile((pinfer-pcrit)/pinfer,na.rm = T,.5),
                   relativ_dis25=quantile((pinfer-pcrit)/pinfer,na.rm = T,.25),
                   relativ_dis75=quantile((pinfer-pcrit)/pinfer,na.rm = T,.75),
                   Size_tipping50=quantile(Size_tipping,na.rm = T,.5),
                   Size_tipping25=quantile(Size_tipping,na.rm = T,.25),
                   Size_tipping75=quantile(Size_tipping,na.rm = T,.75))%>%
  filter(., Site %in% keep_sites)%>%
  add_column(.,Cover=d_biocom$Cover[keep_sites],
             p=d2$mean_p,q=d2$mean_q)



p1=ggplot(d_summarized%>%melt(., measure.vars=c("Cover","p","q"))%>%
            mutate(., variable=recode_factor(variable,"Cover"="Vegetation cover",
                                             "p"="Parameter p","q"="Parameter q")))+
  geom_pointrange(aes(value,abs_dis50,
                      ymin=abs_dis25,
                      ymax=abs_dis75),fill="#B471D8",
                  color="#6A3388",alpha=.5,shape=21,size=.7)+
  facet_wrap(.~variable,ncol=3,scales="free")+
  the_theme+
  labs(x="Value",y="Predicted abslute distance \n to desert stat")+
  theme(strip.text.x = element_text(size=13,face="italic"))

p2=ggplot(d_summarized%>%melt(., measure.vars=c("Cover","p","q"))%>%
            mutate(., variable=recode_factor(variable,"Cover"="Vegetation cover",
                                             "p"="Parameter p","q"="Parameter q")))+
  geom_pointrange(aes(value,Size_tipping50,
                      ymin=Size_tipping25,
                      ymax=Size_tipping75),fill="#A6D67E",
                  color="#6F964F",alpha=.5,shape=21,size=.7)+
  facet_wrap(.~variable,ncol=3,scales="free")+
  the_theme+
  labs(x="Value",y="Predicted tipping size")+
  theme(strip.text.x = element_blank())

p3=ggplot(d_summarized%>%melt(., measure.vars=c("Cover","p","q"))%>%
            mutate(., variable=recode_factor(variable,"Cover"="Vegetation cover",
                                             "p"="Parameter p","q"="Parameter q")))+
  geom_pointrange(aes(value,relativ_dis50,
                      ymin=relativ_dis25,
                      ymax=relativ_dis75),fill="#7095D2",
                  color="#7095D2",alpha=.5,shape=21,size=.7)+
  facet_wrap(.~variable,ncol=3,scales="free")+
  the_theme+
  labs(x="Value",y="Predicted relative distance \n to desert stat")+
  theme(strip.text.x = element_blank())


ggsave("../Figures/Final_figs/Cover_param_association_resilience_metrics.pdf",
       ggarrange(p1,p3,p2,nrow=3,labels=letters[1:3],heights = c(1.2,1,1)),
       width=9,height = 9)



## >> Maping vulnerability ----

melting_d=melting_d%>%
  filter(., variable=="Aridity RCP 8.5")

melting_d$long=d_biocom$Longitude[melting_d$Site]
melting_d$lat=d_biocom$Lattitude[melting_d$Site]


coord_rect=tibble(xmin=c(-10,-85,-125),xmax=c(15,-35,-107),
                  ymin=c(25,-45,30),ymax=c(50,0,45),label=c("(1)","(2)","(3)"))

world_map <- map_data("world")
p2=ggplot(NULL) +
  geom_polygon(data=world_map, aes(x = long, y = lat, group = group),
               fill="lightgray", colour = "white")+
  geom_jitter(data=melting_d,
              aes(x=long,y=lat,color=as.character(ID)),size=1,width = .25,height = .25)+
  geom_rect(data=coord_rect,
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="transparent",color="black")+
  geom_text(data=coord_rect,aes(x=xmax+6,y=ymax+6,label=label))+
  the_theme+
  scale_color_manual(values=cols)+
  labs(x="",y="")+
  theme(axis.title = element_blank())+
  theme(legend.position = "none",axis.line = element_blank())

p3_1=ggplot(NULL) +
  geom_polygon(data=world_map%>%filter(., long<coord_rect$xmax[1] & long>coord_rect$xmin[1],
                                       lat>coord_rect$ymin[1] & lat<coord_rect$ymax[1]), aes(x = long, y = lat, group = group),
               fill="lightgray", colour = "white")+
  geom_jitter(data=melting_d%>%filter(., long<coord_rect$xmax[1] & long>coord_rect$xmin[1],
                                      lat>coord_rect$ymin[1] & lat<coord_rect$ymax[1])%>%arrange(., relativ_dis50),
              aes(x=long,y=lat,color=as.character(ID),size=relativ_dis50),width = .25,height = .25)+
  scale_color_manual(values=cols)+scale_size_continuous(range=c(.3,.7))+
  labs(x="",y="")+
  theme_transparent()+
  theme(legend.position ="none" )

p3_2=ggplot(NULL) +
  geom_polygon(data=world_map%>%filter(., long<coord_rect$xmax[2] & long>coord_rect$xmin[2],
                                       lat>coord_rect$ymin[2] & lat<coord_rect$ymax[2]), aes(x = long, y = lat, group = group),
               fill="lightgray", colour = "white")+
  geom_jitter(data=melting_d%>%filter(., long<coord_rect$xmax[2] & long>coord_rect$xmin[2],
                                      lat>coord_rect$ymin[2] & lat<coord_rect$ymax[2])%>%arrange(., relativ_dis50),
              aes(x=long,y=lat,color=as.character(ID),size=relativ_dis50),width = .25,height = .25)+
  scale_color_manual(values=cols)+scale_size_continuous(range=c(.3,.7))+
  labs(x="",y="")+
  theme_transparent()+
  theme(legend.position ="none" )


p3_3=ggplot(NULL) +
  geom_polygon(data=world_map%>%filter(., long<coord_rect$xmax[3]+10 & long>coord_rect$xmin[3],
                                       lat>coord_rect$ymin[3] & lat<coord_rect$ymax[3]+5), aes(x = long, y = lat, group = group),
               fill="lightgray", colour = "white")+
  geom_jitter(data=melting_d%>%filter(., long<coord_rect$xmax[3] & long>coord_rect$xmin[3],
                                      lat>coord_rect$ymin[3] & lat<coord_rect$ymax[3])%>%arrange(., relativ_dis50),
              aes(x=long,y=lat,color=as.character(ID),size=relativ_dis50),width = .25,height = .25)+
  scale_color_manual(values=cols)+scale_size_continuous(range=c(.3,.7))+
  labs(x="",y="")+
  theme_transparent()+
  theme(legend.position ="none" )


p_tot=ggarrange(p1,
                ggarrange(p2,
                          ggarrange(p3_1,p3_2,p3_3,ncol=3,labels = coord_rect$label),nrow=2,heights = c(1, .8),labels=c("c","")),
                nrow=2,heights = c(1,1.3),labels = c("a",""))

ggsave("../Figures/Final_figs/Climatic_trend.pdf",p_tot,width = 7,height = 10)



# ---------------------------- SI figures ------------------------------

# >> 1) Characteristics of empirical data ----
## Resolution, geographical repartition

d_biocom=read.table("../Data_new/biocom_data.csv",sep=";")

# spatial resolution

p=ggplot(d_biocom)+
  geom_histogram(aes(x=Nbpixels),fill=alpha("blue",.5))+
  the_theme+
  labs(x='# of pixels',y="Count")

ggsave("../Figures/Final_figs/SI/Spatial_resolution_data.pdf",p,width = 6,height = 3)


# Distribution of empirical data, map, aridity and sand cover
d_biocom=read.table("../Data_new/biocom_data.csv",sep=";")

world_map <- map_data("world")
p=ggplot(NULL) +
  geom_polygon(data=world_map, aes(x = long, y = lat, group = group),
               fill="lightgray", colour = "white")+
  geom_point(data=d_biocom,aes(x=Longitude,y=Lattitude,color=Aridity),size=3)+
  the_theme+
  scale_color_gradientn(colors = my_pal(4))+
  labs(x="Longitude",y="Lattitue")

ggsave("../Figures/Final_figs/SI/Map_empirical_sites.pdf",p,width = 6,height = 4)




## Density of summary statistics: ecosystem type, type patterns

#Pair correlation metrics used

p=ggpairs(d_biocom%>%add_column(., for_color=1)%>%dplyr::select(., -Cover)%>%
            dplyr::rename(., SDR=Spectral_ratio,Clustering=clustering,Cover=rho_p,
                          "# neigh"=nb_neigh,Skewness=skewness,Variance=variance,"Autocorr."=moran_I,
                          "Exponent p.l."=PL_expo,"CV PSD"=cv_psd,Fmax=fmax_psd),
          columns = c(13:23),
          mapping = ggplot2::aes(color = as.factor(for_color),size=as.factor(for_color)),
          upper = "blank",
          diag = NULL)+
  scale_color_manual(values=c("#96C3DC"))+
  scale_size_manual(values=.4)+
  the_theme

ggsave(paste0("../Figures/Final_figs/SI/Pair_corr_metrics.pdf"),p,width = 12,height = 12)




#coupling with the one made biogeo consideration

classif_biogeo=read.table("../Data_new/Veg_type_biogeo.csv",sep=";")
d_plant_type=d_biocom%>%add_column(., Type_vege=classif_biogeo$type)

p=ggplot(d_plant_type%>%melt(., measure.vars=colnames(d_plant_type)[14:24])%>%
           mutate(., variable=recode_factor(variable,
                                            "rho_p"="Cover","nb_neigh"="# neighbors","clustering"= "Clustering",
                                            "skewness"="Skewness","variance"="Variance",
                                            "moran_I"="Autocorrelation","Spectral_ratio"="SDR",
                                            "cv_psd"="CV PSD","fmax_psd"="Frac. max","PL_expo"="Exponent p.l." ,
           )))+
  geom_density(aes(x=value,fill=Type_vege),alpha=.5)+
  the_theme+
  facet_wrap(.~variable,scales = "free")+
  labs(x="Value",y="Density",fill="")+
  scale_fill_manual(values=c("#85AB61","#B55960","#ECC570"))

ggsave(paste0("../Figures/Final_figs/SI/Density_type_vege.pdf"),p,width = 8,height = 6)




# Type of vegetation patterns: regular versus irregular

p=ggplot(d_biocom%>%
           melt(., measure.vars=colnames(d_biocom)[14:ncol(d_biocom)])%>%
           mutate(., variable=recode_factor(variable,
                                            "rho_p"="Cover","nb_neigh"="# neighbors","clustering"= "Clustering",
                                            "skewness"="Skewness","variance"="Variance",
                                            "moran_I"="Autocorrelation","Spectral_ratio"="SDR",
                                            "cv_psd"="CV PSD","fmax_psd"="Frac. max","PL_expo"="Exponent p.l." ,
           ),Regular_berdugo=recode_factor(Regular_berdugo,"0"="Irregular","1"="Regular")))+
  geom_density(aes(x=value,fill=as.factor(Regular_berdugo)),alpha=.8)+
  the_theme+
  labs(x="Value",y="Density",fill="Type of pattern  ")+
  facet_wrap(.~variable,scales = "free",nrow=3)+
  scale_fill_manual(values=c("#CCB6EA","#BADCA1"))

ggsave("../Figures/Final_figs/SI/Density_empirical_type_patterns.pdf",p+theme(strip.background.x = element_blank()),width = 9,height = 6)







# >> 2) Optimizing ABC ----
## Optimization of the ABC method: pre- and post-processing


#Simple rejection algorithm
d=read.table("../Data_new/NRMSE/RMSE_param_BoxCox_rejection_optim_lambda_yes_N1_1000.csv",sep=";")

mean_rmse_rej=d%>%
  melt(.)%>%
  group_by(.,variable)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value,na.rm=T))%>%
  dplyr::rename(., Parameter=variable)%>%
  add_column(., Method="BoxCox & Rejection")

p2=ggplot(d%>%melt(.)%>%add_column(., Method="BoxCox & Rejection")%>%
           dplyr::rename(., "Parameter"="variable"))+
  geom_jitter(aes(x=Method,y=value),color="gray",alpha=.5,width =.05,height=0)+
  geom_point(data=mean_rmse_rej,aes(x=Method,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="",y="NRMSE",color="")+
  facet_grid(.~Parameter)+
  the_theme+
  theme(strip.text.x = element_text(size=10),axis.text.x = element_text(angle=60,hjust=1))+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  theme(legend.position = "none",strip.background.x = element_blank())




#With post sampling adjustment

all_sim=expand.grid(N1=c(1000,3000),
                    lambda=c("yes"),
                    Preproc=c("BoxCox","None"),
                    postproc=c("loclinear","neuralnet"))

d=tibble()
for (i in 1:nrow(all_sim)){
  d=rbind(d,read.table(paste0("../Data_new/NRMSE/RMSE_param_",all_sim$Preproc[i],"_",all_sim$postproc[i],"_optim_lambda_",
                              all_sim$lambda[i],"_N1_",all_sim$N1[i],".csv"),sep=";")%>%
            add_column(., N1=all_sim$N1[i],optim_lambda=all_sim$lambda[i],Post=all_sim$postproc[i],Pre=all_sim$Preproc[i]))
}


mean_rmse=d%>%
  melt(., id.vars=c("N1","optim_lambda","Post","Pre"))%>%
  mutate(., Post=recode_factor(Post,"loclinear"="Linear regression","neuralnet"="Non-linear regression"))%>%
  mutate(., Pre=recode_factor(Pre,"None"="No BoxCox","BoxCox"="Box-Cox"))%>%
  add_column(., Treatment=paste0(.$Pre," & \n ",.$Post))%>%
  group_by(.,variable,N1,optim_lambda,Post,Pre,Treatment)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value,na.rm=T))%>%
  dplyr::rename(., Parameter=variable)

p1=ggplot(d%>%melt(., id.vars=c("N1","optim_lambda","Post","Pre"))%>%
            dplyr::rename(., "Parameter"="variable")%>%
            mutate(., Post=recode_factor(Post,"loclinear"="Linear regression","neuralnet"="Non-linear regression"))%>%
            mutate(., Pre=recode_factor(Pre,"None"="No BoxCox","BoxCox"="Box-Cox"))%>%
            add_column(., Treatment=paste0(.$Pre," & \n ",.$Post)))+
  geom_jitter(aes(x=Treatment,y=value,color=interaction(Pre)),
              position = position_jitterdodge(jitter.width = 0.3,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=Treatment,y=mean_rmse,shape=Pre),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="",y="NRMSE",color="")+
  facet_grid(N1~Parameter,labeller = label_bquote(cols="Parameter"==.(as.character(Parameter)),rows=N[1]==.(N1)))+
  the_theme+
  geom_hline(data=tibble(Parameter=c("p","q"),hpos=mean_rmse_rej$mean_rmse[1:2]),aes(yintercept = hpos),color="gray50")+
  theme(strip.text.x = element_text(size=10),axis.text.x = element_text(angle=60,hjust=1))+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values=c("#C46FC5","#80BD5C"))+
  theme(legend.position = "none",axis.title.y = element_blank())

p=ggarrange(ggarrange(ggplot()+theme_void(),p2,ggplot()+theme_void(),nrow=3,heights = c(1,3,1)),p1,ncol=2,widths = c(1,3),labels = LETTERS[1:2])
ggsave(paste0("../Figures/Final_figs/SI/Optimization_inference_preprocessing.pdf"),p,width = 8,height = 6)


## Optimization of the ABC method: PLS versus no-PLS




d=rbind(read.table(paste0("../Data_new/NRMSE/RMSE_hidden_preprocessing_PLS_10_Nnet_10.csv"),sep=";")%>%
          add_column(., PLS="Yes"),
        read.table(paste0("../Data_new/NRMSE/RMSE_hidden_preprocessing_NoPLS_10_Nnet_10.csv"),sep=";")%>%
          add_column(., PLS="No"))



mean_rmse=d%>%
  melt(., id.vars=c("PLS"))%>%
  group_by(.,variable,PLS)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  dplyr::rename(., "Parameter"="variable")


p=ggplot(d%>%melt(., id.vars=c("PLS"))%>%
           dplyr::rename(., "Parameter"="variable"))+
  geom_jitter(aes(x=PLS,y=value,color=as.factor(PLS)),
              position = position_jitterdodge(jitter.width = 0.3,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=PLS,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="Using PLS during pre-processing",y="NRMSE",color="")+
  facet_grid(.~Parameter,labeller = label_bquote(cols="Parameter"==.(as.character(Parameter))))+
  the_theme+
  
  theme(strip.text.x = element_text(size=10),legend.position = "none")+
  scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))

ggsave(paste0("../Figures/Final_figs/SI/Optimization_PLS.pdf"),p,width = 6,height = 3)








## Optimization of the ABC method: neural-network


all_sim=expand.grid(rep_network=seq(10,30,by=10),N_hidden=seq(5,25,by=5))
d=tibble()

for (i in 1:nrow(all_sim)){
  d=rbind(d,read.table(paste0("../Data_new/NRMSE/RMSE_hidden_preprocessing_NoPLS_",
                              all_sim$N_hidden[i],"_Nnet_",all_sim$rep_network[i],".csv"),sep=";")%>%
            add_column(., N_hidden=all_sim$N_hidden[i],N_rep_net=all_sim$rep_network[i]))
}


mean_rmse=d%>%
  melt(., id.vars=c("N_hidden","N_rep_net"))%>%
  group_by(.,variable,N_rep_net,N_hidden)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  dplyr::rename(., "Parameter"="variable")%>%
  mutate(., N_hidden=as.character(N_hidden))


p=ggplot(d%>%melt(., id.vars=c("N_hidden","N_rep_net"))%>%
           dplyr::rename(., "Parameter"="variable")%>%
           mutate(., N_hidden=as.character(N_hidden)))+
  geom_jitter(aes(x=factor(N_hidden,level=c("5","10",'15',"20",'25')),y=value,color=as.factor(N_hidden)),
              position = position_jitterdodge(jitter.width = 0.3,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=N_hidden,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="Number hidden neurons",y="NRMSE",color="")+
  facet_grid(Parameter~N_rep_net,labeller = label_bquote(rows="Parameter"==.(as.character(Parameter)),cols="# evaluation NN"==.(N_rep_net)))+
  the_theme+
  theme(strip.text.x = element_text(size=10),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"),breaks=c('5', '10', '15',"20","25"))

ggsave(paste0("../Figures/Final_figs/SI/Optimization_NN.pdf"),
       p,width = 7,height = 4)






## Number of simulations kept

d=tibble()
list_f=list.files("../Data_new/NRMSE/","NA")[-grep(pattern = "rej",list.files("../Data_new/NRMSE/","NA"))]
for (k in list_f){
  d=rbind(d,read.table(paste0("../Data_new/NRMSE/",k),sep=";")%>%
            add_column(., Nkept=as.numeric(gsub(".csv","",strsplit(k,"_")[[1]][3]))))
}



mean_rmse=d%>%
  melt(., id.vars=c("Nkept"))%>%
  group_by(.,variable,Nkept)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value,na.rm=T))%>%
  dplyr::rename(., Parameter=variable)

p1=ggplot(d%>%melt(., id.vars=c("Nkept"))%>%dplyr::rename(., Parameter=variable))+
  geom_jitter(aes(x=Nkept,y=value,color=Parameter),
              position = position_jitterdodge(jitter.width = 10,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=Nkept,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="Number of simulations kept",y="NRMSE",color="")+
  facet_grid(.~Parameter,labeller = label_bquote(cols="Parameter"==.(as.character(Parameter))))+
  the_theme+
  
  theme(strip.text.x = element_text(size=10),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values=c("#C46FC5","#80BD5C"))+
  theme(legend.position = "none")




d=tibble()
list_f=list.files("../Data_new/NRMSE/","NA")[grep(pattern = "rej",list.files("../Data_new/NRMSE/","NA"))]
for (k in list_f){
  d=rbind(d,read.table(paste0("../Data_new/NRMSE/",k),sep=";")%>%
            add_column(., Nkept=as.numeric(gsub(".csv","",strsplit(k,"_")[[1]][4]))))
}



mean_rmse=d%>%
  melt(., id.vars=c("Nkept"))%>%
  group_by(.,variable,Nkept)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value,na.rm=T))%>%
  dplyr::rename(., Parameter=variable)

p2=ggplot(d%>%melt(., id.vars=c("Nkept"))%>%dplyr::rename(., Parameter=variable))+
  geom_jitter(aes(x=Nkept,y=value,color=Parameter),
              position = position_jitterdodge(jitter.width = 10,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=Nkept,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="Number of simulations kept",y="NRMSE",color="")+
  facet_grid(.~Parameter,labeller = label_bquote(cols="Parameter"==.(as.character(Parameter))))+
  the_theme+
  
  theme(strip.text.x = element_text(size=10),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values=c("#C46FC5","#80BD5C"))+
  theme(legend.position = "none")

ggsave(paste0("../Figures/Final_figs/SI/N_sim_kept.pdf"),
       ggarrange(p1,p2,nrow=2,labels = c("A, non-linear regression","B, no post-sampling adjustment"),vjust=c(4,4),hjust=c(-.25,-.22))
       ,width = 7,height = 8)


## Best summary statistics


d=tibble()
list_f=list.files("../Data_new/Best_sumstat")
all_name=c("All","No PLR","No Exponent p.l.","No PLR & \n Exponent p.l.","No CV PSD","No Frac. max",
           "No CV PSD & \n Frac. max","No CV PSD & \n Exponent p.l. ","No CV PSD, PLR & \n Exponent p.l.")

for (i in 1:length(list_f)){
  d=rbind(d,read.table(paste0("../Data_new/Best_sumstat/",list_f[i]),sep=";")%>%
            add_column(.,Name=all_name[i]))
}

mean_rmse=d%>%
  melt(., id.vars=c("Name"))%>%
  group_by(.,variable,Name)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value,na.rm=T))%>%
  dplyr::rename(., Parameter=variable)

p=ggplot(d%>%
           melt(., id.vars=c("Name")))+
  geom_jitter(aes(x=Name,y=value,color=interaction(Name)),
              position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=Name,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="",y="NRMSE",color="")+
  facet_wrap(.~Parameter,labeller = label_bquote(cols="Parameter"==.(as.character(Parameter))))+
  the_theme+
  
  theme(strip.text.x = element_text(size=10),axis.text.x = element_text(angle=60,hjust=1),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  #scale_color_manual(values=my_pal(9))+
  scale_color_manual(values=colorRampPalette(colors=c("#C46FC5","#80BD5C"))(9))+
  theme(legend.position = "none")+
  ylim(0,.2)

ggsave(paste0("../Figures/Final_figs/SI/Combination_sumstats.pdf"),p,width = 9,height = 6)





d=tibble()
list_f=list.files("../Data_new/Inferrence/","NRMSE_param_rej")
all_name=c("No CV PSD","No Exponent p.l.","No PLR & \n Exponent p.l.","No PLR",
           "No CV PSD, PLR & \n Exponent p.l.","No CV PSD, PLR \n Exponent p.l. & PLR")

d_all_sumstat=read.table(paste0("../Data_new/Inferrence/NRMSE_param_rej_all.csv"),sep=";")
for (i in 1:(length(list_f))){
  post=read.table(paste0("../Data_new/Inferrence/",list_f[i+1]),sep=";")
  d=rbind(d,data.frame(Name=all_name[i],p=apply(post[1:345],2,median)-apply(d_all_sumstat[1:345],2,median),
                       q=apply(post[346:690],2,median)-apply(d_all_sumstat[346:690],2,median),Site=1:345))
}


p=ggplot(d%>%
           melt(., id.vars=c("Name","Site"))%>%
           dplyr::rename(., Parameter=variable))+
  geom_line(aes(x=Name,y=value,group=Site),
             color="gray",lwd=.3,alpha=.4)+
  geom_violin(aes(x=Name,y=value,color=interaction(Name)),width=3)+
  labs(x="",y="Change in parameter values compared \n with all spatial statistics",color="")+
  facet_wrap(.~Parameter,labeller = label_bquote(cols="Parameter"==.(as.character(Parameter))))+
  the_theme+
  theme(strip.text.x = element_text(size=10),axis.text.x = element_text(angle=60,hjust=1),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values=colorRampPalette(colors=c("#C46FC5","#80BD5C"))(7))+
  theme(legend.position = "none")


ggsave(paste0("../Figures/Final_figs/SI/Combination_sumstats_on_p_q.pdf"),p,width = 10,height = 5)


# >> 3) Spatial resolution & system size ---- 
## Change spatial stats with resolution

stat_sim=read.table("../Data_new/All_new_sim2.csv",sep=";")%>%
  mutate(., Id_sim=rep(1:(nrow(.)/5),each=5))%>%
  mutate(., Pooling=as.character(Pooling),rho_p=round(rho_p,5),variable=as.character(variable))%>%
  melt(., id.vars=c("Pooling","Id_sim"))%>%
  mutate(., variable=recode_factor(variable,
                                   "rho_p"="Cover","nb_neigh"="# neighbors","clustering"= "Clustering","skewness"="Skewness","variance"="Variance",
                                   "moran_I"="Autocorrelation","Spectral_ratio"="SDR"
  ))

set.seed(123)
p=ggplot(NULL)+
  geom_line(data=stat_sim%>%filter(., Id_sim %in% sample(unique(.$Id_sim),50)),
            aes(x=Pooling,y=value,group=Id_sim),color="gray50",alpha=.3,lwd=.3)+
  geom_line(data=stat_sim%>%
           group_by(., Pooling,variable)%>%
             dplyr::summarise(., .groups = "keep",mean_value=mean(value,na.rm = T)),
           aes(x=Pooling,y=mean_value,group=interaction(variable)),lwd=1,color="red")+
  geom_point(data=stat_sim%>%
               group_by(., Pooling,variable)%>%
               dplyr::summarise(., .groups = "keep",mean_value=mean(value,na.rm = T)),
             aes(x=Pooling,y=mean_value),color="red",fill="white",shape=21,size=2.5)+
  facet_wrap(.~variable,scales = "free",nrow = 2)+
  labs(x="Change in resolution",y="Mean value across all simulations")+
  scale_x_discrete(labels=c("No change","x2","x3","x4","x5"))+
  the_theme+theme(axis.text.x = element_text(hjust=1,angle=60))+
  guides(color = guide_legend(override.aes = list(size = 2)),fill="none")

ggsave(paste0("../Figures/Final_figs/SI/Change_metrics_spatial_resolution_model.pdf"),p,width = 8,height = 5)

## Robustness inference with spatial resolution

d_RMSE_param=read.table("../Data_new/Scale_obs_indentifiability/Retrieving_parameters_different_resolution_RMSE_param.csv",sep=";")
p1=ggplot(d_RMSE_param%>% #we remove the scale of observation
            mutate(., Scale_obs=as.character(Scale_obs))%>%
            filter(., Method=="NeuralNet")%>%
            melt(., id.vars=c("Site_ID","Method","Scale_obs"))%>%
            filter(., variable=="Pooling")%>%
            mutate(., variable=recode_factor(variable, "Pooling"="Scale obs.")))+
  geom_jitter(aes(x=Scale_obs,y=value,color=Scale_obs),alpha=.5,size=1,position = position_jitter(height = 0,width = .1))+
  labs(x="Change in spatial resolution",y="NRMSE",color="")+
  the_theme+
  guides(Method="none")+
  scale_color_manual(values=my_pal(5))+
  scale_x_discrete(labels = c("No change","x2","x3","x4","x5"))+
  guides(color="none")

ggsave("../Figures/Final_figs/SI/NMRSE_consistency_inference_param_scale.pdf",p1,width = 6,height = 3)


x_y_param=read.table("../Data_new/Scale_obs_indentifiability/Retrieving_parameters_different_resolution_x_y.csv",sep=";")
p2=ggplot(x_y_param%>% #we remove the scale of observation
            filter(., Method=="NeuralNet")%>%
            melt(., id.vars=c("Site_ID","Method","Scale_obs","Type"))%>%
            mutate(., variable=recode_factor(variable, "Pooling"="Scale obs."))%>%
            filter(., Type=="Sim",variable != "Scale obs."))+
  geom_line(aes(x=Scale_obs,y=value,group=Site_ID),alpha=.5,lwd=.5,color="gray")+
  facet_wrap(.~variable,labeller = label_bquote(cols = Parameter==.(as.character(variable))))+
  labs(x=substitute(paste("Scale of observation (",eta,")")),y="Median of the posterior distribution",color="")+
  the_theme+
  guides(Method="none")

ggsave("../Figures/Final_figs/SI/Consistency_inference_param_scale.pdf",p2,width = 7,height = 4)

## System size and summary statistics


list_f=list.files("../Data_new/System_size")
d=tibble()
for (k in list_f){ #to simplify, we average the replicates
  d2=read.table(paste0("../Data_new/System_size/",k),sep=",")%>%
    filter(., V3>0)
  colnames(d2)=c("p","q","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                 "Spectral_ratio","PLR","PL_expo","cv_psd","fmax_psd")
  
  if (length(unique(d2$p,6))==nrow(d2)/60){
    for (z in 1:(nrow(d2)/15)){
      d=rbind(d,as_tibble(t(colMeans(d2[(1+(z-1)*15):(z*15),1:13],na.rm = T))))
    }
  }
}
colnames(d)=c("p","q","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
              "Spectral_ratio","PLR","PL_expo","cv_psd","fmax_psd")

d$System_size=rep(seq(75,225,50),nrow(d)/4)
d$Id_sim=rep(1:(nrow(d)/4),each=4)
stat_sim=d%>%
  melt(., id.vars=c("System_size","Id_sim"))%>%
  filter(., variable %in% c("rho_p","nb_neigh","clustering","skewness","variance","moran_I","Spectral_ratio","cv_psd","fmax_psd"))%>%
  mutate(., variable=recode_factor(variable,
                                   "rho_p"="Cover","nb_neigh"="# neighbors","clustering"= "Clustering","skewness"="Skewness","variance"="Variance",
                                   "moran_I"="Autocorrelation","Spectral_ratio"="SDR","cv_psd"="CV PSD","fmax_psd"="Frac. max"
  ))

set.seed(123)

p=ggplot(NULL)+
  geom_line(data=stat_sim,
            aes(x=System_size,y=value,group=Id_sim),lwd=.3,color="gray",alpha=.4)+
  
  geom_line(data=stat_sim%>%
               group_by(., System_size,variable)%>%
               dplyr::summarise(., .groups = "keep",mean_value=mean(value,na.rm = T)),
             aes(x=System_size,y=mean_value),color="red",lwd=1)+
  
  geom_point(data=stat_sim%>%
               group_by(., System_size,variable)%>%
               dplyr::summarise(., .groups = "keep",mean_value=mean(value,na.rm = T)),
             aes(x=System_size,y=mean_value),color="red",fill="white",shape=21,size=2.5)+
  facet_wrap(.~variable,scales = "free",nrow = 3)+
  labs(x="System size",y="Mean value across all simulations")+
  the_theme+
  guides(color = guide_legend(override.aes = list(size = 2)),fill="none")

ggsave("../Figures/Final_figs/SI/Change_statistics_system_size.pdf",p,width = 8,height = 6)




# >> 4) Comparison data-model : PCA and densities ----


# Density data & model
stat_sim=read.table("../Data_new/All_new_sim2.csv",sep=";")%>%
  sample_n(., 400000)%>%
  mutate(., Pooling=recode_factor(Pooling,"1"="Model, no change",
                                  "2" = "Model, x2","3" = "Model, x3","4" = "Model, x4","5" = "Model, x5"))
d_biocom=read.table("../Data_new/biocom_data.csv",sep=";")

all_d=rbind(stat_sim[,-c(1:2)],
            d_biocom[,c(14:ncol(d_biocom))]%>%
              add_column(.,Pooling='Data'))

p=ggplot(all_d%>%
         melt(., id.vars=c("Pooling"))%>%
           mutate(., variable=recode_factor(variable,
                                            "rho_p"="Cover","nb_neigh"="# neighbors","clustering"= "Clustering",
                                            "skewness"="Skewness","variance"="Variance",
                                            "moran_I"="Autocorrelation","Spectral_ratio"="SDR",
                                            "cv_psd"="CV PSD","fmax_psd"="Frac. max","PL_expo"="Exponent p.l." ,
           )))+
  geom_density(aes(x=value,fill=Pooling),alpha=.5)+
  the_theme+
  labs(x="Value",y="Density",fill="")+
  facet_wrap(.~variable,scales = "free",nrow=3)+
  scale_fill_manual(values=c(my_pal(5),"red"))

ggsave("../Figures/Final_figs/SI/Density_model_versus_data.pdf",p,width = 9,height = 7)

# PCA: raw model and data 

d=all_d%>%
  filter(., Pooling %in% c("Model, no change","Data"))
sumstat_name=colnames(d)[1:11]
res.comp=imputePCA(d[,which(colnames(d) %in% sumstat_name)],ncp=3,scale = T) 

if ("completeObs" %in% names(res.comp)){
  res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
}else {
  res.pca=PCA(res.comp, ncp = 3,  graph=F)
}

axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

for (i in 1:3){
  assign(paste0("p",i),
         d%>%
           mutate(., Pooling=recode_factor(Pooling,"Model, no change"="Model"))%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_point(aes(x = PC1, y = PC2, color = Pooling,fill=Pooling),alpha=.5)+
           scale_color_manual(values=c("gray","red"))+
           scale_fill_manual(values=c("gray","red"))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="",fill="")+
           ggtitle("")+guides(shape="none")+
           theme_classic()+theme(legend.position = "bottom")+
           guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")
  )
}

p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                      p2+theme(legend.position = "none"),
                      p3+theme(legend.position = "none"),
                      ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.1))

ggsave(paste0("../Figures/Final_figs/SI/PCA_raw_model_data.pdf"),p,width = 11,height = 5)



#PCA vegetion type versus model

d_biocom=read.table("../Data_new/biocom_data.csv",sep=";")

#coupling with the one made biogeo consideration
classif_biogeo=read.table("../Data_new/Veg_type_biogeo.csv",sep=";")



d=rbind(stat_sim[,-c(1:2,15)]%>%dplyr::sample_n(., 400000)%>%
          add_column(.,own_classif="Simulations",biogeo="Simulations"),
        d_biocom[,c(14:ncol(d_biocom))]%>%
          add_column(.,Pooling='Data')%>%
          add_column(.,own_classif=classif_own$Type,biogeo=classif_biogeo$type))
sumstat_name=colnames(d)[1:11]
res.comp=imputePCA(d[,which(colnames(d) %in% sumstat_name)],ncp=3,scale = T) 

if ("completeObs" %in% names(res.comp)){
  res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
}else {
  res.pca=PCA(res.comp, ncp = 3,  graph=F)
}

axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

for (i in 1:3){
  assign(paste0("p",i),
         d%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_point(aes(x = PC1, y = PC2, color = Pooling,fill=Pooling,size=Pooling,shape=own_classif))+
           scale_size_manual(values=c(rep(.5,5),1.5))+
           scale_shape_manual(values=c(6,11,3,19))+
           scale_color_manual(values=c(my_pal(5),"black"))+
           scale_fill_manual(values=c(my_pal(5),"black"))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),shape="Type of vegetation",
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="Change in resolution",fill="")+
           ggtitle("")+guides()+
           theme_classic()+theme(legend.position = "bottom")+
           guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")+
           theme(legend.box = "vertical")
  )
}

p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                      p2+theme(legend.position = "none"),
                      p3+theme(legend.position = "none"),
                      ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.2))


ggsave(paste0("../Figs/SI/PCA_spatial_resolution_model_and_data_own_classif_veg.pdf"),p,width = 11,height = 6)



# >> 5) ABC-Posteriors ----

## Comparison x-y obs sim with regular/irregular patterns


keeping_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
x_y_stat=read.table(paste0("../Data_new/Inferrence/x_y_stat_all.csv"),sep=";")
x_y_stat=filter(x_y_stat,Site_ID %in% keeping_sites)%>%
  add_column(., Regular=sapply(1:nrow(.),function(x){
    return(d_biocom$Regular_berdugo[which(d_biocom$ID==.$Site_ID[x])])
  }
 )
)

list_plots=list()
name_plot=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
            "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")
index=1
for (i in c(1:11)){
  d_fil=cbind(filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type","Regular")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Sim")%>%dplyr::rename(., value_sim=value),
              filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type","Regular")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Obs")%>%dplyr::rename(., value_obs=value)%>%dplyr::select(., value_obs))
  
  list_plots[[index]]=ggplot(d_fil)+
    geom_point(aes(x=value_obs,y=value_sim,color=as.factor(Regular)),alpha=.75)+the_theme+
    labs(x="",y="",color="Is regular ?")+
    geom_abline(slope=1,intercept = 0,color="black")+
    ggtitle(name_plot[i])+
    theme(title = element_text(size=10))
  
  index=index+1
}

p=annotate_figure(ggarrange(plotlist=list_plots,ncol = 4,nrow = 3,common.legend = T,legend="bottom"),
                  left=text_grob("Closest simulations",rot=90,color="black",size=15,face ="bold",vjust=1,family = "NewCenturySchoolbook"),
                  bottom = text_grob("Observed spatial statistic",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"))

ggsave("../Figures/Final_figs/SI/Inference_stats_regular_irregular.pdf",p,width = 10,height = 8)


## NMRSE for the summary statistics ----


mean_rmse=d%>%
  melt(., id.vars=c("Type"))%>%
  group_by(.,variable,Type)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value,na.rm=T))

p=ggplot(d%>%melt(., id.vars=c("Type")))+
  geom_jitter(aes(x=Type,y=value,color=Type),
              position = position_jitterdodge(jitter.width = 0.3,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=Type,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="",y="NRMSE",color="")+
  facet_wrap(.~variable,scales = "free")+
  the_theme+
  theme(strip.text.x = element_text(size=10),legend.position = "bottom",
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  geom_hline(yintercept = 1)

ggsave("../Figures/Final_figs/SI/NRMSE_sumstats.pdf",p,width = 10,height = 8)


## Correlation parameters and drivers ----
# 
# 
# d=read.table("../Data_new/posterior_param.csv",sep=";",header=T)
# d=tibble(Site=1:345,p_50=apply(d[,1:345],2,median),p_25=apply(d[,1:345],2,quantile,.25),p_75=apply(d[,1:345],2,quantile,.75),
#          q_50=apply(d[,346:690],2,median),q_25=apply(d[,346:690],2,quantile,.25),q_75=apply(d[,346:690],2,quantile,.75))
# d=d%>%add_column(., Sand=d_biocom$Sand,Aridity=d_biocom$Aridity,MF=d_biocom$MF,Cover=d_biocom$Cover)
# keeping_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
# 
# 
# p1=ggplot(d%>%filter(., Site %in% keeping_sites))+
#   geom_pointrange(aes(x=Cover,y=p_50,ymin=p_25,ymax=p_75),shape=21,fill="gray",color="black")+
#   geom_smooth(aes(x=Cover,y=p_50),fill="#D38DEF",color="#782898",alpha=.3)+
#   the_theme+ggtitle("p")+
#   labs(x="Vegetation cover",y="Posterior parameter value")
#   
# p2=ggplot(d%>%filter(., Site %in% keeping_sites))+
#     geom_pointrange(aes(x=Cover,y=q_50,ymin=q_25,ymax=q_75),shape=21,fill="gray",color="black")+
#     geom_smooth(aes(x=Cover,y=q_50),fill="#D38DEF",color="#782898",alpha=.3)+
#     the_theme+ggtitle("q")+
#     labs(x="Vegetation cover",y="Posterior parameter value")
# 
# ggsave("../Figures/Final_figs/SI/Correlation_param_cover.pdf",
#        ggarrange(p1,p2,ncol=2),
#        width = 9,height = 4)
#   
## Building linear mixed models effects with cover, aridity, MF, SR, sand content and facilitation ----

d=read.table("../Data_new/posterior_param.csv",sep=";",header=T)
d_raw=readxl::read_xlsx("../Data_new/biocom_raw.xlsx")
keep_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
d=tibble(Site=1:345,mean_p=apply(d[,1:345],2,mean),sd_p=apply(d[,1:345],2,sd),
         mean_q=apply(d[,346:690],2,mean),sd_q=apply(d[,346:690],2,sd),
         Plot_n=d_biocom$Plot_n)%>%
  filter(., Site %in% keep_sites)%>%
  add_column(., 
             Facilitation=sapply(1:nrow(.),function(x){return(d_raw$Facil[which(d_raw$plotn==.$Plot_n[x])])}),
             SR=sapply(1:nrow(.),function(x){return(d_raw$SR[which(d_raw$plotn==.$Plot_n[x])])})
  )

d2=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")%>%
  group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  filter(., Site %in% keep_sites)

d=cbind(d,d2,tibble(Sand=d_biocom$Sand,Aridity=d_biocom$Aridity,
                    MF=d_biocom$MF,Cover=d_biocom$Cover)%>%
          filter(., c(1:345) %in% keep_sites))

#as the parameters have uncertainty -> Monte Carlo approach

nsim=1000
d_mod=tibble()
for (k in 1:nsim){
  
  d2=tibble(p=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_p[x],d$sd_p[x]))}))[,1],
            q=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_q[x],d$sd_q[x]))}))[,1],
            Size_tipping=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$Size_mean[x],d$Size_sd[x]))}))[,1],
            abs_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$abs_mean[x],d$abs_sd[x]))}))[,1],
            rela_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$relativ_mean[x],d$relativ_sd[x]))}))[,1],
            Sand=(d$Sand-mean(d$Sand,na.rm=T))/sd(d$Sand,na.rm = T),
            Site=d$Site,
            MF=(d$MF-mean(d$MF,na.rm=T))/sd(d$MF,na.rm = T),
            SR=(d$SR-mean(d$SR,na.rm=T))/sd(d$SR,na.rm = T),
            Facilitation=(d$Facilitation-mean(d$Facilitation,na.rm=T))/sd(d$Facilitation,na.rm = T),
            Cover=(d$Cover-mean(d$Cover,na.rm=T))/sd(d$Cover,na.rm = T),
            Aridity=(d$Aridity-mean(d$Aridity,na.rm=T))/sd(d$Aridity,na.rm = T),
            Plot_n=d$Plot_n)
  
  model_p=summary(d2%>% lmer(formula = paste("p ~ Aridity + MF + Sand + Cover + SR + Facilitation + ( 1 | Plot_n)"), data = .))
  model_q=summary(d2%>% lmer(formula = paste("q ~ Aridity + MF + Sand + Cover + SR + Facilitation + ( 1 | Plot_n)"), data = .))
  model_abs=summary(d2%>% lmer(formula = paste("abs_dist ~ Aridity + MF + Sand + Cover + SR + Facilitation + ( 1 | Plot_n)"), data = .))
  model_rela=summary(d2%>% lmer(formula = paste("rela_dist ~ Aridity + MF + Sand + Cover + SR + Facilitation + ( 1 | Plot_n)"), data = .))
  model_size=summary(d2%>% lmer(formula = paste("Size_tipping ~ Aridity + MF + Sand + Cover + SR + Facilitation + ( 1 | Plot_n)"), data = .))
  
  d_mod=rbind(d_mod,tibble(Aridity=c(model_p$coefficients[2,1],model_q$coefficients[2,1],model_abs$coefficients[2,1],
                                     model_rela$coefficients[2,1],model_size$coefficients[2,1]),
                           Multifunctionality=c(model_p$coefficients[3,1],model_q$coefficients[3,1],model_abs$coefficients[3,1],
                                                model_rela$coefficients[3,1],model_size$coefficients[3,1]),
                           Sand=c(model_p$coefficients[4,1],model_q$coefficients[4,1],model_abs$coefficients[4,1],
                                  model_rela$coefficients[4,1],model_size$coefficients[4,1]),
                           Cover=c(model_p$coefficients[5,1],model_q$coefficients[5,1],model_abs$coefficients[5,1],
                                   model_rela$coefficients[5,1],model_size$coefficients[5,1]),
                           SR=c(model_p$coefficients[6,1],model_q$coefficients[6,1],model_abs$coefficients[6,1],
                                model_rela$coefficients[6,1],model_size$coefficients[6,1]),
                           Facil=c(model_p$coefficients[7,1],model_q$coefficients[7,1],model_abs$coefficients[7,1],
                                   model_rela$coefficients[7,1],model_size$coefficients[7,1]),
                           Param=c("p","q","Absolute distance","Relative distance","Size tipping")))  
  
}

p1=ggplot(d_mod%>%melt(., id.vars=c("Param"))%>%filter(., Param %in% c("p","q")))+
  geom_violin(aes(y=variable,x=value,fill=variable),alpha=.8,color="transparent")+
  geom_boxplot(aes(y=variable,x=value,fill=variable),width=.3,alpha=.8)+
  the_theme+
  facet_wrap(.~Param,scales="free",nrow = 2)+
  labs(y="",x="Fixed effect value",fill="Fixed effect  ")+
  geom_vline(xintercept = 0,linetype=9)+
  scale_fill_manual(values=c("lightgrey",  "#c0f176","#51a8afff","#d29d48ff"))+
  guides(fill="none")

p2=ggplot(d_mod%>%melt(., id.vars=c("Param"))%>%filter(., Param %!in% c("p","q")))+
  geom_violin(aes(y=variable,x=value,fill=variable),alpha=.8,color="transparent")+
  geom_boxplot(aes(y=variable,x=value,fill=variable),width=.3,alpha=.8)+
  the_theme+
  facet_wrap(.~Param,scales="free",nrow = 5)+
  labs(y="",x="Fixed effect value",fill="Fixed effect  ")+
  geom_vline(xintercept = 0,linetype=9)+
  scale_fill_manual(values=c("lightgrey",  "#c0f176","#51a8afff","#d29d48ff"))+
  guides(fill="none")

p_tot=ggarrange(p1,p2,nrow=2,labels=letters[1:2],heights = c(1,1.5))

ggsave("../Figures/Final_figs/SI/LME_param_drivers.pdf",p_tot,width=5,height = 11)


## Is spatial structure or cover good predictors of tipping point distance and height ?  ----

d=read.table("../Data_new/posterior_param.csv",sep=";",header=T)
keep_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
d=tibble(Site=1:345,mean_p=apply(d[,1:345],2,mean),sd_p=apply(d[,1:345],2,sd),
         mean_q=apply(d[,346:690],2,mean),sd_q=apply(d[,346:690],2,sd),
         Plot_n=d_biocom$Plot_n,Cover=d_biocom$Cover)%>%
  filter(., Site %in% keep_sites)

d2=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")%>%
  group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  filter(., Site %in% keep_sites)

d=cbind(d,d2,tibble(Sand=d_biocom$Sand,Aridity=d_biocom$Aridity,
                    MF=d_biocom$MF,Cover=d_biocom$Cover)%>%
          filter(., c(1:345) %in% keep_sites))

#as the parameters have uncertainty -> Monte Carlo approach

nsim=1000
d_mod=tibble()
for (k in 1:nsim){
  
  d2=tibble(p=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_p[x],d$sd_p[x]))}))[,1],
                  q=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_q[x],d$sd_q[x]))}))[,1],
                  Size_tipping=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$Size_mean[x],d$Size_sd[x]))}))[,1],
            abs_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$abs_mean[x],d$abs_sd[x]))}))[,1],
            rela_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$relativ_mean[x],d$relativ_sd[x]))}))[,1],
            Site=d$Site,
            Cover=scale(d$Cover)[,1],
            Plot_n=d$Plot_n)
  
  #Comparing cover + spatial structure with only cover to see whether spatial structure helps to indicate distance or
  #size of tipping point
  
  model_abs=(d2%>% lm(formula = paste("abs_dist ~ Cover + p + q"), data = .))
  model_rela=(d2%>% lm(formula = paste("rela_dist ~ Cover + p + q"), data = .))
  model_size=(d2%>% lm(formula = paste("Size_tipping ~ Cover + p + q"), data = .))
  
  model_abs_0=(d2%>% lm(formula = paste("abs_dist ~ Cover"), data = .))
  model_rela_0=(d2%>% lm(formula = paste("rela_dist ~ Cover"), data = .))
  model_size_0=(d2%>% lm(formula = paste("Size_tipping ~ Cover"), data = .))
  
  AIC_abs=AIC(model_abs,model_abs_0)
  AIC_rela=AIC(model_rela,model_rela_0)
  AIC_size=AIC(model_size,model_size_0)
  
  anova_abs=anova(model_abs,model_abs_0)
  anova_rela=anova(model_rela,model_rela_0)
  anova_size=anova(model_size,model_size_0)
  
  d_mod=rbind(d_mod,tibble(AIC_0=c(AIC_abs$AIC[2],AIC_rela$AIC[2],AIC_size$AIC[2]),
                           AIC=c(AIC_abs$AIC[1],AIC_rela$AIC[1],AIC_size$AIC[1]),
                           F_stat=c(anova_abs$F[2],anova_rela$F[2],anova_size$F[2]),
                           pvalue=c(anova_abs$`Pr(>F)`[2],anova_rela$`Pr(>F)`[2],anova_size$`Pr(>F)`[2]),
                           Param=c("Absolute distance","Relative distance","Size tipping")))  
  
}

p1=ggplot(d_mod)+
  geom_histogram(aes(x=AIC_0-AIC),alpha=.5,color="black",fill="gray")+
  the_theme+
  facet_wrap(.~Param,scales="free",nrow = 2)+
  labs(y="",x="Fixed effect value",fill="Fixed effect  ")+
  geom_vline(xintercept = 0,linetype=9)+
  scale_fill_manual(values=c("lightgrey",  "#c0f176","#51a8afff","#d29d48ff"))+
  guides(fill="none")

## Comparison posterior of landscapes from the same site ----

keep_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
d=read.table("../Data_new/posterior_param.csv",sep=";") #posterior
d2=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";") #predicted metrics
d_biocom=read.table("../Data_new/biocom_data.csv",sep=";")%>%
  filter(., c(1:nrow(.)) %in% keep_sites)

d_summarized=cbind(d2%>%
  group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   mean_rela=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   sd_rela=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   mean_abs=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   sd_abs=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   mean_size=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   sd_size=sd((pinfer-pcrit)/pinfer,na.rm = T))%>%
  filter(., Site %in% keep_sites),
  tibble(p_=colMeans(d)[keep_sites],
         q_=colMeans(d)[keep_sites+345],
         p_sd=apply(d,2,sd,na.rm=T)[keep_sites],
         q_sd=apply(d,2,sd,na.rm=T)[keep_sites+345])
  )%>%
  add_column(., Plot_n=d_biocom$Plot_n[.$Site])

d2=tibble()
for (k in unique(d_biocom$Plot_n)){ #each site
    
  for (rep_id in 1:100){
    
    
    #First parameters 
    
    
    p_site=sapply(which(d_biocom$Plot_n==k),function(x){
      return(rnorm(1,mean=d_summarized$p_[x],
                   sd = d_summarized$p_sd[x]))
    }) #get the values of p
      
    q_site=sapply(which(d_biocom$Plot_n==k),function(x){
      return(rnorm(1,mean=d_summarized$q_[x],
                   sd = d_summarized$q_sd[x]))
    }) #get the values of p
    
    p_other=sapply(which(d_biocom$Plot_n!=k),function(x){
      return(rnorm(1,mean=d_summarized$p_[x],sd=d_summarized$p_sd[x]))
    }) #get the values of p aside this site
    q_other=sapply(which(d_biocom$Plot_n!=k),function(x){
      return(rnorm(1,mean=d_summarized$q_[x],sd=d_summarized$q_sd[x]))
    }) #get the values of q aside this site
    
    RMSE_p_within = mean(sapply(1:length(p_site),function(x){
      return(sqrt(sum((p_site-p_site[x])**2,na.rm = T)/length(p_site) ))})) #NRMSE of p within
    
    RMSE_q_within = mean(sapply(1:length(q_site),function(x){
      return(sqrt(sum((q_site-q_site[x])**2,na.rm = T)/length(q_site) ))})) #NRMSE of q within
    
    RMSE_p_all = mean(sapply(1:length(p_site),function(x){
      return(sqrt(sum((p_other-p_site[x])**2,na.rm = T)/length(p_other) ))})) #NRMSE of p between
    RMSE_q_all = mean(sapply(1:length(q_site),function(x){
      return(sqrt(sum((q_other-q_site[x])**2,na.rm = T)/length(q_other) ))})) #NRMSE of q between
    
    NRMSE_p=RMSE_p_within/RMSE_p_all
    NRMSE_q=RMSE_q_within/RMSE_q_all
    
    
    
    #second for resilience related predictions
    
    #1) relative distance
    
    dist_site=sapply(which(d_biocom$Plot_n==k),function(x){
      return(rnorm(1,mean=d_summarized$mean_rela[x],
                   sd = d_summarized$sd_rela[x]))
    }) #get the values of distance to desert state within the site
    
    dist_other=sapply(which(d_biocom$Plot_n!=k),function(x){
      return(rnorm(1,mean=d_summarized$mean_rela[x],sd=d_summarized$sd_rela[x]))
    }) #get the values of distance aside this site
    
    if (any (dist_other<0)){
      dist_other=dist_other[-which(dist_other<0)] #can happen
    }
    
    RMSE_dist_within = mean(sapply(1:length(dist_site),function(x){
      return(sqrt(sum((dist_site-dist_site[x])**2,na.rm = T)/length(dist_site) ))})) #NRMSE of p within
    
    RMSE_dist_all = mean(sapply(1:length(dist_site),function(x){
      return(sqrt(sum((dist_other-dist_site[x])**2,na.rm = T)/length(dist_other) ))})) #NRMSE of p between
    
    NRMSE_rela=RMSE_dist_within/RMSE_dist_all
    
    #2) absolute distance
    
    dist_site=sapply(which(d_biocom$Plot_n==k),function(x){
      return(rnorm(1,mean=d_summarized$mean_abs[x],
                   sd = d_summarized$sd_abs[x]))
    }) #get the values of distance to desert state within the site
    
    dist_other=sapply(which(d_biocom$Plot_n!=k),function(x){
      return(rnorm(1,mean=d_summarized$mean_abs[x],sd=d_summarized$sd_abs[x]))
    }) #get the values of distance aside this site
    
    if (any (dist_other<0)){
      dist_other=dist_other[-which(dist_other<0)] #can happen
    }
    
    RMSE_dist_within = mean(sapply(1:length(dist_site),function(x){
      return(sqrt(sum((dist_site-dist_site[x])**2,na.rm = T)/length(dist_site) ))})) #NRMSE of p within
    
    RMSE_dist_all = mean(sapply(1:length(dist_site),function(x){
      return(sqrt(sum((dist_other-dist_site[x])**2,na.rm = T)/length(dist_other) ))})) #NRMSE of p between
    
    NRMSE_abs=RMSE_dist_within/RMSE_dist_all
    
    
    #3) size tipping
    
    dist_site=sapply(which(d_biocom$Plot_n==k),function(x){
      return(rnorm(1,mean=d_summarized$mean_size[x],
                   sd = d_summarized$sd_size[x]))
    }) #get the values of distance to desert state within the site
    
    dist_other=sapply(which(d_biocom$Plot_n!=k),function(x){
      return(rnorm(1,mean=d_summarized$mean_size[x],sd=d_summarized$sd_size[x]))
    }) #get the values of distance aside this site
    
    if (any (dist_other<0)){
      dist_other=dist_other[-which(dist_other<0)] #can happen
    }
    
    RMSE_dist_within = mean(sapply(1:length(dist_site),function(x){
      return(sqrt(sum((dist_site-dist_site[x])**2,na.rm = T)/length(dist_site) ))})) #NRMSE of p within
    
    RMSE_dist_all = mean(sapply(1:length(dist_site),function(x){
      return(sqrt(sum((dist_other-dist_site[x])**2,na.rm = T)/length(dist_other) ))})) #NRMSE of p between
    
    NRMSE_size=RMSE_dist_within/RMSE_dist_all
    
    
    d2=rbind(d2,tibble(NRMSE_p=NRMSE_p,NRMSE_q=NRMSE_q,
                       NRMSE_size=NRMSE_size,
                       NRMSE_abs=NRMSE_abs,
                       NRMSE_rela=NRMSE_rela,
                       Plot_n=k,ID_rep=rep_id))
    
  }
}

p=d2%>%
  group_by(.,Plot_n)%>%
  dplyr::summarise(., .groups = "keep",
                   q_median=quantile(NRMSE_p,.5),p_median=quantile(NRMSE_q,.5),
                   rela_median=quantile(NRMSE_rela,.5),abs_median=quantile(NRMSE_abs,.5),
                   size_median=quantile(NRMSE_size,.5))%>%
  melt(., id.vars=c("Plot_n"))%>%
  mutate(., variable=recode_factor(variable,
                                   "q_median"="q",
                                   "p_median"="p",
                                   "rela_median"="Relative dist.",
                                   "abs_median"="Absolute dist.",
                                   "size_median"="Size tipping"))%>%
  ggplot(.)+
  geom_violin(aes(x=variable,y=value,fill=variable),
              color="transparent",alpha=.5,width=.3)+
  geom_boxplot(aes(x=variable,y=value,fill=variable),
               width=.1,alpha=.5,color="black")+
  the_theme+
  labs(x="",y="RMSE within / between sites")+
  geom_hline(yintercept = 1,color="black")+
  scale_color_manual(values=c("#E46767","#7B3B3B","#A6D67E","#B471D8","#7095D2"))+
  scale_fill_manual(values=c("#E46767","#7B3B3B","#A6D67E","#B471D8","#7095D2"))+
  theme(legend.position = "none")

ggsave("../Figures/Final_figs/SI/Comparison_within_between_sites.pdf",p,width = 6,height = 3)




## Observed versus simulated spatial statistics ----

keeping_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
x_y_stat=read.table(paste0("../Data_new/Inferrence/x_y_stat_all.csv"),sep=";")
x_y_stat=filter(x_y_stat,Site_ID %in% keeping_sites)

list_plots=list()
name_plot=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
            "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")
index=1
for (i in c(1:11)){
  d_fil=cbind(filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Sim")%>%dplyr::rename(., value_sim=value),
              filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Obs")%>%dplyr::rename(., value_obs=value)%>%dplyr::select(., value_obs))
  
  list_plots[[index]]=ggplot(d_fil)+
    geom_point(aes(x=value_obs,y=value_sim),color="#96C3DC",alpha=.75)+the_theme+
    labs(x="",y="")+
    geom_abline(slope=1,intercept = 0,color="black")+
    ggtitle(name_plot[i])+
    theme(title = element_text(size=10))
  
  index=index+1
}

p=annotate_figure(ggarrange(plotlist=list_plots,ncol = 4,nrow = 3),
                  left=text_grob("Closest simulations",rot=90,color="black",size=15,face ="bold",vjust=1,family = "NewCenturySchoolbook"),
                  bottom = text_grob("Observed spatial statistic",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"))
ggsave("../Figures/Final_figs/SI/Inference_stats.pdf",p,width = 10,height = 8)



# >> 6) Bimodality in posterior sites ----

site=15
post=read.table("../Data_new/posterior_param.csv",sep=";")[,c(site,site+345)]

p_land1=ggplot(Get_empirical_site(site)%>%melt(.))+
  geom_tile(aes(Var1,Var2,fill=as.factor(value)))+
  scale_fill_manual(values=c("1"="black","0"="white"))+
  theme_transparent()+
  theme(legend.position="none")

colnames(post)=c("p","q")
p_post1=ggplot(post%>%melt(.)%>%
                mutate(., variable=as.character(variable)))+
  geom_histogram(aes(value),color="black",fill="gray")+
  labs(x="Parameter value",y="Count")+
  facet_wrap(.~variable,scales = "free",labeller = label_bquote(cols = Parameters == .(variable) ))+
  the_theme+
  theme(strip.text.x = element_text(size=12),strip.background = element_rect(fill = "#CCE8D8"))



site=20
post=read.table("../Data_new/posterior_param.csv",sep=";")[,c(site,site+345)]

p_land2=ggplot(Get_empirical_site(site)%>%melt(.))+
  geom_tile(aes(Var1,Var2,fill=as.factor(value)))+
  scale_fill_manual(values=c("1"="black","0"="white"))+
  theme_transparent()+
  theme(legend.position="none")

colnames(post)=c("p","q")
p_post2=ggplot(post%>%melt(.)%>%
                mutate(., variable=as.character(variable)))+
  geom_histogram(aes(value),color="black",fill="gray")+
  labs(x="Parameter value",y="Count")+
  facet_wrap(.~variable,scales = "free",labeller = label_bquote(cols = Parameters == .(variable) ))+
  the_theme+
  theme(strip.text.x = element_text(size=12),strip.background = element_rect(fill = "#CCE8D8"))


ggsave("../Figures/Final_figs/SI/Example_bimodal_distrib.pdf",
       ggarrange(ggarrange(
         ggarrange(labels = c("","A1",""),
                   ggplot()+theme_void(),
                   p_land1,ggplot()+theme_void(),ncol=3,widths =c(.3,1,.3)),
         p_post1,labels = c("","B1"),nrow = 2,heights = c(1,1.3)),
         ggarrange(
           ggarrange(labels = c("","A2",""),hjust = 1,
                     ggplot()+theme_void(),
                     p_land2,ggplot()+theme_void(),ncol=3,widths =c(.3,1,.3)),
           p_post2,labels = c("","B2"),nrow = 2,heights = c(1,1.3)),ncol=2),
         width = 12,height = 6)


# >> 7) Validating predictions using Kefi model ----

#First the distance predicted by the kefi model


dist_kefi=tibble()
for (site in list.files("../Data_new/Confirm_kefi/Dist_kefi","Dist")){
  
  site_id=as.numeric(gsub(".csv","",strsplit(site,"_")[[1]][3]))
  
  pred=read.table(paste0("../Data_new/Confirm_kefi/Dist_kefi/",site),sep=",")%>%
    filter(., V1>0)
  colnames(pred)=c("f","b","delta","cover")
  pred$cover[pred$cover<.01]=0

  if (max(pred$cover)>.05){
    dist_kefi=rbind(dist_kefi,tibble(Site=site_id,
                                     size_tipping=min(pred$cover[pred$cover !=0]),
                                     abs_dist=max(pred$b)-max(pred$b[pred$cover==0]),
                                     relativ_dist=(max(pred$b)-max(pred$b[pred$cover==0]))/(max(pred$b[pred$cover==0])),
                                     f=unique(pred$f),delta=unique(pred$delta),
                                     Cover=max(pred$cover)))
  }
}

dist_kefi=dist_kefi%>%add_column(.,Bistab=ifelse(.$Site>81,"yes","no"))

dist_eby=tibble();step_size=0.005
for (site in list.files("../Data_new/Confirm_kefi/Dist_Eby","Dist")){
  site_id=as.numeric(gsub("Eby","",gsub(".csv","",strsplit(site,"_")[[1]][3])))
  
  pred=read.table(paste0("../Data_new/Confirm_kefi/Dist_Eby/",site),sep=",")%>%
    filter(., V2>0)
  colnames(pred)=c("p","q","cover")
  pred$cover[pred$cover<0.01] = 0
  
  if (max(pred$cover)>.05){
    index=0;pred$ID_sim=NA
    for (x in 1:nrow(pred)){
      pred$ID_sim[x]=index
      if (pred$p[x]==0){
        index=index+1
      }
    }
    
    p_desert=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      if (any(d_fil$cover>0)){
        return(d_fil$p[max(which(d_fil$cover !=0))]-step_size)
      }else {
        return(NA)
      }
    }) 
    
    p_infer=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      return(d_fil$p[1])
    }) 
    
    q_infer=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      return(d_fil$q[1])
    }) 
    
    max_cover=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      return(d_fil$cover[1])
    }) 
    
    dist_eby=rbind(dist_eby,tibble(Site=site_id,ID_sim=1:length(p_desert),
                                   abs_dist=p_infer-p_desert,
                                   relativ_dist=(p_infer-p_desert)/p_desert,
                                   Cover=max_cover))
  }else{
    dist_eby=rbind(dist_eby,tibble(Site=site_id,ID_sim=0,
                                   abs_dist=0,
                                   relativ_dist=0,
                                   Cover=0))
  }
  
}

#summarize by quantiles, mean and sd
dist_eby=dist_eby%>%arrange(., Site)

d_eby=dist_eby%>%
  dplyr::group_by(., Site)%>%
  dplyr::summarize(., .groups = "keep",
                   cover_q1=quantile(Cover,.25,na.rm=T),cover_q3=quantile(Cover,.75,na.rm=T),
                   Cover=mean(Cover,na.rm=T),
                   abs_q1=quantile(abs_dist,.25,na.rm=T),abs_q3=quantile(abs_dist,.75,na.rm=T),
                   rela_q1=quantile(relativ_dist,.25,na.rm=T),rela_q3=quantile(relativ_dist,.75,na.rm=T),
                   abs_dist=median(abs_dist,na.rm=T),mean_abs_dist=mean(abs_dist,na.rm=T),
                   mean_rela_dist=mean(relativ_dist,na.rm=T),
                   sd_rela_dist=sd(relativ_dist,na.rm = T),
                   relativ_dist=median(relativ_dist,na.rm=T))%>%
  add_column(., Model="Eby")%>%
  arrange(., Site)%>%
  add_column(.,Bistab=ifelse(.$Site>81,"yes","no"))

sd_abs_dist=sapply(unique(dist_eby$Site),function(x){
  return(sd(dist_eby$abs_dist[which(dist_eby$Site==x)],na.rm = T))
})
d_eby$sd_abs_dist=sd_abs_dist

d_kefi=dist_kefi%>%dplyr::select(., -f,-delta)%>%
  add_column(., abs_q1=0,abs_q3=0,rela_q1=0,rela_q3=0,Model="Kefi")%>%
  arrange(., Site)


stats_param_kefi=read.table("../Data_new/Confirm_kefi/Stats_kefi.csv",sep=",")
colnames(stats_param_kefi)=c("r","d","f","m","b","c","delta","rho_p","nb_neigh","clustering","skewness","variance","moran_I","Spectral_ratio","PLR","PL_expo","cv_psd","fmax_psd")

d_eby=d_eby%>%
  filter(., !(Site %in%  c(1:162)[-unique(d_kefi$Site)]))

d_kefi=d_kefi%>%
  filter(., Site %in% unique(d_eby$Site))


# As there is uncertainty around each inferred distance to the desertification point, 
# we assume that the each inference ~ follow a normal distrib with mean = mean posterior samples
# and sd = sd posterior samples. We then compute the spearman correlation and its associated p-value



#catastrophic shift parameters
d_kefi_cat=filter(d_kefi,Bistab=="yes",Site %in% c(1:nrow(stats_param_kefi))[which(stats_param_kefi$delta==.1)])
d_eby_cat=filter(d_eby,Bistab=="yes",Site %in% c(1:nrow(stats_param_kefi))[which(stats_param_kefi$delta==.1)])

#gradual shift parameters
d_kefi_nocat=filter(d_kefi,Bistab=="no",Site %in% c(1:nrow(stats_param_kefi))[which(stats_param_kefi$delta %in% c(.1))])
d_eby_nocat=filter(d_eby,Bistab=="no",Site %in% c(1:nrow(stats_param_kefi))[which(stats_param_kefi$delta %in% c(.1))])


n=1000
d_spearman=tibble()
for (x in 1:n){
    
  #catastrophic shift params --- absolute distance
  
  abs_eby=sapply(1: nrow(d_eby_cat),function(x){
    return(rnorm(1,mean=d_eby_cat$mean_abs_dist[x],sd=d_eby_cat$sd_abs_dist[x]))
  })
  r_spearman=cor.test(d_kefi_cat$abs_dist,abs_eby,method = "spearman",exact = F)
  d_spearman=rbind(d_spearman,tibble(Stat=r_spearman$estimate,Pval=r_spearman$p.value,Type_dist="Abs",Type_bifu="Cat"))
  
  #catastrophic shift params --- relative distance
  rela_eby=sapply(1: nrow(d_eby_cat),function(x){
    return(rnorm(1,mean=d_eby_cat$mean_rela_dist[x],sd=d_eby_cat$sd_rela_dist[x]))
  })
  r_spearman=cor.test(d_kefi_cat$relativ_dist,rela_eby,method = "spearman",exact = F)
  d_spearman=rbind(d_spearman,tibble(Stat=r_spearman$estimate,Pval=r_spearman$p.value,Type_dist="Rela",Type_bifu="Cat"))
  
  #gradual shift params --- absolute distance
  abs_eby=sapply(1: nrow(d_eby_nocat),function(x){
    return(rnorm(1,mean=d_eby_nocat$mean_abs_dist[x],sd=d_eby_nocat$sd_abs_dist[x]))
  })
  r_spearman=cor.test(d_kefi_nocat$abs_dist,abs_eby,method = "spearman",exact = F)
  d_spearman=rbind(d_spearman,tibble(Stat=r_spearman$estimate,Pval=r_spearman$p.value,Type_dist="Abs",Type_bifu="Grad"))
  
  #gradual shift params --- relative distance
  rela_eby=sapply(1: nrow(d_eby_nocat),function(x){
    return(rnorm(1,mean=d_eby_nocat$mean_rela_dist[x],sd=d_eby_nocat$sd_rela_dist[x]))
  })
  r_spearman=cor.test(d_kefi_nocat$relativ_dist,rela_eby,method = "spearman",exact = F)
  d_spearman=rbind(d_spearman,tibble(Stat=r_spearman$estimate,Pval=r_spearman$p.value,Type_dist="Rela",Type_bifu="Grad"))
  
}



pdf("../Figures/Final_figs/SI/Kefi_C_dist_tipping.pdf",width = 8,height = 8)
{
  par(mfrow=c(2,2))
  
  corr_sp=filter(d_spearman, Type_dist=="Abs",Type_bifu=="Cat")
  
  plot(d_kefi_cat$abs_dist,d_eby_cat$abs_dist,col="#82A6E2",lwd=7,
       main="Absolute distance",xlab="Mean abs. dist from Kefi model",
       ylab="Mean abs. dist from Eby model")
  errbar(d_kefi_cat$abs_dist,d_eby_cat$abs_dist,d_eby_cat$abs_q3,d_eby_cat$abs_q1,
         col="#82A6E2",add=T)
  text(.14,0.28,paste0("r = ",round(median(corr_sp$Stat),2),
                       " (",round(quantile(corr_sp$Stat,.05),2),
                        ", ",round(quantile(corr_sp$Stat,.95),2),
                       ")"))
  
  corr_sp=filter(d_spearman, Type_dist=="Rela",Type_bifu=="Cat")
  
  plot(d_kefi_cat$relativ_dist,d_eby_cat$relativ_dist,col="#82A6E2",lwd=7,pch=19,
       main="Relative distance",
       xlab="Mean relat. dist from Kefi model",
       ylab="Mean relat. dist from Eby model")
  errbar(d_kefi_cat$relativ_dist,d_eby_cat$relativ_dist,
         d_eby_cat$rela_q3,d_eby_cat$rela_q1,
         col="#82A6E2",add=T)
  text(.35,0.45,paste0("r = ",round(median(corr_sp$Stat),2),
                      " (",round(quantile(corr_sp$Stat,.05),2),
                      ", ",round(quantile(corr_sp$Stat,.95),2),
                      ")"))
  
  corr_sp=filter(d_spearman, Type_dist=="Abs",Type_bifu=="Grad")
  
  plot(d_kefi_nocat$abs_dist,d_eby_nocat$abs_dist,col="#EAA96C",lwd=7,
       xlab="Mean abs. dist from Kefi model",
       ylab="Mean abs. dist from Eby model",
       main=NULL)
  errbar(d_kefi_nocat$abs_dist,d_eby_nocat$abs_dist,d_eby_nocat$abs_q3,d_eby_nocat$abs_q1,
         col="#EAA96C",add=T)
  text(.13,0.35,paste0("r = ",round(median(corr_sp$Stat),2),
                       " (",round(quantile(corr_sp$Stat,.05),2),
                       ", ",round(quantile(corr_sp$Stat,.95),2),
                       ")"))
  
  corr_sp=filter(d_spearman, Type_dist=="Rela",Type_bifu=="Grad")
  
  plot(d_kefi_nocat$relativ_dist,d_eby_nocat$relativ_dist,col="#EAA96C",lwd=7,
       xlab="Mean relat. dist from Kefi model",
       ylab="Mean relat. dist from Eby model",
       main=NULL)
  errbar(d_kefi_nocat$relativ_dist,d_eby_nocat$relativ_dist,
         d_eby_nocat$rela_q3,d_eby_nocat$rela_q1,
         col="#EAA96C",add=T)
  text(.35,0.62,paste0("r = ",round(median(corr_sp$Stat),2),
                      " (",round(quantile(corr_sp$Stat,.05),2),
                      ", ",round(quantile(corr_sp$Stat,.95),2),
                      ")"))
  
}
dev.off()



#secondly, the correlation in parameters values

param_rej=colMeans(read.table("../Data_new/Confirm_kefi/param_rej.csv",sep=";")[,-1])
stats_param_kefi=read.table("../Data_new/Confirm_kefi/Stats_kefi.csv",sep=",")
colnames(stats_param_kefi)=c("r","d","f","m","b","c","delta","rho_p","nb_neigh","clustering","skewness","variance","moran_I","Spectral_ratio","PLR","PL_expo","cv_psd","fmax_psd")

p=stats_param_kefi%>%
  add_column(., 
             p=param_rej[1:nrow(stats_param_kefi)],
             q=param_rej[(nrow(stats_param_kefi)+1):(2*nrow(stats_param_kefi))])%>%
  melt(., measure.vars=c("b","c","f","delta"))%>%
  dplyr::rename(., "Param_kefi"="variable","Value_kefi"="value")%>%
  mutate(., Value_kefi=as.character(Value_kefi))%>%
  melt(., measure.vars=c("p","q"))%>%
  dplyr::group_by(., Value_kefi,variable,Param_kefi)%>%
  dplyr::summarise(.,.groups="keep",q2=quantile(value,.5,na.rm=T),
                   q1=quantile(value,.25,na.rm=T),q3=quantile(value,.75,na.rm=T))%>%
  ggplot(.)+
  geom_pointrange(aes(x=Value_kefi,y=q2,ymin=q1,ymax=q3))+
  facet_grid(variable~Param_kefi,scales="free",labeller = label_parsed)+the_theme+
  labs(x="Parameter from the Kefi model",y="Inferred parameter p and q")+
  theme(strip.text.x = element_text(size=14),strip.text.y = element_text(size=14))

ggsave("../Figures/Final_figs/SI/Kefi_A_comparizon_params.pdf",p,width = 7,height = 4)



#first comparizon on how much we fitted the statistics from the Kefi model

x_y_stat=read.table("../Data_new/Confirm_kefi/x_y_stat.csv",sep=";")
x_y_stat=x_y_stat%>%filter(., Site_ID %in% d_kefi$Site)%>%
  add_column(., Low_cov=sapply(1:nrow(.), function(x){
    if (x%%2==0){
      if(.$rho_p[x-1]<.11){
        return(T)
      }else {
        return(F)
      }
    }else{
      if(.$rho_p[x]<.11){
        return(T)
      }else {
        return(F)
      }
    }
  }))%>%
  filter(., Low_cov==F)%>%
  add_column(., Bistab=sapply(1:nrow(.),function(x){
    return(d_kefi$Bistab[which(d_kefi$Site==.$Site_ID[x])])
  }))


list_plots=list()
name_plot=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
            "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")
index=1
for (i in c(1:11)){
  d_fil=cbind(filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type","Bistab")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Sim")%>%dplyr::rename(., value_sim=value),
              filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type","Bistab")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Obs")%>%dplyr::rename(., value_obs=value)%>%dplyr::select(., value_obs))
  
  list_plots[[index]]=ggplot(d_fil)+
    geom_point(aes(x=value_obs,y=value_sim,color=Bistab),alpha=.75,size=2)+the_theme+
    labs(x="",y="",color="")+
    geom_abline(slope=1,intercept = 0,color="black")+
    ggtitle(name_plot[i])+
    theme(title = element_text(size=10))+
    scale_color_manual(values=c("#EAA96C","#82A6E2"),labels=c("Gradual change","Abrupt change"))+
    guides(color = guide_legend(override.aes = list(size = 3)))+
    theme(legend.text = element_text(size=13))
  
  index=index+1
}

p=annotate_figure(ggarrange(plotlist=list_plots,ncol = 4,nrow = 3,common.legend = T,legend="bottom"),
                  left=text_grob("Stats in Eby model (closest selected simulations)",rot=90,color="black",size=15,face ="bold",vjust=1,family = "NewCenturySchoolbook"),
                  bottom = text_grob("Stats in Kefi model (virtual observation)",color="black",size=15,face="bold",vjust=-4,family = "NewCenturySchoolbook"))
ggsave("../Figures/Final_figs/SI/Kefi_B_xystats.pdf",p,width = 10,height = 8)





# >> 9) Comparing with classical Eby model ----

#x and y observed stats
keeping_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
x_y_stat=read.table("../Data_new/Eby/x_y_stat.csv",sep=";")
x_y_stat=filter(x_y_stat,Site_ID %in% keeping_sites)

list_plots=list()
name_plot=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
            "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")
index=1
for (i in c(1:11)){
  d_fil=cbind(filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Sim")%>%dplyr::rename(., value_sim=value),
              filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Obs")%>%dplyr::rename(., value_obs=value)%>%dplyr::select(., value_obs))
  
  list_plots[[index]]=ggplot(d_fil)+
    geom_point(aes(x=value_obs,y=value_sim),color="#96C3DC",alpha=.75)+the_theme+
    labs(x="",y="")+
    geom_abline(slope=1,intercept = 0,color="black")+
    ggtitle(name_plot[i])+
    theme(title = element_text(size=10))
  
  index=index+1
}

p=annotate_figure(ggarrange(plotlist=list_plots,ncol = 4,nrow = 3),
                  left=text_grob("Closest simulations",rot=90,color="black",size=15,face ="bold",vjust=1,family = "NewCenturySchoolbook"),
                  bottom = text_grob("Observed spatial statistic",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"))
ggsave("../Figures/Final_figs/SI/x_y_stat_Eby_classic.pdf",p,width = 10,height = 8)

d=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")
keep_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x

# all sites
d_summarized=d%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",abs_dis50=quantile(pinfer-pcrit,na.rm = T,.5),
                   abs_dis25=quantile(pinfer-pcrit,na.rm = T,.25),
                   abs_dis75=quantile(pinfer-pcrit,na.rm = T,.75),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_dis50=quantile((pinfer-pcrit)/pinfer,na.rm = T,.5),
                   relativ_dis25=quantile((pinfer-pcrit)/pinfer,na.rm = T,.25),
                   relativ_dis75=quantile((pinfer-pcrit)/pinfer,na.rm = T,.75),
                   Size_tipping_sd=sd(Size_tipping,na.rm = T),
                   Size_tipping_mean=mean(Size_tipping,na.rm = T),
                   Size_tipping50=quantile(Size_tipping,na.rm = T,.5),
                   Size_tipping25=quantile(Size_tipping,na.rm = T,.25),
                   Size_tipping75=quantile(Size_tipping,na.rm = T,.75))%>%
  filter(., Site %in% keep_sites)

d2=read.table("../Data_new/Eby/Raw_stability_metrics.csv",sep=";")

d_summarized2=d2%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",abs_dis50=quantile(pinfer-pcrit,na.rm = T,.5),
                   abs_dis25=quantile(pinfer-pcrit,na.rm = T,.25),
                   abs_dis75=quantile(pinfer-pcrit,na.rm = T,.75),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_dis50=quantile((pinfer-pcrit)/pinfer,na.rm = T,.5),
                   relativ_dis25=quantile((pinfer-pcrit)/pinfer,na.rm = T,.25),
                   relativ_dis75=quantile((pinfer-pcrit)/pinfer,na.rm = T,.75),
                   Size_tipping_sd=sd(Size_tipping,na.rm = T),
                   Size_tipping_mean=mean(Size_tipping,na.rm = T),
                   Size_tipping50=quantile(Size_tipping,na.rm = T,.5),
                   Size_tipping25=quantile(Size_tipping,na.rm = T,.25),
                   Size_tipping75=quantile(Size_tipping,na.rm = T,.75))%>%
  filter(., Site %in% keep_sites)

d_summarized=filter(d_summarized,Site %in% unique(d_summarized2$Site))


#estimating rank correlation with the uncertainty in x and y

N_sim=3000
d_cor=tibble()
for (sim in 1:N_sim){
  
  #absolute distance
  
  abs_x=sapply(1: nrow(d_summarized),function(x){
    return(rnorm(1,mean=d_summarized$abs_mean[x],sd=d_summarized$abs_sd[x]))
  })
  abs_y=sapply(1: nrow(d_summarized2),function(x){
    return(rnorm(1,mean=d_summarized2$abs_mean[x],sd=d_summarized2$abs_sd[x]))
  })
  r_spearman=cor.test(abs_x,abs_y,method = "spearman",exact = F)
  d_cor=rbind(d_cor,tibble(Stat=r_spearman$estimate,Pval=r_spearman$p.value,Type_dist="Abs"))
  
  
  #relative distance
  
  abs_x=sapply(1: nrow(d_summarized),function(x){
    return(rnorm(1,mean=d_summarized$relativ_mean[x],sd=d_summarized$relativ_sd[x]))
  })
  abs_y=sapply(1: nrow(d_summarized2),function(x){
    return(rnorm(1,mean=d_summarized2$relativ_mean[x],sd=d_summarized2$relativ_sd[x]))
  })
  r_spearman=cor.test(abs_x,abs_y,method = "spearman",exact = F)
  d_cor=rbind(d_cor,tibble(Stat=r_spearman$estimate,Pval=r_spearman$p.value,Type_dist="Rela"))
  
  
  #size tipping
  abs_x=sapply(1: nrow(d_summarized),function(x){
    return(rnorm(1,mean=d_summarized$Size_tipping_mean[x],sd=d_summarized$Size_tipping_sd[x]))
  })
  abs_y=sapply(1: nrow(d_summarized2),function(x){
    return(rnorm(1,mean=d_summarized2$Size_tipping_mean[x],sd=d_summarized2$Size_tipping_sd[x]))
  })
  r_spearman=cor.test(abs_x,abs_y,method = "spearman",exact = F)
  d_cor=rbind(d_cor,tibble(Stat=r_spearman$estimate,Pval=r_spearman$p.value,Type_dist="Size"))
  
}



d_cor%>%
  dplyr::group_by(., Type_dist)%>%
  dplyr::summarise(., .groups = "keep",
                   mean_cor=mean(Stat,na.rm=T),sd_cor=sd(Stat,na.rm = T))


pdf("../Figures/Verification/Comparing_Eby_true_modified.pdf",width = 7,height = 4)


print(ggplot(NULL)+
        geom_pointrange(aes(x=d_summarized$relativ_dis50,xmin=d_summarized$relativ_dis25,xmax=d_summarized$relativ_dis75,
                            y=d_summarized2$relativ_dis50,ymin=d_summarized2$relativ_dis25,ymax=d_summarized2$relativ_dis75),
                        color=alpha("blue",.5))+
        geom_pointrange(aes(xmax = d_summarized$relativ_dis75, xmin =d_summarized$relativ_dis25,
                            x=d_summarized$relativ_dis50,y=d_summarized2$relativ_dis50),
                        color=alpha("blue",.5))+
        theme_classic()+
        geom_abline(slope = 1,intercept = 0)+
        labs(x="Eby model (modified)",y="Eby model (true)")+
        ggtitle("Relative distance"))



print(ggplot(NULL)+
        geom_pointrange(aes(x=d_summarized$abs_dis50,xmin=d_summarized$abs_dis25,xmax=d_summarized$abs_dis75,
                            y=d_summarized2$abs_dis50,ymin=d_summarized2$abs_dis25,ymax=d_summarized2$abs_dis75),
                        color=alpha("blue",.5))+
        geom_pointrange(aes(xmax = d_summarized$abs_dis75, xmin =d_summarized$abs_dis25,
                            x=d_summarized$abs_dis50,y=d_summarized2$abs_dis50),
                        color=alpha("blue",.5))+
        
        theme_classic()+
        geom_abline(slope = 1,intercept = 0)+
        labs(x="Eby model (modified)",y="Eby model (true)")+
        ggtitle("Absolute distance"))


print(ggplot(NULL)+
        geom_pointrange(aes(x=d_summarized$Size_tipping50,xmin=d_summarized$Size_tipping25,xmax=d_summarized$Size_tipping75,
                            y=d_summarized2$Size_tipping50,ymin=d_summarized2$Size_tipping25,ymax=d_summarized2$Size_tipping75),
                        color=alpha("blue",.5))+
        geom_pointrange(aes(xmax = d_summarized$Size_tipping75, xmin =d_summarized$Size_tipping25,
                            x=d_summarized$Size_tipping50,y=d_summarized2$Size_tipping50),
                        color=alpha("blue",.5))+
        theme_classic()+
        geom_abline(slope = 1,intercept = 0)+
        labs(x="Eby model (modified)",y="Eby model (true)")+
        ggtitle("Size tipping")+
        xlim(0,.3))
dev.off()






# >> 10) Climatic projections with RCP 4.5 scenario ----

keep_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
d=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")
list_clim=list.files("../Data_new/Climatic_data",pattern = ".csv")#[-grep("mean",list.files("../Data_new/Climatic_data",pattern = ".csv"))]

proj_clim=as.data.frame(matrix(NA,length(keep_sites),length(list_clim)))
for (x in 1:ncol(proj_clim)){
  proj_clim[,x]=read.table(paste0("../Data_new/Climatic_data/",list_clim[x]),sep=";")[keep_sites,1]
}
colnames(proj_clim)=c("Aridity RCP 4.5","Aridity RCP 8.5","Temperature RCP 4.5","Temperature RCP 8.5")

# summarizing information in each site
d_summarized=d%>%
  group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",abs_dis50=quantile(pinfer-pcrit,na.rm = T,.5),
                   abs_dis25=quantile(pinfer-pcrit,na.rm = T,.25),
                   abs_dis75=quantile(pinfer-pcrit,na.rm = T,.75),
                   relativ_dis50=quantile((pinfer-pcrit)/pinfer,na.rm = T,.5),
                   relativ_dis25=quantile((pinfer-pcrit)/pinfer,na.rm = T,.25),
                   relativ_dis75=quantile((pinfer-pcrit)/pinfer,na.rm = T,.75),
                   Size_tipping50=quantile(Size_tipping,na.rm = T,.5),
                   Size_tipping25=quantile(Size_tipping,na.rm = T,.25),
                   Size_tipping75=quantile(Size_tipping,na.rm = T,.75))%>%
  filter(., Site %in% keep_sites)%>%
  add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[keep_sites])

d_summarized=cbind(d_summarized,proj_clim)


#first pair correlation between size tipping and distance to desert state


p1=ggplot(d_summarized)+
  geom_pointrange(aes(x=Size_tipping50,ymin=relativ_dis25,ymax=relativ_dis75,y=relativ_dis50,
                      fill=Cover,color=Cover),
                  shape=21)+
  geom_pointrange(aes(x=Size_tipping50,xmin=Size_tipping25,xmax=Size_tipping75,y=relativ_dis50,
                      fill=Cover,color=Cover),
                  shape=21)+
  the_theme+
  scale_fill_gradientn(colours = colorRampPalette(c("#D3EFD3","#90C390","#47A747","#126312"))(100))+
  scale_color_gradientn(colours = colorRampPalette(c("#D3EFD3","#90C390","#47A747","#126312"))(100))+
  labs(y="Distance to desertification",x="Height of the tipping point",fill="Vegetation cover")+
  guides(color="none")+
  theme(legend.position = c(.8, .5),legend.key.size = unit(.5, 'cm'))




#Relative distance to desert state

#defining the colors
mat_tmp=Get_color_classif()
d_summarized$sort_relativ_dist=unlist(sapply(d_summarized$relativ_dis50,function(x){
  return(which(sort(d_summarized$relativ_dis50)==x)[1])
}))
d_summarized$sort_aridity_85=unlist(sapply(d_summarized$`Aridity RCP 4.5`,function(x){
  return(which(sort(d_summarized$`Aridity RCP 4.5`)==x)[1])
}))
d_summarized$color=unlist(sapply(1:nrow(d_summarized),function(x){
  return(mat_tmp[d_summarized$sort_aridity_85[x],d_summarized$sort_relativ_dist[x]])
}))

melting_d=d_summarized%>%
  melt(., measure.vars=colnames(d_summarized)[15:(ncol(d_summarized)-3)])
cols=melting_d$color[1:nrow(d_summarized)];names(cols)=1:nrow(d_summarized)


p2=ggplot(melting_d%>%
            filter(., variable %in% c("Aridity RCP 4.5"))%>%arrange(., relativ_dis50))+
  geom_pointrange(aes(value,ymin=relativ_dis25,ymax=relativ_dis75,y=relativ_dis50,
                      fill=as.character(ID),color=as.character(ID)),
                  shape=21)+
  the_theme+
  theme(legend.position = "none")+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(y="Distance to desertification",x="Vulnerability to climate change")


mat_tmp=Get_color_classif(50)

cols=as.character(mat_tmp);names(cols)=paste0(1:length(cols))

p0_1=ggplot(melt(mat_tmp)%>%add_column(., Val_id=1:nrow(.)))+
  geom_tile(aes(Var1,Var2,fill=as.factor(Val_id)))+
  scale_fill_manual(values = cols)+
  theme_transparent()+
  theme(legend.position = "none",text = element_text(face = "plain"))

p2=p2+annotation_custom(grob=ggplotGrob(p0_1),
                        xmin = -0.0002, xmax=0.00035, ymin=.65, ymax=1.05)+
  geom_richtext(aes(x=0.000075,y=.66,label="Vulnerability"),
                label.size = NA,size=2.5,family = "NewCenturySchoolbook")+
  geom_richtext(aes(x=0.000075,y=.60,label="to climate change"),
                label.size = NA,size=2.5,family = "NewCenturySchoolbook")+
  geom_richtext(aes(x=-0.00018,y=.83,label="Dist. to desert state (b)",angle=90),
                label.size = NA,size=2.5,family = "NewCenturySchoolbook")+
  geom_richtext(aes(x=-0.00028,y=.85,label="Height tip. point (c)",angle=90),
                label.size = NA,size=2.5,family = "NewCenturySchoolbook")

# 
# p0_2=ggplot(d_summarized)+
#   geom_histogram(aes(relativ_dis50),color="black",fill="gray",alpha=.5)+
#   the_theme+theme(legend.position = "none")+
#   labs(x="Distance to desertification",y="Count")
# 
# p1=p1+annotation_custom(grob=ggplotGrob(p0_2),
#                         xmin = 0.00107, xmax=0.00177, ymin=.65, ymax=1.05)+
#   geom_richtext(aes(x=0.0011,y=1,label="b"),
#                 label.size = NA,size=5,fontface="bold")
# 


#Size of the tipping point
#defining the colors
mat_tmp=Get_color_classif()
d_summarized$sort_size_tipping=unlist(sapply(d_summarized$Size_tipping50,function(x){
  return(which(sort(d_summarized$Size_tipping50)==x)[1])
}))
d_summarized$sort_aridity_85=unlist(sapply(d_summarized$`Aridity RCP 4.5`,function(x){
  return(which(sort(d_summarized$`Aridity RCP 4.5`)==x)[1])
}))
d_summarized$color=unlist(sapply(1:nrow(d_summarized),function(x){
  return(mat_tmp[d_summarized$sort_aridity_85[x],d_summarized$sort_size_tipping[x]])
}))
melting_d=d_summarized%>%
  melt(., measure.vars=colnames(d_summarized)[15:(ncol(d_summarized)-3)])
cols=melting_d$color[1:nrow(d_summarized)];names(cols)=1:nrow(d_summarized)


p3=ggplot(melting_d%>%
            filter(., variable %in% c("Aridity RCP 4.5"))%>%arrange(., relativ_dis50))+
  geom_pointrange(aes(value,ymin=Size_tipping25,ymax=Size_tipping75,y=Size_tipping50,
                      fill=as.character(ID),color=as.character(ID)),
                  shape=21)+
  the_theme+
  theme(legend.position = "none")+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  labs(y="Height of the tipping point",x="Vulnerability to climate change")


mat_tmp=Get_color_classif(50)

cols=as.character(mat_tmp);names(cols)=paste0(1:length(cols))

p0_1=ggplot(melt(mat_tmp)%>%add_column(., Val_id=1:nrow(.)))+
  geom_tile(aes(Var1,Var2,fill=as.factor(Val_id)))+
  scale_fill_manual(values = cols)+
  theme_transparent()+
  theme(legend.position = "none",text = element_text(face = "plain"))

# p3=p3+annotation_custom(grob=ggplotGrob(p0_1),
#                         xmin = 0, xmax=0.00055, ymin=.65, ymax=1.05)+
#   geom_richtext(aes(x=0.000275,y=.65,label="Vulnerability to \n climate change"),
#                 label.size = NA,size=.65,family = "NewCenturySchoolbook")+
#   geom_richtext(aes(x=-0.00001,y=.85,label="Height tipping point",angle=90),
#                 label.size = NA,size=3,family = "NewCenturySchoolbook")
# 
# p0_2=ggplot(d_summarized)+
#   geom_histogram(aes(Size_tipping50),color="black",fill="gray",alpha=.5)+
#   the_theme+theme(legend.position = "none")+
#   labs(x="Height of the tipping point",y="Count")
# 
# p2=p2+annotation_custom(grob=ggplotGrob(p0_2),
#                         xmin = 0.00107, xmax=0.00177, ymin=.3, ymax=.6)+
#   geom_richtext(aes(x=0.0011,y=1,label="b"),
#                 label.size = NA,size=5,fontface="bold")



p_tot=ggarrange(p1,p2,p3,nrow=3,labels=letters[1:3])
ggsave("../Figures/Final_figs/SI/Predicting_stability_RCP4_5.pdf",p_tot,width = 5,height = 9)




## --------------------------------------------OTHER) Influence of the number of simulations kept ----

### x_y obs-sim ----

list_f=list.files("../Data_new/N_sim_kept","x_y")
d=data.frame()
for (k in list_f){
  d=rbind(d,read.table(paste0("../Data_new/N_sim_kept/",k),sep=";")%>%
            add_column(., Nkeep=as.numeric(gsub(".csv","",strsplit(k,split = "_")[[1]][5]))))
}


pdf("../Figures/Final_figs/SI/x_y_Nkeep_sim_obs.pdf",width = 10,height = 8)
for (x in unique(d$Nkeep)){
  
  x_y_stat=filter(d,Nkeep==x)
  
  list_plots=list()
  name_plot=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
              "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")
  
  for (i in 1:11){
    d_fil=cbind(filter(x_y_stat%>%
                         melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
                  filter(., Type=="Sim")%>%dplyr::rename(., value_sim=value),
                filter(x_y_stat%>%
                         melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
                  filter(., Type=="Obs")%>%dplyr::rename(., value_obs=value)%>%dplyr::select(., value_obs))
    
    list_plots[[i]]=ggplot(d_fil)+
      geom_point(aes(x=value_obs,y=value_sim),color="#96C3DC",alpha=.75)+the_theme+
      labs(x="",y="")+
      geom_abline(slope=1,intercept = 0,color="black")+
      ggtitle(name_plot[i])+
      theme(title = element_text(size=10))
  }
  
  p=annotate_figure(ggarrange(plotlist=list_plots,ncol = 4,nrow = 3),
                    left=text_grob("Closest simulations",rot=90,color="black",size=15,face ="bold",vjust=1,family = "NewCenturySchoolbook"),
                    bottom = text_grob("Observed spatial statistic",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"),
                    top = text_grob(paste0("Nkeep = ",unique(x_y_stat$Nkeep)),color="black",size=15,face="bold",vjust=1,family = "NewCenturySchoolbook"))
  print(p)
  
}
dev.off()


### posterior of p,q ----

list_f=list.files("../Data_new/N_sim_kept","rej")
d=lapply(list.files("../Data_new/N_sim_kept","rej"), function(x){
  return(read.table(paste0("../Data_new/N_sim_kept/",x),sep=";"))
})

pdf("../Figures/Posterior_sites_Nkeep.pdf",width = 20,height = 4)
for (x in 1:345){
  par(mfrow=c(2,5),mar=c(3,3,3,3))
  
  for (k in 1:length(d)){
    nb_keep=gsub(".csv","",strsplit(list_f[k],"_")[[1]][5])
    hist(d[[k]][,x],main=paste0("p, ",nb_keep),col=alpha("blue",.4))
  }
  for (k in 1:length(d)){
    nb_keep=gsub(".csv","",strsplit(list_f[k],"_")[[1]][5])
    hist(d[[k]][,x+345],main=paste0("q, ",nb_keep),col=alpha("green",.4))
  }
}
dev.off()

#Seeing it on a simple plot



d=tibble()
list_f=list.files("../Data_new/N_sim_kept","rej")[-1]
d_all_sumstat=read.table(paste0("../Data_new/N_sim_kept/NRMSE_param_rej_all_100.csv"),sep=";")
for (i in 1:(length(list_f))){
  post=read.table(paste0("../Data_new/N_sim_kept/",list_f[i]),sep=";")
  d=rbind(d,data.frame(N_keep=gsub(".csv","",strsplit(list_f[i],"_")[[1]][5]),
                       p=apply(post[1:345],2,median)-apply(d_all_sumstat[1:345],2,median),
                       q=apply(post[346:690],2,median)-apply(d_all_sumstat[346:690],2,median),
                       p25=apply(post[1:345],2,quantile,.25)-apply(d_all_sumstat[1:345],2,quantile,.25),
                       q25=apply(post[346:690],2,quantile,.25)-apply(d_all_sumstat[346:690],2,quantile,.25),
                       p75=apply(post[1:345],2,quantile,.75)-apply(d_all_sumstat[1:345],2,quantile,.75),
                       q75=apply(post[346:690],2,quantile,.75)-apply(d_all_sumstat[346:690],2,quantile,.75),
                       Site=1:345))
}


p=ggplot(d%>%
           melt(., id.vars=c("N_keep","Site"))%>%
           add_column(., Parameter=rep(rep(rep(c("p","q"),each=345*4),3)))%>%
           add_column(., Quantile=rep(rep(c(50,25,75),each=345*8)))%>%
           mutate(., N_keep=as.numeric(N_keep)))+
  geom_line(aes(x=N_keep,y=value,group=Site),
            color="gray",lwd=.3,alpha=.4)+
  geom_violin(aes(x=N_keep,y=value,color=interaction(N_keep)),width=50)+
  labs(x="",y="Change in parameter values compared \n with 50 simulations kept",color="")+
  facet_grid(Quantile~Parameter,labeller = label_bquote(cols="Parameter"==.(as.character(Parameter)),
                                                        rows="Quantile, "==.(Quantile)))+
  the_theme+
  theme(strip.text.x = element_text(size=10),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values=colorRampPalette(colors=c("#C46FC5","#80BD5C"))(4))+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = c(100,150,200,250))


ggsave(paste0("../Figures/Final_figs/SI/Change_median_p_q_N_sim_kept.pdf"),p,width = 8,height = 10)
























rm(list=ls())
source("./ABC_drylands_function.R")

dir.create("./Figures",showWarnings=F)
dir.create("./Figures/SI",showWarnings=F)

# ---------------------------- Main figures ------------------------------

## >> Fig. 2 Validation models ----

all_d_kefi=readRDS("./Data/Model_confirmation_Kefi/d_for_figure.rds")
for (k in 1:length(all_d_kefi)){assign(names(all_d_kefi)[k],all_d_kefi[[k]])}
param_kefi=read.table("./Data/Model_confirmation_Kefi/Parameters_kefi.csv",sep=";")%>%
  add_column(., Site=1:nrow(.))
colnames(param_kefi)=c("d","r","c","delta","f","b","m","Site")
param_kefi=filter(param_kefi, Site%in% d_kefi$Site)

corr_sp=filter(d_spearman, Type_dist=="Rela",ID_sim==1)

d_fig=d_eby%>%add_column(., 
                         true_dist_abs=d_kefi$abs_dist,
                         true_size_tipping=d_kefi$size_tipping,
                         true_dist_rela=d_kefi$relativ_dist)%>%
  filter(.,f ==unique(.$f)[1])

p1=ggplot(d_fig)+
  geom_pointrange(aes(x=true_dist_rela,y=relativ_dist,ymax=rela_q3,ymin=rela_q1),
                  color="black",fill="#DEC8EE",shape=23,lwd=.8,size=1)+
  the_theme+
  theme(strip.text.x = element_blank())+
  labs(x="Distance in the dryland model",y="Distance by the inference approach")+
  ggtitle(paste0("r (Spearman) = ",round(median(corr_sp$Stat),2),
                 " (",round(quantile(corr_sp$Stat,.05),2),
                 ", ",round(quantile(corr_sp$Stat,.95),2),
                 ")"))+
  theme(title = element_text(size=10),
        axis.title.x = element_text(colour = "#FF3B3B",size=12),
        axis.title.y = element_text(colour = "#2670B5",size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))


all_d_guichard=readRDS("./Data/Model_confirmation_Guichard/d_for_figure.rds")
for (k in 1:length(all_d_guichard)){assign(names(all_d_guichard)[k],all_d_guichard[[k]])}



corr_sp=filter(d_spearman,ID_sim==2)
d_fig=d_eby%>%add_column(.,
                         true_dist_abs=d_guichard$abs_dist,
                         true_dist_rela=d_guichard$relativ_dist)%>%
  filter(.,a0 == .2)

p2=ggplot(d_fig)+
  geom_pointrange(aes(x=true_dist_rela,y=relativ_dist,ymax=rela_q3,ymin=rela_q1),
                  color="black",fill="#DEC8EE",shape=23,lwd=.8,size=1)+
  the_theme+
  labs(x="Distance in the mussel model",y="Distance by the inference approach")+
  ggtitle(paste0("r (Spearman) = ",round(median(corr_sp$Stat),2),
                 " (",round(quantile(corr_sp$Stat,.05),2),
                 ", ",round(quantile(corr_sp$Stat,.95),2),
                 ")"))+
  theme(title = element_text(size=10),
        axis.title.x = element_text(colour = "#FF3B3B",size=12),
        axis.title.y = element_text(colour = "#2670B5",size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))


p_tot=ggarrange(p1,p2,nrow=2,ncol=1)

ggsave("./Figures/Validation_models_examples.pdf",
       p_tot,
       width = 4,height = 7)



plot(param_kefi$b,all_d_kefi$d_eby$abs_dist)


## >> Fig. 3 Drivers of the resilience of drylands ----


d_partial_res=rbind(
  readRDS("./Data/Drivers_stability_metrics_data_uncertainty_with_facilitation.rds")$Partial_res_data,
  readRDS("./Data/Drivers_stability_metrics_data_uncertainty_without_facilitation.rds")$Partial_res_data
)


d_boot_CI=rbind(
  readRDS("./Data/Drivers_stability_metrics_data_uncertainty_with_facilitation.rds")$Boot_effects,
  readRDS("./Data/Drivers_stability_metrics_data_uncertainty_without_facilitation.rds")$Boot_effects
)

d_boot_CI=d_boot_CI%>%
  dplyr::group_by(., Driver_name,Response)%>%
  dplyr::summarise(., .groups = "keep",q2=median(Slopes),q1=quantile(Slopes,.025),q3=quantile(Slopes,.975),pval=get_bootstrapped_pval(Slopes))


id=1
x_vector=c(-1.5,-1.2,2)
y_vector=c(-1.1,-1.2,-1.15)

id_annotate=1
for (k in c("Multifunctionality","Aridity","Facilitation")){
  
  assign(paste0("p1_",id),
         ggplot(d_partial_res%>%filter(., Driver_name==k,Response=="abs_dist"))+
           geom_point(aes(x=Driver_value,Resids),shape=21,color="grey20",fill="#C9ACDE")+the_theme+
           geom_smooth(aes(x=Driver_value,Resids),method = "lm",se = T,color="black",fill="#C9ACDE")+ 
           annotate("text",x=x_vector[id_annotate],y=y_vector[id_annotate],label=paste0("n = ",nrow(d_partial_res%>%
                                                                                                      filter(., Driver_name==k,Response=="abs_dist"))))+
           labs(x=k,y=title_distance)+
           theme(axis.title = element_text(size=13)))
  id_annotate=id_annotate+1
  id=id+1
}

# p0=ggarrange(ggplot()+theme_void(),p0,ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3))
p1=ggarrange(p1_1,p1_2,p1_3,ncol=3,align = "v",labels = LETTERS[1:3],vjust = -.5,
             font.label = list(size = 23))
p1=ggarrange(ggplot()+theme_void(),p1,nrow=2,heights = c(.1,1))
p_tot=ggarrange(p1,ggplot()+theme_void(),nrow=2,labels=c("","D"),heights = c(1,1),font.label = list(size = 20))
ggsave("./Figures/Drivers_resilience_drylands.pdf",
       p_tot,
       height = 8,width = 10)

#For total and direct effects:
Direct_sem=read.table("./Data/Direct_effects_without_facilitation.csv",sep=";")
Total_sem=read.table("./Data/Total_effects_without_facilitation.csv",sep=";")

## >> Fig. 4 Climatic projections ----

keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
d=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")
clim_trend=read.table("./Data/Climatic_data/mean_aridity_trend.csv",sep=";")

# summarizing information in each site
d_summarized=d%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_dis50_log=log(quantile(pinfer-pcrit,na.rm = T,.5)),
                   abs_dis50=(quantile(pinfer-pcrit,na.rm = T,.5)))%>%
  add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[keep_sites])

d_summarized=d_summarized%>%
  add_column(., Proj_aridity=clim_trend$mean_trend[which(clim_trend$RCP=="rcp85")])

set.seed(123)
#kmeans with 5 clusters to better interpretation of the clusters
kmean_sites = kmeans(scale(d_summarized[,c("abs_dis50_log","Proj_aridity")]), 5)
d_summarized=d_summarized%>%add_column(., cluster_id=as.character(kmean_sites$cluster))

p=ggplot(d_summarized)+
  geom_point(aes(x=abs_dis50_log,y=Proj_aridity,fill=cluster_id,
                 size=Site%in%c(91,23,253,2),shape=Site%in%c(91,23,253,2)),
             color="black")+
  scale_fill_manual(values=c("#FFB77C","#FF707B","#BC8DFF","#FDE7BB","#9CECE5"))+
  scale_size_manual(values=c(1.6,3.7))+
  scale_shape_manual(values=c(21,22))+
  annotate("text",
           x=c(-4.8,-4.8,-3,-2.5),
           y=c(1.5e-3,1e-4,1.5e-3,1e-4),
           label=c("High risk","Ecological risk","Climatic risk","Low risk"),
           color=c("#FF707B","#BC8DFF","#FFB77C","#41D8B9"),family="NewCenturySchoolbook")+
  the_theme+
  theme(legend.position="none")+
  labs(x=expression(paste("Distance to desertification",italic(" (Dist"),", log)")),
       y="Projected aridity change (RCP 8.5)")

p=ggarrange(ggplot()+theme_void(),p,ggplot()+theme_void(),nrow=3,heights =  c(.5,2,.1))

ggsave("./Figures/Mapping_vulnerability.pdf",p,width = 5,height = 5)

# ---------------------------- SI figures ------------------------------

# >> 1) Correlation between predictors ----

d=read.table("./Data/posterior_param.csv",sep=";",header=T)
keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
d=tibble(Site=1:345,mean_p=apply(d[,1:345],2,mean),sd_p=apply(d[,1:345],2,sd),median_p=apply(d[,1:345],2,median),
         mean_q=apply(d[,346:690],2,mean),sd_q=apply(d[,346:690],2,sd),median_q=apply(d[,346:690],2,median),
         Plot_n=d_biocom$Plot_n,
         Aridity=d_biocom$Aridity,
         Sand=d_biocom$Sand,
         MF=d_biocom$MF,
         Soil_A=d_biocom$Soil_A,
         Slope=d_biocom$Slope,
         Facilitation=d_biocom$Facilitation,
         SR=d_biocom$SR,
         Long_cos=d_biocom$Long_cos,
         Long_sin=d_biocom$Long_sin,
         Lat=d_biocom$Lat,
         Elevation=d_biocom$Elevation,
         Cover=d_biocom$Cover,
         Groups=d_biocom$Grp_Kefi)%>%
  dplyr::filter(., Site %in% keep_sites)

d2=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_median=median(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pcrit,na.rm = T),
                   relativ_median=median((pinfer-pcrit)/pcrit,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pcrit,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  dplyr::filter(., Site %in% keep_sites)

d=cbind(d,d2)


d2=tibble(p=logit(d$median_p),
          q=logit(d$median_q),
          abs_dist=scale(log(d$abs_median))[,1],
          rela_dist=scale(d$relativ_median)[,1],
          Sand=(d$Sand-mean(d$Sand,na.rm=T))/sd(d$Sand,na.rm = T),
          Site=d$Site,
          moran=scale(d_biocom$moran_I[d$Site])[,1],
          MF=(d$MF-mean(d$MF,na.rm=T))/sd(d$MF,na.rm = T),
          SR=(d$SR-mean(d$SR,na.rm=T))/sd(d$SR,na.rm = T),
          Soil_A=(d$Soil_A-mean(d$Soil_A,na.rm=T))/sd(d$Soil_A,na.rm = T),
          Facilitation=(d$Facilitation-mean(d$Facilitation,na.rm=T))/sd(d$Facilitation,na.rm = T),
          Cover=(d$Cover-mean(d$Cover,na.rm=T))/sd(d$Cover,na.rm = T),
          Aridity=(d$Aridity-mean(d$Aridity,na.rm=T))/sd(d$Aridity,na.rm = T),
          Lat=(d$Lat-mean(d$Lat,na.rm=T))/sd(d$Lat,na.rm = T),
          Long_cos=(d$Long_cos-mean(d$Long_cos,na.rm=T))/sd(d$Long_cos,na.rm = T),
          Long_sin=(d$Long_sin-mean(d$Long_sin,na.rm=T))/sd(d$Long_sin,na.rm = T),
          Elevation=(d$Elevation-mean(d$Elevation,na.rm=T))/sd(d$Elevation,na.rm = T),
          Slope=(d$Slope-mean(d$Slope,na.rm=T))/sd(d$Slope,na.rm = T),
          Plot_n=d$Plot_n)


corr_pred=corr.test(d2[,c(3,5,7,8,11:18)],use = "na.or.complete",adjust = "none")

colnames(corr_pred$r)=rownames(corr_pred$r)=c("Dist. to desertif. point","Sand","Spatial aggregation (Moran I)",
                                              "Multifunctionality","Facilitation","Cover","Aridity",
                                              "Lattitude","Long (cos)","Long (sin)","Elevation","Slope")

corr_pred$r=round(corr_pred$r,2)
corr_pred$r[lower.tri(corr_pred$r)]=NA
diag(corr_pred$r)=NA

p=ggplot(corr_pred$r%>%
           melt(.)%>%
           add_column(., pval=sapply(1:nrow(.), function(x){
             if (is.na(.$value[x])){
               return(NA)
             }else{return(melt(corr_pred$p)$value[x])}
           })))+
  geom_tile(aes(x=Var1,Var2,fill=value))+
  the_theme+
  geom_text(aes(x=Var1,Var2,label=ifelse(pval<.1,"","X")))+
  scale_fill_gradient2(low="red",mid="white",high = "blue",midpoint = 0,na.value = "white")+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="",fill="")

ggsave("./Figures/SI/Correlation_predictors.pdf",p,width = 6.5,height = 6.5)




# >> 2) Examples of the three sites ----

d=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")
keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1

set.seed(5)
# site=c(244,82,168)
#site=c(126,245,86)
site=c(126,86,252)

pred=tibble();index=0
for (k in site){
  d=read.table(paste0("./Data/Prediction/Dist_tipping_",k,".csv"),sep=",")%>%
    filter(., V1 != 0)%>%add_column(., Site=k)%>%filter(., V2==median(V2)) #selecting median of q
  colnames(d)=c("p","q","cover","Site")
  d$ID_sim=NA
  
  for (x in 1:nrow(d)){
    if (d$p[x]==0.005){
      index=index+1
    }
    d$ID_sim[x]=index
  }
  d=d%>%filter(., ID_sim==unique(.$ID_sim)[1])
  
  pred=rbind(pred,d)
}

pred$color="simulated"
for (k in unique(pred$ID_sim)){
  # pred$p[which(pred$ID_sim==k)]=(pred$p[which(pred$ID_sim==k)]-max(pred$p[which(pred$cover==0 & pred$ID_sim==k)]))/max(pred$p[which(pred$cover==0 & pred$ID_sim==k)])
  pred$p[which(pred$ID_sim==k)]=(pred$p[which(pred$ID_sim==k)]-max(pred$p[which(pred$cover==0 & pred$ID_sim==k)]))
  pred$color[which(pred$ID_sim==k & pred$p==max(pred$p[which(pred$ID_sim==k)]))]="Accepted"
}


p1=ggplot(NULL)+
  geom_point(data=pred,aes(p,y=cover,shape=as.factor(Site),color=color),size=2,shape=1)+
  the_theme+
  scale_shape_manual(values=c(11,10,8))+
  scale_color_manual(values=c("#60BB59","black"))+
  labs(y="Vegetation cover",x="Distance to desertification point (p-pc)",shape="")+
  guides(color = "none",size="none")+
  theme(legend.position = "none")+
  xlim(-0.02,.15)+
  theme(axis.text = element_text(size=13))

ggsave("./Figures/SI/Example_inference_3_sites.pdf",p1,width = 5,height = 3)


# >> 3) Relation between cover, q and distance to tipping point ----

d=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")
post_param=read.table("./Data/posterior_param.csv",sep=";")
# summarizing information in each site
d_summarized=d%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",abs_dis50=quantile(pinfer-pcrit,na.rm = T,.5),
                   abs_dis25=quantile(pinfer-pcrit,na.rm = T,.25),
                   abs_dis75=quantile(pinfer-pcrit,na.rm = T,.75),
                   relativ_dis50=quantile((pinfer-pcrit)/pcrit,na.rm = T,.5),
                   relativ_dis25=quantile((pinfer-pcrit)/pcrit,na.rm = T,.25),
                   relativ_dis75=quantile((pinfer-pcrit)/pcrit,na.rm = T,.75),
                   Size_tipping50=quantile(Size_tipping,na.rm = T,.5),
                   Size_tipping25=quantile(Size_tipping,na.rm = T,.25),
                   Size_tipping75=quantile(Size_tipping,na.rm = T,.75),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pcrit,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pcrit,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T)
  )%>%
  arrange(., Site)%>%
  add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[.$Site],
             Plot_n=d_biocom$Plot_n[.$Site],
             mean_p = colMeans(post_param[,.$Site]),
             mean_q = colMeans(post_param[,.$Site+345]),
             sd_p = apply(post_param[,.$Site],2,sd),
             sd_q = apply(post_param[,.$Site+345],2,sd),
             median_p = apply(post_param[,.$Site],2,median),
             median_q = apply(post_param[,.$Site+345],2,median),
             q_25 = apply(post_param[,.$Site+345],2,quantile,.25),
             q_75 = apply(post_param[,.$Site+345],2,quantile,.75))


p1=ggplot(d_summarized)+
  geom_pointrange(aes(x=Cover,ymin=abs_dis25,ymax=abs_dis75,y=abs_dis50,
                      fill=median_q,color=median_q),
                  shape=21)+
  the_theme+
  scale_fill_gradientn(colours = colorRampPalette(c("#D3EFD3","#90C390","#47A747","#126312"))(100))+
  scale_color_gradientn(colours = colorRampPalette(c("#D3EFD3","#90C390","#47A747","#126312"))(100))+
  labs(y="Distance to the desertification point (Dist)",x="Vegetation cover",fill="Median of q")+
  guides(color="none")+ggtitle("Absolute distance")+
  theme(legend.position = c(.2, .7),legend.key.size = unit(.5, 'cm'),axis.text = element_text(size=13),
        axis.title = element_text(size=14))

p2=ggplot(d_summarized)+
  geom_pointrange(aes(x=Cover,ymin=relativ_dis25,ymax=relativ_dis75,y=relativ_dis50,
                      fill=median_q,color=median_q),
                  shape=21)+
  the_theme+
  scale_fill_gradientn(colours = colorRampPalette(c("#D3EFD3","#90C390","#47A747","#126312"))(100))+
  scale_color_gradientn(colours = colorRampPalette(c("#D3EFD3","#90C390","#47A747","#126312"))(100))+
  labs(y="Distance to the desertification point \n (relative distance)",x="Vegetation cover",fill="Median of q")+
  guides(color="none")+ggtitle("Relative distance")+
  theme(legend.position = c(.2, .7),legend.key.size = unit(.5, 'cm'),axis.text = element_text(size=13),
        axis.title = element_text(size=14))

ggsave("./Figures/SI/Relation_cover_q_resilience.pdf",
       ggarrange(p1,p2,nrow=2,labels = LETTERS[1:2],align = "v"),width = 7,height = 8)

# >> 4) Correlation estimated parameters and estimated distance to the tipping point ----

keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
d=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")

# summarizing information in each site
d_summarized=d%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",abs_dis50=quantile(pinfer-pcrit,na.rm = T,.5),
                   abs_dis25=quantile(pinfer-pcrit,na.rm = T,.25),
                   abs_dis75=quantile(pinfer-pcrit,na.rm = T,.75),
                   relativ_dis50=quantile((pinfer-pcrit)/pcrit,na.rm = T,.5),
                   relativ_dis25=quantile((pinfer-pcrit)/pcrit,na.rm = T,.25),
                   relativ_dis75=quantile((pinfer-pcrit)/pcrit,na.rm = T,.75),
                   Size_tipping50=quantile(Size_tipping,na.rm = T,.5),
                   Size_tipping25=quantile(Size_tipping,na.rm = T,.25),
                   Size_tipping75=quantile(Size_tipping,na.rm = T,.75),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pcrit,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pcrit,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[.$Site],
             Plot_n=d_biocom$Plot_n[.$Site])

#then, predictors of stability metrics
d2=read.table("./Data/posterior_param.csv",sep=";",header=T)
d2=tibble(Site=1:345,mean_p=apply(d2[,1:345],2,mean),sd_p=apply(d2[,1:345],2,sd),
          mean_q=apply(d2[,346:690],2,mean),sd_q=apply(d2[,346:690],2,sd),
          median_p=apply(d2[,1:345],2,median),median_q=apply(d2[,346:690],2,median),
          Plot_n=d_biocom$Plot_n)%>%
  filter(., Site %in% d_summarized$Site)

d_summarized=d_summarized%>%
  add_column(.,p=d2$mean_p,q=d2$mean_q,
             mean_p=d2$mean_p,mean_q=d2$mean_q,
             median_p=d2$median_p,median_q=d2$median_q,
             sd_p=d2$sd_p,sd_q=d2$sd_q,
             moran_I=d_biocom$moran_I[keep_sites])


#ploting the relationships between parameters/cover and predictions

p1_1=ggplot(d_summarized)+
  geom_pointrange(aes(median_q,abs_dis50,
                      ymin=abs_dis25,
                      ymax=abs_dis75),fill="#313131",
                  color="black",alpha=.5,shape=21,size=.7)+
  the_theme+
  labs(x="Posterior median of parameter q",y="Distance to the \n desertification point (Dist)")

p1_2=ggplot(d_summarized)+
  geom_pointrange(aes(median_p,abs_dis50,
                      ymin=abs_dis25,
                      ymax=abs_dis75),fill="#313131",
                  color="black",alpha=.5,shape=21,size=.7)+
  the_theme+
  labs(x="Posterior median of parameter p",y="Distance to the \n desertification point (Dist)")

p1_3=ggplot(d_summarized)+
  geom_pointrange(aes(moran_I,abs_dis50,
                      ymin=abs_dis25,
                      ymax=abs_dis75),fill="#313131",
                  color="black",alpha=.5,shape=21,size=.7)+
  the_theme+
  labs(x="Moran I (Spatial auto-correlation)",y="Distance to the \n desertification point (Dist)")

p2_1=ggplot(d_summarized)+
  geom_point(aes(median_q,Cover),fill="#313131",
             color="black",alpha=.5,shape=21,size=3)+
  the_theme+
  labs(x="Posterior median of parameter q",y="Vegetation cover")

p2_2=ggplot(d_summarized)+
  geom_point(aes(median_p,Cover),fill="#313131",
             color="black",alpha=.5,shape=21,size=3)+
  the_theme+
  labs(x="Posterior median of parameter p",y="Vegetation cover")

p2_3=ggplot(d_summarized)+
  geom_point(aes(moran_I,Cover),fill="#313131",
             color="black",alpha=.5,shape=21,size=3)+
  the_theme+
  labs(x="Moran I (Spatial auto-correlation)",y="Vegetation cover")

p_tot=ggarrange(p1_2,p1_1,p1_3,p2_2,p2_1,p2_3,nrow=2,ncol=3,labels = c("A","","","B","",""))

ggsave("./Figures/SI/Predictors_stability.pdf",p_tot,width = 11,height = 6)

# >> 5) Bootstrapped AIC (q, cover, q+cover) ----

d_mod=read.table("./Data/Cover_vs_spatial_structure_data_uncertainty.csv",sep=";")%>%
  mutate(.,  Predictor=recode_factor( Predictor,
                                      "Cover + Spatial \n structure"="Cover + Moran I",
                                      "Spatial \n structure"="Moran I"))

p1=ggplot(d_mod%>%
            melt(., measure.vars=c("AIC"))%>%
            dplyr::filter(., Stability=="Absolute distance"))+
  geom_density(aes(x=value,fill=Predictor,group=Predictor),alpha=.7)+
  geom_vline(data=d_mod%>%
               dplyr::group_by(., Predictor,Stability)%>%
               dplyr::summarise(., .groups = "keep",median_AIC=median(AIC,na.rm = T))%>%
               filter(., Stability=="Absolute distance"),
             aes(xintercept=median_AIC),color="red",lwd=1)+
  scale_fill_manual(values=c("#7B7B7B","#C1BEBE","#313131"))+
  the_theme+
  guides(fill = guide_legend(override.aes = list(size = 2.5)))+
  labs(x=paste0("Model boostraped AIC (Distance to the desertification point)"),y="Count",fill="")+
  theme(legend.text = element_text(size=7.5),
        legend.title = element_blank(),
        axis.text = element_text(size=13),
        legend.position = c(.15,.8))

ggsave("./Figures/SI/Bootstraped_AIC_q_cover.pdf",p1,width = 6,height = 3)

# >> 6) Replicating figure 3 with p or q ----

#With p or q

d_partial_res=rbind(
  readRDS("./Data/Drivers_stability_metrics_data_uncertainty_with_facilitation.rds")$Partial_res_data,
  readRDS("./Data/Drivers_stability_metrics_data_uncertainty_without_facilitation.rds")$Partial_res_data
)

id=1
for (k in c("Multifunctionality","Facilitation","Aridity")){
  
  assign(paste0("p1_",id),
         ggplot(d_partial_res%>%filter(., Driver_name==k,Response=="q"))+
           geom_point(aes(x=Driver_value,Resids),shape=21,color="grey20",fill="#DEC8EE")+the_theme+
           geom_smooth(aes(x=Driver_value,Resids),method = "lm",se = T,color="black",fill="#DEC8EE")+
           labs(x=k,y="Median of posterior q")+theme(axis.title = element_text(size=13)))

  assign(paste0("p2_",id),
         ggplot(d_partial_res%>%filter(., Driver_name==k,Response=="p"))+
           geom_point(aes(x=Driver_value,Resids),shape=21,color="grey20",fill="#DEC8EE")+the_theme+
           geom_smooth(aes(x=Driver_value,Resids),method = "lm",se = T,color="black",fill="#DEC8EE")+
           labs(x=k,y="Median of posterior p")+theme(axis.title = element_text(size=13)))
  
  id=id+1
}

p1=ggarrange(ggarrange(p1_1,p1_2+theme(axis.title.y = element_blank()),
                       p1_3+theme(axis.title.y = element_blank()),ncol=3),
             ggarrange(p2_1,p2_2+theme(axis.title.y = element_blank()),
                       p2_3+theme(axis.title.y = element_blank()),ncol=3),
             nrow=2,labels=letters[1:2])

ggsave("./Figures/SI/Drivers_resilience_drylands_parameters.pdf",
       p1,
       height = 8,width = 8)


# >> 7) Posteriors: median and bimodality ----

#Bimodality posterior

site=88
post=read.table("./Data/posterior_param.csv",sep=";")[,c(site,site+345)]

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



ggsave("./Figures/SI/Example_bimodal_distrib.pdf",
       ggarrange(
         ggarrange(labels = c("","A",""),
                   ggplot()+theme_void(),
                   p_land1,ggplot()+theme_void(),ncol=3,widths =c(.3,1,.3)),
         p_post1,labels = c("","B"),nrow = 2,heights = c(1,1.3)),
       width = 6,height = 6)


#Posterior median
d=read.table("./Data/posterior_param.csv",sep=";",header=T)
keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1

d=tibble(Site=1:345,mean_p=apply(d[,1:345],2,mean),sd_p=apply(d[,1:345],2,sd),median_p=apply(d[,1:345],2,median),
         mean_q=apply(d[,346:690],2,mean),sd_q=apply(d[,346:690],2,sd),median_q=apply(d[,346:690],2,median),
         median_eta=apply(d[,691:1035],2,median),
         Cover=d_biocom$Cover)%>%
  dplyr::filter(., Site %in% keep_sites)%>%
  dplyr::select(.,-Site)

d2=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_median=median(pinfer-pcrit,na.rm = T),
                   relativ_median=median((pinfer-pcrit)/pcrit,na.rm = T))%>%
  dplyr::filter(., Site %in% keep_sites)

d=cbind(d,d2)


p4=ggplot(d)+
  geom_point(aes(x=median_p,y=median_q,color=(log(abs_median)),fill=(log(abs_median))),shape=21)+
  scale_color_viridis_c(option = "A",labels=c("  Close","Far"),breaks=c(-5,-2))+
  scale_fill_viridis_c(option = "A",labels=c("  Close","Far"),breaks=c(-5,-2))+
  the_theme+
  labs(x="Median of estimated p",y="Median of estimated q",
       fill="Distance to the \n desertification point (Dist)  ",
       color="Distance to the \n desertification point (Dist)  ")+
  theme(axis.title = element_text(size=13),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))

p=ggarrange(
  ggplot(d)+
    geom_histogram(aes(x=median_p))+
    the_theme+
    labs(x="Median of p",y="Number of sites"),
  ggplot(d)+
    geom_histogram(aes(x=median_q))+
    the_theme+
    labs(x="Median of q",y="Number of sites"),
  ggplot(d)+
    geom_histogram(aes(x=median_eta))+
    the_theme+
    labs(x="Median of the scale parameter",y="Number of sites"),
  p4,
  ncol = 2,nrow=2,
  labels = LETTERS[1:4]
)

ggsave("./Figures/SI/Posterior_median.pdf",p,width = 8,height = 7)


#Example two sites

site=179
post=read.table("./Data/posterior_param.csv",sep=";")[,c(site,site+345)]

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



site=33
post=read.table("./Data/posterior_param.csv",sep=";")[,c(site,site+345)]

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


ggsave("./Figures/SI/Example_sites_distrib.pdf",
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




# >> 8) Characteristics of empirical data and PCA on simulations ----
## Resolution, geographical repartition

d_biocom=read.table("./Data/data_sites.csv",sep=";")

# spatial resolution

p=ggplot(d_biocom)+
  geom_histogram(aes(x=Nbpixels),fill="lightblue")+
  the_theme+
  labs(x='# of pixels',y="Count")

ggsave("./Figures/SI/Spatial_resolution_data.pdf",p,width = 6,height = 3)


# Distribution of empirical data, map, aridity and sand cover
d_biocom=read.table("./Data/data_sites.csv",sep=";")

world_map = map_data("world")
p=ggplot(NULL) +
  geom_polygon(data=world_map, aes(x = long, y = lat, group = group),
               fill="lightgray", colour = "white")+
  geom_point(data=d_biocom,aes(x=Longitude,y=Lattitude,color=Aridity),size=3)+
  the_theme+
  scale_color_gradientn(colors = my_pal(4))+
  labs(x="Longitude",y="Lattitue",color="Aridity level")

ggsave("./Figures/SI/Map_empirical_sites.pdf",p,width = 6,height = 4)




## Density of summary statistics

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

ggsave(paste0("./Figures/SI/Pair_corr_metrics.pdf"),p,width = 12,height = 12)



# PCA on simulations by coloring by the parameter value p, q

d=read.table("./Data/Simulations.csv",sep=";",header=T)%>%
  filter(., Pooling==1)
colnames(d)[3:(ncol(d)-1)]=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation","SDR","PLR","Exponent p.l.","CV PSD","Frac. max")

sumstat_name=colnames(d)[3:13]
res.comp=imputePCA(d[,which(colnames(d) %in% sumstat_name)],ncp=3,scale = T) 

if ("completeObs" %in% names(res.comp)){
  res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
}else {
  res.pca=PCA(res.comp, ncp = 3,  graph=F)
}


vPC1 = res.pca$var$coord[,1]
vPC2 = res.pca$var$coord[,2]
vPC3 = res.pca$var$coord[,3]
vlabs = rownames(res.pca$var$coord)
vPCs = data.frame(cbind(vPC1,vPC2,vPC3),res.pca$var$contrib)#[c(1,2,3,4,5,6,8,10,11),]
rownames(vPCs) = vlabs
colnames(vPCs) = c("PC1","PC2","PC3","Contrib1","Contrib2","Contrib3")
random_vertical_distrib=runif(nrow(vPCs),-.5,.5)

save_vPC=vPCs

axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))


for (i in 1:3){
  assign(paste0("p1_",i),
         d%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_point(aes(x = PC1, y = PC2, color = p,fill=p),alpha=.5)+
           geom_text_repel(data=vPCs%>%mutate(., PC1=save_vPC[,axes_for_plot$x[i]],PC2=save_vPC[,axes_for_plot$y[i]]),
                           aes(x=PC1*9,y=PC2*8+random_vertical_distrib,label=rownames(vPCs)), size=4)+
           geom_segment(data=vPCs%>%mutate(., PC1=save_vPC[,axes_for_plot$x[i]],PC2=save_vPC[,axes_for_plot$y[i]]),
                        aes(x = 0, y = 0, xend = PC1*8, yend = PC2*8+random_vertical_distrib), arrow = arrow(length = unit(1/2, 'picas')), color = "grey30")+
           scale_color_viridis_c()+
           scale_fill_viridis_c()+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="p",fill="p")+
           ggtitle("")+guides(shape="none")+
           theme_classic()+theme(legend.position = "bottom")+
           guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")+
           xlim(9*min(vPCs[,paste0("PC",axes_for_plot$x[i])])-2,2+9*max(vPCs[,paste0("PC",axes_for_plot$x[i])]))+
           ylim(9*min(vPCs[,paste0("PC",axes_for_plot$y[i])])-2,2+9*max(vPCs[,paste0("PC",axes_for_plot$y[i])]))
         
         
  )
}

for (i in 1:3){
  assign(paste0("p2_",i),
         d%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_point(aes(x = PC1, y = PC2, color = q,fill=q),alpha=.5)+
           geom_segment(data=vPCs%>%mutate(., PC1=save_vPC[,axes_for_plot$x[i]],PC2=save_vPC[,axes_for_plot$y[i]]),
                        aes(x = 0, y = 0, xend = PC1*6, yend = PC2*6+random_vertical_distrib), arrow = arrow(length = unit(1/2, 'picas')), color = "grey30")+
           geom_text_repel(data=vPCs%>%mutate(., PC1=save_vPC[,axes_for_plot$x[i]],PC2=save_vPC[,axes_for_plot$y[i]]),
                           aes(x=PC1*7,y=PC2*6+random_vertical_distrib,label=rownames(vPCs)), size=4)+
           scale_color_viridis_c()+
           scale_fill_viridis_c()+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="q",fill="q")+
           ggtitle("")+guides(shape="none")+
           theme_classic()+theme(legend.position = "bottom")+
           guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")+
           xlim(9*min(vPCs[,paste0("PC",axes_for_plot$x[i])])-2,2+9*max(vPCs[,paste0("PC",axes_for_plot$x[i])]))+
           ylim(9*min(vPCs[,paste0("PC",axes_for_plot$y[i])])-2,2+9*max(vPCs[,paste0("PC",axes_for_plot$y[i])]))  
  )
}

for (i in 1:3){
  assign(paste0("p3_",i),
         d%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           # geom_point(aes(x = PC1, y = PC2, color = q,fill=q),alpha=.5)+
           geom_segment(data=vPCs%>%mutate(., PC1=save_vPC[,axes_for_plot$x[i]],PC2=save_vPC[,axes_for_plot$y[i]])%>%
                          add_column(., color_segment=save_vPC[,paste0("Contrib",i)]),
                        aes(x = 0, y = 0, xend = PC1*6, yend = PC2*6+random_vertical_distrib,color=color_segment), arrow = arrow(length = unit(1/2, 'picas')),lwd=1)+
           geom_text_repel(data=vPCs%>%mutate(., PC1=save_vPC[,axes_for_plot$x[i]],PC2=save_vPC[,axes_for_plot$y[i]]),
                           aes(x=PC1*7,y=PC2*6+random_vertical_distrib,label=rownames(vPCs)), size=4)+
           scale_color_gradientn(colors = c("lightblue","yellow","orange"))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color=paste0("Contribution to PC",i))+
           ggtitle("")+guides(shape="none")+
           theme_classic()+theme(legend.position = "bottom")+
           guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")+
           xlim(9*min(vPCs[,paste0("PC",axes_for_plot$x[i])])-2,2+9*max(vPCs[,paste0("PC",axes_for_plot$x[i])]))+
           ylim(9*min(vPCs[,paste0("PC",axes_for_plot$y[i])])-2,2+9*max(vPCs[,paste0("PC",axes_for_plot$y[i])]))  
  )
}


p_tot=ggarrange(ggarrange(p1_1+ggtitle("Parameter p"),p1_2,p1_3,ncol=3,common.legend = T,legend = "bottom"),
                ggarrange(p2_1+ggtitle("Parameter q"),p2_2,p2_3,ncol=3,common.legend = T,legend = "bottom"),
                ggarrange(p3_1+ggtitle("Contribution of spatial statistics \n to principal components"),p3_2,p3_3,ncol=3,legend="bottom"),
                labels = letters[1:3],nrow=3)

ggsave("./Figures/SI/PCA_on_simulations_p_q.pdf",p_tot,width=13,height=13)



# >> 9) Optimizing ABC ----

## Number of simulations kept

d=tibble()
list_f=list.files("./Data/NRMSE/","NA")[-grep(pattern = "rej",list.files("./Data/NRMSE/","NA"))]
for (k in list_f){
  d=rbind(d,read.table(paste0("./Data/NRMSE/",k),sep=";")%>%
            add_column(., Nkept=as.numeric(gsub(".csv","",strsplit(k,"_")[[1]][3]))))
}

mean_rmse=d%>%
  melt(., id.vars=c("Nkept"))%>%
  dplyr::group_by(.,variable,Nkept)%>%
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
list_f=list.files("./Data/NRMSE/","NA")[grep(pattern = "rej",list.files("./Data/NRMSE/","NA"))]
for (k in list_f){
  d=rbind(d,read.table(paste0("./Data/NRMSE/",k),sep=";")%>%
            add_column(., Nkept=as.numeric(gsub(".csv","",strsplit(k,"_")[[1]][4]))))
}

mean_rmse=d%>%
  melt(., id.vars=c("Nkept"))%>%
  dplyr::group_by(.,variable,Nkept)%>%
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

ggsave(paste0("./Figures/SI/N_sim_kept.pdf"),
       p2,width = 7,height = 4)


## Best summary statistics

d=tibble()
list_f=list.files("./Data/Best_sumstat","RMSE")
all_name=c("All","No PLR","No Exponent p.l.","No PLR & \n Exponent p.l.","No CV PSD","No Frac. max",
           "No CV PSD & \n Frac. max","No CV PSD & \n Exponent p.l. ","No CV PSD, PLR & \n Exponent p.l.",
           "All not sensitive \n to spatial res.","Using PLS")

for (i in 1:length(list_f)){
  d=rbind(d,read.table(paste0("./Data/Best_sumstat/",list_f[i]),sep=";")%>%
            add_column(.,Name=all_name[i]))
}

d=d%>%add_column(., Color_id=sapply(1:nrow(.),function(x){
  if (.$Name[x] == "All"){
    return("ID_1")
  }else if (.$Name[x]=='All not sensitive \n to spatial res.'){
    return("ID_2")
  }else if (.$Name[x]=='Using PLS'){
    return("ID_3")
  }else{
    return("ID_4")
  }
}))

mean_rmse=d%>%
  melt(., id.vars=c("Name","Color_id"))%>%
  dplyr::group_by(.,variable,Name)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value,na.rm=T))%>%
  dplyr::rename(., Parameter=variable)

pA=ggplot(d%>%
           melt(., id.vars=c("Name","Color_id")))+
  geom_jitter(aes(x=Name,y=value,color=Color_id),
              position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=Name,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="",y="NRMSE",color="")+
  facet_wrap(.~Parameter,labeller = label_bquote(cols="Parameter"==.(as.character(Parameter))))+
  the_theme+
  theme(strip.text.x = element_text(size=10),axis.text.x = element_text(angle=60,hjust=1),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  #scale_color_manual(values=my_pal(9))+
  scale_color_manual(values=c("red","lightgreen","lightblue","grey"))+
  theme(legend.position = "none")


# 
# #Having an index of the summary statistics importance using cov-var
# 
# imp=read.table("./Data/Best_sumstat/Importance_stats.csv",sep=";")
# var_stat=read.table("./Data/Best_sumstat/Variance_scaled.csv",sep=";")
# 
# order_var=imp%>%melt(.)%>%arrange(value)%>%add_column(., Order=nrow(.):1)
# 
# pB=ggplot(rbind(imp,var_stat)%>%
#          add_column(., Type=c("% explained var. accounting for covariance",
#                               "prop. of spatial stat.' variance"))%>%
#          melt(., id.vars="Type")%>%
#            add_column(., Order_var=sapply(1:nrow(.),
#                                           function(x){
#                                             return(order_var$Order[which(.$variable[x]==order_var$variable)])}))%>%
#            mutate(., variable=recode_factor(variable,
#                                             "rho_p"="Cover","nb_neigh"="# neighbors","clustering"= "Clustering","skewness"="Skewness","variance"="Variance",
#                                             "moran_I"="Autocorrelation","Spectral_ratio"="SDR","cv_psd"="CV PSD","fmax_psd"="Frac. max",
#                                             "PL_expo"="Exponent PL fit"
#            ))%>%
#            mutate(., variable=fct_reorder(variable,Order_var)))+
#   geom_bar(aes(x=variable,y=value,fill=Type),
#            position="stack", stat="identity",width = .5)+
#   the_theme+
#   scale_fill_manual(values=c("lightgreen","red"))+
#   theme(axis.text.x = element_text(angle=60,hjust=1))+
#   labs(y="Explained variance",x="",fill="")

ggsave("./Figures/SI/Combination_sumstats.pdf",
       # ggarrange(pA,pB,nrow=2,labels=letters[1:2]),
       pA,
       width = 7,height = 4)#8)



# Removing some stats and impact on the posterior distribution

d=tibble()
list_f=list.files("./Data/Inferrence/","NRMSE_param_rej")
all_name=c("No CV PSD","No Exponent p.l.","No PLR & \n Exponent p.l.","No PLR",
           "No CV PSD, PLR & \n Exponent p.l.","No CV PSD, PLR \n Exponent p.l. & PLR")

d_all_sumstat=read.table(paste0("./Data/Inferrence/NRMSE_param_rej_all.csv"),sep=";")
for (i in 1:(length(list_f))){
  post=read.table(paste0("./Data/Inferrence/",list_f[i+1]),sep=";")
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


ggsave(paste0("./Figures/SI/Combination_sumstats_on_p_q.pdf"),p,width = 10,height = 5)


# >> 10) Spatial resolution & system size ---- 
## Change spatial stats with resolution

stat_sim=read.table("./Data/Simulations.csv",sep=";",header=T)%>%
  dplyr::mutate(., Id_sim=rep(1:(nrow(.)/5),each=5))%>%
  dplyr::mutate(., Pooling=as.character(Pooling),rho_p=round(rho_p,5))%>%
  melt(., id.vars=c("Pooling","Id_sim"))%>%
  filter(., variable %!in% c("p","q"))%>%
  dplyr::mutate(., variable=recode_factor(variable,
                                          "rho_p"="Cover","nb_neigh"="# neighbors","clustering"= "Clustering","skewness"="Skewness","variance"="Variance",
                                          "moran_I"="Autocorrelation","Spectral_ratio"="SDR"
  ))

set.seed(123)
p=ggplot(NULL)+
  geom_line(data=stat_sim%>%filter(., Id_sim %in% sample(unique(.$Id_sim),50)),
            aes(x=Pooling,y=value,group=Id_sim),color="gray50",alpha=.3,lwd=.3)+
  geom_line(data=stat_sim%>%
              dplyr::group_by(., Pooling,variable)%>%
              dplyr::summarise(., .groups = "keep",mean_value=mean(value,na.rm = T)),
            aes(x=Pooling,y=mean_value,group=interaction(variable)),lwd=1,color="red")+
  geom_point(data=stat_sim%>%
               dplyr::group_by(., Pooling,variable)%>%
               dplyr::summarise(., .groups = "keep",mean_value=mean(value,na.rm = T)),
             aes(x=Pooling,y=mean_value),color="red",fill="white",shape=21,size=2.5)+
  facet_wrap(.~variable,scales = "free",nrow = 4)+
  labs(x="Change in resolution",y="Mean value across all simulations")+
  scale_x_discrete(labels=c("No change","x2","x3","x4","x5"))+
  the_theme+theme(axis.text.x = element_text(hjust=1,angle=60))+
  guides(color = guide_legend(override.aes = list(size = 2)),fill="none")

ggsave(paste0("./Figures/SI/Change_metrics_spatial_resolution_model.pdf"),p,width = 8,height = 9)

## Robustness inference with spatial resolution

d_RMSE_param=read.table("./Data/Scale_obs_indentifiability/Retrieving_parameters_different_resolution_RMSE_param.csv",sep=";")

mean_rmse_rej=d_RMSE_param%>% #we remove the scale of observation
  mutate(., Scale_obs=as.character(Scale_obs))%>%
  filter(., Method=="Rejection")%>%
  melt(., id.vars=c("Site_ID","Method","Scale_obs"))%>%
  dplyr::group_by(.,Method,Scale_obs,variable)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value,na.rm=T))


for (k in 1:3){
  assign(paste0("p_",k),
         ggplot(d_RMSE_param%>% #we remove the scale of observation
                  mutate(., Scale_obs=as.character(Scale_obs))%>%
                  filter(., Method=="Rejection")%>%
                  melt(., id.vars=c("Site_ID","Method","Scale_obs"))%>%
                  filter(., variable==unique(.$variable)[k])%>%
                  mutate(., variable=recode_factor(variable, "Pooling"="Scale obs.")))+
           geom_jitter(aes(x=Scale_obs,y=value,color=Scale_obs),alpha=.5,size=1,position = position_jitter(height = 0,width = .1))+
           geom_point(data=mean_rmse_rej%>%filter(., variable==c("p","q","Pooling")[k]),
                      aes(x=Scale_obs,y=mean_rmse),
                      color="white",fill="black",shape=24,size=3)+
           labs(x="",y="NRMSE",color="")+
           the_theme+
           guides(Method="none")+
           scale_color_manual(values=my_pal(5))+
           scale_x_discrete(labels = c(expression(paste(eta," = 1")),expression(paste(eta," = 2")),expression(paste(eta," = 3")),
                                       expression(paste(eta," = 4")),expression(paste(eta," = 5"))))+
           guides(color="none"))
}

ggsave("./Figures/SI/NMRSE_consistency_inference_param_scale.pdf",
       ggarrange(ggplot()+theme_void(),
                 ggarrange(p_1+ggtitle(TeX("A) Parameter p")),
                           p_2+ggtitle(TeX("B) Parameter q")),
                           p_3+ggtitle(TeX("C) Parameter \\eta")),nrow=3),nrow=2,heights = c(1,10)),
       width = 4,height = 7)


#Does estimated parameters change with the model of observation ?

x_y_param=read.table("./Data/Scale_obs_indentifiability/Retrieving_parameters_different_resolution_x_y.csv",sep=";")
p1=ggplot(x_y_param%>%dplyr::select(.,p,q,Site_ID,Method,Type,Scale_obs)%>%
            melt(., id.vars=c("Site_ID","Scale_obs","Type","Method"))%>%mutate(., value=as.numeric(value))%>%
            mutate(., variable=recode_factor(variable, "Pooling"="Scale obs."))%>%
            filter(., Type=="Sim",variable != "Scale obs."))+
  geom_line(aes(x=Scale_obs,y=value,group=Site_ID),alpha=.5,lwd=.5,color="gray")+
  facet_wrap(.~variable,labeller = label_bquote(cols = Parameter==.(as.character(variable))))+
  labs(x=substitute(paste("Scale of observation (",eta,")")),y="Median of the posterior distribution",color="")+
  the_theme


for (i in 1:3){
  d_fil=cbind(filter(x_y_param%>%
                       melt(., id.vars=c("Site_ID","Method", "Type","Scale_obs")),variable==colnames(x_y_param)[i])%>%
                filter(., Type=="Sim")%>%dplyr::rename(., value_sim=value),
              filter(x_y_param%>%
                       melt(., id.vars=c("Site_ID","Method", "Type","Scale_obs")),variable==colnames(x_y_param)[i])%>%
                filter(., Type=="Obs")%>%dplyr::rename(., value_obs=value)%>%dplyr::select(., value_obs))
  
  assign(paste0("p_",i),ggplot(d_fil)+
           geom_point(aes(x=value_obs,y=value_sim,color=as.factor(Scale_obs)),alpha=.75)+the_theme+
           labs(x="True parameter value",y="Estimated parameter value",color=TeX("$\\eta"))+
           scale_color_manual(values=my_pal(5))+
           geom_abline(slope=1,intercept = 0,color="black")+
           theme(title = element_text(size=10)))
  
}
p_tot=ggarrange(p1,ggarrange(p_1+ggtitle("Parameter p"),
                             p_2+ggtitle("Parameter q"),
                             p_3+ggtitle(TeX("$\\eta$")),common.legend = T,ncol=3),
                nrow=2,labels=LETTERS[1:2])

ggsave("./Figures/SI/Consistency_inference_param_scale.pdf",p_tot,width = 9,height = 7)


## System size and summary statistics

list_f=list.files("./Data/System_size")
d=tibble()
for (k in list_f){ #to simplify, we average the replicates
  d2=read.table(paste0("./Data/System_size/",k),sep=",")%>%
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
  filter(., variable %in% c("rho_p","nb_neigh","clustering","skewness",
                            "variance","moran_I","Spectral_ratio","cv_psd","fmax_psd",
                            "PLR","PL_expo"))%>%
  mutate(., variable=recode_factor(variable,
                                   "rho_p"="Cover","nb_neigh"="# neighbors","clustering"= "Clustering","skewness"="Skewness","variance"="Variance",
                                   "moran_I"="Autocorrelation","Spectral_ratio"="SDR","cv_psd"="CV PSD","fmax_psd"="Frac. max",
                                   "PL_expo"="Exponent PL fit"
  ))

set.seed(123)

p=ggplot(NULL)+
  geom_line(data=stat_sim,
            aes(x=System_size,y=value,group=Id_sim),lwd=.3,color="gray",alpha=.4)+
  
  geom_line(data=stat_sim%>%
              dplyr::group_by(., System_size,variable)%>%
              dplyr::summarise(., .groups = "keep",mean_value=mean(value,na.rm = T)),
            aes(x=System_size,y=mean_value),color="red",lwd=1)+
  
  geom_point(data=stat_sim%>%
               dplyr::group_by(., System_size,variable)%>%
               dplyr::summarise(., .groups = "keep",mean_value=mean(value,na.rm = T)),
             aes(x=System_size,y=mean_value),color="red",fill="white",shape=21,size=2.5)+
  facet_wrap(.~variable,scales = "free",nrow = 3)+
  labs(x="System size",y="Mean value across all simulations")+
  the_theme+
  guides(color = guide_legend(override.aes = list(size = 2)),fill="none")

ggsave("./Figures/SI/Change_statistics_system_size.pdf",p,width = 8,height = 6)




# >> 11) Comparison data-model : PCA and densities ----

# Density data & model
stat_sim=read.table("./Data/Simulations.csv",sep=";",header=T)%>%
  sample_n(., 30000)%>%
  mutate(., Pooling=recode_factor(Pooling,"1"="Model, no change",
                                  "2" = "Model, x2","3" = "Model, x3","4" = "Model, x4","5" = "Model, x5"))
d_biocom=read.table("./Data/data_sites.csv",sep=";")

all_d=rbind(stat_sim[,-c(1:2)],
            d_biocom[,c(14:24)]%>%
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

ggsave("./Figures/SI/Density_model_versus_data.pdf",p,width = 9,height = 7)

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

ggsave(paste0("./Figures/SI/PCA_raw_model_data.pdf"),p,width = 11,height = 5)


# >> 12) PCA levels of aggregation data ----


d_biocom=read.table("./Data/data_sites.csv",sep=";")
d_sim=read.table("./Data/Simulations.csv",sep=";",header=T)%>%
  mutate(., Pooling=recode_factor(Pooling,"1"="Model, no change",
                                  "2" = "Model, x2","3" = "Model, x3","4" = "Model, x4","5" = "Model, x5"))

set.seed(123)
d=rbind(d_sim[,-c(1:2,15)],
        d_biocom[,c(14:24)]%>%
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

ggsave(paste0("./Figures/SI/PCA_spatial_resolution_model_and_data.pdf"),p,width = 11,height = 5)


# >> 13) ABC-Posteriors ----


## Observed versus simulated spatial statistics ----

keeping_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
x_y_stat=read.table(paste0("./Data/x_y_stat.csv"),sep=";")
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
                  left=text_grob("Selected simulations",rot=90,color="black",size=15,face ="bold",vjust=1,family = "NewCenturySchoolbook"),
                  bottom = text_grob("Observed spatial statistic",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"))
ggsave("./Figures/SI/Inference_stats.pdf",p,width = 10,height = 8)





## Empirical priors ----
d=read.table("./Data/Simulations.csv",sep=";",header = T)

p1=ggplot(d,aes(x=p,y=q))+geom_density_2d_filled(alpha=.7)+labs(x="Parameter p",y="Parameter q")+the_theme+theme(legend.position = "none")

p21=ggplot(d)+geom_density(aes(p),fill="gray30")+the_theme+labs(x="Parameter p",y="Density")+
  ggtitle("Parameter p")
p22=ggplot(d)+geom_density(aes(q),fill="gray30")+the_theme+labs(x="Parameter q",y="Density")+
  ggtitle("Parameter q")

p_tot=ggarrange(p1,ggarrange(p21,p22,nrow=2),labels = LETTERS[1:2],widths = c(1,.7))
ggsave("./Figures/SI/Empirical_priors.pdf",p_tot,width = 9,height = 5)

# >> 14) Validating ---- 
## Validation predictions using Kefi dryland vegetation model ----

all_d_kefi=readRDS("./Data/Model_confirmation_Kefi/d_for_figure.rds")
for (k in 1:length(all_d_kefi)){assign(names(all_d_kefi)[k],all_d_kefi[[k]])}


pos_x=.25
pos_y1=.7

for (k in unique(d_spearman$ID_sim)){
  
  corr_sp=filter(d_spearman, Type_dist=="Rela",ID_sim==k)
  
  d_fig=d_eby%>%add_column(., 
                           true_dist_abs=d_kefi$abs_dist,
                           true_size_tipping=d_kefi$size_tipping,
                           true_dist_rela=d_kefi$relativ_dist)%>%
    filter(.,f ==unique(.$f)[k])
  
  assign(paste0("p1_",k),ggplot(d_fig)+
           geom_pointrange(aes(x=true_dist_rela,y=relativ_dist,ymax=rela_q3,ymin=rela_q1),
                           color="black",fill="white",shape=23,lwd=.8,size=1)+
           the_theme+
           theme(strip.text.x = element_blank())+
           labs(x="Distance in the dryland model",y="Relative distance \n to the degradation point")+
           ggtitle(paste0("r (Spearman) = ",round(median(corr_sp$Stat),2),
                          " (",round(quantile(corr_sp$Stat,.05),2),
                          ", ",round(quantile(corr_sp$Stat,.95),2),
                          ")"))+
           theme(title = element_text(size=10),
                 axis.title.x = element_text(colour = "#FF3B3B",size=12),
                 axis.title.y = element_text(colour = "#8F45C7",size=12),
                 axis.text.x = element_text(size=12),
                 axis.text.y = element_text(size=12))+
           annotate("text",label=ifelse(k==1,"Low facilitation",ifelse(k==2,"Medium facilitation","High facilitation")),
                    x= min(d_fig$true_dist_rela) + pos_x * diff(range(d_fig$true_dist_rela)),
                    y=min(d_fig$relativ_dist) + pos_y1 * diff(range(d_fig$relativ_dist))))
  
  
  
  corr_sp=filter(d_spearman, Type_dist=="Abs",ID_sim==k)
  
  d_fig=d_eby%>%add_column(., 
                           true_dist_abs=d_kefi$abs_dist,
                           true_size_tipping=d_kefi$size_tipping,
                           true_dist_rela=d_kefi$relativ_dist)%>%
    filter(.,f ==unique(.$f)[k])
  
  assign(paste0("p2_",k),ggplot(d_fig)+
           geom_pointrange(aes(x=true_dist_abs,y=abs_dist,ymax=abs_q3,ymin=abs_q1),
                           color="black",fill="white",shape=23,lwd=.8,size=1)+
           the_theme+
           theme(strip.text.x = element_blank())+
           labs(x="Distance in the dryland model",y=expression(paste("Distance to the degradation point", italic(" (Dist)"))))+
           ggtitle(paste0("r (Spearman) = ",round(median(corr_sp$Stat),2),
                          " (",round(quantile(corr_sp$Stat,.05),2),
                          ", ",round(quantile(corr_sp$Stat,.95),2),
                          ")"))+
           theme(title = element_text(size=10),
                 axis.title.x = element_text(colour = "#FF3B3B",size=12),
                 axis.title.y = element_text(colour = "#8F45C7",size=12),
                 axis.text.x = element_text(size=12),
                 axis.text.y = element_text(size=12))+
           annotate("text",label=ifelse(k==1,"Low facilitation","High facilitation"),
                    x= min(d_fig$true_dist_abs) + pos_x * diff(range(d_fig$true_dist_abs)),
                    y=min(d_fig$abs_dist) + pos_y1 * diff(range(d_fig$abs_dist))))
  
  
}

p_tot=ggarrange(ggarrange(p2_1,p2_2,p2_3,ncol=3),
                ggarrange(p1_1,p1_2,p1_3,ncol=3),nrow=2,labels = LETTERS[1:2],align = "hv")

ggsave("./Figures/SI/Dryland_model_dist_tipping.pdf",p_tot,width = 10,height = 7)




#first comparizon on how much we fitted the statistics from the Kefi model

x_y_stat=read.table("./Data/Model_confirmation_Kefi/x_y_stat.csv",sep=";")
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
  filter(., Low_cov==F)


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
    geom_point(aes(x=value_obs,y=value_sim),alpha=.75,size=2,color="#82A6E2")+the_theme+
    labs(x="",y="",color="")+
    geom_abline(slope=1,intercept = 0,color="black")+
    ggtitle(name_plot[i])+
    theme(title = element_text(size=10))+
    guides(color = guide_legend(override.aes = list(size = 3)))+
    theme(legend.text = element_text(size=13))
  
  index=index+1
}

p=annotate_figure(ggarrange(plotlist=list_plots,ncol = 4,nrow = 3,common.legend = T,legend="bottom"),
                  left=text_grob("Stats in Contact process (selected simulations)",rot=90,color="black",size=15,face ="bold",vjust=1,family = "NewCenturySchoolbook"),
                  bottom = text_grob("Stats in the dryland model (virtual observation)",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"))
ggsave("./Figures/SI/Dryland_model_xystat.pdf",p,width = 10,height = 8)


## Validating predictions using Guichard mussel bed model ----

all_d_guichard=readRDS("./Data/Model_confirmation_Guichard/d_for_figure.rds")
for (k in 1:length(all_d_guichard)){assign(names(all_d_guichard)[k],all_d_guichard[[k]])}

model_comb=expand.grid(a0=unique(d_eby$a0),a2=unique(d_eby$a2))

pos_x=.4
pos_x2=.3
pos_y11=.9
pos_y21=1
pos_y12=.9
pos_y22=1

for (k in unique(d_spearman$ID_sim)){
  
  corr_sp=filter(d_spearman, Type_dist=="Rela",ID_sim==k)
  
  d_fig=d_eby%>%add_column(., 
                           true_dist_abs=d_guichard$abs_dist,
                           true_size_tipping=d_guichard$size_tipping,
                           true_dist_rela=d_guichard$relativ_dist)%>%
    filter(.,a0 ==model_comb$a0[k],a2 ==model_comb$a2[k])
  
  assign(paste0("p1_",k),ggplot(d_fig)+
           geom_pointrange(aes(x=true_dist_rela,y=relativ_dist,ymax=rela_q3,ymin=rela_q1),
                           color="black",fill="white",shape=23,lwd=.8,size=1)+
           the_theme+
           theme(strip.text.x = element_blank())+
           labs(x="Distance in the mussel model",y="Relative distance \n to the degradation point")+
           ggtitle(paste0("r (Spearman) = ",round(median(corr_sp$Stat),2),
                          " (",round(quantile(corr_sp$Stat,.05),2),
                          ", ",round(quantile(corr_sp$Stat,.95),2),
                          ")"))+
           theme(title = element_text(size=10),
                 axis.title.x = element_text(colour = "#FF3B3B",size=12),
                 axis.title.y = element_text(colour = "#8F45C7",size=12),
                 axis.text.x = element_text(size=12),
                 axis.text.y = element_text(size=12))+
           annotate("text",label=ifelse(k%in%c(1,2),"Low disturbance",
                                        ifelse(k%in%c(3,4),"Medium disturbance","High disturbance")),
                    x= min(d_fig$true_dist_rela) + pos_x * diff(range(d_fig$true_dist_rela)),
                    y=min(d_fig$relativ_dist) + pos_y11 * diff(range(d_fig$relativ_dist)))+
           annotate("text",label=ifelse(unique(d_fig$a2)==unique(d_eby$a2)[2],"High recovery","Low recovery"),
                    x= min(d_fig$true_dist_rela) + pos_x * diff(range(d_fig$true_dist_rela)),
                    y=min(d_fig$relativ_dist) + pos_y21 * diff(range(d_fig$relativ_dist))))
  
  
  
  corr_sp=filter(d_spearman, Type_dist=="Abs",ID_sim==k)
  
  d_fig=d_eby%>%add_column(., 
                           true_dist_abs=d_guichard$abs_dist,
                           true_size_tipping=d_guichard$size_tipping,
                           true_dist_rela=d_guichard$relativ_dist)%>%
    filter(.,a0 ==model_comb$a0[k],a2 ==model_comb$a2[k])
  
  assign(paste0("p2_",k),ggplot(d_fig)+
           geom_pointrange(aes(x=true_dist_abs,y=abs_dist,ymax=abs_q3,ymin=abs_q1),
                           color="black",fill="white",shape=23,lwd=.8,size=1)+
           the_theme+
           theme(strip.text.x = element_blank())+
           labs(x="Distance in the mussel model",y=expression(paste("Distance to the degradation point", italic(" (Dist)"))))+
           ggtitle(paste0("r (Spearman) = ",round(median(corr_sp$Stat),2),
                          " (",round(quantile(corr_sp$Stat,.05),2),
                          ", ",round(quantile(corr_sp$Stat,.95),2),
                          ")"))+
           theme(title = element_text(size=10),
                 axis.title.x = element_text(colour = "#FF3B3B",size=12),
                 axis.title.y = element_text(colour = "#8F45C7",size=12),
                 axis.text.x = element_text(size=12),
                 axis.text.y = element_text(size=12))+
           annotate("text",label=ifelse(k%in%c(1,2),"Low disturbance",
                                        ifelse(k%in%c(3,4),"Medium disturbance","High disturbance")),
                    x= min(d_fig$true_dist_abs) + pos_x2 * diff(range(d_fig$true_dist_abs)),
                    y=min(d_fig$abs_dist) + pos_y12 * diff(range(d_fig$abs_dist)))+
           annotate("text",label=ifelse(unique(d_fig$a2)==unique(d_eby$a2)[2],"High recovery","Low recovery"),
                    x= min(d_fig$true_dist_abs) + pos_x2 * diff(range(d_fig$true_dist_abs)),
                    y=min(d_fig$abs_dist) + pos_y22 * diff(range(d_fig$abs_dist))))
  
  
}

p_tot=ggarrange(ggarrange(p2_1,p2_2,p2_3,p2_4,p2_5,p2_6,ncol=3,nrow=2),
                ggarrange(p1_1,p1_2,p1_3,p1_4,p1_5,p1_6,ncol=3,nrow=2),nrow=2,labels = LETTERS[1:2])

ggsave("./Figures/SI/Mussel_model_dist_tipping.pdf",p_tot,width = 10,height = 12)



# comparizon on how much we fitted the statistics from the Guichard model

x_y_stat=read.table("./Data/Model_confirmation_Guichard/x_y_stat.csv",sep=";")%>%
  filter(., Site_ID %in% d_guichard$Site) #to avoid full covered landscapes

list_plots=list()
name_plot=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
            "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")
index=1
for (i in c(1:11)){
  d_fil=cbind(filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method","Type")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Sim")%>%dplyr::rename(., value_sim=value),
              filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method","Type")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Obs")%>%dplyr::rename(., value_obs=value)%>%dplyr::select(., value_obs))
  
  list_plots[[index]]=ggplot(d_fil)+
    geom_point(aes(x=value_obs,y=value_sim),alpha=.75,size=3,color="gray")+the_theme+
    labs(x="",y="",color="")+
    geom_abline(slope=1,intercept = 0,color="black")+
    ggtitle(name_plot[i])+
    theme(title = element_text(size=10))+
    guides(color = guide_legend(override.aes = list(size = 3)))+
    theme(legend.text = element_text(size=13))
  
  index=index+1
}

p=annotate_figure(ggarrange(plotlist=list_plots,ncol = 4,nrow = 3,common.legend = T,legend="bottom"),
                  left=text_grob("Stats in Contact process (closest selected simulations)",rot=90,color="black",size=15,face ="bold",vjust=1,family = "NewCenturySchoolbook"),
                  bottom = text_grob("Stats in the mussel model (virtual observation)",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"))

ggsave("./Figures/SI/Mussel_model_xystats.pdf",p,width = 10,height = 8)


## Relating parameters in the two models with infered distance ----

all_d_kefi=readRDS("./Data/Model_confirmation_Kefi/d_for_figure.rds")
for (k in 1:length(all_d_kefi)){assign(names(all_d_kefi)[k],all_d_kefi[[k]])}
param_kefi=read.table("./Data/Model_confirmation_Kefi/Parameters_kefi.csv",sep=";")%>%
  add_column(., Site=1:nrow(.))
colnames(param_kefi)=c("d","r","c","delta","f","b","m","Site")
param_kefi=filter(param_kefi, Site%in% d_kefi$Site)

all_d_guichard=readRDS("./Data/Model_confirmation_Guichard/d_for_figure.rds")
for (k in 1:length(all_d_guichard)){assign(names(all_d_guichard)[k],all_d_guichard[[k]])}
param_guichard=read.table("./Data/Model_confirmation_Guichard/Parameters_guichard.csv",sep=";")%>%
  add_column(., Site=1:nrow(.))%>%
  filter(., Site%in%d_eby$Site)
colnames(param_guichard)=c("d","a0","a2","Site")

p1=ggplot(NULL)+
  geom_line(aes(x=param_kefi$b,y=all_d_kefi$d_eby$abs_dist,
                group=as.factor(param_kefi$f)),color="grey",lwd=.5)+
  geom_point(aes(x=param_kefi$b,y=all_d_kefi$d_eby$abs_dist,
                 fill=as.factor(param_kefi$f),color=as.factor(param_kefi$f)),shape=21,size=3)+
  the_theme+
  scale_fill_manual(values=c("#AAECEB","#66B38C","#2211F9"))+
  scale_color_manual(values=c("#AAECEB","#66B38C","#2211F9"))+
  labs(x="Recruitment of vegetation (b)",y="Distance to the \n desertification point (Dist)",
       fill="Facilitation parameter (f)",color="Facilitation parameter (f)")

p2=ggplot(NULL)+
  geom_line(aes(x=param_guichard$d,y=all_d_guichard$d_eby$abs_dist,
                group=interaction(as.factor(param_guichard$a0),as.factor(param_guichard$a2))),color="grey",lwd=.5)+
  geom_point(aes(x=param_guichard$d,y=all_d_guichard$d_eby$abs_dist,
                 color=as.factor(param_guichard$a0),
                 fill=as.factor(param_guichard$a0),
                 shape=as.factor(param_guichard$a2)),size=3)+
  the_theme+
  scale_fill_manual(values=c("#AAECEB","#66B38C","#2211F9"))+
  scale_color_manual(values=c("#AAECEB","#66B38C","#2211F9"))+
  scale_shape_manual(values=c(21,22))+
  labs(x="Mortality rate of mussels (d)",y="Distance to the \n desertification point (Dist)",
       color="Density-dependent mortality (a0)",fill="Density-dependent mortality (a0)",
       shape="Density-dependent colonization (a2)")+
  theme(legend.box = "vertical")

p=ggarrange(p1,p2,nrow=2,labels = LETTERS[1:2],heights = c(1,1.15))

ggsave("./Figures/SI/Models_parameters_estimated_distance.pdf",p,width =6,height = 8)

## Relating parameters in the two models with estimated parameters in the minimal model ----

all_d_kefi=readRDS("./Data/Model_confirmation_Kefi/d_for_figure.rds")
for (k in 1:length(all_d_kefi)){assign(names(all_d_kefi)[k],all_d_kefi[[k]])}
param_kefi=read.table("./Data/Model_confirmation_Kefi/Parameters_kefi.csv",sep=";")%>%
  add_column(., Site=1:nrow(.))
colnames(param_kefi)=c("d","r","c","delta","f","b","m","Site")
param_kefi=filter(param_kefi, Site%in% d_kefi$Site)
param_eby_kefi=read.table("./Data/Model_confirmation_Kefi/param_rej.csv",sep=";")
param_eby_kefi=tibble(Median_p=apply(param_eby_kefi[,1:60], 2,median),
                      Median_q=apply(param_eby_kefi[,61:120], 2,median),
                      Site=1:60)%>%
  filter(., Site%in% d_kefi$Site)%>%
  add_column(., b=param_kefi$b,f=param_kefi$f)%>%
  filter(., Median_p!=0)




all_d_guichard=readRDS("./Data/Model_confirmation_Guichard/d_for_figure.rds")
for (k in 1:length(all_d_guichard)){assign(names(all_d_guichard)[k],all_d_guichard[[k]])}
param_guichard=read.table("./Data/Model_confirmation_Guichard/Parameters_guichard.csv",sep=";")%>%
  add_column(., Site=1:nrow(.))%>%
  filter(., Site%in%d_eby$Site)
colnames(param_guichard)=c("d","a0","a2","Site")
param_eby_guichard=read.table("./Data/Model_confirmation_Guichard/param_rej.csv",sep=";")
param_eby_guichard=tibble(Median_p=apply(param_eby_guichard[,1:60], 2,median),
                          Median_q=apply(param_eby_guichard[,61:120], 2,median),
                          Site=1:60)%>%
  filter(., Site%in% d_guichard$Site)%>%
  add_column(., d=param_guichard$d,
             a0=param_guichard$a0,
             a2=param_guichard$a2)%>%
  filter(., Median_p!=0)

p1=ggplot(param_eby_kefi)+
  geom_line(aes(x=b,y=Median_p,
                group=as.factor(f)),color="grey",lwd=.5)+
  geom_point(aes(x=b,y=Median_p,
                 fill=as.factor(f),color=as.factor(f)),shape=21,size=3)+
  the_theme+
  scale_fill_manual(values=c("#AAECEB","#66B38C","#2211F9"))+
  scale_color_manual(values=c("#AAECEB","#66B38C","#2211F9"))+
  labs(x="Recruitment of vegetation (b)",y="Median of posterior \n distribution of p",
       fill="Facilitation parameter (f)",color="Facilitation parameter (f)")

p2=ggplot(param_eby_guichard)+
  geom_line(aes(x=d,y=Median_p,
                group=interaction(as.factor(a0),as.factor(a2))),color="grey",lwd=.5)+
  geom_point(aes(x=d,y=Median_p,
                 color=as.factor(a0),
                 fill=as.factor(a0),
                 shape=as.factor(a2)),size=3)+
  the_theme+
  scale_fill_manual(values=c("#AAECEB","#66B38C","#2211F9"))+
  scale_color_manual(values=c("#AAECEB","#66B38C","#2211F9"))+
  scale_shape_manual(values=c(21,22))+
  labs(x="Mortality rate of mussels (d)",y="Median of posterior \n distribution of p",
       color="Density-dependent mortality (a0)",fill="Density-dependent mortality (a0)",
       shape="Density-dependent colonization (a2)")+
  theme(legend.box = "vertical")

p=ggarrange(p1,p2,nrow=2,labels = LETTERS[1:2],heights = c(1,1.15))

ggsave("./Figures/SI/Models_parameters_estimated_parameters.pdf",p,width =6,height = 8)




# >> 15) Characteristics of the climatic clusters ----


keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
d=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")
clim_trend=read.table("./Data/Climatic_data/mean_aridity_trend.csv",sep=";")

# summarizing information in each site
d_summarized=d%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_dis50_log=log(quantile(pinfer-pcrit,na.rm = T,.5)),
                   abs_dis50=(quantile(pinfer-pcrit,na.rm = T,.5)))%>%
  add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[keep_sites],Facilitation=d_biocom$Facilitation[keep_sites])

d_summarized=d_summarized%>%
  add_column(., Proj_aridity=clim_trend$mean_trend[which(clim_trend$RCP=="rcp85")])

set.seed(123)
kmean_sites = kmeans(scale(d_summarized[,c("abs_dis50_log","Proj_aridity")]), 5)


name_driver=c("Multifunctionality","Current aridity level","Facilitation","Distance to the desertification point \n (Dist, log)",
              "Yearly trend in aridity aridity \n (1950-2100)","Vegetation cover")
id_plot=1
for (k in c("MF","aridity","Facilitation","abs_dis50_log","Proj_aridity","Cover","legend")){
  
  if (k =="legend"){
    
    p_legend=ggplot(d_summarized%>%
             add_column(., cluster_id=as.character(kmean_sites$cluster))%>%
             mutate(., cluster_id=recode_factor(cluster_id,
                                                "1"="Climatic risk, \n medium ecological risk",
                                                "2"= "High risk",
                                                "3"="Ecological risk",
                                                '4'='Climatic risk, \n low ecological risk',
                                                "5"="Low risk"))%>%
             melt(., measure.vars="MF"))+
      geom_violin(aes(x=cluster_id,y=value,group=cluster_id,fill=as.character(cluster_id)),alpha=0.8,trim=T)+
      the_theme+
      labs(x="",y="")+
      theme(axis.text.x = element_blank(),legend.position = "bottom",
            axis.ticks.x = element_blank(),axis.title.x = element_blank())+  
      scale_fill_manual(values=c("#FDE7BB","#FFB77C","#BC8DFF","#FF707B","#9CECE5"))
  } else{
    
    assign(paste0("p_",id_plot),ggplot(d_summarized%>%
                                         add_column(., cluster_id=as.character(kmean_sites$cluster))%>%
                                         mutate(., cluster_id=recode_factor(cluster_id,
                                                                            "1"="Climatic risk, \n medium ecological risk",
                                                                            "2"= "High risk",
                                                                            "3"="Ecological risk",
                                                                            '4'='Climatic risk, \n low ecological risk',
                                                                            "5"="Low risk"))%>%
                                         melt(., measure.vars=k))+
             # geom_boxplot(aes(x=cluster_id,y=value,group=cluster_id,fill=as.character(cluster_id)))+
             # geom_jitter(aes(x=cluster_id,y=value,fill=as.character(cluster_id)),height = 0,width = .2,shape=21)+
             geom_violin(aes(x=cluster_id,y=value,group=cluster_id,fill=as.character(cluster_id)),alpha=0.8,trim=T)+
             geom_boxplot(aes(x=cluster_id,y=value,fill=as.character(cluster_id)),width=0.2,outlier.shape = NA,color="black")+
             the_theme+
             # labs(x="",y=name_driver[id_plot],)+
             labs(x="",y=name_driver[id_plot])+
             # theme(axis.text.x = element_text(angle=60,hjust=1),legend.position = "none")+  
             theme(axis.text.x = element_blank(),legend.position = "none",
                   axis.ticks.x = element_blank(),axis.title.x = element_blank())+  
             scale_fill_manual(values=c("#FDE7BB","#FFB77C","#BC8DFF","#FF707B","#9CECE5")))
  }
  id_plot=id_plot+1
}


p_tot=ggarrange(ggarrange(p_4,p_5,p_6,p_1,p_2,p_3,ncol=3,nrow=2),
                get_legend(p_legend+labs(fill="")),nrow=2,heights = c(1,.1))
ggsave("./Figures/SI/Clusters_characteristics.pdf",p_tot,width = 9,height = 6)



# >> 16) Changes in climatic projection with models and RCP scenario ----

keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
d=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")
clim_trend=read.table("./Data/Climatic_data/mean_aridity_trend.csv",sep=";")

# summarizing information in each site
d_summarized=d%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_dis50_log=log(quantile(pinfer-pcrit,na.rm = T,.5)),
                   abs_dis50=(quantile(pinfer-pcrit,na.rm = T,.5)))%>%
  add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[keep_sites],Facilitation=d_biocom$Facilitation[keep_sites])

d_summarized=d_summarized%>%
  add_column(., Proj_aridity=clim_trend$mean_trend[which(clim_trend$RCP=="rcp85")])

#First; optimal number of clusters

Optimal_n_clust=mclust::mclustBIC(scale(d_summarized[,c("abs_dis50_log","Proj_aridity")]))
Optimal_n_clust=as.data.frame(Optimal_n_clust[1:9,])

p=Optimal_n_clust%>%
  add_column(., Nclust=1:9)%>%
  melt(.,id.vars = "Nclust")%>%
  ggplot(.)+
  geom_line(aes(x=Nclust,y=value,group=variable),color="grey",lwd=.5,alpha=.5)+
  geom_line(data=.%>%dplyr::group_by(.,Nclust)%>%dplyr::summarise(., .groups = "keep",mean_val=mean(value)),
            aes(x=Nclust,y=mean_val),color="black",lwd=2)+
  the_theme+
  labs(x="Number of clusters",y="Bayesian Information Criterion (BIC)")
ggsave("./Figures/SI/Optimal_number_clusters.pdf",p,width = 5,height = 4)



#SEcond; Change position clusters with the climatic model 

d_clim=read.table("./Data/Climatic_data/aridity_trend_all_models.csv",sep=";")
d_clim_averaged=read.table("./Data/Climatic_data/mean_aridity_trend.csv",sep=";")%>%
  dplyr::rename(., trend=mean_trend)%>%
  arrange(., RCP)%>%
  add_column(., model=rep((max(d_clim$model)+1):(max(d_clim$model)+2),each=293))

d_clim=rbind(d_clim,
             d_clim_averaged[,c(3,4,1,2)])
d_center=tibble()

for (model_id in unique(d_clim$model)){
  
  keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
  d=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")
  clim_trend=d_clim%>%dplyr::filter(., model==model_id)
  
  # summarizing information in each site
  d_summarized=d%>%
    dplyr::group_by(., Site,MF,aridity,Sand)%>%
    dplyr::summarise(., .groups = "keep",
                     abs_dis50_log=log(quantile(pinfer-pcrit,na.rm = T,.5)),
                     abs_dis50=(quantile(pinfer-pcrit,na.rm = T,.5)))%>%
    add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[keep_sites])
  
  d_summarized=d_summarized%>%
    add_column(., Proj_aridity=clim_trend$trend)
  
  set.seed(123)
  #kmeans with 5 clusters to better interpretation of the clusters
  kmean_sites = kmeans(scale(d_summarized[,c("abs_dis50_log","Proj_aridity")]), 5)
  d_summarized=d_summarized%>%
    add_column(., cluster_id=as.character(kmean_sites$cluster))
  
  d_center=rbind(d_center,kmean_sites$centers%>%
                   as_tibble(.)%>%
                   add_column(., 
                              Model=model_id,
                              RCP=unique(clim_trend$RCP),
                              Cluster=1:5,
                              sd_abs=sd(d_summarized$abs_dis50_log),
                              mean_abs=mean(d_summarized$abs_dis50_log),
                              sd_proj=sd(d_summarized$Proj_aridity),
                              mean_proj=mean(d_summarized$Proj_aridity)))
  
}


d_center$Average_models="No"
d_center$Average_models[which(d_center$Model %in% 19:20)] = "Yes"

p1=ggplot(d_center%>%
         mutate(., Cluster=as.factor(Cluster),
                RCP=recode_factor(RCP,"rcp45"="RCP 4.5 Scenario",
                                  "rcp85"="RCP 8.5 Scenario"))%>%
         mutate(., Cluster=recode_factor(Cluster,
                                            "1"="Climatic risk, \n medium ecological risk",
                                            "2"= "High risk",
                                            "3"="Ecological risk",
                                            '4'='Climatic risk, \n low ecological risk',
                                            "5"="Low risk"))%>%
         mutate(., 
                abs_dis50_log=abs_dis50_log*sd_abs+mean_abs,
                Proj_aridity=Proj_aridity*sd_proj+mean_proj))+
  geom_point(aes(x=abs_dis50_log,y=Proj_aridity,fill=Cluster,color=Cluster,
                 shape=Average_models,size=Average_models))+
  facet_wrap(.~RCP)+
  the_theme2+
  guides(shape="none",size="none",fill="none")+
  scale_size_manual(values=c(2,5))+
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=c("#FFB77C","#FF707B","#BC8DFF","#FDE7BB","#9CECE5"))+
  scale_color_manual(values=c("#FFB77C","#FF707B","#BC8DFF","#FDE7BB","#9CECE5"))+
  labs(x=expression(paste("Distance to the tipping point",italic(" (Dist"),", log)")),
       y="Projected aridity change (RCP 8.5)",color="",fill="")


ggsave("./Figures/SI/Sensitivity_clusters_climatic_data.pdf",p1,width = 8,height = 4)



p1=ggplot(d_clim_averaged)+geom_line(aes(x=RCP,y=trend,group=Site_ID))+the_theme2+
  labs(x="RCP scenario",y="Temporal trend in aridity")

p2=ggplot(NULL)+geom_point(aes(x=d_clim_averaged$trend[which(d_clim_averaged$RCP=="rcp45")],
                               y=d_clim_averaged$trend[which(d_clim_averaged$RCP=="rcp85")]))+the_theme2+
  labs(x="Temporal trend in aridity with the RCP 4.5 scenario",y="Temporal trend in aridity with the RCP 8.5 scenario")

# summarizing information in each site
d_summarized=d%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_dis50_log=log(quantile(pinfer-pcrit,na.rm = T,.5)),
                   abs_dis50=(quantile(pinfer-pcrit,na.rm = T,.5)))%>%
  add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[keep_sites],Facilitation=d_biocom$Facilitation[keep_sites])

d_summarized=d_summarized%>%
  add_column(., 
             Proj_aridity45=d_clim_averaged$trend[which(d_clim_averaged$RCP=="rcp45")],
             Proj_aridity85=d_clim_averaged$trend[which(d_clim_averaged$RCP=="rcp85")])
set.seed(123)
#kmeans with 5 clusters to better interpretation of the clusters
kmean_sites = kmeans(scale(d_summarized[,c("abs_dis50_log","Proj_aridity45")]), 5)
kmean_sites2 = kmeans(scale(d_summarized[,c("abs_dis50_log","Proj_aridity85")]), 5)
d_summarized=d_summarized%>%add_column(., cluster_id=as.character(kmean_sites$cluster))

p3=ggplot(d_summarized)+
  geom_point(aes(x=abs_dis50_log,y=Proj_aridity45,fill=cluster_id),
             color="black")+
  scale_fill_manual(values=c("#FFB77C","#FF707B","#BC8DFF","#FDE7BB","#9CECE5"))+
  annotate("text",
           x=c(-4.8,-4.8,-3,-2.5),
           y=c(9e-4,1e-4,9e-4,1e-4),
           label=c("High risk","Ecological risk","Climatic risk","Low risk"),
           color=c("#FF707B","#BC8DFF","#FFB77C","#41D8B9"),family="NewCenturySchoolbook")+
  the_theme+
  theme(legend.position="none")+
  labs(x=expression(paste("Distance to the desertification point",italic(" (Dist"),", log)")),
       y="Projected aridity change (RCP 8.5)")

p4=ggplot(melt(table(tibble(RCP45=kmean_sites$cluster,
                            RCP85=kmean_sites2$cluster)))%>%
            mutate(., RCP45=recode_factor(RCP45,"1"="Climatic risk, \n medium ecological risk",
                                        "2"= "High risk",
                                        "3"="Ecological risk",
                                        '4'='Climatic risk, \n low ecological risk',
                                        "5"="Low risk"))%>%
            mutate(., RCP85=recode_factor(RCP85,"4"="Climatic risk, \n medium ecological risk",
                                         "2"= "High risk",
                                         "1"="Ecological risk",
                                         '5'='Climatic risk, \n low ecological risk',
                                         "3"="Low risk")))+
  geom_tile(aes(x=RCP45,RCP85,fill=value))+the_theme+scale_fill_viridis_c(option = "A")+
  labs(x="Vulnerability groups using the \n RCP 4.5 scenario",
       y="Vulnerability groups using the \n RCP 8.5 scenario",
       fill="")+guides(color="none",fill="none")+
  geom_text(aes(x=RCP45,y=RCP85,label=value,color=value>30))+
  scale_color_manual(values=c("white","black"))


ggsave("./Figures/SI/Scenario_RCP_trends_aridity.pdf",ggarrange(ggarrange(p1,p2,ncol=2,labels = letters[1:2]),
                                                                ggarrange(p3,p4,ncol=2,labels = letters[3:4]),nrow=2),width = 11,height = 9)

# then getting an index of similarity

d_similarity=tibble()

for (rcp_scena in 1:2){
  for (model_id1 in unique(d_clim$model)[(1+(rcp_scena-1)*9):(9*rcp_scena)]){
    
    keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
    d=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")
    clim_trend=d_clim%>%
      dplyr::filter(., model==model_id1)
    
    # summarizing information in each site
    d_summarized=d%>%
      dplyr::group_by(., Site,MF,aridity,Sand)%>%
      dplyr::summarise(., .groups = "keep",
                       abs_dis50_log=log(quantile(pinfer-pcrit,na.rm = T,.5)),
                       abs_dis50=(quantile(pinfer-pcrit,na.rm = T,.5)))%>%
      add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[keep_sites])
    
    d_summarized=d_summarized%>%
      add_column(., Proj_aridity=clim_trend$trend)
    
    set.seed(123)
    #kmeans with 5 clusters to better interpretation of the clusters
    kmean_sites1 = kmeans(scale(d_summarized[,c("abs_dis50_log","Proj_aridity")]), 5)
    
    
    for (model_id2 in unique(d_clim$model)[(1+(rcp_scena-1)*9):(9*rcp_scena)]){
      
      if (model_id2!=model_id1){
        
        keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
        d=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")
        clim_trend=d_clim%>%
          dplyr::filter(., model==model_id2)
        
        # summarizing information in each site
        d_summarized=d%>%
          dplyr::group_by(., Site,MF,aridity,Sand)%>%
          dplyr::summarise(., .groups = "keep",
                           abs_dis50_log=log(quantile(pinfer-pcrit,na.rm = T,.5)),
                           abs_dis50=(quantile(pinfer-pcrit,na.rm = T,.5)))%>%
          add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[keep_sites])
        
        d_summarized=d_summarized%>%
          add_column(., Proj_aridity=clim_trend$trend)
        
        set.seed(123)
        #kmeans with 5 clusters to better interpretation of the clusters
        kmean_sites2 = kmeans(scale(d_summarized[,c("abs_dis50_log","Proj_aridity")]), 5)
        
        d_similarity=rbind(d_similarity,tibble(Model1=model_id1,
                                               Model2=model_id2,
                                               RCP_scena=unique(clim_trend$RCP),
                                               Similarity=adjustedRandIndex(kmean_sites1$cluster,kmean_sites2$cluster)))
        
      }
    }
  }
}

p2=ggplot(d_similarity)+
  geom_boxplot(aes(y=Similarity,x=RCP_scena))+
  the_theme+
  labs(x="RCP scenarios",y="Similarity of the K-means partitions")+
  geom_hline(yintercept = c(0,1))


ggsave("./Figures/SI/Sensitivity_clusters_climatic_data.pdf",ggarrange(p1,p2,nrow=2,labels = letters[1:2]),width = 8,height = 8)




#Same but at the site level: does site classification changes with the climatic model used

d_center2=d_center%>%dplyr::filter(.,Average_models=="Yes")
d_single=as.data.frame(table(d_center2[,c(1,3)]))%>%
  dplyr::filter(., Freq>0)%>%
  mutate(., value=recode_factor(value,
                                "1"="Climatic risk, \n medium ecological risk",
                                "2"= "High risk",
                                "3"="Ecological risk",
                                '4'='Climatic risk, \n low ecological risk',
                                "5"="Low risk"))%>%
  arrange(., value)%>%
  mutate(., ID=fct_reorder(ID, desc(value)))



d_center=d_center%>%dplyr::filter(.,RCP=="rcp85")
ggplot(as.data.frame(table(d_center[,c(1,3)]))%>%
         mutate(., value=recode_factor(value,
                                       "1"="Climatic risk, \n medium ecological risk",
                                       "2"= "High risk",
                                       "3"="Ecological risk",
                                       '4'='Climatic risk, \n low ecological risk',
                                       "5"="Low risk"))%>%
         add_column(., Pos=sapply(1:nrow(.),function(x){
           return(d_single$Pos[which(d_single$ID==.$ID[x])])
         }))%>%
         mutate(., ID=fct_reorder(ID, Pos)))+
         geom_bar(aes(x=ID,y=Freq,fill=value),position = "fill",stat = "identity")+
  scale_fill_manual(values=c("#FFB77C","#FF707B","#BC8DFF","#FDE7BB","#9CECE5"))


# >> 17) Mean-field model ----
# 
# d=read.table("./Data/posterior_param.csv",sep=";",header=T)
# keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
# 
# d=tibble(Site=1:345,mean_p=apply(d[,1:345],2,mean),sd_p=apply(d[,1:345],2,sd),median_p=apply(d[,1:345],2,median),
#          mean_q=apply(d[,346:690],2,mean),sd_q=apply(d[,346:690],2,sd),median_q=apply(d[,346:690],2,median),
#          median_eta=apply(d[,691:1035],2,median),
#          Cover=d_biocom$Cover)%>%
#   dplyr::filter(., Site %in% keep_sites)%>%
#   dplyr::select(.,-Site)
# 
# d2=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")%>%
#   dplyr::group_by(., Site,MF,aridity,Sand)%>%
#   dplyr::summarise(., .groups = "keep",
#                    abs_median=median(pinfer-pcrit,na.rm = T),
#                    relativ_median=median((pinfer-pcrit)/pcrit,na.rm = T))%>%
#   dplyr::filter(., Site %in% keep_sites)
# 
# d=cbind(d,d2)
d_MF=read.table("./Data/Mean_field_data.csv",sep=";")

p=ggplot(NULL)+
  geom_tile(data=d_MF%>%filter(.,Branch!="backward"),
            aes(x=p,y=q,fill=Cover_bis))+
  the_theme+
  scale_fill_viridis_c(na.value = "#D8C0A5",option = "G")+
  labs(fill="Vegetation cover",x="Parameter p",
       y="Parameter q",color="Distance to the \n tipping point (Dist)  ")+
  xlim(.2,1.01)+
  annotate("text",x=.35,y=.4,label="No vegetation")+
  geom_hline(yintercept = unique(d_MF$q)[c(1,100,140)])+
  annotate("text",x=1.01,y=unique(d_MF$q)[1]+.025,label="d")+
  annotate("text",x=1.01,y=unique(d_MF$q)[100]+.025,label="c")+
  annotate("text",x=1.01,y=unique(d_MF$q)[140]+.025,label="b")
  
id=1
for (k in unique(d_MF$q)[c(1,100,140)]){
  
  assign(paste0("p_",id),
         ggplot(d_MF%>%filter(., q==k))+
           geom_point(aes(x=p,y=Cover),shape=21,fill="white",color="black")+
           the_theme+
           labs(y="Vegetation cover",x="Parameter p")
         )
  
  id=id+1
}

p_bifu=ggarrange(p_3,p_2,p_1,nrow=3,labels = letters[2:4])

ggsave("./Figures/SI/Phase_diagram.pdf",
       ggarrange(p,p_bifu,ncol=2,widths = c(1,.4)),
       width = 9,height = 6)

# >> 18) Sensitivity to spatial resolution ----

#Change in estimated parameters with coarse graining

spatial_res=c("/2","/3","No change")
list_f=list.files("./Data/Resolution/","posterior")
d_param=tibble()
for (k in list_f){
  d=read.table(paste0("./Data/Resolution/",k),sep=";")
  d_param=rbind(d_param,tibble(median_p=apply(d[,1:345],2,median),
                               median_q=apply(d[,346:690],2,median),
                               median_eta=apply(d[,691:1035],2,median),
                               ID=1:345,
                               Resolution=spatial_res[as.numeric(gsub(".csv","",strsplit(k,"_")[[1]][3]))]))
}

d_param=d_param%>%
  filter(.,ID %in% keep_sites)%>%
  mutate(., Resolution=recode_factor(Resolution,"/2"="Resolution /2","/3"="Resolution /3"))%>%
  mutate(., Resolution=factor(Resolution, 
                              level = c('No change', 'Resolution /2', 'Resolution /3')))%>%
  melt(.,id.vars = c("Resolution","ID"))%>%
  mutate(., variable=recode_factor(variable,
                                   "median_p"="Median of parameter p",
                                   "median_q"="Median of parameter q",
                                   "median_eta"="Median of the scale parameter"))
  


p_param=ggplot(d_param)+
  geom_line(aes(x=Resolution,y=value,group=ID),lwd=.1,color="grey",alpha=.2)+
  geom_point(data=d_param%>%
               dplyr::group_by(., Resolution,variable)%>%
               dplyr::summarise(., .groups = "keep",mean_value=mean(value,na.rm = T)),
             aes(x=Resolution,y=mean_value),color="red",fill="white",shape=21,size=2.5)+
  geom_line(data=d_param%>%
              dplyr::group_by(., Resolution,variable)%>%
              dplyr::summarise(., .groups = "keep",mean_value=mean(value,na.rm = T)),
            aes(x=Resolution,y=mean_value,group=interaction(variable)),lwd=1,color="red")+
  the_theme+
  facet_wrap(.~variable,scales="free")+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="",fill="")


#Change in distance between observation and closest simulations with coarse graining
spatial_res=c("/2","/3","No change")
list_f=list.files("./Data/Resolution/","sumstat")
d_sumstat=tibble()
for (k in list_f){
  d=read.table(paste0("./Data/Resolution/",k),sep=";")
  d_sumstat=rbind(d_sumstat,d%>%
                    add_column(.,Resolution=spatial_res[as.numeric(gsub(".csv","",
                                                                        strsplit(k,"_")[[1]][3]))],
                               ID=1:345))
}
colnames(d_sumstat)[1:11]=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
                            "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")
d_sumstat=d_sumstat%>%
  mutate(., Resolution=recode_factor(Resolution,"/2"="Resolution /2","/3"="Resolution /3"))%>%
  mutate(., Resolution=factor(Resolution, 
                              level = c('No change', 'Resolution /2', 'Resolution /3')))
mean_rmse=d_sumstat%>%
  melt(., id.vars=c("Resolution","ID"))%>%
  dplyr::group_by(.,variable,Resolution)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value,na.rm=T))

p_nrmse=d_sumstat%>%
  filter(.,ID %in% keep_sites)%>%
  melt(.,id.vars = c("Resolution","ID"))%>%
  ggplot(.)+
  geom_jitter(aes(x=Resolution,y=value,color=Resolution),width = .1,alpha=.4)+
  geom_point(data=mean_rmse,aes(x=Resolution,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  the_theme+
  facet_wrap(.~variable,scales="free")+
  labs(color="",y="NRMSE")+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())


p_tot=ggarrange(p_nrmse,p_param,nrow=2,labels = LETTERS[1:2],heights = c(1,.5))

ggsave("./Figures/SI/Coarse_graining_observations.pdf",p_tot,width = 8,height = 9)

# >> 19) Comparing the estimations when varying p with when varying p & q ----

keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
d_with_p=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_median=median(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pcrit,na.rm = T),
                   relativ_median=median((pinfer-pcrit)/pcrit,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pcrit,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  dplyr::filter(., Site %in% keep_sites)


d_with_pq=read.table("./Data/Resilience_metrics_with_p_q.csv",sep=";")%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_median=median(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pcrit,na.rm = T),
                   relativ_median=median((pinfer-pcrit)/pcrit,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pcrit,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  dplyr::filter(., Site %in% keep_sites)


p1=ggplot(NULL)+
  geom_point(aes(x=d_with_p$abs_median,d_with_pq$abs_median),shape=21,color="black",fill="#E0B6FF")+
  the_theme+
  labs(x="Distance to the desertification point (along p)",y="Distance to the desertification point (along p & q)")+
  geom_abline(slope = 1)+
  annotate("text",x=.1,y=.2,label=paste0("Spearman = ",round(cor(d_with_p$abs_median,d_with_pq$abs_median,method="spearman"),3)))

p2=ggplot(NULL)+
  geom_point(aes(x=log(d_with_p$abs_median),log(d_with_pq$abs_median)),shape=21,color="black",fill="#E0B6FF")+
  the_theme+
  labs(x="Distance to the desertification point (along p, log)",y="Distance to the desertification point (along p & q, log)")+
  geom_abline(slope = 1)+
  annotate("text",x=-4,y=-2,label=paste0("Spearman = ",round(cor(log(d_with_p$abs_median),log(d_with_pq$abs_median),method="spearman"),3)))

p_cor=ggarrange(p1,p2,ncol=2)


boot_function_lm = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}

d=read.table("./Data/posterior_param.csv",sep=";",header=T)
keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
d=tibble(Site=1:345,mean_p=apply(d[,1:345],2,mean),sd_p=apply(d[,1:345],2,sd),median_p=apply(d[,1:345],2,median),
         mean_q=apply(d[,346:690],2,mean),sd_q=apply(d[,346:690],2,sd),median_q=apply(d[,346:690],2,median),
         Plot_n=d_biocom$Plot_n,
         Aridity=d_biocom$Aridity,
         Sand=d_biocom$Sand,
         MF=d_biocom$MF,
         Soil_A=d_biocom$Soil_A,
         Slope=d_biocom$Slope,
         Facilitation=d_biocom$Facilitation,
         SR=d_biocom$SR,
         Long_cos=d_biocom$Long_cos,
         Long=d_biocom$Longitude,
         Long_sin=d_biocom$Long_sin,
         Lat=d_biocom$Lat,
         Elevation=d_biocom$Elevation,
         Cover=d_biocom$Cover)%>%
  filter(., Site %in% keep_sites)

d2=read.table("./Data/Resilience_metrics_with_p_q.csv",sep=";")%>%
  group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_median=median(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_median=median((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  filter(., Site %in% keep_sites)

d=cbind(d,d2)


d2=tibble(p=logit(d$median_p),
          q=logit(d$median_q),
          Size_tipping=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$Size_mean[x],d$Size_sd[x]))}))[,1],
          abs_dist=scale(log(d$abs_median))[,1],
          rela_dist=scale(log(d$relativ_median))[,1],
          Sand=(d$Sand-mean(d$Sand,na.rm=T))/sd(d$Sand,na.rm = T),
          Site=d$Site,
          MF=(d$MF-mean(d$MF,na.rm=T))/sd(d$MF,na.rm = T),
          SR=(d$SR-mean(d$SR,na.rm=T))/sd(d$SR,na.rm = T),
          Facilitation=(d$Facilitation-mean(d$Facilitation,na.rm=T))/sd(d$Facilitation,na.rm = T),
          Cover=(d$Cover-mean(d$Cover,na.rm=T))/sd(d$Cover,na.rm = T),
          Aridity=(d$Aridity-mean(d$Aridity,na.rm=T))/sd(d$Aridity,na.rm = T),
          Lat=(d$Lat-mean(d$Lat,na.rm=T))/sd(d$Lat,na.rm = T),
          Long_cos=(d$Long_cos-mean(d$Long_cos,na.rm=T))/sd(d$Long_cos,na.rm = T),
          Long=d$Long,
          Long_sin=(d$Long_sin-mean(d$Long_sin,na.rm=T))/sd(d$Long_sin,na.rm = T),
          Elevation=(d$Elevation-mean(d$Elevation,na.rm=T))/sd(d$Elevation,na.rm = T),
          Slope=(d$Slope-mean(d$Slope,na.rm=T))/sd(d$Slope,na.rm = T),
          Plot_n=d$Plot_n)

mod_predictors=gsub("\n     ","","Aridity + MF + Sand + 
                    Lat + Long_cos + Long_sin + Slope + Elevation + ( 1 | Plot_n)")
d_mod=list(Boot_effects=tibble(),Partial_res_data=tibble())
d_info_model=list(global_R2=tibble(),R2_partial_res=tibble(),Vif=tibble(),Moran=tibble(),Effects=tibble())

for (response_var in c("abs_dist","rela_dist","q","p")){
  
  model_abs=(lmer(formula = paste(response_var," ~ ",mod_predictors), data = d2)) #fitting the model
  
  mcp.fnc(model_abs) #checking model assumptions
  
  #ARIDITY
  
  resid_mod=visreg::visreg(fit = model_abs,xvar="Aridity",plot=F) 
  boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Aridity)
  
  d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                    tibble(Driver_name="Aridity",
                                           Response=response_var,
                                           R2=rsq(lm(data=resid_mod$res,visregRes~Aridity)),
                                           pval=get_bootstrapped_pval(boot_AI$t[,1])))
  
  d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                               tibble(Driver_value=resid_mod$res$Aridity,
                                      Driver_name="Aridity",
                                      Response=response_var,
                                      Resids=resid_mod$res$visregRes))
  
  d_mod$Boot_effects=rbind(d_mod$Boot_effects,tibble(Slopes=boot_AI$t[,1],Driver_name="Aridity",Response=response_var))
  
  
  # Multifunctionality 
  
  resid_mod=visreg::visreg(fit = model_abs,xvar="MF",plot=F) 
  boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~MF)
  
  d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                    tibble(Driver_name="Multifunctionality",
                                           Response=response_var,
                                           R2=rsq(lm(data=resid_mod$res,visregRes~MF)),
                                           pval=get_bootstrapped_pval(boot_AI$t[,1])))
  
  d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                               tibble(Driver_value=resid_mod$res$MF,
                                      Driver_name="Multifunctionality",
                                      Response=response_var,
                                      Resids=resid_mod$res$visregRes))
  
  d_mod$Boot_effects=rbind(d_mod$Boot_effects,tibble(Slopes=boot_AI$t[,1],Driver_name="Multifunctionality",Response=response_var))
  
  
  #Properties model
  
  #R2  
  d_info_model$global_R2=rbind(d_info_model$global_R2,tibble(Response=response_var,
                                                             R2m=r.squaredGLMM(model_abs)[1],
                                                             R2c=r.squaredGLMM(model_abs)[2]))
  
  #vif
  d_info_model$Vif=rbind(d_info_model$Vif,tibble(Response=response_var,
                                                 Vif=vif(model_abs),
                                                 Name_pred=names(vif(model_abs))))
  
  #Spatial autocorr
  save=d2
  
  save$Resid	= residuals(model_abs)									
  coords = cbind(x=save$Long, y=save$Lat)
  
  for (sp_scale in c(10,30,50)){
    
    col.knn = knearneigh(coords, k=sp_scale)
    nb.col.knn = knn2nb(col.knn)
    Moran_test=moran.test(save$Resid, nb2listw(nb.col.knn))
    d_info_model$Moran=rbind(d_info_model$Moran,tibble(Response=response_var,
                                                       p_val=Moran_test$p.value,
                                                       Sp_scale=sp_scale))
  }
  
  d_info_model$Effects=rbind(d_info_model$Effects,tibble(Response=response_var,
                                                         Estimate=summary(model_abs)$coefficients[-1,1],
                                                         Std_error=summary(model_abs)$coefficients[-1,2],
                                                         Name_pred=rownames(summary(model_abs)$coefficients)[-1]))
}
save_d_mod=d_mod



d2=d2%>%
  filter(., !is.na(Facilitation))

mod_predictors=gsub("\n     ","","Aridity + MF + Sand + Facilitation +
                    Lat + Long_cos + Long_sin + Slope + Elevation + ( 1 | Plot_n)")

model_abs=(lmer(formula = paste("abs_dist ~ ",mod_predictors), data = d2)) #fitting the model
model_rela=(lmer(formula = paste("rela_dist ~ ",mod_predictors), data = d2)) #fitting the model
model_q=(lmer(formula = paste("q ~ ",mod_predictors), data = d2)) #fitting the model
model_p=(lmer(formula = paste("p ~ ",mod_predictors), data = d2)) #fitting the model
model_cover=(lmer(formula = paste("Cover ~ ",mod_predictors), data = d2)) #fitting the model

d_mod=list(Boot_effects=tibble(),Partial_res_data=tibble())

# Facilitation: distance

resid_mod=visreg::visreg(fit = model_abs,xvar="Facilitation",plot=F) 
boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Facilitation)

d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                  tibble(Driver_name="Facilitation",
                                         Response="abs_dist",
                                         R2=rsq(lm(data=resid_mod$res,visregRes~Facilitation)),
                                         pval=get_bootstrapped_pval(boot_AI$t[,1])))

d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                             tibble(Driver_value=resid_mod$res$Facilitation,
                                    Driver_name="Facilitation",
                                    Response="abs_dist",
                                    Resids=resid_mod$res$visregRes))

# Facilitation: relative distance

resid_mod=visreg::visreg(fit = model_rela,xvar="Facilitation",plot=F) 
boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Facilitation)

d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                  tibble(Driver_name="Facilitation",
                                         Response="rela_dist",
                                         R2=rsq(lm(data=resid_mod$res,visregRes~Facilitation)),
                                         pval=get_bootstrapped_pval(boot_AI$t[,1])))

d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                             tibble(Driver_value=resid_mod$res$Facilitation,
                                    Driver_name="Facilitation",
                                    Response="rela_dist",
                                    Resids=resid_mod$res$visregRes))

# Facilitation: q

resid_mod=visreg::visreg(fit = model_q,xvar="Facilitation",plot=F) 
boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Facilitation)

d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                  tibble(Driver_name="Facilitation",
                                         Response="q",
                                         R2=rsq(lm(data=resid_mod$res,visregRes~Facilitation)),
                                         pval=get_bootstrapped_pval(boot_AI$t[,1])))

d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                             tibble(Driver_value=resid_mod$res$Facilitation,
                                    Driver_name="Facilitation",
                                    Response="q",
                                    Resids=resid_mod$res$visregRes))

# Facilitation: p

resid_mod=visreg::visreg(fit = model_p,xvar="Facilitation",plot=F) 
boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Facilitation)

d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                  tibble(Driver_name="Facilitation",
                                         Response="p",
                                         R2=rsq(lm(data=resid_mod$res,visregRes~Facilitation)),
                                         pval=get_bootstrapped_pval(boot_AI$t[,1])))

d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                             tibble(Driver_value=resid_mod$res$Facilitation,
                                    Driver_name="Facilitation",
                                    Response="p",
                                    Resids=resid_mod$res$visregRes))


d_partial_res=rbind(save_d_mod$Partial_res_data,d_mod$Partial_res_data)


d_boot_CI=rbind(save_d_mod$Boot_effects,d_mod$Boot_effects)

d_boot_CI=d_boot_CI%>%
  dplyr::group_by(., Driver_name,Response)%>%
  dplyr::summarise(., .groups = "keep",q2=median(Slopes),q1=quantile(Slopes,.025),q3=quantile(Slopes,.975),pval=get_bootstrapped_pval(Slopes))


id=1
x_vector=c(-1.5,-1.2,2)
y_vector=c(-1.1,-1.2,-1.15)

id_annotate=1
for (k in c("Multifunctionality","Aridity","Facilitation")){
  
  assign(paste0("p1_",id),
         ggplot(d_partial_res%>%filter(., Driver_name==k,Response=="abs_dist"))+
           geom_point(aes(x=Driver_value,Resids),shape=21,color="grey20",fill="#C9ACDE")+the_theme+
           geom_smooth(aes(x=Driver_value,Resids),method = "lm",se = T,color="black",fill="#C9ACDE")+ 
           annotate("text",x=x_vector[id_annotate],y=y_vector[id_annotate],label=paste0("n = ",nrow(d_partial_res%>%
                                                                                                      filter(., Driver_name==k,Response=="abs_dist"))))+
           labs(x=k,y=title_distance)+
           theme(axis.title = element_text(size=13)))
  id_annotate=id_annotate+1
  id=id+1
}

# p0=ggarrange(ggplot()+theme_void(),p0,ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3))
p_partial=ggarrange(p1_1,p1_2,p1_3,ncol=3,align = "v")

p_tot=ggarrange(p_cor,p_partial,nrow=2,labels=letters[1:2])

ggsave("./Figures/SI/Distance_to_tipping_pq.pdf",p_tot,width = 9,height = 8)


test=MASS::mvrnorm(100,c(0,0,0),matrix(c(c(1,-.35,0),c(-.35,1,.58),c(0,.58,1)),3,3))
X=test[,1]
Y=test[,2]
Z=test[,3]
par(mfrow=c(1,1))
plot(X,Y)
plot(Z,Y)
plot(Z,X)


# >> 20) Sensitivity all results to resilience definition ----


#Figure 3 
d_partial_res=rbind(
  readRDS("./Data/Drivers_stability_metrics_data_uncertainty_with_facilitation.rds")$Partial_res_data,
  readRDS("./Data/Drivers_stability_metrics_data_uncertainty_without_facilitation.rds")$Partial_res_data
)

id=1
for (k in c("Multifunctionality","Facilitation","Aridity")){
  
  assign(paste0("p1_",id),
         ggplot(d_partial_res%>%filter(., Driver_name==k,Response=="rela_dist"))+
           geom_point(aes(x=Driver_value,Resids),shape=21,color="grey20",fill="#DEC8EE")+the_theme+
           geom_smooth(aes(x=Driver_value,Resids),method = "lm",se = T,color="black",fill="#DEC8EE")+
           labs(x=k,y="Relative distance to \n the desertification point")+theme(axis.title = element_text(size=13)))
  
  id=id+1
}

p1=ggarrange(p1_1,p1_2+theme(axis.title.y = element_blank()),
             p1_3+theme(axis.title.y = element_blank()),ncol=3)


keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
d=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")
clim_trend=read.table("./Data/Climatic_data/mean_aridity_trend.csv",sep=";")

# summarizing information in each site
d_summarized=d%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   rela_dis50_log=log(quantile((pinfer-pcrit)/pcrit,na.rm = T,.5)),
                   abs_dis50_log=log(quantile((pinfer-pcrit),na.rm = T,.5)),
                   abs_dis50=(quantile((pinfer-pcrit)/pcrit,na.rm = T,.5)))%>%
  add_column(.,ID=1:nrow(.),Cover=d_biocom$Cover[keep_sites])

d_summarized=d_summarized%>%
  add_column(., Proj_aridity=clim_trend$mean_trend[which(clim_trend$RCP=="rcp85")])

set.seed(123)
#kmeans with 5 clusters to better interpretation of the clusters
kmean_sites  = kmeans(scale(d_summarized[,c("abs_dis50_log","Proj_aridity")]), 5)
kmean_sites2 = kmeans(scale(d_summarized[,c("rela_dis50_log","Proj_aridity")]), 5)

p2=ggplot(melt(table(tibble(abs=kmean_sites$cluster,
                         rela=kmean_sites2$cluster)))%>%
           mutate(., abs=recode_factor(abs,"1"="Climatic risk, \n medium ecological risk",
                                       "2"= "High risk",
                                       "3"="Ecological risk",
                                       '4'='Climatic risk, \n low ecological risk',
                                       "5"="Low risk"))%>%
           mutate(., rela=recode_factor(rela,"4"="Climatic risk, \n medium ecological risk",
                                       "2"= "High risk",
                                       "1"="Ecological risk",
                                       '5'='Climatic risk, \n low ecological risk',
                                       "3"="Low risk")))+
  geom_tile(aes(x=abs,rela,fill=value))+the_theme+scale_fill_viridis_c(option = "A")+
  labs(x="Vulnerability groups using the \n absolute distance to the desertification point",
       y="Vulnerability groups using the \n relative distance to the desertification point",
       fill="")+guides(color="none",fill="none")+
  geom_text(aes(x=abs,y=rela,label=value,color=value>30))+
  scale_color_manual(values=c("white","black"))


p3=ggplot(d_summarized)+
  geom_point(aes(x=abs_dis50_log,y=rela_dis50_log),fill="lightblue")+
  the_theme+
  labs(x="Absolute distance to \n the desertification point",
       y="Relative distance to \n the desertification point")
  

p_tot=ggarrange(p1,ggarrange(p3,p2,ncol=2,labels=letters[2:3],widths = c(1,1.5)),labels = c(letters[1],""),nrow=2)
ggsave("./Figures/SI/Sensitivity_relative_distance.pdf",p_tot,width = 10,height = 8)

# >> 21) Comparing heights of Kefi/Guichard and Eby ----

all_d_guichard=readRDS("./Data/Model_confirmation_Guichard/d_for_figure.rds")
all_d_kefi=readRDS("./Data/Model_confirmation_Kefi//d_for_figure.rds")


p1=ggplot(tibble(Mean_size_Eby=all_d_kefi$d_eby$mean_size_tipping,
                 Size_Kefi = all_d_kefi$d_kefi$size_tipping))+
  geom_point(aes(x=Size_Kefi,Mean_size_Eby),fill="lightblue",shape=21,color="black")+
  the_theme2+
  labs(x="Amount of vegetation loss before \n desertification (in Kefi's model)",
       y="Amount of vegetation loss before \n desertification (in minimal model)")


p2=ggplot(tibble(Mean_size_Eby=all_d_guichard$d_eby$mean_size_tipping,
                 Size_Guichard = all_d_guichard$d_guichard$size_tipping))+
  geom_point(aes(x=Size_Guichard,Mean_size_Eby),fill="lightblue",shape=21,color="black")+
  the_theme2+
  labs(x="Amount of vegetation loss before \n desertification (in Guichard's model)",
       y="Amount of vegetation loss before \n desertification (in minimal model)")

ggsave("./Figures/SI/Predicting_height_vegetation_loss.pdf",ggarrange(p1,p2,ncol=2,labels = letters[1:2]),width = 8,height = 4)


# >> 22) Distance by pairs of sites ----

all_d_guichard=readRDS("./Data/Model_confirmation_Guichard/d_for_figure.rds")
all_d_kefi=readRDS("./Data/Model_confirmation_Kefi//d_for_figure.rds")
for (k in 1:length(all_d_kefi)){assign(names(all_d_kefi)[k],all_d_kefi[[k]])}


for (k in unique(d_spearman$ID_sim)){
  
  corr_sp=filter(d_spearman, Type_dist=="Rela",ID_sim==k)
  
  d_fig=d_eby%>%add_column(., 
                           true_dist_abs=d_kefi$abs_dist,
                           true_size_tipping=d_kefi$size_tipping,
                           true_dist_rela=d_kefi$relativ_dist)%>%
    filter(.,f ==unique(.$f)[k])
  
  assign(paste0("p1_",k),ggplot(NULL)+
           geom_point(aes(x=as.numeric(dist(cbind(d_fig$mean_rela_dist,d_fig$mean_rela_dist))),
                          y=as.numeric(dist(cbind(d_fig$true_dist_rela,d_fig$true_dist_rela)))),
                           color="black",fill="white",shape=23,lwd=.8,size=1)+
           the_theme+
           theme(strip.text.x = element_blank())+
           labs(x="Distance between pairs of sites in the Kefi model \n (using the relative Dist)",
                y="Distance between pairs of sites in the Eby model \n (using the relative Dist)"))
  
  
  
  corr_sp=filter(d_spearman, Type_dist=="Abs",ID_sim==k)
  
  d_fig=d_eby%>%add_column(., 
                           true_dist_abs=d_kefi$abs_dist,
                           true_size_tipping=d_kefi$size_tipping,
                           true_dist_rela=d_kefi$relativ_dist)%>%
    filter(.,f ==unique(.$f)[k])
  
  assign(paste0("p2_",k),ggplot(NULL)+
           geom_point(aes(x=as.numeric(dist(cbind(d_fig$mean_rela_dist,d_fig$mean_rela_dist))),
                          y=as.numeric(dist(cbind(d_fig$true_dist_rela,d_fig$true_dist_rela)))),
                      color="black",fill="white",shape=23,lwd=.8,size=1)+
           the_theme+
           theme(strip.text.x = element_blank())+
           labs(x="Distance between pairs of sites in the Kefi model \n (using the absolute Dist)",
                y="Distance between pairs of sites in the Eby model \n (using the absolute Dist)"))
  
  
}

p_tot=ggarrange(ggarrange(p2_1,p2_2,p2_3,ncol=3),
                ggarrange(p1_1,p1_2,p1_3,ncol=3),nrow=2,labels = LETTERS[1:2],align = "hv")


ggsave("./Figures/SI/Pair_wise_distance_sites_Kefi.pdf",p_tot,width = 12,height = 8)

# >> 24) Rankings using spatial indic in kefi and guichard models ----

all_d_kefi=readRDS("./Data/Model_confirmation_Kefi//d_for_figure.rds")
d_kefi=read.table("./Data/Model_confirmation_Kefi/Stats_kefi.csv",sep=",")[,-c(1:7)]
colnames(d_kefi)=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
                   "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")

d_kefi=d_kefi%>%
  add_column(., ID=1:60)%>%
  dplyr::filter(.,ID%in%all_d_kefi$d_eby$Site)%>%
  add_column(., 
             "Dist in Kefi model"=all_d_kefi$d_kefi$abs_dist,
             "Estimated Dist"=all_d_kefi$d_eby$mean_abs_dist)


all_ranks_kefi=lapply(c("Variance","Autocorrelation","PLR","Exponent p.l.","Dist in Kefi model",
                   "Estimated Dist"),function(x){
                     rank_x=tibble(Rank=(d_kefi[,x]))
                     colnames(rank_x)=x
                     return(rank_x)
                   })%>%bind_cols(.)

all_d_kefi=readRDS("./Data/Model_confirmation_Guichard/d_for_figure.rds")
d_kefi=read.table("./Data/Model_confirmation_Guichard/Stats_guichard.csv",sep=",")[,-c(1:3,15:18)]
colnames(d_kefi)=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
                   "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")

d_kefi=d_kefi%>%
  add_column(., ID=1:60)%>%
  dplyr::filter(.,ID%in%all_d_kefi$d_eby$Site)%>%
  add_column(., 
             "Dist in Guichard model"=all_d_kefi$d_guichard$abs_dist,
             "Estimated Dist"=all_d_kefi$d_eby$mean_abs_dist)



all_ranks_guichard=lapply(c("Variance","Autocorrelation","PLR","Exponent p.l.","Dist in Guichard model",
                   "Estimated Dist"),function(x){
                     rank_x=tibble(Rank=(d_kefi[,x]))
                     colnames(rank_x)=x
                     return(rank_x)
                   })%>%bind_cols(.)



d=read.table("./Data/data_sites.csv",sep=";")
keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1

d2=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_median=median(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pcrit,na.rm = T),
                   relativ_median=median((pinfer-pcrit)/pcrit,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pcrit,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  dplyr::filter(., Site %in% keep_sites)

d=cbind(d[keep_sites,],d2)%>%
  mutate(., PL_expo=PL_expo,abs_mean=abs_mean,cv_psd=cv_psd,fmax_psd=fmax_psd) #to have same trends of variables

all_ranks=lapply(c("moran_I","Spectral_ratio","PL_expo","fmax_psd","abs_mean"),function(x){
  rank_x=tibble(Rank=(d[,x]))
  colnames(rank_x)=x
  return(rank_x)
})%>%bind_cols(.)

colnames(all_ranks)=c("Autocorrelation","SDR","Frac. max","Exponent PL fit","Dist. to desert.")



cor_rank_data=cor(all_ranks,method = "spearman",use = "na.or.complete")
diag(cor_rank_data)=NA
cor_rank_data[upper.tri(cor_rank_data)]=NA

cor_rank_kefi=cor(all_ranks_kefi,method = "spearman")
diag(cor_rank_kefi)=NA
cor_rank_kefi[upper.tri(cor_rank_kefi)]=NA


cor_rank_guichard=cor(all_ranks_guichard,method = "spearman")
diag(cor_rank_guichard)=NA
cor_rank_guichard[upper.tri(cor_rank_guichard)]=NA

p1=ggplot(melt(cor_rank_kefi))+
  geom_tile(aes(x=Var1,Var2,fill=round(abs(value),2)))+
  geom_text(aes(x=Var1,Var2,label=round(abs(value),2)))+
  scale_fill_gradient2()+
  the_theme2+
  labs(x="",y="",fill="Spearman correlation between metrics")+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))

p2=ggplot(melt(cor_rank_guichard))+
  geom_tile(aes(x=Var1,Var2,fill=round(abs(value),2)))+
  geom_text(aes(x=Var1,Var2,label=round(abs(value),2)))+
  scale_fill_gradient2()+
  the_theme2+
  labs(x="",y="",fill="Spearman correlation between metrics")+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))

p3=ggplot(melt(cor_rank_data))+
  geom_tile(aes(x=Var1,Var2,fill=round(abs(value),2)))+
  geom_text(aes(x=Var1,Var2,label=round(abs(value),2)))+
  scale_fill_gradient2(na.value = "white")+
  the_theme2+
  labs(x="",y="",fill="Spearman correlation between metrics")+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))

ggsave("./Figures/SI/Performance_EWS_Dist_models_data.pdf",
       ggarrange(p1+ggtitle("Using Kéfi model"),
                 p2+ggtitle("Using Guichard model"),
                 p3+ggtitle("Using observed ecosystem landscapes"),
                 nrow=3,common.legend = T,legend = "bottom",labels = letters[1:3]),width = 6,height = 14)


rm(list=ls())
source("./ABC_drylands_function.R")

# ---------------------------- Main figures ------------------------------

## >> Fig. 2 Validation models ----

all_d_kefi=readRDS("./Data/Model_confirmation_Kefi/d_for_figure.rds")
for (k in 1:length(all_d_kefi)){assign(names(all_d_kefi)[k],all_d_kefi[[k]])}

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






## >> Fig. 3 Drivers of the resilience of drylands ----

d_mod2=readRDS("./Data/Drivers_stability_metrics_data_uncertainty_without_facilitation.rds")
id=1
for (k in c("Multifunctionality","Aridity","Sand","Soil amelioration")){
  assign(paste0("p_",id),
         ggplot(d_mod2$Partial_res_data%>%filter(., Driver_name==k,Response_var=="abs_dist"))+
           geom_point(aes(x=Driver_value,Resids),shape=21,color="black",fill="#DEC8EE")+the_theme+
           geom_smooth(aes(x=Driver_value,Resids),method = "lm",se = T,color="black",fill="#DEC8EE")+
           labs(x=k,y="Distance to the tipping point"))
  id=id+1
}

p_val=d_mod2$Boot_effects%>%filter(., Response_var=="abs_dist")%>%
  dplyr::group_by(., Driver_name,Response_var)%>%
  dplyr::summarise(., p_val=get_bootstrapped_pval(Slopes))%>%
  mutate(., p_val=ifelse(.$p_val<.001,"",.$p_val))

p_5=ggplot(NULL)+
  geom_violin(data=d_mod2$Boot_effects%>%filter(., Response_var=="abs_dist"),
              aes(x=Driver_name,y=Slopes),fill="#DEC8EE",width=.5,draw_quantiles =c(.5))+
  geom_text(data=NULL,aes(x=1:4,y=max(d_mod2$Boot_effects$Slopes)+.02,label=p_val$p_val),size=3.5)+
  the_theme+
  geom_hline(yintercept = 0,linetype=9)+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="Bootstrapped slopes")


ggsave("./Figures/Drivers_resilience_drylands.pdf",ggarrange(ggarrange(p_1,p_4,p_3,p_2,nrow=2,ncol=2),
                 ggarrange(ggplot()+theme_void(),p_5,ggplot()+theme_void(),nrow=3,heights = c(.3,1,.3)),
                 ncol = 2,labels = letters[1:2],widths = c(1,.6)),height = 6,width = 9)

# ---------------------------- SI figures ------------------------------

# >> 0) Correlation between predictors ----

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


corr_pred=corr.test(d2[,c(1:3,5,7,9,10,11:17)],use = "na.or.complete",adjust = "none")

colnames(corr_pred$r)=rownames(corr_pred$r)=c("Parameter p","Parameter q", "Dist. to tipping point","Sand",
                                              "Multifunctionality","Soil amelioration","Facilitation","Cover","Aridity",
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










# >> 1) Examples of the three sites ----

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
  geom_point(data=pred,aes(p,y=cover,shape=as.factor(Site),color=color),size=2)+
  the_theme+
  scale_shape_manual(values=c(11,10,8))+
  scale_color_manual(values=c("#60BB59","black"))+
  labs(y="Vegetation cover",x="Distance to desertification point (p-pc)",shape="")+
  guides(color = "none",size="none")+
  theme(legend.position = "none")+
  xlim(-0.02,.15)+
  theme(axis.text = element_text(size=13))

ggsave("./Figures/SI/Example_inference_3_sites.pdf",p1,width = 5,height = 3)


# >> 2) Relation between cover, q and distance to tipping point ----

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
  labs(y="Distance to the tipping point (p-pc)",x="Vegetation cover",fill="Median of q")+
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
  labs(y="Distance to the tipping point (p-pc)/pc",x="Vegetation cover",fill="Median of q")+
  guides(color="none")+ggtitle("Relative distance")+
  theme(legend.position = c(.2, .7),legend.key.size = unit(.5, 'cm'),axis.text = element_text(size=13),
        axis.title = element_text(size=14))

ggsave("./Figures/SI/Relation_cover_q_resilience.pdf",
       ggarrange(p1,p2,nrow=2,labels = letters[1:2]),width = 7,height = 8)

# >> 3) Predictors of stability metrics ----

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
             sd_p=d2$sd_p,sd_q=d2$sd_q)


#ploting the relationships between parameters/cover and predictions

p1_1=ggplot(d_summarized)+
  geom_pointrange(aes(median_q,abs_dis50,
                      ymin=abs_dis25,
                      ymax=abs_dis75),fill="#313131",
                  color="black",alpha=.5,shape=21,size=.7)+
  the_theme+
  labs(x="Posterior median of parameter q",y="Distance to the tipping point (p-pc)")

p1_2=ggplot(d_summarized)+
  geom_pointrange(aes(median_p,abs_dis50,
                      ymin=abs_dis25,
                      ymax=abs_dis75),fill="#313131",
                  color="black",alpha=.5,shape=21,size=.7)+
  the_theme+
  labs(x="Posterior median of parameter p",y="Distance to the tipping point (p-pc)")


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

p_tot=ggarrange(p1_2,p1_1,p2_2,p2_1,nrow=2,ncol=2,labels=c(letters[1],"",letters[2],""))

ggsave("./Figures/SI/Predictors_stability.pdf",p_tot,width = 8,height = 6)

# >> 4) Bootstrapped AIC ----
d_mod=read.table("./Data/Cover_vs_spatial_structure_data_uncertainty.csv",sep=";")%>%
  mutate(.,  Predictor=recode_factor( Predictor,
                                      "Cover + Spatial \n structure"="Cover + parameter q",
                                      "Spatial \n structure"="Parameter q"))

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
  the_theme+theme(legend.position = c(.55,.7))+
  guides(fill = guide_legend(override.aes = list(size = 2.5)))+
  labs(x=paste0("Model boostraped AIC (Distance to the tipping point)"),y="Count",fill="")+
  theme(legend.text = element_text(size=7.5),
        legend.title = element_blank(),
        axis.text = element_text(size=13),
        legend.position = c(.25,.8))

ggsave("./Figures/SI/Bootstraped_AIC_q_cover.pdf",p1,width = 6,height = 3)

# >> 5) Replicating figure 3 with facilitation ----

d_mod2=readRDS("./Data/Drivers_stability_metrics_data_uncertainty_with_facilitation.rds")
id=1
for (k in unique(d_mod2$Partial_res_data$Response)){
  assign(paste0("p_",id),
         ggplot(d_mod2$Partial_res_data%>%filter(., Response==k))+
           geom_point(aes(x=Driver_value,Resids),shape=21,color="black",fill="#DEC8EE")+the_theme+
           geom_smooth(aes(x=Driver_value,Resids),method = "lm",se = T,color="black",fill="#DEC8EE")+
           labs(x="Facilitation by nurses",y=k))
  id=id+1
}

p_val=d_mod2$Boot_effects%>%
  dplyr::group_by(., Driver_name,Response)%>%
  dplyr::summarise(., p_val=get_bootstrapped_pval(Slopes))%>%
  mutate(., p_val=ifelse(.$p_val<.001,"",.$p_val))

p_4=ggplot(NULL)+
  geom_violin(data=d_mod2$Boot_effects%>%filter(., Driver_name !="Species richness"),
              aes(x=Response,y=Slopes,group=Response),fill="#DEC8EE",width=.5,draw_quantiles =c(.5))+
  geom_text(data=NULL,aes(x=1:3,y=min(d_mod2$Boot_effects$Slopes)-.02,label=p_val$p_val),size=3.5)+
  the_theme+
  geom_hline(yintercept = 0,linetype=9)+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="Bootstrapped slopes")



p_1=ggarrange(ggplot()+theme_void(),p_1,ggplot()+theme_void(),ncol=3,widths = c(.5,1,.5),labels = "a")
p_4=ggarrange(ggplot()+theme_void(),p_4,ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3),labels = "c")
p_2=ggarrange(p_2,p_3,ncol = 2,labels = "b")
ggsave("./Figures/SI/Drivers_resilience_drylands_with_facilitation.pdf",
       ggarrange(p_1,p_2,p_4,nrow= 3,heights = c(1,1,1.3)),
         height = 10,width = 6)

# >> 6) Replicating figure 3 with relative distance, q or vegetation cover ----

d_mod2=readRDS("./Data/Drivers_stability_metrics_data_uncertainty_without_facilitation.rds")
id=1
for (k in c("Multifunctionality","Aridity","Sand","Soil amelioration")){
  assign(paste0("p_",id),
         ggplot(d_mod2$Partial_res_data%>%filter(., Driver_name==k,Response_var=="rela_dist"))+
           geom_point(aes(x=Driver_value,Resids),shape=21,color="black",fill="#DEC8EE")+the_theme+
           geom_smooth(aes(x=Driver_value,Resids),method = "lm",se = T,color="black",fill="#DEC8EE")+
           labs(x=k,y="Distance to the tipping point (relative)"))
  id=id+1
}

p_val=d_mod2$Boot_effects%>%filter(., Response_var=="rela_dist")%>%
  dplyr::group_by(., Driver_name,Response_var)%>%
  dplyr::summarise(., p_val=get_bootstrapped_pval(Slopes))%>%
  mutate(., p_val=ifelse(.$p_val<.001,"",.$p_val))

p_5=ggplot(NULL)+
  geom_violin(data=d_mod2$Boot_effects%>%filter(., Response_var=="rela_dist"),
              aes(x=Driver_name,y=Slopes),fill="#DEC8EE",width=.5,draw_quantiles =c(.5))+
  geom_text(data=NULL,aes(x=1:4,y=min(d_mod2$Boot_effects$Slopes)-.02,label=p_val$p_val),size=3.5)+
  the_theme+
  geom_hline(yintercept = 0,linetype=9)+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="Bootstrapped slopes")


ggsave("./Figures/SI/Drivers_resilience_drylands_relative_distance.pdf",
       ggarrange(ggarrange(p_1,p_4,p_3,p_2,nrow=2,ncol=2),
                 ggarrange(ggplot()+theme_void(),p_5,ggplot()+theme_void(),nrow=3,heights = c(.3,1,.3)),
                 ncol = 2,labels = letters[1:2],widths = c(1,.6)),height = 6,width = 9)


#With aggregation parameter (q)
d_mod2=readRDS("./Data/Drivers_stability_metrics_data_uncertainty_without_facilitation.rds")
id=1
for (k in c("Multifunctionality","Aridity","Sand","Soil amelioration")){
  assign(paste0("p_",id),
         ggplot(d_mod2$Partial_res_data%>%filter(., Driver_name==k,Response_var=="q"))+
           geom_point(aes(x=Driver_value,Resids),shape=21,color="black",fill="#DEC8EE")+the_theme+
           geom_smooth(aes(x=Driver_value,Resids),method = "lm",se = T,color="black",fill="#DEC8EE")+
           labs(x=k,y="Aggregation parameter (q)"))
  id=id+1
}

p_val=d_mod2$Boot_effects%>%filter(., Response_var=="q")%>%
  dplyr::group_by(., Driver_name,Response_var)%>%
  dplyr::summarise(., p_val=get_bootstrapped_pval(Slopes))%>%
  mutate(., p_val=ifelse(.$p_val<.001,"",.$p_val))

p_5=ggplot(NULL)+
  geom_violin(data=d_mod2$Boot_effects%>%filter(., Response_var=="q"),
              aes(x=Driver_name,y=Slopes),fill="#DEC8EE",width=.5,draw_quantiles =c(.5))+
  geom_text(data=NULL,aes(x=1:4,y=min(d_mod2$Boot_effects$Slopes)-.02,label=p_val$p_val),size=3.5)+
  the_theme+
  geom_hline(yintercept = 0,linetype=9)+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="Bootstrapped slopes")


ggsave("./Figures/SI/Drivers_resilience_drylands_with_q.pdf",
       ggarrange(ggarrange(p_1,p_4,p_3,p_2,nrow=2,ncol=2),
                 ggarrange(ggplot()+theme_void(),p_5,ggplot()+theme_void(),nrow=3,heights = c(.3,1,.3)),
                 ncol = 2,labels = letters[1:2],widths = c(1,.6)),height = 6,width = 9)


#with vegetation cover 

d_mod2=readRDS("./Data/Drivers_stability_metrics_data_uncertainty_without_facilitation.rds")
id=1
for (k in c("Multifunctionality","Aridity","Sand","Soil amelioration")){
  assign(paste0("p_",id),
         ggplot(d_mod2$Partial_res_data%>%filter(., Driver_name==k,Response_var=="Cover"))+
           geom_point(aes(x=Driver_value,Resids),shape=21,color="black",fill="#DEC8EE")+the_theme+
           geom_smooth(aes(x=Driver_value,Resids),method = "lm",se = T,color="black",fill="#DEC8EE")+
           labs(x=k,y="Vegetation cover"))
  id=id+1
}

p_val=d_mod2$Boot_effects%>%filter(., Response_var=="Cover")%>%
  dplyr::group_by(., Driver_name,Response_var)%>%
  dplyr::summarise(., p_val=get_bootstrapped_pval(Slopes))%>%
  mutate(., p_val=ifelse(.$p_val<.001,"",.$p_val))

p_5=ggplot(NULL)+
  geom_violin(data=d_mod2$Boot_effects%>%filter(., Response_var=="Cover"),
              aes(x=Driver_name,y=Slopes),fill="#DEC8EE",width=.5,draw_quantiles =c(.5))+
  geom_text(data=NULL,aes(x=1:4,y=min(d_mod2$Boot_effects$Slopes)-.02,label=p_val$p_val),size=3.5)+
  the_theme+
  geom_hline(yintercept = 0,linetype=9)+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="Bootstrapped slopes")


ggsave("./Figures/SI/Drivers_resilience_drylands_with_cover.pdf",
       ggarrange(ggarrange(p_1,p_4,p_3,p_2,nrow=2,ncol=2),
                 ggarrange(ggplot()+theme_void(),p_5,ggplot()+theme_void(),nrow=3,heights = c(.3,1,.3)),
                 ncol = 2,labels = letters[1:2],widths = c(1,.6)),height = 6,width = 9)


# >> 7) Posteriors: median and bimodality ----

#Bimodality posterior

site=11
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
       width = 12,height = 6)


#Posterior median
d=read.table("./Data/posterior_param.csv",sep=";",header=T)
keep_sites=read.table("./Data/Keeping_sites.csv",sep=";")$V1
d=tibble(Site=keep_sites, 
         Cover=d_biocom$Cover[keep_sites],
         Median_p=apply(d[,keep_sites],2,median),
         Median_q=apply(d[,keep_sites+345],2,median),
         Median_eta=apply(d[,keep_sites+690],2,median))

p=ggarrange(
  ggplot(d)+
    geom_histogram(aes(x=Median_p))+
    the_theme+
    labs(x="Median of p",y="Number of sites"),
  ggplot(d)+
    geom_histogram(aes(x=Median_q))+
    the_theme+
    labs(x="Median of q",y="Number of sites"),
  ggplot(d)+
    geom_histogram(aes(x=Median_eta))+
    the_theme+
    labs(x="Median of the scale parameter",y="Number of sites"),
  ggplot(d)+
    geom_point(aes(x=Median_p,Median_q,color=Cover,fill=Cover))+
    the_theme+
    scale_fill_gradientn(colours = colorRampPalette(c("#D3EFD3","#90C390","#47A747","#126312"))(100))+
    scale_color_gradientn(colours = colorRampPalette(c("#D3EFD3","#90C390","#47A747","#126312"))(100))+
    labs(x="Median of p",y="Median of q"),ncol = 2,nrow=2
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
  labs(x="Longitude",y="Lattitue")

ggsave("./Figures/SI/Map_empirical_sites.pdf",p,width = 6,height = 4)




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

ggsave(paste0("./Figures/SI/Pair_corr_metrics.pdf"),p,width = 12,height = 12)


# PCA on simulations by coloring by the parameter value p, q

d=read.table("./Data/Simulations.csv",sep=";",header=T)%>%
  dplyr::sample_n(., 20000)

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
vPCs = data.frame(cbind(vPC1,vPC2,vPC3))#[c(1,2,3,4,5,6,8,10,11),]
rownames(vPCs) = vlabs
colnames(vPCs) = c("PC1","PC2","PC3")
random_vertical_distrib=runif(nrow(vPCs),-1,1)

save_vPC=vPCs

axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

for (i in 1:3){
  assign(paste0("p1_",i),
         d%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_text(data=vPCs%>%mutate(., PC1=save_vPC[,axes_for_plot$x[i]],PC2=save_vPC[,axes_for_plot$y[i]]),
                     aes(x=PC1*9,y=PC2*9,label=rownames(vPCs)), size=4)+
           geom_segment(data=vPCs%>%mutate(., PC1=save_vPC[,axes_for_plot$x[i]],PC2=save_vPC[,axes_for_plot$y[i]]),
                        aes(x = 0, y = 0, xend = PC1*8, yend = PC2*8+random_vertical_distrib), arrow = arrow(length = unit(1/2, 'picas')), color = "grey30")+
           geom_point(aes(x = PC1, y = PC2, color = p,fill=p),alpha=.5)+
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
           geom_text(data=vPCs%>%mutate(., PC1=save_vPC[,axes_for_plot$x[i]],PC2=save_vPC[,axes_for_plot$y[i]]),
                     aes(x=PC1*7,y=PC2*7,label=rownames(vPCs)), size=4)+
           geom_segment(data=vPCs%>%mutate(., PC1=save_vPC[,axes_for_plot$x[i]],PC2=save_vPC[,axes_for_plot$y[i]]),
                        aes(x = 0, y = 0, xend = PC1*6, yend = PC2*6+random_vertical_distrib), arrow = arrow(length = unit(1/2, 'picas')), color = "grey30")+
           geom_point(aes(x = PC1, y = PC2, color = q,fill=q),alpha=.5)+
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


p_tot=ggarrange(ggarrange(p1_1+ggtitle("Parameter p"),p1_2,p1_3,ncol=3,common.legend = T,legend = "bottom"),
                ggarrange(p2_1+ggtitle("Parameter q"),p2_2,p2_3,ncol=3,common.legend = T,legend = "bottom"),labels = letters[1:2],nrow=2)

ggsave("./Figures/SI/PCA_on_simulations_p_q.pdf",p_tot,width=10,height=7)



# >> 9) Optimizing ABC ----
## Optimization of the ABC method: pre- and post-processing


#Simple rejection algorithm
d=read.table("./Data/NRMSE/RMSE_param_BoxCox_rejection_optim_lambda_yes_N1_1000.csv",sep=";")

mean_rmse_rej=d%>%
  melt(.)%>%
 dplyr::group_by(.,variable)%>%
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

all_sim=expand.grid(N1=c(1000),
                    lambda=c("yes"),
                    Preproc=c("BoxCox","None"),
                    postproc=c("loclinear","neuralnet"))

d=tibble()
for (i in 1:nrow(all_sim)){
  d=rbind(d,read.table(paste0("./Data/NRMSE/RMSE_param_",all_sim$Preproc[i],"_",all_sim$postproc[i],"_optim_lambda_",
                              all_sim$lambda[i],"_N1_",all_sim$N1[i],".csv"),sep=";")%>%
            add_column(., N1=all_sim$N1[i],optim_lambda=all_sim$lambda[i],Post=all_sim$postproc[i],Pre=all_sim$Preproc[i]))
}


mean_rmse=d%>%
  melt(., id.vars=c("N1","optim_lambda","Post","Pre"))%>%
  mutate(., Post=recode_factor(Post,"loclinear"="Linear regression","neuralnet"="Non-linear regression"))%>%
  mutate(., Pre=recode_factor(Pre,"None"="No BoxCox","BoxCox"="Box-Cox"))%>%
  add_column(., Treatment=paste0(.$Pre," & \n ",.$Post))%>%
  dplyr::group_by(.,variable,N1,optim_lambda,Post,Pre,Treatment)%>%
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
  facet_grid(.~Parameter,labeller = label_bquote(cols="Parameter"==.(as.character(Parameter))))+
  the_theme+
  geom_hline(data=tibble(Parameter=c("p","q"),hpos=mean_rmse_rej$mean_rmse[1:2]),aes(yintercept = hpos),color="gray50")+
  theme(strip.text.x = element_text(size=10),axis.text.x = element_text(angle=60,hjust=1))+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values=c("#C46FC5","#80BD5C"))+
  theme(legend.position = "none",axis.title.y = element_blank())

p=ggarrange(p2,p1,ncol=2,widths = c(1,3),labels = LETTERS[1:2])
ggsave(paste0("./Figures/SI/Optimization_inference_preprocessing.pdf"),p,width = 8,height = 6)








## Optimization of the ABC method: neural-network


all_sim=expand.grid(rep_network=seq(10,30,by=10),N_hidden=seq(5,25,by=5))
d=tibble()

for (i in 1:nrow(all_sim)){
  d=rbind(d,read.table(paste0("./Data/NRMSE/RMSE_hidden_preprocessing_NoPLS_",
                              all_sim$N_hidden[i],"_Nnet_",all_sim$rep_network[i],".csv"),sep=";")%>%
            add_column(., N_hidden=all_sim$N_hidden[i],N_rep_net=all_sim$rep_network[i]))
}


mean_rmse=d%>%
  melt(., id.vars=c("N_hidden","N_rep_net"))%>%
 dplyr::group_by(.,variable,N_rep_net,N_hidden)%>%
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

ggsave(paste0("./Figures/SI/Optimization_NN.pdf"),
       p,width = 7,height = 4)






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
list_f=list.files("./Data/Best_sumstat")
all_name=c("All","No PLR","No Exponent p.l.","No PLR & \n Exponent p.l.","No CV PSD","No Frac. max",
           "No CV PSD & \n Frac. max","No CV PSD & \n Exponent p.l. ","No CV PSD, PLR & \n Exponent p.l.")

for (i in 1:length(list_f)){
  d=rbind(d,read.table(paste0("./Data/Best_sumstat/",list_f[i]),sep=";")%>%
            add_column(.,Name=all_name[i]))
}

mean_rmse=d%>%
  melt(., id.vars=c("Name"))%>%
 dplyr::group_by(.,variable,Name)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value,na.rm=T))%>%
  dplyr::rename(., Parameter=variable)

p=ggplot(d%>%
           melt(., id.vars=c("Name")))+
  geom_jitter(aes(x=Name,y=value,color=Name=="All"),
              position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=Name,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="",y="NRMSE",color="")+
  facet_wrap(.~Parameter,labeller = label_bquote(cols="Parameter"==.(as.character(Parameter))))+
  the_theme+
  theme(strip.text.x = element_text(size=10),axis.text.x = element_text(angle=60,hjust=1),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  #scale_color_manual(values=my_pal(9))+
  scale_color_manual(values=c("gray","red"))+
  theme(legend.position = "none")

ggsave(paste0("./Figures/SI/Combination_sumstats.pdf"),p,width = 7,height = 4)





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
           labs(x="Change in spatial resolution",y="NRMSE",color="")+
           the_theme+
           guides(Method="none")+
           scale_color_manual(values=my_pal(5))+
           scale_x_discrete(labels = c("No change","x2","x3","x4","x5"))+
           guides(color="none"))
}

ggsave("./Figures/SI/NMRSE_consistency_inference_param_scale.pdf",
       ggarrange(p_1+ggtitle(TeX("A) Parameter p")),
                 p_2+ggtitle(TeX("B) Parameter q")),
                 p_3+ggtitle(TeX("C) Parameter \\eta")),nrow=3),
       width = 5,height = 7)


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
                nrow=2,labels=letters[1:2])

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
stat_sim=read.table("./Data/Simulations.csv",sep=";",header = T)%>%
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
d_sim=read.table("./Data/Simulations.csv",sep=";",header = T)%>%
  mutate(., Pooling=recode_factor(Pooling,"1"="Model, no change",
                                  "2" = "Model, x2","3" = "Model, x3","4" = "Model, x4","5" = "Model, x5"))

set.seed(123)
d=rbind(stat_sim[,-c(1:2,15)],
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

p21=ggplot(d)+geom_density(aes(p),fill="grey")+the_theme+labs(x="Parameter p",y="Density")+
  ggtitle("Parameter p")
p22=ggplot(d)+geom_density(aes(q),fill="grey")+the_theme+labs(x="Parameter q",y="Density")+
  ggtitle("Parameter q")

p_tot=ggarrange(p1,ggarrange(p21,p22,nrow=2),labels = letters[1:2],widths = c(1,.7))
ggsave("./Figures/SI/Empirical_priors.pdf",p_tot,width = 9,height = 5)

# >> 14) Validating ---- 
## Validation predictions using Kefi dryland vegetation model ----
  
all_d_kefi=readRDS("./Data/Model_confirmation_Kefi/d_for_figure.rds")
for (k in 1:length(all_d_kefi)){assign(names(all_d_kefi)[k],all_d_kefi[[k]])}
param_kefi=read.table("./Data/Model_confirmation_Kefi/Parameters_kefi.csv",sep=";")


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
    labs(x="Distance in the dryland model",y="Distance by the inference approach")+
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
           labs(x="Distance in the dryland model",y="Distance by the inference approach")+
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
                ggarrange(p1_1,p1_2,p1_3,ncol=3),nrow=2,labels = letters[1:2])

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
param_guichard=read.table("./Data/Model_confirmation_Guichard/Parameters_guichard.csv",sep=";")

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
           labs(x="Distance in the mussel model",y="Distance by the inference approach")+
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
           labs(x="Distance in the mussel model",y="Distance by the inference approach")+
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
                ggarrange(p1_1,p1_2,p1_3,p1_4,p1_5,p1_6,ncol=3,nrow=2),nrow=2,labels = letters[1:2])

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


## Validating predictions using null model without spatial structure ----

# comparizon on how much we fitted the statistics from the landscapes without spatial structure

x_y_stat=read.table("./Data/Models_confirmation/Confirm_null/x_y_stat.csv",sep=";")

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
                  bottom = text_grob("Stats in the null vegatation landscapes (virtual observation)",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"))

ggsave("./Figures/SI/Null_xystats.pdf",p,width = 10,height = 8)





d_median=read.table("./Data/Similarity_within_sites_median.csv",sep=";")
d_ABC=read.table("./Data/Similarity_within_sites_ABC_uncertainty.csv",sep=";")

p1=d_median%>%
  melt(., id.vars=c("Plot_n"))%>%
  mutate(., variable=recode_factor(variable,
                                   "NRMSE_p"="Parameter q",
                                   "NRMSE_q"="Parameter p",
                                   "NRMSE_rela"="Relative distance",
                                   "NRMSE_abs"="Absolute distance"))%>%
  ggplot(.)+
  geom_violin(aes(x=variable,y=value,fill=variable),
              color="transparent",alpha=.5,width=.3)+
  geom_boxplot(aes(x=variable,y=value,fill=variable),
               width=.1,alpha=.5,color="black")+
  the_theme+
  labs(x="",y="RMSE within / between sites")+
  geom_hline(yintercept = 1,color="black")+
  scale_color_manual(values=c("#A6D67E","#237F2E","#D0A3E8","#8C26C3"))+
  scale_fill_manual(values=c("#A6D67E","#237F2E","#D0A3E8","#8C26C3"))+
  theme(legend.position = "none")

p2=d_ABC%>%
  dplyr::group_by(.,Plot_n)%>%
  melt(., id.vars=c("Plot_n","ID_rep"))%>%
  mutate(., variable=recode_factor(variable,
                                   "NRMSE_p"="Parameter q",
                                   "NRMSE_q"="Parameter p",
                                   "NRMSE_rela"="Relative distance",
                                   "NRMSE_abs"="Absolute distance"))%>%
  ggplot(.)+
  geom_violin(aes(x=variable,y=value,fill=variable),
              color="transparent",alpha=.5,width=.3)+
  geom_boxplot(aes(x=variable,y=value,fill=variable),
               width=.1,alpha=.5,color="black",outlier.shape = NA)+
  the_theme+
  labs(x="",y="RMSE within / between sites")+
  geom_hline(yintercept = 1,color="black")+
  scale_color_manual(values=c("#A6D67E","#237F2E","#D0A3E8","#8C26C3"))+
  scale_fill_manual(values=c("#A6D67E","#237F2E","#D0A3E8","#8C26C3"))+
  theme(legend.position = "none")

p=ggarrange(
  p1+ggtitle("With posterior median"),
  p2+ggtitle("With posterior uncertainty"),
  labels = letters[1:2],nrow=2
)
ggsave("./Figures/SI/Comparison_within_between_sites.pdf",p,width = 6,height = 6)






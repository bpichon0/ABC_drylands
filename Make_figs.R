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







## >> Observed versus simulated spatial statistics ----


x_y_stat=read.table(paste0("../Data_new/Inferrence/x_y_stat_all.csv"),sep=";")
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
ggsave("../Figures/Final_figs/Inference_stats.pdf",p,width = 10,height = 8)


## >> Prediction distance tipping point ----


d=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")

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
  geom_line(data=pred%>%filter(., abs(q-median(q)) < quantile(abs(q-median(q)),.5),p>.45)%>%
              add_column(., median_or_not=sapply(1:nrow(.),function(x){
                if (.$q[x]==median(.$q)){return("Yes")
                }else{return("No")}
              })),
            aes(p,y=cover,group=ID_sim,color=median_or_not,size=median_or_not))+
  geom_point(data=pinfered%>%filter(., abs(q-median(q)) < quantile(abs(q-median(q)),.5)),
             aes(x=pinfered,y=cover),color=alpha("red",.3))+
  scale_color_manual(values=(c("gray","#A564C3")))+
  scale_size_manual(values=c(.5,1.3))+
  xlim(0.46,.62)+
  the_theme+
  labs(y="Vegetation cover",x="Reproduction parameter (p)",color="Median of q ? ")+
  guides(color = guide_legend(override.aes = list(size = 1.5)),size="none")

landscape_vege=ggplot(Get_empirical_site(site)%>%melt(.))+
  geom_tile(aes(x=Var1,y=Var2,fill=as.factor(value)))+
  scale_fill_manual(values=c("black","white"))+
  theme_transparent()+
  theme(legend.position = "none")

p1=p1 +
  annotation_custom(grob=ggplotGrob(landscape_vege),
                    ymin = 0, ymax=0.27, xmin=.55, xmax=.63)

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
                   Size_tipping75=quantile(Size_tipping,na.rm = T,.75))

p2=ggplot(d_summarized%>%melt(., measure.vars=c("Sand","MF","aridity")))+
  geom_pointrange(aes(value,ymin=relativ_dis25,ymax=relativ_dis75,y=relativ_dis50),shape=21,color="black",fill="#77ADE6")+
  facet_wrap(.~variable,scales = "free")+
  the_theme+
  labs(y="Relative dist. to desert state",x="Driver value")

p3=ggplot(d_summarized%>%melt(., measure.vars=c("Sand","MF","aridity")))+
  geom_pointrange(aes(value,ymin=Size_tipping25,ymax=Size_tipping75,y=Size_tipping50),shape=21,color="black",fill="#77ADE6")+
  facet_wrap(.~variable,scales = "free")+
  the_theme+
  labs(y="Size shift",x="Driver value")

ggsave("../Figures/Final_figs/Dist_tipping.pdf",
       ggarrange(ggarrange(ggplot()+theme_void(),
                           p1,
                           ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3),labels = c("",LETTERS[1],""))
                 ,p2,p3,nrow=3,labels=c("",LETTERS[2:3])),width = 7,height = 10)



#map

world_map <- map_data("world")
p=ggplot(NULL) +
  geom_polygon(data=world_map, aes(x = long, y = lat, group = group),
               fill="lightgray", colour = "white")+
  geom_point(data=d_biocom[unique(d$Site),],aes(x=Longitude,y=Lattitude,color=d_summarized$relativ_dis50),size=3)+
  the_theme+
  scale_color_gradientn(colors = my_pal(4))+
  labs(x="Longitude",y="Lattitue")







# ---------------------------- SI figures ------------------------------

## >> 0) Characteristics of empirical data ----
### Resolution, geographical repartition ----

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




### Density of summary statistics: ecosystem type, type patterns ----

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
           ),Regular2=recode_factor(Regular2,"0"="Irregular","1"="Regular")))+
  geom_density(aes(x=value,fill=as.factor(Regular2)),alpha=.8)+
  the_theme+
  labs(x="Value",y="Density",fill="Type of pattern  ")+
  facet_wrap(.~variable,scales = "free",nrow=3)+
  scale_fill_manual(values=c("#CCB6EA","#BADCA1"))

ggsave("../Figures/Final_figs/SI/Density_empirical_type_patterns.pdf",p+theme(strip.background.x = element_blank()),width = 9,height = 6)







## >> 1) Optimizing ABC ----
### Optimization of the ABC method: pre- and post-processing ----


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


### Optimization of the ABC method: PLS versus no-PLS ----




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








### Optimization of the ABC method: neural-network ----


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






### Number of simulations kept ----

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


### Best summary statistics ----


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


## >> 2) Spatial resolution & system size ---- 
### Change spatial stats with resolution ----

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

### Robustness inference with spatial resolution ----

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
  theme(axis.text.x = element_text(angle=60,hjust=1))+
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
  labs(x="Method ABC (post-processing)",y="Parameter value",color="")+
  the_theme+
  guides(Method="none")

ggsave("../Figures/Final_figs/SI/Consistency_inference_param_scale.pdf",p2,width = 7,height = 4)

### System size and summary statistics ----


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




## >> 3) Comparison data-model : PCA and densities ----


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



## >> 4) ABC-Posteriors ----


# NMRSE for the different combination of summary statistics

list_f=list.files(paste0("../Data_new/Inferrence/"),"NRMSE_sumstat")
d=tibble()
for (k in 1:length(list_f)){
  name_cols=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
                        "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")
  if (k==2){
    d2=read.table(paste0("../Data_new/Inferrence/",list_f[k]),sep=";")
    colnames(d2)=name_cols[-c(10)] #no CV
    
    d2=d2%>%add_column(.,"CV PSD"=NA)%>%
      relocate(., "CV PSD",.before="Frac. max")%>%
      add_column(., Type="No CV PSD")
    
  }else if (k==3){
    
    d2=read.table(paste0("../Data_new/Inferrence/",list_f[k]),sep=";")
    colnames(d2)=name_cols[-c(9)] #No PL
    d2=d2%>%add_column(.,"Exponent p.l."=NA)%>%
      relocate(., "Exponent p.l.",.after="PLR")%>%
      add_column(., Type="No Exponent p.l.")
      
    
  }else if (k==4){ #no pl plr
    
    d2=read.table(paste0("../Data_new/Inferrence/",list_f[k]),sep=";")
    colnames(d2)=name_cols[-c(8,9)]
    d2=d2%>%add_column(.,"PLR"=NA,"Exponent p.l."=NA)%>%
      add_column(., Type="No PLR & Exponent p.l.")%>%
      relocate(., "PLR",.after="SDR")%>%
      relocate(., "Exponent p.l.",.after="PLR")
      
      
  }  else if (k==5){  #no PLR
    
    d2=read.table(paste0("../Data_new/Inferrence/",list_f[k]),sep=";")
    colnames(d2)=name_cols[-c(8)] #no PL
    d2=d2%>%add_column(.,"PLR"=NA)%>%
      add_column(., Type="No PLR")%>%
      relocate(., "PLR",.after="SDR")
    
  }else if (k==6){ #no the 4
    d2=read.table(paste0("../Data_new/Inferrence/",list_f[k]),sep=";")
    colnames(d2)=name_cols[-c(8:11)]
    d2=d2%>%add_column(.,"PLR"=NA,"Exponent p.l."=NA,"CV PSD"=NA,"Frac. max"=NA)%>%
      add_column(., Type="No PLR, Exponent p.l. & \n CV PSD & Frac. max")
    
  }else if (k==7){

    d2=read.table(paste0("../Data_new/Inferrence/",list_f[k]),sep=";")
    colnames(d2)=name_cols[-c(8:10)]    
    d2=d2%>%add_column(.,"PLR"=NA,"Exponent p.l."=NA,"CV PSD"=NA)%>%
      add_column(., Type="No PLR, Exponent p.l. & \n CV PSD")%>%
      relocate(., "Frac. max",.after="CV PSD")
    
  }else {
    d2=read.table(paste0("../Data_new/Inferrence/",list_f[k]),sep=";")
    colnames(d2)=name_cols
    d2=d2%>%add_column(., Type="All")
  }

  
  d=rbind(d,d2)
}


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







# Inference for all combination of summary statistics

list_f=list.files("../Data_new/Inferrence/","x_y")

pdf("../Figures/Final_figs/SI/All_inference.pdf",width = 7,height = 8)
for (k in 1:length(list_f)){
  
  x_y_stat=read.table(paste0("../Data_new/Inferrence/",list_f[k]),sep=";")
  if (k==1) save=x_y_stat
  
  par(mfrow=c(4,3),mar=rep(2,4))
  for (i in 1:(length(colnames(save))-3)){
    if (colnames(save)[i] %in% colnames(x_y_stat)){
      plot(y=filter(x_y_stat,Type=="Sim")[,colnames(save)[i]],xlab="Sim",ylab="Obs",main=colnames(save)[i],
           x=filter(x_y_stat,Type=="Obs")[,colnames(save)[i]],col="gray")
      abline(a=0,b=1)
      
    } else{
      image(matrix(1,1,1,1),col="white")
    }
    
  }
  
  
  x_y_stat=read.table(paste0("../Data_new/Inferrence/",list_f[k]),sep=";")%>%
    filter(., Site_ID %in% which(d_biocom$Nbpixels<80000))
  if (k==1) save=x_y_stat
  
  par(mfrow=c(4,3),mar=rep(2,4))
  for (i in 1:(length(colnames(save))-3)){
    if (colnames(save)[i] %in% colnames(x_y_stat)){
      plot(y=filter(x_y_stat,Type=="Sim")[,colnames(save)[i]],xlab="Sim",ylab="Obs",main=colnames(save)[i],
           x=filter(x_y_stat,Type=="Obs")[,colnames(save)[i]],col="#DC7575")
      abline(a=0,b=1)
      
    } else{
      image(matrix(1,1,1,1),col="white")
    }
    
  }
  
}
dev.off()

##  Correlation parameters and drivers


d=read.table("../Data_new/Inferrence/Posterior_modes_each_sites.csv",sep=";",header=T)[,-1]
colnames(d)=c("Site","p_50","p_25","p_75","q_50","q_25","q_75","Scale")
d=d%>%add_column(., Sand=d_biocom$Sand,Aridity=d_biocom$Aridity,MF=d_biocom$MF)

d_melt=d%>%
  melt(., measure.vars=c("p_50","q_50"),variable.name="Parameter_var")%>%
  dplyr::rename(., Parameter_val=value)%>%
  mutate(., Parameter_var=recode_factor(Parameter_var,"p_50"="Median of p posterior","q_50"="Median of q posterior"))%>%
  melt(., measure.vars=c("MF","Aridity","Sand"),variable.name="Driver_var")%>%
  mutate(., Driver_var=recode_factor(Driver_var,"MF"="Multifunctionality"))%>%
  dplyr::rename(., Driver_val=value)

list_plots=list()
index=1
for (i in unique(d_melt$Parameter_var)){
  for (j in unique(d_melt$Driver_var)){
    list_plots[[index]]=ggplot(d_melt%>%filter(., Driver_var==j,Parameter_var==i))+
      geom_point(aes(x=Driver_val,y=Parameter_val),color="#96C3DC",alpha=.75)+the_theme+
      labs(x=j,y=i)+the_theme
    index=index+1
  }
}

p=annotate_figure(ggarrange(plotlist=list_plots,nrow = 3,ncol = 2),
                  left=text_grob("Median of posterior parameters",rot=90,color="black",size=15,face ="bold",family = "NewCenturySchoolbook"),
                  bottom = text_grob("Drivers",color="black",size=15,face="bold",family = "NewCenturySchoolbook"))
ggsave("../Figures/Final_figs/Correlation_drivers_parameters.pdf",p,width = 10,height = 8)




## >> 5) Influence of the number of simulations kept ----

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






















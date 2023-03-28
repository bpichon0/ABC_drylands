rm(list=ls())
source("./ABC_drylands_function.R")

# ---------------------------- Main figures ------------------------------

## >> PCA levels of aggregation data ----


d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models.csv",sep=";"),
             Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
             PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
  filter(., rho_p>0.03)%>%
  arrange(.,Model)

for (i in 1:(nrow(d_sim)/4)){
  d_sim$PL_expo[(4*(i-1)+1):(4*(i-1)+3)]=d_sim$PL_expo[(4*(i-1)+4)]
}
d_biocom$PL_expo[d_biocom$PL_expo==0]=NA

d_sim_data=rbind(d_sim%>%filter(., Model=="Eby_feedback"),
                 d_biocom[which(d_biocom$Nbpixels<80000),17:ncol(d_biocom)]%>% #keeping only the data with relatively low quality
                   add_column(., Pooling="Data",Model=NA)%>%
                   mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering)))

sumstat_name=colnames(d_sim_data)[1:9]
res.comp=imputePCA(d_sim_data[,which(colnames(d_sim_data) %in% sumstat_name)],ncp=3,scale = T) 

if ("completeObs" %in% names(res.comp)){
  res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
}else {
  res.pca=PCA(res.comp, ncp = 3,  graph=F)
}

axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

for (i in 1:3){
  assign(paste0("p",i),
         d_sim_data%>%
           mutate(., Pooling=recode_factor(Pooling,"1/4"='x4',"1/3"="x3","1/2"="x2","1"="No change"))%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_point(aes(x = PC1, y = PC2, color = Pooling,fill=Pooling,size=Pooling),alpha=.5)+
           scale_size_manual(values=c(rep(.5,4),1.5))+
           scale_color_manual(values=c(brewer_pal(palette = "BrBG")(4),"black"))+
           scale_fill_manual(values=c(brewer_pal(palette = "BrBG")(4),"black"))+
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





# ---------------------------- SI figures ------------------------------

## >> 1) Optimizing ABC ----
### Optimization of the ABC method: pre- and post-processing ----

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

p=ggplot(d%>%melt(., id.vars=c("N1","optim_lambda","Post","Pre"))%>%
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
  
  theme(strip.text.x = element_text(size=10),axis.text.x = element_text(angle=60,hjust=1),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values=c("#C46FC5","#80BD5C"))+
  theme(legend.position = "none")

ggsave(paste0("../Figures/Final_figs/SI/Optimization_inference_preprocessing.pdf"),p,width = 6,height = 7)






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
list_f=list.files("../Data_new/NRMSE/","NA")
for (k in list_f){
  d=rbind(d,read.table(paste0("../Data_new/NRMSE/",k),sep=";")%>%
            add_column(., Nkept=as.numeric(gsub(".csv","",strsplit(k,"_")[[1]][3]))))
}



mean_rmse=d%>%
  melt(., id.vars=c("Nkept"))%>%
  group_by(.,variable,Nkept)%>%
  dplyr::summarise(., .groups = "keep",mean_rmse=mean(value,na.rm=T))%>%
  dplyr::rename(., Parameter=variable)

p=ggplot(d%>%melt(., id.vars=c("Nkept"))%>%dplyr::rename(., Parameter=variable))+
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


ggsave(paste0("../Figures/Final_figs/SI/N_sim_kept.pdf"),p,width = 7,height = 4)



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









## >> 2) Spatial resolution ---- 
### Change spatial stats with resolution ----

stat_sim=read.table("../Data_new/All_new_sim.csv",sep=";")%>%
  mutate(., Id_sim=rep(1:(nrow(.)/5),each=5))%>%
  mutate(., Pooling=as.character(Pooling),rho_p=round(rho_p,5),variable=as.character(variable))%>%
  melt(., id.vars=c("Pooling","Id_sim"))%>%
  filter(., variable %in% c("rho_p","nb_neigh","clustering","skewness","variance","moran_I","Spectral_ratio","cv_psd","fmax_psd"))%>%
  mutate(., variable=recode_factor(variable,
                                   "rho_p"="Cover","nb_neigh"="# neighbors","clustering"= "Clustering","skewness"="Skewness","variance"="Variance",
                                   "moran_I"="Autocorrelation","Spectral_ratio"="SDR","cv_psd"="CV PSD","fmax_psd"="Frac. max"
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

## >> 3) Resolution, densities & characteristics of empirical data ----


d_biocom=read.table("../Data_new/biocom_data.csv",sep=";")

# spatial resolution

p=ggplot(d_biocom)+
  geom_histogram(aes(x=Nbpixels),fill=alpha("blue",.5))+
  the_theme+
  labs(x='# of pixels',y="Count")

ggsave("../Figures/Final_figs/SI/Spatial_resolution_data.pdf",p,width = 6,height = 3)


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


## >> 4) Comparison data-model : PCA and densities ----


# Density data & model
stat_sim=read.table("../Data_new/All_new_sim.csv",sep=";")%>%
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



# PCA model at different observation scale and data 


d=all_d
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
           scale_size_manual(values=c(rep(.7,5),1.5))+
           scale_color_manual(values=c(my_pal(5),"black"))+
           scale_fill_manual(values=c(my_pal(5),"black"))+
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

ggsave(paste0("../Figures/Final_figs/PCA_scale_observation_model_data.pdf"),p,width = 11,height = 5)


## >> 5) ABC-Posteriors ----


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
    colnames(d2)=name_cols[-c(10:11)] #No cv fmax
    d2=d2%>%add_column(.,"CV PSD"=NA,"Frac. max"=NA)%>%
      add_column(., Type="No Frac. max & CV PSD")
      
    
  }else if (k==4){

    d2=read.table(paste0("../Data_new/Inferrence/",list_f[k]),sep=";")
    colnames(d2)=name_cols[-c(9)] #no PL
    d2=d2%>%add_column(.,"Exponent p.l."=NA)%>%
      add_column(., Type="No Exponent p.l.")%>%
      relocate(., "Exponent p.l.",.after="PLR")
    
  }else if (k==5){ #no pl plr
    
    d2=read.table(paste0("../Data_new/Inferrence/",list_f[k]),sep=";")
    colnames(d2)=name_cols[-c(8,9)]
    d2=d2%>%add_column(.,"PLR"=NA,"Exponent p.l."=NA)%>%
      add_column(., Type="No PLR & Exponent p.l.")%>%
      relocate(., "PLR",.after="SDR")%>%
      relocate(., "Exponent p.l.",.after="PLR")
      
      
    
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




# x-y obs-simulated for each spatial statistic

x_y_stat=read.table(paste0("../Data_new/Inferrence/x_y_stat_all.csv"),sep=";")
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
                bottom = text_grob("Observed spatial statistic",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"))
ggsave("../Figures/Final_figs/Inference_stats.pdf",p,width = 10,height = 8)


# 
# #With the ranking of sites
# x_y_stat=read.table(paste0("../Data_new/Inferrence/x_y_stat_all.csv"),sep=";")
# list_plots=list()
# name_plot=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
#             "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")
# 
# for (i in 1:11){
#   d_fil=cbind(filter(x_y_stat%>%
#                        melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
#                 filter(., Type=="Sim")%>%dplyr::rename(., value_sim=value),
#               filter(x_y_stat%>%
#                        melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
#                 filter(., Type=="Obs")%>%dplyr::rename(., value_obs=value)%>%dplyr::select(., value_obs))%>%
#     add_column(., Closest=final_rank$Rank)
#   
#   list_plots[[i]]=ggplot(d_fil)+
#     geom_point(aes(x=value_obs,y=value_sim,color=Closest),alpha=.75)+the_theme+
#     labs(x="",y="")+
#     geom_abline(slope=1,intercept = 0,color="black")+
#     ggtitle(name_plot[i])+
#     theme(title = element_text(size=10))+
#     scale_color_gradientn(colours = viridis(345))
# }
# 
# p=annotate_figure(ggarrange(plotlist=list_plots,ncol = 4,nrow = 3),
#                   left=text_grob("Closest simulations",rot=90,color="black",size=15,face ="bold",vjust=1,family = "NewCenturySchoolbook"),
#                   bottom = text_grob("Observed spatial statistic",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"))
# 





## >> 6) Prediction ----

d=read.table("../Data_new/Inferrence/NRMSE_param_NN_all.csv",sep=";")
d2=read.table("../Data_new/Inferrence/NRMSE_param_rej_all.csv",sep=";")

site_p_NN=d[,1:345]
site_q_NN=d[,346:690]
site_p_rej=d2[,1:345]
site_q_rej=d2[,346:690]



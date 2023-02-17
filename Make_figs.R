rm(list=ls())
source("./ABC_drylands_function.R")



# ---------------------------- SI figures ------------------------------

## >> Optimization of the ABC method: pre- and post-processing ----

all_sim=expand.grid(N1=c(1000,3000),
                    lambda=c("yes"),
                    Preproc=c("BoxCox","None"),
                    postproc=c("loclinear","neuralnet"))

d=tibble()
for (i in 1:nrow(all_sim)){
  d=rbind(d,read.table(paste0("../Data/Step6_Optimizing_inferrence/Pre_post/RMSE_param_",all_sim$Preproc[i],"_",all_sim$postproc[i],"_optim_lambda_",
                              all_sim$lambda[i],"_N1_",all_sim$N1[i],".csv"),sep=";")%>%
            add_column(., N1=all_sim$N1[i],optim_lambda=all_sim$lambda[i],Post=all_sim$postproc[i],Pre=all_sim$Preproc[i]))
}


mean_rmse=d%>%
  melt(., id.vars=c("N1","optim_lambda","Post","Pre"))%>%
  mutate(., Post=recode_factor(Post,"loclinear"="Linear regression","neuralnet"="Non-linear regression"))%>%
  mutate(., Pre=recode_factor(Pre,"None"="No BoxCox","BoxCox"="Box-Cox"))%>%
  add_column(., Treatment=paste0(.$Pre," & \n ",.$Post))%>%
  group_by(.,variable,N1,optim_lambda,Post,Pre,Treatment)%>%
  summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  rename(., "Parameter"="variable")

p=ggplot(d%>%melt(., id.vars=c("N1","optim_lambda","Post","Pre"))%>%
           rename(., "Parameter"="variable")%>%
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
  ylim(0,.5)+
  theme(strip.text.x = element_text(size=10),axis.text.x = element_text(angle=60,hjust=1),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898AA4","#E4C035"))

ggsave(paste0("../Figures/Final_figs/SI/Optimization_inference_preprocessing.pdf"),p,width = 6,height = 5)






## >> Optimization of the ABC method: PLS versus no-PLS ----




d=rbind(read.table(paste0("../Data/Step6_Optimizing_inferrence/Neural_net/RMSE_hidden_preprocessing_PLS_10_Nnet_10.csv"),sep=";")%>%
          add_column(., PLS="Yes"),
        read.table(paste0("../Data/Step6_Optimizing_inferrence/Neural_net/RMSE_hidden_preprocessing_NoPLS_10_Nnet_10.csv"),sep=";")%>%
          add_column(., PLS="No"))



mean_rmse=d%>%
  melt(., id.vars=c("PLS"))%>%
  group_by(.,variable,PLS)%>%
  summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  rename(., "Parameter"="variable")


p=ggplot(d%>%melt(., id.vars=c("PLS"))%>%
           rename(., "Parameter"="variable"))+
  geom_jitter(aes(x=PLS,y=value,color=as.factor(PLS)),
              position = position_jitterdodge(jitter.width = 0.3,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=PLS,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="Using PLS during pre-processing",y="NRMSE",color="")+
  facet_grid(.~Parameter,labeller = label_bquote(cols="Parameter"==.(as.character(Parameter))))+
  the_theme+
  ylim(0,.5)+
  theme(strip.text.x = element_text(size=10),legend.position = "none")+
  scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))

ggsave(paste0("../Figures/Final_figs/SI/Optimization_PLS.pdf"),p,width = 6,height = 3)








## >> Optimization of the ABC method: neural-network ----


all_sim=expand.grid(rep_network=seq(10,30,by=10),N_hidden=seq(5,25,by=5))

d=tibble()

for (i in 1:nrow(all_sim)){
  d=rbind(d,read.table(paste0("../Data/Step6_Optimizing_inferrence/Neural_net/RMSE_hidden_preprocessing_NoPLS_",
                              all_sim$N_hidden[i],"_Nnet_",all_sim$rep_network[i],".csv"),sep=";")%>%
            add_column(., N_hidden=all_sim$N_hidden[i],N_rep_net=all_sim$rep_network[i]))
}


mean_rmse=d%>%
  melt(., id.vars=c("N_hidden","N_rep_net"))%>%
  group_by(.,variable,N_rep_net,N_hidden)%>%
  summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  rename(., "Parameter"="variable")%>%
  mutate(., N_hidden=as.character(N_hidden))


p=ggplot(d%>%melt(., id.vars=c("N_hidden","N_rep_net"))%>%
           rename(., "Parameter"="variable")%>%
           mutate(., N_hidden=as.character(N_hidden)))+
  geom_jitter(aes(x=factor(N_hidden,level=c("5","10",'15',"20",'25')),y=value,color=as.factor(N_hidden)),
              position = position_jitterdodge(jitter.width = 0.3,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=N_hidden,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="Number hidden neurons",y="NRMSE",color="")+
  facet_grid(Parameter~N_rep_net,labeller = label_bquote(rows="Parameter"==.(as.character(Parameter)),cols="# evaluation NN"==.(N_rep_net)))+
  the_theme+
  ylim(0,.5)+
  theme(strip.text.x = element_text(size=10),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"),breaks=c('5', '10', '15',"20","25"))

ggsave(paste0("../Figures/Final_figs/SI/Optimization_NN.pdf"),
       p,width = 7,height = 4)






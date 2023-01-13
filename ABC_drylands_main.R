rm(list=ls())
source("./ABC_drylands_function.R")


# Step 1 : Generating parameters for sensitivity analysis ----

dir.create("../Data/Step1_sensitivity",showWarnings = F)

## One dimensional analysis ----


range_params=data.frame(min = c(0, 0.02, 0, 0.005,0,0,0),
                        max = c(1, 1, 1, 1,1,1,1))
rownames(range_params)=c("r", "d", "f", "m","b","c","delta")

based_param=c(.05,.1,.9,.1,.5,.2,.1)
N_sim_param=40

sensi_table_param=tibble()
for (id_param in 1:nrow(range_params)){
  
  sensi_table_param=rbind(sensi_table_param,
                          expand_grid(based_param[1],based_param[2],based_param[3],based_param[4],
                                      based_param[5],based_param[6],based_param[7],
                                      seq(range_params$min[id_param],range_params$max[id_param],
                                          length.out=N_sim_param)
                                      ))
  
}
colnames(sensi_table_param)=c("r", "d", "f", "m","b","c","delta","range")

for (i in 1:(ncol(sensi_table_param)-1)){
  sensi_table_param[(N_sim_param*(i-1)+1):(i*N_sim_param),i]=sensi_table_param$range[(N_sim_param*(i-1)+1):(i*N_sim_param)]
}
sensi_table_param=sensi_table_param%>%
  select(., -range)

write.table(sensi_table_param,"../Data/Step1_sensitivity/sensitivity_1D_param.csv",sep=";",row.names = F)



#plotting metrics against parameters

results1D=read.table("../Data/Step1_sensitivity/All_data_1D.csv",sep=",")
colnames(results1D)=c("r", "d", "f", "m","b","c","delta",
                      "rho_p","nb_neigh","clustering","skewness","variance","moran_I","Spectral_ratio",
                      "PLR","PL_expo")

n_param=7;nsim=40;n_metric=9

pdf("../Figures/Sensitivity_1D.pdf",width = 6,height = 7)
for (i in 1:n_param){
  
  par(mar = c(4,2,.2,1))
  
  par(mfrow=c(5,2))
  
  d=results1D[((i-1)*nsim+1):(i*nsim),]
  
  col_metric=c("green","blue","red","pink","black","orange","grey","lightblue","purple")
  name_metric=c("Cover","# neighbors","Clustering","Skeness","Variance","Moran I",
                "Spectral ratio","PLR","PL expo")
  
  for (metric in 1:n_metric){
    if (metric==3){ #for clustering we constrain the values as they go up with low cover
      plot(x=d[,i],y=d[,colnames(results1D)[metric+n_param]],col=col_metric[metric],
           xlab="",ylab=name_metric[metric],ylim=c(0,10))
    } else{
      plot(x=d[,i],y=d[,colnames(results1D)[metric+n_param]],col=col_metric[metric],
           xlab="",ylab=name_metric[metric])
    }
    
    legend(x = "topleft",legend=c(name_metric[metric]),fill = c(col_metric[metric]))
  }
  mtext(paste0("Parameter = ",colnames(results1D)[i]),line=-52,outer=TRUE)
  
}
dev.off()






# Getting the trend for each metric and plotting it

results1D=read.table("../Data/Step1_sensitivity/All_data_1D.csv",sep=",")
colnames(results1D)=c("r", "d", "f", "m","b","c","delta",
                      "rho_p","nb_neigh","clustering","skewness","variance","moran_I","Spectral_ratio",
                      "PLR","PL_expo")

n_param=7;n_metric=9;nsim=40

trends_param=trends_with_scaling=tibble()
for (metric in 1:n_metric){
  for (i in 1:n_param){
    
    d=results1D[((i-1)*nsim+1):(i*nsim),c(i,8,metric+n_param)] #we keep the parameter, the cover and the focal metric
    if (any(which(d[,2]<.1))) d=d[-which(d[,2]<.2),] #to avoid problems of linear regression fitting 
    
    reg=lm(d[,3]~d[,1])$coefficients[2]
    trends_param=rbind(trends_param,tibble(Coeff=reg,
                                           Param=as.character(colnames(results1D)[i]),
                                           Metric=as.character(colnames(results1D)[metric+n_param])))
    
    trends_with_scaling=rbind(trends_with_scaling,
                              tibble(Coeff=reg,
                                     Param=as.character(colnames(results1D)[i]),
                                     Metric=as.character(colnames(results1D)[metric+n_param])))
  }
  #for each metric, we normalize the trends across parameters
  trends_with_scaling$Coeff[((metric-1)*n_param+1):(metric*n_param)] = as.numeric(scale(trends_with_scaling$Coeff[((metric-1)*n_param+1):(metric*n_param)]))
}

p1=ggplot(trends_param)+
  geom_tile(aes(x=Param,y=Metric,fill=Coeff),color="black")+
  theme_classic()+
  scale_fill_gradient2(low="#D82828",mid="white",high="#0E5DCE")

p2=ggplot(trends_with_scaling)+
  geom_tile(aes(x=Param,y=Metric,fill=Coeff),color="black")+
  theme_classic()+
  scale_fill_gradient2(low="#D82828",mid="white",high="#0E5DCE")

ggsave("../Figures/Trends_1D_param.pdf",ggarrange(p1,p2,ncol=2,labels=letters[1:2]),
       width = 10,height = 4)



## Two dimensional analysis for parameters pairs ----

range_params=data.frame(min = c(0, 0.02, 0, 0.005,0,0,0),
                        max = c(1, 1, 1, 1,1,1,1))
rownames(range_params)=c("r", "d", "f", "m","b","c","delta")

based_param=c(.05,.1,.9,.1,.5,.2,.1)
N_sim_param=20 #we reduce the number of simulations

sensi_table_param=tibble()
for (id_param1 in 1:(nrow(range_params)-1)){ # for each pair of parameters
  for (id_param2 in (id_param1+1):nrow(range_params)){
  
    sensi_table_param=rbind(sensi_table_param,
                            expand_grid(based_param[1],based_param[2],based_param[3],based_param[4],
                                        based_param[5],based_param[6],based_param[7],
                                        seq(range_params$min[id_param1],range_params$max[id_param1],
                                            length.out=N_sim_param),
                                        seq(range_params$min[id_param2],range_params$max[id_param2],
                                            length.out=N_sim_param)
                            ))
  }
  
}
colnames(sensi_table_param)=c("r", "d", "f", "m","b","c","delta","range1","range2")

index=1
for (id_param1 in 1:(nrow(range_params)-1)){ # for each pair of parameters
  for (id_param2 in (id_param1+1):nrow(range_params)){
    sensi_table_param[((N_sim_param**2)*(index-1)+1):(index*(N_sim_param**2)),id_param1]=
      sensi_table_param$range1[((N_sim_param**2)*(index-1)+1):(index*(N_sim_param**2))]
    
    sensi_table_param[((N_sim_param**2)*(index-1)+1):(index*(N_sim_param**2)),id_param2]=
      sensi_table_param$range2[((N_sim_param**2)*(index-1)+1):(index*(N_sim_param**2))]
    index=index+1
  }
}
sensi_table_param=sensi_table_param%>%
  select(., -range1,-range2)

write.table(sensi_table_param,"../Data/Step1_sensitivity/sensitivity_2D_param.csv",sep=";",row.names = F)

#plotting
results2D=read.table("../Data/Step1_sensitivity/All_data_2D.csv",sep=",")
colnames(results2D)=c("r", "d", "f", "m","b","c","delta",
                      "rho_p","nb_neigh","clustering","skewness","variance","moran_I")

n_param=7;nsim=400
pdf("../Figures/Sensitivity_2D.pdf",width = 10,height = 6)
index=1
for (i in 1:(n_param-1)){
  for (j in (i+1):n_param){
  
    d=results2D[((index-1)*nsim+1):(index*nsim),]
    d$Var1=d[,i]
    d$Var2=d[,j]
      
    for (metric in 1:6){
      which_metric=c("rho_p","nb_neigh","clustering","skewness","variance","moran_I")[metric]
      
      if (which_metric != "rho_p"){
        
        assign(paste0("p_",metric),d%>%
                 melt(., measure.vars=which_metric)%>%
                 mutate(value = ifelse(value>10 | value==0, NA, value))%>%
                 
                 ggplot(.)+
                 geom_tile(aes(x=Var1,y=Var2,
                               fill=value))+
                 theme_classic()+
                 labs(x=colnames(d)[i],y=colnames(d)[j],fill=which_metric)+
                 scale_fill_gradientn(colors = colorRampPalette(c("#F1DB71","#9DE284","#84AEE2"))(100)))
        
      } else{
        
        assign(paste0("p_",metric),d%>%
                 melt(., measure.vars=which_metric)%>%
                 mutate(value = ifelse(value>10 | value==0, NA, value))%>%
                 
                 ggplot(.)+
                 geom_tile(aes(x=Var1,y=Var2,
                               fill=value))+
                 theme_classic()+
                 labs(x=colnames(d)[i],y=colnames(d)[j],fill=which_metric)+
                 scale_fill_gradient2(low = "white",mid = "#6FDE47",
                                                          high = "#218647",midpoint = .5))
        
      }
    }    
    
    print(ggarrange(p_1,p_2,p_3,p_4,p_5,p_6,ncol = 3,nrow = 2 ))
  
    index=index+1
  }
}
dev.off()

#example for overleaf with f and c

index=14
d=results2D[((index-1)*nsim+1):(index*nsim),]
d$Var1=d[,3]
d$Var2=d[,6]


for (metric in 1:6){
  which_metric=c("rho_p","nb_neigh","clustering","skewness","variance","moran_I")[metric]
  
  if (which_metric != "rho_p"){
    
    assign(paste0("p_",metric),d%>%
             melt(., measure.vars=which_metric)%>%
             mutate(value = ifelse(value>10 | value==0, NA, value))%>%
             
             ggplot(.)+
             geom_tile(aes(x=Var1,y=Var2,
                           fill=value))+
             theme_classic()+
             labs(x=colnames(d)[i],y=colnames(d)[j],fill=which_metric)+
             scale_fill_gradientn(colors = colorRampPalette(c("#F1DB71","#9DE284","#84AEE2"))(100)))
    
  } else{
    
    assign(paste0("p_",metric),d%>%
             melt(., measure.vars=which_metric)%>%
             mutate(value = ifelse(value>10 | value==0, NA, value))%>%
             
             ggplot(.)+
             geom_tile(aes(x=Var1,y=Var2,
                           fill=value))+
             theme_classic()+
             labs(x=colnames(d)[i],y=colnames(d)[j],fill=which_metric)+
             scale_fill_gradient2(low = "white",mid = "#6FDE47",
                                  high = "#218647",midpoint = .5))
    
  }
}    

p=ggarrange(p_1,p_2,p_3,p_4,p_5,p_6,ncol = 3,nrow = 2 )
ggsave("../Figures/Example_interaction_f_c.pdf",p,width = 10,height = 6)



# Step 2 : Generating pseudo data-sets for cross verification -----

set.seed(123)
#defining the priors
range_priors=data.frame(min = c(0.02, 0, 0.005,0,0,0),
                        max = c(1, 1, 1,1,1,1))
rownames(range_priors)=c("d", "f", "m","b","c","delta")

# Latin hypercube sampling on the priors
pseudo_param=as.data.frame(Latinhyper(range_priors, 5e4))
pseudo_param$r=.05
pseudo_param[, "m"]=qlnorm(pseudo_param[,"m"], -2.25,.7)


test=Latinhyper(range_priors, 1e3)
test[,"m"]=qlnorm(test[,"m"], -2.25,.7)
pairs(test) #as example

write.table(pseudo_param,'../Data/Pseudo_parameters.csv',sep=";",row.names = F)

























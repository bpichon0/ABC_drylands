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

pdf("../Figures/Sensitivity/Sensitivity_1D.pdf",width = 6,height = 7)
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

ggsave("../Figures/Sensitivity/Trends_1D_param.pdf",ggarrange(p1,p2,ncol=2,labels=letters[1:2]),
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
pdf("../Figures/Sensitivity/Sensitivity_2D.pdf",width = 10,height = 6)
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
ggsave("../Figures/Sensitivity/Example_interaction_f_c.pdf",p,width = 10,height = 6)



# Step 2 : Generating pseudo data-sets for cross verification -----

set.seed(123)
#defining the priors
range_priors=data.frame(min = c(0, 0, 0.005,0,0,0),
                        max = c(1, 1, 1,1,1,1))
rownames(range_priors)=c("d", "f", "m","b","c","delta")

# Latin hypercube sampling on the priors
pseudo_param=as.data.frame(Latinhyper(range_priors, 2.5e5))
pseudo_param$r=.05
pseudo_param[, "m"]=qlnorm(pseudo_param[,"m"], -2.25,.7)

#reorganizing
pseudo_param=pseudo_param[,c(7,1:6)]
write.table(pseudo_param,'../Data/Pseudo_parameters.csv',sep=";",row.names = F)


#Merging the julia outputs

list_simu=list.files('../Data/Step2_cross_validation',pattern = ".csv")

d_all=tibble()
for (file_simu in list_simu){
  
  d=read.table(paste0("../Data/Step2_cross_validation/",file_simu),sep=",")
  colnames(d)= c("r","d", "f", "m","b","c","delta",
                 "rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                 "Spectral_ratio","PLR","PL_expo")
  
  d_all=rbind(d_all,d)
}

d_all=d_all[-which(is.nan(d_all$PLR) | is.nan(d_all$PL_expo)),]

matrix_param=d_all[,1:7]
matrix_sumstat=d_all[,8:(ncol(d_all))]



#cross-validation



N_for_cross_validation = 100

  
for (method_abc in c("rejection","neuralnet","loclinear")){
  
  mat_cor_param=array(0,c(6,6,100)) #correlation matrix for parameters
  
  pdf(paste0("../Figures/Cross_validation/Cross_validation_n",N_for_cross_validation,"_",method_abc,".pdf"),width = 8,height = 4)
  d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
  
  for (n in 1:N_for_cross_validation){
    
    n_cross=ifelse(n==55,n+1,n)
    
    #for each virtual data, we perform ABC rejection algorithm with linear regression adjustment for posterior
    cross_valid=abc(target = matrix_sumstat[n_cross,],
                    param = matrix_param,sumstat = matrix_sumstat,
                    tol = .006,method = method_abc)
    
    if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
    
    
    #Matrix of correlation between parameters
    mat_cor_param[,,n]=cor(cross_valid$adj.values[,-1])
    
    
    #We plot the differences in posterior distribution/true parameter
    d_melt=as.data.frame(cross_valid$adj.values)%>%
      melt(., id.vars="r")
    
    p1=ggplot(d_melt)+
      geom_density(aes(x=value))+
      facet_grid(.~variable,scales = "free_y")+
      theme_classic()+xlim(0,1)+
      theme(legend.position = "bottom")+
      labs(x="Value",y="Density",col="")
    
    
    for (i in colnames(matrix_param)[-1]){
      p1=p1+geom_vline(data=tibble(x=c(matrix_param[n_cross,i],colMeans(cross_valid$adj.values)[i]),col=c('True param',"Mean posterior"),variable=i),
                   aes(xintercept =x ,color=col))
      
    }
      
    print(p1)
    
    d_melt=as.data.frame(cross_valid$ss)%>%
      melt(.)
    
    id=1
    par(mfrow=c(2,5))
    # for (i in colnames(matrix_sumstat)){
    #   assign(paste0("p_",id),ggplot(filter(d_melt,variable==i))+
    #            geom_density(aes(x=value))+
    #            geom_vline(data=tibble(x=c(matrix_sumstat[n_cross,i],colMeans(cross_valid$ss)[i]),col=c('True param',"Mean posterior"),variable=i),
    #                       aes(xintercept =x ,color=col))+
    #            theme_classic()+
    #            theme(axis.title.y = element_blank())+
    #            theme(legend.position = "bottom")+
    #            labs(x="Value",y="Density",col="")+
    #            ggtitle(i))
    #          
    #   id=id+1
    # }
    for (i in colnames(matrix_sumstat)){
      plot(density(d_melt$value[which(d_melt$variable==i)]),main=i,xlab="Value")
      abline(v = matrix_sumstat[n_cross,i],col="blue")
      abline(v = colMeans(cross_valid$ss)[i],col="red")
    }
    
    
    # p2=ggarrange(p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,ncol=9,legend = "none")
    # 
    # print(ggarrange(p1,p2,nrow=2,common.legend = F,heights = c(1,1.2)))
    
    #we save the mean posterior distribution for each and the true observed parameters
    d_cross_param=rbind(d_cross_param,as_tibble(t(colMeans(cross_valid$adj.values)))%>%add_column(., Type="Sim"))
    d_cross_param=rbind(d_cross_param,as_tibble((matrix_param[n_cross,]))%>%add_column(., Type="Obs"))
    
    
    #As we work with virtual data, we do the same for the summary stats we save the mean posterior distribution for each and the true observed parameters
    d_cross_sumstat=rbind(d_cross_sumstat,as_tibble(t(colMeans(cross_valid$ss)))%>%add_column(., Type="Sim"))
    d_cross_sumstat=rbind(d_cross_sumstat,as_tibble((matrix_sumstat[n_cross,]))%>%add_column(., Type="Obs"))
    
    
    #We compute the mean squared error (RMSE) 
    RMSE = sapply(1:ncol(cross_valid$adj.values),function(x){
      sqrt(sum((cross_valid$adj.values[,x]-matrix_param[n_cross,x])**2)/nrow(cross_valid$adj.values) )
      }
    )
    
    #normalize it by the RMSE under the prior distribution
    RMSE_prior=sapply(1:ncol(matrix_param),function(x){
      sqrt(sum((matrix_param[,x]-matrix_param[n_cross,x])**2)/nrow(matrix_param) )
    }
    )
    NRMSE = RMSE/RMSE_prior
    
    d_NRMSE_param=rbind(d_NRMSE_param,as_tibble(t(NRMSE)))
    
    
    #We repeat the same for the summary statistics observed
    RMSE = sapply(1:ncol(cross_valid$ss),function(x){
      sqrt(sum((cross_valid$ss[,x]-matrix_sumstat[n_cross,x])**2)/nrow(cross_valid$ss) )
    }
    )
    
    RMSE_prior=sapply(1:ncol(matrix_sumstat),function(x){
      sqrt(sum((matrix_sumstat[,x]-matrix_sumstat[n_cross,x])**2)/nrow(matrix_sumstat) )
    }
    )
    NRMSE = RMSE/RMSE_prior
    
    d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,as_tibble(t(NRMSE)))
    
  } #end loop Nvirtual data
  
  colnames(d_NRMSE_param)=colnames(d_cross_param)[-length(colnames(d_cross_param))]
  colnames(d_NRMSE_sumstat)=colnames(d_cross_sumstat)
  dev.off()
  
  get_upper_tri=function(mat){
    mat[lower.tri(mat)]= NA
    diag(mat)=NA
    return(mat)
  }
  
  #Ploting the correlation between parameters
  colnames(mat_cor_param)=rownames(mat_cor_param)=colnames(d_NRMSE_param)[-1]
  
  p=ggplot(get_upper_tri(rowMeans(mat_cor_param, dims = 2))%>%
           melt(.)) + 
    geom_tile(aes(Var2, Var1,fill=value), color = "white")+
    geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4) +
    labs(x="",y="",fill="")+
    theme_classic()+
    theme(legend.position = "bottom")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",limit = c(-1,1),na.value = "white")
  
  ggsave(paste0("../Figures/Cross_validation/Correlation_parameters_",method_abc,".pdf"),width = 6,height = 5)
  
  
  
  #Ploting f(x,y) with x=True parameter/summary stat, y=simulated
  ##d-melting the tibble
  
  pdf(paste0("../Figures/Cross_validation/x_y_obs_true_param_",method_abc,".pdf"),width = 8,height = 5)
  par(mfrow=c(2,3))
  for (i in colnames(matrix_param)[-1]){
    d=d_cross_param[,c(which(colnames(d_cross_param)==i),ncol(d_cross_param))]%>%
      mutate(., idx=rep(1:(nrow(d_cross_param)/2),each=2))%>%
      dcast(., idx ~ Type,value.var=i)%>%
      mutate(., Parameter=i)
    
    plot(d$Obs,d$Sim,xlab="True parameter",ylab="Simulated parameter",col=alpha("blue",.5),main=i)
    abline(coef = c(0,1),col="red")
  }
  dev.off()
  
  
  
  pdf(paste0("../Figures/Cross_validation/x_y_obs_true_summarystat_",method_abc,".pdf"),width = 10,height = 7)
  
  par(mfrow=c(3,3))
  for (i in colnames(matrix_sumstat)[-1]){
    d=d_cross_sumstat[,c(which(colnames(d_cross_sumstat)==i),ncol(d_cross_sumstat))]%>%
      mutate(., idx=rep(1:(nrow(d_cross_sumstat)/2),each=2))%>%
      dcast(., idx ~ Type,value.var=i)%>%
      mutate(., Parameter=i)
    
    plot(d$Obs,d$Sim,xlab="True parameter",ylab="Simulated parameter",col=alpha("blue",.5),main=i)
    abline(coef = c(0,1),col="red")
  }
  dev.off()
  
  p=ggplot(d_NRMSE_param%>%
           melt(., id.vars="r"))+
    geom_jitter(aes(x=variable,y=value,color=variable),
                position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0), alpha=.5)+
    geom_point(data=as_tibble(t(colMeans(d_NRMSE_param)))%>%
                 melt(.,id.vars="r"),
               aes(x=variable,y=value),color="black",size=3,shape=18)+
    labs(x="Parameter",y="NRMSE",color="")+
    geom_hline(yintercept = 1)+
    ylim(0,ifelse(max(d_NRMSE_param[,-1])>3,3,max(d_NRMSE_param[,-1])))+
    theme_classic()+
    theme(legend.position = "none")
  
  ggsave(paste0("../Figures/Cross_validation/NRMSE_param_",method_abc,".pdf"),p,width = 8,height = 5)
  
  
  p=ggplot(d_NRMSE_sumstat%>%
             melt(.))+
    geom_jitter(aes(x=variable,y=value,color=variable),
                position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0), alpha=.5)+
    geom_point(data=as_tibble(t(colMeans(d_NRMSE_sumstat)))%>%
                melt(.),
              aes(x=variable,y=value),color="black",size=3,shape=18)+
    labs(x="Parameter",y="NRMSE",color="")+
    geom_hline(yintercept = 1)+
    theme_classic()+
    theme(legend.position = "none")
  
  ggsave(paste0("../Figures/Cross_validation/NRMSE_summarystat_",method_abc,".pdf"),p,width = 8,height = 5)
  
}# loop over the two method for ABC




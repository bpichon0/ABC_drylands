rm(list=ls())
source("./ABC_drylands_function.R")

#***********************************************************

# ---------------------- Step 1 : Generating parameters for sensitivity analysis ----------

#***********************************************************

dir.create("../Data/Step1_sensitivity",showWarnings = F)

## >> 1) One dimensional analysis ----


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

trends_param=trends_with_scaling=corr_param_metric=tibble()
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
    
    corr_param_metric=rbind(corr_param_metric,
                            tibble(Coeff=cor(d[,3],d[,1]),
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

p3=ggplot(corr_param_metric)+
  geom_tile(aes(x=Param,y=Metric,fill=Coeff),color="black")+
  theme_classic()+
  scale_fill_gradient2(low="#D82828",mid="white",high="#0E5DCE")

ggsave("../Figures/Sensitivity/Trends_1D_param.pdf",ggarrange(ggarrange(p1,p2,ncol=2,labels=letters[1:2]),
                                                              ggarrange(ggplot()+theme_void(),p3,ggplot()+theme_void(), ncol=3,widths = c(.5,1,.5)),
                                                              nrow=2,labels=c("",letters[3])),
       width = 10,height = 7)



## >> 2) Two dimensional analysis for parameters pairs ----

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
                      "rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                      "Spectral_ratio","PLR","PL_expo")

n_param=7;nsim=400
pdf("../Figures/Sensitivity/Sensitivity_2D.pdf",width = 10,height = 8)
index=1
for (i in 1:(n_param-1)){
  for (j in (i+1):n_param){
  
    d=results2D[((index-1)*nsim+1):(index*nsim),]
    d$Var1=d[,i]
    d$Var2=d[,j]
    
    for (metric in 1:9){
      which_metric=c("rho_p","nb_neigh","clustering","skewness","variance","moran_I","Spectral_ratio","PLR","PL_expo")[metric]
      
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
    
    print(ggarrange(p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,ncol = 3,nrow = 3 ))
  
    
    #Performing ACP
    if (any(is.nan(d$PL_expo) | is.nan(d$PLR)))  d=d[-which(is.nan(d$PL_expo) | is.nan(d$PLR)),]
    
    res.pca=PCA(d[,8:ncol(d)], scale.unit = T, ncp = 3,  graph=F)
    
    p1=fviz_pca_biplot(res.pca, geom.ind = "point", 
                 axes=c(1,2), col.ind = d$Var1,col.var="black",
                 label = "var", repel = T)+
      scale_color_gradientn(colours = 
                              colorRampPalette(colors = 
                                                 c("#960261","#C5218B","#E266B6",
                                                   "#9F94EF","#2F48DA","#2F48DA"))(100),
      )+
      labs(x=paste0("PC 1 (",round(res.pca$eig[1,2], 1)," %)"),
           y=paste0("PC 2 (",round(res.pca$eig[2,2], 1)," %)"),color=paste0(colnames(d)[i],"  "))+
      ggtitle(paste0(colnames(d)[i],"  "))+
      theme_classic()+theme(legend.position = "bottom",plot.title = element_text(size=25))
    
    
    
    p2=fviz_pca_biplot(res.pca, geom.ind = "point", 
                    axes=c(1,2), col.ind = d$Var2,col.var="black",
                    label = "var", repel = T)+
      scale_color_gradientn(colours = 
                              colorRampPalette(colors = 
                                                 c("#960261","#C5218B","#E266B6",
                                                   "#9F94EF","#2F48DA","#2F48DA"))(100),
      )+
      labs(x=paste0("PC 1 (",round(res.pca$eig[1,2], 1)," %)"),
           y=paste0("PC 2 (",round(res.pca$eig[2,2], 1)," %)"),color=paste0(colnames(d)[j],"  "))+
      ggtitle(paste0(colnames(d)[j],"  "))+
      theme_classic()+theme(legend.position = "bottom",plot.title = element_text(size=25))
    
    
    print(ggarrange(p1,p2))
    
    
    
    
    
    index=index+1
  }
}
dev.off()

#example for overleaf with f and c

index=14
d=results2D[((index-1)*nsim+1):(index*nsim),]
d$Var1=d[,3]
d$Var2=d[,6]


for (metric in 1:9){
  which_metric=c("rho_p","nb_neigh","clustering","skewness","variance","moran_I","Spectral_ratio","PLR","PL_expo")[metric]
  
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

p=ggarrange(p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,ncol = 3,nrow = 3 )
ggsave("../Figures/Sensitivity/Example_interaction_f_c.pdf",p,width = 10,height = 6)


#***********************************************************

# ---------------------- Step 2 : Cross verification ---------------------------

#***********************************************************

## >> 1) Generating pseudo-parameters using latin hypercube sampling ----

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


## >> 2) Analysis ----

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
write.table(d_all,"../Data/All_sim_ABC.csv",sep=";")



# Analyzing outputs using linear models

d_all=read.table("../Data/All_sim_ABC.csv",sep=";")

n_metric=9;n_param=7
d=tibble()
for (metric in 1:n_metric){
  
  reg0=lm(as.formula(paste(colnames(d_all)[n_param+metric],"~1")),data=d_all[,-c((n_param+1):ncol(d_all))[-which(c((n_param+1):ncol(d_all)) == metric+n_param)]])
  select_param=step(reg0,scope=as.formula(paste("~",paste(colnames(d_all)[1:n_param],collapse = "+"))),direction="both")
  d=rbind(d,as_tibble(t(select_param$coefficients[-1]))%>%
            add_column(., Y=colnames(d_all)[n_param+metric]))
}

d





#correlation summary stat


matrix_param=d_all[,1:7]
matrix_sumstat=d_all[,8:(ncol(d_all))]



mat_cor_sumstat=cor(matrix_sumstat)

p=ggplot(get_upper_tri(mat_cor_sumstat)%>%
           melt(.)) + 
  geom_tile(aes(Var2, Var1,fill=value), color = "white")+
  geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4) +
  labs(x="",y="",fill="")+
  theme_classic()+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 60,vjust=.7))+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white",limit = c(-1,1),na.value = "white")

ggsave(paste0("../Figures/Kefi_inferrence/Cross_validation/Correlation_sumstats.pdf"),width = 5,height = 5)


#Dimentionality of EWS
mat_cor_sumstat_scaled =2**((mat_cor_sumstat - min(mat_cor_sumstat)) / diff(range(mat_cor_sumstat)))
diag(mat_cor_sumstat)=NA


graph_cor=graph_from_adjacency_matrix(adjmatrix = as.matrix(mat_cor_sumstat_scaled),mode = "undirected",weighted = T)
plot(graph_cor)
modules=igraph::edge.betweenness.community(graph=graph_cor)
modules$membership #each seem to bring some information

t(matrix_sumstat)%>% scale %>% dist %>% 
  hclust %>% as.dendrogram %>%
  ggdendrogram




#cross-validation
d_all=read.table("../Data/All_sim_ABC.csv",sep=";")

condition_cover=which(d_all$rho_p < 0.1  | d_all$rho_p>.8)
matrix_param=d_all[-condition_cover,2:7]
matrix_sumstat=d_all[-condition_cover,8:(ncol(d_all))]
N_for_cross_validation = 100
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)

for (method_abc in c("rejection","neuralnet","loclinear")){
  
  mat_cor_param=array(0,c(6,6,N_for_cross_validation)) #correlation matrix for parameters
  
  pdf(paste0("../Figures/Kefi_inferrence/Cross_validation/Cross_validation_n",N_for_cross_validation,"_",method_abc,".pdf"),width = 8,height = 4)
  d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
  
  for (n in 1:N_for_cross_validation){
    
    n_cross=nrow_for_sample[n]
    
    #for each virtual data, we perform ABC rejection algorithm with linear regression adjustment for posterior
    cross_valid=abc(target = matrix_sumstat[n_cross,],
                    param = matrix_param[-n_cross,],sumstat = matrix_sumstat[-n_cross,], #removing the target data
                    tol = 100/nrow(matrix_param),method = method_abc) #we keep the 100 closest simulations
    
    if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
    
    cross_valid$adj.values=cross_valid$adj.values
    #Matrix of correlation between parameters & sumstats
    mat_cor_param[,,n]=cor(cross_valid$adj.values)

    
    #We plot the differences in posterior distribution/true parameter
    
    
    par(mfrow=c(2,4))
    for (i in colnames(matrix_param)){
      plot(density(cross_valid$adj.values[,i]),main=i,xlab="Value")
      abline(v = matrix_param[n_cross,i],col="blue")
      abline(v = colMeans(cross_valid$adj.values)[i],col="red")
    }
  
    d_melt=as.data.frame(cross_valid$ss)%>%
      melt(.)
    
    par(mfrow=c(2,4))
    for (i in 1:(length(colnames(cross_valid$adj.values))-1)){
      for (j in (i+1):length(colnames(cross_valid$adj.values))){
        plot(x=cross_valid$adj.values[,i],y=cross_valid$adj.values[,j],
             xlab=colnames(cross_valid$adj.values)[i],ylab=colnames(cross_valid$adj.values)[j],
             col=alpha("blue",.8))
      }
    }

    par(mfrow=c(2,5))
    for (i in colnames(matrix_sumstat)){
      plot(density(d_melt$value[which(d_melt$variable==i)]),main=i,xlab="Value")
      abline(v = matrix_sumstat[n_cross,i],col="blue")
      abline(v = colMeans(cross_valid$ss)[i],col="red")
    }
    
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
    RMSE_prior=sapply(1:(ncol(matrix_param)),function(x){
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
  
  
  #Ploting the correlation between parameters 
  colnames(mat_cor_param)=rownames(mat_cor_param)=colnames(d_NRMSE_param)
  
  p=ggplot(get_upper_tri(rowMeans(mat_cor_param, dims = 2,na.rm = T))%>%
           melt(.)) + 
    geom_tile(aes(Var2, Var1,fill=value), color = "white")+
    geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4) +
    labs(x="",y="",fill="")+
    theme_classic()+
    theme(legend.position = "bottom")+
    scale_fill_gradient2(low = "red", high = "blue", mid = "white",limit = c(-1,1),na.value = "white")
  
  ggsave(paste0("../Figures/Kefi_inferrence/Cross_validation/Correlation_parameters_",method_abc,".pdf"),width = 6,height = 5)
  
  
  #Ploting f(x,y) with x=True parameter/summary stat, y=simulated
  ##d-melting the tibble
  
  pdf(paste0("../Figures/Kefi_inferrence/Cross_validation/x_y_obs_true_param_",method_abc,".pdf"),width = 8,height = 5)
  par(mfrow=c(2,3))
  for (i in colnames(matrix_param)){
    d=d_cross_param[,c(which(colnames(d_cross_param)==i),ncol(d_cross_param))]%>%
      mutate(., idx=rep(1:(nrow(d_cross_param)/2),each=2))%>%
      dcast(., idx ~ Type,value.var=i)%>%
      mutate(., Parameter=i)
    
    plot(d$Obs,d$Sim,xlab="True parameter",ylab="Simulated parameter",col=alpha("blue",.5),main=i)
    abline(coef = c(0,1),col="red")
  }
  dev.off()
  
  
  
  pdf(paste0("../Figures/Kefi_inferrence/Cross_validation/x_y_obs_true_summarystat_",method_abc,".pdf"),width = 10,height = 7)
  
  par(mfrow=c(3,3))
  for (i in colnames(matrix_sumstat)){
    d=d_cross_sumstat[,c(which(colnames(d_cross_sumstat)==i),ncol(d_cross_sumstat))]%>%
      mutate(., idx=rep(1:(nrow(d_cross_sumstat)/2),each=2))%>%
      dcast(., idx ~ Type,value.var=i)%>%
      mutate(., Parameter=i)
    
    plot(d$Obs,d$Sim,xlab="True parameter",ylab="Simulated parameter",col=alpha("blue",.5),main=i)
    abline(coef = c(0,1),col="red")
  }
  dev.off()
  
  p=ggplot(d_NRMSE_param%>%
           melt(.))+
    geom_jitter(aes(x=variable,y=value,color=variable),
                position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0), alpha=.5)+
    geom_point(data=as_tibble(t(colMeans(d_NRMSE_param)))%>%
                 melt(.),
               aes(x=variable,y=value),color="black",size=3,shape=18)+
    labs(x="Parameter",y="NRMSE",color="")+
    geom_hline(yintercept = 1)+
    ylim(0,ifelse(max(d_NRMSE_param[,-1])>3,3,max(d_NRMSE_param[,-1])))+
    theme_classic()+
    theme(legend.position = "none")
  
  ggsave(paste0("../Figures/Kefi_inferrence/Cross_validation/NRMSE_param_",method_abc,".pdf"),p,width = 8,height = 5)
  
  
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
  
  ggsave(paste0("../Figures/Kefi_inferrence/Cross_validation/NRMSE_summarystat_",method_abc,".pdf"),p,width = 8,height = 5)
  
}# loop over the two method for ABC



## >> 3) PCA and variables ----


d_all=read.table("../Data/All_sim_ABC.csv",sep=";")

sample_row=sample(1:nrow(d_all),size=10000,replace = F)
res.pca=PCA(d_all[sample_row,8:ncol(d_all)], scale.unit = T, ncp = 3,  graph=F)

pdf("../Figures/Space_param_EWS/Space_EWS_parameters.pdf",width = 12,height = 4)
for (param in 1:7){
  axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))
  for (i in 1:3){
    assign(paste0("p",i),
           
           fviz_pca_biplot(res.pca, geom.ind = "point", 
                           axes=c(axes_for_plot$x[i],axes_for_plot$y[i]), col.ind = d_all[sample_row,param],col.var="black",
                           label = "var", repel = T)+
             scale_color_gradientn(colours = 
                                     colorRampPalette(colors = 
                                                        c("#960261","#C5218B","#E266B6",
                                                          "#9F94EF","#2F48DA","#2F48DA"))(100),
             )+
             labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color=colnames(d_all)[param])+
             ggtitle(paste0(colnames(d_all)[param]))+
             theme_classic()+theme(legend.position = "bottom",plot.title = element_text(size=25))
           
    )
  }
  print(ggarrange(p1,p2,p3,ncol=3))
  
}
dev.off()
ggsave("../Figures/Space_param_EWS/ACP_delta.pdf",
       ggarrange(p1+ggtitle(TeX("$\\delta$")),
                 p2+ggtitle(TeX("$\\delta$")),
                 p3+ggtitle(TeX("$\\delta$")),ncol=3),width = 12,height = 4)




#***********************************************************

# ---------------------- Step 3: Influence of the number of photos to average ----

#***********************************************************

nb_picture = c(3,5,15,25,35)
d_distance=tibble()

for (nb_pic in nb_picture){
  d=tibble()
  list_simu=list.files(paste0("../Data/Step3_nbr_pictures/",nb_pic,"_pic"))
  
  for (simu_file in list_simu){
    d=rbind(d,read.table(paste0("../Data/Step3_nbr_pictures/",nb_pic,"_pic/",simu_file),sep=","))
  }
  colnames(d)= c("r","d", "f", "m","b","c","delta",
                 "rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                 "Spectral_ratio","PLR","PL_expo")
  
  write.table(d,paste0("../Data/Step3_nbr_pictures/Simu_",nb_pic,"_pictures.csv"),sep=";")

  condition_cover=which(d$rho_p < 0.1  | d$rho_p>.8)
  matrix_param=d[-condition_cover,2:7]
  matrix_sumstat=d[-condition_cover,8:(ncol(d))]
  N_for_cross_validation = 150
  set.seed(123)
  nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)
  
  for (method_abc in c("rejection")){
    
    mat_cor_param=array(0,c(6,6,N_for_cross_validation)) #correlation matrix for parameters
    
    d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
    
    for (n in 1:N_for_cross_validation){
      
      n_cross=nrow_for_sample[n]
      
      #for each virtual data, we perform ABC rejection algorithm with linear regression adjustment for posterior
      cross_valid=abc(target = matrix_sumstat[n_cross,],
                      param = matrix_param[-n_cross,],sumstat = matrix_sumstat[-n_cross,], #removing the target data
                      tol = 100/nrow(matrix_param),method = method_abc) #we keep the 100 closest simulations
      
      if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
      
      cross_valid$adj.values=cross_valid$adj.values
      #Matrix of correlation between parameters & sumstats
      mat_cor_param[,,n]=cor(cross_valid$adj.values)
      
      
      
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
      RMSE_prior=sapply(1:(ncol(matrix_param)),function(x){
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
  
    
    #Ploting f(x,y) with x=True parameter/summary stat, y=simulated
    ##d-melting the tibble
    
    pdf(paste0("../Figures/Kefi_inferrence/Number_pictures/x_y_obs_true_param_",method_abc,"_nbpic_",nb_pic,".pdf"),width = 8,height = 5)
    par(mfrow=c(2,3))
    for (i in colnames(matrix_param)){
      d=d_cross_param[,c(which(colnames(d_cross_param)==i),ncol(d_cross_param))]%>%
        mutate(., idx=rep(1:(nrow(d_cross_param)/2),each=2))%>%
        dcast(., idx ~ Type,value.var=i)%>%
        mutate(., Parameter=i)
      
      plot(d$Obs,d$Sim,xlab="True parameter",ylab="Simulated parameter",col=alpha("blue",.5),main=i)
      abline(coef = c(0,1),col="red")
    }
    dev.off()
    
    
    
    pdf(paste0("../Figures/Kefi_inferrence/Number_pictures/x_y_obs_true_summarystat_",method_abc,"_nbpic_",nb_pic,".pdf"),width = 10,height = 7)
    
    par(mfrow=c(3,3))
    for (i in colnames(matrix_sumstat)){
      d=d_cross_sumstat[,c(which(colnames(d_cross_sumstat)==i),ncol(d_cross_sumstat))]%>%
        mutate(., idx=rep(1:(nrow(d_cross_sumstat)/2),each=2))%>%
        dcast(., idx ~ Type,value.var=i)%>%
        mutate(., Parameter=i)
      
      plot(d$Obs,d$Sim,xlab="True parameter",ylab="Simulated parameter",col=alpha("blue",.5),main=i)
      abline(coef = c(0,1),col="red")
    }
    dev.off()
    
    p=ggplot(d_NRMSE_param%>%
               melt(.))+
      geom_jitter(aes(x=variable,y=value,color=variable),
                  position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0), alpha=.5)+
      geom_point(data=as_tibble(t(colMeans(d_NRMSE_param,na.rm = T)))%>%
                   melt(.),
                 aes(x=variable,y=value),color="black",size=3,shape=18)+
      labs(x="Parameter",y="NRMSE",color="")+
      geom_hline(yintercept = 1)+
      ylim(0,ifelse(max(d_NRMSE_param[,-1])>3,3,max(d_NRMSE_param[,-1])))+
      theme_classic()+
      theme(legend.position = "none")
    
    ggsave(paste0("../Figures/Kefi_inferrence/Number_pictures/NRMSE_param_",method_abc,"_nbpic_",nb_pic,".pdf"),p,width = 8,height = 5)
    
    
    #keeping the distance of NRMSE to 1 for each 
    
    d_distance=rbind(d_distance,as_tibble(t(colMeans(d_NRMSE_param,na.rm = T)))%>%add_column(., nb_picture=nb_pic))
    
    p=ggplot(d_NRMSE_sumstat%>%
               melt(.))+
      geom_jitter(aes(x=variable,y=value,color=variable),
                  position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0), alpha=.5)+
      geom_point(data=as_tibble(t(colMeans(d_NRMSE_sumstat,na.rm = T)))%>%
                   melt(.),
                 aes(x=variable,y=value),color="black",size=3,shape=18)+
      labs(x="Parameter",y="NRMSE",color="")+
      geom_hline(yintercept = 1)+
      theme_classic()+
      theme(legend.position = "none")
    
    ggsave(paste0("../Figures/Kefi_inferrence/Number_pictures/NRMSE_summarystat_",method_abc,"_nbpic_",nb_pic,".pdf"),p,width = 8,height = 5)
    
  }# loop over the two method for ABC
  
}


p=ggplot(d_distance%>%melt(., id.vars=c("nb_picture")))+
  geom_line(aes(x=nb_picture,y=value,color=variable),lwd=1)+
  geom_point(aes(x=nb_picture,y=value,color=variable),size=3,shape=21,fill="white")+
  geom_point(data=tibble(mean_NRMSE=rowMeans(d_distance[,-ncol(d_distance)]),nb_picture=d_distance$nb_picture),
             aes(x=nb_picture,y=mean_NRMSE),shape=15)+
  geom_line(data=tibble(mean_NRMSE=rowMeans(d_distance[,-ncol(d_distance)]),nb_picture=d_distance$nb_picture),
             aes(x=nb_picture,y=mean_NRMSE),linetype=9)+
  the_theme+
  geom_hline(yintercept = 1)+
  labs(x="# of pictures",y="Mean NRMSE (150 virtual data \n and keeping the best 100 datasets)",color="")+
  scale_x_continuous(breaks = d_distance$nb_picture)
  
ggsave("../Figures/Kefi_inferrence/Number_pictures/Influence_#_pictures.pdf",p,width = 6,height = 4)








#***********************************************************

# ---------------------- Step 4: Fixing some parameters, varying others ----

#***********************************************************


## >> 1) Pseudo-parameters ----
# 14 combinations: we fix b or d, or r or c and their combination.

set.seed(123)


experience_parameters=matrix(0,nrow = 14,7) #14 different combinations. As maximum: 3 parameters fixed
colnames(experience_parameters)=c('r',"d", "f", "m","b","c","delta")
experience_parameters[,1]=c(0,1,0,1,1,0,1,rep(0,7))
experience_parameters[,2]=c(rep(0,7),0,1,0,1,1,0,1)
experience_parameters[,3]=c(rep(0,14))
experience_parameters[,4]=c(rep(0,14))
experience_parameters[,5]=rep(c(1,0,0,0,1,1,1),2)
experience_parameters[,6]=rep(c(0,0,1,1,0,1,1),2)
experience_parameters[,7]=c(rep(0,14))
based_param=c('r'=.05,"d"=.1, "f"=.9, "m"=.1,"b"=.8,"c"=.2,"delta"=.1)

d_prior_all=tibble()
for (i in 1:nrow(experience_parameters)){
  
  range_priors=data.frame(min = c(0,0, 0, 0.005,0,0,0),
                          max = c(1,1, 1, 1,1,1,1))
  rownames(range_priors)=c('r',"d", "f", "m","b","c","delta")
  
  # Latin hypercube sampling on the priors
  pseudo_param=as.data.frame(Latinhyper(range_priors, 4e4))
  pseudo_param[, "m"]=qlnorm(pseudo_param[,"m"], -2.25,.7)
  
  
  for (j in 1:length(names(which(experience_parameters[i,]!=0)))){
    pseudo_param[,names(which(experience_parameters[i,]!=0))[j]] = based_param[names(which(experience_parameters[i,]!=0))[j]] 
  }

  #reorganizing
  d_prior_all=rbind(d_prior_all,pseudo_param)
}

write.table(d_prior_all,paste0("../Data/Pseudo_param_all_combinations.csv"),sep=";",row.names = F)


## >> 2) Analysis ----


for (virtual_exp in 1:14){ #for each combination of parameters that were knocked-out, we do the analysis
  
  d=tibble()
  for (simu_id in ((virtual_exp-1)*400+1):(virtual_exp*400)){
    d=rbind(d,read.table(paste0("../Data/Step4_combination_param/Simulation_ABC_number_",simu_id,".csv"),sep=","))
  }
  
  colnames(d)= c("r","d", "f", "m","b","c","delta",
                 "rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                 "Spectral_ratio","PLR","PL_expo")
  
  knocked_param=which(experience_parameters[virtual_exp,] == 1)
  
  condition_cover=which(d$rho_p < 0.1  | d$rho_p>.8)
  matrix_param=d[-condition_cover,1:7] #removing very low cover and very high cover where PL can't be fitter 
  
  print(nrow(matrix_param))
  
  matrix_param=matrix_param[,-knocked_param] #removing the fixed parameters
  
  matrix_sumstat=d[-condition_cover,8:(ncol(d))]
  N_for_cross_validation = 100
  set.seed(123)
  nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F) # the virtual data sampled
  
  for (method_abc in c("fixed_nb")){ #old 20/01, does not change anything for (method_abc in c("fixed_tol","fixed_nb")){
    
    mat_cor_param=array(0,c(7-length(knocked_param),
                            7-length(knocked_param),
                            N_for_cross_validation)) #correlation matrix for parameters
    
    d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
    
    for (n in 1:N_for_cross_validation){
      
      n_cross=nrow_for_sample[n]
      
      #for each virtual data, we perform ABC rejection algorithm 
      cross_valid=abc(target = matrix_sumstat[n_cross,],
                      param = matrix_param[-n_cross,],sumstat = matrix_sumstat[-n_cross,], #removing the target data
                      tol = ifelse(method_abc=="fixed_tol",.002,100/nrow(matrix_param)),method = "rejection") #we keep the 100 closest simulations
      
      if (names(cross_valid)[1]=="unadj.values") names(cross_valid)[1] = "adj.values"
      
      cross_valid$adj.values=cross_valid$adj.values
      
      #Matrix of correlation between parameters & sumstats
      mat_cor_param[,,n]=cor(cross_valid$adj.values)
      
      #we save the mean posterior distribution for each and the true observed parameters
      d_cross_param=rbind(d_cross_param,as_tibble(t(colMeans(cross_valid$adj.values)))%>%add_column(., Type="Sim"))
      d_cross_param=rbind(d_cross_param,as_tibble((matrix_param[n_cross,]))%>%add_column(., Type="Obs"))
      
      #As we work with virtual data, we do the same for the summary stats we save the mean posterior distribution for each and the true observed parameters
      d_cross_sumstat=rbind(d_cross_sumstat,as_tibble(t(colMeans(cross_valid$ss)))%>%add_column(., Type="Sim"))
      d_cross_sumstat=rbind(d_cross_sumstat,as_tibble((matrix_sumstat[n_cross,]))%>%add_column(., Type="Obs"))
      
      
      #We compute the mean squared error (RMSE) 
      RMSE = sapply(1:ncol(cross_valid$adj.values),function(x){
        sqrt(sum((cross_valid$adj.values[,x]-matrix_param[n_cross,x])**2)/nrow(cross_valid$adj.values) )}
      )
      
      #normalize it by the RMSE under the prior distribution
      RMSE_prior=sapply(1:(ncol(matrix_param)),function(x){
        sqrt(sum((matrix_param[,x]-matrix_param[n_cross,x])**2)/nrow(matrix_param) )})
      
      NRMSE = RMSE/RMSE_prior
      d_NRMSE_param=rbind(d_NRMSE_param,as_tibble(t(NRMSE)))
      
      
      #We repeat the same for the summary statistics observed
      RMSE = sapply(1:ncol(cross_valid$ss),function(x){
        sqrt(sum((cross_valid$ss[,x]-matrix_sumstat[n_cross,x])**2)/nrow(cross_valid$ss) )
      }
      )
      
      RMSE_prior=sapply(1:ncol(matrix_sumstat),function(x){
        sqrt(sum((matrix_sumstat[,x]-matrix_sumstat[n_cross,x])**2)/nrow(matrix_sumstat) )})
      
      NRMSE = RMSE/RMSE_prior
      d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,as_tibble(t(NRMSE)))
      
    } #end loop Nvirtual data
    
    colnames(d_NRMSE_param)=colnames(d_cross_param)[-length(colnames(d_cross_param))]
    colnames(d_NRMSE_sumstat)=colnames(d_cross_sumstat)
    
    
    
    p=ggplot(d_NRMSE_param%>%
               melt(.))+
      geom_jitter(aes(x=variable,y=value,color=variable),
                  position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0), alpha=.5)+
      geom_point(data=as_tibble(t(colMeans(d_NRMSE_param)))%>%
                   melt(.),
                 aes(x=variable,y=value),color="black",size=3,shape=18)+
      labs(x="Parameter",y="NRMSE",color="")+
      geom_hline(yintercept = 1)+
      ylim(0,ifelse(max(d_NRMSE_param[,-1])>3,3,max(d_NRMSE_param[,-1])))+
      theme_classic()+
      theme(legend.position = "none")+
      ggtitle(paste0(c("r","d", "f", "m","b","c","delta")[knocked_param],collapse = ", "))
    
    ggsave(paste0("../Figures/Kefi_inferrence/Combination_param/NRMSE_param_",virtual_exp,"_",method_abc,".pdf"),p,width = 8,height = 5)
    
    
    
    pdf(paste0("../Figures/Kefi_inferrence/Combination_param/x_y_obs_true_param_",virtual_exp,"_",method_abc,".pdf"),width = 8,height = 5)
    par(mfrow=c(2,4))
    for (i in colnames(matrix_param)){
      d=d_cross_param[,c(which(colnames(d_cross_param)==i),ncol(d_cross_param))]%>%
        mutate(., idx=rep(1:(nrow(d_cross_param)/2),each=2))%>%
        dcast(., idx ~ Type,value.var=i)%>%
        mutate(., Parameter=i)
      
      plot(d$Obs,d$Sim,xlab="True parameter",ylab="Simulated parameter",col=alpha("blue",.5),main=i)
      abline(coef = c(0,1),col="red")
    }
    dev.off()
    
    

    
    #Ploting the correlation between parameters 
    colnames(mat_cor_param)=rownames(mat_cor_param)=colnames(d_NRMSE_param)
    
    p=ggplot(get_upper_tri(rowMeans(mat_cor_param, dims = 2,na.rm = T))%>%
               melt(.)) + 
      geom_tile(aes(Var2, Var1,fill=value), color = "white")+
      geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4) +
      labs(x="",y="",fill="")+
      theme_classic()+
      theme(legend.position = "bottom")+
      scale_fill_gradient2(low = "red", high = "blue", mid = "white",limit = c(-1,1),na.value = "white")
    
    ggsave(paste0("../Figures/Kefi_inferrence/Combination_param/Correlation_parameters_",virtual_exp,"_",method_abc,".pdf"),width = 6,height = 5)
    
    
  }
  
}












#***********************************************************

# ---------------------- Step 5: Testing two step procedure proposed by Siren et al., 2019 ----

#***********************************************************

#cross-validation
d_all=read.table("../Data/All_sim_ABC.csv",sep=";")
condition_cover=which(d_all$rho_p < 0.1  | d_all$rho_p>.8)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
matrix_param=d_all[,2:7]
matrix_sumstat=d_all[,8:(ncol(d_all))]
N_for_cross_validation = 50
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)

for (two_step in c(T,F)){
  
  mat_cor_param=array(0,c(6,6,N_for_cross_validation)) #correlation matrix for parameters
  
  d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
  
  for (n in 1:N_for_cross_validation){
    
    matrix_param=d_all[,2:7]
    matrix_sumstat=d_all[,8:(ncol(d_all))]
    save_sumstat=matrix_sumstat
    
    n_cross=nrow_for_sample[n]
    
    if (two_step){ #Applying the two step procedure used in Siren MEE paper : Don't know whether it make sense in our case. TO discuss Monday

      #First box cox transformation of variables to that they approach normality
      #As we have negative values, we used the transformation coined by Manly in 1971
      for (x in 1:ncol(matrix_sumstat)) if (x %in% c(4,6)){
        matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(.5)) -1)/(.5)
      }else {matrix_sumstat[,x] = (matrix_sumstat[,x]^(.5) -1)/(.5)}
      
      #Second we scale
      for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)

      #and finally, we perform the first PLS
  

      pls_1=plsr(f + d + m + b + c + delta~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
           data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
      
      
      n_comp_pls=selectNcomp(pls_1,method = "onesigma")
      
      
      if (n_comp_pls > 1){
        mat_sumstat_pls=pls_1$scores[,1:n_comp_pls] # selecting # components
      } else if (n_comp_pls==1){ #otherwise we take the whole components
        mat_sumstat_pls=matrix(pls_1$scores[,1:n_comp_pls],ncol=1)
      } else {mat_sumstat_pls=pls_1$scores[,1:ncol(pls_1$scores)]}
      
      
      cross_valid=abc(target = mat_sumstat_pls[n_cross,],
                      param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                      tol = 2000/nrow(matrix_param),method = "rejection") #we keep the 2000 closest simulations for the first step
      
      #Keeping 2000 simulations and doing the same steps again: normality, scaling and PLS
      
      mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid$ss)),8:(ncol(d_all))] #we keep information with the true values
      mat_sumstat_step1=rbind(mat_sumstat_step1,d_all[n_cross,8:(ncol(d_all))])
      
      #again, first box cox
      for (x in 1:ncol(mat_sumstat_step1)) if (x %in% c(4,6)){
        mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(.5)) -1)/(.5)
      }else {mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(.5) -1)/(.5)}
      
      #and normalization
      for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
      
      pls_2=plsr(f + d + m + b + c + delta~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                          mat_sumstat_step1)), scale=TRUE, validation="CV")
      
      
      n_comp_pls=selectNcomp(pls_2,method = "onesigma")
      

      if (n_comp_pls > 1){
        mat_sumstat_pls2=pls_2$scores[,1:n_comp_pls] #pls 2 selecting # components
      } else if (n_comp_pls==1){ #otherwise we take the whole components
        mat_sumstat_pls2=matrix(pls_2$scores[,1:n_comp_pls],ncol=1)
      } else {mat_sumstat_pls2=pls_2$scores[,1:ncol(pls_2$scores)]}
      
      cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),], #removing the target data
                      tol = 100/nrow(mat_sumstat_pls2),method = "rejection") #we keep the 100 closest simulations
      
      cross_valid$ss=d_all[as.numeric(rownames(cross_valid$ss)),8:(ncol(d_all))] #we keep information with the true values
      
      
      
      
    } else {
      
      #for each virtual data, we perform ABC rejection algorithm with linear regression adjustment for posterior
      cross_valid=abc(target = matrix_sumstat[n_cross,],
                      param = matrix_param[-n_cross,],sumstat = matrix_sumstat[-n_cross,], #removing the target data
                      tol = 100/nrow(matrix_param),method = "rejection") #we keep the 100 closest simulations
    }
    
    
    matrix_sumstat=save_sumstat
    
    
    if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
    
    cross_valid$adj.values=cross_valid$adj.values
    mat_cor_param[,,n]=cor(cross_valid$adj.values)
    
    
    #We plot the differences in posterior distribution/true parameter

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
    RMSE_prior=sapply(1:(ncol(matrix_param)),function(x){
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

  
  #Ploting the correlation between parameters 
  colnames(mat_cor_param)=rownames(mat_cor_param)=colnames(d_NRMSE_param)
  
  p=ggplot(get_upper_tri(rowMeans(mat_cor_param, dims = 2,na.rm = T))%>%
             melt(.)) + 
    geom_tile(aes(Var2, Var1,fill=value), color = "white")+
    geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4) +
    labs(x="",y="",fill="")+
    theme_classic()+
    theme(legend.position = "bottom")+
    scale_fill_gradient2(low = "red", high = "blue", mid = "white",limit = c(-1,1),na.value = "white")
  
  ggsave(paste0("../Figures/Kefi_inferrence/Siren_two_steps/Correlation_parameters_",ifelse(two_step,"twostep","classic"),".pdf"),width = 6,height = 5)
  
  
  #Ploting f(x,y) with x=True parameter/summary stat, y=simulated
  ##d-melting the tibble
  

  p=ggplot(d_NRMSE_param%>%
             melt(.))+
    geom_jitter(aes(x=variable,y=value,color=variable),
                position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0), alpha=.5)+
    geom_point(data=as_tibble(t(colMeans(d_NRMSE_param)))%>%
                 melt(.),
               aes(x=variable,y=value),color="black",size=3,shape=18)+
    labs(x="Parameter",y="NRMSE",color="")+
    geom_hline(yintercept = 1)+
    ylim(0,ifelse(max(d_NRMSE_param[,-1])>3,3,max(d_NRMSE_param[,-1])))+
    theme_classic()+
    theme(legend.position = "none")
  
  ggsave(paste0("../Figures/Kefi_inferrence/Siren_two_steps/NRMSE_param_",ifelse(two_step,"twostep","classic"),".pdf"),p,width = 8,height = 5)
  
  
}# loop over the two method for ABC





#***********************************************************

# ---------------------- Step 6: Adding post-processing with linear regression ----

#***********************************************************




#cross-validation
d_all=read.table("../Data/All_sim_ABC.csv",sep=";")
condition_cover=which(d_all$rho_p < 0.1  | d_all$rho_p>.8)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
matrix_param=d_all[,2:7]
matrix_sumstat=d_all[,8:(ncol(d_all))]
N_for_cross_validation = 50
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)

for (method_abc in c("rejection","loclinear","neuralnet")){
  
  mat_cor_param=array(0,c(6,6,N_for_cross_validation)) #correlation matrix for parameters
  d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
  
  for (n in 1:N_for_cross_validation){
    
    matrix_param=d_all[,2:7]
    matrix_sumstat=d_all[,8:(ncol(d_all))]
    save_sumstat=matrix_sumstat

    n_cross=nrow_for_sample[n]
    
    
    #First box cox transformation of variables to that they approach normality
    #As we have negative values, we used the transformation coined by Manly in 1971
    for (x in 1:ncol(matrix_sumstat)) if (x %in% c(4,6)){
      matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(.5)) -1)/(.5)
    }else {matrix_sumstat[,x] = (matrix_sumstat[,x]^(.5) -1)/(.5)}
    
    #Second we scale
    for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
    
    #and finally, we perform the first PLS
    
    
    pls_1=plsr(f + d + m + b + c + delta~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
               data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
    
    
    n_comp_pls=selectNcomp(pls_1,method = "onesigma")
    
    
    if (n_comp_pls > 1){
      mat_sumstat_pls=pls_1$scores[,1:n_comp_pls] # selecting # components
    } else if (n_comp_pls==1){ #otherwise we take the whole components
      mat_sumstat_pls=matrix(pls_1$scores[,1:n_comp_pls],ncol=1)
    } else {mat_sumstat_pls=pls_1$scores[,1:ncol(pls_1$scores)]}
    
    
    cross_valid=abc(target = mat_sumstat_pls[n_cross,],
                    param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                    tol = 2000/nrow(matrix_param),method = "rejection") #we keep the 2000 closest simulations for the first step
    
    #Keeping 2000 simulations and doing the same steps again: normality, scaling and PLS
    
    mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid$ss)),8:(ncol(d_all))] #we keep information with the true values
    mat_sumstat_step1=rbind(mat_sumstat_step1,d_all[n_cross,8:(ncol(d_all))])
    
    #again, first box cox
    for (x in 1:ncol(mat_sumstat_step1)) if (x %in% c(4,6)){
      mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(.5)) -1)/(.5)
    }else {mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(.5) -1)/(.5)}
    
    #and normalization
    for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
    
    pls_2=plsr(f + d + m + b + c + delta~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
               data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                        mat_sumstat_step1)), scale=TRUE, validation="CV")
    
    
    n_comp_pls=selectNcomp(pls_2,method = "onesigma")
    
    
    if (n_comp_pls > 1){
      mat_sumstat_pls2=pls_2$scores[,1:n_comp_pls] #pls 2 selecting # components
    } else if (n_comp_pls==1){ #otherwise we take the whole components
      mat_sumstat_pls2=matrix(pls_2$scores[,1:n_comp_pls],ncol=1)
    } else {mat_sumstat_pls2=pls_2$scores[,1:ncol(pls_2$scores)]}
    
    cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),],
                    param = cross_valid$unadj.values,
                    sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),], #removing the target data
                    tol = 100/nrow(mat_sumstat_pls2),method = method_abc,transf = rep("logit",6), #as parameters are proba, we perform logit regression
                    logit.bounds = matrix(c(0,1),6,2,byrow = T))
    
    cross_valid$ss=d_all[as.numeric(rownames(cross_valid$ss)),8:(ncol(d_all))] #we keep information with the true values
    if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
    
  
    mat_cor_param[,,n]=cor(cross_valid$adj.values)
    
    
    #We plot the differences in posterior distribution/true parameter
    
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
    RMSE_prior=sapply(1:(ncol(matrix_param)),function(x){
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
  matrix_sumstat=save_sumstat
  
  
  #Ploting the correlation between parameters 
  colnames(mat_cor_param)=rownames(mat_cor_param)=colnames(d_NRMSE_param)
  
  p=ggplot(get_upper_tri(rowMeans(mat_cor_param, dims = 2,na.rm = T))%>%
             melt(.)) + 
    geom_tile(aes(Var2, Var1,fill=value), color = "white")+
    geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4) +
    labs(x="",y="",fill="")+
    theme_classic()+
    theme(legend.position = "bottom")+
    scale_fill_gradient2(low = "red", high = "blue", mid = "white",limit = c(-1,1),na.value = "white")
  
  ggsave(paste0("../Figures/Kefi_inferrence/Adding_postprocessing/Correlation_parameters_",method_abc,"post_processing.pdf"),width = 6,height = 5)
  
  
  #Ploting f(x,y) with x=True parameter/summary stat, y=simulated
  ##d-melting the tibble
  
  
  p=ggplot(d_NRMSE_param%>%
             melt(.))+
    geom_jitter(aes(x=variable,y=value,color=variable),
                position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0), alpha=.5)+
    geom_point(data=as_tibble(t(colMeans(d_NRMSE_param)))%>%
                 melt(.),
               aes(x=variable,y=value),color="black",size=3,shape=18)+
    labs(x="Parameter",y="NRMSE",color="")+
    geom_hline(yintercept = 1)+
    ylim(0,ifelse(max(d_NRMSE_param[,-1])>3,3,max(d_NRMSE_param[,-1])))+
    theme_classic()+
    theme(legend.position = "none")
  
  ggsave(paste0("../Figures/Kefi_inferrence/Adding_postprocessing/NRMSE_param_",method_abc,".pdf"),p,width = 8,height = 5)
  
}







#***********************************************************

# ---------------------- Step 7: Sensitivity analysis on the landscape size (50, 75, 100, 125) ----

#***********************************************************
#XXX do figures & sims



#***********************************************************

# ---------------------- Step 8: Testing with the Eby model  ----

#***********************************************************

## >> 1) Pseudo-parameters ----


set.seed(123)
range_priors=data.frame(min = c(0, 0),
                        max = c(1, 1))
rownames(range_priors)=c("p","q")

# Latin hypercube sampling on the priors
pseudo_param=as.data.frame(Latinhyper(range_priors, 100000))
write.table(pseudo_param,'../Data/Pseudo_parameters_Eby.csv',sep=";",row.names = F)









## >> 2) Analysis ----

list_simu=list.files('../Data/Step5_Eby_model',pattern = ".csv")

d_all=tibble()
for (file_simu in list_simu){
  
  d=read.table(paste0("../Data/Step5_Eby_model/",file_simu),sep=",")
  colnames(d)= c("p","q",
                 "rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                 "Spectral_ratio","PLR","PL_expo")
  
  d_all=rbind(d_all,d)
}
d_all=d_all[-which(is.nan(d_all$PLR) | is.nan(d_all$PL_expo) | sign(d_all$PL_exp)== -1),]
write.table(d_all,"../Data/All_sim_Eby.csv",sep=";")



# Analyzing outputs using linear models


d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
matrix_param=d_all[,1:2]
matrix_sumstat=d_all[,3:(ncol(d_all))]

N_for_cross_validation = 100
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)

for (method_abc in c("rejection","loclinear","neuralnet")){
  
  for (two_step in c(T,F)){
    
    mat_cor_param=array(0,c(2,2,N_for_cross_validation)) #correlation matrix for parameters
    
    pdf(paste0("../Figures/Eby_model/First_ABC/Cross_validation_n",N_for_cross_validation,"_",
               ifelse(two_step,"twostep","classic"),"_",method_abc,".pdf"),width = 8,height = 4)
    
    d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
    
    for (n in 1:N_for_cross_validation){
      
      matrix_param=d_all[,1:2]
      matrix_sumstat=d_all[,3:(ncol(d_all))]
      save_sumstat=matrix_sumstat
      
      n_cross=nrow_for_sample[n]
      
      if (two_step){ #Applying the two step procedure used in Siren MEE paper : Don't know whether it make sense in our case. TO discuss Monday
        
        #First box cox transformation of variables to that they approach normality
        #As we have negative values, we used the transformation coined by Manly in 1971
        for (x in 1:ncol(matrix_sumstat)) if (x %in% c(4,6)){
          matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(.5)) -1)/(.5)
        }else {matrix_sumstat[,x] = (matrix_sumstat[,x]^(.5) -1)/(.5)}
        
        #Second we scale
        for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
        
        #and finally, we perform the first PLS
        
        pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                   data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
        
        
        n_comp_pls=selectNcomp(pls_1,method = "onesigma")
        
        
        if (n_comp_pls > 1){
          mat_sumstat_pls=pls_1$scores[,1:n_comp_pls] # selecting # components
        } else if (n_comp_pls==1){ #otherwise we take the whole components
          mat_sumstat_pls=matrix(pls_1$scores[,1:n_comp_pls],ncol=1)
        } else {mat_sumstat_pls=pls_1$scores[,1:ncol(pls_1$scores)]}
        
        
        cross_valid=abc(target = mat_sumstat_pls[n_cross,],
                        param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                        tol = 2000/nrow(matrix_param),method = "rejection") #we keep the 2000 closest simulations for the first step
        
        #Keeping 2000 simulations and doing the same steps again: normality, scaling and PLS
        
        mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid$ss)),3:(ncol(d_all))] #we keep information with the true values
        mat_sumstat_step1=rbind(mat_sumstat_step1,d_all[n_cross,3:(ncol(d_all))])
        
        #again, first box cox
        for (x in 1:ncol(mat_sumstat_step1)) if (x %in% c(4,6)){
          mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(.5)) -1)/(.5)
        }else {mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(.5) -1)/(.5)}
        
        #and normalization
        for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
        
        pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                   data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                            mat_sumstat_step1)), scale=TRUE, validation="CV")
        
        
        n_comp_pls=selectNcomp(pls_2,method = "onesigma")
        
        
        if (n_comp_pls > 1){
          mat_sumstat_pls2=pls_2$scores[,1:n_comp_pls] #pls 2 selecting # components
        } else if (n_comp_pls==1){ #otherwise we take the whole components
          mat_sumstat_pls2=matrix(pls_2$scores[,1:n_comp_pls],ncol=1)
        } else {mat_sumstat_pls2=pls_2$scores[,1:ncol(pls_2$scores)]}
        
        cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),],
                        param = cross_valid$unadj.values,
                        sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),], #removing the target data
                        tol = 100/nrow(mat_sumstat_pls2),method = method_abc,transf = rep("logit",2), #as parameters are proba, we perform logit regression
                        logit.bounds = matrix(c(0,1),2,2,byrow = T)) 
        
        cross_valid$ss=d_all[as.numeric(rownames(cross_valid$ss)),3:(ncol(d_all))] #we keep information with the true values
        
        
        
        
      } else {
        
        #for each virtual data, we perform ABC rejection algorithm with linear regression adjustment for posterior
        cross_valid=abc(target = matrix_sumstat[n_cross,],
                        param = matrix_param[-n_cross,],sumstat = matrix_sumstat[-n_cross,], #removing the target data
                        tol = 100/nrow(matrix_param),method = method_abc,transf = rep("logit",2), #as parameters are proba, we perform logit regression
                        logit.bounds = matrix(c(0,1),2,2,byrow = T))
      }    
      
      
      matrix_sumstat=save_sumstat
      
      
      if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
      
      cross_valid$adj.values=cross_valid$adj.values
      #Matrix of correlation between parameters & sumstats
      mat_cor_param[,,n]=cor(cross_valid$adj.values)
      
      
      #We plot the differences in posterior distribution/true parameter
      
      
      par(mfrow=c(1,2))
      for (i in colnames(matrix_param)){
        plot(density(cross_valid$adj.values[,i]),main=i,xlab="Value")
        abline(v = matrix_param[n_cross,i],col="blue")
        abline(v = colMeans(cross_valid$adj.values)[i],col="red")
      }
      
      
      par(mfrow=c(1,1))
      for (i in 1:(length(colnames(cross_valid$adj.values))-1)){
        for (j in (i+1):length(colnames(cross_valid$adj.values))){
          plot(x=cross_valid$adj.values[,i],y=cross_valid$adj.values[,j],
               xlab=colnames(cross_valid$adj.values)[i],ylab=colnames(cross_valid$adj.values)[j],
               col=alpha("blue",.8))
        }
      }
      
      
      d_melt=as.data.frame(cross_valid$ss)%>%
        melt(.)
      
      par(mfrow=c(2,5))
      for (i in colnames(matrix_sumstat)){
        plot(density(d_melt$value[which(d_melt$variable==i)]),main=i,xlab="Value")
        abline(v = matrix_sumstat[n_cross,i],col="blue")
        abline(v = colMeans(cross_valid$ss)[i],col="red")
      }
      
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
      RMSE_prior=sapply(1:(ncol(matrix_param)),function(x){
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
    
    write.table(d_NRMSE_param,paste0("../Data/Step5_Eby_model/RMSE/RMSE_",
                                     ifelse(two_step,"twostep","classic"),"_",method_abc,".csv"),sep=";")
    
    
    
    colnames(d_NRMSE_sumstat)=colnames(d_cross_sumstat)
    dev.off()
    
    
    
    
    #Ploting the correlation between parameters 
    colnames(mat_cor_param)=rownames(mat_cor_param)=colnames(d_NRMSE_param)
    
    p=ggplot(get_upper_tri(rowMeans(mat_cor_param, dims = 2,na.rm = T))%>%
               melt(.)) + 
      geom_tile(aes(Var2, Var1,fill=value), color = "white")+
      geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4) +
      labs(x="",y="",fill="")+
      theme_classic()+
      theme(legend.position = "bottom")+
      scale_fill_gradient2(low = "red", high = "blue", mid = "white",limit = c(-1,1),na.value = "white")
    
    ggsave(paste0("../Figures/Eby_model/First_ABC/Correlation_parameters_",ifelse(two_step,"twostep","classic"),"_",method_abc,".pdf"),width = 6,height = 5)
    
    
    #Ploting f(x,y) with x=True parameter/summary stat, y=simulated
    ##d-melting the tibble
    
    pdf(paste0("../Figures/Eby_model/First_ABC/x_y_obs_true_param_",ifelse(two_step,"twostep","classic"),"_",method_abc,".pdf"),width = 6,height = 4)
    par(mfrow=c(1,2))
    for (i in colnames(matrix_param)){
      d=d_cross_param[,c(which(colnames(d_cross_param)==i),ncol(d_cross_param))]%>%
        mutate(., idx=rep(1:(nrow(d_cross_param)/2),each=2))%>%
        dcast(., idx ~ Type,value.var=i)%>%
        mutate(., Parameter=i)
      
      plot(d$Obs,d$Sim,xlab="True parameter",ylab="Simulated parameter",col=alpha("blue",.5),main=i)
      abline(coef = c(0,1),col="red")
    }
    dev.off()
    
    
    
    p=ggplot(d_NRMSE_param%>%
               melt(.))+
      geom_jitter(aes(x=variable,y=value,color=variable),
                  position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0), alpha=.5)+
      geom_point(data=as_tibble(t(colMeans(d_NRMSE_param)))%>%
                   melt(.),
                 aes(x=variable,y=value),color="black",size=3,shape=18)+
      labs(x="Parameter",y="NRMSE",color="")+
      geom_hline(yintercept = 1)+
      ylim(0,ifelse(max(d_NRMSE_param[,-1])>3,3,max(d_NRMSE_param[,-1])))+
      theme_classic()+
      theme(legend.position = "none")
    
    ggsave(paste0("../Figures/Eby_model/First_ABC/NRMSE_param_",ifelse(two_step,"twostep","classic"),"_",method_abc,".pdf"),p,width = 6,height = 3)
    
  }
}


# Clean figure pre and post-processing 

d=tibble()
for (pre in c("classic","twostep")){
  for (post in c("rejection",'loclinear',"neuralnet")){
    d=rbind(d,read.table(paste0("../Data/Step5_Eby_model/RMSE/RMSE_",pre,"_",post,".csv"),sep=";")%>%
              add_column(., Method = paste0(ifelse(pre=="classic","Pre","No_pre"),"_",post)))
  }
}

mean_rmse=d%>%
  melt(., id.vars=c("Method"))%>%
  group_by(.,Method,variable)%>%
  summarise(., .groups = "keep",mean_rmse=mean(value))

p=ggplot(d%>%melt(., id.vars=c("Method")))+
  geom_jitter(aes(x=Method,y=value,color=Method),position = position_jitterdodge(jitter.width = 0.3,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=Method,y=mean_rmse),color="white",fill="black",shape=24,size=2.5)+
  labs(x="",y="NRMSE",color="")+
  facet_wrap(.~variable,labeller = label_bquote(cols= Parameter == .(as.character(variable))))+
  the_theme+
  theme(strip.text.x = element_text(size=12),axis.text.x = element_text(angle=60,hjust=1),legend.position = "none")

ggsave('../Figures/Eby_model/First_ABC/NRMSE_post_and_pre_processing.pdf',p,width = 8,height = 5)


#what values of parameters are better inferred 
d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
d_virtual_param=d_all[1:100,1:2]

d=d%>%
  melt(., id.vars=c("Method"))%>%
  add_column(.,virtual_param=c(rep(d_virtual_param$p,6),rep(d_virtual_param$q,6)))

p=ggplot(d)+
  geom_point(aes(x=virtual_param,y=value,color=Method),alpha=.5)+
  labs(x="",y="NRMSE",color="")+
  facet_wrap(.~variable,labeller = label_bquote(cols= Parameter == .(as.character(variable))))+
  the_theme+
  theme(strip.text.x = element_text(size=12),axis.text.x = element_text(angle=60,hjust=1))

ggsave('../Figures/Eby_model/First_ABC/Quality_inference_NRMSE_virtual_data.pdf',p,width = 8,height = 5)


## >> 3) PCA and variables ----

# PCA parameter and summary statistics in Eby model
d_all=read.table("../Data/All_sim_Eby.csv",sep=";")

sample_row=sample(1:nrow(d_all),size=10000,replace = F)
res.pca=PCA(d_all[sample_row,3:ncol(d_all)], scale.unit = T, ncp = 3,  graph=F)

for (param in 1:2){
  axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))
  for (i in 1:3){
    assign(paste0("p",i,"_",param),
           
           fviz_pca_biplot(res.pca, geom.ind = "point", 
                           axes=c(axes_for_plot$x[i],axes_for_plot$y[i]), col.ind = d_all[sample_row,param],col.var="black",
                           label = "var", repel = T)+
             scale_color_gradientn(colours = 
                                     colorRampPalette(colors = 
                                                        c("#960261","#C5218B","#E266B6",
                                                          "#9F94EF","#2F48DA","#2F48DA"))(100),
             )+
             labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color=colnames(d_all)[param])+
             ggtitle(paste0(colnames(d_all)[param]))+
             theme_classic()+theme(legend.position = "bottom",plot.title = element_text(size=25))
           
    )
  }
  
}

ggsave("../Figures/Eby_model/ACP_and_model_description/ACP_param_Eby.pdf",
       ggarrange(
       ggarrange(p1_1+ggtitle(TeX("$\\p$")),
                 p2_1+ggtitle(TeX("$\\p$")),
                 p3_1+ggtitle(TeX("$\\p$")),ncol=3),
       ggarrange(p1_2+ggtitle(TeX("$\\q$")),
                 p2_2+ggtitle(TeX("$\\q$")),
                 p3_2+ggtitle(TeX("$\\q$")),ncol=3),nrow=2),width = 12,height = 7)


# Comparing coverage property of Eby model and Kefi one

d_all_eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%
  filter(., rho_p !=0)
d_all_kefi=read.table("../Data/All_sim_ABC.csv",sep=";")%>%
  filter(., rho_p > .05 & rho_p < .9)

n_keep=1000
d_tot=rbind(d_all_eby[sample(1:nrow(d_all_eby),size=n_keep,replace = F),3:11]%>%mutate(., model="Eby model"),
            d_all_kefi[sample(1:nrow(d_all_kefi),size=n_keep,replace = F),8:16]%>%mutate(., model="Kefi model"))

res.pca=PCA(d_tot[,-ncol(d_tot)], scale.unit = T, ncp = 3,  graph=F)
axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

for (i in 1:3){
  
  assign(paste0("p",i),
         
         fviz_pca_biplot(res.pca, geom.ind = "point", 
                         axes=c(axes_for_plot$x[i],axes_for_plot$y[i]), col.ind = d_tot$model,col.var="black",alpha=.5,
                         label = "var", repel = T,pointsize =1)+
           scale_color_manual(values=c("#A078DA","#68B77A"))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="")+
           ggtitle("")+guides(shape="none")+
           theme_classic()+theme(legend.position = "bottom",plot.title = element_text(size=25))
         
  )
}

p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
            p2+theme(legend.position = "none"),
            p3+theme(legend.position = "none"),
            ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.1))

ggsave("../Figures/Eby_model/ACP_and_model_description/Coverage_Eby_Kefi_models.pdf",
       p,width=9,height = 4)






#***********************************************************


#***********************************************************

# ---------------------- Step 9: Improving inference ----

#***********************************************************

## >> 1) Optimizing pre and post-processing methods ----
#we play on both the lambda of the boxcox method and the size of the sample for
#stage 1 of the two step preprocessing procedure

d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0 | d_all$PLR ==1) #we remove very healthy sites for which PLR = 1 to avoid errors.
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
matrix_param=d_all[,1:2]
matrix_sumstat=d_all[,3:(ncol(d_all))]

N_for_cross_validation = 100
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)

for (optim_lambda in c(T,F)){

  for (size_step1 in c(1000,3000)){
  
    for (method_abc in c("loclinear","neuralnet")){
      
      for (preprocessing in c("PLS_BoxCox","BoxCox","None")){
        
        mat_cor_param=array(0,c(2,2,N_for_cross_validation)) #correlation matrix for parameters
        
        pdf(paste0("../Figures/Eby_model/Optimizing_inferrence/Pre_post/Cross_validation_n",N_for_cross_validation,"_",preprocessing,
                   "_",method_abc,"_optim_lambda_",ifelse(optim_lambda,"yes","no"),"_N1_",size_step1,".pdf"),width = 8,height = 4)
        
        d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
        
        for (n in 1:N_for_cross_validation){
          
          matrix_param=d_all[,1:2]
          matrix_sumstat=d_all[,3:(ncol(d_all))]
          save_sumstat=matrix_sumstat
          
          n_cross=nrow_for_sample[n]
          
          if (preprocessing %in% c("BoxCox", "PLS_BoxCox")){ #Applying the two step procedure used in Siren MEE paper : Don't know whether it make sense in our case. TO discuss Monday
            
            if (optim_lambda==T){
              for (x in 1:ncol(matrix_sumstat)) if (x %in% c(4,6)){
                
                b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
                lambda_x=b$x[which.max(b$y)]
                if (lambda_x !=0){ #to avoid errors
                  matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
                }
                
              }else {
                b=boxcox(lm(matrix_sumstat[,x] ~ 1),plotit = F,eps = .05)    
                lambda_x=b$x[which.max(b$y)]
                if (lambda_x !=0){ #to avoid errors
                  matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
                }
              }

            } else {
              for (x in 1:ncol(matrix_sumstat)) if (x %in% c(4,6)){
                matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(.5)) -1)/(.5)
              }else {matrix_sumstat[,x] = (matrix_sumstat[,x]^(.5) -1)/(.5)}
            }
            
            #Second we scale
            for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
            
            if (preprocessing=="PLS_BoxCox"){
              #and finally, we perform the first PLS
              
              pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                         data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
              
              
              n_comp_pls=selectNcomp(pls_1,method = "onesigma")
              
              if (n_comp_pls > 1){
                mat_sumstat_pls=pls_1$scores[,1:n_comp_pls] # selecting # components
              } else if (n_comp_pls==1){ #otherwise we take the whole components
                mat_sumstat_pls=as.data.frame(matrix(pls_1$scores[,1:n_comp_pls],ncol=1))
              } else {mat_sumstat_pls=pls_1$scores[,1:ncol(pls_1$scores)]}
              
              
            } else {mat_sumstat_pls=matrix_sumstat}
            
            cross_valid1=abc(target = mat_sumstat_pls[n_cross,],
                            param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                            tol = size_step1/nrow(matrix_param),method = "rejection") #we keep the size_step1 closest simulations for the first step
            
            #Keeping size_step1 simulations and doing the same steps again: normality, scaling and PLS
            
            mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid1$ss)),3:(ncol(d_all))] #we keep information with the true values
            mat_sumstat_step1=rbind(mat_sumstat_step1,d_all[n_cross,3:(ncol(d_all))])
            
            #again, first box cox
            
            
            if (optim_lambda==T){
              for (x in 1:ncol(mat_sumstat_step1)) if (x %in% c(4,6)){
                
                b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
                lambda_x=b$x[which.max(b$y)]
                
                if (lambda_x !=0){ #to avoid errors
                  mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
                }
                
                
              }else {
                b=boxcox(lm(mat_sumstat_step1[,x] ~ 1),plotit = F,eps = .05)    
                lambda_x=b$x[which.max(b$y)]
                if (lambda_x !=0){ #to avoid errors
                  mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
                }
                
              }
              
            } else {
              for (x in 1:ncol(mat_sumstat_step1)) if (x %in% c(4,6)){
                mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(.5)) -1)/(.5)
              }else {mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(.5) -1)/(.5)}
            }
            
            #and normalization
            for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
            
            
            if (preprocessing=="PLS_BoxCox"){
              
              pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                         data=as.data.frame(cbind(rbind(cross_valid1$unadj.values,matrix_param[n_cross,]),
                                                  mat_sumstat_step1)), scale=TRUE, validation="CV")
              n_comp_pls=selectNcomp(pls_2,method = "onesigma")
              
              
              if (n_comp_pls > 1){
                mat_sumstat_pls2=pls_2$scores[,1:n_comp_pls] #pls 2 selecting # components
              } else if (n_comp_pls==1){ #otherwise we take the whole components
                mat_sumstat_pls2=matrix(pls_2$scores[,1:n_comp_pls],ncol=1)
              } else {mat_sumstat_pls2=pls_2$scores[,1:ncol(pls_2$scores)]}
              
            } else {
              
              mat_sumstat_pls2=mat_sumstat_step1
              
            }
            
          
            cross_valid2=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),],
                            param = cross_valid1$unadj.values,
                            sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),], #removing the target data
                            tol = 100/nrow(mat_sumstat_pls2),method = method_abc,transf = rep("logit",2), #as parameters are proba, we perform logit regression
                            logit.bounds = matrix(c(0,1),2,2,byrow = T)) 
            
            cross_valid2$ss=d_all[as.numeric(rownames(cross_valid2$ss)),3:(ncol(d_all))] #we keep information with the true values
            
            
          } else {

            #for each virtual data, we perform ABC rejection algorithm with linear regression adjustment for posterior
            cross_valid2=abc(target = matrix_sumstat[n_cross,],
                            param = matrix_param[-n_cross,],sumstat = matrix_sumstat[-n_cross,], #removing the target data
                            tol = 100/nrow(matrix_param),method = method_abc,transf = rep("logit",2), #as parameters are proba, we perform logit regression
                            logit.bounds = matrix(c(0,1),2,2,byrow = T))
          }    
          
          
          matrix_sumstat=save_sumstat
          
          
          if (names(cross_valid2)[1]=="unadj.values")names(cross_valid2)[1] = "adj.values"
          
          cross_valid2$adj.values=cross_valid2$adj.values
          #Matrix of correlation between parameters & sumstats
          mat_cor_param[,,n]=cor(cross_valid2$adj.values)
          
          
          #We plot the differences in posterior distribution/true parameter
          
          
          par(mfrow=c(1,2))
          for (i in colnames(matrix_param)){
            plot(density(cross_valid2$adj.values[,i]),main=i,xlab="Value")
            abline(v = matrix_param[n_cross,i],col="blue")
            abline(v = colMeans(cross_valid2$adj.values)[i],col="red")
          }
          
          
          par(mfrow=c(1,1))
          for (i in 1:(length(colnames(cross_valid2$adj.values))-1)){
            for (j in (i+1):length(colnames(cross_valid2$adj.values))){
              plot(x=cross_valid2$adj.values[,i],y=cross_valid2$adj.values[,j],
                   xlab=colnames(cross_valid2$adj.values)[i],ylab=colnames(cross_valid2$adj.values)[j],
                   col=alpha("blue",.8))
            }
          }
          
          
          d_melt=as.data.frame(cross_valid2$ss)%>%
            melt(.)
          
          par(mfrow=c(2,5))
          for (i in colnames(matrix_sumstat)){
            plot(density(d_melt$value[which(d_melt$variable==i)]),main=i,xlab="Value")
            abline(v = matrix_sumstat[n_cross,i],col="blue")
            abline(v = colMeans(cross_valid2$ss)[i],col="red")
          }
          
          #we save the mean posterior distribution for each and the true observed parameters
          d_cross_param=rbind(d_cross_param,as_tibble(t(colMeans(cross_valid2$adj.values)))%>%add_column(., Type="Sim"))
          d_cross_param=rbind(d_cross_param,as_tibble((matrix_param[n_cross,]))%>%add_column(., Type="Obs"))
          
          #As we work with virtual data, we do the same for the summary stats we save the mean posterior distribution for each and the true observed parameters
          d_cross_sumstat=rbind(d_cross_sumstat,as_tibble(t(colMeans(cross_valid2$ss)))%>%add_column(., Type="Sim"))
          d_cross_sumstat=rbind(d_cross_sumstat,as_tibble((matrix_sumstat[n_cross,]))%>%add_column(., Type="Obs"))
          
          
          #We compute the mean squared error (RMSE) 
          RMSE = sapply(1:ncol(cross_valid2$adj.values),function(x){
            sqrt(sum((cross_valid2$adj.values[,x]-matrix_param[n_cross,x])**2)/nrow(cross_valid2$adj.values) )
          }
          )
          
          #normalize it by the RMSE under the prior distribution
          RMSE_prior=sapply(1:(ncol(matrix_param)),function(x){
            sqrt(sum((matrix_param[,x]-matrix_param[n_cross,x])**2)/nrow(matrix_param) )
          }
          )
          NRMSE = RMSE/RMSE_prior
          
          d_NRMSE_param=rbind(d_NRMSE_param,as_tibble(t(NRMSE)))
          
          
          #We repeat the same for the summary statistics observed
          RMSE = sapply(1:ncol(cross_valid2$ss),function(x){
            sqrt(sum((cross_valid2$ss[,x]-matrix_sumstat[n_cross,x])**2)/nrow(cross_valid2$ss) )
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
        
        write.table(d_NRMSE_param,paste0("../Data/Step6_Optimizing_inferrence/Pre_post/RMSE_param_",preprocessing,
                                         "_",method_abc,"_optim_lambda_",ifelse(optim_lambda,"yes","no"),"_N1_",size_step1,".csv"),sep=";")
        
        colnames(d_NRMSE_sumstat)=colnames(d_cross_sumstat)
        
        write.table(d_NRMSE_sumstat,paste0("../Data/Step6_Optimizing_inferrence/Pre_post/RMSE_sumstat_",preprocessing,
                                         "_",method_abc,"_optim_lambda_",ifelse(optim_lambda,"yes","no"),"_N1_",size_step1,".csv"),sep=";")
        

        dev.off()
        
        

        
      }
    }
  }
}


# Ploting figure size

all_sim=expand.grid(N1=c(1000,3000),
                    lambda=c("yes","no"),
                    Preproc=c("PLS_BoxCox","BoxCox","None"),
                    postproc=c("loclinear","neuralnet"))

d=tibble()
for (i in 1:nrow(all_sim)){
  d=rbind(d,read.table(paste0("../Data/Step6_Optimizing_inferrence/Pre_post/RMSE_param_",all_sim$Preproc[i],"_",all_sim$postproc[i],"_optim_lambda_",
                              all_sim$lambda[i],"_N1_",all_sim$N1[i],".csv"),sep=";")%>%
            add_column(., N1=all_sim$N1[i],optim_lambda=all_sim$lambda[i],Post=all_sim$postproc[i],Pre=all_sim$Preproc[i]))
}


mean_rmse=d%>%
  melt(., id.vars=c("N1","optim_lambda","Post","Pre"))%>%
  group_by(.,variable,N1,optim_lambda,Post,Pre)%>%
  summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  rename(., "Parameter"="variable")

p=ggplot(d%>%melt(., id.vars=c("N1","optim_lambda","Post","Pre"))%>%
           rename(., "Parameter"="variable"))+
  geom_jitter(aes(x=interaction(Pre,Post),y=value,color=interaction(Pre,Post)),
              position = position_jitterdodge(jitter.width = 0.3,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=interaction(Pre,Post),y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="",y="NRMSE",color="")+
  facet_grid(N1+Parameter~optim_lambda,labeller = label_both)+
  the_theme+
  ylim(0,.5)+
  theme(strip.text.x = element_text(size=10),axis.text.x = element_text(angle=60,hjust=1),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898AA4","#E4C035"))

ggsave(paste0("../Figures/Eby_model/Optimizing_inferrence/Pre_post/Optimization_inference_all.pdf"),p,width = 7,height = 9)


p=d%>%
  melt(., id.vars=c("N1","optim_lambda","Post","Pre"))%>%
  group_by(.,variable,N1,optim_lambda,Post,Pre)%>%
  summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  ggplot(.)+
  geom_line(aes(x=N1,y=mean_rmse,color=interaction(Pre,Post),group=interaction(Pre,Post,optim_lambda),linetype=optim_lambda),lwd=.8)+
  geom_point(aes(x=N1,y=mean_rmse,fill=interaction(Pre,Post)),color="black",shape=21,size=1)+
  facet_wrap(.~variable,labeller = label_bquote(cols="Parameter = "==.(as.character(variable))))+
  labs(x="# virtual data kept stage 1 pre-processing",color="",linetype=TeX(r"(MLE for $\lambda$)"),y="Mean NRMSE")+
  scale_x_continuous(breaks = c(1000,3000))+
  the_theme+
  guides(color = guide_legend(override.aes = list(size = 2)),fill="none")+
  theme(legend.box = "vertical")+
  scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898AA4","#E4C035"))


ggsave("../Figures/Eby_model/Optimizing_inferrence/Pre_post/Optimization_inference_size_step1_preproc.pdf",p,width = 7,height = 4)



 ###XXXX do the same for sumstats





## >> 2) Optimizing the structure of the neural-network ----



d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
matrix_param=d_all[,1:2]
matrix_sumstat=d_all[,3:(ncol(d_all))]

N_for_cross_validation = 100
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)

for (method_pre in c("PLS","NoPLS")){
  
  for (size_hidden in seq(5,25,by=5)){
    
    for (rep_network in c(10,20,30)){
      
      mat_cor_param=array(0,c(2,2,N_for_cross_validation)) #correlation matrix for parameters
      
      
      d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
      
      for (n in 1:N_for_cross_validation){
        
        matrix_param=d_all[,1:2]
        matrix_sumstat=d_all[,3:(ncol(d_all))]
        save_sumstat=matrix_sumstat
        
        n_cross=nrow_for_sample[n]
        
        for (x in 1:ncol(matrix_sumstat)) if (x %in% c(4,6)){
          
          b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
          }
          
        }else {
          b=boxcox(lm(matrix_sumstat[,x] ~ 1),plotit = F,eps = .05)
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
          }
        }
        
        
        #Second we scale
        for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
        
        #and finally, we perform the first PLS
        if (method_pre=="PLS"){
          
          pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                     data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
          n_comp_pls=selectNcomp(pls_1,method = "onesigma")
          
          
          if (n_comp_pls > 1){
            mat_sumstat_pls=pls_1$scores[,1:n_comp_pls] # selecting # components
          } else if (n_comp_pls==1){ #otherwise we take the whole components
            mat_sumstat_pls=matrix(pls_1$scores[,1:n_comp_pls],ncol=1)
          } else {mat_sumstat_pls=pls_1$scores[,1:ncol(pls_1$scores)]}
          
          
        } else {
          mat_sumstat_pls=matrix_sumstat
        }
        
        
        
        cross_valid=abc(target = mat_sumstat_pls[n_cross,],
                        param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                        tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
        
        #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
        
        mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid$ss)),3:(ncol(d_all))] #we keep information with the true values
        mat_sumstat_step1=rbind(mat_sumstat_step1,d_all[n_cross,3:(ncol(d_all))])
        
        #again, first box cox
        
        for (x in 1:ncol(mat_sumstat_step1)) if (x %in% c(4,6)){
          
          b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
          lambda_x=b$x[which.max(b$y)]
          
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
          }
          
          
        }else {
          b=boxcox(lm(mat_sumstat_step1[,x] ~ 1),plotit = F,eps = .05)
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
          }
          
        }
        
        #and normalization
        for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
        
        
        
        if (method_pre=="PLS"){
          
          pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                     data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                              mat_sumstat_step1)), scale=TRUE, validation="CV")
          
          
          n_comp_pls=selectNcomp(pls_2,method = "onesigma")
          
          
          if (n_comp_pls > 1){
            mat_sumstat_pls2=pls_2$scores[,1:n_comp_pls] #pls 2 selecting # components
          } else if (n_comp_pls==1){ #otherwise we take the whole components
            mat_sumstat_pls2=matrix(pls_2$scores[,1:n_comp_pls],ncol=1)
          } else {mat_sumstat_pls2=pls_2$scores[,1:ncol(pls_2$scores)]}
          
        } else {
          mat_sumstat_pls2=mat_sumstat_step1
        }
        
        
        cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),],
                        param = cross_valid$unadj.values,
                        sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),], #removing the target data
                        tol = 100/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                        logit.bounds = matrix(c(0,1),2,2,byrow = T),
                        numnet = rep_network,sizenet = size_hidden)
        
        cross_valid$ss=d_all[as.numeric(rownames(cross_valid$ss)),3:(ncol(d_all))] #we keep information with the true values
        
        matrix_sumstat=save_sumstat
        
        if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
        
        cross_valid$adj.values=cross_valid$adj.values
        #Matrix of correlation between parameters & sumstats
        mat_cor_param[,,n]=cor(cross_valid$adj.values)
        
        
        #We plot the differences in posterior distribution/true parameter
        
        if (any( is.nan(cross_valid$adj.values[,1]) | is.nan(cross_valid$adj.values[,2]))){
          cross_valid$adj.values = cross_valid$adj.values[-which(is.nan(cross_valid$adj.values[,1])
                                                                 | is.nan(cross_valid$adj.values[,2])),]
          
        }
        
        
        
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
        RMSE_prior=sapply(1:(ncol(matrix_param)),function(x){
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
      
      write.table(d_NRMSE_param,paste0("../Data/Step6_Optimizing_inferrence/Neural_net/RMSE_hidden_preprocessing_",method_pre,"_",
                                       size_hidden,"_Nnet_",rep_network,".csv"),sep=";")
      
      
      
    }
  }
}





all_sim=expand.grid(rep_network=seq(10,30,by=10),N_hidden=seq(5,25,by=5))

for (plot_id in 1:2){
  d=tibble()
  
  for (i in 1:nrow(all_sim)){
    d=rbind(d,read.table(paste0("../Data/Step6_Optimizing_inferrence/Neural_net/RMSE_hidden_preprocessing_",c("NoPLS","PLS")[plot_id],"_",
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
    facet_grid(Parameter~N_rep_net,labeller = label_both)+
    the_theme+
    ylim(0,.5)+
    theme(strip.text.x = element_text(size=10),legend.position = "bottom")+
    guides(color = guide_legend(override.aes = list(size = 2)))+
    scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"),breaks=c('5', '10', '15',"20","25"))
  
  ggsave(paste0("../Figures/Eby_model/Optimizing_inferrence/Neural_net/Optimization_NN_all_",c("NoPLS","PLS")[plot_id],".pdf"),
         p,width = 7,height = 4)
  
  
  p=d%>%
    melt(., id.vars=c("N_hidden","N_rep_net"))%>%
    group_by(.,variable,N_hidden,N_rep_net)%>%
    summarise(., .groups = "keep",mean_rmse=mean(value))%>%
    mutate(., N_rep_net=as.character(N_rep_net))%>%
    ggplot(.)+
    geom_line(aes(x=N_hidden,y=mean_rmse,color=N_rep_net,group=N_rep_net),lwd=.8)+
    geom_point(aes(x=N_hidden,y=mean_rmse,fill=N_rep_net),color="black",shape=21,size=1)+
    facet_wrap(.~variable,labeller = label_bquote(cols="Parameter = "==.(as.character(variable))),scales = "free")+
    labs(x="Number hidden neurons",color="# neural networks  ",y="Mean NRMSE")+
    scale_x_continuous(breaks = seq(5,25,by=5))+
    the_theme+
    guides(color = guide_legend(override.aes = list(size = 2)),fill="none")+
    scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5"))+
    scale_fill_manual(values=c("#C46FC5","#80BD5C","#568DC5"))
  
  ggsave(paste0("../Figures/Eby_model/Optimizing_inferrence/Neural_net/Optimization_NN_mean_",c("NoPLS","PLS")[plot_id],".pdf"),
         p,width = 7,height = 4)
  
}



## >> 3) Optimizing the number of principal components ----



d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
matrix_param=d_all[,1:2]
matrix_sumstat=d_all[,3:(ncol(d_all))]

N_for_cross_validation = 100
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)

for (n1_PLS in 1:9){
  
  for (n2_PLS in 1:9){
    
    mat_cor_param=array(0,c(2,2,N_for_cross_validation)) #correlation matrix for parameters
    
    pdf(paste0("../Figures/Eby_model/Optimizing_inferrence/Pls_comp/Cross_validation_n1_PLS_",n1_PLS,
               "_n2_PLS_",n2_PLS,".pdf"),width = 8,height = 4)
    
    d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
    
    for (n in 1:N_for_cross_validation){
      
      matrix_param=d_all[,1:2]
      matrix_sumstat=d_all[,3:(ncol(d_all))]
      save_sumstat=matrix_sumstat
      
      n_cross=nrow_for_sample[n]
      
      for (x in 1:ncol(matrix_sumstat)) if (x %in% c(4,6)){
        
        b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
        lambda_x=b$x[which.max(b$y)]
        if (lambda_x !=0){ #to avoid errors
          matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
        }
        
      }else {
        b=boxcox(lm(matrix_sumstat[,x] ~ 1),plotit = F,eps = .05)    
        lambda_x=b$x[which.max(b$y)]
        if (lambda_x !=0){ #to avoid errors
          matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
        }
      }
      
      
      #Second we scale
      for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
      
      #and finally, we perform the first PLS
      
      pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
      
      
      n_comp_pls=n1_PLS
      
      
      if (n_comp_pls > 1){
        mat_sumstat_pls=pls_1$scores[,1:n_comp_pls] # selecting # components
      } else if (n_comp_pls==1){ #otherwise we take the whole components
        mat_sumstat_pls=matrix(pls_1$scores[,1:n_comp_pls],ncol=1)
      } else {mat_sumstat_pls=pls_1$scores[,1:ncol(pls_1$scores)]}
      
      
      cross_valid=abc(target = mat_sumstat_pls[n_cross,],
                      param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
      
      #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
      
      mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid$ss)),3:(ncol(d_all))] #we keep information with the true values
      mat_sumstat_step1=rbind(mat_sumstat_step1,d_all[n_cross,3:(ncol(d_all))])
      
      #again, first box cox
      
      for (x in 1:ncol(mat_sumstat_step1)) if (x %in% c(4,6)){
        
        b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
        lambda_x=b$x[which.max(b$y)]
        
        if (lambda_x !=0){ #to avoid errors
          mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
        }
        
        
      }else {
        b=boxcox(lm(mat_sumstat_step1[,x] ~ 1),plotit = F,eps = .05)    
        lambda_x=b$x[which.max(b$y)]
        if (lambda_x !=0){ #to avoid errors
          mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
        }
        
      }
      
      #and normalization
      for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
      
      pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                          mat_sumstat_step1)), scale=TRUE, validation="CV")
      
      
      n_comp_pls=n2_PLS
      
      
      if (n_comp_pls > 1){
        mat_sumstat_pls2=pls_2$scores[,1:n_comp_pls] #pls 2 selecting # components
      } else if (n_comp_pls==1){ #otherwise we take the whole components
        mat_sumstat_pls2=matrix(pls_2$scores[,1:n_comp_pls],ncol=1)
      } else {mat_sumstat_pls2=pls_2$scores[,1:ncol(pls_2$scores)]}
      
      cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),], #removing the target data
                      tol = 100/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),2,2,byrow = T),
                      numnet = 10,sizenet = 10) 
      
      cross_valid$ss=d_all[as.numeric(rownames(cross_valid$ss)),3:(ncol(d_all))] #we keep information with the true values
      
      matrix_sumstat=save_sumstat
      
      if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
      
      cross_valid$adj.values=cross_valid$adj.values
      #Matrix of correlation between parameters & sumstats
      mat_cor_param[,,n]=cor(cross_valid$adj.values)
      
      
      #We plot the differences in posterior distribution/true parameter
      
      if (any( is.nan(cross_valid$adj.values[,1]) | is.nan(cross_valid$adj.values[,2]))){
        cross_valid$adj.values = cross_valid$adj.values[-which(is.nan(cross_valid$adj.values[,1])
                                                               | is.nan(cross_valid$adj.values[,2])),]
      }      
     
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
      RMSE_prior=sapply(1:(ncol(matrix_param)),function(x){
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
    
    write.table(d_NRMSE_param,paste0("../Data/Step6_Optimizing_inferrence/Pls_comp/RMSE_n1_PLS_",n1_PLS,
                                     "_n2_PLS_",n2_PLS,".csv"),sep=";")
    

    
  }
}


#loading the case where we optimize the number of pls comp
d_optim=read.table("../Data/Step6_Optimizing_inferrence/Neural_net/RMSE_hidden_10_Nnet_10.csv",sep=";")

all_sim=expand.grid(n1_PLS=1:9,n2_PLS=1:9)

d=tibble()
for (i in 1:nrow(all_sim)){
  d=rbind(d,read.table(paste0("../Data/Step6_Optimizing_inferrence/Pls_comp/RMSE_n1_PLS_",
                              all_sim$n1_PLS[i],"_n2_PLS_",all_sim$n2_PLS[i],".csv"),sep=";")%>%
            mutate(., p=(p-d_optim$p),q=(q-d_optim$q))%>%
            add_column(., n1_PLS=all_sim$n1_PLS[i],n2_PLS=all_sim$n2_PLS[i]))
}



mean_rmse=d%>%
  melt(., id.vars=c("n1_PLS","n2_PLS"))%>%
  group_by(.,variable,n1_PLS,n2_PLS)%>%
  summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  rename(., "Parameter"="variable")
  

p=ggplot(mean_rmse)+
  geom_tile(aes(x=n1_PLS,y=n2_PLS,fill=mean_rmse))+
  labs(x="n1 comp PLS",y="n2 comp PLS",color="",fill="Difference with optimized case  ")+
  facet_grid(.~Parameter,labeller = label_both)+
  the_theme+
  theme(strip.text.x = element_text(size=10),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_fill_gradientn(colours = colorRampPalette(c("blue","red"))(100))


ggsave(paste0("../Figures/Eby_model/Optimizing_inferrence/Pls_comp/PLS_comp_all.pdf"),p,width = 7,height = 4)


## >> 4) Sensibility to summary statistics removal ----


d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
matrix_param=d_all[,1:2]
matrix_sumstat=d_all[,3:(ncol(d_all))]

N_for_cross_validation = 100
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)


  
for (type_removal in c("no_removal","PLR","Expo_PL","PLR_explo_PL","Spectral_ratio","Spec_PLR_expo","clustering","nb_neigh","neigh_clust")){
  
  if (type_removal=="no_removal"){
    sumstat_name=colnames(d_all[,3:(ncol(d_all))])
  } else if (type_removal=="PLR"){
    sumstat_name=colnames(d_all[,3:(ncol(d_all))])[-8]    
  } else if (type_removal=="Expo_PL"){
    sumstat_name=colnames(d_all[,3:(ncol(d_all))])[-9]    
  } else if (type_removal=="Spectral_ratio"){
    sumstat_name=colnames(d_all[,3:(ncol(d_all))])[-7]    
  } else if (type_removal=="PLR_explo_PL"){
    sumstat_name=colnames(d_all[,3:(ncol(d_all))])[-c(8:9)]    
  }else if (type_removal=="clustering"){
    sumstat_name=colnames(d_all[,3:(ncol(d_all))])[-c(3)]    
  }else if (type_removal=="nb_neigh"){
    sumstat_name=colnames(d_all[,3:(ncol(d_all))])[-c(2)]    
  }else if (type_removal=="neigh_clust"){
    sumstat_name=colnames(d_all[,3:(ncol(d_all))])[-c(2:3)]    
  } else {
    sumstat_name=colnames(d_all[,3:(ncol(d_all))])[-c(7:9)]  
  }
    
  mat_cor_param=array(0,c(2,2,N_for_cross_validation)) #correlation matrix for parameters
  d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
  
  for (n in 1:N_for_cross_validation){
    
    matrix_param=d_all[,1:2]
    matrix_sumstat=d_all[,which(colnames(d_all) %in% sumstat_name)]
    save_sumstat=matrix_sumstat
    n_cross=nrow_for_sample[n]
    
    for (x in 1:ncol(matrix_sumstat)) if (x %in% which(colnames(matrix_sumstat) %in% c("skewness","moran_I"))){
      
      b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
      }
    }else {
      b=boxcox(lm(matrix_sumstat[,x] ~ 1),plotit = F,eps = .05)    
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
      }
    }
    
    #Second we scale
    for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
    
    #and finally, we perform the first PLS
    if (type_removal=="no_removal"){
      pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
    } else if (type_removal=="PLR"){
      pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PL_expo,
                 data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
    } else if (type_removal=="Expo_PL"){
      pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR,
                 data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
    } else if (type_removal=="Spectral_ratio"){
      pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+PL_expo+PLR,
                 data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
    } else if (type_removal=="clustering"){
      pls_1=plsr(p + q~rho_p+nb_neigh+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
    } else if (type_removal=="nb_neigh"){
      pls_1=plsr(p + q~rho_p+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
    } else if (type_removal=="neigh_clust"){
      pls_1=plsr(p + q~rho_p+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
    } else if (type_removal=="Spec_PLR_expo"){
      pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I,
                 data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
    } else {
      pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio,
                 data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
    }

    n_comp_pls=selectNcomp(pls_1,method = "onesigma")
    
    if (n_comp_pls > 1){
      mat_sumstat_pls=pls_1$scores[,1:n_comp_pls] # selecting # components
    } else if (n_comp_pls==1){ #otherwise we take the whole components
      mat_sumstat_pls=matrix(pls_1$scores[,1:n_comp_pls],ncol=1)
    } else {mat_sumstat_pls=pls_1$scores[,1:ncol(pls_1$scores)]}
    
    cross_valid=abc(target = mat_sumstat_pls[n_cross,],
                    param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                    tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
    
    #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
    mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid$ss)),which(colnames(d_all) %in% sumstat_name)] #we keep information with the true values
    mat_sumstat_step1=rbind(mat_sumstat_step1,d_all[n_cross,which(colnames(d_all) %in% sumstat_name)])
    
    #again, first box cox
    for (x in 1:ncol(mat_sumstat_step1)) if (x %in% which(colnames(matrix_sumstat) %in% c("skewness","moran_I"))){
      
      b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
      }

    }else {
      b=boxcox(lm(mat_sumstat_step1[,x] ~ 1),plotit = F,eps = .05)    
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
      }
    }
    
    #and normalization
    for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
    
    if (type_removal=="no_removal"){
      pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                          mat_sumstat_step1)), scale=TRUE, validation="CV")
    } else if (type_removal=="PLR"){
      pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PL_expo,
                 data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                          mat_sumstat_step1)), scale=TRUE, validation="CV")
    } else if (type_removal=="Expo_PL"){
      pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR,
                 data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                          mat_sumstat_step1)), scale=TRUE, validation="CV")
    } else if (type_removal=="Spectral_ratio"){
      pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+PLR+PL_expo,
                 data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                          mat_sumstat_step1)), scale=TRUE, validation="CV")
    } else if (type_removal=="Spec_PLR_expo"){
      pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I,
                 data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                          mat_sumstat_step1)), scale=TRUE, validation="CV")
    } else if (type_removal=="clustering"){
      pls_2=plsr(p + q~rho_p+nb_neigh+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                          mat_sumstat_step1)), scale=TRUE, validation="CV")
    } else if (type_removal=="nb_neigh"){
      pls_2=plsr(p + q~rho_p+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                          mat_sumstat_step1)), scale=TRUE, validation="CV")
    } else if (type_removal=="neigh_clust"){
      pls_2=plsr(p + q~rho_p+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                          mat_sumstat_step1)), scale=TRUE, validation="CV")
    } else {
      pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio,
                 data=as.data.frame(cbind(rbind(cross_valid$unadj.values,matrix_param[n_cross,]),
                                          mat_sumstat_step1)), scale=TRUE, validation="CV")
    }
    
    n_comp_pls=selectNcomp(pls_1,method = "onesigma")

    if (n_comp_pls > 1){
      mat_sumstat_pls2=pls_2$scores[,1:n_comp_pls] #pls 2 selecting # components
    } else if (n_comp_pls==1){ #otherwise we take the whole components
      mat_sumstat_pls2=matrix(pls_2$scores[,1:n_comp_pls],ncol=1)
    } else {mat_sumstat_pls2=pls_2$scores[,1:ncol(pls_2$scores)]}
    
    cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),],
                    param = cross_valid$unadj.values,
                    sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),], #removing the target data
                    tol = 100/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                    logit.bounds = matrix(c(0,1),2,2,byrow = T),
                    numnet = 10,sizenet = 10) 
    
    cross_valid$ss=d_all[as.numeric(rownames(cross_valid$ss)),which(colnames(d_all) %in% sumstat_name)] #we keep information with the true values
    
    matrix_sumstat=save_sumstat
    
    if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
    
    #Matrix of correlation between parameters & sumstats
    mat_cor_param[,,n]=cor(cross_valid$adj.values)
    
    if (any( is.nan(cross_valid$adj.values[,1]) | is.nan(cross_valid$adj.values[,2]))){
      cross_valid$adj.values = cross_valid$adj.values[-which(is.nan(cross_valid$adj.values[,1])
                                                             | is.nan(cross_valid$adj.values[,2])),]
    }      
    
    #we save the mean posterior distribution for each and the true observed parameters
    d_cross_param=rbind(d_cross_param,as_tibble(t(colMeans(cross_valid$adj.values)))%>%add_column(., Type="Sim"))
    d_cross_param=rbind(d_cross_param,as_tibble((matrix_param[n_cross,]))%>%add_column(., Type="Obs"))
    
    #As we work with virtual data, we do the same for the summary stats we save the mean posterior distribution for each and the true observed parameters
    d_cross_sumstat=rbind(d_cross_sumstat,as_tibble(t(colMeans(cross_valid$ss)))%>%add_column(., Type="Sim"))
    d_cross_sumstat=rbind(d_cross_sumstat,as_tibble((matrix_sumstat[n_cross,]))%>%add_column(., Type="Obs"))
    
    
    #We compute the mean squared error (RMSE) 
    RMSE = sapply(1:ncol(cross_valid$adj.values),function(x){
      sqrt(sum((cross_valid$adj.values[,x]-matrix_param[n_cross,x])**2)/nrow(cross_valid$adj.values))})
    
    #normalize it by the RMSE under the prior distribution
    RMSE_prior=sapply(1:(ncol(matrix_param)),function(x){
      sqrt(sum((matrix_param[,x]-matrix_param[n_cross,x])**2)/nrow(matrix_param) )})
    
    NRMSE = RMSE/RMSE_prior
    d_NRMSE_param=rbind(d_NRMSE_param,as_tibble(t(NRMSE)))
    
    #We repeat the same for the summary statistics observed
    RMSE = sapply(1:ncol(cross_valid$ss),function(x){
      sqrt(sum((cross_valid$ss[,x]-matrix_sumstat[n_cross,x])**2)/nrow(cross_valid$ss) )})
    
    RMSE_prior=sapply(1:ncol(matrix_sumstat),function(x){
      sqrt(sum((matrix_sumstat[,x]-matrix_sumstat[n_cross,x])**2)/nrow(matrix_sumstat) )})
    NRMSE = RMSE/RMSE_prior
    
    d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,as_tibble(t(NRMSE)))
    
  } #end loop Nvirtual data
  
  colnames(d_NRMSE_param)=colnames(d_cross_param)[-length(colnames(d_cross_param))]
  
  write.table(d_NRMSE_param,paste0("../Data/Step6_Optimizing_inferrence/Sensi_removal/RMSE_",type_removal,".csv"),sep=";")

}




d=tibble()
for (i in list.files("../Data/Step6_Optimizing_inferrence/Sensi_removal/")){
  d=rbind(d,read.table(paste0("../Data/Step6_Optimizing_inferrence/Sensi_removal/",i),sep=";")%>%
            add_column(.,Type=gsub(gsub(pattern = ".csv",replacement = "",x=i),pattern = 'RMSE_',replacement = '')))
}



mean_rmse=d%>%
  melt(., id.vars=c("Type"))%>%
  group_by(.,variable,Type)%>%
  summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  rename(., "Parameter"="variable")


p=ggplot(d%>%melt(., id.vars=c("Type"))%>%
           rename(., "Parameter"="variable"))+
  geom_jitter(aes(x=Type,y=value,color=as.factor(Type)),
              position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=Type,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="",y="NRMSE",color="")+
  facet_grid(.~Parameter,labeller = label_both)+
  the_theme+
  ylim(0,.5)+
  theme(strip.text.x = element_text(size=10),legend.position = "bottom",axis.text.x = element_text(angle=60,hjust=1))+
  guides(color = guide_legend(override.aes = list(size = 2)))

ggsave(paste0("../Figures/Eby_model/Optimizing_inferrence/Sensi_removal/Sensi_removal.pdf"),p,width = 7,height = 4)








# ---------------------- Step 10: Empirical data ----

## >> 0) Summary statistics of the data & PCA model/data ----

for (id in 1:5){
  pdf(paste0("../Figures/Empirical_data/Landscapes/all_landscapes_",id,".pdf"),width = 6,height = 6)
  for (land in (((id-1)*69)+1):(id *69 )){
    Plot_empirical(land)
  }
  dev.off()
}




d_biocom=as_tibble(readxl::read_xlsx("../Data/Data_Biocom/biocom-infos.xlsx"))
d_sumstat=read.table("../Data/Data_Biocom/Summary_stats_data.csv",sep=",")
colnames(d_sumstat)=c("rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                      "Spectral_ratio","PLR","PL_expo")

d_biocom=d_biocom%>%
  
  dplyr::select(., file,id, plotn, MF,TypeMF, Aridity,Sand,lat,long,imgcover,nbpixels)%>%
  
  dplyr::rename(., File_ID=file, ID=id, Plot_n=plotn, Lattitude=lat,Longitude=long,Cover=imgcover,Nbpixels=nbpixels)%>%
  
  bind_cols(., d_sumstat) %>%
  
  as.data.frame(.)

write.table(d_biocom,"../Data/Data_Biocom/biocom_data.csv",sep=";")


#Distribution of summary statistics

pdf("../Figures/Empirical_data/Comparizon_models_data/Distrib_sum_stat_data.pdf",width = 7,height = 6)
par(mfrow=c(3,3))
for (sumstat in 1:9){
  hist(d_biocom[,10+sumstat],xlab=colnames(d_biocom)[10+sumstat],main="",col=alpha("blue",.3))
}
dev.off()

#and comparizon with Eby model
d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)
p=ggplot(rbind(d_biocom[,11:ncol(d_biocom)]%>%add_column(.,Type="Data"),d_Eby[,3:ncol(d_Eby)]%>%add_column(., Type="Model"))%>%
           melt(., id.vars=c("Type")))+
    geom_density(aes(x=value,fill=Type),alpha=.7)+
    the_theme+
    facet_wrap(.~variable,scales = "free")+
    scale_fill_manual(values=c("#C0CEAA","#5D589C"))

ggsave("../Figures/Empirical_data/Comparizon_models_data/Density_Eby_and_data.pdf",p,width = 7,height = 6)


# PCA for comparing data and model
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)
d_Kefi=read.table("../Data/All_sim_ABC.csv",sep=";")%>%filter(., rho_p > 0.1  && rho_p < 0.8)
d_tot=rbind(d_Eby[sample(1:nrow(d_Eby),15000,F),3:ncol(d_Eby)]%>%add_column(., Type="Eby model"),
            d_biocom[,11:ncol(d_biocom)]%>%add_column(., Type="zEmpirical sites"), #to plot the empirical sites above simulations
            d_Kefi[sample(1:nrow(d_Kefi),15000,F),8:ncol(d_Kefi)]%>%add_column(., Type="Kefi model"))%>%
  arrange(., Type)

#First raw PCA

sumstat_name=colnames(d_Eby)[3:ncol(d_Eby)]
res.pca=PCA(d_tot[,which(colnames(d_tot) %in% sumstat_name)], scale.unit = T, ncp = 3,  graph=F)
axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

#getting the centroids
centroids=tibble(Type=d_tot$Type)%>%
  mutate(., Type=recode_factor(Type,"zEmpirical sites"="Empirical sites"))%>%
  add_column(., PC1=res.pca$ind$coord[,1],PC2=res.pca$ind$coord[,2],PC3=res.pca$ind$coord[,3])%>%
  group_by(., Type)%>%
  summarise(., .groups = "keep",C1=mean(PC1),C2=mean(PC2),C3=mean(PC3))%>%
  as.data.frame()
  
  
for (i in 1:3){
  assign(paste0("p",i),
         d_tot%>%
           mutate(., Type=recode_factor(Type,"zEmpirical sites"="Empirical sites"))%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_point(aes(x = PC1, y = PC2, color = Type,fill=Type,size=Type),alpha=.5)+
           scale_size_manual(values=c(.5,1,.5))+
           scale_color_manual(values=c("#FD4848",alpha("#C0CEAA",.5),alpha("#A078DA",.5)))+
           scale_fill_manual(values=c("#FD4848",alpha("#C0CEAA",.5),alpha("#A078DA",.5)))+
           geom_point(data=centroids,aes(x=centroids[,axes_for_plot$x[i]+1],y=centroids[,axes_for_plot$y[i]+1]),shape=24,fill="black",color="white")+
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
ggsave("../Figures/Empirical_data/Comparizon_models_data/PCA_coverage_models_data.pdf",p, width=9,height = 4)


#Second, PCA on Eby model or Kefi model only with data
for (model in c("Kefi","Eby")){
  
  sumstat_name=colnames(d_Eby)[3:ncol(d_Eby)]
  dat=as.data.frame(d_tot[-which(d_tot$Type==paste0(model," model")),])
  
  res.pca=PCA(as.data.frame(dat[,which(colnames(dat) %in% sumstat_name)]),
              scale.unit = T, ncp = 3,  graph=F)
  axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))
  
  for (i in 1:3){
    assign(paste0("p",i),
           dat%>%
             mutate(., Type=recode_factor(Type,"zEmpirical sites"="Empirical sites"))%>%
             add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = Type,fill=Type),size=.5,alpha=.5)+
             
             scale_color_manual(values=c("#FD4848",ifelse(model=="Kefi",alpha("#C0CEAA",.5),alpha("#A078DA",.5))))+
             scale_fill_manual(values=c("#FD4848",ifelse(model=="Kefi",alpha("#C0CEAA",.5),alpha("#A078DA",.5))))+
             labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="",fill="")+
             ggtitle("")+guides(shape="none")+
             theme_classic()+theme(legend.position = "bottom")+
             guides(color = guide_legend(override.aes = list(size = 3)),fill="none")  
    )
  }
  
  p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                        p2+theme(legend.position = "none"),
                        p3+theme(legend.position = "none"),
                        ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
              nrow=2,heights = c(1,.1))
  
  ggsave(paste0("../Figures/Empirical_data/Comparizon_models_data/PCA_",ifelse(model=="Kefi","Eby","Kefi"),"_empirical_data.pdf"),p, width=9,height = 4)
  
}




#Knock-out of some summary statistics and projection on PCA simulations + data
pdf("../Figures/Empirical_data/Comparizon_models_data/Knock_out_PCA_Eby_data.pdf", width=9,height = 4)

for (scena_sumstat_removal in c("None","PLR","Expo_PL","PLR_explo_PL","Spectral_ratio","Spec_PLR_expo","clustering","nb_neigh","neigh_clust")){
  
  if (scena_sumstat_removal=="None"){
    sumstat_name=colnames(d_Eby)[c(3:ncol(d_Eby))] 
  } else if (scena_sumstat_removal=="PLR"){
    sumstat_name=colnames(d_Eby)[c(3:ncol(d_Eby))[-8]] 
  } else if (scena_sumstat_removal=="Expo_PL"){
    sumstat_name=colnames(d_Eby)[c(3:ncol(d_Eby))[-9]]    
  } else if (scena_sumstat_removal=="PLR_explo_PL"){
    sumstat_name=colnames(d_Eby)[c(3:ncol(d_Eby))[-c(8:9)]]    
  } else if (scena_sumstat_removal=="Spectral_ratio"){
    sumstat_name=colnames(d_Eby)[c(3:ncol(d_Eby))[-7]]    
  } else if (scena_sumstat_removal=="clustering"){
    sumstat_name=colnames(d_Eby)[c(3:ncol(d_Eby))[-3]]    
  } else if (scena_sumstat_removal=="Spec_PLR_expo"){
    sumstat_name=colnames(d_Eby)[c(3:ncol(d_Eby))[-c(7:9)]]    
  } else if (scena_sumstat_removal=="nb_neigh"){
    sumstat_name=colnames(d_Eby)[c(3:ncol(d_Eby))[-2]]    
  } else {
    sumstat_name=colnames(d_Eby)[c(3:ncol(d_Eby))[-c(2:3)]]    
  }
  
  dat=as.data.frame(d_tot[-which(d_tot$Type=="Kefi model"),])
  

  res.pca=PCA(as.data.frame(dat[,which(colnames(dat) %in% sumstat_name)]),
              scale.unit = T, ncp = 3,  graph=F)
  axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))
  
  for (i in 1:3){
    assign(paste0("p",i),
           dat%>%
             add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = Type,fill=Type),size=.5)+
             
             scale_color_manual(values=c(alpha("#C0CEAA",.5),"#FD4848"))+
             scale_fill_manual(values=c(alpha("#C0CEAA",.5),"#FD4848"))+
             labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="",fill="")+
             ggtitle("")+guides(shape="none")+
             theme_classic()+theme(legend.position = "bottom")+
             guides(color = guide_legend(override.aes = list(size = 3)),fill="none")  
    )
  }
  
  p=ggarrange(ggarrange(p1+theme(legend.position = "none",plot.title=element_blank()),
                        p2+theme(legend.position = "none",plot.title=element_blank()),
                        p3+theme(legend.position = "none",plot.title=element_blank()),
                        ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
              nrow=2,heights = c(1,.1))
  
  print(annotate_figure(p, top = text_grob(paste0("Knock-out = ",scena_sumstat_removal), 
                                              face = "bold", size = 14))
  )

}
dev.off()



## >> 1) Naive ABC on all empirical sites ----

d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()


d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)

for (method_data in c("Raw","NoPLS")){
  
  d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
  
  d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
  d_posterior=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")[as.numeric(which(closest_sites)),]
  

  index=1
  for (site_ID in 1:nrow(d_biocom)){
    
    #observed summary statistics in the site
    observed_sumstat=d_biocom[site_ID,which(colnames(d_biocom) %in% c("rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                                                                      "Spectral_ratio","PLR","PL_expo"))]
    
    if (d_biocom$PL_expo[site_ID] ==0) observed_sumstat=observed_sumstat[-9] #if we cannot fit a PL or tPL, we remove the PL_expo
    
    matrix_param=d_all[,1:2]
    matrix_sumstat=d_all[,which(colnames(d_all) %in% names(observed_sumstat))]
    save_sumstat=matrix_sumstat
    matrix_sumstat=rbind(matrix_sumstat,observed_sumstat)
    

    for (x in 1:ncol(matrix_sumstat)) if (x %in% which(colnames(matrix_sumstat) %in% c("skewness","moran_I"))){
      
      b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
      }
    }else {
      b=boxcox(lm(matrix_sumstat[,x] ~ 1),plotit = F,eps = .05)    
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
      }
    }
    
    #Second we scale
    for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
    
    #and finally, we perform the first PLS (excluding empirical sdata)
    if (method_data != "NoPLS"){
      
      if (d_biocom$PL_expo[site_ID] ==0){
        pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR,
                   data=cbind(matrix_param,matrix_sumstat[-nrow(matrix_sumstat),]), scale=TRUE, validation="CV")
      } else {
        pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                   data=cbind(matrix_param,matrix_sumstat[-nrow(matrix_sumstat),]), scale=TRUE, validation="CV")
      }
      
      n_comp_pls=selectNcomp(pls_1,method = "onesigma")
      
      if (n_comp_pls > 1){
        mat_sumstat_pls=pls_1$scores[,1:n_comp_pls] # selecting # components
      } else if (n_comp_pls==1){ #otherwise we take the whole components
        mat_sumstat_pls=matrix(pls_1$scores[,1:n_comp_pls],ncol=1)
      } else {mat_sumstat_pls=pls_1$scores[,1:ncol(pls_1$scores)]}
      
      #predicting the values of observed summary statistics under the pls model
      observed_sumstat_pls1=predict(pls_1, matrix_sumstat[nrow(matrix_sumstat),], ncomp = 1:n_comp_pls, type = "scores")
      
    } else {
      
      mat_sumstat_pls=matrix_sumstat[-nrow(matrix_sumstat),]
      observed_sumstat_pls1=matrix_sumstat[nrow(matrix_sumstat),]
    }
    

    cross_valid1=abc(target = observed_sumstat_pls1,
                     param = matrix_param,sumstat = mat_sumstat_pls, #removing the target data
                     tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
    
    #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
    mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid1$ss)),which(colnames(d_all) %in% names(observed_sumstat))] #we keep information with the true values
    mat_sumstat_step1=rbind(mat_sumstat_step1,observed_sumstat)
    
    
    #again, first box cox
    for (x in 1:ncol(mat_sumstat_step1)) if (x %in% which(colnames(matrix_sumstat) %in% c("skewness","moran_I"))){
      
      b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
      }
      
    }else {
      b=boxcox(lm(mat_sumstat_step1[,x] ~ 1),plotit = F,eps = .05)    
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
      }
    }
    
    #and normalization
    for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
    
    if (method_data != "NoPLS"){
      
      if (d_biocom$PL_expo[site_ID] ==0){
        pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR,
                   data=as.data.frame(cbind(rbind(cross_valid1$unadj.values),
                                            mat_sumstat_step1[-nrow(mat_sumstat_step1),])), scale=TRUE, validation="CV")  
        
      } else {
        pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                   data=as.data.frame(cbind(rbind(cross_valid1$unadj.values),
                                            mat_sumstat_step1[-nrow(mat_sumstat_step1),])), scale=TRUE, validation="CV")  
      }
      
      n_comp_pls=selectNcomp(pls_2,method = "onesigma")
      
      #predicting the values of observed summary statistics under the pls model
      observed_sumstat_pls2=predict(pls_2, matrix_sumstat[nrow(matrix_sumstat),], ncomp = 1:n_comp_pls, type = "scores")
      
      
      if (n_comp_pls > 1){
        mat_sumstat_pls2=pls_2$scores[,1:n_comp_pls] #pls 2 selecting # components
      } else if (n_comp_pls==1){ #otherwise we take the whole components
        mat_sumstat_pls2=as.data.frame(matrix(pls_2$scores[,1:n_comp_pls],ncol=1))
      } else {mat_sumstat_pls2=pls_2$scores[,1:ncol(pls_2$scores)]}
      
      rownames(mat_sumstat_pls2)=rownames(pls_2$scores)
      
      
      
    } else {
      mat_sumstat_pls2=mat_sumstat_step1[-nrow(mat_sumstat_step1),]
      observed_sumstat_pls2=mat_sumstat_step1[nrow(mat_sumstat_step1),]
    }
    
    
    cross_valid2=abc(target = observed_sumstat_pls2,
                     param = cross_valid1$unadj.values,
                     sumstat = mat_sumstat_pls2, #removing the target data
                     tol = 75/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                     logit.bounds = matrix(c(0,1),2,2,byrow = T),
                     numnet = 10,sizenet = 10) 
    
    cross_valid2$ss=as.data.frame(cross_valid2$ss)
    cross_valid2$ss=d_all[as.numeric(rownames(cross_valid2$ss)),which(colnames(d_all) %in% names(observed_sumstat))] #we keep information with the true values
    

    
    matrix_sumstat=save_sumstat
    
    if (names(cross_valid2)[1]=="unadj.values")names(cross_valid2)[1] = "adj.values"
    
    
    if (any( is.nan(cross_valid2$adj.values[,1]) | is.nan(cross_valid2$adj.values[,2]))){
      cross_valid2$adj.values = cross_valid2$adj.values[-which(is.nan(cross_valid2$adj.values[,1])
                                                               | is.nan(cross_valid2$adj.values[,2])),]
    }      
    

    if (d_biocom$PL_expo[site_ID] ==0){
      d_cross_sumstat=rbind(d_cross_sumstat,as_tibble(t(colMeans(cross_valid2$ss)))%>%
                              add_column(.,PL_expo=NA, Type="Sim",Site=d_biocom$File_ID[site_ID]))
      d_cross_sumstat=rbind(d_cross_sumstat,as_tibble((observed_sumstat))%>%
                              add_column(., PL_expo=NA,Type="Obs",Site=d_biocom$File_ID[site_ID]))
    } else {
      d_cross_sumstat=rbind(d_cross_sumstat,as_tibble(t(colMeans(cross_valid2$ss)))%>%
                              add_column(.,Type="Sim",Site=d_biocom$File_ID[site_ID]))
      d_cross_sumstat=rbind(d_cross_sumstat,as_tibble((observed_sumstat))%>%
                              add_column(., Type="Obs",Site=d_biocom$File_ID[site_ID]))
      
    }
    
    
    
    #We compute the RMSE for the summary statistics observed
    RMSE = sapply(1:ncol(cross_valid2$ss),function(x){
      sqrt(sum((cross_valid2$ss[,x]-as.numeric(observed_sumstat[x]))**2)/nrow(cross_valid2$ss) )})
    
    RMSE_prior=sapply(1:ncol(matrix_sumstat),function(x){
      sqrt(sum((matrix_sumstat[,x]-as.numeric(observed_sumstat[x]))**2)/nrow(matrix_sumstat) )})
    NRMSE = RMSE/RMSE_prior
    names(NRMSE)= colnames(matrix_sumstat)
    
    if (d_biocom$PL_expo[site_ID] ==0){
      d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,as_tibble(t(NRMSE))%>%
                              add_column(., PL_expo = NA,Site=site_ID))
    } else {
      d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,as_tibble(t(NRMSE))%>%
                              add_column(.,Site=site_ID))
    }
    
    
    
    d_posterior[index,20]=mean(cross_valid2$adj.values[,1]) #mean posteriors
    d_posterior[index,21]=mean(cross_valid2$adj.values[,2]) #mean posteriors
    index=index+1
    
  }
  dev.off()
  
  
  d_posterior=d_posterior%>%
    rename(., p=V20,q=V21)
  
  write.table(d_posterior,paste0("../Data/Step7_Empirical_data/ABC_all_sites/Posteriors_sites_",method_data,".csv"),sep=";")
  write.table(d_NRMSE_sumstat,paste0("../Data/Step7_Empirical_data/ABC_all_sites/NRMSE_sumstat_",method_data,".csv"),sep=";")
  write.table(d_cross_sumstat,paste0("../Data/Step7_Empirical_data/ABC_all_sites/Cross_valid_sumstat_",method_data,".csv"),sep=";")
  
}


#NRMSE sumstat
d_NRMSE_sumstat=rbind(
  read.table("../Data/Step7_Empirical_data/ABC_all_sites/NRMSE_sumstat_NoPLS.csv",sep=";")%>%add_column(., Method="No PLS"),
  read.table("../Data/Step7_Empirical_data/ABC_all_sites/NRMSE_sumstat_PLS.csv",sep=";")%>%add_column(., Method="PLS")
  )
p=d_NRMSE_sumstat%>%
  melt(., id.vars=c("Site","Method"))%>%
  ggplot(.,aes(x=Method,y=value))+
  geom_line(aes(group=Site),lwd=.2,color=alpha("gray",.3))+
  geom_violin(aes(fill=variable),alpha=.4)+
  facet_wrap(.~variable,scales = "free")+
  stat_summary(fun = "median", geom = "point", size = 3,shape=24,fill="black",color="white")+
  the_theme+
  labs(x="",y="NRMSE",color="",fill="")+
  geom_hline(yintercept = 1)+
  guides(fill="none") 

ggsave("../Figures/Empirical_data/ABC/NRMSE_sumstat_data.pdf",p,width = 9,height = 6)


#mean simultion summary stat versus observed sumstat
d_cross_sumstat=rbind(
  read.table("../Data/Step7_Empirical_data/ABC_all_sites/Cross_valid_sumstat_NoPLS.csv",sep=";")%>%add_column(., Method="No PLS"),
  read.table("../Data/Step7_Empirical_data/ABC_all_sites/Cross_valid_sumstat_PLS.csv",sep=";")%>%add_column(., Method="PLS")
)

pdf("../Figures/Empirical_data/ABC/x_y_obs_sum_sumstat_data.pdf",width = 7,height = 7)
d_sim=as.data.frame(filter(d_cross_sumstat,Type=="Sim"))
d_obs=as.data.frame(filter(d_cross_sumstat,Type=="Obs"))
par(mfrow=c(3,3))
for (i in 1:9){
  plot(as.numeric(d_sim[,i]),
       as.numeric(d_obs[,i]),col=c(alpha("#C46FC5",.3),alpha("#80BD5C",.3))[as.factor(d_sim$Method)],pch=21,xlab="Sim",ylab="Obs",
       main=colnames(d_cross_sumstat)[i])
  abline(a=0,b=1)
}
dev.off()



d_posterior=tibble()
for (x_file in list.files("../Data/Step7_Empirical_data/ABC_all_sites/",pattern = "Posterior")){
  d_posterior = rbind(
    d_posterior,
    read.table(paste0("../Data/Step7_Empirical_data/ABC_all_sites/",x_file),sep=";")%>%
      add_column(., Type=gsub(".csv","",strsplit(x_file,"_")[[1]][3]))
  )
}

p=ggplot(d_posterior%>%
           melt(., measure.vars=c("p","q"))%>%
           melt(., measure.vars=c("Aridity","MF","Sand"),variable.name="Driver",value.name = "Driver_value"))+
  geom_point(aes(x=Driver_value,y=value,color=Type),alpha=.3)+
  geom_smooth(method = "lm",aes(x=Driver_value,y=value,color=Type),se=F)+
  facet_grid(variable~Driver,scales = "free")+
  labs(x="Value of the driver",y="Inferred parameter value",color="Type of pre-processing")+
  the_theme+
  scale_color_manual(values=c("#C46FC5","#80BD5C"))
ggsave("../Figures/Empirical_data/ABC/Driver_inferred_parameters.pdf",p,width = 7,height = 6)




## >> 2) Linking quality estimation with distance in PCA  ----

d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)

d_tot=rbind(d_Eby[,3:ncol(d_Eby)]%>%add_column(., Type="Eby model"),
            d_biocom[,11:ncol(d_biocom)]%>%add_column(., Type="Empirical sites"))

res.pca=PCA(d_tot[,-ncol(d_tot)], scale.unit = T, ncp = 3,  graph=F)

#centroids in PCA of simulations
centroid_sim=colMeans(res.pca$ind$coord[1:nrow(d_Eby),])

#getting euclidean distance from centroid for each empirical site
dist_empirical = sapply((nrow(d_Eby)+1):(nrow(d_tot)),function(x){
  return( sqrt(sum((res.pca$ind$coord[x,] - centroid_sim)^2)) )
})
write.table(dist_empirical,"../Data/Step7_Empirical_data/Closest_sim/Dist_centroid.csv",sep=";")

#Then we plot NRMSE = f(distance to centroid) 
d_NRMSE_sumstat=read.table("../Data/Step7_Empirical_data/ABC_all_sites/NRMSE_sumstat_NoPLS.csv",sep=";")

p=ggplot(d_NRMSE_sumstat%>%
         melt(., id.vars=c("Site"))%>%
         add_column(., Euclid_dist=rep(dist_empirical,9)))+
  geom_point(aes(x=Euclid_dist,y=value,color=variable),alpha=.4)+
  facet_wrap(.~variable,scales = "free")+
  the_theme+
  theme(legend.position = "none")+
  labs(x="Distance to centroid in PCA",y="NRMSE")

ggsave("../Figures/Empirical_data/ABC/Distance_closest_simu/Dist_centroid_NMRSE_fit.pdf",p,width = 7,height = 6)

#and we do the same with the distance to the x=y line 
d_cross_sumstat=read.table("../Data/Step7_Empirical_data/ABC_all_sites/Cross_valid_sumstat_NoPLS.csv",sep=";")%>%
  melt(., id.vars=c("Type","Site"))

error_abc=sapply(1:(nrow(d_cross_sumstat)/2),function(x){
  return(
    abs(d_cross_sumstat$value[2*(x-1)+1]-d_cross_sumstat$value[2*x] ) #absolute distance between true and fitted 
  )
})

p=d_cross_sumstat%>%
  filter(., Type=="Obs")%>%
  dplyr::select(., -Type)%>%
  add_column(., Error_abc=error_abc,Euclid_dist=rep(dist_empirical,9))%>%
  ggplot(.)+
  geom_point(aes(x=Euclid_dist,y=Error_abc,color=variable),alpha=.4)+
  facet_wrap(.~variable,scales = "free")+
  the_theme+
  theme(legend.position = "none")+
  labs(x="Distance to centroid in PCA",y="Absolute difference between fitted and observed ")

ggsave("../Figures/Empirical_data/ABC/Distance_closest_simu/Dist_centroid_ABC_sumstat_fitted.pdf",p,width = 7,height = 6)


# Similarly for the link between aridity and inferred p values

#doing the same with the smallest distance to a simulation point
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)

d_tot=rbind(d_Eby[,3:ncol(d_Eby)]%>%add_column(., Type="Eby model"),
            d_biocom[,11:ncol(d_biocom)]%>%add_column(., Type="Empirical sites"))

res.pca=PCA(d_tot[,-ncol(d_tot)], scale.unit = T, ncp = 3,  graph=F)

distances=apply(res.pca$ind$coord[(nrow(res.pca$ind$coord)-344):(nrow(res.pca$ind$coord)),],1, #for each empirical data
                   function(x) min(apply(res.pca$ind$coord[1:(nrow(res.pca$ind$coord)-345),],1,
                                         function(y) dist(rbind(x, y))))) #we get the min distance to simulated data


write.table(distances,"../Data/Step7_Empirical_data/Closest_sim/Distance_closest_point.csv",sep=";")

#Then we plot NRMSE = f(distance to centroid) 
distances=read.table("../Data/Step7_Empirical_data/Closest_sim/Distance_closest_point.csv",sep=";")$x
d_NRMSE_sumstat=read.table("../Data/Step7_Empirical_data/ABC_all_sites/NRMSE_sumstat_NoPLS.csv",sep=";")

p=ggplot(d_NRMSE_sumstat%>%
           melt(., id.vars=c("Site"))%>%
           add_column(., Euclid_dist=rep(distances,9)))+
  geom_point(aes(x=Euclid_dist,y=value,color=variable),alpha=.4)+
  facet_wrap(.~variable,scales = "free")+
  the_theme+
  theme(legend.position = "none")+
  labs(x="Distance closest simulation point",y="NRMSE")

ggsave("../Figures/Empirical_data/ABC/Distance_closest_simu/Dist_simupoint_NMRSE_fit.pdf",p,width = 7,height = 6)

#and we do the same with the distance to the x=y line 
d_cross_sumstat=read.table("../Data/Step7_Empirical_data/ABC_all_sites/Cross_valid_sumstat_NoPLS.csv",sep=";")%>%
  melt(., id.vars=c("Type","Site"))

error_abc=sapply(1:(nrow(d_cross_sumstat)/2),function(x){
  return(
    abs(d_cross_sumstat$value[2*(x-1)+1]-d_cross_sumstat$value[2*x] ) #absolute distance between true and fitted 
  )
})

p=d_cross_sumstat%>%
  filter(., Type=="Obs")%>%
  dplyr::select(., -Type)%>%
  add_column(., Error_abc=error_abc,Euclid_dist=rep(distances,9))%>%
  ggplot(.)+
  geom_point(aes(x=Euclid_dist,y=Error_abc,color=variable),alpha=.4)+
  facet_wrap(.~variable,scales = "free")+
  the_theme+
  theme(legend.position = "none")+
  labs(x="Distance closest simulation point",y="Absolute difference between fitted and observed ")

ggsave("../Figures/Empirical_data/ABC/Distance_closest_simu/Dist_simupoint_ABC_sumstat_fitted.pdf",p,width = 7,height = 6)









#example of sites for which Eby model cannot do much
pdf("../Figures/Empirical_data/landscapes/Example_for_PDE.pdf",width = 20,height = 4)
par(mfrow=c(1,5))
Plot_empirical(10)
Plot_empirical(12)
Plot_empirical(283)
Plot_empirical(296)
Plot_empirical(345)
dev.off()



## >> 3) ABC with some summary stat removal ----


d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)

for (preprocessing in c("NOPLS","PLS")){

  for (remove_ID in c("neighB","clustering","spectralR","moranI","variance")){
    
    d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
    d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
  
    
    if (remove_ID == "neighB"){
      sumstat_kept=colnames(d_all)[3:ncol(d_all)][-2]
    } else if (remove_ID == "clustering"){
      sumstat_kept=colnames(d_all)[3:ncol(d_all)][-3]
    } else if (remove_ID == "spectralR"){
      sumstat_kept=colnames(d_all)[3:ncol(d_all)][-7]
    } else if (remove_ID == "moranI"){
      sumstat_kept=colnames(d_all)[3:ncol(d_all)][-6]
    } else {
      sumstat_kept=colnames(d_all)[3:ncol(d_all)][-5]
    }
    
    for (site_ID in 1:nrow(d_biocom)){
    
      #observed summary statistics in the site
      observed_sumstat=d_biocom[site_ID,which(colnames(d_biocom) %in% sumstat_kept)]
  
      if (d_biocom$PL_expo[site_ID] ==0) observed_sumstat=observed_sumstat[-ncol(observed_sumstat)] #if we cannot fit a PL or tPL, we remove the PL_expo
      
      matrix_param=d_all[,1:2]
      matrix_sumstat=d_all[,which(colnames(d_all) %in% names(observed_sumstat))]
      save_sumstat=matrix_sumstat
      matrix_sumstat=rbind(matrix_sumstat,observed_sumstat)
      
      for (x in 1:ncol(matrix_sumstat)) if (x %in% which(colnames(matrix_sumstat) %in% c("skewness","moran_I"))){
        
        b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
        lambda_x=b$x[which.max(b$y)]
        if (lambda_x !=0){ #to avoid errors
          matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
        }
      }else {
        b=boxcox(lm(matrix_sumstat[,x]+.01 ~ 1),plotit = F,eps = .05)    
        lambda_x=b$x[which.max(b$y)]
        if (lambda_x !=0){ #to avoid errors
          matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
        }
      }
      
      #Second we scale
      for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
      
      #and finally, we perform the first PLS (excluding empirical sdata)
      print(as.formula(paste(paste(colnames(d_all)[1:2],collapse =  " + "),
                             " ~ ",
                             paste(sumstat_kept,collapse =  " + "))))
      
      if (preprocessing=="PLS"){
        
        if (d_biocom$PL_expo[site_ID] ==0){
          
          pls_1=plsr(as.formula(paste(paste(colnames(d_all)[1:2],collapse =  " + "),
                                      " ~ ",
                                      paste(sumstat_kept[-which(sumstat_kept == "PL_expo")],collapse =  " + "))),
                     data=cbind(matrix_param,matrix_sumstat[-nrow(matrix_sumstat),]), scale=TRUE, validation="CV")
        } else {
          pls_1=plsr(as.formula(paste(paste(colnames(d_all)[1:2],collapse =  " + "),
                                      " ~ ",
                                      paste(sumstat_kept,collapse =  " + "))),
                     data=cbind(matrix_param,matrix_sumstat[-nrow(matrix_sumstat),]), scale=TRUE, validation="CV")
        }
        
        n_comp_pls=selectNcomp(pls_1,method = "onesigma")
        
        if (n_comp_pls > 1){
          mat_sumstat_pls=pls_1$scores[,1:n_comp_pls] # selecting # components
        } else if (n_comp_pls==1){ #otherwise we take the whole components
          mat_sumstat_pls=as.data.frame(matrix(pls_1$scores[,1:n_comp_pls],ncol=1))
        } else {mat_sumstat_pls=pls_1$scores[,1:ncol(pls_1$scores)]}
        
        #predicting the values of observed summary statistics under the pls model
        observed_sumstat_pls1=predict(pls_1, matrix_sumstat[nrow(matrix_sumstat),], ncomp = 1:n_comp_pls, type = "scores")
        
      } else {
        
        mat_sumstat_pls=matrix_sumstat[-nrow(matrix_sumstat),]
        observed_sumstat_pls1=matrix_sumstat[nrow(matrix_sumstat),]
      }
  
      cross_valid1=abc(target = observed_sumstat_pls1,
                      param = matrix_param,sumstat = mat_sumstat_pls, #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
      
      #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
      mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid1$ss)),which(colnames(d_all) %in% names(observed_sumstat))] #we keep information with the true values
      mat_sumstat_step1=rbind(mat_sumstat_step1,observed_sumstat)
      
      #again, first box cox
      for (x in 1:ncol(mat_sumstat_step1)) if (x %in% which(colnames(matrix_sumstat) %in% c("skewness","moran_I"))){
        
        b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
        lambda_x=b$x[which.max(b$y)]
        if (lambda_x !=0){ #to avoid errors
          mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
        }
        
      }else {
        b=boxcox(lm(mat_sumstat_step1[,x]+.01 ~ 1),plotit = F,eps = .05)    
        lambda_x=b$x[which.max(b$y)]
        if (lambda_x !=0){ #to avoid errors
          mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
        }
      }
      
      #and normalization
      for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
      
      if (preprocessing=="PLS"){
        
        if (d_biocom$PL_expo[site_ID] ==0){
          pls_2=plsr(as.formula(paste(paste(colnames(d_all)[1:2],collapse =  " + "),
                                      " ~ ",
                                      paste(sumstat_kept[-which(sumstat_kept == "PL_expo")],collapse =  " + "))),
                     data=as.data.frame(cbind(rbind(cross_valid1$unadj.values),
                                              mat_sumstat_step1[-nrow(mat_sumstat_step1),])), scale=TRUE, validation="CV")
        } else {
          pls_2=plsr(as.formula(paste(paste(colnames(d_all)[1:2],collapse =  " + "),
                                      " ~ ",
                                      paste(sumstat_kept,collapse =  " + "))),
                     data=as.data.frame(cbind(rbind(cross_valid1$unadj.values),
                                              mat_sumstat_step1[-nrow(mat_sumstat_step1),])), scale=TRUE, validation="CV")
        }
        
        n_comp_pls=selectNcomp(pls_2,method = "onesigma")
        
        #predicting the values of observed summary statistics under the pls model
        observed_sumstat_pls2=predict(pls_2, matrix_sumstat[nrow(matrix_sumstat),], ncomp = 1:n_comp_pls, type = "scores")
        
        if (n_comp_pls > 1){
          mat_sumstat_pls2=pls_2$scores[,1:n_comp_pls] #pls 2 selecting # components
        } else if (n_comp_pls==1){ #otherwise we take the whole components
          mat_sumstat_pls2=as.data.frame(matrix(pls_2$scores[,1:n_comp_pls],ncol=1))
        } else {mat_sumstat_pls2=pls_2$scores[,1:ncol(pls_2$scores)]}
        
        rownames(mat_sumstat_pls2)=rownames(pls_2$scores)
        
      } else {
        mat_sumstat_pls2=mat_sumstat_step1[-nrow(mat_sumstat_step1),]
        observed_sumstat_pls2=mat_sumstat_step1[nrow(mat_sumstat_step1),]
      }
      
      cross_valid2=abc(target = observed_sumstat_pls2,
                      param = cross_valid1$unadj.values,
                      sumstat = mat_sumstat_pls2, #removing the target data
                      tol = 75/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),2,2,byrow = T),
                      numnet = 10,sizenet = 10) 
      
      cross_valid2$ss=d_all[as.numeric(rownames(cross_valid2$ss)),which(colnames(d_all) %in% names(observed_sumstat))] #we keep information with the true values
      
      
      matrix_sumstat=save_sumstat
      
      if (names(cross_valid2)[1]=="unadj.values")names(cross_valid2)[1] = "adj.values"
      
      
      if (any( is.nan(cross_valid2$adj.values[,1]) | is.nan(cross_valid2$adj.values[,2]))){
        cross_valid2$adj.values = cross_valid2$adj.values[-which(is.nan(cross_valid2$adj.values[,1])
                                                               | is.nan(cross_valid2$adj.values[,2])),]
      }      
      
      if (d_biocom$PL_expo[site_ID] ==0){
        d_cross_sumstat=rbind(d_cross_sumstat,as_tibble(t(colMeans(cross_valid2$ss)))%>%
                                add_column(.,PL_expo=NA, Type="Sim",Site=d_biocom$File_ID[site_ID]))
        d_cross_sumstat=rbind(d_cross_sumstat,as_tibble((observed_sumstat))%>%
                                add_column(., PL_expo=NA,Type="Obs",Site=d_biocom$File_ID[site_ID]))
      } else {
        d_cross_sumstat=rbind(d_cross_sumstat,as_tibble(t(colMeans(cross_valid2$ss)))%>%
                                add_column(.,Type="Sim",Site=d_biocom$File_ID[site_ID]))
        d_cross_sumstat=rbind(d_cross_sumstat,as_tibble((observed_sumstat))%>%
                                add_column(., Type="Obs",Site=d_biocom$File_ID[site_ID]))
      }
      
      
      
      #We compute the RMSE for the summary statistics observed
      RMSE = sapply(1:ncol(cross_valid2$ss),function(x){
        sqrt(sum((cross_valid2$ss[,x]-as.numeric(observed_sumstat[x]))**2)/nrow(cross_valid2$ss) )})
      
      RMSE_prior=sapply(1:ncol(matrix_sumstat),function(x){
        sqrt(sum((matrix_sumstat[,x]-as.numeric(observed_sumstat[x]))**2)/nrow(matrix_sumstat) )})
      NRMSE = RMSE/RMSE_prior
      names(NRMSE)= colnames(matrix_sumstat)
      
      if (d_biocom$PL_expo[site_ID] ==0){
        d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,as_tibble(t(NRMSE))%>%
                                add_column(., PL_expo = NA,Site=site_ID))
      } else {
        d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,as_tibble(t(NRMSE))%>%
                                add_column(.,Site=site_ID))
      }
      
      
      d_biocom[site_ID,20]=mean(cross_valid2$adj.values[,1]) #mean posteriors
      d_biocom[site_ID,21]=mean(cross_valid2$adj.values[,2]) #mean posteriors
      
      
    }

    d_biocom=d_biocom%>%
      rename(., p=V20,q=V21)
    
    write.table(d_biocom,paste0("../Data/Step7_Empirical_data/Removal_sumary_stat/Posteriors_sites_",remove_ID,"_",preprocessing,".csv"),sep=";")
    write.table(d_NRMSE_sumstat,paste0("../Data/Step7_Empirical_data/Removal_sumary_stat/NRMSE_sumstat_",remove_ID,"_",preprocessing,".csv"),sep=";")
    write.table(d_cross_sumstat,paste0("../Data/Step7_Empirical_data/Removal_sumary_stat/Cross_valid_sumstat_",remove_ID,"_",preprocessing,".csv"),sep=";")
  }
}

for (method_data in c("PLS","NOPLS")){

  d_removal=tibble()
  for (x_file in list.files("../Data/Step7_Empirical_data/Removal_sumary_stat/",pattern = "NRMSE")[
    grep(paste0("_",method_data),list.files("../Data/Step7_Empirical_data/Removal_sumary_stat/",pattern = "NRMSE"))]){ #only selecting NoPLS or PLS
    
    table=read.table(paste0("../Data/Step7_Empirical_data/Removal_sumary_stat/",x_file),sep=";")%>%
      add_column(., Removed =strsplit(split = "_",x_file)[[1]][3],
                 Other_colum=NA)
    
    if (strsplit(split = "_",x_file)[[1]][3]=="moranI") {colnames(table)[ncol(table)]="moran_I"
      } else if (strsplit(split = "_",x_file)[[1]][3]=="neighB"){ colnames(table)[ncol(table)]="nb_neigh"
    } else if (strsplit(split = "_",x_file)[[1]][3]=="spectralR"){ colnames(table)[ncol(table)]="Spectral_ratio"
    } else {colnames(table)[ncol(table)]=strsplit(split = "_",x_file)[[1]][3]}
    
    #all tables have the same structure despite variable removal
    table=table%>%
      dplyr::select(.,rho_p,nb_neigh,clustering,skewness,moran_I,variance,
             Spectral_ratio,PLR,PL_expo,Site,Removed)
    
    d_removal=rbind(d_removal,table)
  }
  d_removal=rbind(d_removal,
                  read.table(paste0("../Data/Step7_Empirical_data/ABC_all_sites/NRMSE_sumstat_",method_data,".csv"),sep=";")%>%
                    add_column(., Removed="None")
                    )
  
  p=d_removal%>%
      melt(., id.vars=c("Site","Removed"))%>%
      ggplot(.,aes(x=Removed,y=value))+
      geom_line(aes(group=Site),lwd=.2,color=alpha("gray",.3))+
      geom_violin(aes(fill=Removed),alpha=.4)+
      facet_wrap(.~variable,scales="free")+
      stat_summary(fun = "median", geom = "point", size = 3,shape=24,fill="black",color="white")+
      the_theme+
      labs(x="",y="NRMSE",color="",fill="")+
      geom_hline(yintercept = 1)+
      guides(fill="none")+
    theme(axis.text.x = element_text(hjust=1,angle=60))
  
  ggsave(paste0("../Figures/Empirical_data/ABC/Sensi_removal/NRMSE_removal_sumstat_",method_data,".pdf"),p,width =14,height = 10 )
}



for (method_data in c("NOPLS")){
  
  pdf(paste0("../Figures/Empirical_data/ABC/Sensi_removal/Observed_retained_x_y_",method_data,".pdf"),width = 7,height = 6)

  d_removal=tibble()
  for (x_file in list.files("../Data/Step7_Empirical_data/Removal_sumary_stat/",pattern = "Cross_valid")[
    grep(paste0("_",method_data),list.files("../Data/Step7_Empirical_data/Removal_sumary_stat/",pattern = "Cross_valid"))
  ]){
    
    table=read.table(paste0("../Data/Step7_Empirical_data/Removal_sumary_stat/",x_file),sep=";")%>%
      add_column(., Removed =strsplit(split = "_",x_file)[[1]][4],
                 Other_colum=NA)
    
    if (strsplit(split = "_",x_file)[[1]][4]=="moranI") {colnames(table)[ncol(table)]="moran_I"
    } else if (strsplit(split = "_",x_file)[[1]][4]=="neighB"){ colnames(table)[ncol(table)]="nb_neigh"
    } else if (strsplit(split = "_",x_file)[[1]][4]=="spectralR"){ colnames(table)[ncol(table)]="Spectral_ratio"
    } else {colnames(table)[ncol(table)]=strsplit(split = "_",x_file)[[1]][4]}
    
    #all tables have the same structure despite variable removal
    table=table%>%
      dplyr::select(.,rho_p,nb_neigh,clustering,skewness,moran_I,variance,
                    Spectral_ratio,PLR,PL_expo,Site,Removed,Type)
    
    d_removal=rbind(d_removal,table)
  }
  d_removal=rbind(d_removal,
                  read.table(paste0("../Data/Step7_Empirical_data/ABC_all_sites/Cross_valid_sumstat_",method_data,".csv"),sep=";")%>%
                    add_column(., Removed="None"))%>%
    melt(.,id.vars=c("Site","Removed","Type"))
  

  sapply(1:length(unique(d_removal$variable)),function(y){
    par(mfrow=c(2,3))
    sapply(1:length(unique(d_removal$Removed)),function(x){
      if (any(filter(d_removal,Type=="Sim",!is.na(value),Removed==unique(d_removal$Removed)[x],variable==unique(d_removal$variable)[y])$value)){
        plot(x=filter(d_removal,Type=="Sim",!is.na(value),Removed==unique(d_removal$Removed)[x],variable==unique(d_removal$variable)[y])$value,
             y=filter(d_removal,Type=="Obs",!is.na(value),Removed==unique(d_removal$Removed)[x],variable==unique(d_removal$variable)[y])$value,
             col=c(alpha("blue",.3)),pch=21,xlab="Sim",ylab="Obs",
             main=paste0("Removed = ",unique(d_removal$Removed)[x]))
        abline(a=0,b=1)
        
      }
    })    } )         
  
  
  dev.off()
}





## >> 4) Filtering the landscapes which are closer to the simulations in PCA ----

#first we visualize whether the metrics we use are relevant in the PCA space: this is good

distances=read.table("../Data/Step7_Empirical_data/Closest_sim/Distance_closest_point.csv",sep=";")$x
dist_empirical=read.table("../Data/Step7_Empirical_data/Closest_sim/Dist_centroid.csv",sep=";")$x
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)%>%
  dplyr::select(., -p, -q)
d_tot=rbind(d_biocom[,11:ncol(d_biocom)],d_Eby)

res.pca=PCA(d_tot, scale.unit = T, ncp = 3,  graph=F)
projec_space=res.pca$ind$coord[1:nrow(d_biocom),]%>%
  as.data.frame(.)%>%
  add_column(., Dist_to_sim=distances+dist_empirical)

projec_space$Closest_sim =  projec_space$Dist_to_sim <=  quantile(projec_space$Dist_to_sim,probs = .05) 

axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

pdf("../Figures/Empirical_data/ABC/Distance_closest_simu/Linking_PCA_and_distance.pdf",width = 4,height = 3)
for (i in 1:3){
  print(
    projec_space%>%
      ggplot(.) +
      geom_hline(yintercept = 0, lty = 2) +
      geom_vline(xintercept = 0, lty = 2) +
      geom_point(aes(x = projec_space[,axes_for_plot$x[i]], y = projec_space[,axes_for_plot$y[i]], 
                     color = Closest_sim,
                     fill=Closest_sim))+
      scale_color_viridis_d()+
      scale_fill_viridis_d()+
      labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
           y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="Closest sims",fill="")+
      ggtitle("")+guides(shape="none")+
      theme_classic()+theme(legend.position = "bottom")
  )
}

dev.off()


# Then we keep the simulations which are the closer to simulations and centroids

#First we look at their characteristics
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")

p=d_biocom%>%
  add_column(., 
             Distance_centroids=distances,
             Dist_point_sim=dist_empirical,
             Sum_dist=dist_empirical+distances)%>%
  add_column(., 
             closest_sites = sapply(1:nrow(.),function(x){return(.$Sum_dist[x]  <=  quantile(.$Sum_dist,probs = .05) )}) )%>%
  melt(., id.vars=c("File_ID","ID","Plot_n","Nbpixels",
                    "Distance_centroids","Dist_point_sim","Sum_dist","closest_sites"))%>%
  ggplot(.)+
  labs(fill="Closest site ?",x="")+
  geom_density(aes(x=value,fill=closest_sites),alpha=.3)+
  the_theme+
  facet_wrap(.~variable,scales = "free")

ggsave("../Figures/Empirical_data/ABC/Filtering_data/Density_closest_vs_others.pdf",width = 7,height = 6)


d_NRMSE_sumstat=read.table("../Data/Step7_Empirical_data/ABC_all_sites/NRMSE_sumstat_NoPLS.csv",sep=";")
p=d_NRMSE_sumstat%>%
  add_column(., Distance_centroids=distances,Dist_point_sim=dist_empirical,Sum_dist=dist_empirical+distances)%>%
  add_column(., closest_sites = sapply(1:nrow(.),function(x){return(.$Sum_dist[x]  <=  quantile(.$Sum_dist,probs = .05) )}) )%>% #we keep the 5% closest with the simulations.
  melt(., id.vars=c("closest_sites","Site","Distance_centroids","Sum_dist","Dist_point_sim"))%>%
  ggplot(.,aes(x=variable,y=value))+
  geom_line(aes(group=Site,color=closest_sites),lwd=.2)+
  geom_violin(aes(fill=variable),alpha=.4)+
  stat_summary(fun = "median", geom = "point", size = 3,shape=24,fill="black",color="white")+
  the_theme+
  labs(x="",y="NRMSE",color="Closest sites ?",fill="")+
  geom_hline(yintercept = 1)+
  scale_color_manual(values=c(alpha("gray",.3),"red"))+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  theme(legend.box = "vertical")

ggsave("../Figures/Empirical_data/ABC/Filtering_data/NRMSE_sumstat_data_filtered.pdf",p,width = 9,height = 6)




d_cross_sumstat=read.table("../Data/Step7_Empirical_data/ABC_all_sites/Cross_valid_sumstat_NoPLS.csv",sep=";")
pdf("../Figures/Empirical_data/ABC/Filtering_data/x_y_obs_sum_sumstat_data_filtered.pdf",width = 7,height = 7)

#Getting the closest sites compared to simulations
closest_sites = sapply(1:nrow(d_biocom),function(x){return((distances+dist_empirical)[x]  <=  quantile(distances+dist_empirical,probs = .1) )})
d_sim=as.data.frame(filter(d_cross_sumstat,Type=="Sim"))
d_obs=as.data.frame(filter(d_cross_sumstat,Type=="Obs"))
par(mfrow=c(3,3))

for (i in 1:9){
  plot(as.numeric(d_sim[,i]),bg=c(alpha("gray",.4),"red")[as.factor(closest_sites)],col ="transparent",
       as.numeric(d_obs[,i]),pch=21,xlab="Sim",ylab="Obs",cex=c(.7,1)[as.factor(closest_sites)],
       main=colnames(d_cross_sumstat)[i])
  abline(a=0,b=1)
}
dev.off()



#images of the closest landscapes
pdf("../Figures/Empirical_data/ABC/Filtering_data/Closest_landscapes.pdf",width = 5,height = 5)
for (i in which(closest_sites)){
  Plot_empirical(i)
}
dev.off()



d_posterior=read.table(paste0("../Data/Step7_Empirical_data/ABC_all_sites/Posteriors_sites_NoPLS.csv"),sep=";")
p=ggplot(d_posterior%>%
           add_column(., 
                      Distance_centroids=distances,
                      Dist_point_sim=dist_empirical,
                      Sum_dist=dist_empirical+distances)%>%
           add_column(., 
                      closest_sites = sapply(1:nrow(.),function(x){return(.$Sum_dist[x]  <=  quantile(.$Sum_dist,probs = .1) )}) )%>%
           melt(., measure.vars=c("p","q"))%>%
           melt(., measure.vars=c("Aridity","MF","Sand"),variable.name="Driver",value.name = "Driver_value"))+
  geom_point(aes(x=Driver_value,y=value,color=closest_sites),alpha=.3)+
  geom_smooth(method = "lm",aes(x=Driver_value,y=value,color=closest_sites),se=F)+
  facet_grid(variable~Driver,scales = "free")+
  labs(x="Value of the driver",y="Inferred parameter value",color="Closest sites ?")+
  the_theme+
  scale_color_manual(values=c("#ADAAAA","#E45050"))
ggsave("../Figures/Empirical_data/ABC/Filtering_data/Driver_inferred_parameters_filtered.pdf",p,width = 7,height = 6)


## >> 5) Filtering the sites which perform the best ----

d_NRMSE_sumstat=read.table("../Data/Step7_Empirical_data/ABC_all_sites/NRMSE_sumstat_NoPLS.csv",sep=";")%>%
  add_column(., Sum_NRMSE = sapply(1:nrow(.),function(x){
    return(rowSums(.[x,1:9],na.rm = T))
  }))%>%
  arrange(., Sum_NRMSE,decreasing=T)
closest_sites = sapply(1:nrow(d_biocom),function(x){return((distances+dist_empirical)[x]  <=  quantile(distances+dist_empirical,probs = .1) )})

as.numeric(which(closest_sites))
best_sites=sort(d_NRMSE_sumstat$Site[1:length(as.numeric(which(closest_sites)))])
best_sites
#Not exactly the same, let's do the same analysis with this filtered sites.



d_NRMSE_sumstat=read.table("../Data/Step7_Empirical_data/ABC_all_sites/NRMSE_sumstat_NoPLS.csv",sep=";")
p=d_NRMSE_sumstat%>%
  add_column(., Best_sites = sapply(1:nrow(.),function(x){return(ifelse (x %in% best_sites,T,F))}) )%>% 
  melt(., id.vars=c("Best_sites","Site"))%>%
  ggplot(.,aes(x=variable,y=value))+
  geom_line(aes(group=Site,color=Best_sites),lwd=.2)+
  geom_violin(aes(fill=variable),alpha=.4)+
  stat_summary(fun = "median", geom = "point", size = 3,shape=24,fill="black",color="white")+
  the_theme+
  labs(x="",y="NRMSE",color="Best_sites ?",fill="")+
  geom_hline(yintercept = 1)+
  scale_color_manual(values=c(alpha("gray",.3),"red"))+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  theme(legend.box = "vertical")

ggsave("../Figures/Empirical_data/ABC/Filtering_data/NRMSE_sumstat_data_best_sites.pdf",p,width = 9,height = 6)



# ---------------------- Step 11: Using posteriors to predict the distance to the tipping point --------

## >> 1) Getting the posterior characteristics of the best and closest sites ----
#We recompute briefly the posteriror of parameters while keeping an IC around the mean posterior

d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
d_posterior=tibble()


d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
closest_sites=read.table("../Data/Step7_Empirical_data/Closest_sites.csv",sep=";")


  
index=1

#we keep the sites that fit the best and that are the closest to the simulations
for (site_ID in which(1:nrow(d_biocom) %in% closest_sites$Closest & 1:nrow(d_biocom) %in% closest_sites$Best)){
  
  #observed summary statistics in the site
  observed_sumstat=d_biocom[site_ID,which(colnames(d_biocom) %in% c("rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                                                                    "Spectral_ratio","PLR","PL_expo"))]
  
  if (d_biocom$PL_expo[site_ID] ==0) observed_sumstat=observed_sumstat[-9] #if we cannot fit a PL or tPL, we remove the PL_expo
  
  matrix_param=d_all[,1:2]
  matrix_sumstat=d_all[,which(colnames(d_all) %in% names(observed_sumstat))]
  save_sumstat=matrix_sumstat
  matrix_sumstat=rbind(matrix_sumstat,observed_sumstat)

  for (x in 1:ncol(matrix_sumstat)) if (x %in% which(colnames(matrix_sumstat) %in% c("skewness","moran_I"))){
    
    b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
    lambda_x=b$x[which.max(b$y)]
    if (lambda_x !=0){ #to avoid errors
      matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
    }
  }else {
    b=boxcox(lm(matrix_sumstat[,x] ~ 1),plotit = F,eps = .05)    
    lambda_x=b$x[which.max(b$y)]
    if (lambda_x !=0){ #to avoid errors
      matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
    }
  }
  
  #Second we scale
  for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
  
  mat_sumstat_pls=matrix_sumstat[-nrow(matrix_sumstat),]
  observed_sumstat_pls1=matrix_sumstat[nrow(matrix_sumstat),]
  
  cross_valid1=abc(target = observed_sumstat_pls1,
                   param = matrix_param,sumstat = mat_sumstat_pls, #removing the target data
                   tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
  
  #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
  mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid1$ss)),which(colnames(d_all) %in% names(observed_sumstat))] #we keep information with the true values
  mat_sumstat_step1=rbind(mat_sumstat_step1,observed_sumstat)

  #again, first box cox
  for (x in 1:ncol(mat_sumstat_step1)) if (x %in% which(colnames(matrix_sumstat) %in% c("skewness","moran_I"))){
    
    b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
    lambda_x=b$x[which.max(b$y)]
    if (lambda_x !=0){ #to avoid errors
      mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
    }
    
  }else {
    b=boxcox(lm(mat_sumstat_step1[,x] ~ 1),plotit = F,eps = .05)    
    lambda_x=b$x[which.max(b$y)]
    if (lambda_x !=0){ #to avoid errors
      mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
    }
  }
  #and normalization
  for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
  
  mat_sumstat_pls2=mat_sumstat_step1[-nrow(mat_sumstat_step1),]
  observed_sumstat_pls2=mat_sumstat_step1[nrow(mat_sumstat_step1),]

  cross_valid2=abc(target = observed_sumstat_pls2,
                   param = cross_valid1$unadj.values,
                   sumstat = mat_sumstat_pls2, #removing the target data
                   tol = 75/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                   logit.bounds = matrix(c(0,1),2,2,byrow = T),
                   numnet = 10,sizenet = 10) 
  
  cross_valid2$ss=as.data.frame(cross_valid2$ss)
  cross_valid2$ss=d_all[as.numeric(rownames(cross_valid2$ss)),which(colnames(d_all) %in% names(observed_sumstat))] #we keep information with the true values
  
  matrix_sumstat=save_sumstat
  
  if (names(cross_valid2)[1]=="unadj.values")names(cross_valid2)[1] = "adj.values"
  
  
  if (any( is.nan(cross_valid2$adj.values[,1]) | is.nan(cross_valid2$adj.values[,2]))){
    cross_valid2$adj.values = cross_valid2$adj.values[-which(is.nan(cross_valid2$adj.values[,1])
                                                             | is.nan(cross_valid2$adj.values[,2])),]
  }      
  
  #we keep the mean as well as the quantiles at 5 and 95% for IC
  d_posterior=rbind(d_posterior,tibble(Mean_p=mean(cross_valid2$adj.values[,1]),Mean_q=mean(cross_valid2$adj.values[,2]),
                                           p_q05=quantile(cross_valid2$adj.values[,1],probs = .05),q_q05=quantile(cross_valid2$adj.values[,2],probs = .05),
                                           p_q95=quantile(cross_valid2$adj.values[,1],probs = .95),q_q95=quantile(cross_valid2$adj.values[,2],probs = .95)))
  index=index+1
}


write.table(d_posterior,"../Data/Step7_Empirical_data/Tipping_point/Posterior_best_sites.csv",sep=";")


## >> 2) Computing a distance to the tipping point ----











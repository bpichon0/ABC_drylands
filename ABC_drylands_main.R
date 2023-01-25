rm(list=ls())
source("./ABC_drylands_function.R")

#***********************************************************

# Step 1 : Generating parameters for sensitivity analysis ----

#***********************************************************

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

# Step 2 : Cross verification -----

#***********************************************************

## 1) Generating pseudo-parameters using latin hypercube sampling ----

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


## 2) Analysis ----

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

ggsave(paste0("../Figures/Cross_validation/Correlation_sumstats.pdf"),width = 5,height = 5)


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
  
  pdf(paste0("../Figures/Cross_validation/Cross_validation_n",N_for_cross_validation,"_",method_abc,".pdf"),width = 8,height = 4)
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
  
  ggsave(paste0("../Figures/Cross_validation/Correlation_parameters_",method_abc,".pdf"),width = 6,height = 5)
  
  
  #Ploting f(x,y) with x=True parameter/summary stat, y=simulated
  ##d-melting the tibble
  
  pdf(paste0("../Figures/Cross_validation/x_y_obs_true_param_",method_abc,".pdf"),width = 8,height = 5)
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
  
  
  
  pdf(paste0("../Figures/Cross_validation/x_y_obs_true_summarystat_",method_abc,".pdf"),width = 10,height = 7)
  
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



## 3) PCA and variables ----


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
             labs(x=paste0("PC 1 (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC 2 (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color=colnames(d_all)[param])+
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

# Step 3: Influence of the number of photos to average ----

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
    
    pdf(paste0("../Figures/Number_pictures/x_y_obs_true_param_",method_abc,"_nbpic_",nb_pic,".pdf"),width = 8,height = 5)
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
    
    
    
    pdf(paste0("../Figures/Number_pictures/x_y_obs_true_summarystat_",method_abc,"_nbpic_",nb_pic,".pdf"),width = 10,height = 7)
    
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
    
    ggsave(paste0("../Figures/Number_pictures/NRMSE_param_",method_abc,"_nbpic_",nb_pic,".pdf"),p,width = 8,height = 5)
    
    
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
    
    ggsave(paste0("../Figures/Number_pictures/NRMSE_summarystat_",method_abc,"_nbpic_",nb_pic,".pdf"),p,width = 8,height = 5)
    
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
  
ggsave("../Figures/Number_pictures/Influence_#_pictures.pdf",p,width = 6,height = 4)








#***********************************************************

# Step 4: Fixing some parameters, varying others ----

#***********************************************************


## Pseudo-parameters ----
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


## Analysis ----


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
    
    ggsave(paste0("../Figures/Combination_param/NRMSE_param_",virtual_exp,"_",method_abc,".pdf"),p,width = 8,height = 5)
    
    
    
    pdf(paste0("../Figures/Combination_param/x_y_obs_true_param_",virtual_exp,"_",method_abc,".pdf"),width = 8,height = 5)
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
    
    ggsave(paste0("../Figures/Combination_param/Correlation_parameters_",virtual_exp,"_",method_abc,".pdf"),width = 6,height = 5)
    
    
  }
  
}












#***********************************************************

# Step 5: Testing two step procedure proposed by Siren et al., 2019 ----

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
        mat_sumstat_pls=pls_1$Yscores[,1:n_comp_pls] # selecting # components
      } else if (n_comp_pls==1){ #otherwise we take the whole components
        mat_sumstat_pls=matrix(pls_1$Yscores[,1:n_comp_pls],ncol=1)
      } else {mat_sumstat_pls=pls_1$Yscores[,1:ncol(pls_1$Yscores)]}
      
      
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
        mat_sumstat_pls2=pls_2$Yscores[,1:n_comp_pls] #pls 2 selecting # components
      } else if (n_comp_pls==1){ #otherwise we take the whole components
        mat_sumstat_pls2=matrix(pls_2$Yscores[,1:n_comp_pls],ncol=1)
      } else {mat_sumstat_pls2=pls_2$Yscores[,1:ncol(pls_2$Yscores)]}
      
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
  
  ggsave(paste0("../Figures/Siren_two_steps/Correlation_parameters_",ifelse(two_step,"twostep","classic"),".pdf"),width = 6,height = 5)
  
  
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
  
  ggsave(paste0("../Figures/Siren_two_steps/NRMSE_param_",ifelse(two_step,"twostep","classic"),".pdf"),p,width = 8,height = 5)
  
  
}# loop over the two method for ABC





#***********************************************************

# Step 6: Adding post-processing with linear regression ----

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
      mat_sumstat_pls=pls_1$Yscores[,1:n_comp_pls] # selecting # components
    } else if (n_comp_pls==1){ #otherwise we take the whole components
      mat_sumstat_pls=matrix(pls_1$Yscores[,1:n_comp_pls],ncol=1)
    } else {mat_sumstat_pls=pls_1$Yscores[,1:ncol(pls_1$Yscores)]}
    
    
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
      mat_sumstat_pls2=pls_2$Yscores[,1:n_comp_pls] #pls 2 selecting # components
    } else if (n_comp_pls==1){ #otherwise we take the whole components
      mat_sumstat_pls2=matrix(pls_2$Yscores[,1:n_comp_pls],ncol=1)
    } else {mat_sumstat_pls2=pls_2$Yscores[,1:ncol(pls_2$Yscores)]}
    
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
  
  ggsave(paste0("../Figures/Adding_postprocessing/Correlation_parameters_",method_abc,"post_processing.pdf"),width = 6,height = 5)
  
  
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
  
  ggsave(paste0("../Figures/Adding_postprocessing/NRMSE_param_",method_abc,".pdf"),p,width = 8,height = 5)
  
}







#***********************************************************

# Step 7: Sensitivity analysis on the landscape size (50, 75, 100, 125) ----

#***********************************************************
#XXX do figures & sims



#***********************************************************

# Step 8: Testing with the Eby model  ----

#***********************************************************

## 1) Pseudo-parameters ----


set.seed(123)
range_priors=data.frame(min = c(0, 0),
                        max = c(1, 1))
rownames(range_priors)=c("p","q")

# Latin hypercube sampling on the priors
pseudo_param=as.data.frame(Latinhyper(range_priors, 100000))
write.table(pseudo_param,'../Data/Pseudo_parameters_Eby.csv',sep=";",row.names = F)









## 2) Analysis ----

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
    
    pdf(paste0("../Figures/Eby_model/Cross_validation_n",N_for_cross_validation,"_",
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
          mat_sumstat_pls=pls_1$Yscores[,1:n_comp_pls] # selecting # components
        } else if (n_comp_pls==1){ #otherwise we take the whole components
          mat_sumstat_pls=matrix(pls_1$Yscores[,1:n_comp_pls],ncol=1)
        } else {mat_sumstat_pls=pls_1$Yscores[,1:ncol(pls_1$Yscores)]}
        
        
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
          mat_sumstat_pls2=pls_2$Yscores[,1:n_comp_pls] #pls 2 selecting # components
        } else if (n_comp_pls==1){ #otherwise we take the whole components
          mat_sumstat_pls2=matrix(pls_2$Yscores[,1:n_comp_pls],ncol=1)
        } else {mat_sumstat_pls2=pls_2$Yscores[,1:ncol(pls_2$Yscores)]}
        
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
    
    ggsave(paste0("../Figures/Eby_model/Correlation_parameters_",ifelse(two_step,"twostep","classic"),"_",method_abc,".pdf"),width = 6,height = 5)
    
    
    #Ploting f(x,y) with x=True parameter/summary stat, y=simulated
    ##d-melting the tibble
    
    pdf(paste0("../Figures/Eby_model/x_y_obs_true_param_",ifelse(two_step,"twostep","classic"),"_",method_abc,".pdf"),width = 6,height = 4)
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
    
    
    
    # pdf(paste0("../Figures/Eby_model/x_y_obs_true_summarystat_",ifelse(two_step,"twostep","classic"),"_",method_abc,".pdf"),width = 10,height = 7)
    # 
    # par(mfrow=c(3,3))
    # for (i in colnames(matrix_sumstat)){
    #   d=d_cross_sumstat[,c(which(colnames(d_cross_sumstat)==i),ncol(d_cross_sumstat))]%>%
    #     mutate(., idx=rep(1:(nrow(d_cross_sumstat)/2),each=2))%>%
    #     dcast(., idx ~ Type,value.var=i)%>%
    #     mutate(., Parameter=i)
    #   
    #   plot(d$Obs,d$Sim,xlab="True parameter",ylab="Simulated parameter",col=alpha("blue",.5),main=i)
    #   abline(coef = c(0,1),col="red")
    # }
    # dev.off()
    
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
    
    ggsave(paste0("../Figures/Eby_model/NRMSE_param_",ifelse(two_step,"twostep","classic"),"_",method_abc,".pdf"),p,width = 6,height = 3)
    
    # 
    # p=ggplot(d_NRMSE_sumstat%>%
    #            melt(.))+
    #   geom_jitter(aes(x=variable,y=value,color=variable),
    #               position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0), alpha=.5)+
    #   geom_point(data=as_tibble(t(colMeans(d_NRMSE_sumstat)))%>%
    #                melt(.),
    #              aes(x=variable,y=value),color="black",size=3,shape=18)+
    #   labs(x="Parameter",y="NRMSE",color="")+
    #   geom_hline(yintercept = 1)+
    #   theme_classic()+
    #   theme(legend.position = "none")
    # 
    # ggsave(paste0("../Figures/Eby_model/NRMSE_summarystat_",ifelse(two_step,"twostep","classic"),"_",method_abc,".pdf"),p,width = 8,height = 5)
    # 
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

ggsave('../Figures/Eby_model/NRMSE_post_and_pre_processing.pdf',p,width = 8,height = 5)


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

ggsave('../Figures/Eby_model/Quality_inference_NRMSE_virtual_data.pdf',p,width = 8,height = 5)


## 3) PCA and variables ----

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
             labs(x=paste0("PC 1 (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC 2 (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color=colnames(d_all)[param])+
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
           labs(x=paste0("PC 1 (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC 2 (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="")+
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

# Step 9: Improving inference ----

#***********************************************************

## 1) Number of simulations kept in phase 1, optimization for boxcox parameter and post-processing methods ----
#we play on both the lambda of the boxcox method and the size of the sample for
#stage 1 of the two step preprocessing procedure

d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
matrix_param=d_all[,1:2]
matrix_sumstat=d_all[,3:(ncol(d_all))]

N_for_cross_validation = 100
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)

for (optim_lambda in c(T,F)){

  for (size_step1 in c(1000,2000,3000)){
  
    for (method_abc in c("loclinear","neuralnet")){
      
      for (two_step in c(T,F)){
        
        mat_cor_param=array(0,c(2,2,N_for_cross_validation)) #correlation matrix for parameters
        
        pdf(paste0("../Figures/Optimizing_inferrence/Pre_post/Cross_validation_n",N_for_cross_validation,"_",ifelse(two_step,"twostep","classic"),
                   "_",method_abc,"_optim_lambda_",ifelse(optim_lambda,"yes","no"),"_N1_",size_step1,".pdf"),width = 8,height = 4)
        
        d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
        
        for (n in 1:N_for_cross_validation){
          
          matrix_param=d_all[,1:2]
          matrix_sumstat=d_all[,3:(ncol(d_all))]
          save_sumstat=matrix_sumstat
          
          n_cross=nrow_for_sample[n]
          
          if (two_step){ #Applying the two step procedure used in Siren MEE paper : Don't know whether it make sense in our case. TO discuss Monday
            
            if (optim_lambda==T){
              for (x in 1:ncol(matrix_sumstat)) if (x %in% c(4,6)){
                
                b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
                lambda_x=b$x[which.max(b$y)]
                matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
                
              }else {
                b=boxcox(lm(matrix_sumstat[,x] ~ 1),plotit = F,eps = .05)    
                lambda_x=b$x[which.max(b$y)]
                matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
              }

            } else {
              for (x in 1:ncol(matrix_sumstat)) if (x %in% c(4,6)){
                matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(.5)) -1)/(.5)
              }else {matrix_sumstat[,x] = (matrix_sumstat[,x]^(.5) -1)/(.5)}
            }
            
            #Second we scale
            for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
            
            #and finally, we perform the first PLS
            
            pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                       data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
            
            
            n_comp_pls=selectNcomp(pls_1,method = "onesigma")
            
            
            if (n_comp_pls > 1){
              mat_sumstat_pls=pls_1$Yscores[,1:n_comp_pls] # selecting # components
            } else if (n_comp_pls==1){ #otherwise we take the whole components
              mat_sumstat_pls=matrix(pls_1$Yscores[,1:n_comp_pls],ncol=1)
            } else {mat_sumstat_pls=pls_1$Yscores[,1:ncol(pls_1$Yscores)]}
            
            
            cross_valid=abc(target = mat_sumstat_pls[n_cross,],
                            param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                            tol = size_step1/nrow(matrix_param),method = "rejection") #we keep the size_step1 closest simulations for the first step
            
            #Keeping size_step1 simulations and doing the same steps again: normality, scaling and PLS
            
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
              mat_sumstat_pls2=pls_2$Yscores[,1:n_comp_pls] #pls 2 selecting # components
            } else if (n_comp_pls==1){ #otherwise we take the whole components
              mat_sumstat_pls2=matrix(pls_2$Yscores[,1:n_comp_pls],ncol=1)
            } else {mat_sumstat_pls2=pls_2$Yscores[,1:ncol(pls_2$Yscores)]}
            
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
        
        write.table(d_NRMSE_param,paste0("../Data/Step6_Optimizing_inferrence/Pre_post/RMSE_",ifelse(two_step,"twostep","classic"),
                                         "_",method_abc,"_optim_lambda_",ifelse(optim_lambda,"yes","no"),"_N1_",size_step1,".csv"),sep=";")
        
        colnames(d_NRMSE_sumstat)=colnames(d_cross_sumstat)
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
        
        ggsave(paste0("../Figures/Optimizing_inferrence/Pre_post/NRMSE_param_",ifelse(two_step,"twostep","classic"),
                      "_",method_abc,"_optim_lambda_",ifelse(optim_lambda,"yes","no"),"_N1_",size_step1,".pdf"),p,width = 6,height = 3)
        
      }
    }
  }
}


# Ploting figure size

all_sim=expand.grid(N1=c(1000,2000,3000),lambda=c("yes","no"),Preproc=c("classic","twostep"),postproc=c("loclinear","neuralnet"))

d=tibble()

for (i in 1:nrow(all_sim)){
  d=rbind(d,read.table(paste0("../Data/Step6_Optimizing_inferrence/Pre_post/RMSE_",all_sim$Preproc[i],"_",all_sim$postproc[i],"_optim_lambda_",
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
  theme(strip.text.x = element_text(size=12),axis.text.x = element_text(angle=60,hjust=1),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_color_viridis_d()

ggsave(paste0("../Figures/Optimizing_inferrence/Pre_post/Optimization_inference_all.pdf"),p,width = 7,height = 9)


p=d%>%
  melt(., id.vars=c("N1","optim_lambda","Post","Pre"))%>%
  group_by(.,variable,N1,optim_lambda,Post,Pre)%>%
  summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  ggplot(.)+
  geom_line(aes(x=N1,y=mean_rmse,color=interaction(Pre,Post),group=interaction(Pre,Post,optim_lambda),linetype=optim_lambda),lwd=.8)+
  geom_point(aes(x=N1,y=mean_rmse,fill=interaction(Pre,Post)),color="black",shape=21,size=1)+
  facet_wrap(.~variable,labeller = label_bquote(cols="Parameter = "==.(as.character(variable))))+
  labs(x="# virtual data kept stage 1 pre-processing",color="",linetype=TeX(r"(MLE for $\lambda$)"),y="Mean NRMSE")+
  scale_x_continuous(breaks = c(1000,2000,3000))+
  the_theme+
  guides(color = guide_legend(override.aes = list(size = 2)),fill="none")+
  theme(legend.box = "vertical")+
  scale_color_viridis_d()


ggsave("../Figures/Optimizing_inferrence/Pre_post/Optimization_inference_size_step1_preproc.pdf",p,width = 9,height = 4)



## 2) Structure of the neural-network: number neurons and of neural networks ----



d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)
matrix_param=d_all[,1:2]
matrix_sumstat=d_all[,3:(ncol(d_all))]

N_for_cross_validation = 100
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)

for (size_hidden in seq(5,25,by=5)){
  
  for (rep_network in c(10,20,30)){

    mat_cor_param=array(0,c(2,2,N_for_cross_validation)) #correlation matrix for parameters
    
    pdf(paste0("../Figures/Optimizing_inferrence/Neural_net/Cross_validation_n",N_for_cross_validation,
               "_hidden_",size_hidden,"_Nnet_",rep_network,".pdf"),width = 8,height = 4)
    
    d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
    
    for (n in 1:N_for_cross_validation){
      
      matrix_param=d_all[,1:2]
      matrix_sumstat=d_all[,3:(ncol(d_all))]
      save_sumstat=matrix_sumstat
      
      n_cross=nrow_for_sample[n]
      
      for (x in 1:ncol(matrix_sumstat)) if (x %in% c(4,6)){
        
        b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
        lambda_x=b$x[which.max(b$y)]
        matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
        
      }else {
        b=boxcox(lm(matrix_sumstat[,x] ~ 1),plotit = F,eps = .05)    
        lambda_x=b$x[which.max(b$y)]
        matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
      }
      
      
      #Second we scale
      for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
      
      #and finally, we perform the first PLS
      
      pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo,
                 data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
      
      
      n_comp_pls=selectNcomp(pls_1,method = "onesigma")
      
      
      if (n_comp_pls > 1){
        mat_sumstat_pls=pls_1$Yscores[,1:n_comp_pls] # selecting # components
      } else if (n_comp_pls==1){ #otherwise we take the whole components
        mat_sumstat_pls=matrix(pls_1$Yscores[,1:n_comp_pls],ncol=1)
      } else {mat_sumstat_pls=pls_1$Yscores[,1:ncol(pls_1$Yscores)]}
      
      
      cross_valid=abc(target = mat_sumstat_pls[n_cross,],
                      param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
      
      #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
      
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
        mat_sumstat_pls2=pls_2$Yscores[,1:n_comp_pls] #pls 2 selecting # components
      } else if (n_comp_pls==1){ #otherwise we take the whole components
        mat_sumstat_pls2=matrix(pls_2$Yscores[,1:n_comp_pls],ncol=1)
      } else {mat_sumstat_pls2=pls_2$Yscores[,1:ncol(pls_2$Yscores)]}
      
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
    
    write.table(d_NRMSE_param,paste0("../Data/Step6_Optimizing_inferrence/Neural_net/RMSE_hidden_",
                                     size_hidden,"_Nnet_",rep_network,".csv"),sep=";")
    
    colnames(d_NRMSE_sumstat)=colnames(d_cross_sumstat)
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
    
    ggsave(paste0("../Figures/Optimizing_inferrence/Neural_net/NRMSE_param_hidden_",
                  size_hidden,"_Nnet_",rep_network,".pdf"),p,width = 6,height = 3)
    
    
  }
}


all_sim=expand.grid(N_hidden=seq(5,25,by=5),rep_network=seq(10,30,by=10))

d=tibble()

for (i in 1:nrow(all_sim)){
  d=rbind(d,read.table(paste0("../Data/Step6_Optimizing_inferrence/Neural_net/RMSE_hidden_",
                              all_sim$N_hidden[i],"_Nnet_",all_sim$rep_network[i],".csv"),sep=";")%>%
            add_column(., N_hidden=all_sim$N_hidden[i],N_rep_net=all_sim$rep_network[i]))
}


mean_rmse=d%>%
  melt(., id.vars=c("N_hidden","N_rep_net"))%>%
  group_by(.,variable,N_rep_net,N_hidden)%>%
  summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  rename(., "Parameter"="variable")

p=ggplot(d%>%melt(., id.vars=c("N_hidden","N_rep_net"))%>%
           rename(., "Parameter"="variable"))+
  geom_jitter(aes(x=N_hidden,y=value,color=N_hidden),
              position = position_jitterdodge(jitter.width = 0.3,jitter.height = 0),alpha=.5)+
  geom_point(data=mean_rmse,aes(x=N_hidden,y=mean_rmse),
             color="white",fill="black",shape=24,size=2.5)+
  labs(x="Number hidden neurons",y="NRMSE",color="")+
  facet_grid(Parameter~N_hidden,labeller = label_both)+
  the_theme+
  ylim(0,.5)+
  theme(strip.text.x = element_text(size=12),axis.text.x = element_text(angle=60,hjust=1),legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_color_viridis_d()

ggsave(paste0("../Figures/Optimizing_inferrence/Neural_net/Optimization_NN_all.pdf"),p,width = 7,height = 9)


p=d%>%
  melt(., id.vars=c("N_hidden","N_rep_net"))%>%
  group_by(.,variable,N_hidden,N_rep_net)%>%
  summarise(., .groups = "keep",mean_rmse=mean(value))%>%
  ggplot(.)+
  geom_line(aes(x=N_hidden,y=mean_rmse,color=N_rep_net,group=N_rep_net),lwd=.8)+
  geom_point(aes(x=N_hidden,y=mean_rmse,fill=N_rep_net),color="black",shape=21,size=1)+
  facet_wrap(.~variable,labeller = label_bquote(cols="Parameter = "==.(as.character(variable))))+
  labs(x="Number hidden neurons",color="",y="Mean NRMSE")+
  scale_x_continuous(breaks = seq(5,25,by=5))+
  the_theme+
  guides(color = guide_legend(override.aes = list(size = 2)),fill="none")+
  scale_color_viridis_d()


ggsave("../Figures/Optimizing_inferrence/Neural_net/Optimization_NN_mean.pdf",p,width = 9,height = 4)






# Other: Correlation and EWS ----

#***********************************************************


## Correlation between EWS in a factorial analysis (b, m) ----

d=read.table("../Data/Correlation_EWS_data.csv",sep=",")
colnames(d)=c("r","d", "f", "m","b","c","delta","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
              "Spectral_ratio","PLR","PL_expo")

color_for_m=colorRampPalette(colors = c("#960261","#C5218B","#E266B6","#9F94EF","#2F48DA","#2F48DA"))(30)

pdf("../Figures/Space_param_EWS/mortality_recruitment_EWS.pdf",width = 6,height = 4)
par(mfrow=c(1,1))
for (metric in 1:9){
  for (m_id in seq(1,length(d$m),by=5)){
    
    d_fil=d[which(d$m==unique(d$m)[m_id]),]
    
    if (any (which(is.nan(d_fil[,7+metric])))) d_fil=d_fil[-which(is.nan(d_fil[,7+metric])),]

    
    if (m_id==1 & metric!=3 & metric!=9 & metric !=8){
      plot(1-d_fil[,5],d_fil[,metric+7],col=color_for_m[m_id],xlab="b",ylab=colnames(d)[metric+7],ylim=c(min(d[,metric+7]),max(d[,metric+7])),
           type="p",lwd=2)
    }else if (m_id ==1 & metric==3) {
      plot(1-d_fil[,5],d_fil[,metric+7],col=color_for_m[m_id],xlab="b",ylab=colnames(d)[metric+7],
           type="p",lwd=2,ylim=c(0,7))
    }else if (m_id ==1 & (metric==8 | metric==9)) {
      plot(1-d_fil[,5],d_fil[,metric+7],col=color_for_m[m_id],xlab="b",ylab=colnames(d)[metric+7],
           type="p",lwd=2,ylim=c(0,7))
    } else {
      points(1-d_fil[,5],d_fil[,metric+7],col=color_for_m[m_id],lwd=2)
    }
  }
}
dev.off()









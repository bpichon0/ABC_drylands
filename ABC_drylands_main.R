rm(list=ls())
source("./ABC_drylands_function.R")




# ---------------------- Step 1: Testing with the Eby model  ----

#***********************************************************

## >> 1) Pseudo-parameters ----


set.seed(123)
range_priors=data.frame(min = c(0, 0),
                        max = c(1, 1))
rownames(range_priors)=c("p","q")

# Latin hypercube sampling on the priors
pseudo_param=as.data.frame(Latinhyper(range_priors, 100000))
write.table(pseudo_param,'../Data/Pseudo_parameters_Eby.csv',sep=";",row.names = F)




range_priors=data.frame(min = c(0.01, 0.0001,.25),
                        max = c(.2, .007,.6))

rownames(range_priors)=c("d","m","b")

pseudo_param=as.data.frame(Latinhyper(range_priors, 8000))

write.table(pseudo_param,'./Param_schneider.csv',sep=";",row.names = F)



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
condition_cover=which(d_all$rho_p ==0 | d_all$rho_p==1 | d_all$PLR==1)
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
d_all_kefi=read.table("../Data/All_sim_Kefi.csv",sep=";")%>%
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
           scale_color_manual(values=c("#68B77A","#D4A261"))+
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

# ---------------------- Step 2: Improving inference ----

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






#***********************************************************

# ---------------------- Step 3: Empirical data ----


#***********************************************************
## >> 0) Summary statistics of the data: plus adding regularity of patterns and characteristics of patches ----

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



##         Adding the characteristics of the vegetation patterns: regular or irregular

mat_site=list.files("../Data/Data_Biocom/landscapes")

pdf("../Figures/Empirical_data/Landscapes/r_spectrum_sites.pdf",width = 12,height = 6)
for (i in 1:length(mat_site)){
  landscape=Get_empirical_site(i)
  landscape=as.matrix(landscape>0)
  p1=plot_spectrum(spectral_sews(landscape))+
    ggtitle(i)
  
  print(p1)
  
  p2=plot_spectrum(spectral_sews(landscape[1:round(nrow(landscape)/2),1:round(nrow(landscape)/2)]))+ggtitle("S ",1)+
    theme(axis.title = element_blank(),axis.text = element_text(size=8))
  p3=plot_spectrum(spectral_sews(landscape[(round(nrow(landscape)/2)+1):nrow(landscape),1:round(nrow(landscape)/2)]))+
    ggtitle("S ",2)+
    theme(axis.title = element_blank(),axis.text = element_text(size=8))
  p4=plot_spectrum(spectral_sews(landscape[1:round(nrow(landscape)/2),(round(nrow(landscape)/2)+1):nrow(landscape)]))+
    ggtitle("S ",3)+
    theme(axis.title = element_blank(),axis.text = element_text(size=8))
  p5=plot_spectrum(spectral_sews(landscape[(round(nrow(landscape)/2)+1):nrow(landscape),(round(nrow(landscape)/2)+1):nrow(landscape)]))+
    ggtitle("S ",4)+
    theme(axis.title = element_blank(),axis.text = element_text(size=8))
  Plot_empirical(i)
  print(ggarrange(ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.5,1,.5)),ggarrange(p2,p3,p4,p5,ncol=4),nrow = 2))
}
dev.off()

#script for classifying regular and irregular patterns following Berdugo advices in 2019 paper
Classif_regular_irregular=function(n_begin=1,n_end=345){

  d_subplots=tibble()
  
  for (i in n_begin:n_end){
    d_biocom=read.table("../Data/Step7_Empirical_data/Patterns/biocom_data_patterns.csv",sep=";")
    
    
    
    plot_i=readline(paste0("Is site ",i,"regular (1) or irregular (2) ?, (plot ",d_biocom$File_ID[i],")") ) 
    
    x1=x2=x3=x4=NA
    
    if (plot_i !=1 & plot_i!=2){ #if there is a doubt, we divide into 4 landscapes
      
      x1=readline("Is subplot 1 regular (1) or irregular (2) ?")
      x2=readline("Is subplot 2 regular (1) or irregular (2) ?")  
      x3=readline("Is subplot 3 regular (1) or irregular (2) ?")  
      x4=readline("Is subplot 4 regular (1) or irregular (2) ?")
      print(paste0("Number of regular subsites = ",x1+x2+x3+x4))
      plot_i=readline("Therefore how should it be classified ?")  
      
    }
    d_biocom$Regular[i] = plot_i
    d_subplots=rbind(d_subplots,tibble(Plot=plot_i,subplot1=x1,subplot2=x2,subplot3=x3,subplot4=x4))
    write.table(d_biocom,"../Data/Step7_Empirical_data/Patterns/biocom_data_patterns.csv",sep=";")
    
  }
  
  return(list(d_biocom,d_subplots))
  
}

classif=Classif_regular_irregular(1,345)

#Now we compare the classif with the one in Berdugo et al 2019
d_regular=read.table("../Data/Step7_Empirical_data/Patterns/biocom_data_patterns.csv",sep=";")
d_berdugo=read.table("../Data/Step7_Empirical_data/Patterns/Classif_berdugo.csv",sep=";")

d_regular$Regular=NA
#applying for each site the following rule: if 2 subplot are classified as regular, the site is classified as regular
for (i in seq(1,nrow(d_regular),by=3)){
  d_regular$Regular[i:(i+2)]=ifelse(sum(d_regular$Regular[i:(i+2)])>4,0,1)
}

#adding Berdugo classification
d_regular$Regular_berdugo=sapply(1:nrow(d_regular),function(x){
  return(d_berdugo$ReguConsensus[which(d_berdugo$plotn == gsub("-a","",gsub("-b","",gsub("-c","",d_regular$File_ID[x]))))])
})

table(d_regular$Regular_berdugo,d_regular$Regular)

#We keep the classification in Berdugo et al 2019

write.table(d_regular[,c(1:11,21:22,12:20)],"../Data/Data_Biocom/biocom_data.csv",sep=";")


##           Adding the characteristics of the patches: mean, max, min, sd
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")

d_patches=tibble()
for (i in 1:345){
  land=Get_empirical_site(i)>0
  psd=sort(patchsizes(land>0)) 
  d_patches=rbind(d_patches,tibble(max_psd=max(psd),mean_psd=mean(psd),sd_psd=sd(psd)))
}

d_biocom=cbind(d_biocom,d_patches)
write.table(d_biocom[,c(1:13,23:25,14:22)],"../Data/Data_Biocom/biocom_data.csv",sep=";")






##           Distribution of summary statistics

pdf("../Figures/Empirical_data/Comparizon_models_data/Distrib_sum_stat_data.pdf",width = 7,height = 6)
par(mfrow=c(3,3))
for (sumstat in 1:9){
  hist(d_biocom[,17+sumstat],xlab=colnames(d_biocom)[14+sumstat],main="",col=alpha("blue",.3))
}
dev.off()

#and comparizon with Eby model
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)
p=ggplot(rbind(d_biocom[,17:ncol(d_biocom)]%>%
                 add_column(.,Type="Data")%>%
                 mutate(., Spectral_ratio=log(Spectral_ratio),
                        clustering=log(clustering)),
               d_Eby[,3:ncol(d_Eby)]%>%
                 add_column(., Type="Model")%>%
                 mutate(., Spectral_ratio=log(Spectral_ratio),
                        clustering=log(clustering)))%>%
           melt(., id.vars=c("Type")))+
    geom_density(aes(x=value,fill=Type),alpha=.7)+
    the_theme+
    facet_wrap(.~variable,scales = "free")+
    scale_fill_manual(values=c("#C0CEAA","#5D589C"))

ggsave("../Figures/Empirical_data/Comparizon_models_data/Density_Eby_and_data.pdf",p,width = 7,height = 6)




# Seems to be an effect of the number of pixels. Let's separate the data set into two groups depending on their number of pixels
# And we look at the summary stat of these two group

p1=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")%>%
  ggplot(.)+
  geom_histogram(aes(x=Nbpixels),alpha=.3,fill="blue")+
  geom_vline(xintercept = 80000,col="red")+
  the_theme+
  labs(x="Number of pixels",y="")

p2=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")%>%
  mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering))%>%
  mutate(., Nbpixels=.$Nbpixels>80000)%>%
  melt(., measure.vars=colnames(.)[17:ncol(.)])%>%
  ggplot(.)+
  geom_density(aes(x=value,fill=Nbpixels),alpha=.7)+
  the_theme+
  facet_wrap(.~variable,scales = "free")+
  labs(x="",y="",fill="Large image")+
  scale_fill_manual(values=c("#E4C88D","#AA91CE"))

ggsave("../Figures/Empirical_data/Comparizon_models_data/Understanding_differences/Consequence_image_site_sumstat.pdf",
       ggarrange(ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),p2,nrow=2,labels=letters[1:2],heights = c(.7,1.5)),
       width = 7,height = 8)

#Finally we cross this information with the data resolution of each site that can be accessed in Berdugo et al 2017 NEE data

d_berdugo=as_tibble(readxl::read_xlsx("../Data//BerdugoetalData.xlsx"))
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_biocom$Resolution=sapply(1:nrow(d_biocom),function(x){
  return(d_berdugo$Resolution[which(d_berdugo$plotn==d_biocom$Plot_n[x])])
})


p=d_biocom%>%
  mutate(., Nbpixels=.$Nbpixels>80000)%>%
  ggplot(.)+
  geom_density(aes(x=Resolution,fill=Nbpixels),alpha=.7)+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(x="Image resolution",y="",fill="Large image")+
  scale_fill_manual(values=c("#E4C88D","#AA91CE"))

ggsave("../Figures/Empirical_data/Comparizon_models_data/Understanding_differences/Image_resolution.pdf",p, width = 6,height = 4)


test=d_biocom%>%
  mutate(., Nbpixels=Nbpixels>80000)%>%
  arrange(., Nbpixels)

pdf("./test.pdf",width = 5,height = 5)
for (i in seq(1,nrow(test),by=20)){
  Plot_empirical(i) 
  mtext(paste0(test$Nbpixels[i]))
}
dev.off()



#Finally let's understand the correlations between all summary statistics in the data
p=ggpairs(d_biocom[,c(14:25)]%>%
            mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
                   max_psd=log(max_psd),mean_psd=log(mean_psd),sd_psd=log(sd_psd)))+
  the_theme

ggsave("../Figures/Empirical_data/Pair_correlation_stats.pdf",p,width = 15,height = 15)



## >> 1) Separating each large landscape into 4 of equal size ----

## Let's try to separate each large landscape into 4 of equal size 
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_biocom$Nbpixels=d_biocom$Nbpixels>80000

d_divided=tibble()
for (x in 1:nrow(d_biocom)){
  
  if (d_biocom$Nbpixels[x]){ #if large landscape
    
    for (sub_site in 1:4){
      
      landscape=Get_empirical_site(x)
      
      if (sub_site==1){
        landscape=landscape[1:round(nrow(landscape)/2),1:round(nrow(landscape)/2)] #top left
      } else if (sub_site==2){
        landscape=landscape[(round(nrow(landscape)/2)+1):nrow(landscape),1:round(nrow(landscape)/2)] #bottom left
      } else if (sub_site==3){
        landscape=landscape[1:round(nrow(landscape)/2),(round(nrow(landscape)/2)+1):nrow(landscape)] #top right
      } else {
        landscape=landscape[(round(nrow(landscape)/2)+1):nrow(landscape),(round(nrow(landscape)/2)+1):nrow(landscape)] #bottom right
      }
      
      cover = sum(landscape) / (dim(landscape)[1]**2)
      
      # number of neighbors
      #vegetation clustering
      neighbors_mat = simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
      mean_nb_neigh = mean(neighbors_mat[which(landscape == 1)]) #mean number of plant neighbors
      mean_clustering = mean_nb_neigh / cover
      spatial_ews = generic_sews(landscape>0,4,moranI_coarse_grain = T)$value

      spectral_ratio = as.data.frame(spectral_sews(landscape>0,quiet=T))$value
      
      psd=spatialwarnings::patchdistr_sews(landscape>0)
      max_patchsize=max(psd$psd_obs)
      PLR=spatialwarnings::raw_plrange(landscape>0)
      if (nrow(psd$psd_type)==1){ 
            alpha_exp=NA        
        } else {alpha_exp = psd$psd_type$plexpo[which(psd$psd_type$best==T)]} #i.e., when there is no good fit, return NA


      d_divided=rbind(d_divided,d_biocom[x,c(1:3)]%>%
                        add_column(.,
                                   rho_p=cover,
                                   nb_neigh=mean_nb_neigh,clustering=mean_clustering,
                                   skewness=spatial_ews[2],variance=spatial_ews[1],moran_I=spatial_ews[3],
                                   Spectral_ratio=spectral_ratio,PLR=PLR,PL_expo=alpha_exp,
                                   Id_subsite=sub_site,
                                   Type_landscape="High_res"))
      

      
    }
    
  }else {
    d_divided=rbind(d_divided,d_biocom[x,c(1:3,17:25)]%>%
                      add_column(.,
                                 Id_subsite=1,
                                 Type_landscape="Low_res"))
  }
  print(x)
  
}

# write.table(d_divided,"../Data/Step7_Empirical_data/Dividing_data.csv",sep=";")

d_divided=read.table("../Data/Step7_Empirical_data/Dividing_data.csv",sep=";")
ggplot(d_divided%>%filter(., Type_landscape=="High_res")%>% 
         melt(.,measure.vars=colnames(d_divided)[4:12]))+
  geom_histogram(aes(x=value,fill=as.factor(Id_subsite)),alpha=.7)+
  the_theme+
  labs(x="",y="",fill="Sub-landscape id")+
  facet_wrap(.~variable,scales = "free")

#not so much variation across sublandscapes


ggplot(d_divided%>%
         melt(.,measure.vars=colnames(d_divided)[4:12]))+
  geom_density(aes(x=value,fill=Type_landscape),alpha=.7)+
  the_theme+
  facet_wrap(.~variable,scales = "free")+
  scale_fill_manual(values=c("#C0CEAA","#5D589C"))


d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)
d_tot=rbind(d_Eby[sample(1:nrow(d_Eby),15000,F),3:ncol(d_Eby)]%>%add_column(., Type_landscape="Eby model"),
            d_divided[,c(4:12,14)])%>%
  arrange(., Type_landscape)


p=ggplot(d_tot%>%
           mutate(., Spectral_ratio=log(Spectral_ratio),
                  clustering=log(clustering))%>%
         melt(.,measure.vars=colnames(d_tot)[1:(ncol(d_tot)-1)]))+
  geom_density(aes(x=value,fill=Type_landscape),alpha=.5)+
  the_theme+
  facet_wrap(.~variable,scales = "free")+
  scale_fill_manual(values=c("#7B3636","#FD4848","#C0CEAA"))+
  labs(x="",y="",fill="Type landscape")


ggsave(paste0("../Figures/Empirical_data/Comparizon_models_data/Understanding_differences/",
       "Density_after_spliting_landscapes.pdf"),
       width = 7,height = 5)



d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")

p=ggplot(rbind(d_biocom[,c(11,17:25)]%>%
                 filter(.,Nbpixels>80000)%>%
                 dplyr::select(., -Nbpixels)%>%
                 add_column(., Type_landscape="Before"),
               d_divided[,c(4:12,14)]%>%filter(., Type_landscape=="High_res"))%>%
           mutate(., Type_landscape=recode_factor(Type_landscape,"High_res"="After"))%>%
           mutate(., Spectral_ratio=log(Spectral_ratio),
                  clustering=log(clustering))%>%
           melt(.,measure.vars=colnames(d_tot)[1:(ncol(d_tot)-1)]))+
  geom_density(aes(x=value,fill=Type_landscape),alpha=.5)+
  the_theme+
  facet_wrap(.~variable,scales = "free")+
  scale_fill_manual(values=c("#7B3636","#FD4848","#C0CEAA"))+
  labs(x="",y="",fill="Type landscape")

ggsave(paste0("../Figures/Empirical_data/Comparizon_models_data/",
       "Understanding_differences/Comparizon_before_after_dividing_high_res_landscapes.pdf"),
       width = 7,height = 5)


## >> 2) PCA comparing data and model ---- 
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)
d_Kefi=read.table("../Data/All_sim_Kefi.csv",sep=";")%>%filter(., rho_p !=0)
d_tot=rbind(d_Eby[sample(1:nrow(d_Eby),15000,F),3:ncol(d_Eby)]%>%add_column(., Type="Eby model"),
            d_biocom[,17:ncol(d_biocom)]%>%add_column(., Type=paste0("z",d_biocom$Regular_berdugo)), #to plot the empirical sites above simulations
            d_Kefi[sample(1:nrow(d_Kefi),15000,F),8:ncol(d_Kefi)]%>%add_column(., Type="Kefi model"))%>%
  arrange(., Type)

#First raw PCA

sumstat_name=colnames(d_Eby)[3:ncol(d_Eby)]
res.pca=PCA(d_tot[,which(colnames(d_tot) %in% sumstat_name)], scale.unit = T, ncp = 3,  graph=F)
axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

#getting the centroids
centroids=tibble(Type=d_tot$Type)%>%
  mutate(., Type=recode_factor(Type,"z0"="Empirical sites, irregular","z1"="Empirical sites, regular"))%>%
  add_column(., PC1=res.pca$ind$coord[,1],PC2=res.pca$ind$coord[,2],PC3=res.pca$ind$coord[,3])%>%
  group_by(., Type)%>%
  summarise(., .groups = "keep",C1=mean(PC1),C2=mean(PC2),C3=mean(PC3))%>%
  as.data.frame()
  
  
for (i in 1:3){
  assign(paste0("p",i),
         d_tot%>%
           mutate(., Type=recode_factor(Type,"z0"="Empirical sites, irregular","z1"="Empirical sites, regular"))%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_point(aes(x = PC1, y = PC2, color = Type,fill=Type,size=Type),alpha=.5)+
           scale_size_manual(values=c(1,1,.5,.5))+
           scale_color_manual(values=c("#7B3636","#FD4848",alpha("#C0CEAA",.5),alpha("#ADDDE2",.5)))+
           scale_fill_manual(values=c("#7B3636","#FD4848",alpha("#C0CEAA",.5),alpha("#ADDDE2",.5)))+
           # geom_point(data=centroids,aes(x=centroids[,axes_for_plot$x[i]+1],y=centroids[,axes_for_plot$y[i]+1]),shape=24,fill="black",color="white")+
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


#Second, PCA on Eby model, Kefi, Schneider model only with data

list_model=c("Kefi","Eby","Schneider")

for (model in 1:length(list_model)){
  
  color_model=c(alpha("#C0CEAA",.5),alpha("#ADDDE2",.5),alpha("#3C82E4",.2))
  
  d_sim=read.table(paste0("../Data/All_sim_",list_model[model],".csv"),sep=";")
  d_sim=d_sim[,which(colnames(d_sim) %in% sumstat_name)]
  
  sumstat_name=colnames(d_Eby)[3:ncol(d_Eby)]
  dat=rbind(d_sim%>%
      filter(., rho_p !=0)%>%
      sample_n(., 25000)%>%add_column(., Type=paste0(list_model[model]," model")),
      d_biocom[,17:ncol(d_biocom)]%>%add_column(., Type=paste0("z",d_biocom$Regular_berdugo)))
  
  res.comp=imputePCA(dat[,-ncol(dat)], scale = T, ncp = 3,  graph=F)
  
  if ("completeObs" %in% names(res.comp)){
    res.pca=PCA(res.comp$completeObs,
                scale.unit = T, ncp = 3,  graph=F)
  }else {
    res.pca=PCA(res.comp,
                scale.unit = T, ncp = 3,  graph=F)
  }
  
  axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))
  
  for (i in 1:3){
    assign(paste0("p",i),
           dat%>%
             mutate(., Type=recode_factor(Type,"z0"="Empirical sites, irregular","z1"="Empirical sites, regular"))%>%
             add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = Type,fill=Type,size=Type),alpha=.4)+
             scale_size_manual(values=c(1,1,.5))+
             scale_color_manual(values=c("#7B3636","#FD4848",color_model[model]))+
             scale_fill_manual(values=c("#7B3636","#FD4848",color_model[model]))+
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
  
  ggsave(paste0("../Figures/Empirical_data/Comparizon_models_data/PCA_",list_model[model],"_empirical_data.pdf"),p, width=9,height = 4)
  
  if (model == "Eby"){
    dat=dat%>%
      mutate(., Type=recode_factor(Type,"z0"="Empirical sites, irregular","z1"="Empirical sites, regular"))
    
    for (i in 1:3){
      
      assign(paste0("p",i),
             fviz_pca_biplot(res.pca, geom.ind = "point", 
                             axes=c(axes_for_plot$x[i],axes_for_plot$y[i]), col.ind =dat$rho_p ,col.var="black",
                             label = "var", repel = T)+
               scale_color_viridis_c()+
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
    
    ggsave(paste0("../Figures/Empirical_data/Comparizon_models_data/Understanding_differences/PCA_COVER_Eby_empirical_data.pdf"),p, width=9,height = 4)
  }
}


## PCA with empirical sites colored by the characteristics of the psd distribution

for (metric_psd in c("mean_psd","max_psd",'sd_psd',"Nbpixels")){
  
  d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
  d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)
  d_Kefi=read.table("../Data/All_sim_Kefi.csv",sep=";")%>%filter(., rho_p > 0.1  && rho_p < 0.8)
  d_tot=rbind(d_Eby[sample(1:nrow(d_Eby),15000,F),3:ncol(d_Eby)]%>%add_column(., Type=NA),
              d_Kefi[sample(1:nrow(d_Kefi),15000,F),8:ncol(d_Kefi)]%>%add_column(., Type=NA),
              d_biocom[,17:ncol(d_biocom)]%>%add_column(., Type=d_biocom[,metric_psd]))
  #First raw PCA
  
  sumstat_name=colnames(d_Eby)[3:ncol(d_Eby)]
  res.pca=PCA(d_tot[,which(colnames(d_tot) %in% sumstat_name)], scale.unit = T, ncp = 3,  graph=F)
  axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))
  
  for (i in 1:3){
    assign(paste0("p",i),
           d_tot%>%
             filter(., !is.na(Type))%>%
             add_column(., PC1=res.pca$ind$coord[!is.na(d_tot$Type),axes_for_plot$x[i]],PC2=res.pca$ind$coord[!is.na(d_tot$Type),axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = log(Type),fill=log(Type)),alpha=.5)+
             scale_color_viridis_c(option = "C")+
             scale_fill_viridis_c(option = "C")+
             labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color=metric_psd,fill="")+
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
  ggsave(paste0("../Figures/Empirical_data/Comparizon_models_data/Understanding_differences/PCA_",metric_psd,".pdf"),p, width=9,height = 4)
  
  
}



##    Knock-out of some summary statistics and projection on PCA simulations + data
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



## >> 3) Naive ABC on all empirical sites ----

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
  d_posterior=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
  

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
    
    
    
    d_posterior[index,26]=mean(cross_valid2$adj.values[,1]) #mean posteriors
    d_posterior[index,27]=mean(cross_valid2$adj.values[,2]) #mean posteriors
    index=index+1
    
  }
  dev.off()
  
  
  d_posterior=d_posterior%>%
    rename(., p=V26,q=V27)
  
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




## >> 4) Linking quality estimation with distance in PCA  ----

d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)

d_tot=rbind(d_Eby[,3:ncol(d_Eby)]%>%add_column(., Type="Eby model"),
            d_biocom[,17:ncol(d_biocom)]%>%add_column(., Type="Empirical sites"))

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
            d_biocom[,17:ncol(d_biocom)]%>%add_column(., Type="Empirical sites"))

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



## >> 5) ABC with some summary stat removal ----


d_all=read.table("../Data/All_sim_Eby.csv",sep=";")
condition_cover=which(d_all$rho_p ==0)
d_all=d_all[-condition_cover,]
rownames(d_all)=1:nrow(d_all)

for (preprocessing in c("NOPLS","PLS")[1]){

  for (remove_ID in c("neighB","clustering","spectralR","moranI","variance","neighclust")[6]){
    
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
    } else if (remove_ID == "neighclust"){
      sumstat_kept=colnames(d_all)[3:ncol(d_all)][-c(2:3)]
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
      
      
      d_biocom[site_ID,23]=mean(cross_valid2$adj.values[,1]) #mean posteriors
      d_biocom[site_ID,24]=mean(cross_valid2$adj.values[,2]) #mean posteriors
      
      
    }

    d_biocom=d_biocom%>%
      rename(., p=V26,q=V27)
    
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



#Checking the posterior of each parameters depending on the sumstat that is removed. 

  

d_removal=tibble()
for (x_file in list.files("../Data/Step7_Empirical_data/Removal_sumary_stat/",pattern = "Posteriors")[
  grep(paste0("_NOPLS"),list.files("../Data/Step7_Empirical_data/Removal_sumary_stat/",pattern = "Posteriors"))
]){
  d_removal=rbind(d_removal,read.table(paste0("../Data/Step7_Empirical_data/Removal_sumary_stat/",x_file),sep=";")%>%
                    dplyr::select(., p, q,File_ID)%>%
                    add_column(., Removed =strsplit(split = "_",x_file)[[1]][3]))
}
d_removal=rbind(d_removal,
                read.table(paste0("../Data/Step7_Empirical_data/ABC_all_sites/Posteriors_sites_NOPLS.csv"),sep=";")%>%
                  dplyr::select(., p, q,File_ID)%>%
                  add_column(., Removed="None"))%>%
  melt(.,id.vars=c("Removed","File_ID"))

p=d_removal%>%
  ggplot(.,aes(x=Removed,y=value))+
  geom_line(aes(group=File_ID),lwd=.2,color=alpha("gray",.3))+
  geom_violin(aes(fill=Removed),alpha=.4)+
  facet_wrap(.~variable,scales="free")+
  stat_summary(fun = "median", geom = "point", size = 3,shape=24,fill="black",color="white")+
  the_theme+
  labs(x="",y="NRMSE",color="",fill="")+
  geom_hline(yintercept = 1)+
  guides(fill="none")+
  theme(axis.text.x = element_text(hjust=1,angle=60))
ggsave("../Figures/Empirical_data/ABC/Sensi_removal/Posterior_parameters_NoPLS.pdf",p,width = 8,height = 5)






## >> 6) Filtering the landscapes which are closer to the simulations in PCA ----

#first we visualize whether the metrics we use are relevant in the PCA space: this is good

distances=read.table("../Data/Step7_Empirical_data/Closest_sim/Distance_closest_point.csv",sep=";")$x
dist_empirical=read.table("../Data/Step7_Empirical_data/Closest_sim/Dist_centroid.csv",sep=";")$x
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)%>%
  dplyr::select(., -p, -q)
d_tot=rbind(d_biocom[,17:ncol(d_biocom)],d_Eby)

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


# >> Filtering the sites which perform the best

d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
distances=read.table("../Data/Step7_Empirical_data/Closest_sim/Distance_closest_point.csv",sep=";")$x
dist_empirical=read.table("../Data/Step7_Empirical_data/Closest_sim/Dist_centroid.csv",sep=";")$x

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
  add_column(., Sum_NRMSE = sapply(1:nrow(.),function(x){return(rowSums(.[x,1:9],na.rm = T))}) )%>% 
  add_column(., Best_sites2 = sapply(1:nrow(.),function(x){return(.$Sum_NRMSE[x] <= quantile(.$Sum_NRMSE,.1))}) )%>% 
  melt(., id.vars=c("Best_sites","Site","Sum_NRMSE","Best_sites2"))%>%
  ggplot(.,aes(x=variable,y=value))+
  geom_line(aes(group=Site,color=Best_sites2),lwd=.2)+
  geom_violin(aes(fill=variable),alpha=.4)+
  stat_summary(fun = "median", geom = "point", size = 3,shape=24,fill="black",color="white")+
  the_theme+
  labs(x="",y="NRMSE",color="Best sites ?",fill="")+
  geom_hline(yintercept = 1)+
  scale_color_manual(values=c(alpha("gray",.3),"red"))+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  theme(legend.box = "vertical")

ggsave("../Figures/Empirical_data/ABC/Filtering_data/NRMSE_sumstat_data_best_sites.pdf",p,width = 9,height = 6)







#***********************************************************

# ---------------------- Step 4: Bringing data and models closer ----

#***********************************************************
## >> 1) PCA of models and data for space cover of different models ----

all_dat=list.files("../Data/Step8_Solutions_gap/Data_Eby_feedbacks/")
all_dat=all_dat[-grep(".csv",all_dat)]

d_merged_sim=tibble()
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")

for (type_feedback in all_dat){
  
  list_simu=list.files(paste0('../Data/Step8_Solutions_gap/Data_Eby_feedbacks/',type_feedback),pattern = ".csv")
  
  d_all=tibble()
  for (file_simu in list_simu){
    
    d=read.table(paste0("../Data/Step8_Solutions_gap/Data_Eby_feedbacks/",type_feedback,"/",file_simu),sep=",")
    colnames(d)= c("p","q",
                   "rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                   "Spectral_ratio","PLR","PL_expo")
    
    d_all=rbind(d_all,d)
  }
  
  d_all=filter(d_all,rho_p!=0 & PL_expo >= 0 | is.na(PL_expo))
  rownames(d_all)=1:nrow(d_all)
  write.table(d_all,paste0("../Data/Step8_Solutions_gap/Data_Eby_feedbacks/",type_feedback,".csv"),sep=";")
  
  
  #classic one colored by simulation type
  
  d_tot=rbind(d_biocom[,17:ncol(d_biocom)]%>%add_column(., Type=paste0("z",d_biocom$Regular_berdugo),NBpixels=d_biocom$Nbpixels), #to plot the empirical sites above simulations
              d_all[sample(1:nrow(d_all),20000,F),3:ncol(d_all)]%>%add_column(., Type=type_feedback,NBpixels=NA))%>%
    arrange(., Type)
  
  #First raw PCA
  
  sumstat_name=colnames(d_all)[3:ncol(d_all)]
  res.comp=imputePCA(d_tot[,which(colnames(d_tot) %in% sumstat_name)],ncp=3,scale = T) 
  res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
  axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))
  
  
  for (i in 1:3){
    assign(paste0("p",i),
           d_tot%>%
             mutate(., Type=recode_factor(Type,"z0"="Empirical sites, irregular","z1"="Empirical sites, regular"))%>%
             add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = Type,fill=Type,size=Type),alpha=.5)+
             scale_size_manual(values=c(1,1,.5,.5,.5))+
             scale_color_manual(values=c("#7B3636","#FD4848",alpha("#C0CEAA",.5)))+
             scale_fill_manual(values=c("#7B3636","#FD4848",alpha("#C0CEAA",.5)))+
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
  ggsave(paste0("../Figures/Eby_model/With_feedbacks/PCA_coverage_",type_feedback,".pdf"),p, width=9,height = 4)
  
  
  
  #colored by pixel size  simulations
  
  sumstat_name=colnames(d_all)[3:ncol(d_all)]
  res.comp=imputePCA(d_tot[,which(colnames(d_tot) %in% sumstat_name)],ncp=3,scale = T) 
  res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
  axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))
  
  
  for (i in 1:3){
    assign(paste0("p",i),
           d_tot%>%
             mutate(., Type=recode_factor(Type,"z0"="Empirical sites, irregular","z1"="Empirical sites, regular"))%>%
             add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = log(NBpixels),fill=log(NBpixels),size=(Type)),alpha=.5)+
             scale_size_manual(values=c(1,1,.5,.5,.5))+
             scale_color_viridis_c()+
             scale_fill_viridis_c()+
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
  ggsave(paste0("../Figures/Eby_model/With_feedbacks/PCA_coverage_PIXELS_",type_feedback,".pdf"),p, width=9,height = 4)
  
  
  
  
  # ONly small pixel images
  
  sumstat_name=colnames(d_all)[3:ncol(d_all)]
  d_tot=filter(d_tot,NBpixels<80000 | is.na(NBpixels))
  res.comp=imputePCA(d_tot[,which(colnames(d_tot) %in% sumstat_name)],ncp=3,scale = T) 
  res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
  axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))
  
  
  for (i in 1:3){
    assign(paste0("p",i),
           d_tot%>%
             mutate(., Type=recode_factor(Type,"z0"="Empirical sites, irregular","z1"="Empirical sites, regular"))%>%
             add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = log(NBpixels),fill=log(NBpixels),size=(Type)),alpha=.5)+
             scale_size_manual(values=c(1,1,.5,.5,.5))+
             scale_color_viridis_c()+
             scale_fill_viridis_c()+
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
  
  ggsave(paste0("../Figures/Eby_model/With_feedbacks/PCA_coverage_only_small_pixels_",type_feedback,".pdf"),p, width=9,height = 4)
  
  d_merged_sim=rbind(d_merged_sim,d_all%>%add_column(., Type=type_feedback))
  
}

p=ggplot(d_merged_sim%>%
           mutate(., clustering=log(clustering),Spectral_ratio=log(Spectral_ratio))%>%
           melt(., id.vars=c("Type","p","q")))+
  geom_density(aes(x=value,fill=Type),alpha=.5)+
  the_theme+
  facet_wrap(.~variable,scales = "free")

ggsave("../Figures/Eby_model/With_feedbacks/Merged_densities.pdf",p, width=7,height = 7)





## >> 2) ABC to compare the models ----


d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")

#we only save the sumstats of the closest simulations. for posteriors we need to calibrate the inference
d_cross_sumstat=d_NRMSE_sumstat=tibble()

all_dat=list.files("../Data/Step8_Solutions_gap/Data_Eby_feedbacks/")


for (type_feedback in all_dat){
  
  list_simu=list.files(paste0('../Data/Step8_Solutions_gap/Data_Eby_feedbacks/',type_feedback),pattern = ".csv")
  
  d_all=tibble()
  for (file_simu in list_simu){
    
    d=read.table(paste0("../Data/Step8_Solutions_gap/Data_Eby_feedbacks/",type_feedback,"/",file_simu),sep=",")
    colnames(d)= c("p","q",
                   "rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                   "Spectral_ratio","PLR","PL_expo")
    d_all=rbind(d_all,d)
  }
  d_all=filter(d_all,rho_p!=0,!is.na(PL_expo),PL_expo>0)
  rownames(d_all)=1:nrow(d_all)
  
  d_cross_param=d_NRMSE_param=d_cross_sumstat=d_NRMSE_sumstat=tibble()
  
  d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")


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
    
    #and finally, we perform the first PLS (excluding empirical data)
    
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
                     tol = 75/nrow(mat_sumstat_pls2),method = "rejection") # we only perform simple rejection here as we are not yet interested in the posterior
    
    cross_valid2$ss=as.data.frame(cross_valid2$ss)
    cross_valid2$ss=d_all[as.numeric(rownames(cross_valid2$ss)),which(colnames(d_all) %in% names(observed_sumstat))] #we keep information with the true values
    
    matrix_sumstat=save_sumstat
    
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
    
    
  }

  write.table(d_NRMSE_sumstat,paste0("../Data/Step8_Solutions_gap/Data_Eby_feedbacks//NRMSE_sumstat_",type_feedback,".csv"),sep=";")
  write.table(d_cross_sumstat,paste0("../Data/Step8_Solutions_gap/Data_Eby_feedbacks/ABC_empirical/Cross_valid_sumstat_",type_feedback,".csv"),sep=";")
}

d_NRMSE_sumstat=tibble()
for (i in list.files("../Data/Step8_Solutions_gap/Data_Eby_feedbacks/ABC_empirical/",pattern = "NRMSE")){
  d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,
                        read.table(paste0("../Data/Step8_Solutions_gap/Data_Eby_feedbacks/ABC_empirical/",i),sep=";")%>%
                          add_column(.,Type_feedback=gsub(".csv","",gsub("NRMSE_sumstat_4NN_","",i)))
  )
}

p=d_NRMSE_sumstat%>%
  melt(., id.vars=c("Type_feedback","Site"))%>%
  ggplot(.,aes(x=Type_feedback,y=value))+
  geom_line(aes(group=Site),lwd=.2,color=alpha("gray",.3))+
  geom_violin(aes(fill=variable),alpha=.4)+
  facet_wrap(.~variable,scales = "free")+
  stat_summary(fun = "median", geom = "point", size = 3,shape=24,fill="black",color="white")+
  the_theme+
  labs(x="",y="NRMSE",color="",fill="")+
  geom_hline(yintercept = 1)+
  guides(fill="none") +
  theme(axis.text.x = element_text(angle=60,hjust=1))

ggsave("../Figures/Eby_model/With_feedbacks/Comparing_summary_stat_ABC_feedback.pdf",p,width = 10,height = 9)








## >> 3) Pooling pixels in empirical landscapes ----

d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")

pdf("../Figures/Empirical_data/Pooling/Pooling_landscapes.pdf",width = 9,height = 3)
par(mfrow=c(1,2),mar=c(1,1,1,1))
example_landcapes=which(d_biocom$Nbpixels>80000)
d_pooled=d_double_pool=tibble()
for (id in example_landcapes){
  image(Get_empirical_site(id),col=c("white","black"),axes=F)
  
  simple_pool=pooling(Get_empirical_site(id),2)>.5
  
  image(simple_pool,col=c("white","black"),axes=F)
  d_pooled=rbind(d_pooled,Get_sumstat(simple_pool))
  
  douple_pool=pooling(pooling(Get_empirical_site(id),2)>.3,2)>.3
  image(douple_pool,col=c("white","black"),axes=F)
  d_double_pool=rbind(d_double_pool,Get_sumstat(douple_pool))
}
dev.off()


#for the example

pdf("../Figures/Empirical_data/Pooling/Example_simple_pooling_landscapes.pdf",width = 6,height = 6)
par(mfrow=c(2,2),mar=c(1,1,1,1))
example_landcapes=which(d_biocom$Nbpixels>80000)
for (id in example_landcapes[c(10,50)]){
  image(Get_empirical_site(id),col=c("white","black"),axes=F)
  
  simple_pool=pooling(Get_empirical_site(id),2)>.5
  image(simple_pool,col=c("white","black"),axes=F)
  
}
dev.off()


pdf("../Figures/Empirical_data/Pooling/Example_double_pooling_landscapes.pdf",width = 9,height = 3)
par(mfrow=c(1,3),mar=c(1,1,1,1))
example_landcapes=which(d_biocom$Nbpixels>80000)
for (id in example_landcapes[50]){
  image(Get_empirical_site(id),col=c("white","black"),axes=F)
  
  simple_pool=pooling(Get_empirical_site(id),2)>.5
  image(simple_pool,col=c("white","black"),axes=F)

  douple_pool=pooling(pooling(Get_empirical_site(id),2)>.3,2)>.3
  image(douple_pool,col=c("white","black"),axes=F)
}
dev.off()


#comparing sumstat

d_plot=rbind(d_pooled%>%
                 add_column(., Site=rep(example_landcapes,1)),
               d_biocom[example_landcapes,17:25]%>%
                 add_column(., Site=example_landcapes))%>%
  add_column(., Pooled=rep(c("yes, 2","no"),each=length(example_landcapes)))


p=ggplot(d_plot%>%
           mutate(., clustering=log(clustering),Spectral_ratio=log(Spectral_ratio))%>%
           melt(., id.vars=c("Pooled","Site")))+
  geom_line(aes(x=Pooled,y=value,group=Site),color="gray")+
  facet_wrap(.~variable,scales="free")+
  the_theme+
  labs(x="Pooled landscapes ?",y="")
ggsave("../Figures/Empirical_data/Pooling/Change_in_sumstat_pooling.pdf",p,width = 7,height=6)

write.table(d_pooled,"../Data/Step8_Solutions_gap/Pooling_pixels/Pooled_sites.csv",sep=";")
write.table(d_double_pool,"../Data/Step8_Solutions_gap/Pooling_pixels/Pooled_sites_double.csv",sep=";")


## Density of sumstat when compared to empirical data and low res landscapes
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_pooled=read.table("../Data/Step8_Solutions_gap/Pooling_pixels/Pooled_sites.csv",sep=";")
d_biocom[which(d_biocom$Nbpixels>80000),17:25]=d_pooled[,1:9]


p=d_biocom%>%
  mutate(., Nbpixels=.$Nbpixels>80000)%>%
  mutate(.,Spectral_ratio=log(Spectral_ratio),clustering=log(clustering))%>%
  melt(., measure.vars=colnames(.)[17:ncol(.)])%>%
  ggplot(.)+
  geom_density(aes(x=value,fill=Nbpixels),alpha=.7)+
  the_theme+
  facet_wrap(.~variable,scales = "free")+
  labs(x="",y="",fill="Large image")+
  scale_fill_manual(values=c("#E4C88D","#AA91CE"))

ggsave("../Figures/Empirical_data/Pooling/Sumstat_pooling_landscapes.pdf",p,width = 7,height=6)



## Seeing how these sites are changed in the PCA
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%filter(., rho_p !=0)
d_tot=rbind(d_Eby[,3:ncol(d_Eby)]%>%add_column(., Type="Eby model"),
            d_biocom[,17:ncol(d_biocom)]%>%add_column(., Type=paste0("z",d_biocom$Regular_berdugo)))%>%
  arrange(., Type)

#First raw PCA

sumstat_name=colnames(d_Eby)[3:ncol(d_Eby)]
res.pca=PCA(d_tot[,which(colnames(d_tot) %in% sumstat_name)], ncp = 3,  graph=F)
axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))


d_pooled=read.table("../Data/Step8_Solutions_gap/Pooling_pixels/Pooled_sites.csv",sep=";")
d_pooled$PL_expo[is.na(d_pooled$PL_expo)]=0
PCA_push=predict.PCA(res.pca,d_pooled[,1:9])#predicting the projection
d_all=rbind(res.pca$ind$coord[(nrow(d_Eby)+1):nrow(res.pca$ind$coord),][which(d_biocom$Nbpixels>80000),],PCA_push$coord)%>%
  as.data.frame(.)%>%
  add_column(., ID=rep(1:length(which(d_biocom$Nbpixels>80000)),2))


for (i in 1:3){
  assign(paste0("p",i),
         tibble(PC1=d_all[,axes_for_plot$x[i]],
                PC2=d_all[,axes_for_plot$y[i]],
                ID=d_all$ID,group=rep(c("Before pooling","After pooling"),each=length(which(d_biocom$Nbpixels>80000))))%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_line(aes(x = PC1, y = PC2, group=ID),color="gray",alpha=.7,lwd=.5)+
           geom_point(aes(x = PC1, y = PC2, color = group,fill=group),size=2,alpha=.5)+
           scale_color_manual(values=c("#B32DA6","#64C568"))+
           scale_fill_manual(values=c("#B32DA6","#64C568"))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="",fill="")+
           ggtitle("")+guides(shape="none")+
           theme_classic()+theme(legend.position = "bottom")+
           guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")  
  )
}


p=ggarrange(ggarrange(p1+theme(legend.position = "none",plot.title=element_blank()),
                      p2+theme(legend.position = "none",plot.title=element_blank()),
                      p3+theme(legend.position = "none",plot.title=element_blank()),
                      ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.1))


ggsave("../Figures/Empirical_data/Pooling/Change_in_PCA_space_simple_pooling.pdf",p,width = 9,height=4)



d_pooled=read.table("../Data/Step8_Solutions_gap/Pooling_pixels/Pooled_sites_double.csv",sep=";")
d_pooled$PL_expo[is.na(d_pooled$PL_expo)]=0
PCA_push=predict.PCA(res.pca,d_pooled[,1:9])#predicting the projection
d_all=rbind(res.pca$ind$coord[(nrow(d_Eby)+1):nrow(res.pca$ind$coord),][which(d_biocom$Nbpixels>80000),],PCA_push$coord)%>%
  as.data.frame(.)%>%
  add_column(., ID=rep(1:length(which(d_biocom$Nbpixels>80000)),2))


for (i in 1:3){
  assign(paste0("p",i),
         tibble(PC1=d_all[,axes_for_plot$x[i]],
                PC2=d_all[,axes_for_plot$y[i]],
                ID=d_all$ID,group=rep(c("Before pooling","After pooling"),each=length(which(d_biocom$Nbpixels>80000))))%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_line(aes(x = PC1, y = PC2, group=ID),color="gray",alpha=.7,lwd=.5)+
           geom_point(aes(x = PC1, y = PC2, color = group,fill=group),size=2,alpha=.5)+
           scale_color_manual(values=c("#B32DA6","#64C568"))+
           scale_fill_manual(values=c("#B32DA6","#64C568"))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="",fill="")+
           ggtitle("")+guides(shape="none")+
           theme_classic()+theme(legend.position = "bottom")+
           guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")  
  )
}


p=ggarrange(ggarrange(p1+theme(legend.position = "none",plot.title=element_blank()),
                      p2+theme(legend.position = "none",plot.title=element_blank()),
                      p3+theme(legend.position = "none",plot.title=element_blank()),
                      ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.1))


ggsave("../Figures/Empirical_data/Pooling/Change_in_PCA_space_double_pooling.pdf",p,width = 9,height=4)






# ---------------------- Step 5: Increasing/decreasing spatial resolution artificially ----


## 1) >> Doing it on empirical data ----

d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")


d_pooled=tibble()
for (id in 1:345){
  landscape=Get_empirical_site(id)
  for (pooling_coeff in c(.25,.33,.5,2,3,4)){
    d_pooled=rbind(d_pooled,Get_sumstat(pooling(landscape,pooling_coeff))%>%
                     add_column(., Pooling=pooling_coeff))
    print(id)
  }
}


write.table(d_pooled,"../Data/Step9_Spatial_resolution/d_pooled_empirical_data.csv",sep=";")




# Ploting the change in summary statistics values against the resolution


d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_pooled=read.table("../Data/Step9_Spatial_resolution/d_pooled_empirical_data.csv",sep=";")

dat=rbind(d_pooled%>%add_column(., Site=rep(1:345,each=length(unique(d_pooled$Pooling)))),
          d_biocom[,17:ncol(d_biocom)]%>%add_column(., Pooling=1, Site=1:345))%>%
  mutate(., clustering=log(clustering),Spectral_ratio=log(Spectral_ratio))%>%
  group_by(., Site)%>%arrange(.,Site,Pooling)

for (i in 1:345){
  dat$PL_expo[(7*(i-1)+1):(7*(i-1)+3)]=dat$PL_expo[(7*(i-1)+4)]
}

p=ggplot(dat%>%
         filter(., Pooling<3)%>%
         melt(., id.vars=c("Pooling","Site"))%>%
         mutate(., Pooling=as.numeric(as.factor(Pooling))),
       aes(x=Pooling,y=value))+
  geom_line(aes(group=Site),color="gray70",alpha=.2)+
  stat_summary(fun = "mean", geom = "line", lwd = 1,fill="red",color="red")+
  the_theme+
  facet_wrap(.~variable,scales = "free")+
  labs(x="Change in resolution",y="Mean value")+
  scale_x_reverse(labels = c("x4","x3","x2","1","/2","/3","/4"),
                     breaks = 1:length(unique(dat$Pooling)))

ggsave("../Figures/Spatial_resolution/Data/Change_metrics_spatial_resolution_data.pdf",p,width = 7,height = 6)





#doing a PCA colored by spatial resolution

sumstat_name=colnames(dat)[1:9]
res.comp=imputePCA(dat[,which(colnames(dat) %in% sumstat_name)],ncp=3,scale = T) 
res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))


for (i in 1:3){
  assign(paste0("p",i),
         fviz_pca_biplot(res.pca, geom.ind = "point", 
                         axes=c(axes_for_plot$x[i],axes_for_plot$y[i]), col.ind = log(dat$Pooling),col.var="black",
                         label = "var", repel = T,alpha=.3)+
           scale_color_viridis_c(option = "C",breaks=c(-1,1),labels=c("High res","Low res"))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="")+
           ggtitle("")+
           theme_classic()+theme(legend.position = "bottom",legend.text = element_text(size=12))
  )
}

p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                      p2+theme(legend.position = "none"),
                      p3+theme(legend.position = "none"),
                      ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.2))

ggsave("../Figures/Spatial_resolution/Data/PCA_spatial_resolution_data.pdf",p,width = 11,height = 5)








# How does it impacts the EWS trends
condition_res=which(d_biocom$Nbpixels>80000)

d_pooled=filter(d_pooled, Pooling<3) #otherwise too much, we don't keep the structure

dat=rbind(d_pooled%>%add_column(., ID=rep(1:345,each=length(unique(d_pooled$Pooling))))%>%
            filter(., ID %!in% condition_res), #we exclude data that are with initial high cover
  
          d_biocom[-condition_res,17:ncol(d_biocom)]%>%
            add_column(., Pooling=1, ID=c(1:345)[-condition_res]))%>%
  
  mutate(., clustering=log(clustering),Spectral_ratio=log(Spectral_ratio),
         Aridity=c(rep(d_biocom$Aridity[-condition_res],each=length(unique(d_pooled$Pooling))),d_biocom$Aridity[-condition_res]))


d_trends=tibble()
dat_melt=melt(dat, id.vars=c("Pooling","Aridity","ID"))

for (pool in unique(dat_melt$Pooling)){
  for (sumstat in unique(dat_melt$variable)){
    
    #fitting linear mixed effect model with site as random factor
    
    mod=lmer(value~Aridity+(1|Site),data = filter(dat_melt,variable==sumstat,Pooling==pool)%>%
               add_column(., Site=rep(1:66,each=3))%>%
               filter(., !is.na(value))%>%
               mutate(., value=(value-mean(value))/sd(value))) #standardizing
    
    # boot_mod=bootMer(mod,FUN=function(x) predict(x, newdat, re.form=NA),nsim = 1000)
    
    
    d_trends=rbind(d_trends,tibble(Sumstat=sumstat,Pooling=pool,Trend=mod@beta[2]))
    
                                   # Low_quant=apply(boot_mod$t, 2, quantile, 0.025),
                                   # High_quant=apply(boot_mod$t, 2, quantile, 0.975))) #getting the trend
    
  }
}



p=d_trends%>%filter(., Sumstat %!in% c("rho_p","nb_neigh","clustering"))%>%
  mutate(.,Sign_trend=.$Trend>0,
         Value_pooling=sapply(1:nrow(.),function(x){
           if (.$Pooling[x]>1){return("Decreasing")
           }else if (.$Pooling[x]==1){return("Reference")
           }else {
             return("Increasing")
           }
         }))%>%
  mutate(., Pooling_intensity=sapply(1:nrow(.),function(x){
    if (.$Pooling[x] %in% c(.25,4)){return("High")
    }else if (.$Pooling[x] %in% c(.5,2)){return("Low")
    }else {
      return("Medium")
    }
  }))%>%
  ggplot(.)+
  geom_point(aes(x=Trend,y=Sumstat,color=Sign_trend,size=as.character(Pooling_intensity),shape=Value_pooling))+
  geom_vline(xintercept = 0,linetype=9)+
  the_theme+
  scale_color_manual(values=c("red","blue"))+
  guides(color="none")+
  labs(shape="Pooling",x="Trend along aridity gradient",y="",shape="Change in resolution",size="Scale of change")+
  scale_y_discrete(labels=c("Moran I","PL exponent","PLR","Skewness","Spectral ratio","Variance"))+
  scale_shape_manual(values=c(11,8,7))+
  scale_size_manual(values=c(5,2,3.5))+
  theme(legend.box = "vertical")

ggsave("../Figures/Spatial_resolution/Model/Change_of_trends_spatial_resolution_data.pdf",p,width = 7,height = 6)


p=dat%>%
  melt(., id.vars=c("Aridity","ID","Pooling"))%>%
  mutate(., Pooling=as.character(round(Pooling,2)))%>%
  mutate(.,Pooling=recode_factor(Pooling,"0.25"="x4","0.33"="x3","0.5"="x2",'1'="Original","2"="/2"))%>%
  mutate(.,variable=recode_factor(variable,"Spectral_ratio"="SDR","moran_I"="Moran I"))%>%
  filter(., variable %!in% c("rho_p","nb_neigh","clustering"))%>%
  ggplot(.)+
  geom_point(aes(x=Aridity,y=value,color=variable),alpha=.5,size=1)+
  the_theme+
  facet_grid(variable~Pooling,scales="free",labeller = label_bquote(cols="Res "==.(as.character(Pooling))))+
  labs(x="Aridity",y="",color="")+
  guides(color = guide_legend(override.aes = list(size = 3)))
         


ggsave("../Figures/Spatial_resolution/Model/Explicit_change_sumstat_aridity_data.pdf",p,width = 7,height = 6)





## 2) >> Doing it on simulations ----

#For that we first select 25000 simulations per model

d_Eby=read.table("../Data/All_sim_Eby.csv",sep=";")%>%
  filter(., rho_p > 0.01 & !is.na(PLR) & !is.na(PL_expo))%>%
  sample_n(., 25000)
d_Eby_feedback=read.table("../Data/All_sim_Eby_feedback.csv",sep=";")%>%
  filter(., rho_p > 0.01 & !is.na(PLR) & !is.na(PL_expo))%>%
  sample_n(., 25000)
d_Kefi=read.table("../Data/All_sim_Kefi.csv",sep=";")%>%
  filter(., rho_p > 0.01 & !is.na(PLR) & !is.na(PL_expo))%>%
  sample_n(., 25000)
d_Schneider=read.table("../Data/All_sim_Schneider.csv",sep=";")%>%
  filter(., rho_p > 0.01 & !is.na(PLR) & !is.na(PL_expo))%>%
  sample_n(., 25000)


write.table(d_Eby[,-c(3:ncol(d_Eby))],
            '../Data/Step9_Spatial_resolution/Simu/Parameter_Eby.csv',sep=";",
            row.names = F)
write.table(d_Eby_feedback[,-c(3:ncol(d_Eby_feedback))],
            '../Data/Step9_Spatial_resolution/Simu/Parameter_Eby_feedback.csv',sep=";",
            row.names = F)
write.table(d_Kefi[,-c(8:ncol(d_Kefi))],
            '../Data/Step9_Spatial_resolution/Simu/Parameter_Kefi.csv',sep=";",
            row.names = F)
write.table(d_Schneider[,-c(9:ncol(d_Schneider))],
            '../Data/Step9_Spatial_resolution/Simu/Parameter_Schneider.csv',sep=";",
            row.names = F)




# Once simulations are made in Julia, we extract and analyse the results

list_folders=list.dirs("../Data/Step9_Spatial_resolution/Simu",recursive = F)[-5]
d_simu=tibble()

for (k in list_folders){
  
  if (tail(strsplit(k,"/")[[1]],1) %in% c("Eby","Eby_feedback")){
    n_param=2
  }else if (tail(strsplit(k,"/")[[1]],1)=="Kefi"){
    n_param=7
  }else {
    n_param=8
  }
  
  list_sim=list.files(k,".csv")
  
  d=tibble()
  for (sim in list_sim){
    d=rbind(d,read.table(paste0(k,"/",sim),sep=",")[,(n_param+1):(n_param+9)])
  }
  
  colnames(d)=c("rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                "Spectral_ratio","PLR","PL_expo")

  d=add_column(d,Pooling=rep(c("1/4","1/3",'1/2',"1"),nrow(d)/4))%>%
    add_column(., Model=tail(strsplit(k,"/")[[1]],1))%>%
    filter(., rho_p != 0)
  

  d_simu=rbind(d_simu,d)
}


write.table(d_simu,"../Data/Step9_Spatial_resolution/All_sims_models.csv",sep=";")


#adding Schneider aggregated simulations

d_simu2=list.files("../Data/Step9_Spatial_resolution/Simu/Schneider_Aggregated/",".csv")
n_param=3
d=tibble()
for (sim in d_simu2){
  d=rbind(d,read.table(paste0("../Data/Step9_Spatial_resolution/Simu/Schneider_Aggregated/",sim),sep=",")[,(n_param+1):(n_param+9)])
}

colnames(d)=c("rho_p","nb_neigh","clustering","skewness","variance","moran_I",
              "Spectral_ratio","PLR","PL_expo")

d=add_column(d,Pooling=rep(c("1"),nrow(d)))%>%
  add_column(., Model="Schneider")%>%
  filter(., rho_p > 0.03)

d_simu=rbind(d_simu,d)

write.table(d_simu,"../Data/Step9_Spatial_resolution/All_sims_models_with_Schn_aggre.csv",sep=";")


#Doing PCA and coloring by pooling intensity for each model

d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models.csv",sep=";"),
             Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
             PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
  filter(., rho_p>0.03)%>%
  arrange(.,Model)

for (i in 1:(nrow(d_sim)/4)){
  d_sim$PL_expo[(4*(i-1)+1):(4*(i-1)+3)]=d_sim$PL_expo[(4*(i-1)+4)]
}

#d_sim$PL_expo[is.na(d_sim$PL_expo)]=0
d_biocom$PL_expo[d_biocom$PL_expo==0]=NA
index=1

for (model in unique(d_sim$Model)){
  
  
  #simulations alone
  sumstat_name=colnames(d_sim)[1:9]
  res.comp=imputePCA(d_sim[which(d_sim$Model==model),which(colnames(d_sim) %in% sumstat_name)],ncp=3,scale = T) 

  if ("completeObs" %in% names(res.comp)){
    res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
  }else {
    res.pca=PCA(res.comp, ncp = 3,  graph=F)
  }
  
  axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))
  
  
  for (i in 1:3){
    assign(paste0("p",i),
           d_sim%>%
             filter(., Model==model)%>%
             mutate(., Pooling=recode_factor(Pooling,"1/4"='x4',"1/3"="x3","1/2"="x2","1"="No change"))%>%
             add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = Pooling,fill=Pooling),alpha=.3)+
             scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))+
             scale_fill_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))+
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
  
  ggsave(paste0("../Figures/Spatial_resolution/Model/PCA_spatial_resolution_",model,"_model.pdf"),p,width = 11,height = 5)
  
  #Simulations and empirical data
  
  d_sim_data=rbind(d_sim%>%filter(., Model==model),
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
             geom_point(aes(x = PC1, y = PC2, color = Pooling,fill=Pooling),alpha=.5)+
             scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","black"))+
             scale_fill_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","black"))+
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
  
  ggsave(paste0("../Figures/Spatial_resolution/Model/PCA_spatial_resolution_",model,"_model_and_data.pdf"),p,width = 11,height = 5)
  

  d_sim_data=rbind(d_sim%>%filter(., Model==model)%>%add_column(.,Vege=NA),
                   d_biocom[which(d_biocom$Nbpixels<80000),17:ncol(d_biocom)]%>% #keeping only the data with relatively low quality
                     add_column(., Pooling="Data",Model=NA,Vege=.$rho_p)%>%
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
             geom_point(aes(x = PC1, y = PC2, color = Vege,fill=Vege))+
             scale_color_viridis_c(alpha = 1,na.value = alpha("gray",.3))+
             scale_fill_viridis_c(alpha = 1,na.value = alpha("gray",.3))+
             labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="Vegetation cover",fill="")+
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
  
  ggsave(paste0("../Figures/Spatial_resolution/Model/PCA_spatial_resolution_",model,"_model_and_data_COVER.pdf"),p,width = 11,height = 5)
  
  
  
  # Finally, we compare the density of simulations versus data
  
  
  d_sim_data=rbind(d_sim%>%filter(., Model==model,Pooling != "1"),
                   d_biocom[which(d_biocom$Nbpixels<80000),17:ncol(d_biocom)]%>% #keeping only the data with relatively low quality
                     add_column(., Pooling="Data",Model="zEmp")%>%
                     mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering)))%>%
    arrange(., Model)
  
  
  assign(paste0("p",index),ggplot(d_sim_data%>%melt(., id.vars=c("Pooling","Model")))+
    geom_density(aes(x=value,fill=Model),alpha=.3)+
    facet_wrap(.~variable,scales = "free")+
    theme_classic()+theme(legend.position = "bottom")+
    scale_fill_manual(values=c("blue","forestgreen"),labels=c("Sim","Emp")))
  
  
  index=index+1
}

p=ggarrange(p1+ggtitle("Eby"),
            p2+ggtitle("Eby feedback"),
            p3+ggtitle("Kefi"),
            p4+ggtitle("Schneider"),nrow=2,ncol=2,common.legend = T,legend="bottom")


ggsave(paste0("../Figures/Spatial_resolution/Model/Density_emp_sim_all_models.pdf"),p,width = 12,height = 10)



#Varying each metric along pooling intensity gradient

d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models.csv",sep=";"),
             Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
             PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
  filter(., rho_p>.05)

for (i in 1:(nrow(d_sim)/4)){
  d_sim$PL_expo[(4*(i-1)+1):(4*(i-1)+3)]=d_sim$PL_expo[(4*(i-1)+4)]
}
#d_sim$PL_expo[is.na(d_sim$PL_expo)]=0

p=d_sim%>%
  mutate(., Id_sim=rep(1:(nrow(.)/4),each=4))%>%
  melt(., id.vars=c("Model","Pooling","Id_sim"))%>%
  mutate(., Model=recode_factor(Model,"Eby_feedback"="Eby with feedbacks"))%>%
  group_by(., Model,Pooling,variable)%>%
  summarise(., .groups = "keep",mean_value=mean(value,na.rm = T))%>%
  ggplot(.)+
  geom_line(aes(x=Pooling,y=mean_value,group=interaction(Model),color=Model),lwd=1)+
  geom_point(aes(x=Pooling,y=mean_value,color=Model),fill="white",shape=21,size=2.5)+
  facet_wrap(.~variable,scales = "free")+
  labs(x="Change in resolution \n Landscape size",y="Mean value",color="Model")+
  scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))+
  scale_fill_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))+
  scale_x_discrete(labels=c("No change \n 50","x2 \n 100","x3 \n 150","x4 \n 200"))+
  the_theme+
  guides(color = guide_legend(override.aes = list(size = 2)),fill="none")
  
ggsave(paste0("../Figures/Spatial_resolution/Model/Change_metrics_spatial_resolution_models.pdf"),p,width = 8,height = 7)



for (i in c(.2, .3, .4)){
  
  d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models.csv",sep=";"),
               Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
               PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
    filter(., rho_p>i & rho_p<i+0.01)
  
  assign(paste0("p_",which(i==c(.2,.3,.4))),d_sim%>%
    mutate(., Id_sim=rep(1:(nrow(.)/4),each=4))%>%
    melt(., id.vars=c("Model","Pooling","Id_sim"))%>%
    mutate(., Model=recode_factor(Model,"Eby_feedback"="Eby with feedbacks"))%>%
    ggplot(.)+
    geom_line(aes(x=Pooling,y=value,group=interaction(Model,Id_sim),color=Model),lwd=.5,alpha=.5)+
    facet_wrap(.~variable,scales = "free")+
    labs(x="Change in resolution \n Landscape size",y="Mean value",color="Model")+
    scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))+
    scale_fill_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))+
    scale_x_discrete(labels=c("No change \n 50","x2 \n 100","x3 \n 150","x4 \n 200"))+
    the_theme+
    guides(color = guide_legend(override.aes = list(size = 3)),fill="none")
  )

}
ggsave(paste0("../Figures/Spatial_resolution/Model/Change_metrics_resolution_example_cover.pdf"),
       ggarrange(p_1,p_2,p_3,nrow=3),
       width = 8,height = 17)






#Trends of each metric along the stressor gradient for each model 

list_sim=list.files("../Data/Step9_Spatial_resolution/Trends/")
d_trends=tibble()

for (sim in list_sim){
  
  d_sim=read.table(paste0("../Data/Step9_Spatial_resolution/Trends/",sim),sep=",")%>%
    filter(., V1>.1) #10% of cover at least
  colnames(d_sim)=c("rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                    "Spectral_ratio","PLR","PL_expo","driver","Pooling")
  
  # print(ggplot(d_sim%>%melt(., id.vars=c("Pooling","driver"))%>%filter(., variable %!in% c("rho_p","nb_neigh","clustering")))+
  #   geom_point(aes(x=driver,y=value,color=as.character(Pooling)))+
  #   labs(x="Driver",y="")+
  #   facet_wrap(.~variable,scales="free")+
  #   the_theme)
    
  dat_melt=melt(d_sim, id.vars=c("Pooling","driver"))
  
  for (pool in unique(dat_melt$Pooling)){
    for (sumstat in unique(dat_melt$variable)){
      
      #fitting linear mixed effect model with site as random factor
      
      mod=lm(value~driver,data = filter(dat_melt,variable==sumstat,Pooling==pool)%>%
                 filter(., !is.na(value))%>%
                 mutate(., value=(value-mean(value))/sd(value))) #standardizing
      
      d_trends=rbind(d_trends,tibble(Sumstat=sumstat,Pooling=pool,Trend=mod$coefficients[2],  #getting the trend
                                     Model=gsub("Trends_","",gsub(".csv","",sim))))
      
    }
  }
}

p=d_trends%>%filter(., Sumstat %!in% c("rho_p","nb_neigh","clustering"))%>%
  mutate(.,Sign_trend=.$Trend>0,
         Value_pooling=sapply(1:nrow(.),function(x){
           if (.$Pooling[x]>1){return("Decreasing")
           }else if (.$Pooling[x]==1){return("Reference")
           }else {
             return("Increasing")
           }
         }))%>%
  mutate(., Pooling_intensity=sapply(1:nrow(.),function(x){
    if (.$Pooling[x] %in% c(.25,4)){return("High")
    }else if (.$Pooling[x] %in% c(.5,2)){return("Low")
    }else {
      return("Medium")
    }
  }))%>%
  ggplot(.)+
  geom_point(aes(x=Trend,y=Sumstat,color=Sign_trend,size=as.character(Pooling_intensity),shape=Value_pooling))+
  geom_vline(xintercept = 0,linetype=9)+
  the_theme+
  facet_wrap(.~Model)+
  scale_color_manual(values=c("red","blue"))+
  guides(color="none")+
  labs(x="Trend along aridity gradient",y="",shape="Change in resolution",size="Scale of change")+
  scale_y_discrete(labels=c("Moran I","PL exponent","PLR","Skewness","Spectral ratio","Variance"))+
  scale_shape_manual(values=c(11,8,7))+
  scale_size_manual(values=c(5,2,3.5))+
  theme(legend.box = "vertical")

ggsave("../Figures/Spatial_resolution/Model/Change_of_trends_spatial_resolution_model.pdf",p,width = 9,height = 8)



d_sim_all=tibble()
list_sim=list.files("../Data/Step9_Spatial_resolution/Trends/")

for (sim in list_sim){
  d_sim=read.table(paste0("../Data/Step9_Spatial_resolution/Trends/",sim),sep=",")%>%
    filter(., V1>.1)%>%add_column(.,Model=gsub("Trends_","",gsub(".csv","",sim))) #10% of cover at least
  colnames(d_sim)=c("rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                    "Spectral_ratio","PLR","PL_expo","driver","Pooling","Model")
  d_sim_all=rbind(d_sim_all,d_sim)
}


p=d_sim_all%>%
  mutate(.,Spectral_ratio=log(Spectral_ratio),clustering=log(clustering))%>%
  melt(., id.vars=c("driver","Pooling","Model"))%>%
  filter(., variable %!in% c("nb_neigh","clustering"))%>%
  ggplot(.)+
  geom_point(aes(x=driver,y=value,color=as.character(Pooling)),alpha=.5,size=2)+
  facet_grid(variable~Model,scales="free")+
  the_theme+
  labs(x="1 - driver",y="Summary statistics value",color="")+
  scale_colour_hue(labels=c("x4","x3","x2","Original","/2"))+
  guides(color = guide_legend(override.aes = list(size = 2)),fill="none")


ggsave("../Figures/Spatial_resolution/Model/Explicit_change_sumstat_aridity_model.pdf",p,width = 7,height = 8)







## Same but starting with landscapes of size 100 
## and doing both decrease/increase in pixel size

list_folders=list.dirs("../Data/Step9_Spatial_resolution/Simu2",recursive = F)
d_simu=tibble()

for (k in list_folders){
  
  if (tail(strsplit(k,"/")[[1]],1) %in% c("Eby","Eby_feedback")){
    n_param=2
  }else if (tail(strsplit(k,"/")[[1]],1)=="Kefi"){
    n_param=7
  }else {
    n_param=8
  }
  
  list_sim=list.files(k,".csv")
  
  d=tibble()
  for (sim in list_sim){
    d=rbind(d,read.table(paste0(k,"/",sim),sep=",")[,(n_param+1):(n_param+9)])
  }
  
  colnames(d)=c("rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                "Spectral_ratio","PLR","PL_expo")
  
  d=add_column(d,Pooling=rep(c("2",'1',"1/2"),nrow(d)/3))%>%
    add_column(., Model=tail(strsplit(k,"/")[[1]],1))%>%
    filter(., rho_p != 0)
  
  
  d_simu=rbind(d_simu,d)
}



d_sim=mutate(d_simu,
             Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
             PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%filter(., Pooling != "2")%>%
  add_column(., ID=rep(1:(nrow(.)/2),each=2))

n_keep=999

d_sim%>%
  sample_n(n_keep)%>%
  melt(., id.vars=c("Model","Pooling","ID"))%>%
  mutate(., Model=recode_factor(Model,"Eby_feedback"="Eby with feedbacks"))%>%
  ggplot(.)+
  geom_line(aes(x=Pooling,y=value,group=interaction(Model,ID),color=Model),lwd=.3,alpha=.5)+
  # geom_point(aes(x=Pooling,y=value,color=Model),fill="white",shape=21,size=2.5)+
  facet_wrap(.~variable,scales = "free")+
  labs(x="Change in resolution \n Landscape size",y="",color="Model")+
  scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))+
  scale_fill_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))+
  the_theme+
  guides(color = guide_legend(override.aes = list(size = 2)),fill="none")







## 3) >> Which are the outsiders in the PCA ? ----

d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")%>%
  filter(., Nbpixels<80000)

d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models.csv",sep=";"),
             Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
             PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
  filter(., rho_p>0.05)%>%filter(., Model=="Eby_feedback",Pooling=="1/2")

d_sim_data=rbind(d_sim,
                 d_biocom[,17:ncol(d_biocom)]%>% #keeping only the data with relatively low quality
                   add_column(., Pooling="Data",Model=NA)%>%
                   mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering)))

sumstat_name=colnames(d_sim_data)[1:9]
res.comp=imputePCA(d_sim_data[,which(colnames(d_sim_data) %in% sumstat_name)],ncp=3,scale = T) 

if ("completeObs" %in% names(res.comp)){
  res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
}else {
  res.pca=PCA(res.comp, ncp = 3,  graph=F)
}

#centroids in PCA of simulations
centroid_sim=colMeans(as.data.frame(res.pca$ind$coord)[1:nrow(d_sim),],na.rm = T)

threshold=.8
#getting euclidean distance from centroid for each empirical site
dist_empirical = sapply(1:nrow(d_biocom),function(x){
  return( sqrt(sum((res.pca$ind$coord[x+nrow(d_sim),] - centroid_sim)^2)) )
})

far_sites=as.numeric(rownames(d_biocom)[which(dist_empirical>quantile(dist_empirical,threshold))])


d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")%>%
  add_column(., Far=sapply(1:nrow(.),function(x){return(ifelse(x %in% far_sites,"T","F"))}))

ggplot(d_biocom%>%mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering))%>%
         melt(., measure.vars=colnames(d_biocom)[17:25]))+
  geom_density(aes(x=value,fill=Far),alpha=.4)+
  facet_wrap(.~variable,scales = "free")+
  theme_classic()



## 4) >> Example change fit PL exponent with resolution ----


id=86
plot_distr(spatialwarnings::patchdistr_sews(Get_empirical_site(id)==1,xmin = 1))
plot_distr(spatialwarnings::patchdistr_sews(pooling(Get_empirical_site(id),.5)==1,xmin = 4))
plot_distr(spatialwarnings::patchdistr_sews(pooling(Get_empirical_site(id),1/3)==1,xmin = 9))












# ---------------------- Step 6: Distinguishing between the different types of vegetation ----

# 0) Processing images


#We need to distinguish herbaceous from shrubs


#transforming into png images

list_raw=list.files("../Data/Data_Biocom/Raw_images/")
for (img_id in list_raw){
  img=image_read(paste0("../Data/Data_Biocom/Raw_images/",img_id))
  img=image_convert(img, format = "png")
  image_write(img, path = paste0("../Data/Data_Biocom/png_img/",gsub(".tif","",img_id),".png"))
}


d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")

#transforming each image into a 2, 3 or 4 clusters to differentiate vegetation

for (img_id in list.files("../Data/Data_Biocom/png_img/")){
}



## 1) >> Classification herb vs shrubs ----

d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")


for (i in 1:345){
  d=read.table("../Data/Data_Biocom/type_vege.csv",sep=";")
  
  d=rbind(d,tibble(Name=i,Type=Analyse_image(i)))
  
  write.table(d,"../Data/Data_Biocom/type_vege.csv",sep=";")
}




for (mod in c("Schneider","Eby_feedback")){
  
  d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
  d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models.csv",sep=";"),
               Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
               PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
    filter(., rho_p>0.03,Model==mod)
  
  for (i in 1:(nrow(d_sim)/4)){
    d_sim$PL_expo[(4*(i-1)+1):(4*(i-1)+3)]=d_sim$PL_expo[(4*(i-1)+4)]
  }
  
  d=read.table("../Data/Data_Biocom/type_vege.csv",sep=";")
  
  d_biocom$PL_expo[d_biocom$PL_expo==0]=NA
  
  d_sim_data=rbind(d_sim,
                   d_biocom[,17:ncol(d_biocom)]%>% #keeping only the data with relatively low quality
                     add_column(., Pooling="Data",Model=d$Type)%>%
                     filter(., d_biocom$Nbpixels<80000)%>%
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
             geom_point(aes(x = PC1, y = PC2, color = Pooling,fill=Pooling,shape=Model,size=Model),alpha=.5)+
             scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","black"))+
             scale_fill_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","black"))+
             scale_size_manual(values=c(.5,2,2,2,2))+
             labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="Change in resolution",fill="")+
             ggtitle("")+
             theme_classic()+theme(legend.position = "bottom",legend.box = "vertical")+
             guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none",
                    shape = guide_legend(override.aes = list(size = 5)))
    )
  }
  
  p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                        p2+theme(legend.position = "none"),
                        p3+theme(legend.position = "none"),
                        ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
              nrow=2,heights = c(1,.2))
  
  ggsave(paste0("../Figures/Plant_type/PCA_shrub_herb_coex_",mod,".pdf"),p,width = 11,height = 6)
  
}




d_sim_data=rbind(d_sim%>%mutate(., Model=NA, Pooling=NA),
                 d_biocom[,17:ncol(d_biocom)]%>% #keeping only the data with relatively low quality
                   add_column(., Pooling="Data",Model=log(d_biocom$Nbpixels))%>%
                   filter(d_biocom$Nbpixels<80000)%>%
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
           geom_point(aes(x = PC1, y = PC2, color = Model,fill=Model))+
           scale_color_viridis_c(alpha = 1,na.value = alpha("gray",.3))+
           scale_fill_viridis_c(alpha = 1,na.value = alpha("gray",.3))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="# Pixels",fill="")+
           ggtitle("")+guides(shape="none")+
           theme_classic()+theme(legend.position = "bottom")+
           guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")
  )
}

p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                      p2+theme(legend.position = "none"),
                      p3+theme(legend.position = "none"),
                      ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.2))

ggsave(paste0("../Figures/Plant_type/PCA_NB_pixels.pdf"),p,width = 11,height = 6)





d_sim_data=rbind(d_sim%>%mutate(., Model=NA, Pooling=NA),
                 d_biocom[,17:ncol(d_biocom)]%>% #keeping only the data with relatively low quality
                   add_column(., Pooling="Data",Model=log(d_biocom$clustering))%>%
                   filter(d_biocom$Nbpixels<80000)%>%
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
           geom_point(aes(x = PC1, y = PC2, color = Model,fill=Model))+
           scale_color_viridis_c(alpha = 1,na.value = alpha("gray",.3))+
           scale_fill_viridis_c(alpha = 1,na.value = alpha("gray",.3))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="Clustering",fill="")+
           ggtitle("")+guides(shape="none")+
           theme_classic()+theme(legend.position = "bottom")+
           guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")
  )
}

p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                      p2+theme(legend.position = "none"),
                      p3+theme(legend.position = "none"),
                      ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.2))

ggsave(paste0("../Figures/Plant_type/PCA_clustering.pdf"),p,width = 11,height = 6)






d_sim_data=rbind(d_sim%>%mutate(., Model=NA,Pooling=NA),
                 d_biocom[,17:ncol(d_biocom)]%>% #keeping only the data with relatively low quality
                   add_column(., Pooling="Data",Model=(d_biocom$nb_neigh))%>%
                   filter(., d_biocom$Nbpixels<80000)%>%
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
           geom_point(aes(x = PC1, y = PC2, color = Model,fill=Model))+
           scale_color_viridis_c(alpha = 1,na.value = alpha("gray",.3))+
           scale_fill_viridis_c(alpha = 1,na.value = alpha("gray",.3))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),color="# neighbors",fill="")+
           ggtitle("")+guides(shape="none")+
           theme_classic()+theme(legend.position = "bottom")+
           guides(color = guide_legend(override.aes = list(size = 3)),fill="none",size="none")
  )
}

p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                      p2+theme(legend.position = "none"),
                      p3+theme(legend.position = "none"),
                      ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.2))

ggsave(paste0("../Figures/Plant_type/PCA_nb_neighbors.pdf"),p,width = 11,height = 6)



## 2) >> Statistics on shrubs -----
library(spatstat)

list_landscape=list.files("../Data/Data_Biocom/type_landscape/",pattern = "type2") #we keep shrubs
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")

d=tibble()

for (i in 1:length(list_landscape)){
  
  mat=as.matrix(read.table(paste0("../Data/Data_Biocom/type_landscape/",list_landscape[i]),sep=";"))
  # Plot_landscape(mat)
  mat_filter=Filtering_small_patches(mat,cutoff = 40)
  # Plot_landscape(mat_filter)
  d2=Centroid_patches(mat_filter) #getting the centroïd and size of each patch
  d=rbind(d,d2%>%add_column(., Id_landscape=gsub(".csv","",strsplit(list_landscape[i],"_")[[1]][3]))) 
  
}


ggplot(d)+
  geom_histogram(aes(x=log(Size)),fill=alpha("blue",.3),color="blue")+
  theme_classic()


pdf("../Figures/Plant_type/Spatial_distrib_big_patches.pdf",width = 6,height = 4)
for (i in unique(d$Id_landscape)){
  
  d_fil=filter(d,Id_landscape==i)
  dat=as.ppp(data.frame(x=d_fil$centroid_x,y=d_fil$centroid_y),c(0,unique(d_fil$window),0,unique(d_fil$window)))

  Kdot=envelope(dat, Kest, nsim = 99)
  plot(Kdot)
  
}
dev.off()



# ---------------------- Step 7: Testing ABC with the observation parameter: can we retrieve parameters despite the scale of observation ----

## 0) >> Aggregating the simulations ----

list_sim=list.files("../Data/Step10_ABC_scale/Sim",pattern = ".csv")
d=tibble()

for (sim in list_sim){

  d2=read.table(paste0("../Data/Step10_ABC_scale/Sim/",sim),sep=",")
  colnames(d2)= c("p","q","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                  "Spectral_ratio","PLR","PL_expo")
  d=rbind(d,filter(d2,rho_p != 0)%>%
            add_column(.,Pooling=rep(c("x2","x3"),nrow(.)/2)))
}
write.table(d,"../Data/Step10_ABC_scale/All_sim_Eby.csv",sep=";")



## 1) >> ABC with multiple scale of observation ----

d_simu=tibble()
list_sim=list.files("../Data/Step9_Spatial_resolution/Simu/Eby_feedback",".csv")
n_param=2

d=tibble()
for (sim in list_sim){
  d=rbind(d,read.table(paste0("../Data/Step9_Spatial_resolution/Simu/Eby_feedback","/",sim),sep=",")[,1:(n_param+9)])
}

colnames(d)=c("p","q","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
              "Spectral_ratio","PLR","PL_expo")

d_sim=add_column(d,Pooling=rep(c("1/4","1/3",'1/2',"1"),nrow(d)/4))%>%
  add_column(., Model="Eby_feedback")%>%
  filter(., rho_p != 0)%>%
  filter(., Pooling !="1")%>%
  mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
         PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
  filter(., rho_p>0.03)
  
list_sample=sample(1:nrow(d_sim),300)
d_RMSE_param=x_y_param=tibble()

for (sample_id in list_sample){
  
  target=d_sim[sample_id,3:11]
  
  if (any(is.na(target))){
    
    which_na=which(is.na(target))
    abc_sim=abc(target=target[-which_na],
                param = d_sim[-sample_id,1:2],sumstat = d_sim[-sample_id,c(3:11)[-which_na]],
                tol = 50/nrow(d_sim),method="rejection")
    
  }else {
    abc_sim=abc(target=target,
                param = d_sim[-sample_id,1:2],sumstat = d_sim[-sample_id,3:11],
                tol = 50/nrow(d_sim),method="rejection")
  }
  
  RMSE = sapply(1:ncol(abc_sim$unadj.values),function(x){
    sqrt(sum((abc_sim$unadj.values[,x]-as.numeric(d_sim[sample_id,1:2])[x])**2)/nrow(abc_sim$unadj.values) )})
  
  RMSE_prior=sapply(1:ncol(d_sim[,1:2]),function(x){
    sqrt(sum((d_sim[,x]-as.numeric(d_sim[sample_id,1:2])[x])**2)/nrow(d_sim) )})
  
  NMRSE=RMSE/RMSE_prior

  d_RMSE_param=rbind(d_RMSE_param,as_tibble(t(NMRSE))%>%add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection"))
  
  x_y_param=rbind(x_y_param,as_tibble(t(colMeans(abc_sim$unadj.values)))%>%
                    add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection",Type="Sim"))
  
  x_y_param=rbind(x_y_param,d_sim[sample_id,1:2]%>%
                    add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection",Type="Obs"))
  
  
  if (any(is.na(target))){
    
    which_na=which(is.na(target))
    abc_sim=abc(target=target[-which_na],
                param = d_sim[-sample_id,1:2],sumstat = d_sim[-sample_id,c(3:11)[-which_na]],
                tol = 50/nrow(d_sim),method="neuralnet",transf = rep("logit",2),
                logit.bounds = matrix(c(0,1),2,2,byrow = T),
                numnet = 10,sizenet = 10)
    
  }else {
    abc_sim=abc(target=target,
                param = d_sim[-sample_id,1:2],sumstat = d_sim[-sample_id,3:11],
                tol = 50/nrow(d_sim),method="neuralnet",transf = rep("logit",2),
                logit.bounds = matrix(c(0,1),2,2,byrow = T),
                numnet = 10,sizenet = 10)
  }
  
  RMSE = sapply(1:ncol(abc_sim$adj.values),function(x){
    sqrt(sum((abc_sim$adj.values[,x]-as.numeric(d_sim[sample_id,1:2])[x])**2)/nrow(abc_sim$adj.values) )})
  
  RMSE_prior=sapply(1:ncol(d_sim[,1:2]),function(x){
    sqrt(sum((d_sim[,x]-as.numeric(d_sim[sample_id,1:2])[x])**2)/nrow(d_sim) )})
  
  NMRSE=RMSE/RMSE_prior
  
  d_RMSE_param=rbind(d_RMSE_param,as_tibble(t(NMRSE))%>%add_column(., Site_ID=which(list_sample == sample_id),Method="NeuralNet"))
  
  x_y_param=rbind(x_y_param,as_tibble(t(colMeans(abc_sim$adj.values)))%>%
    add_column(., Site_ID=which(list_sample == sample_id),Method="NeuralNet",Type="Sim"))
  
  x_y_param=rbind(x_y_param,d_sim[sample_id,1:2]%>%
    add_column(., Site_ID=which(list_sample == sample_id),Method="NeuralNet",Type="Obs"))
  
}

# write.table(d_RMSE_param,"../Data/Step10_ABC_scale/RMSE_stat_Eby_feedback_merged_scales.csv",sep=";")
# write.table(x_y_param,"../Data/Step10_ABC_scale/x_y_obs_sim_param.csv",sep=";")

d_RMSE_param=read.table("../Data/Step10_ABC_scale/RMSE_stat_Eby_feedback_merged_scales.csv",sep=";")

p=ggplot(d_RMSE_param%>%rename(., p=V1,q=V2)%>%melt(., id.vars=c("Site_ID","Method")))+
  geom_jitter(aes(x=Method,y=value,color=Method),alpha=.5,size=.5,position = position_jitter(height = 0,width = .1))+
  facet_wrap(.~variable,labeller = label_bquote(cols = Parameter==.(as.character(variable))))+
  labs(x="Method ABC (post-processing)",y="NRMSE",color="")+
  the_theme+
  guides(Method="none")+
  scale_color_manual(values=c("purple","green"))

ggsave("../Figures/ABC_scale/ABC_param_all_scales.pdf",p,width = 6,height = 4)


x_y_param=read.table("../Data/Step10_ABC_scale/x_y_obs_sim_param.csv",sep=";")

par(mfrow=c(2,2))
plot(x=filter(x_y_param,Method=="Rejection",Type=="Sim")[,1],xlab="Sim",ylab="Obs",main="Rejection, p",
     y=filter(x_y_param,Method=="Rejection",Type=="Obs")[,1],col="gray")
abline(a=0,b=1)

plot(x=filter(x_y_param,Method=="Rejection",Type=="Sim")[,2],xlab="Sim",ylab="Obs",main="Rejection, q",
      y=filter(x_y_param,Method=="Rejection",Type=="Obs")[,2],col="gray")
abline(a=0,b=1)
plot(x=filter(x_y_param,Method=="NeuralNet",Type=="Sim")[,1],xlab="Sim",ylab="Obs",main="NeuralNet, p",
      y=filter(x_y_param,Method=="NeuralNet",Type=="Obs")[,1],col="gray")
abline(a=0,b=1)
plot(x=filter(x_y_param,Method=="NeuralNet",Type=="Sim")[,2],xlab="Sim",ylab="Obs",main="NeuralNet, q",
      y=filter(x_y_param,Method=="NeuralNet",Type=="Obs")[,2],col="gray")
abline(a=0,b=1)


## 2) >> ABC with only x2, x3, x4 as scale of observation (separately) ----

d_simu=tibble()
list_sim=list.files("../Data/Step9_Spatial_resolution/Simu/Eby_feedback",".csv")
n_param=2

d=tibble()
for (sim in list_sim){
  d=rbind(d,read.table(paste0("../Data/Step9_Spatial_resolution/Simu/Eby_feedback","/",sim),sep=",")[,1:(n_param+9)])
}

colnames(d)=c("p","q","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
              "Spectral_ratio","PLR","PL_expo")
for (i in 1:(nrow(d)/4)){
  d$PL_expo[(4*(i-1)+1):(4*(i-1)+3)]=d$PL_expo[(4*(i-1)+4)]
}


d_RMSE_param=x_y_param=d_RMSE_stat=tibble()
set.seed(123)
list_sample=sample(1:5559,100)

for (scale_obs in c("1/2","1/3","1/4")){

  d_sim=add_column(d,Pooling=rep(c("1/4","1/3",'1/2',"1"),nrow(d)/4))%>%
    add_column(., Model="Eby_feedback")%>%
    filter(., rho_p != 0)%>%
    filter(., Pooling ==scale_obs)%>%
    mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
           PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
    filter(., rho_p>0.03)
  
  
  
  for (sample_id in list_sample){ #we use the same landscape to compare what we are inferring
    
    target=d_sim[sample_id,3:11]
    
    if (any(is.na(target))){
      which_na=which(is.na(target))
      abc_sim=abc(target=target[-which_na],
                  param = d_sim[-sample_id,1:2],sumstat = d_sim[-sample_id,c(3:11)[-which_na]],
                  tol = 50/nrow(d_sim),method="rejection")
    } else {
      abc_sim=abc(target=target,
                  param = d_sim[-sample_id,1:2],sumstat = d_sim[-sample_id,3:11],
                  tol = 50/nrow(d_sim),method="rejection")
    }
    
    RMSE = sapply(1:ncol(abc_sim$unadj.values),function(x){
      sqrt(sum((abc_sim$unadj.values[,x]-as.numeric(d_sim[sample_id,1:2])[x])**2)/nrow(abc_sim$unadj.values) )})
    
    RMSE_prior=sapply(1:ncol(d_sim[,1:2]),function(x){
      sqrt(sum((d_sim[,x]-as.numeric(d_sim[sample_id,1:2])[x])**2)/nrow(d_sim) )})
    
    NMRSE=RMSE/RMSE_prior
    
    d_RMSE_param=rbind(d_RMSE_param,as_tibble(t(NMRSE))%>%
                         add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection",Pooling=scale_obs))
    x_y_param=rbind(x_y_param,d_sim[sample_id,1:2]%>%
                      add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection",Pooling=scale_obs,Type="Obs"))
    x_y_param=rbind(x_y_param,as_tibble(t(colMeans(abc_sim$unadj.values)))%>%
                      add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection",Pooling=scale_obs,Type="Sim"))
    
    RMSE = sapply(1:ncol(abc_sim$ss),function(x){
      sqrt(sum((abc_sim$ss[,x]-as.numeric(target)[x])**2)/nrow(abc_sim$ss) )})
    
    RMSE_prior=sapply(1:9,function(x){
      sqrt(sum((d_sim[,2+x]-as.numeric(target)[x])**2)/nrow(d_sim) )})
    
    NMRSE=RMSE/RMSE_prior
    
    d_RMSE_stat=rbind(d_RMSE_stat,as_tibble(t(NMRSE))%>%
                        add_column(., Site_ID=which(list_sample == sample_id),Method="NeuralNet",Pooling=scale_obs))
    
    
    if (any(is.na(target))){
      
      which_na=which(is.na(target))
      abc_sim=abc(target=target[-which_na],
                  param = d_sim[-sample_id,1:2],sumstat = d_sim[-sample_id,c(3:11)[-which_na]],
                  tol = 50/nrow(d_sim),method="neuralnet",transf = rep("logit",2),
                  logit.bounds = matrix(c(0,1),2,2,byrow = T),
                  numnet = 10,sizenet = 10)
      
    }else {
      abc_sim=abc(target=target,
                  param = d_sim[-sample_id,1:2],sumstat = d_sim[-sample_id,3:11],
                  tol = 50/nrow(d_sim),method="neuralnet",transf = rep("logit",2),
                  logit.bounds = matrix(c(0,1),2,2,byrow = T),
                  numnet = 10,sizenet = 10)
    }
    
    #param
    RMSE = sapply(1:ncol(abc_sim$adj.values),function(x){
      sqrt(sum((abc_sim$adj.values[,x]-as.numeric(d_sim[sample_id,1:2])[x])**2)/nrow(abc_sim$adj.values) )})
    
    RMSE_prior=sapply(1:ncol(d_sim[,1:2]),function(x){
      sqrt(sum((d_sim[,x]-as.numeric(d_sim[sample_id,1:2])[x])**2)/nrow(d_sim) )})
    
    NMRSE=RMSE/RMSE_prior
    
    d_RMSE_param=rbind(d_RMSE_param,as_tibble(t(NMRSE))%>%
                         add_column(., Site_ID=which(list_sample == sample_id),Method="NeuralNet",Pooling=scale_obs))
  }
}

# write.table(d_RMSE_param,"../Data/Step10_ABC_scale/RMSE_param_Eby_feedback_scales.csv",sep=";")
# write.table(d_RMSE_stat,"../Data/Step10_ABC_scale/RMSE_stat_Eby_feedback_scales.csv",sep=";")


d_RMSE_stat=read.table("../Data/Step10_ABC_scale/RMSE_stat_Eby_feedback_scales.csv",sep=";")

colnames(d_RMSE_stat)[1:9]=c("rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                        "Spectral_ratio","PLR","PL_expo")

p1=ggplot(d_RMSE_stat%>%melt(., id.vars=c("Site_ID","Method","Pooling"))%>%
         mutate(., Pooling=recode_factor(Pooling,"1/2"="x2","1/3"="x3","1/4"="x4")))+
  geom_jitter(aes(x=Pooling,y=value,color=Pooling),alpha=.3,position = position_jitter(height = 0,width = .1))+
  facet_wrap(.~variable,scales = "free")+
  labs(x="Scale of observation",y="NRMSE",color="")+
  the_theme+
  guides(Method="none")+
  scale_color_manual(values=c("purple","green","orange"))

ggsave("../Figures/Spatial_resolution/ABC_scale/ABC_stat_Eby_feedback.pdf",p1,width = 7,height = 6)


d_RMSE_param=read.table("../Data/Step10_ABC_scale/RMSE_param_Eby_feedback_scales.csv",sep=";")
colnames(d_RMSE_param)[1:2]=c("p","q")

p2=ggplot(d_RMSE_param%>%melt(., id.vars=c("Site_ID","Method","Pooling"))%>%
         mutate(., Pooling=recode_factor(Pooling,"1/2"="x2","1/3"="x3","1/4"="x4")))+
  geom_jitter(aes(x=Pooling,y=value,color=Pooling),alpha=.3,position = position_jitter(height = 0,width = .1))+
  facet_wrap(.~variable,scales = "free")+
  labs(x="Scale of observation",y="NRMSE",color="")+
  the_theme+
  guides(Method="none")+
  scale_color_manual(values=c("purple","green","orange"))

ggsave("../Figures/Spatial_resolution/ABC_scale/ABC_param_Eby_feedback.pdf",p2,width = 7,height = 6)


## 3) >> Quantifying the robustness of inference across scales of observation ----
# 
# d_simu=tibble()
# list_sim=list.files("../Data/Step9_Spatial_resolution/Simu/Eby_feedback",".csv")
# n_param=2
# 
# d=tibble()
# for (sim in list_sim){
#   d=rbind(d,read.table(paste0("../Data/Step9_Spatial_resolution/Simu/Eby_feedback","/",sim),sep=",")[,1:(n_param+9)])
# }
# colnames(d)=c("p","q","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
#               "Spectral_ratio","PLR","PL_expo")
# 
# for (i in 1:(nrow(d)/4)){
#   d$PL_expo[(4*(i-1)+1):(4*(i-1)+3)]=d$PL_expo[(4*(i-1)+4)]
# }
# 
# set.seed(123)
# 
# list_virtual=d%>%add_column(.,Pooling=rep(c("1/4","1/3",'1/2',"1"),nrow(d)/4))%>%
#   mutate(., ID_row=1:nrow(.))%>%
#   filter(., rho_p != 0,Pooling !="1")%>%
#   mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
#          PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
#   filter(., rho_p>0.03)%>%
#   sample_n(., 100)
# 
# d_RMSE_param=x_y_param=d_RMSE_stat=tibble()
# 
# for (scale_obs in c("1/2","1/3","1/4","all")){
#   
#   for (sample_id in 1:nrow(list_virtual)){ #we use the same landscape to compare what we are inferring
#     
#     if (scale_obs=="all"){scale_obs=c("1/2","1/3","1/4")}
#     
#     d_sim=add_column(d,Pooling=rep(c("1/4","1/3",'1/2',"1"),nrow(d)/4))%>%
#       add_column(., Model="Eby_feedback",ID_row=1:nrow(.))%>%
#       filter(., rho_p != 0)%>%
#       filter(., Pooling %in% scale_obs)%>%
#       mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
#              PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
#       filter(., rho_p>0.03)
#     
#     
#     target=list_virtual[sample_id,3:11]
#     target_param=list_virtual[sample_id,1:2]
#       
#     if (any(which(d_sim$ID_row==list_virtual$ID_row[sample_id]))){ #in case virtual data is in the simulation df, we remove it
#       d_sim=d_sim[-which(d_sim$ID_row==list_virtual$ID_row[sample_id]),]
#     }
#     
#     
#     if (any(is.na(target))){
#       which_na=which(is.na(target))
#       abc_sim=abc(target=target[-which_na],
#                   param = d_sim[,1:2],
#                   sumstat = d_sim[,c(3:11)[-which_na]],
#                   tol = 50/nrow(d_sim),method="rejection")
#     } else {
#       abc_sim=abc(target=target,
#                   param = d_sim[,1:2],
#                   sumstat = d_sim[,3:11],
#                   tol = 50/nrow(d_sim),method="rejection")
#     }
#     
# 
#     if (length(scale_obs)==3){scale_obs="all"}
#     
#     x_y_param=rbind(x_y_param,list_virtual[sample_id,1:2]%>%
#                       add_column(., Site_ID=sample_id,Method="Rejection",Pooling=scale_obs,Type="Obs"))
#     x_y_param=rbind(x_y_param,as_tibble(t(colMeans(abc_sim$unadj.values)))%>%
#                       add_column(., Site_ID=sample_id,Method="Rejection",Pooling=scale_obs,Type="Sim"))
#     
#     if (any(is.na(target))){
#       
#       which_na=which(is.na(target))
#       abc_sim=abc(target=target[-which_na],
#                   param = d_sim[,1:2],sumstat = d_sim[,c(3:11)[-which_na]],
#                   tol = 50/nrow(d_sim),method="neuralnet",transf = rep("logit",2),
#                   logit.bounds = matrix(c(0,1),2,2,byrow = T),
#                   numnet = 10,sizenet = 10)
#       
#     }else {
#       abc_sim=abc(target=target,
#                   param = d_sim[,1:2],sumstat = d_sim[,3:11],
#                   tol = 50/nrow(d_sim),method="neuralnet",transf = rep("logit",2),
#                   logit.bounds = matrix(c(0,1),2,2,byrow = T),
#                   numnet = 10,sizenet = 10)
#     }
#     
# 
#     x_y_param=rbind(x_y_param,list_virtual[sample_id,1:2]%>%
#                       add_column(., Site_ID=sample_id,Method="NeuralNet",Pooling=scale_obs,Type="Obs"))
#     x_y_param=rbind(x_y_param,as_tibble(t(colMeans(abc_sim$adj.values)))%>%
#                       add_column(., Site_ID=sample_id,Method="NeuralNet",Pooling=scale_obs,Type="Sim"))
#     
#   }
# }
# 
# write.table(x_y_param,"../Data/Step10_ABC_scale/Quality_inference_scale.csv",sep=";")
# 
# 
# 
# x_y_param=read.table("../Data/Step10_ABC_scale/Quality_inference_scale.csv",sep=";")
# 
# ggplot(x_y_param%>%melt(., measure.vars=c("p","q"))%>%
#          filter(., Method=="NeuralNet"))+
#   geom_line(aes(x=Pooling,y=value,group=Site_ID),color="gray")+
#   the_theme+
#   facet_wrap(.~variable)
# 
# par(mfrow=c(3,3))
# for (i in unique(x_y_param$Method)){
#   for (j in unique(x_y_param$Pooling)){
#     plot(x=filter(x_y_param,Method==i,Type=="Sim",Pooling==j)[,1],xlab="Sim",ylab="Obs",main="Rejection, p",
#          y=filter(x_y_param,Method==i,Type=="Obs",Pooling==j)[,1],col="gray")
#     abline(a=0,b=1)
#   }
# }

## 4) >> Comparing the inference for similar landscapes but observed at different scales ----

d_simu=tibble()
list_sim=list.files("../Data/Step9_Spatial_resolution/Simu/Eby_feedback",".csv")
n_param=2

d=tibble()
for (sim in list_sim){
  d=rbind(d,read.table(paste0("../Data/Step9_Spatial_resolution/Simu/Eby_feedback","/",sim),sep=",")[,1:(n_param+9)])
}

colnames(d)=c("p","q","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
              "Spectral_ratio","PLR","PL_expo")

d_sim=add_column(d,Pooling=rep(c("1/4","1/3",'1/2',"1"),nrow(d)/4))%>%
  add_column(., Model="Eby_feedback",ID_sim=rep(1:(nrow(.)/4),each=4))%>%
  filter(., rho_p != 0)%>%
  filter(., Pooling !="1")%>%
  mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
         PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
  filter(., rho_p>0.03)

set.seed(123)
list_sample=sample(unique(d_sim$ID_sim),100) #for each of these 100 parameter sets, we need to verify whether we infer the same parameters of not
d_RMSE_param=x_y_param=tibble()

for (sample_id in list_sample){
  
  for (scale_obs in c("1/4","1/3","1/2")){
    
    target=d_sim%>%filter(., ID_sim==sample_id,Pooling==scale_obs)%>%
      dplyr::select(., -Pooling,-Model,-ID_sim)
    
    all_sim=d_sim[-which(d_sim$Pooling==scale_obs & d_sim$ID_sim==sample_id),]%>%
                   dplyr::select(., -Pooling,-Model,-ID_sim)
    
    if (any(is.na(target))){
      
      which_na=which(is.na(target))
      abc_sim=abc(target=target[-which_na,3:11],
                  param = all_sim[,1:2],sumstat = all_sim[,c(3:11)[-which_na]],
                  tol = 50/nrow(all_sim),method="rejection")
      
    }else {
      abc_sim=abc(target=target[,3:11],
                  param = all_sim[,1:2],sumstat = all_sim[,3:11],
                  tol = 50/nrow(all_sim),method="rejection")
    }
    
    RMSE = sapply(1:ncol(abc_sim$unadj.values),function(x){
      sqrt(sum((abc_sim$unadj.values[,x]-target[,x])**2)/nrow(abc_sim$unadj.values) )})
    
    RMSE_prior=sapply(1:2,function(x){
      sqrt(sum(((all_sim[,x]-target[,x])**2)/nrow(all_sim) ))})
    
    NMRSE=RMSE/RMSE_prior
    
    d_RMSE_param=rbind(d_RMSE_param,data.frame((t(NMRSE)))%>%
                         add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection",Pooling=scale_obs))
    
    x_y_param=rbind(x_y_param,data.frame((t(colMeans(abc_sim$unadj.values))))%>%
                      add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection",Type="Sim",Pooling=scale_obs))
    
    x_y_param=rbind(x_y_param,as.data.frame(target[,1:2])%>%
                      add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection",Type="Obs",Pooling=scale_obs))
    
    
    if (any(is.na(target))){
      
      which_na=which(is.na(target))
      abc_sim=abc(target=target[-which_na,3:11],
                  param = all_sim[,1:2],sumstat = all_sim[,c(3:11)[-which_na]],
                  tol = 50/nrow(all_sim),method="neuralnet",transf = rep("logit",2),
                  logit.bounds = matrix(c(0,1),2,2,byrow = T),
                  numnet = 10,sizenet = 10)
      
    }else {
      abc_sim=abc(target=target[,3:11],
                  param = all_sim[,1:2],sumstat = all_sim[,3:11],
                  tol = 50/nrow(all_sim),method="neuralnet",transf = rep("logit",2),
                  logit.bounds = matrix(c(0,1),2,2,byrow = T),
                  numnet = 10,sizenet = 10)
    }
    
    RMSE = sapply(1:ncol(abc_sim$unadj.values),function(x){
      sqrt(sum((abc_sim$unadj.values[,x]-target[,x])**2)/nrow(abc_sim$unadj.values) )})
    
    RMSE_prior=sapply(1:2,function(x){
      sqrt(sum(((all_sim[,x]-target[,x])**2)/nrow(all_sim) ))})
    
    NMRSE=RMSE/RMSE_prior
    
    d_RMSE_param=rbind(d_RMSE_param,data.frame((t(NMRSE)))%>%
                         add_column(., Site_ID=which(list_sample == sample_id),Method="NeuralNet",Pooling=scale_obs))
    
    x_y_param=rbind(x_y_param,data.frame(t(colMeans(abc_sim$adj.values)))%>%
                      add_column(., Site_ID=which(list_sample == sample_id),Method="NeuralNet",Type="Sim",Pooling=scale_obs))
    
    x_y_param=rbind(x_y_param,target[,1:2]%>%
                      add_column(., Site_ID=which(list_sample == sample_id),Method="NeuralNet",Type="Obs",Pooling=scale_obs))
  }
}

#write.table(d_RMSE_param,"../Data/Step10_ABC_scale/Retrieving_parameters_different_resolution_RMSE_param.csv",sep=";")
#write.table(x_y_param,"../Data/Step10_ABC_scale/Retrieving_parameters_different_resolution_x_y.csv",sep=";")

p1=ggplot(d_RMSE_param%>%
           melt(., id.vars=c("Site_ID","Method","Pooling")))+
  geom_jitter(aes(x=Pooling,y=value,color=Pooling),alpha=.5,size=.5,position = position_jitter(height = 0,width = .1))+
  facet_wrap(Method~variable,labeller = label_bquote(cols = Parameter==.(as.character(variable))))+
  labs(x="Method ABC (post-processing)",y="NRMSE",color="")+
  the_theme+
  guides(Method="none")+
  scale_color_manual(values=c("purple","green","red"))



p2=ggplot(x_y_param%>%
           melt(., id.vars=c("Site_ID","Method","Pooling","Type"))%>%
         filter(., Type=="Sim"))+
  geom_line(aes(x=Pooling,y=value,group=Site_ID),alpha=.5,lwd=.5,color="gray")+
  facet_grid(Method~variable,labeller = label_bquote(cols = Parameter==.(as.character(variable))))+
  labs(x="Method ABC (post-processing)",y="NRMSE",color="")+
  the_theme+
  guides(Method="none")+
  scale_color_manual(values=c("purple","green"))


ggsave("../Figures/ABC_scale/Consistency_inference_param_scale.pdf",
       #ggarrange(p1,p2,heights =c(1,1.1) ,labels=letters[1:2],nrow=2),
       p2,
       width = 7,height = 5)






# ---------------------- Step 8: Trying to infer the stability of empirical sites ------------------



for (model in c("Kefi","Eby_feedback","Eby","Schneider")){
  
  
  n_param=ifelse(model %in% c("Eby","Eby_feedback"),2,ifelse(model=="Kefi",7,8))
  
  
  d_simu=tibble()
  list_sim=list.files(paste0("../Data/Step9_Spatial_resolution/Simu/",model),".csv")
  
  d=tibble()
  for (sim in list_sim){
    d=rbind(d,read.table(paste0("../Data/Step9_Spatial_resolution/Simu/",model,"/",sim),sep=",")[,1:(n_param+9)])
  }
  
  if (model %in% c("Eby","Eby_feedback")){
    colnames(d)=c("p","q","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                  "Spectral_ratio","PLR","PL_expo")
  }else if(model =="Kefi"){
    colnames(d)=c("r","d","f","m","b","c","delta","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                  "Spectral_ratio","PLR","PL_expo")
  }else{
    colnames(d)=c("r","d","f","m","b","c","delta","g0","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
                  "Spectral_ratio","PLR","PL_expo")
  }
  
  d_sim=add_column(d,Pooling=rep(c("1/4","1/3",'1/2',"1"),nrow(d)/4))%>%
    add_column(., Model=model)%>%
    filter(., rho_p != 0)%>%
    mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
           PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
    filter(., rho_p>0.03)
  
  
  d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")%>%
    mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering))
  
  d_biocom$PL_expo[d_biocom$PL_expo==0]=NA
  
  
  
  d_param_infer=x_y_stat=tibble()
  
  for (empirical_id in 1:nrow(d_biocom)){
    
    target=d_biocom[empirical_id,17:25]
    
    if (any(is.na(target))){
      
      which_na=which(is.na(target))
      abc_sim=abc(target=target[-which_na],
                  param = d_sim[,1:n_param],sumstat = d_sim[,c((n_param+1):(n_param+9))[-which_na]],
                  tol = 50/nrow(d_sim),method="rejection")
      
      x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(abc_sim$ss)))%>%
                       add_column(.,PL_expo=NA, Site_ID=empirical_id,Method="Rejection",Type="Sim"))
      
      x_y_stat=rbind(x_y_stat,target%>%
                       add_column(., Site_ID=empirical_id,Method="Rejection",Type="Obs"))
      
      
      
      
    }else {
      abc_sim=abc(target=target,
                  param = d_sim[,1:n_param],sumstat = d_sim[,c((n_param+1):(n_param+9))],
                  tol = 50/nrow(d_sim),method="rejection")
      
      x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(abc_sim$ss)))%>%
                       add_column(., Site_ID=empirical_id,Method="Rejection",Type="Sim"))
      
      x_y_stat=rbind(x_y_stat,target%>%
                       add_column(., Site_ID=empirical_id,Method="Rejection",Type="Obs"))
      
    }
    
    
    d_param_infer=rbind(d_param_infer,as_tibble(t(colMeans(abc_sim$unadj.values)))%>%
                          add_column(., Site_ID=empirical_id,Method="Rejection"))
    
    
    
    
    if (any(is.na(target))){
      
      which_na=which(is.na(target))
      abc_sim=abc(target=target[-which_na],
                  param = d_sim[,1:n_param],sumstat = d_sim[,c((n_param+1):(n_param+9))[-which_na]],
                  tol = 50/nrow(d_sim),method="neuralnet",transf = rep("logit",n_param),
                  logit.bounds = matrix(c(0,1),n_param,2,byrow = T),
                  numnet = 10,sizenet = 10)
      
      
      x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(abc_sim$ss)))%>%
                       add_column(.,PL_expo=NA, Site_ID=empirical_id,Method="NeuralNet",Type="Sim"))
      
      x_y_stat=rbind(x_y_stat,target%>%
                       add_column(., Site_ID=empirical_id,Method="NeuralNet",Type="Obs"))
      
    }else {
      abc_sim=abc(target=target,
                  param = d_sim[,1:n_param],sumstat = d_sim[,c((n_param+1):(n_param+9))],
                  tol = 50/nrow(d_sim),method="neuralnet",transf = rep("logit",n_param),
                  logit.bounds = matrix(c(0,1),n_param,2,byrow = T),
                  numnet = 10,sizenet = 10)
      
      x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(abc_sim$ss)))%>%
                       add_column(., Site_ID=empirical_id,Method="NeuralNet",Type="Sim"))
      
      x_y_stat=rbind(x_y_stat,target%>%
                       add_column(., Site_ID=empirical_id,Method="NeuralNet",Type="Obs"))
      
    }
    
    d_param_infer=rbind(d_param_infer,as_tibble(t(colMeans(abc_sim$unadj.values)))%>%
                          add_column(., Site_ID=empirical_id,Method="NeuralNet"))
    
    
    
  }
  
  write.table(d_param_infer,paste0("../Data/Step11_Inferrence/Inferred_param_sites_",model,".csv"),sep=";")
  write.table(x_y_stat,paste0("../Data/Step11_Inferrence/x_y_obs_sim_stat_",model,".csv"),sep=";")
  
  
  pdf(paste0("../Figures/ABC_scale/obs_sim_all_models/obs_sim_sumstat_empirical_",model,".pdf"),width = 7,height = 7)
  par(mfrow=c(3,3))
  for (i in 1:9){
    plot(x=filter(x_y_stat,Type=="Sim",Method=="Rejection")[,i],xlab="Sim",ylab="Obs",main=colnames(x_y_stat)[i],
         y=filter(x_y_stat,Type=="Obs",Method=="Rejection")[,i],col="gray")
    abline(a=0,b=1)
  }
  dev.off()
  
}




x_y_stat=read.table("../Data/Step11_Inferrence/x_y_obs_sim_stat_Eby_feedback.csv",sep=";")

pdf("../Figures/ABC_scale/obs_sim_sumstat_empirical_all.pdf",width = 7,height = 7)
par(mfrow=c(3,3))
for (i in 1:9){
      plot(x=filter(x_y_stat,Type=="Sim",Method=="Rejection")[,i],xlab="Sim",ylab="Obs",main=colnames(x_y_stat)[i],
           y=filter(x_y_stat,Type=="Obs",Method=="Rejection")[,i],col="gray")
      abline(a=0,b=1)
}
dev.off()

x_y_stat=read.table("../Data/Step11_Inferrence/x_y_obs_sim_stat_Eby_feedback.csv",sep=";")
d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")%>%
  mutate(., Spectral_ratio=log(Spectral_ratio),clustering=log(clustering))

pdf("../Figures/ABC_scale/obs_sim_sumstat_empirical_low_res.pdf",width = 7,height = 7)
par(mfrow=c(3,3))
for (i in 1:9){
  d_sim=filter(x_y_stat,Type=="Sim",Method=="Rejection")%>%
    filter(.,d_biocom$Nbpixels<80000)
  d_obs=filter(x_y_stat,Type=="Obs",Method=="Rejection")%>%
    filter(.,d_biocom$Nbpixels<80000)
  plot(x=d_sim[,i],xlab="Sim",ylab="Obs",main=colnames(x_y_stat)[i],
       y=d_obs[,i],col="gray")
  abline(a=0,b=1)
}
dev.off()









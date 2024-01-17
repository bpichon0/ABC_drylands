rm(list=ls())
source("./ABC_drylands_function.R")

# ---------------------- Step 0: Merging simulations ----
"
Once simulations are made (using the Sim_ABC_main.jl file),
please run this script to merge all simulations into a single dataframe.
"
d_simu=tibble()
n_param=2

list_sim=list.files("./Data/Simulations/",".csv")

d=tibble()
for (sim in list_sim){
  d=rbind(d,read.table(paste0("./Data/Simulations/","/",sim),sep=",")[,(1):(n_param+14)])%>%
    filter(., V3>0.03)
}

colnames(d)=c("p","q","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
              "Spectral_ratio","PLR","PL_expo","cv_psd","median_psd","mean_psd","sd_psd","fmax_psd")


d=add_column(d,Pooling=rep(1:5,nrow(d)/5))
for (i in 1:(nrow(d)/5)){
  d$PL_expo[(5*(i-1)+2):(5*(i-1)+5)]=d$PL_expo[(5*(i-1)+1)]
  d$PLR[(5*(i-1)+2):(5*(i-1)+5)]=d$PLR[(5*(i-1)+1)]
  d$median_psd[(5*(i-1)+2):(5*(i-1)+5)]=d$median_psd[(5*(i-1)+1)]*d$Pooling[(5*(i-1)+2):(5*(i-1)+5)]
  d$mean_psd[(5*(i-1)+2):(5*(i-1)+5)]=d$mean_psd[(5*(i-1)+1)]*d$Pooling[(5*(i-1)+2):(5*(i-1)+5)]
  d$sd_psd[(5*(i-1)+2):(5*(i-1)+5)]=d$sd_psd[(5*(i-1)+1)]*d$Pooling[(5*(i-1)+2):(5*(i-1)+5)]
  d$cv_psd[(5*(i-1)+2):(5*(i-1)+5)]=d$cv_psd[(5*(i-1)+1)]
  d$fmax_psd[(5*(i-1)+2):(5*(i-1)+5)]=d$fmax_psd[(5*(i-1)+1)]
}

write.table(d,"./Data/Simulations.csv",sep=";")

# ---------------------- Step 1: Inference, running ABC ----
## >> 1.1) Inference parameters ----


Running_ABC=function(id,n_sim_kept=100){
  method_abc="rejection"
  `%!in%` = Negate(`%in%`)
  n_param=3
  
  d_biocom=read.table("./Data/data_sites.csv",sep=";")
  d_sim=read.table("./Data/Simulations.csv",sep=";",header=T)%>%
    dplyr::relocate(., Pooling,.after =q )%>%
    filter(., PL_expo>0)
  
  rownames(d_sim)=1:nrow(d_sim)
  
  
  d_param_infer_NN=array(0,c(n_sim_kept,nrow(d_biocom),3))
  d_param_infer_rej=array(0,c(n_sim_kept,nrow(d_biocom),3))
  
  d_NRMSE_sumstat=x_y_stat=tibble()
  
  
  sumstat_kept=4:14
  
  
  
  for (empirical_id in 1:345){
    
    target=d_biocom[empirical_id,10+sumstat_kept]
    matrix_param=d_sim[,1:3]
    
    mat_sumstat=rbind(d_sim[,sumstat_kept],target)
    
    #1) Boxcox
    
    for (x in 1:ncol(mat_sumstat)){
      if (any(is.na(target)) & x %!in% which(is.na(target))){
        if (colnames(mat_sumstat)[x] %in% c("skewness","moran_I","fmax_psd")){
          
          
          
          b=boxcox(lm(mat_sumstat[,x]+abs(min(mat_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat[,x] = (exp(mat_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
          }
          
        }else {
          b=boxcox(lm(mat_sumstat[,x]+.5 ~ 1),plotit = F,eps = .05)
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat[,x] = (mat_sumstat[,x]^(lambda_x) -1)/(lambda_x)
          }
        }
      }
    }
    
    
    #2) Scaling
    
    for (x in 1:ncol(mat_sumstat)) mat_sumstat[,x] = (mat_sumstat[,x]-mean(mat_sumstat[,x],na.rm = T))/sd(mat_sumstat[,x],na.rm = T)
    
    if (any(is.na(mat_sumstat[nrow(mat_sumstat),]))){
      
      which_na=which(is.na(mat_sumstat[nrow(mat_sumstat),]))
      
      cross_valid=abc(target = mat_sumstat[nrow(mat_sumstat),-which_na],
                      param = matrix_param[-nrow(mat_sumstat),],sumstat = mat_sumstat[-nrow(mat_sumstat),-which_na], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
      
    }else {
      cross_valid=abc(target = mat_sumstat[nrow(mat_sumstat),],
                      param = matrix_param[-nrow(mat_sumstat),],sumstat = mat_sumstat[-nrow(mat_sumstat),], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
    }
    
    
    #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
    
    mat_sumstat_step1=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
    mat_sumstat_step1=rbind(mat_sumstat_step1,target)
    
    #again, first box cox
    which_na=which(is.na(mat_sumstat[nrow(mat_sumstat),]))
    
    for (x in 1:ncol(mat_sumstat_step1)){
      
      if (any(which_na) & x %!in% which_na){
        
        if (colnames(mat_sumstat_step1)[x] %in% c("skewness","moran_I","fmax_psd")){
          
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
      }
    }
    
    #and normalization
    for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
    
    if (any(is.na(mat_sumstat_step1[nrow(mat_sumstat_step1),]))){
      
      which_na=which(is.na(mat_sumstat_step1[nrow(mat_sumstat_step1),]))
      
      cross_valid=abc(target = mat_sumstat_step1[nrow(mat_sumstat_step1),-which_na],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_step1[-nrow(mat_sumstat_step1),-which_na], #removing the target data
                      tol = n_sim_kept/nrow(mat_sumstat_step1),method = method_abc,transf = c(rep("logit",2),"none"), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),3,2,byrow = T),
                      numnet = 10,sizenet = 15)
      
      x_y_stat=rbind(x_y_stat,target[-which_na]%>%
                       add_column(., PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Obs")%>%
                       relocate(., PL_expo,.after =PLR ))
      
      x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)))%>%
                       add_column(.,PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Sim")%>%
                       relocate(., PL_expo,.after =PLR ))
      
      
    }else {
      cross_valid=abc(target = mat_sumstat_step1[nrow(mat_sumstat_step1),],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_step1[-nrow(mat_sumstat_step1),], #removing the target data
                      tol = n_sim_kept/nrow(mat_sumstat_step1),method = method_abc,transf = c(rep("logit",2),"none"), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),3,2,byrow = T),
                      numnet = 10,sizenet = 15)
      
      x_y_stat=rbind(x_y_stat,target%>%
                       add_column(., Site_ID=empirical_id,Method="rejection",Type="Obs"))
      
      x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)))%>%
                       add_column(.,Site_ID=empirical_id,Method="rejection",Type="Sim"))
      
    }
    
    
    cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
    
    mat_sumstat=d_sim[,sumstat_kept]
    
    if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
    
    cross_valid$adj.values=cross_valid$adj.values
    
    #NRMSE for the summary statistics observed
    RMSE = sapply(1:ncol(cross_valid$ss),function(x){
      sqrt(sum((cross_valid$ss[,x]-target[,x])**2,na.rm = T)/nrow(cross_valid$ss) )
    }
    )
    
    RMSE_prior=sapply(1:ncol(mat_sumstat),function(x){
      sqrt(sum((mat_sumstat[,x]-target[,x])**2,na.rm = T)/nrow(mat_sumstat) )
    }
    )
    NRMSE = RMSE/RMSE_prior
    
    d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,as_tibble(t(NRMSE)))
    
    d_param_infer_rej[,empirical_id,1]=cross_valid$unadj.values[,1] # we keep the whole distribution for p
    d_param_infer_rej[,empirical_id,2]=cross_valid$unadj.values[,2] # for q
    d_param_infer_rej[,empirical_id,3]=cross_valid$unadj.values[,3] # for the scale of observation
  }
  
  write.table(d_NRMSE_sumstat,paste0("./Data/NRMSE_sumstat.csv"),sep=";")
  write.table(x_y_stat,paste0("./Data/x_y_stat.csv"),sep=";")
  write.table(d_param_infer_rej,paste0("./Data/posterior_param.csv"),sep=";")

}

library(parallel)
mclapply(1:8,Running_ABC,mc.cores = 8)

## >> 1.2) Selecting relevant empirical data ----

#Testing for bimodality in the posterior distribution of parameters
library(diptest)
post_param=read.table("./Data/posterior_param.csv",sep=";")

d_biocom$bimod=sapply(1:nrow(d_biocom),function(x){
  if (dip.test(post_param[,x])$p.value<.05 | dip.test(post_param[,x+345])$p.value<.05){
    return("bimod")
  }else {
    return("unimod")
  }
})
write.table(which(d_biocom$bimod!="bimod"),"./Data/Keeping_sites.csv",sep=";")


## >> 1.3) Inference distance to the tipping point: merging bifurcation diagrams ----

# First Step 2 of Sim_ABC_main.jl called "Computing the distance to a tipping point"

d=tibble();step_size=0.005

for (site in list.files("./Data/Prediction/","Dist")){
  
  site_id=as.numeric(gsub(".csv","",strsplit(site,"_")[[1]][3]))
  
  pred=read.table(paste0("./Data/Prediction/",site),sep=",")%>%
    filter(., V1 != 0)
  colnames(pred)=c("p","q","cover")
  
  index=0;pred$ID_sim=NA
  for (x in 1:nrow(pred)){
    if (pred$p[x]==0.005){
      index=index+1
    }
    pred$ID_sim[x]=index
  }
  p_desert=sapply(unique(pred$ID_sim),function(x){
    d_fil=filter(pred,ID_sim==x,cover>0)
    if (any(d_fil$cover>0)){
      return(d_fil$p[1])
    }else {
      return(NA)
    }
  }) 
  
  p_infer=sapply(unique(pred$ID_sim),function(x){
    d_fil=filter(pred,ID_sim==x)
    return(d_fil$p[nrow(d_fil)])
  }) 
  
  q_infer=sapply(unique(pred$ID_sim),function(x){
    d_fil=filter(pred,ID_sim==x)
    return(d_fil$q[nrow(d_fil)])
  }) 
  
  size_tipping=sapply(unique(pred$ID_sim),function(x){
    d_fil=filter(pred,ID_sim==x)
    if (any(d_fil$cover>0)){
      return(d_fil$cover[1])
    }else {
      return(NA)
    }
    
  }) 
  
  d=rbind(d,tibble(ID_sim=1:length(p_desert),pcrit=p_desert,pinfer=p_infer,qinfer=q_infer,Size_tipping=size_tipping,
                   Site=site_id,
                   aridity=d_biocom$Aridity[site_id],Sand=d_biocom$Sand[site_id],
                   MF=d_biocom$MF[site_id]))
  
  #displaying the distribution
  
  d2=tibble(abs_dist=p_infer-p_desert,relativ_dist=(p_infer-p_desert)/(p_desert),Size_tipping=size_tipping)
  print(site)
  
  
}

write.table(d,"./Data/Resilience_metrics_1_neigh.csv",sep=";")


# ---------------------- Step 2: Methodology around inference ----
## >> 2.1) Optimizing pre and post-processing methods ----
#we play on both the lambda of the boxcox method and the size of the sample for
#stage 1 of the two step pre-processing procedure

dir.create("./Data/NRMSE",showWarnings = F)

d_all=read.table("./Data/Simulations.csv",sep=";",header=T)%>%
  filter(., Pooling==1)%>%
  dplyr::select(., -Pooling)%>%
  filter(., !is.na(PLR), !is.na(PL_expo)) #to avoid problems. this correspond to very high cover landscapes

rownames(d_all)=1:nrow(d_all)

N_for_cross_validation = 100
set.seed(123)
nrow_for_sample=sample(c(1:nrow(d_all)),N_for_cross_validation,replace = F)

for (optim_lambda in c(T,F)){
  
  for (size_step1 in c(1000,3000)){
    
    for (method_abc in c("loclinear","neuralnet")){
      
      for (preprocessing in c("BoxCox","None")){#c("PLS_BoxCox","BoxCox","None")){
        
        mat_cor_param=array(0,c(2,2,N_for_cross_validation)) #correlation matrix for parameters
        

        d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
        
        for (n in 1:N_for_cross_validation){
          
          matrix_param=d_all[,1:2]
          matrix_sumstat=d_all[,3:(ncol(d_all))]
          save_sumstat=matrix_sumstat
          
          n_cross=nrow_for_sample[n]
          
          if (preprocessing %in% c("BoxCox", "PLS_BoxCox")){ #Applying the two step procedure used in Siren MEE paper : Don't know whether it make sense in our case. TO discuss Monday
            
            if (optim_lambda==T){
              for (x in 1:ncol(matrix_sumstat)) if (colnames(matrix_sumstat)[x] %in% c("skewness","moran_I","fmax_psd")){
                
                b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
                lambda_x=b$x[which.max(b$y)]
                if (lambda_x !=0){ #to avoid errors
                  matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
                }
                
              }else {
                b=boxcox(lm(matrix_sumstat[,x]+.5 ~ 1),plotit = F,eps = .05)    
                lambda_x=b$x[which.max(b$y)]
                if (lambda_x !=0){ #to avoid errors
                  matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
                }
              }
              
            } else {
              for (x in 1:ncol(matrix_sumstat)) if (colnames(matrix_sumstat)[x] %in% c("skewness","moran_I","fmax_psd")){
                matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(.5)) -1)/(.5)
              }else {matrix_sumstat[,x] = (matrix_sumstat[,x]^(.5) -1)/(.5)}
            }
            
            #Second we scale
            for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
            
            if (preprocessing=="PLS_BoxCox"){
              #and finally, we perform the first PLS
              
              pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo+cv_psd+fmax_psd,
                         data=cbind(matrix_param,matrix_sumstat), scale=TRUE, validation="CV")
              
              
              n_comp_pls=selectNcomp(pls_1,method = "onesigma")
              
              if (n_comp_pls > 1){
                mat_sumstat_pls=pls_1$scores[,1:n_comp_pls] # selecting # components
              } else if (n_comp_pls==1){ #otherwise we take the whole components
                mat_sumstat_pls=as.data.frame(matrix(pls_1$scores[,1:n_comp_pls],ncol=1))
              } else {mat_sumstat_pls=pls_1$scores[,1:ncol(pls_1$scores)]}
              
              
            } else {mat_sumstat_pls=matrix_sumstat}
            
            if (any(is.na(mat_sumstat_pls[n_cross,]))){
              
              which_na=which(is.na(mat_sumstat_pls[n_cross,]))
              
              cross_valid1=abc(target = mat_sumstat_pls[n_cross,-which_na],
                              param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,-which_na], #removing the target data
                              tol = size_step1/(nrow(mat_sumstat_pls)),method = "rejection") #we keep the 1000 closest simulations for the first step
              
            }else {
              cross_valid1=abc(target = mat_sumstat_pls[n_cross,],
                              param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                              tol = size_step1/(nrow(mat_sumstat_pls)),method = "rejection") #we keep the 1000 closest simulations for the first step
            }
            
            if (nrow(cross_valid1$unadj.values) > size_step1){
              row_addi=sample(1:nrow(cross_valid1$unadj.values),abs(size_step1-nrow(cross_valid1$unadj.values)))
              cross_valid1$unadj.values=cross_valid1$unadj.values[-row_addi,]
              cross_valid1$ss=cross_valid1$unadj.values[-row_addi,]
            }
            
            #Keeping size_step1 simulations and doing the same steps again: normality, scaling and PLS
            
            mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid1$ss)),3:(ncol(d_all))] #we keep information with the true values
            mat_sumstat_step1=rbind(mat_sumstat_step1,d_all[n_cross,3:(ncol(d_all))])
            
            #again, first box cox
            
            
            if (optim_lambda==T){
              for (x in 1:ncol(mat_sumstat_step1)) if (colnames(mat_sumstat_step1)[x] %in% c("skewness","moran_I","fmax_psd")){
                
                b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
                lambda_x=b$x[which.max(b$y)]
                
                if (lambda_x !=0){ #to avoid errors
                  mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
                }
                
                
              }else {
                b=boxcox(lm(mat_sumstat_step1[,x]+.5 ~ 1),plotit = F,eps = .05)    
                lambda_x=b$x[which.max(b$y)]
                if (lambda_x !=0){ #to avoid errors
                  mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
                }
                
              }
              
            } else {
              for (x in 1:ncol(mat_sumstat_step1)) if (colnames(mat_sumstat_step1)[x] %in% c("skewness","moran_I","fmax_psd")){
                mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(.5)) -1)/(.5)
              }else {mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(.5) -1)/(.5)}
            }
            
            #and normalization
            for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
            
            
            if (preprocessing=="PLS_BoxCox"){
              
              pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo+cv_psd+fmax_psd,
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
            
            
            if (any(is.na(mat_sumstat_pls2[nrow(mat_sumstat_pls2),]))){
              
              which_na=which(is.na(mat_sumstat_pls2[nrow(mat_sumstat_pls2),]))
              
              cross_valid2=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),-which_na],
                               param = cross_valid1$unadj.values,
                               sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),-which_na], #removing the target data
                               tol = 100/(nrow(mat_sumstat_pls2)),method = method_abc,transf = rep("logit",2), #as parameters are proba, we perform logit regression
                               logit.bounds = matrix(c(0,1),2,2,byrow = T)) 
              
            }else {
              cross_valid2=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),],
                               param = cross_valid1$unadj.values,
                               sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),], #removing the target data
                               tol = 100/(nrow(mat_sumstat_pls2)),method = method_abc,transf = rep("logit",2), #as parameters are proba, we perform logit regression
                               logit.bounds = matrix(c(0,1),2,2,byrow = T)) 
            }
            
            
            
            cross_valid2$ss=d_all[as.numeric(rownames(cross_valid2$ss)),3:(ncol(d_all))] #we keep information with the true values
            
            
          } else { #no pls and boxcox
            
            if (any(is.na(matrix_sumstat[n_cross,]))){
              
              which_na=which(is.na(matrix_sumstat[n_cross,]))
              
              cross_valid2=abc(target = matrix_sumstat[n_cross,-which_na],
                               param = matrix_param[-n_cross,],
                               sumstat = matrix_sumstat[-n_cross,-which_na], #removing the target data
                               tol = 100/(nrow(matrix_param)),method = method_abc,transf = rep("logit",2), #as parameters are proba, we perform logit regression
                               logit.bounds = matrix(c(0,1),2,2,byrow = T)) 
              
            }else {
              cross_valid2=abc(target = matrix_sumstat[n_cross,],
                               param = matrix_param[-n_cross,],
                               sumstat = matrix_sumstat[-n_cross,], #removing the target data
                               tol = 100/(nrow(matrix_param)),method = method_abc,transf = rep("logit",2), #as parameters are proba, we perform logit regression
                               logit.bounds = matrix(c(0,1),2,2,byrow = T)) 
            }
            
          }    
          
          
          matrix_sumstat=save_sumstat
          
          
          if (names(cross_valid2)[1]=="unadj.values")names(cross_valid2)[1] = "adj.values"
          
          cross_valid2$adj.values=cross_valid2$adj.values
          #Matrix of correlation between parameters & sumstats
          mat_cor_param[,,n]=cor(cross_valid2$adj.values)
          
          
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
        
        write.table(d_NRMSE_param,paste0("./Data/NRMSE/RMSE_param_",preprocessing,
                                         "_",method_abc,"_optim_lambda_",ifelse(optim_lambda,"yes","no"),"_N1_",size_step1,".csv"),sep=";")
        
        colnames(d_NRMSE_sumstat)=colnames(d_cross_sumstat)
        
        write.table(d_NRMSE_sumstat,paste0("./Data/NRMSE/RMSE_sumstat_",preprocessing,
                                           "_",method_abc,"_optim_lambda_",ifelse(optim_lambda,"yes","no"),"_N1_",size_step1,".csv"),sep=";")
        
        

        
        
        
      }
    }
  }
}





## >> 2.2) Optimizing the structure of the neural-network ----

dir.create("./Data/NRMSE",showWarnings = F)


d_all=read.table("./Data/Simulations.csv",sep=";",header = T)%>%
  filter(., Pooling==1)%>%
  dplyr::select(., -Pooling)%>%
  filter(., !is.na(PLR), !is.na(PL_expo)) #to avoid problems. this correspond to very high cover landscapes

rownames(d_all)=1:nrow(d_all)


N_for_cross_validation = 100
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)

for (method_pre in c("NoPLS")){#c("PLS","NoPLS")){
  
  for (size_hidden in seq(5,25,by=5)){
    
    for (rep_network in c(10,20,30)){
      
      mat_cor_param=array(0,c(2,2,N_for_cross_validation)) #correlation matrix for parameters
      
      
      d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
      
      for (n in 1:N_for_cross_validation){
        
        matrix_param=d_all[,1:2]
        matrix_sumstat=d_all[,3:(ncol(d_all))]
        save_sumstat=matrix_sumstat
        
        n_cross=nrow_for_sample[n]
        
        for (x in 1:ncol(matrix_sumstat)) if (colnames(matrix_sumstat)[x] %in% c("skewness","moran_I","fmax_psd")){
          
          b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
          }
          
        }else {
          b=boxcox(lm(matrix_sumstat[,x]+.5 ~ 1),plotit = F,eps = .05)
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
          }
        }
        
        
        #Second we scale
        for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
        
        #and finally, we perform the first PLS
        if (method_pre=="PLS"){
          
          pls_1=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo+cv_psd+fmax_psd,
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
        
        if (any(is.na(mat_sumstat_pls[n_cross,]))){
          
          which_na=which(is.na(mat_sumstat_pls[n_cross,]))
          
          cross_valid=abc(target = mat_sumstat_pls[n_cross,-which_na],
                          param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,-which_na], #removing the target data
                          tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
          
        }else {
          cross_valid=abc(target = mat_sumstat_pls[n_cross,],
                          param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                          tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
        }
        
        
        #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
        
        mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid$ss)),3:(ncol(d_all))] #we keep information with the true values
        mat_sumstat_step1=rbind(mat_sumstat_step1,d_all[n_cross,3:(ncol(d_all))])
        
        #again, first box cox
        
        for (x in 1:ncol(mat_sumstat_step1)) if (colnames(mat_sumstat_step1)[x] %in% c("skewness","moran_I","fmax_psd")){
          
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
          
          pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo+cv_psd+fmax_psd,
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
        
        if (any(is.na(mat_sumstat_pls2[nrow(mat_sumstat_pls2),]))){
          
          which_na=which(is.na(mat_sumstat_pls2[nrow(mat_sumstat_pls2),]))
          
          cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),-which_na],
                          param = cross_valid$unadj.values,
                          sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),-which_na], #removing the target data
                          tol = 100/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                          logit.bounds = matrix(c(0,1),2,2,byrow = T),
                          numnet = rep_network,sizenet = size_hidden)
          
        }else {
          cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),],
                          param = cross_valid$unadj.values,
                          sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),], #removing the target data
                          tol = 100/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                          logit.bounds = matrix(c(0,1),2,2,byrow = T),
                          numnet = rep_network,sizenet = size_hidden)
        }
        
        
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
      
      write.table(d_NRMSE_param,paste0("./Data/NRMSE/RMSE_hidden_preprocessing_",method_pre,"_",
                                       size_hidden,"_Nnet_",rep_network,".csv"),sep=";")
      
      
      
    }
  }
}





## >> 2.3) Influence of the number of simulation kept ----

dir.create("./Data/NRMSE",showWarnings = F)
d_all=read.table("./Data/Simulations.csv",sep=";",header = T)%>%
  filter(., Pooling==1)%>%
  dplyr::select(., -Pooling)%>%
  filter(., !is.na(PLR), !is.na(PL_expo)) #to avoid problems. this correspond to very high cover landscapes

rownames(d_all)=1:nrow(d_all)

N_for_cross_validation = 100
set.seed(123)
nrow_for_sample=sample(c(1:nrow(d_all)),N_for_cross_validation,replace = F)


for (NA_kept in c(50,100,150,200,250)){
  
  
  d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
  
  for (n in 1:N_for_cross_validation){
    
    matrix_param=d_all[,1:2]
    matrix_sumstat=d_all[,3:(ncol(d_all))]
    save_sumstat=matrix_sumstat
    
    n_cross=nrow_for_sample[n]
    
    for (x in 1:ncol(matrix_sumstat)) if (colnames(matrix_sumstat)[x] %in% c("skewness","moran_I","fmax_psd")){
      
      b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
      }
      
    }else {
      b=boxcox(lm(matrix_sumstat[,x]+.5 ~ 1),plotit = F,eps = .05)
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
      }
    }
    
    
    
    for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
    
    
    mat_sumstat_pls=matrix_sumstat
    
    
    if (any(is.na(mat_sumstat_pls[n_cross,]))){
      
      which_na=which(is.na(mat_sumstat_pls[n_cross,]))
      
      cross_valid=abc(target = mat_sumstat_pls[n_cross,-which_na],
                      param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,-which_na], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
      
    }else {
      cross_valid=abc(target = mat_sumstat_pls[n_cross,],
                      param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
    }
    
    
    
    mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid$ss)),3:(ncol(d_all))] #we keep information with the true values
    mat_sumstat_step1=rbind(mat_sumstat_step1,d_all[n_cross,3:(ncol(d_all))])
    
    
    for (x in 1:ncol(mat_sumstat_step1)) if (colnames(mat_sumstat_step1)[x] %in% c("skewness","moran_I","fmax_psd")){
      
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
    
    for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
    
    
    
    if (method_pre=="PLS"){
      
      pls_2=plsr(p + q~rho_p+nb_neigh+clustering+skewness+variance+moran_I+Spectral_ratio+PLR+PL_expo+cv_psd+fmax_psd,
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
    
    if (any(is.na(mat_sumstat_pls2[nrow(mat_sumstat_pls2),]))){
      
      which_na=which(is.na(mat_sumstat_pls2[nrow(mat_sumstat_pls2),]))
      
      cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),-which_na],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),-which_na], #removing the target data
                      tol = NA_kept/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),2,2,byrow = T),
                      numnet = 10,sizenet = 10)
      
    }else {
      cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),], #removing the target data
                      tol = NA_kept/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),2,2,byrow = T),
                      numnet = 10,sizenet = 10)
    }
    
    
    cross_valid$ss=d_all[as.numeric(rownames(cross_valid$ss)),3:(ncol(d_all))] #we keep information with the true values
    
    matrix_sumstat=save_sumstat
    
    if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
    
    cross_valid$adj.values=cross_valid$adj.values
    
    
    
    if (any( is.nan(cross_valid$adj.values[,1]) | is.nan(cross_valid$adj.values[,2]))){
      cross_valid$adj.values = cross_valid$adj.values[-which(is.nan(cross_valid$adj.values[,1])
                                                             | is.nan(cross_valid$adj.values[,2])),]
      
    }
    
    
    d_cross_param=rbind(d_cross_param,as_tibble(t(colMeans(cross_valid$adj.values)))%>%add_column(., Type="Sim"))
    d_cross_param=rbind(d_cross_param,as_tibble((matrix_param[n_cross,]))%>%add_column(., Type="Obs"))
    
    d_cross_sumstat=rbind(d_cross_sumstat,as_tibble(t(colMeans(cross_valid$ss)))%>%add_column(., Type="Sim"))
    d_cross_sumstat=rbind(d_cross_sumstat,as_tibble((matrix_sumstat[n_cross,]))%>%add_column(., Type="Obs"))
    
    
    RMSE = sapply(1:ncol(cross_valid$adj.values),function(x){
      sqrt(sum((cross_valid$adj.values[,x]-matrix_param[n_cross,x])**2)/nrow(cross_valid$adj.values) )
    }
    )
    
    RMSE_prior=sapply(1:(ncol(matrix_param)),function(x){
      sqrt(sum((matrix_param[,x]-matrix_param[n_cross,x])**2)/nrow(matrix_param) )
    }
    )
    NRMSE = RMSE/RMSE_prior
    
    d_NRMSE_param=rbind(d_NRMSE_param,as_tibble(t(NRMSE)))
    
    
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
  
  write.table(d_NRMSE_param,paste0("./Data/NRMSE/RMSE_NAkept_",NA_kept,".csv"),sep=";")
  
}






## >> 2.4) Selecting the best summary statistics ----

dir.create("./Data/Best_sumstat",showWarnings = F)

d_all=read.table("./Data/Simulations.csv",sep=";",header=T)%>%
  filter(., Pooling==1)%>%
  dplyr::select(., -Pooling)

rownames(d_all)=1:nrow(d_all)

N_for_cross_validation = 100
set.seed(123)
nrow_for_sample=sample(c(1:nrow(matrix_param)),N_for_cross_validation,replace = F)

size_hidden=10
rep_network=10
all_name=c("All","no_PLR","no_PL_expo","no_PLR_PL","no_cv","no_fmax","no_fmax_cv","no_PL_cv","no_cv_PL_PLR")
list_sumstat=list(c(1:ncol(matrix_sumstat)),
                  c(1:ncol(matrix_sumstat))[-c(8)],
                  c(1:ncol(matrix_sumstat))[-c(9)],
                  c(1:ncol(matrix_sumstat))[-c(8:9)],
                  c(1:ncol(matrix_sumstat))[-c(10)],
                  c(1:ncol(matrix_sumstat))[-c(11)],
                  c(1:ncol(matrix_sumstat))[-c(10:11)],
                  c(1:ncol(matrix_sumstat))[-c(9:10)],
                  c(1:ncol(matrix_sumstat))[-c(8:10)])


for (which_sumstat in 1:length(all_name)){
  
  name_removal=all_name[which_sumstat]
  
  summary_stat_kept=list_sumstat[[which_sumstat]]
  
  
  mat_cor_param=array(0,c(2,2,N_for_cross_validation)) #correlation matrix for parameters
  
  
  d_cross_param=d_cross_sumstat=d_NRMSE_param=d_NRMSE_sumstat=tibble()
  
  for (n in 1:N_for_cross_validation){
    
    matrix_param=d_all[,1:2]
    matrix_sumstat=d_all[,2+summary_stat_kept]
    save_sumstat=matrix_sumstat
    
    n_cross=nrow_for_sample[n]
    
    for (x in 1:ncol(matrix_sumstat)) if (colnames(matrix_sumstat)[x] %in% c("skewness","moran_I","fmax_psd")){
      
      b=boxcox(lm(matrix_sumstat[,x]+abs(min(matrix_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values, does not change the distribution
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        matrix_sumstat[,x] = (exp(matrix_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
      }
      
    }else {
      b=boxcox(lm(matrix_sumstat[,x]+.5 ~ 1),plotit = F,eps = .05)
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        matrix_sumstat[,x] = (matrix_sumstat[,x]^(lambda_x) -1)/(lambda_x)
      }
    }
    
    
    #Second we scale
    for (x in 1:ncol(matrix_sumstat)) matrix_sumstat[,x] = (matrix_sumstat[,x]-mean(matrix_sumstat[,x],na.rm = T))/sd(matrix_sumstat[,x],na.rm = T)
    
    mat_sumstat_pls=matrix_sumstat
    
    if (any(is.na(mat_sumstat_pls[n_cross,]))){
      
      which_na=which(is.na(mat_sumstat_pls[n_cross,]))
      
      cross_valid=abc(target = mat_sumstat_pls[n_cross,-which_na],
                      param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,-which_na], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
      
    }else {
      cross_valid=abc(target = mat_sumstat_pls[n_cross,],
                      param = matrix_param[-n_cross,],sumstat = mat_sumstat_pls[-n_cross,], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
    }
    
    
    
    #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
    
    mat_sumstat_step1=d_all[as.numeric(rownames(cross_valid$ss)),3:(ncol(d_all))] #we keep information with the true values
    mat_sumstat_step1=rbind(mat_sumstat_step1,d_all[n_cross,3:(ncol(d_all))])
    
    #again, first box cox
    
    for (x in 1:ncol(mat_sumstat_step1)) if (colnames(mat_sumstat_step1)[x] %in% c("skewness","moran_I","fmax_psd")){
      
      b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
      lambda_x=b$x[which.max(b$y)]
      
      if (lambda_x !=0){ #to avoid errors
        mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
      }
      
      
    }else {
      b=boxcox(lm(mat_sumstat_step1[,x]+.5 ~ 1),plotit = F,eps = .05)
      lambda_x=b$x[which.max(b$y)]
      if (lambda_x !=0){ #to avoid errors
        mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
      }
      
    }
    
    #and normalization
    for (x in 1:ncol(mat_sumstat_step1)) mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
    
    mat_sumstat_pls2=mat_sumstat_step1
    
    
    if (any(is.na(mat_sumstat_pls2[nrow(mat_sumstat_pls2),]))){
      
      which_na=which(is.na(mat_sumstat_pls2[nrow(mat_sumstat_pls2),]))
      
      cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),-which_na],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),-which_na], #removing the target data
                      tol = 75/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),2,2,byrow = T),
                      numnet = rep_network,sizenet = size_hidden)
      
    }else {
      cross_valid=abc(target = mat_sumstat_pls2[nrow(mat_sumstat_pls2),],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_pls2[-nrow(mat_sumstat_pls2),], #removing the target data
                      tol = 75/nrow(mat_sumstat_pls2),method = "neuralnet",transf = rep("logit",2), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),2,2,byrow = T),
                      numnet = rep_network,sizenet = size_hidden)
    }
    
    cross_valid$ss=d_all[as.numeric(rownames(cross_valid$ss)),2+summary_stat_kept] #we keep information with the true values
    
    matrix_sumstat=save_sumstat
    
    if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
    
    cross_valid$adj.values=cross_valid$adj.values
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
  
  write.table(d_NRMSE_param,paste0("./Data/Best_sumstat/RMSE_selecting_sumstat_,",name_removal,".csv"),sep=";")

}






## >> 2.5) Can we recover parameters and scale of observation in simulations ----

dir.create("./Data/Scale_obs_indentifiability",showWarnings = F)
d_sim=read.table("./Data/Simulations.csv",sep=";",header=T)%>%
  add_column(.,ID_sim=rep(1:(nrow(.)/5),each=5))%>%
  filter(., rho_p>0.03)%>%
  dplyr::relocate(.,Pooling,.after=q)

set.seed(123)
list_sample=sample(unique(d_sim$ID_sim),100) #for each of these 100 parameter sets, we need to verify whether we infer the same parameters of not
d_RMSE_param=x_y_param=tibble()
n_keep=100

stat_kept=c(4:13,ncol(d_sim)-1)

for (sample_id in list_sample){
  
  for (scale_obs in c(5,4,3,2,1)){ #we include the scale of observation as a parameter of the model
    
    target=d_sim%>%filter(., ID_sim==sample_id,Pooling==scale_obs)%>%
      dplyr::select(., -ID_sim)
    
    all_sim=d_sim[-which(d_sim$Pooling==scale_obs & d_sim$ID_sim==sample_id),]%>%
      dplyr::select(., -ID_sim)
    
    if (any(is.na(target))){
      
      which_na=which(is.na(target))
      abc_sim=abc(target=target[stat_kept[-(which_na-3)]],
                  param = all_sim[,c(1:3)],sumstat = all_sim[,stat_kept[-(which_na-3)]],
                  tol = n_keep/nrow(all_sim),method="rejection")
      
    }else {
      abc_sim=abc(target=target[,stat_kept],
                  param = all_sim[,1:3],sumstat = all_sim[,stat_kept],
                  tol = n_keep/nrow(all_sim),method="rejection")
    }
    
    RMSE = sapply(1:ncol(abc_sim$unadj.values),function(x){
      sqrt(sum((abc_sim$unadj.values[,x]-as.numeric(target[,x]))**2)/nrow(abc_sim$unadj.values) )})
    
    RMSE_prior=sapply(1:3,function(x){
      sqrt(sum(((all_sim[,x]-target[,x])**2)/nrow(all_sim) ))})
    
    NMRSE=RMSE/RMSE_prior
    names(NMRSE)=c("p","q","Pooling")
    
    d_RMSE_param=rbind(d_RMSE_param,data.frame((t(NMRSE)))%>%
                         add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection",Scale_obs=scale_obs))
    
    x_y_param=rbind(x_y_param,data.frame((t(colMeans(abc_sim$unadj.values))))%>%
                      add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection",Type="Sim",Scale_obs=scale_obs))
    
    x_y_param=rbind(x_y_param,as.data.frame(target[,1:3])%>%
                      add_column(., Site_ID=which(list_sample == sample_id),Method="Rejection",Type="Obs",Scale_obs=scale_obs))
    
  }
}

write.table(d_RMSE_param,"./Data/Scale_obs_indentifiability/Retrieving_parameters_different_resolution_RMSE_param.csv",sep=";")
write.table(x_y_param,"./Data/Scale_obs_indentifiability/Retrieving_parameters_different_resolution_x_y.csv",sep=";")


# ---------------------- Step 3: Running mixed effects models ----
## >> 4.1) Without facilitation ----

boot_function_lm = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}

d=read.table("./Data/posterior_param.csv",sep=";",header=T)
keep_sites=read.table("./Data/Keeping_sites_biocom.csv",sep=";")$V1
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

d2=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")%>%
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


d2=tibble(p=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_p[x],d$sd_p[x]))}))[,1],
          q=logit(d$median_q),
          Size_tipping=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$Size_mean[x],d$Size_sd[x]))}))[,1],
          abs_dist=scale(d$abs_median)[,1],
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
          Long=d$Long,
          Long_sin=(d$Long_sin-mean(d$Long_sin,na.rm=T))/sd(d$Long_sin,na.rm = T),
          Elevation=(d$Elevation-mean(d$Elevation,na.rm=T))/sd(d$Elevation,na.rm = T),
          Slope=(d$Slope-mean(d$Slope,na.rm=T))/sd(d$Slope,na.rm = T),
          Plot_n=d$Plot_n)

mod_predictors=gsub("\n     ","","Aridity + MF + Sand + Soil_A +
                    Lat + Long_cos + Long_sin + Slope + Elevation + ( 1 | Plot_n)")
d_mod=list(Boot_effects=tibble(),Partial_res_data=tibble())
d_info_model=list(global_R2=tibble(),R2_partial_res=tibble(),Vif=tibble(),Moran=tibble(),Effects=tibble())

for (response_var in c("Cover","q","abs_dist","rela_dist")){
  
  model_abs=(lmer(formula = paste(response_var," ~ ",mod_predictors), data = d2)) #fitting the model
  
  mcp.fnc(model_abs) #checking model assumptions
  
  #ARIDITY

  resid_mod=visreg::visreg(fit = model_abs,xvar="Aridity",plot=F) 
  boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Aridity)
  
  d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                    tibble(Driver_name="Aridity",
                                           Response_var=response_var,
                                           R2=rsq(lm(data=resid_mod$res,visregRes~Aridity))))
  
  d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                               tibble(Driver_value=resid_mod$res$Aridity,
                                      Driver_name="Aridity",
                                      Response_var=response_var,
                                      Resids=resid_mod$res$visregRes))
  
  d_mod$Boot_effects=rbind(d_mod$Boot_effects,tibble(Slopes=boot_AI$t[,1],Driver_name="Aridity",Response_var=response_var))
  
  
  # Multifunctionality 
  
  resid_mod=visreg::visreg(fit = model_abs,xvar="MF",plot=F) 
  boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~MF)
  
  d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                    tibble(Driver_name="Multifunctionality",
                                           Response_var=response_var,
                                           R2=rsq(lm(data=resid_mod$res,visregRes~MF))))
  
  d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                               tibble(Driver_value=resid_mod$res$MF,
                                      Driver_name="Multifunctionality",
                                      Response_var=response_var,
                                      Resids=resid_mod$res$visregRes))
  
  d_mod$Boot_effects=rbind(d_mod$Boot_effects,tibble(Slopes=boot_AI$t[,1],Driver_name="Multifunctionality",Response_var=response_var))
  
  
  # Sand soil content 
  
  resid_mod=visreg::visreg(fit = model_abs,xvar="Sand",plot=F) 
  boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Sand)
  
  d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                    tibble(Driver_name="Sand",
                                           Response_var=response_var,
                                           R2=rsq(lm(data=resid_mod$res,visregRes~Sand))))
  
  d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                               tibble(Driver_value=resid_mod$res$Sand,
                                      Driver_name="Sand",
                                      Response_var=response_var,
                                      Resids=resid_mod$res$visregRes))
  
  d_mod$Boot_effects=rbind(d_mod$Boot_effects,tibble(Slopes=boot_AI$t[,1],Driver_name="Sand",Response_var=response_var))
  

  
  # Soil amelioration 
  
  resid_mod=visreg::visreg(fit = model_abs,xvar="Soil_A",plot=F) 
  boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Soil_A)
  
  d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                    tibble(Driver_name="Soil amelioration",
                                           Response_var=response_var,
                                           R2=rsq(lm(data=resid_mod$res,visregRes~Soil_A))))
  
  d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                               tibble(Driver_value=resid_mod$res$Soil_A,
                                      Driver_name="Soil amelioration",
                                      Response_var=response_var,
                                      Resids=resid_mod$res$visregRes))
  
  d_mod$Boot_effects=rbind(d_mod$Boot_effects,tibble(Slopes=boot_AI$t[,1],Driver_name="Soil amelioration",Response_var=response_var))
  
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

saveRDS(d_mod,"./Data/Drivers_stability_metrics_data_uncertainty_without_facilitation.rds")



## >> 4.2) With facilitation ----

boot_function_lm = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}

d=read.table("./Data/posterior_param.csv",sep=";",header=T)
keep_sites=read.table("./Data/Keeping_sites_biocom.csv",sep=";")$V1
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
         Cover=d_biocom$Cover)%>%
  filter(., Site %in% keep_sites)

d2=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")%>%
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


d2=tibble(p=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_p[x],d$sd_p[x]))}))[,1],
          q=logit(d$median_q),
          Size_tipping=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$Size_mean[x],d$Size_sd[x]))}))[,1],
          abs_dist=scale(d$abs_median)[,1],
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

mod_predictors=gsub("\n     ","","Aridity + MF + Sand + Soil_A + Facilitation +
                    Lat + Long_cos + Long_sin + Slope + Elevation + ( 1 | Plot_n)")

model_abs=(lmer(formula = paste("abs_dist ~ ",mod_predictors), data = d2)) #fitting the model
model_q=(lmer(formula = paste("q ~ ",mod_predictors), data = d2)) #fitting the model
model_cover=(lmer(formula = paste("Cover ~ ",mod_predictors), data = d2)) #fitting the model

d_mod=list(Boot_effects=tibble(),Partial_res_data=tibble())

# Facilitation: distance

resid_mod=visreg::visreg(fit = model_abs,xvar="Facilitation",plot=F) 
boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Facilitation)

d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                  tibble(Driver_name="Facilitation",
                                         Response_var="abs_dist",
                                         R2=rsq(lm(data=resid_mod$res,visregRes~Facilitation))))

d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                             tibble(Driver_value=resid_mod$res$Facilitation,
                                    Driver_name="Facilitation",
                                    Response="Distance",
                                    Resids=resid_mod$res$visregRes))

# Facilitation: q

d_mod$Boot_effects=rbind(d_mod$Boot_effects,tibble(Slopes=boot_AI$t[,1],Driver_name="Facilitation",Response="Distance"))

resid_mod=visreg::visreg(fit = model_q,xvar="Facilitation",plot=F) 
boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Facilitation)

d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                  tibble(Driver_name="Facilitation",
                                         Response_var="q",
                                         R2=rsq(lm(data=resid_mod$res,visregRes~Facilitation))))

d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                             tibble(Driver_value=resid_mod$res$Facilitation,
                                    Driver_name="Facilitation",
                                    Response="Aggregation \n vegetation (q)",
                                    Resids=resid_mod$res$visregRes))

# Facilitation: cover

d_mod$Boot_effects=rbind(d_mod$Boot_effects,tibble(Slopes=boot_AI$t[,1],Driver_name="Facilitation",Response="Aggregation \n vegetation (q)"))

resid_mod=visreg::visreg(fit = model_cover,xvar="Facilitation",plot=F) 
boot_AI = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Facilitation)

d_info_model$R2_partial_res=rbind(d_info_model$R2_partial_res,
                                  tibble(Driver_name="Facilitation",
                                         Response_var="Cover",
                                         R2=rsq(lm(data=resid_mod$res,visregRes~Facilitation))))

d_mod$Partial_res_data=rbind(d_mod$Partial_res_data,
                             tibble(Driver_value=resid_mod$res$Facilitation,
                                    Driver_name="Facilitation",
                                    Response="Vegetation cover",
                                    Resids=resid_mod$res$visregRes))

d_mod$Boot_effects=rbind(d_mod$Boot_effects,tibble(Slopes=boot_AI$t[,1],Driver_name="Facilitation",Response="Vegetation cover"))

saveRDS(d_mod,"./Data/Drivers_stability_metrics_data_uncertainty_with_facilitation.rds")
saveRDS(d_info_model,"./Data/Properties_models.rds")




# ---------------------- Step 4: Bootstrap AIC q cover ----

#with data uncertainty
d=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")
post_param=read.table("./Data/posterior_param.csv",sep=";")
# summarizing information in each site
d=d%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",abs_dis50=quantile(pinfer-pcrit,na.rm = T,.5),
                   abs_dis25=quantile(pinfer-pcrit,na.rm = T,.25),
                   abs_dis75=quantile(pinfer-pcrit,na.rm = T,.75),
                   relativ_dis50=quantile((pinfer-pcrit)/pinfer,na.rm = T,.5),
                   relativ_dis25=quantile((pinfer-pcrit)/pinfer,na.rm = T,.25),
                   relativ_dis75=quantile((pinfer-pcrit)/pinfer,na.rm = T,.75),
                   Size_tipping50=quantile(Size_tipping,na.rm = T,.5),
                   Size_tipping25=quantile(Size_tipping,na.rm = T,.25),
                   Size_tipping75=quantile(Size_tipping,na.rm = T,.75),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
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

d2=tibble(p=logit(d$median_p),
          q=logit(d$median_q),
          abs_dist=scale(d$abs_dis50)[,1],
          rela_dist=scale(d$relativ_dis50)[,1],
          Site=d$Site,
          Cover=scale(d$Cover)[,1],
          Plot_n=d$Plot_n)

#Comparing cover + spatial structure with only cover to see whether spatial structure helps to indicate distance or
#size of tipping point

n_boot=500
d_AIC=tibble()
for (k in 1:n_boot){
  
  data_sampled=d2[sample(1:nrow(d2),replace = T),]
  
  #models for cover
  model_abs_cover=(data_sampled%>% lmer(formula = paste("abs_dist ~ Cover + ( 1 | Plot_n)"), data = .))
  model_rela_cover=(data_sampled%>% lmer(formula = paste("rela_dist ~ Cover + ( 1 | Plot_n)"), data = .))
  
  #models for q+cover
  model_abs_both=(data_sampled%>% lmer(formula = paste("abs_dist ~ q + Cover + ( 1 | Plot_n)"), data = .))
  model_rela_both=(data_sampled%>% lmer(formula = paste("rela_dist ~ q + Cover + ( 1 | Plot_n)"), data = .))
  
  #models for q
  model_abs_q=(data_sampled%>% lmer(formula = paste("abs_dist ~ q + ( 1 | Plot_n)"), data = .))
  model_rela_q=(data_sampled%>% lmer(formula = paste("rela_dist ~ q + ( 1 | Plot_n)"), data = .))
  
  AIC_abs=AIC(model_abs_cover,model_abs_both,model_abs_q)
  AIC_rela=AIC(model_rela_cover,model_rela_both,model_rela_q)
  d_AIC=rbind(d_AIC,tibble(AIC=c(AIC_abs$AIC,AIC_rela$AIC),
                           Stability=rep(c("Absolute distance","Relative distance"),each=3),
                           Predictor=rep(c("Cover","Cover + Spatial \n structure","Spatial \n structure"),2),
                           ID=k))
  
}

write.table(d_AIC,"./Data/Cover_vs_spatial_structure_data_uncertainty.csv",sep=";")



# ---------------------- Step 5: Similarity inference within site ----

keep_sites=read.table("./Data/Keeping_sites_biocom.csv",sep=";")$V1
d=read.table("./Data/posterior_param.csv",sep=";") #posterior
d2=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";") #predicted metrics
d_biocom=read.table("./Data/biocom_data.csv",sep=";")%>%
  filter(., c(1:nrow(.)) %in% keep_sites)

d_summarized=cbind(d2%>%
                     dplyr::group_by(., Site,MF,aridity,Sand)%>%
                     dplyr::summarise(., .groups = "keep",
                                      mean_rela=mean((pinfer-pcrit)/pcrit,na.rm = T),
                                      sd_rela=sd((pinfer-pcrit)/pcrit,na.rm = T),
                                      mean_abs=mean((pinfer-pcrit)/pcrit,na.rm = T),
                                      sd_abs=sd((pinfer-pcrit)/pcrit,na.rm = T),
                                      median_abs=median((pinfer-pcrit),na.rm = T),
                                      median_rela=median((pinfer-pcrit)/pcrit,na.rm = T),
                                      mean_size=mean((pinfer-pcrit)/pcrit,na.rm = T),
                                      sd_size=sd((pinfer-pcrit)/pcrit,na.rm = T))%>%
                     filter(., Site %in% keep_sites),
                   tibble(mean_p=colMeans(d)[keep_sites],
                          mean_q=colMeans(d)[keep_sites+345],
                          median_p=apply(d,2,median)[keep_sites],
                          median_q=apply(d,2,median)[keep_sites+345],
                          p_sd=apply(d,2,sd,na.rm=T)[keep_sites],
                          q_sd=apply(d,2,sd,na.rm=T)[keep_sites+345]))%>%
  add_column(., Plot_n=d_biocom$Plot_n[.$Site])

d_ABC=d_median=tibble()
for (k in unique(d_biocom$Plot_n)){ #each site
  
  
  #using ABC uncertainty
  for (rep_id in 1:100){
    
    
    #First parameters 
    
    
    p_site=sapply(which(d_biocom$Plot_n==k),function(x){
      return(rnorm(1,mean=d_summarized$mean_p[x],
                   sd = d_summarized$p_sd[x]))
    }) #get the values of p
    
    q_site=sapply(which(d_biocom$Plot_n==k),function(x){
      return(rnorm(1,mean=d_summarized$mean_q[x],
                   sd = d_summarized$q_sd[x]))
    }) #get the values of p
    
    p_other=sapply(which(d_biocom$Plot_n!=k),function(x){
      return(rnorm(1,mean=d_summarized$mean_p[x],sd=d_summarized$p_sd[x]))
    }) #get the values of p aside this site
    q_other=sapply(which(d_biocom$Plot_n!=k),function(x){
      return(rnorm(1,mean=d_summarized$mean_q[x],sd=d_summarized$q_sd[x]))
    }) #get the values of q aside this site
    
    RMSE_p_within = mean(sapply(1:length(p_site),function(x){
      return(sqrt(sum((p_site-p_site[x])**2,na.rm = T)/length(p_site) ))})) #NRMSE of p within
    
    RMSE_q_within = mean(sapply(1:length(q_site),function(x){
      return(sqrt(sum((q_site-q_site[x])**2,na.rm = T)/length(q_site) ))})) #NRMSE of q within
    
    RMSE_p_all = mean(sapply(1:length(p_site),function(x){
      return(sqrt(sum((p_other-p_site[x])**2,na.rm = T)/length(p_other) ))})) #NRMSE of p between
    RMSE_q_all = mean(sapply(1:length(q_site),function(x){
      return(sqrt(sum((q_other-q_site[x])**2,na.rm = T)/length(q_other) ))})) #NRMSE of q between
    
    NRMSE_p=RMSE_p_within/RMSE_p_all
    NRMSE_q=RMSE_q_within/RMSE_q_all
    
    
    
    #second for resilience related predictions
    
    #1) relative distance
    
    dist_site=sapply(which(d_biocom$Plot_n==k),function(x){
      return(rnorm(1,mean=d_summarized$mean_rela[x],
                   sd = d_summarized$sd_rela[x]))
    }) #get the values of distance to desert state within the site
    
    dist_other=sapply(which(d_biocom$Plot_n!=k),function(x){
      return(rnorm(1,mean=d_summarized$mean_rela[x],sd=d_summarized$sd_rela[x]))
    }) #get the values of distance aside this site
    
    if (any (dist_other<0)){
      dist_other=dist_other[-which(dist_other<0)] #can happen
    }
    
    RMSE_dist_within = mean(sapply(1:length(dist_site),function(x){
      return(sqrt(sum((dist_site-dist_site[x])**2,na.rm = T)/length(dist_site) ))})) #NRMSE of p within
    
    RMSE_dist_all = mean(sapply(1:length(dist_site),function(x){
      return(sqrt(sum((dist_other-dist_site[x])**2,na.rm = T)/length(dist_other) ))})) #NRMSE of p between
    
    NRMSE_rela=RMSE_dist_within/RMSE_dist_all
    
    #2) absolute distance
    
    dist_site=sapply(which(d_biocom$Plot_n==k),function(x){
      return(rnorm(1,mean=d_summarized$mean_abs[x],
                   sd = d_summarized$sd_abs[x]))
    }) #get the values of distance to desert state within the site
    
    dist_other=sapply(which(d_biocom$Plot_n!=k),function(x){
      return(rnorm(1,mean=d_summarized$mean_abs[x],sd=d_summarized$sd_abs[x]))
    }) #get the values of distance aside this site
    
    if (any (dist_other<0)){
      dist_other=dist_other[-which(dist_other<0)] #can happen
    }
    
    RMSE_dist_within = mean(sapply(1:length(dist_site),function(x){
      return(sqrt(sum((dist_site-dist_site[x])**2,na.rm = T)/length(dist_site) ))})) #NRMSE of p within
    
    RMSE_dist_all = mean(sapply(1:length(dist_site),function(x){
      return(sqrt(sum((dist_other-dist_site[x])**2,na.rm = T)/length(dist_other) ))})) #NRMSE of p between
    
    NRMSE_abs=RMSE_dist_within/RMSE_dist_all
    
    
    
    d_ABC=rbind(d_ABC,tibble(NRMSE_p=NRMSE_p,NRMSE_q=NRMSE_q,
                             NRMSE_abs=NRMSE_abs,
                             NRMSE_rela=NRMSE_rela,
                             Plot_n=k,ID_rep=rep_id))
    
  }
  
  #Using the posterior median
  #First parameters 
  
  
  p_site=d_summarized$median_p[which(d_summarized$Plot_n==k)]
  q_site=d_summarized$median_q[which(d_summarized$Plot_n==k)]
  
  p_other=d_summarized$median_p[-which(d_summarized$Plot_n==k)]
  q_other=d_summarized$median_q[-which(d_summarized$Plot_n==k)]
  
  RMSE_p_within = mean(sapply(1:length(p_site),function(x){
    return(sqrt(sum((p_site-p_site[x])**2,na.rm = T)/length(p_site) ))})) #NRMSE of p within
  
  RMSE_q_within = mean(sapply(1:length(q_site),function(x){
    return(sqrt(sum((q_site-q_site[x])**2,na.rm = T)/length(q_site) ))})) #NRMSE of q within
  
  RMSE_p_all = mean(sapply(1:length(p_site),function(x){
    return(sqrt(sum((p_other-p_site[x])**2,na.rm = T)/length(p_other) ))})) #NRMSE of p between
  RMSE_q_all = mean(sapply(1:length(q_site),function(x){
    return(sqrt(sum((q_other-q_site[x])**2,na.rm = T)/length(q_other) ))})) #NRMSE of q between
  
  NRMSE_p=RMSE_p_within/RMSE_p_all
  NRMSE_q=RMSE_q_within/RMSE_q_all
  
  
  
  #second for resilience related predictions
  
  #1) relative distance
  
  dist_site=d_summarized$median_rela[which(d_summarized$Plot_n==k)]
  dist_other=d_summarized$median_rela[-which(d_summarized$Plot_n==k)]
  
  RMSE_dist_within = mean(sapply(1:length(dist_site),function(x){
    return(sqrt(sum((dist_site-dist_site[x])**2,na.rm = T)/length(dist_site) ))})) #NRMSE of p within
  
  RMSE_dist_all = mean(sapply(1:length(dist_site),function(x){
    return(sqrt(sum((dist_other-dist_site[x])**2,na.rm = T)/length(dist_other) ))})) #NRMSE of p between
  
  NRMSE_rela=RMSE_dist_within/RMSE_dist_all
  
  #2) absolute distance
  
  dist_site=d_summarized$median_abs[which(d_summarized$Plot_n==k)]
  dist_other=d_summarized$median_abs[-which(d_summarized$Plot_n==k)]
  
  RMSE_dist_within = mean(sapply(1:length(dist_site),function(x){
    return(sqrt(sum((dist_site-dist_site[x])**2,na.rm = T)/length(dist_site) ))})) #NRMSE of p within
  
  RMSE_dist_all = mean(sapply(1:length(dist_site),function(x){
    return(sqrt(sum((dist_other-dist_site[x])**2,na.rm = T)/length(dist_other) ))})) #NRMSE of p between
  
  NRMSE_abs=RMSE_dist_within/RMSE_dist_all
  
  d_median=rbind(d_median,tibble(NRMSE_p=NRMSE_p,NRMSE_q=NRMSE_q,
                                 NRMSE_abs=NRMSE_abs,
                                 NRMSE_rela=NRMSE_rela,
                                 Plot_n=k))
}

write.table(d_median,"./Data/Similarity_within_sites_median.csv",sep=";")
write.table(d_ABC,"./Data/Similarity_within_sites_ABC_uncertainty.csv",sep=";")

# ---------------------- Step 6: Validation using simulations ----
## >> 6.1) Kefi dryland model ----
#First the distance predicted by the kefi model


dist_kefi=tibble()
for (site in list.files("./Data/Model_confirmation_Kefi/Dist_kefi","Dist")){
  
  site_id=as.numeric(gsub(".csv","",strsplit(site,"_")[[1]][3]))
  
  pred=read.table(paste0("./Data/Model_confirmation_Kefi/Dist_kefi/",site),sep=",")%>%
    filter(., V1>0)
  colnames(pred)=c("f","b","delta","cover")
  pred$cover[pred$cover<.01]=0
  
  if (max(pred$cover)>.05){
    dist_kefi=rbind(dist_kefi,tibble(Site=site_id,
                                     size_tipping=min(pred$cover[pred$cover !=0]),
                                     abs_dist=max(pred$b)-min(pred$b[pred$cover!=0]),
                                     relativ_dist=(max(pred$b)-min(pred$b[pred$cover!=0]))/min(pred$b[pred$cover!=0]),
                                     f=unique(pred$f),delta=unique(pred$delta),
                                     Cover=max(pred$cover)))
  }
}


dist_eby=tibble();step_size=0.005
for (site in list.files("./Data/Model_confirmation_Kefi/Dist_Eby","Dist")){
  site_id=as.numeric(gsub("Eby","",gsub(".csv","",strsplit(site,"_")[[1]][3])))
  
  pred=read.table(paste0("./Data/Model_confirmation_Kefi/Dist_Eby/",site),sep=",")%>%
    filter(., V2>0)
  colnames(pred)=c("p","q","cover")
  pred$cover[pred$cover<0.01] = 0
  
  if (max(pred$cover)>.05){
    index=0;pred$ID_sim=NA
    for (x in 1:nrow(pred)){
      pred$ID_sim[x]=index
      if (pred$p[x]==0){
        index=index+1
      }
    }
    
    p_desert=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x,cover>0)
      if (any(d_fil$cover>0)){
        return(d_fil$p[nrow(d_fil)])
      }else {
        return(NA)
      }
    }) 
    
    p_infer=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      return(d_fil$p[1])
    }) 
    
    q_infer=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      return(d_fil$q[1])
    }) 
    
    max_cover=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      return(d_fil$cover[1])
    }) 
    
    size_tipping=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x,cover>0)
      if (any(d_fil$cover>0)){
        return(d_fil$cover[nrow(d_fil)])
      }else {
        return(NA)
      }
    }) 
    
    dist_eby=rbind(dist_eby,tibble(Site=site_id,ID_sim=1:length(p_desert),
                                   abs_dist=p_infer-p_desert,
                                   relativ_dist=(p_infer-p_desert)/p_desert,
                                   size_tipping=size_tipping,
                                   Cover=max_cover))
  }else{
    dist_eby=rbind(dist_eby,tibble(Site=site_id,ID_sim=0,
                                   abs_dist=0,
                                   relativ_dist=0,
                                   size_tipping=0,
                                   Cover=0))
  }
}

#summarize by quantiles, mean and sd
dist_eby=dist_eby%>%dplyr::arrange(., Site)%>%filter(., Cover !=0)%>%
  filter(., Cover <.9)

d_eby=dist_eby%>%
  dplyr::group_by(., Site)%>%
  dplyr::summarize(., .groups = "keep",
                   cover_q1=quantile(Cover,.25,na.rm=T),cover_q3=quantile(Cover,.75,na.rm=T),
                   Cover=mean(Cover,na.rm=T),
                   abs_q1=quantile(abs_dist,.25,na.rm=T),abs_q3=quantile(abs_dist,.75,na.rm=T),
                   rela_q1=quantile(relativ_dist,.25,na.rm=T),rela_q3=quantile(relativ_dist,.75,na.rm=T),
                   abs_dist=median(abs_dist,na.rm=T),mean_abs_dist=mean(abs_dist,na.rm=T),
                   mean_rela_dist=mean(relativ_dist,na.rm=T),mean_size_tipping=mean(size_tipping,na.rm=T),
                   sd_rela_dist=sd(relativ_dist,na.rm = T),
                   relativ_dist=median(relativ_dist,na.rm=T))%>%
  add_column(., Model="Eby")%>%
  arrange(., Site)

sd_abs_dist=sapply(unique(dist_eby$Site),function(x){
  return(sd(dist_eby$abs_dist[which(dist_eby$Site==x)],na.rm = T))
})
d_eby$sd_abs_dist=sd_abs_dist

d_kefi=dist_kefi%>%dplyr::select(., -f,-delta)%>%
  add_column(., abs_q1=0,abs_q3=0,rela_q1=0,rela_q3=0,Model="Kefi")%>%
  arrange(., Site)

d_kefi=d_kefi%>%
  filter(., Site %in% unique(d_eby$Site))

# As there is uncertainty around each inferred distance to the desertification point, 
# we assume that the each inference ~ follow a normal distrib with mean = mean posterior samples
# and sd = sd posterior samples. We then compute the spearman correlation and its associated p-value



param_kefi=read.table("./Data/Model_confirmation_Kefi/Parameters_kefi.csv",sep=";")
colnames(param_kefi)=c("r","d","m","c","f","b","delta")
d_kefi=cbind(d_kefi,param_kefi[d_kefi$Site,])
d_eby=cbind(d_eby,param_kefi[d_kefi$Site,])

d_spearman=tibble()

id=1;  n=200
for (f_ in unique(param_kefi$f)){
  
  for (x in 1:n){
    
    d1=dplyr::filter(d_kefi,f==f_)
    d2=dplyr::filter(d_eby,f==f_)
    
    #absolute distance
    abs_eby=sapply(1: nrow(d2),function(x){
      return(rnorm(1,mean=d2$mean_abs_dist[x],sd=d2$sd_abs_dist[x]))
    })
    r_spearman=cor.test(d1$abs_dist,abs_eby,method = "spearman",exact = F)
    d_spearman=rbind(d_spearman,tibble(Stat=r_spearman$estimate,
                                       Pval=r_spearman$p.value,Type_dist="Abs",
                                       ID_sim=id))
    
    #relative distance
    rela_eby=sapply(1: nrow(d2),function(x){
      return(rnorm(1,mean=d2$mean_rela_dist[x],sd=d2$sd_rela_dist[x]))
    })
    r_spearman=cor.test(d1$relativ_dist,rela_eby,method = "spearman",exact = F)
    d_spearman=rbind(d_spearman,tibble(Stat=r_spearman$estimate,
                                       Pval=r_spearman$p.value,Type_dist="Rela",
                                       ID_sim=id))
    
  }
  id=id+1
  
}

saveRDS(list(d_spearman=d_spearman,
             d_kefi=d_kefi,
             d_eby=d_eby),
        "./Data/Model_confirmation_Kefi/d_for_figure.rds")



## >> 6.2) Guichard mussel-bed model ----


#First the distance predicted by the guichard model


dist_guichard=tibble()
for (site in list.files("./Data/Model_confirmation_Guichard/Dist_guichard","Dist")){
  
  site_id=as.numeric(gsub(".csv","",strsplit(site,"_")[[1]][3]))
  
  pred=read.table(paste0("./Data/Model_confirmation_Guichard/Dist_guichard/",site),sep=",")%>%
    filter(., V4>0)
  colnames(pred)=c("d","a0","a2","cover")
  pred$cover[pred$cover<.01]=0
  
  if (max(pred$cover)>.05){
    
    pred=pred%>%filter(., cover>0)
    
    dist_guichard=rbind(dist_guichard,tibble(Site=site_id,
                                             size_tipping=min(pred$cover[pred$cover !=0]),
                                             abs_dist=abs(diff(range(pred$d))),
                                             relativ_dist=abs(diff(range(pred$d)))/(max(pred$d)),
                                             a0=unique(pred$a0),a2=unique(pred$a2),
                                             Cover=max(pred$cover)))
  }
}


d_guichard=dist_guichard%>%dplyr::select(., -a0,-a2)%>%
  add_column(., abs_q1=0,abs_q3=0,rela_q1=0,rela_q3=0,Model="Guichard")%>%
  arrange(., Site)%>%
  filter(., Cover<.9)


dist_eby=tibble();step_size=0.005
for (site in list.files("./Data/Model_confirmation_Guichard/Dist_Eby","Dist")){
  
  site_id=as.numeric(gsub("Eby","",gsub(".csv","",strsplit(site,"_")[[1]][3])))
  
  pred=read.table(paste0("./Data/Model_confirmation_Guichard/Dist_Eby/",site),sep=",")%>%
    filter(., V2>0)
  colnames(pred)=c("p","q","cover")
  pred$cover[pred$cover<0.01] = 0
  
  if (max(pred$cover)>.05){
    index=0;pred$ID_sim=NA
    for (x in 1:nrow(pred)){
      pred$ID_sim[x]=index
      if (pred$p[x]==0){
        index=index+1
      }
    }
    
    p_desert=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x,cover>0)
      if (any(d_fil$cover>0)){
        return(d_fil$p[nrow(d_fil)])
      }else {
        return(NA)
      }
    }) 
    
    p_infer=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      return(d_fil$p[1])
    }) 
    
    q_infer=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      return(d_fil$q[1])
    }) 
    
    max_cover=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      return(d_fil$cover[1])
    }) 
    
    size_tipping=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x,cover>0)
      if (any(d_fil$cover>0)){
        return(d_fil$cover[nrow(d_fil)])
      }else {
        return(NA)
      }
    }) 
    
    
    dist_eby=rbind(dist_eby,tibble(Site=site_id,ID_sim=1:length(p_desert),
                                   abs_dist=p_infer-p_desert,
                                   relativ_dist=(p_infer-p_desert)/p_desert,
                                   size_tipping=size_tipping,
                                   Cover=max_cover))
  }else{
    dist_eby=rbind(dist_eby,tibble(Site=site_id,ID_sim=0,
                                   abs_dist=0,
                                   relativ_dist=0,
                                   size_tipping=0,
                                   Cover=0))
  }
  
}

#summarize by quantiles, mean and sd
dist_eby=dist_eby%>%dplyr::arrange(., Site)%>%filter(., Cover !=0)%>%
  filter(., Cover <.9)%>%
  filter(., Site %in% unique(d_guichard$Site))

d_eby=dist_eby%>%
  dplyr::group_by(., Site)%>%
  dplyr::summarize(., .groups = "keep",
                   cover_q1=quantile(Cover,.25,na.rm=T),cover_q3=quantile(Cover,.75,na.rm=T),
                   Cover=mean(Cover,na.rm=T),
                   abs_q1=quantile(abs_dist,.25,na.rm=T),abs_q3=quantile(abs_dist,.75,na.rm=T),
                   rela_q1=quantile(relativ_dist,.25,na.rm=T),rela_q3=quantile(relativ_dist,.75,na.rm=T),
                   abs_dist=median(abs_dist,na.rm=T),mean_abs_dist=mean(abs_dist,na.rm=T),
                   mean_rela_dist=mean(relativ_dist,na.rm=T),mean_size_tipping=mean(size_tipping,na.rm=T),
                   sd_rela_dist=sd(relativ_dist,na.rm = T),
                   relativ_dist=median(relativ_dist,na.rm=T))%>%
  add_column(., Model="Eby")%>%
  arrange(., Site)

sd_abs_dist=sapply(unique(dist_eby$Site),function(x){
  return(sd(dist_eby$abs_dist[which(dist_eby$Site==x)],na.rm = T))
})
d_eby$sd_abs_dist=sd_abs_dist





# As there is uncertainty around each inferred distance to the desertification point, 
# we assume that the each inference ~ follow a normal distrib with mean = mean posterior samples
# and sd = sd posterior samples. We then compute the spearman correlation and its associated p-value

param_guichard=read.table("./Data/Model_confirmation_Guichard/Parameters_guichard.csv",sep=";")
colnames(param_guichard)=c("d","a0","a2")
d_guichard=cbind(d_guichard,param_guichard[d_guichard$Site,])%>%
  filter(., Site %in% unique(d_eby$Site))
d_eby=cbind(d_eby,param_guichard[d_guichard$Site,])

d_spearman=d_spearman2=tibble()
id=1;  n=200
for (a0_ in unique(param_guichard$a0)){
  for (a2_ in unique(param_guichard$a2)){
    
    
    for (x in 1:n){
      
      d1=dplyr::filter(d_guichard,a0==a0_,a2==a2_)
      d2=dplyr::filter(d_eby,a0==a0_,a2==a2_)
      
      #absolute distance
      abs_eby=sapply(1: nrow(d2),function(x){
        return(rnorm(1,mean=d2$mean_abs_dist[x],sd=d2$sd_abs_dist[x]))
      })
      r_spearman=cor.test(d1$abs_dist,abs_eby,method = "spearman",exact = F)
      d_spearman=rbind(d_spearman,tibble(Stat=r_spearman$estimate,
                                         Pval=r_spearman$p.value,Type_dist="Abs",
                                         ID_sim=id))
      
      #relative distance
      rela_eby=sapply(1: nrow(d2),function(x){
        return(rnorm(1,mean=d2$mean_rela_dist[x],sd=d2$sd_rela_dist[x]))
      })
      r_spearman=cor.test(d1$relativ_dist,rela_eby,method = "spearman",exact = F)
      d_spearman=rbind(d_spearman,tibble(Stat=r_spearman$estimate,
                                         Pval=r_spearman$p.value,Type_dist="Rela",
                                         ID_sim=id))
      
    }
    id=id+1
  }
}
id=1
for (a0_ in unique(param_guichard$a0)){
  for (x in 1:n){
    
    d1=dplyr::filter(d_guichard,a0==a0_)
    d2=dplyr::filter(d_eby,a0==a0_)
    
    #absolute distance
    abs_eby=sapply(1: nrow(d2),function(x){
      return(rnorm(1,mean=d2$mean_abs_dist[x],sd=d2$sd_abs_dist[x]))
    })
    r_spearman=cor.test(d1$abs_dist,abs_eby,method = "spearman",exact = F)
    d_spearman2=rbind(d_spearman2,tibble(Stat=r_spearman$estimate,
                                         Pval=r_spearman$p.value,Type_dist="Abs",
                                         ID_sim=id))
    
    #relative distance
    rela_eby=sapply(1: nrow(d2),function(x){
      return(rnorm(1,mean=d2$mean_rela_dist[x],sd=d2$sd_rela_dist[x]))
    })
    r_spearman=cor.test(d1$relativ_dist,rela_eby,method = "spearman",exact = F)
    d_spearman2=rbind(d_spearman2,tibble(Stat=r_spearman$estimate,
                                         Pval=r_spearman$p.value,Type_dist="Rela",
                                         ID_sim=id))
    
  }
  id=id+1
}

saveRDS(list(d_spearman=d_spearman,
             d_spearman2=d_spearman2,
             d_guichard=d_guichard,
             d_eby=d_eby),
        "./Data/Model_confirmation_Guichard/d_for_figure.rds")




# ---------------------- Step 7: Bifurcation diagram with q ----

d=read.table("./Data/posterior_param.csv",sep=";",header=T)
keep_sites=read.table("./Data/Keeping_sites_biocom.csv",sep=";")$V1
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
         Cover=d_biocom$Cover)%>%
  dplyr::filter(., Site %in% keep_sites)


d2=read.table("./Data/Resilience_metrics_1_neigh_with_q.csv",sep=";")%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_median=median(qinfer-qcrit,na.rm = T),
                   row_at_desert=min(row_at_desert,na.rm = T)
  )%>% #when decreasing q, we don't go necessarly t
  dplyr::filter(., Site %in% keep_sites,row_at_desert !=1)

d_with_q=cbind(d%>%filter(., Site%in%d2$Site),d2)%>%arrange(., Site)


d2=read.table("./Data/Resilience_metrics_1_neigh.csv",sep=";")%>%
  dplyr::group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_median=median(pinfer-pcrit,na.rm = T))%>%
  dplyr::filter(., Site %in% d_with_q$Site)

d_with_p=cbind(d%>%filter(., Site%in%d2$Site),d2)%>%arrange(., Site)


p=ggplot(NULL)+
  geom_point(aes(x=d_with_p$abs_median,y=d_with_q$abs_median))+
  the_theme+
  labs(x="Absolute distance to the tipping point (decreasing p)",
       y="Absolute distance to the tipping point (decreasing q)")+
  geom_abline(slope = 1,intercept = 0,color="red")




# ---------------------- Step 8: Climatic projections ----

library(stars)     
library(sf)        
library(raster)

for (file_ in list.files("./Data/Climatic_data",pattern = ".nc")[-c(1,2)]){
  
  nc_data_aridity = raster::rotate(raster(paste0("./Data/Climatic_data/aidm_0.5_rcp45_yreof_tsreg_200601-210012.nc"))) 
  dryland_sf = st_as_sf(d_biocom[keep_sites,], coords = c("Longitude", "Lattitude"), crs = 4326)
  
  aridity_values = raster::extract(aridity_years, dryland_sf,method="bilinear")
  
  
  fun=function(x) { 
    if (any(is.na(x))){
      NA
    }else{
      lin_mod=lm(x ~ as.numeric(seq(1950, 2100, by=1)[1:ncol(aridity_years)]));
      summary(lin_mod)$coefficients[2]
    }
  }
  
  aridity_years=extract(aridity_years, SpatialPoints(d_biocom[,c("Longitude","Lattitude")]))
  
  trends_aridity = apply(aridity_years,1,fun)
  
  write.table(trends_aridity,paste0("./Data/Climatic_data/trend_clim_",gsub(".nc","",file_),".csv"),sep=";")

}

test=nc_data_aridity
test_spdf <- as(test, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")



simulated_sites=sapply(1:length(list.files("./Data/Prediction/","Dist")),
            function(x){
              return(as.numeric(gsub(".csv","",strsplit(list.files("./Data/Prediction/","Dist")[x],"_")[[1]][3])))})

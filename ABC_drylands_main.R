rm(list=ls())
source("./ABC_drylands_function.R")



# ---------------------- Step 0: Merging simulations ----


d_simu=tibble()
n_param=2

list_sim=list.files("./Data/Simulations/",".csv")

d=tibble()
for (sim in list_sim){
  d=rbind(d,read.table(paste0("./Data/Simulations/",sim),sep=",")[,(1):(n_param+11)])%>%
    filter(., V3>0.03)
}

colnames(d)=c("p","q","rho_p","nb_neigh","clustering","skewness","variance","moran_I",
              "Spectral_ratio","PLR","PL_expo","cv_psd","median_psd","mean_psd","sd_psd","fmax_psd")


d=add_column(d,Pooling=rep(1:5,nrow(d)/5))
for (i in 1:(nrow(d)/5)){
  d$PL_expo[(5*(i-1)+2):(5*(i-1)+5)]=d$PL_expo[(5*(i-1)+1)]
  d$PLR[(5*(i-1)+2):(5*(i-1)+5)]=d$PLR[(5*(i-1)+1)]
  d$cv_psd[(5*(i-1)+2):(5*(i-1)+5)]=d$cv_psd[(5*(i-1)+1)]
  d$fmax_psd[(5*(i-1)+2):(5*(i-1)+5)]=d$fmax_psd[(5*(i-1)+1)]
}

write.table(d,"./Data/simulations.csv",sep=";")






# ---------------------- Step 1: Improving inference ----
## >> 1) Optimizing pre and post-processing methods ----
#we play on both the lambda of the boxcox method and the size of the sample for
#stage 1 of the two step pre-processing procedure

dir.create("../Data_new/NRMSE",showWarnings = F)

d_all=read.table("../Data_new/All_new_sim.csv",sep=";")%>%
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
        
        write.table(d_NRMSE_param,paste0("../Data_new/NRMSE/RMSE_param_",preprocessing,
                                         "_",method_abc,"_optim_lambda_",ifelse(optim_lambda,"yes","no"),"_N1_",size_step1,".csv"),sep=";")
        
        colnames(d_NRMSE_sumstat)=colnames(d_cross_sumstat)
        
        write.table(d_NRMSE_sumstat,paste0("../Data_new/NRMSE/RMSE_sumstat_",preprocessing,
                                           "_",method_abc,"_optim_lambda_",ifelse(optim_lambda,"yes","no"),"_N1_",size_step1,".csv"),sep=";")
        
        

        
        
        
      }
    }
  }
}





## >> 2) Optimizing the structure of the neural-network ----

dir.create("../Data_new/NRMSE",showWarnings = F)


d_all=read.table("../Data_new/All_new_sim.csv",sep=";")%>%
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
      
      write.table(d_NRMSE_param,paste0("../Data_new/NRMSE/RMSE_hidden_preprocessing_",method_pre,"_",
                                       size_hidden,"_Nnet_",rep_network,".csv"),sep=";")
      
      
      
    }
  }
}





## >> 3) Influence of the number of simulation kept ----

dir.create("../Data_new/NRMSE",showWarnings = F)
d_all=read.table("../Data_new/All_new_sim.csv",sep=";")%>%
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
  
  write.table(d_NRMSE_param,paste0("./Data_new/NRMSE/RMSE_NAkept_",NA_kept,".csv"),sep=";")
  
}






## >> 4) Selecting the best summary statistics ----

dir.create("../Data_new/Best_sumstat",showWarnings = F)

d_all=read.table("../Data_new/All_new_sim.csv",sep=";")%>%
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
  
  write.table(d_NRMSE_param,paste0("../Data_new/Best_sumstat/RMSE_selecting_sumstat_,",name_removal,".csv"),sep=";")

}






# ---------------------- Step 2: Spatial resolution ----
d_sim=read.table("../Data_new/All_new_sim.csv",sep=";")

p=d_sim%>%
  mutate(., Id_sim=rep(1:(nrow(.)/5),each=5))%>%
  add_column(., Model="Eby")%>%
  mutate(., Pooling=as.character(Pooling))%>%
  melt(., id.vars=c("Pooling","Id_sim","Model"))%>%
  filter(., variable %in% c("nb_neigh","clustering","skewness","variance","moran_I","Spectral_ratio"))%>%
  mutate(., Model=recode_factor(Model,"Eby_feedback"="Eby with feedbacks"))%>%
  group_by(., Model,Pooling,variable)%>%
  summarise(., .groups = "keep",mean_value=mean(value,na.rm = T))%>%
  ggplot(.)+
  geom_line(aes(x=Pooling,y=mean_value,group=interaction(Model),color=Model),lwd=1)+
  geom_point(aes(x=Pooling,y=mean_value,color=Model),fill="white",shape=21,size=2.5)+
  facet_wrap(.~variable,scales = "free")+
  labs(x="Change in resolution",y="Mean value",color="Model")+
  scale_color_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))+
  scale_fill_manual(values=c("#C46FC5","#80BD5C","#568DC5","#DE6450","#898C86"))+
  scale_x_discrete(labels=c("No change","x2","x3","x4","x5"))+
  the_theme+theme(axis.text.x = element_text(hjust=1,angle=60))+
  guides(color = guide_legend(override.aes = list(size = 2)),fill="none")


## >> 1) Collecting data and merging data ----


d_simu=tibble()
n_param=2

list_sim=list.files("../Data/Step10_ABC_scale/Sim_new",".csv")

d=tibble()
for (sim in list_sim){
  d=rbind(d,read.table(paste0("../Data/Step10_ABC_scale/Sim_new","/",sim),sep=",")[,(1):(n_param+14)])%>%
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

write.table(d,"../Data_new/All_new_sim.csv",sep=";")


## >> 2) Can we recover parameters and scale of observation in simulations ----

dir.create("../Data_new/Scale_obs_indentifiability",showWarnings = F)
d_sim=read.table("../Data_new/All_new_sim.csv",sep=";")%>%
  add_column(.,ID_sim=rep(1:(nrow(.)/5),each=5))%>%
  filter(., rho_p>0.03)%>%
  dplyr::relocate(.,Pooling,.after=q)

set.seed(123)
list_sample=sample(unique(d_sim$ID_sim),100) #for each of these 100 parameter sets, we need to verify whether we infer the same parameters of not
d_RMSE_param=x_y_param=tibble()
n_keep=75

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
    
    
    if (any(is.na(target))){
      
      which_na=which(is.na(target))
      abc_sim=abc(target=target[,stat_kept[-(which_na-3)]],
                  param = all_sim[,c(1:3)],sumstat = all_sim[,stat_kept[-(which_na-3)]],
                  tol = n_keep/nrow(all_sim),method="neuralnet",transf = c(rep("logit",2),"none"),
                  logit.bounds = matrix(c(0,1),3,2,byrow = T),
                  numnet = 10,sizenet = 10)
      
    }else {
      abc_sim=abc(target=target[,stat_kept],
                  param = all_sim[,1:3],sumstat = all_sim[,stat_kept],
                  tol = n_keep/nrow(all_sim),method="neuralnet",transf = c(rep("logit",2),"none"),
                  logit.bounds = matrix(c(0,1),3,2,byrow = T),
                  numnet = 10,sizenet = 10)
    }
    
    RMSE = sapply(1:ncol(abc_sim$adj.values),function(x){
      sqrt(sum((abc_sim$adj.values[,x]-target[,x])**2)/nrow(abc_sim$adj.values) )})
    
    RMSE_prior=sapply(1:3,function(x){
      sqrt(sum(((all_sim[,x]-target[,x])**2)/nrow(all_sim) ))})
    
    NMRSE=RMSE/RMSE_prior
    names(NMRSE)=c("p","q","Pooling")
    
    d_RMSE_param=rbind(d_RMSE_param,data.frame((t(NMRSE)))%>%
                         add_column(., Site_ID=which(list_sample == sample_id),Method="NeuralNet",Scale_obs=scale_obs))
    
    x_y_param=rbind(x_y_param,data.frame(t(colMeans(abc_sim$adj.values)))%>%
                      add_column(., Site_ID=which(list_sample == sample_id),Method="NeuralNet",Type="Sim",Scale_obs=scale_obs))
    
    x_y_param=rbind(x_y_param,target[,1:3]%>%
                      add_column(., Site_ID=which(list_sample == sample_id),Method="NeuralNet",Type="Obs",Scale_obs=scale_obs))
  }
}


write.table(d_RMSE_param,"../Data_new/Scale_obs_indentifiability/Retrieving_parameters_different_resolution_RMSE_param.csv",sep=";")
write.table(x_y_param,"../Data_new/Scale_obs_indentifiability/Retrieving_parameters_different_resolution_x_y.csv",sep=";")






# ---------------------- Step 3: Selecting relevant empirical data ----

#Testing for bimodality in the posterior distribution of parameters
library(diptest)
post_param=read.table("../Data_new/posterior_param.csv",sep=";")

d_biocom$bimod=sapply(1:nrow(d_biocom),function(x){
  if (dip.test(post_param[,x])$p.value<.05 | dip.test(post_param[,x+345])$p.value<.05){
    return("bimod")
  }else {
    return("unimod")
  }
})
write.table(which(d_biocom$bimod!="bimod"),"../Data_new/Keeping_sites.csv",sep=";")

# ---------------------- Step 4: Inference, what we learn ----



run_abc=function(id,n_sim_kept=100,method_abc="neuralnet"){
  for (id_plot in id){
    `%!in%` = Negate(`%in%`)
    n_param=3
    
    d_biocom=read.table("../Data_new/biocom_data.csv",sep=";")
    d_sim=read.table("../Data_new/All_new_sim.csv",sep=";")%>%
      dplyr::relocate(., Pooling,.after =q )%>%
      filter(., PL_expo>0)
    
    rownames(d_sim)=1:nrow(d_sim)
    
    
    d_param_infer_NN=array(0,c(n_sim_kept,nrow(d_biocom),3))
    d_param_infer_rej=array(0,c(n_sim_kept,nrow(d_biocom),3))
    
    d_NRMSE_sumstat=x_y_stat=tibble()
    
    name_plot=c("all","no_cv","no_cv_fmax","no_PL_PLR","not_the_4","only_fmax","no_PL","no_PLR")
    
    sumstat_kept=4:14
    
    
    if (id_plot==3){
      sumstat_kept=sumstat_kept[-c(10:11)]
    }else if (id_plot==4){
      sumstat_kept=sumstat_kept[-c(8,9)]
    }else if (id_plot==5){
      sumstat_kept=sumstat_kept[-c(8:11)]
    }else if (id_plot==6){
      sumstat_kept=sumstat_kept[-c(8:10)]
    }else if (id_plot==7){
      sumstat_kept=sumstat_kept[-c(9)]
    }else if (id_plot==8){
      sumstat_kept=sumstat_kept[-c(8)]
    }else if (id_plot==2){
      sumstat_kept=sumstat_kept[-c(10)]
    }
    
    
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
      
      d_param_infer_NN[,empirical_id,1]=cross_valid$adj.values[,1] # we keep the whole distribution for p
      d_param_infer_NN[,empirical_id,2]=cross_valid$adj.values[,2] # for q
      d_param_infer_NN[,empirical_id,3]=cross_valid$adj.values[,3] # for the scale of observation
      
      d_param_infer_rej[,empirical_id,1]=cross_valid$unadj.values[,1] # we keep the whole distribution for p
      d_param_infer_rej[,empirical_id,2]=cross_valid$unadj.values[,2] # for q
      d_param_infer_rej[,empirical_id,3]=cross_valid$unadj.values[,3] # for the scale of observation
      
      
    }
    
    write.table(d_NRMSE_sumstat,paste0("../Data_new/Inferrence/NRMSE_sumstat_",name_plot[id_plot],".csv"),sep=";")
    write.table(x_y_stat,paste0("../Data_new/Inferrence/x_y_stat_",name_plot[id_plot],".csv"),sep=";")
    write.table(d_param_infer_NN,paste0("../Data_new/Inferrence/param_",ifelse(method_abc=="neuralnet","NN","loc"),
                                        "_",name_plot[id_plot],"_",n_sim_kept,".csv"),sep=";")
    write.table(d_param_infer_rej,paste0("../Data_new/Inferrence/param_rej_",name_plot[id_plot],".csv"),sep=";")
    
  }
  
}
library(parallel)
mclapply(1:8,run_abc,mc.cores = 8)



#Getting the modes of the posterior of each site

d=tibble()
param_infer_rej=read.table(paste0("../Data_new/posterior_param.csv"),sep=";")
d=rbind(d,tibble(Site=1:345,N_keep=50,
                 p_50=apply(param_infer_rej[1:345],2,quantile,.5),
                 p_25=apply(param_infer_rej[1:345],2,quantile,.25),
                 p_75=apply(param_infer_rej[1:345],2,quantile,.75),
                 q_50=apply(param_infer_rej[346:690],2,quantile,.5),
                 q_25=apply(param_infer_rej[346:690],2,quantile,.25),
                 q_75=apply(param_infer_rej[346:690],2,quantile,.75),
                 Scale=apply(param_infer_rej[691:1035],2,mean)))

write.table(d,"../Data_new/Inferrence/Posterior_modes_each_sites.csv",sep=";",row.names = F,col.names = F)


#Ploting histogram of each posterior distribution of the sites

param_infer_rej=read.table(paste0("../Data_new/posterior_param.csv"),sep=";")

pdf("../Figures/Posterior_sites.pdf",width = 12,height = 5)

par(mfrow=c(1,3))
for (i in 1:345){
  hist(param_infer_rej[,i],main="p",xlab="",ylab="",col=alpha("blue",.4))
  hist(param_infer_rej[,i+345],main="q",xlab="",ylab="",col=alpha("green",.4))
  hist(param_infer_rej[,i+690],main="scale",xlab="",ylab="",col=alpha("red",.4))
  mtext(paste0(mean(param_infer_rej[,i+690])))
}

dev.off()


## >> 1) Running mixed-effect models ----


d=read.table("../Data_new/posterior_param.csv",sep=";",header=T)
d_raw=readxl::read_xlsx("../Data_new/biocom_raw.xlsx")
keep_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
d=tibble(Site=1:345,mean_p=apply(d[,1:345],2,mean),sd_p=apply(d[,1:345],2,sd),
         mean_q=apply(d[,346:690],2,mean),sd_q=apply(d[,346:690],2,sd),
         Plot_n=d_biocom$Plot_n,
         Aridity=d_biocom$Aridity,
         Sand=d_biocom$Sand,
         MF=d_biocom$MF,
         Cover=d_biocom$Cover)%>%
  filter(., Site %in% keep_sites)%>%
  add_column(., 
             Facilitation=sapply(1:nrow(.),function(x){return(d_raw$Facil[which(d_raw$plotn==.$Plot_n[x])])}),
             SR=sapply(1:nrow(.),function(x){return(d_raw$SR[which(d_raw$plotn==.$Plot_n[x])])}),
             Long_cos=sapply(1:nrow(.),function(x){return(cos(d_raw$Longitude[which(d_raw$plotn==.$Plot_n[x])]))}),
             Long_sin=sapply(1:nrow(.),function(x){return(sin(d_raw$Longitude[which(d_raw$plotn==.$Plot_n[x])]))}),
             Lat=sapply(1:nrow(.),function(x){return(d_raw$Latitude[which(d_raw$plotn==.$Plot_n[x])])}),
             Elevation=sapply(1:nrow(.),function(x){return(d_raw$ELE_1[which(d_raw$plotn==.$Plot_n[x])])}),
             Slope=sapply(1:nrow(.),function(x){return(d_raw$SLO[which(d_raw$plotn==.$Plot_n[x])])})
  )

d2=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")%>%
  group_by(., Site,MF,aridity,Sand)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  filter(., Site %in% keep_sites)

d=cbind(d,d2)

#as the parameters have uncertainty -> Monte Carlo approach

nsim=1000
d_mod=tibble()
for (k in 1:nsim){
  
  d2=tibble(p=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_p[x],d$sd_p[x]))}))[,1],
            q=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_q[x],d$sd_q[x]))}))[,1],
            Size_tipping=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$Size_mean[x],d$Size_sd[x]))}))[,1],
            abs_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$abs_mean[x],d$abs_sd[x]))}))[,1],
            rela_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$relativ_mean[x],d$relativ_sd[x]))}))[,1],
            Sand=(d$Sand-mean(d$Sand,na.rm=T))/sd(d$Sand,na.rm = T),
            Site=d$Site,
            MF=(d$MF-mean(d$MF,na.rm=T))/sd(d$MF,na.rm = T),
            SR=(d$SR-mean(d$SR,na.rm=T))/sd(d$SR,na.rm = T),
            Facilitation=(d$Facilitation-mean(d$Facilitation,na.rm=T))/sd(d$Facilitation,na.rm = T),
            Cover=(d$Cover-mean(d$Cover,na.rm=T))/sd(d$Cover,na.rm = T),
            Aridity=(d$Aridity-mean(d$Aridity,na.rm=T))/sd(d$Aridity,na.rm = T),
            Lat=(d$Lat-mean(d$Lat,na.rm=T))/sd(d$Lat,na.rm = T),
            Long_cos=(d$Long_cos-mean(d$Long_cos,na.rm=T))/sd(d$Long_cos,na.rm = T),
            Long_sin=(d$Long_sin-mean(d$Long_sin,na.rm=T))/sd(d$Long_sin,na.rm = T),
            Elevation=(d$Elevation-mean(d$Elevation,na.rm=T))/sd(d$Elevation,na.rm = T),
            Slope=(d$Slope-mean(d$Slope,na.rm=T))/sd(d$Slope,na.rm = T),
            Plot_n=d$Plot_n)
  
  mod_predictors=gsub("\n     ","","Aridity + MF + Sand + Cover + SR + Facilitation +
      Lat + Long_cos + Long_sin + Slope + Elevation + ( 1 | Plot_n)")
  
  model_p=summary(d2%>% lmer(formula = paste("p ~ ",mod_predictors), data = .))
  model_q=summary(d2%>% lmer(formula = paste("q ~ ",mod_predictors), data = .))
  model_abs=summary(d2%>% lmer(formula = paste("abs_dist ~ ",mod_predictors), data = .))
  model_rela=summary(d2%>% lmer(formula = paste("rela_dist ~ ",mod_predictors), data = .))
  model_size=summary(d2%>% lmer(formula = paste("Size_tipping ~ ",mod_predictors), data = .))
  
  d_mod=rbind(d_mod,tibble(Aridity=c(model_p$coefficients[2,1],model_q$coefficients[2,1],model_abs$coefficients[2,1],
                                     model_rela$coefficients[2,1],model_size$coefficients[2,1]),
                           Multifunctionality=c(model_p$coefficients[3,1],model_q$coefficients[3,1],model_abs$coefficients[3,1],
                                                model_rela$coefficients[3,1],model_size$coefficients[3,1]),
                           Sand=c(model_p$coefficients[4,1],model_q$coefficients[4,1],model_abs$coefficients[4,1],
                                  model_rela$coefficients[4,1],model_size$coefficients[4,1]),
                           Cover=c(model_p$coefficients[5,1],model_q$coefficients[5,1],model_abs$coefficients[5,1],
                                   model_rela$coefficients[5,1],model_size$coefficients[5,1]),
                           SR=c(model_p$coefficients[6,1],model_q$coefficients[6,1],model_abs$coefficients[6,1],
                                model_rela$coefficients[6,1],model_size$coefficients[6,1]),
                           Facil=c(model_p$coefficients[7,1],model_q$coefficients[7,1],model_abs$coefficients[7,1],
                                   model_rela$coefficients[7,1],model_size$coefficients[7,1]),
                           Lat=c(model_p$coefficients[8,1],model_q$coefficients[8,1],model_abs$coefficients[8,1],
                                 model_rela$coefficients[8,1],model_size$coefficients[8,1]),
                           Long_cos=c(model_p$coefficients[9,1],model_q$coefficients[9,1],model_abs$coefficients[9,1],
                                      model_rela$coefficients[9,1],model_size$coefficients[9,1]),
                           Long_sin=c(model_p$coefficients[10,1],model_q$coefficients[10,1],model_abs$coefficients[10,1],
                                      model_rela$coefficients[10,1],model_size$coefficients[10,1]),
                           Slope=c(model_p$coefficients[11,1],model_q$coefficients[11,1],model_abs$coefficients[11,1],
                                   model_rela$coefficients[11,1],model_size$coefficients[11,1]),
                           Elevation=c(model_p$coefficients[12,1],model_q$coefficients[12,1],model_abs$coefficients[12,1],
                                       model_rela$coefficients[12,1],model_size$coefficients[12,1]),
                           Param=c("p","q","Absolute distance","Relative distance","Size tipping")))
  
  
  print(k)
  
}

write.table(d_mod,"../Data_new/Prediction/Drivers_stability_metrics.csv",sep=";")



## >> 2) Evaluating importance ----


Run_model_importance=function(id){
  
  list_evalutate_models=expand.grid(Stat_param=c("p","q","rela_dist",
                                                 "Size_tipping","abs_dist"),
                                    Seed_id=1:10)
  
  stat=list_evalutate_models$Stat_param[id] #geting the same of the response var
  Seed_id=list_evalutate_models$Seed_id[id] #getting the seed ID
  
  d=read.table("../Data_new/posterior_param.csv",sep=";",header=T)
  d_raw=readxl::read_xlsx("../Data_new/biocom_raw.xlsx")
  keep_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
  d=tibble(Site=1:345,mean_p=apply(d[,1:345],2,mean),sd_p=apply(d[,1:345],2,sd),
           mean_q=apply(d[,346:690],2,mean),sd_q=apply(d[,346:690],2,sd),
           Plot_n=d_biocom$Plot_n,
           Aridity=d_biocom$Aridity,
           Sand=d_biocom$Sand,
           MF=d_biocom$MF,
           Cover=d_biocom$Cover)%>%
    filter(., Site %in% keep_sites)%>%
    add_column(., 
               Facilitation=sapply(1:nrow(.),function(x){return(d_raw$Facil[which(d_raw$plotn==.$Plot_n[x])])}),
               SR=sapply(1:nrow(.),function(x){return(d_raw$SR[which(d_raw$plotn==.$Plot_n[x])])}),
               Long_cos=sapply(1:nrow(.),function(x){return(cos(d_raw$Longitude[which(d_raw$plotn==.$Plot_n[x])]))}),
               Long_sin=sapply(1:nrow(.),function(x){return(sin(d_raw$Longitude[which(d_raw$plotn==.$Plot_n[x])]))}),
               Lat=sapply(1:nrow(.),function(x){return(d_raw$Latitude[which(d_raw$plotn==.$Plot_n[x])])}),
               Elevation=sapply(1:nrow(.),function(x){return(d_raw$ELE_1[which(d_raw$plotn==.$Plot_n[x])])}),
               Slope=sapply(1:nrow(.),function(x){return(d_raw$SLO[which(d_raw$plotn==.$Plot_n[x])])})
    )
  
  
  d2=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")%>%
    group_by(., Site,MF,aridity,Sand)%>%
    dplyr::summarise(., .groups = "keep",
                     abs_mean=mean(pinfer-pcrit,na.rm = T),
                     abs_sd=sd(pinfer-pcrit,na.rm = T),
                     relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                     relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                     Size_mean=mean(Size_tipping,na.rm = T),
                     Size_sd=sd(Size_tipping,na.rm = T))%>%
    filter(., Site %in% keep_sites)
  
  d=cbind(d,d2)
  
  d_importance=tibble()#for keeping all information
  
  for (seed in ((1+(Seed_id-1)*100):(Seed_id*100))){
    set.seed(seed)
    #as the parameters have uncertainty -> Monte Carlo approach
    d2=tibble(p=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_p[x],d$sd_p[x]))}))[,1],
              q=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_q[x],d$sd_q[x]))}))[,1],
              Size_tipping=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$Size_mean[x],d$Size_sd[x]))}))[,1],
              abs_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$abs_mean[x],d$abs_sd[x]))}))[,1],
              rela_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$relativ_mean[x],d$relativ_sd[x]))}))[,1],
              Sand=(d$Sand-mean(d$Sand,na.rm=T))/sd(d$Sand,na.rm = T),
              Site=d$Site,
              MF=(d$MF-mean(d$MF,na.rm=T))/sd(d$MF,na.rm = T),
              SR=(d$SR-mean(d$SR,na.rm=T))/sd(d$SR,na.rm = T),
              Facilitation=(d$Facilitation-mean(d$Facilitation,na.rm=T))/sd(d$Facilitation,na.rm = T),
              Cover=(d$Cover-mean(d$Cover,na.rm=T))/sd(d$Cover,na.rm = T),
              Aridity=(d$Aridity-mean(d$Aridity,na.rm=T))/sd(d$Aridity,na.rm = T),
              Lat=(d$Lat-mean(d$Lat,na.rm=T))/sd(d$Lat,na.rm = T),
              Long_cos=(d$Long_cos-mean(d$Long_cos,na.rm=T))/sd(d$Long_cos,na.rm = T),
              Long_sin=(d$Long_sin-mean(d$Long_sin,na.rm=T))/sd(d$Long_sin,na.rm = T),
              Elevation=(d$Elevation-mean(d$Elevation,na.rm=T))/sd(d$Elevation,na.rm = T),
              Slope=(d$Slope-mean(d$Slope,na.rm=T))/sd(d$Slope,na.rm = T),
              Plot_n=d$Plot_n)
    
    mod_predictors="Aridity + MF + Sand + Cover + SR + Facilitation +
        Lat + Long_cos + Long_sin + Slope + Elevation + (1|Plot_n)"
    
    
    d_data_mod=d2%>%melt(., measure.vars=stat)%>%
      filter(., !is.na(value), !is.na(Facilitation))
    
    model_predict=lmer(formula(paste("value ~",mod_predictors)),d_data_mod,
                       na.action = na.fail ,REML ="FALSE")
    
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_predict, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    
    model_predict = lmer(formula(paste("value ~",mod_predictors)), data = d_data_out, 
                         na.action = na.fail,REML ="FALSE")
    
    # #do some model selection
    select_model=dredge(model_predict, subset = ~ Slope & Elevation &
                          Long_cos & Long_sin & Lat,
                        options(na.action = "na.fail"),
                        extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                                R2c=function(x) r.squaredGLMM(x)[2])
    )
    
    if (nrow(select_model%>%filter(., AICc<min(AICc)+2))>1){
      R2=select_model%>%filter(., AICc<min(AICc)+2)
      
      #extract the result of model selection
      result_select=model.avg(select_model, subset = delta <2)
      #Get the importance of each metric
      importance_mod=sw(result_select)
      d_importance=rbind(d_importance,cbind(tibble(Response_var=stat,
                                                   N_outliers=rm.outliers$n.removed),
                                            Aggregate_importance(importance_mod))%>%
                           add_column(., Seed_ID=Seed_id,R2m=mean(R2$R2m),R2c=mean(R2$R2c)))
      
    }else{
      result_select=select_model%>%filter(., AICc==min(AICc))
      R2=  select_model%>%filter(., AICc==min(AICc))
      
      d_importance=rbind(d_importance,cbind(tibble(Response_var=stat,
                                                   N_outliers=rm.outliers$n.removed),
                                            Aggregate_importance(result_select,"1mod"))%>%
                           add_column(., Seed_ID=Seed_id,R2m=mean(R2$R2m),R2c=mean(R2$R2c)))
    }
    
    
  }
  
  
  #Saving each DF
  write.table(d_importance,paste0("../Data_new/Prediction/Importance/Importance_",stat,"_",Seed_id,".csv"),sep=";")
  
}

library(parallel)

mclapply(1:50,Run_model_importance,mc.cores = 25)


d_importance=tibble()
for (k in list.files("../Data_new/Prediction/Importance/",".csv")){
  d_importance=rbind(d_importance,read.table(paste0("../Data_new/Prediction/Importance/",k),sep=";"))
}
write.table(d_importance,"../Data_new/Prediction/Importance_drivers_stability_metrics.csv",sep=";")


## >> 3) Indirect SEM effects: aborted ----

Run_model_SEM=function(id){
  

  Seed_id=id #getting the seed ID
  
  d=read.table("../Data_new/posterior_param.csv",sep=";",header=T)
  d_raw=readxl::read_xlsx("../Data_new/biocom_raw.xlsx")
  keep_sites=read.table("../Data_new/Keeping_sites.csv",sep=";")$x
  d=tibble(Site=1:345,mean_p=apply(d[,1:345],2,mean),sd_p=apply(d[,1:345],2,sd),
           mean_q=apply(d[,346:690],2,mean),sd_q=apply(d[,346:690],2,sd),
           Plot_n=d_biocom$Plot_n,
           Aridity=d_biocom$Aridity,
           Sand=d_biocom$Sand,
           MF=d_biocom$MF,
           Cover=d_biocom$Cover)%>%
    filter(., Site %in% keep_sites)%>%
    add_column(., 
               Facilitation=sapply(1:nrow(.),function(x){return(d_raw$Facil[which(d_raw$plotn==.$Plot_n[x])])}),
               SR=sapply(1:nrow(.),function(x){return(d_raw$SR[which(d_raw$plotn==.$Plot_n[x])])}),
               Long_cos=sapply(1:nrow(.),function(x){return(cos(d_raw$Longitude[which(d_raw$plotn==.$Plot_n[x])]))}),
               Long_sin=sapply(1:nrow(.),function(x){return(sin(d_raw$Longitude[which(d_raw$plotn==.$Plot_n[x])]))}),
               Lat=sapply(1:nrow(.),function(x){return(d_raw$Latitude[which(d_raw$plotn==.$Plot_n[x])])}),
               Elevation=sapply(1:nrow(.),function(x){return(d_raw$ELE_1[which(d_raw$plotn==.$Plot_n[x])])}),
               Slope=sapply(1:nrow(.),function(x){return(d_raw$SLO[which(d_raw$plotn==.$Plot_n[x])])})
    )
  
  
  d2=read.table("../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")%>%
    group_by(., Site,MF,aridity,Sand)%>%
    dplyr::summarise(., .groups = "keep",
                     abs_mean=mean(pinfer-pcrit,na.rm = T),
                     abs_sd=sd(pinfer-pcrit,na.rm = T),
                     relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                     relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                     Size_mean=mean(Size_tipping,na.rm = T),
                     Size_sd=sd(Size_tipping,na.rm = T))%>%
    filter(., Site %in% keep_sites)
  
  d=cbind(d,d2)
  
  d_importance=tibble()#for keeping all information
  
  for (seed in ((1+(Seed_id-1)*100):(Seed_id*100))){
    
    set.seed(seed)
    
    #as the parameters have uncertainty -> Monte Carlo approach
    d2=tibble(p=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_p[x],d$sd_p[x]))}))[,1],
              q=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_q[x],d$sd_q[x]))}))[,1],
              Size_tipping=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$Size_mean[x],d$Size_sd[x]))}))[,1],
              abs_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$abs_mean[x],d$abs_sd[x]))}))[,1],
              rela_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$relativ_mean[x],d$relativ_sd[x]))}))[,1],
              Sand=(d$Sand-mean(d$Sand,na.rm=T))/sd(d$Sand,na.rm = T),
              Site=d$Site,
              MF=(d$MF-mean(d$MF,na.rm=T))/sd(d$MF,na.rm = T),
              SR=(d$SR-mean(d$SR,na.rm=T))/sd(d$SR,na.rm = T),
              Facilitation=(d$Facilitation-mean(d$Facilitation,na.rm=T))/sd(d$Facilitation,na.rm = T),
              Cover=(d$Cover-mean(d$Cover,na.rm=T))/sd(d$Cover,na.rm = T),
              Aridity=(d$Aridity-mean(d$Aridity,na.rm=T))/sd(d$Aridity,na.rm = T),
              Lat=(d$Lat-mean(d$Lat,na.rm=T))/sd(d$Lat,na.rm = T),
              Long_cos=(d$Long_cos-mean(d$Long_cos,na.rm=T))/sd(d$Long_cos,na.rm = T),
              Long_sin=(d$Long_sin-mean(d$Long_sin,na.rm=T))/sd(d$Long_sin,na.rm = T),
              Elevation=(d$Elevation-mean(d$Elevation,na.rm=T))/sd(d$Elevation,na.rm = T),
              Slope=(d$Slope-mean(d$Slope,na.rm=T))/sd(d$Slope,na.rm = T),
              Plot_n=d$Plot_n)
    
    
    for (stability_var in c("rela_dist","Size_tipping","q")){ #for each infered value
      
      model_lmer=lm("value ~ Long_cos + Long_sin + Lat + Slope + Elevation + Sand",
                    data = d2%>%melt(., measure.vars=stability_var)%>%
                      filter(., !is.na(value),!is.na(Facilitation)),
                    na.action = na.omit)
      
      #we remove potential outliers
      rm.outliers = romr.fnc(model_lmer, d2%>%melt(., measure.vars=stability_var)%>%
                               filter(., !is.na(value),!is.na(Facilitation)),
                             trim=2.5)
      d_data_out = rm.outliers$data
      
      #First run lm on covariates
      model_lmer=lm("value ~ Long_cos + Long_sin + Lat + Slope + Elevation + Sand",
                    data = d_data_out%>%filter(., !is.na(value)),
                    na.action = na.omit)
      
      #take the residuals
      resid_model=residuals(model_lmer) #extract residuals
      
      save=d_data_out%>%filter(., !is.na(value),!is.na(Facilitation))%>% #add it to the dataframe 
        add_column(., Resid_mod=resid_model)
        
      
      #first running the SEM with all data
      
      d_sem=save
      all_d=summary(psem(
        lmer(Resid_mod ~ (1|Plot_n) + Cover + MF + Aridity + Facilitation + SR, d_sem),
        lmer(Cover ~ (1|Plot_n) + Aridity + Cover , d_sem),
        lmer(MF ~ (1|Plot_n) + Aridity , d_sem),
        lmer(SR ~ (1|Plot_n) + Aridity , d_sem),
        lmer(Facilitation ~ (1|Plot_n) + Aridity, d_sem)
      ))
      
      
      #SEM with low aridity
      
      d_sem=filter(save,Aridity < 0.7)
      low_arid=summary(psem(
        lmer(Resid_mod ~ (1|Plot_n) + Cover + MF + Aridity + Facilitation + SR, d_sem),
        lmer(Cover ~ (1|Plot_n) + Aridity + Cover , d_sem),
        lmer(MF ~ (1|Plot_n) + Aridity , d_sem),
        lmer(SR ~ (1|Plot_n) + Aridity , d_sem),
        lmer(Facilitation ~ (1|Plot_n) + Aridity, d_sem)
      ))
      
      #SEM with high aridity
      
      d_sem=filter(save,Aridity >= 0.7)
      high_arid=summary(psem(
        lmer(Resid_mod ~ (1|Plot_n) + Cover + MF + Aridity + Facilitation + SR, d_sem),
        lmer(Cover ~ (1|Plot_n) + Aridity + Cover , d_sem),
        lmer(MF ~ (1|Plot_n) + Aridity , d_sem),
        lmer(SR ~ (1|Plot_n) + Aridity , d_sem),
        lmer(Facilitation ~ (1|Plot_n) + Aridity, d_sem)
      ))
      
      
      
    }
    
  }
  
}





# ---------------------- Step 5: Prediction ----

d=tibble();step_size=0.005

for (site in list.files("../Data_new/Prediction/","Dist")){
  
  site_id=as.numeric(gsub(".csv","",strsplit(site,"_")[[1]][3]))
  
  pred=read.table(paste0("../Data_new/Prediction/",site),sep=",")%>%
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
    d_fil=filter(pred,ID_sim==x)
    if (any(d_fil$cover>0)){
      return(d_fil$p[min(which(d_fil$cover !=0))]-step_size)
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
      return(d_fil$cover[which(d_fil$p==d_fil$p[min(which(d_fil$cover !=0))])])
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

write.table(d,"../Data_new/Prediction/Raw_stability_metrics.csv",sep=";")







# ---------------------- Step 6: Climatic projections ----

rm(list=ls())
source("./ABC_drylands_function.R")
library(stars)      # To process the raster data
library(sf)         # To work with vector data
library(ggplot2)    # For plotting
library(patchwork)  # To combine different ggplot plots
library(raster)
library(ncdf4)
library(greenbrown)
library(terra)
library(oceanmap)
for (file_ in list.files("../Data_new/Climatic_data",pattern = ".nc")[-c(1,2)]){
  
  # data_aridity = rast(paste0("../Data_new/Climatic_data/",file_))
  # # years_sample=as.numeric(seq(2006, 2100, by=1))
  # names(data_aridity) = as.character(seq(1950, 2100, by=1)[1:length(names(data_aridity))])
  # 
  # fun=function(x) { 
  #   if (any(is.na(x))){
  #     NA
  #   }else{
  #     lin_mod=lm(x[-1] ~ as.character(seq(1950, 2100, by=1)[1:length(names(data_aridity))]));
  #     summary(lin_mod)$coefficients[2]
  #   }
  # }
  # 
  # site_year=raster::extract(data_aridity,
  #                           d_biocom[,(c("Longitude","Lattitude"))],
  #                           method = "bilinear", buffer = 0)
  # 
  # trends_aridity = apply(site_year,1,fun)
  # write.table(trends_aridity,paste0("../Data_new/Climatic_data/trend_clim_",gsub(".nc","",file_),".csv"),sep=";")
  
  nc_data_aridity = ncvar_get(nc_open(paste0("../Data_new/Climatic_data/",file_)),
                              varid=ifelse(any(grep("aridity",file_)),"aridity_annual-mean","BIO01"))
  
  list_years=list()
  for (i in 1:dim(nc_data_aridity)[3]){
    list_years[[i]]=flip(raster(t(as.matrix(nc_data_aridity[,,i])),xmn=-180, xmx=180, ymn=-90, ymx=90),2)
  }
  aridity_years=raster::stack(list_years)
  
  names(aridity_years) = as.character(seq(1950, 2100, by=1)[1:length(names(aridity_years))])
  
  
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
  
  write.table(trends_aridity,paste0("../Data_new/Climatic_data/trend_clim_",gsub(".nc","",file_),".csv"),sep=";")

}



simulated_sites=sapply(1:length(list.files("../Data_new/Prediction/","Dist")),
            function(x){
              return(as.numeric(gsub(".csv","",strsplit(list.files("../Data_new/Prediction/","Dist")[x],"_")[[1]][3])))})



# Human population data

list_years=list()
index=1
for (file_ in list.files("../Data_new/Climatic_data/Population/SSP3",pattern = ".nc")){
  nc_data_pop = ncvar_get(nc_open(paste0("../Data_new/Climatic_data/Population/SSP3/",file_)),
                              varid=gsub(".nc","",file_))
  list_years[[index]]=flip(raster(t(as.matrix(nc_data_pop)),xmn=-180, xmx=180, ymn=-90, ymx=90),2)
  index=index+1
}
  
pop_years=raster::stack(list_years)
names(pop_years) = as.character(seq(2010, 2100, by=10))


fun=function(x) { 
  if (any(is.na(x))){
    NA
  }else{
    lin_mod=lm(x ~ as.numeric(as.character(seq(2010, 2100, by=10))));
    summary(lin_mod)$coefficients[2]
  }
}

pop_years_data=extract(pop_years, SpatialPoints(d_biocom[,c("Longitude","Lattitude")]))

trends_pop = apply(pop_years_data,1,fun)
  
write.table(trends_aridity,paste0("../Data_new/Climatic_data/trend_clim_",gsub(".nc","",file_),".csv"),sep=";")
  









# ---------------------- Step 7: Validation using simulations ----
## >> 1) Kefi dryland model ----
#First the distance predicted by the kefi model


dist_kefi=tibble()
for (site in list.files("../Data_new/Models_confirmation/Confirm_kefi/Dist_kefi","Dist")){
  
  site_id=as.numeric(gsub(".csv","",strsplit(site,"_")[[1]][3]))
  
  pred=read.table(paste0("../Data_new/Models_confirmation/Confirm_kefi/Dist_kefi/",site),sep=",")%>%
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
for (site in list.files("../Data_new/Models_confirmation/Confirm_kefi/Dist_Eby","Dist")){
  site_id=as.numeric(gsub("Eby","",gsub(".csv","",strsplit(site,"_")[[1]][3])))
  
  pred=read.table(paste0("../Data_new/Models_confirmation/Confirm_kefi/Dist_Eby/",site),sep=",")%>%
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



param_kefi=read.table("../Data_new/Models_confirmation/Confirm_kefi/Parameters_kefi.csv",sep=";")
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
        "../Data_new/Models_confirmation/Confirm_kefi/d_for_figure.rds")



## >> 2) Guichard mussel-bed model ----


#First the distance predicted by the guichard model


dist_guichard=tibble()
for (site in list.files("../Data_new/Models_confirmation/Confirm_Guichard/Dist_guichard","Dist")){
  
  site_id=as.numeric(gsub(".csv","",strsplit(site,"_")[[1]][3]))
  
  pred=read.table(paste0("../Data_new/Models_confirmation/Confirm_Guichard/Dist_guichard/",site),sep=",")%>%
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
for (site in list.files("../Data_new/Models_confirmation/Confirm_Guichard/Dist_Eby","Dist")){
  
  site_id=as.numeric(gsub("Eby","",gsub(".csv","",strsplit(site,"_")[[1]][3])))
  
  pred=read.table(paste0("../Data_new/Models_confirmation/Confirm_Guichard/Dist_Eby/",site),sep=",")%>%
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

param_guichard=read.table("../Data_new/Models_confirmation/Confirm_Guichard/Parameters_guichard.csv",sep=";")
colnames(param_guichard)=c("d","a0","a2")
d_guichard=cbind(d_guichard,param_guichard[d_guichard$Site,])%>%
  filter(., Site %in% unique(d_eby$Site))
d_eby=cbind(d_eby,param_guichard[d_guichard$Site,])

d_spearman=tibble()
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

saveRDS(list(d_spearman=d_spearman,
             d_guichard=d_guichard,
             d_eby=d_eby),
        "../Data_new/Models_confirmation/Confirm_Guichard//d_for_figure.rds")

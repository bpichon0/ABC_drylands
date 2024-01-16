rm(list=ls())
source("../../ABC_drylands_function.R")



n_tot=60
n_sim_kept=100
method_abc="rejection"
  `%!in%` = Negate(`%in%`)  
  data_kefi=read.table("./Stats_kefi.csv",sep=",")
  colnames(data_kefi)=c("r","d","f","m","b","c","delta","rho_p","nb_neigh","clustering","skewness","variance","moran_I","Spectral_ratio","PLR","PL_expo","cv_psd","fmax_psd")
  print(colnames(data_kefi))
  print(dim(data_kefi)) 
  d_sim=read.table("../Simulations.csv",sep=";")%>%
        dplyr::relocate(., Pooling,.after =q )%>%
      #filter(., PL_expo>0,!is.na(PLR))%>%
      filter(., PL_expo>0,!is.na(PLR))
  
  rownames(d_sim)=1:nrow(d_sim)
  
  
  d_param_infer_NN=array(0,c(n_sim_kept,nrow(data_kefi),3))
  d_param_infer_rej=array(0,c(n_sim_kept,nrow(data_kefi),3))
  
  d_NRMSE_sumstat=x_y_stat=tibble()
  

  sumstat_kept=4:14
  
  
  
  for (empirical_id in 1:n_tot){
    print(empirical_id)
    target=data_kefi[empirical_id,8:18]
    if (target$rho_p>0.06){
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
       cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values

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
      
       cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
      x_y_stat=rbind(x_y_stat,target%>%
                       add_column(., Site_ID=empirical_id,Method="rejection",Type="Obs"))
      
      x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)))%>%
                       add_column(.,Site_ID=empirical_id,Method="rejection",Type="Sim"))
      
    }
    
    print(cross_valid$ss)
    
    cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
    
    print(cross_valid$ss)

    mat_sumstat=d_sim[,sumstat_kept]
    
    if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
    
    
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
    
    d_param_infer_rej[,empirical_id,1]=cross_valid$adj.values[,1] # we keep the whole distribution for p
    d_param_infer_rej[,empirical_id,2]=cross_valid$adj.values[,2] # for q
    d_param_infer_rej[,empirical_id,3]=cross_valid$adj.values[,3] # for the scale of observation
    }
    
  }
  
 write.table(d_NRMSE_sumstat,paste0("./NRMSE_sumstat.csv"),sep=";")
 write.table(x_y_stat,paste0("./x_y_stat.csv"),sep=";")
 write.table(d_param_infer_rej,paste0("./param_rej.csv"),sep=";",col.names=F)
  





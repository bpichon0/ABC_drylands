x = c("tidyverse", "ggpubr", "latex2exp", "deSolve", "reshape2", "simecol",
      "abc", "spatialwarnings", "FME","phaseR","ggpattern","LMERConvenienceFunctions",
      "ggquiver", "scales","boot","RColorBrewer","ggnewscale","cluster","pls",
      "MuMIn","png","car","ggtext","Hmisc","lme4","car","spdep","psych",
      "factoextra","FactoMineR","missMDA","GGally","diptest","raster","ape","abctools","viridis","rsq",
      "gradientForest","extendedForest","rfPermute","A3","semEff","piecewiseSEM","ggrepel","mclust"
)

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

#install pacakges if not installed already
# install.packages(setdiff(x, rownames(installed.packages())))
lapply(x, require, character.only = TRUE)


# dir.create("./Data/Simulations",showWarnings = F)
# dir.create("./Data/System_size",showWarnings = F)
# dir.create("./Figures/Final_figs",showWarnings = F)
# dir.create("./Figures/Final_figs/SI",showWarnings = F)

d_biocom=read.table("./Data/data_sites.csv",sep=";")

the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill="white",color="white"),
                                strip.text.y = element_text(size = 10, angle = -90, face = "italic"),
                                strip.text.x = element_text(size = 10, face = "italic"),
                                legend.text = element_text(size = 10))


the_theme2 = theme_classic() + theme(
  legend.position = "bottom",
  panel.border = element_rect(colour = "black", fill=NA),
  strip.background = element_rect(fill = "transparent",color="transparent"),
  strip.text.y = element_text(size = 10, angle = -90),
  strip.text.x = element_text(size = 10),title = element_text(size=8),
  axis.title.y=element_text(size = 10),
  axis.title.x=element_text(size = 10),
  #legend.box="vertical",
  legend.text = element_text(size = 10), text = element_text(family = "NewCenturySchoolbook")
)


title_distance=expression(paste("Distance to the desertification point",italic(" (Dist)")))

# Useful functions ----
my_pal=function(n){
  if (n%%2 ==0){
    return(c(brewer_pal(palette = "BrBG")(n)))
  }else{
    return(c(brewer_pal(palette = "BrBG")(n)[1:round(n/2)],"gray",brewer_pal(palette = "BrBG")(n)[(round(n/2)+2):n]))
  }
}

`%!in%` = Negate(`%in%`)

get_bootstrapped_pval=function(x){
  return(ifelse(length(which(x>0))/length(x)>.5,length(which(x<0))/length(x),length(which(x>0))/length(x)))
}

logit=function(x){return(log(x/(1-x)))}

invlogit=function(x){return(exp(x)/(1+exp(x)))}

get_upper_tri=function(mat){
  mat[lower.tri(mat)]= NA
  diag(mat)=NA
  return(mat)
}

distance=function(sim_data, obs_data) { 
  sum(apply(sim_data - obs_data, 2, function(x) {x^2})) 
}

normalize=function(x){
  ((x-min(x))/(max(x)-min(x)) - 0.5 ) *2
}

# Ploting functions ----

Plot_psd=function(id,best=F){
  print(spatialwarnings::plot_distr(spatialwarnings::patchdistr_sews(Get_empirical_site(id)>0),best_only = best))
}

Plot_psd_raw=function(landscape){
  psd_id=spatialwarnings::patchsizes(landscape>0)
  psd_tibble=tibble(patch_size=psd_id)
  
  psd_tibble$freq=sapply(1:nrow(psd_tibble),function(x){
    return(length(which(psd_tibble$patch_size>=psd_tibble$patch_size[x]))/nrow(psd_tibble))
  })
  print(ggplot(psd_tibble)+
    geom_point(aes(x=patch_size,y=freq))+
    the_theme+
    scale_x_log10()+
    scale_y_log10())
  
}

Plot_dynamics=function(d,different_sim=F){
  
  d$Time=1:nrow(d)
  
  d_melt=melt(d,id.vars=c("Time"))
  print(ggplot(d_melt%>%mutate(.,variable=recode_factor(variable,"Rho_V"="Vegetation","Rho_F"="Fertile","Rho_D"="Degraded")))+
          geom_line(aes(x=Time,y=value,color=variable),lwd=1)+scale_color_manual(values=c("#9DD696","#D2C48F","#777777"))+
          labs(x="Time steps",y="Fraction",color="")+the_theme)
  
}

Plot_landscape=function(landscape,txt="",col_img="black"){
  landscape[landscape>1] = 0
  if (txt==""){
    image(t(apply(landscape,2,rev)),xaxt = "n",yaxt ="n",col=(c("white",col_img) ))
  } else {
    image(t(apply(landscape,2,rev)),xaxt = "n",yaxt ="n",col=(c("white",col_img) ),main=txt)
  }
}

Plot_hist_simu_obs=function(simu, obs){
  par(mfrow=c(3,3))
  par(mar=rep(4,4))
  for (i in 1:ncol(simu)){
    hist(simu[,i],main=colnames(simu)[i])
    abline(v=obs[i],col="red",lwd=3)
  }
  par(mfrow=c(1,1))
}

# Getting images & stat functions ----


Get_sumstat=function(landscape,log_=T){
  
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
  cv_patch=sd(psd$psd_obs)/mean(psd$psd_obs)
  PLR=spatialwarnings::raw_plrange(landscape>0)
  if (nrow(psd$psd_type)==1){ 
    alpha_exp=NA        
  } else {alpha_exp = psd$psd_type$plexpo[which(psd$psd_type$best==T)]} #i.e., when there is no good fit, return NA
  
  if (log_){
    mean_clustering=log(mean_clustering)
    spectral_ratio=log(spectral_ratio)
    max_patchsize=log(max_patchsize/length(landscape))
  }
  d=tibble(rho_p=cover,
           nb_neigh=mean_nb_neigh,clustering=mean_clustering,
           skewness=spatial_ews[2],variance=spatial_ews[1],moran_I=spatial_ews[3],
           Spectral_ratio=spectral_ratio,PLR=PLR,PL_expo=alpha_exp,CV_PSD=cv_patch,
           fmax_PSD=max_patchsize)
  
  return(d)
}

pooling=function(mat, submatrix_size) {
  
  pooling_matrix=matrix(nrow = floor(nrow(mat)/submatrix_size), 
                        ncol = floor(ncol(mat)/submatrix_size), 
                        dimnames = NULL)
  
  for (i in 1:nrow(pooling_matrix)) {
    for (j in 1:ncol(pooling_matrix)) {
      start_row=(i-1) * submatrix_size + 1
      end_row=start_row + submatrix_size - 1
      start_col=(j-1) * submatrix_size + 1
      end_col=start_col + submatrix_size - 1
      
      pooling_matrix[i, j]=mean(mat[start_row:end_row, start_col:end_col])
    }
  }
  return(pooling_matrix)
}

pooling2=function(mat, submatrix_size) {
  
  pooling_matrix=matrix(nrow = floor(nrow(mat)/submatrix_size), 
                        ncol = floor(ncol(mat)/submatrix_size), 
                        dimnames = NULL)
  
  for (i in 1:nrow(pooling_matrix)) {
    for (j in 1:ncol(pooling_matrix)) {
      start_row=(i-1) * submatrix_size + 1
      end_row=start_row + submatrix_size - 1
      start_col=(j-1) * submatrix_size + 1
      end_col=start_col + submatrix_size - 1
      
      pooling_matrix[i, j]=mean(mat[start_row:end_row, start_col:end_col])>=.5
    }
  }
  return(pooling_matrix)
}

inverse_pooling=function(mat, submatrix_size) {
  n=nrow(mat)
  m=ncol(mat)
  new_n=n * submatrix_size
  new_m=m * submatrix_size
  new_mat=matrix(0, nrow = new_n, ncol = new_m)
  
  for (i in 1:n) {
    for (j in 1:m) {
      start_row=(i-1) * submatrix_size + 1
      end_row=start_row + submatrix_size - 1
      start_col=(j-1) * submatrix_size + 1
      end_col=start_col + submatrix_size - 1
      
      new_mat[start_row:end_row, start_col:end_col]=kronecker(mat[i,j], matrix(1, nrow = submatrix_size, ncol = submatrix_size))
    }
  }
  
  return(new_mat)
}

Plot_leveraged_landscape=function(id,res=200){
  landscape=Get_empirical_site(id)
  landscape_pooled=pooling(landscape,nrow(landscape)/res)
  par(mfrow=c(1,2))
  Plot_landscape(landscape,paste0("Resolution = ",nrow(landscape),"�"))
  Plot_landscape(landscape_pooled,paste0("Resolution = ",nrow(landscape_pooled),"�"))
  par(mfrow=c(1,1))
}

Plot_distrib_patches=function(mat){
  print(ggplot(tibble(psd=patchsizes(mat==1)))+
          geom_histogram(aes(x=psd),fill=alpha("blue",.3),color="blue")+
          theme_classic()+
          labs(x="Size of patches",y="Count"))
}

Filtering_small_patches=function(mat,cutoff=30){
  
  raster_mat=raster::raster(mat)
  
  mat_clump=clump(raster_mat) #Get an id for each patch
  mat_clump_veg=as.matrix(mat_clump) #transforming into a matrix of number, each number = patch
  mat_clump_veg[is.na(mat_clump_veg)]=0
  
  table_id_patch=table(as.numeric(mat_clump_veg))
  
  for (i in 1:length(table_id_patch)){ #for each patch
    
    if (table_id_patch[i]<cutoff){
      mat_clump_veg[mat_clump_veg==as.numeric(names(table_id_patch)[i])]=0
    }
  }
  
  mat_clump_veg[mat_clump_veg>1]=1
  return (mat_clump_veg)
}

Boxcox_and_scale=function(d){
  
  
  which_neg=unlist(sapply(1:ncol(d),function(x){
    if(any(d[,x] < 0,na.rm = T)){
      return(x)
    }}))
  
  
  for (x in 1:ncol(d)) if (x %in% which_neg){
    
    b=boxcox(lm(d[,x]+abs(min(d[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
    lambda_x=b$x[which.max(b$y)]
    if (lambda_x !=0){ #to avoid errors
      d[,x] = (exp(d[,x]*(lambda_x)) -1)/(lambda_x)
    }
    
  }else {
    b=boxcox(lm(d[,x] ~ 1),plotit = F,eps = .05)    
    lambda_x=b$x[which.max(b$y)]
    if (lambda_x !=0){ #to avoid errors
      d[,x] = (d[,x]^(lambda_x) -1)/(lambda_x)
    }
  }
  
  
  #Second we scale
  for (x in 1:ncol(d)) d[,x] = (d[,x]-mean(d[,x],na.rm = T))/sd(d[,x],na.rm = T)
  
  return(d)
}

Get_psd=function(id){
  psd_mat=sort(spatialwarnings::patchsizes(Get_empirical_site(id)>0))
  print(psd_mat)
  return(psd_mat)
}

Get_id_in_pdf=function(id){
  return(paste0(floor(id/23)+1," & ",23*2*((id/23) -floor(id/23))))
}

rotate = function(x) t(apply(x, 2, rev))

Get_color_classif=function(n=275){
  rotate=function(x) t(apply(x, 2, rev))
  mm=tcrossprod(seq(1,0,length.out = n))
  tmp1=sapply(col2rgb("purple")/255, function(x) 1-mm*(1-x))
  tmp2=sapply(col2rgb("cyan")/255, function(x) 1-rotate(mm)*(1-x))
  tmp3=sapply(col2rgb("orange")/255, function(x) 1-rotate(rotate(mm))*(1-x))
  tmp4=sapply(col2rgb("red")/255, function(x) 1-rotate(rotate(rotate(mm)))*(1-x))
  
  tmp=(tmp1*tmp2*tmp3*tmp4)
  mat_tmp=matrix(sapply(1:nrow(tmp),function(x){return(rgb(tmp[x,1],tmp[x,2],tmp[x,3]))}),n,n)
  return(mat_tmp)
}

Aggregate_importance=function(importance_mod,n_models="nmod"){
  
  
  if (n_models=="1mod"){
    d=tibble(
      Sand =ifelse("Sand" %in% names(importance_mod),1,0),
      Aridity=ifelse("Aridity" %in% names(importance_mod),1,0),
      Cover=ifelse("Cover" %in% names(importance_mod),1,0),
      Facilitation=ifelse("Facilitation" %in% names(importance_mod),1,0),
      SR=ifelse("SR" %in% names(importance_mod),1,0))
    
  }else {
    d=tibble(
      Sand =ifelse("Sand" %in% names(importance_mod),importance_mod["Sand"],0),
      Aridity=ifelse("Aridity" %in% names(importance_mod),importance_mod["Aridity"],0),
      Cover=ifelse("Cover" %in% names(importance_mod),importance_mod["Cover"],0),
      Facilitation=ifelse("Facilitation" %in% names(importance_mod),importance_mod["Facilitation"],0),
      SR=ifelse("SR" %in% names(importance_mod),importance_mod["SR"],0))
    
  }
  return(d)
}

Get_empirical_site=function(id){
  d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
  return(as.matrix(read.table(paste0("../Data/Data_Biocom/landscapes/",d_biocom$File_ID[id],".txt"))))
}

#ODE of the meanfield model
ODE_MF = function(t,y,param){ 
  p=param[1]
  q=param[2]
  dy=y * (1 - y) * p + (y**2) * (1 - y) * q - (1 - p) * y * (1-y) - (1 - q) * (y**2) 
  return(list(c(dy)))
}

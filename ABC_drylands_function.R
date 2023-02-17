library(simecol)
library(tidyverse);library(reshape2);library(latex2exp)
library(animation);library(magick)
library(deSolve);library(rootSolve);library(ggdendro)
library(FME);library(ggpubr);library(spatialwarnings)
library(reshape2);library(abc);library(igraph);library(cluster)
library(FactoMineR) ;library(factoextra);library(pls)
library(missMDA);library(GGally);library(scales)



the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90, face = "italic"),
                                strip.text.x = element_text(size = 10, face = "italic"),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))



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


fourneighbors = function(landscape, state = 1, bounds = 1) {
  
  neighborhood = matrix(
    c(0, 1, 0,
      1, 0, 1,
      0, 1, 0),
    nrow = 3)
  nb_neighbors = simecol::neighbors(
    x = landscape,
    state = state,
    wdist = neighborhood,
    bounds = bounds)
  
}


Plot_dynamics=function(d,different_sim=F){
  
  d$Time=1:nrow(d)
  
  d_melt=melt(d,id.vars=c("Time"))
  print(ggplot(d_melt%>%mutate(.,variable=recode_factor(variable,"Rho_V"="Vegetation","Rho_F"="Fertile","Rho_D"="Degraded")))+
    geom_line(aes(x=Time,y=value,color=variable),lwd=1)+scale_color_manual(values=c("#9DD696","#D2C48F","#777777"))+
    labs(x="Time steps",y="Fraction",color="")+the_theme)
  
}

Plot_landscape=function(landscape){
  landscape[landscape>1] = 0
  image(landscape,xaxt = "n",yaxt ="n",col=(c("white","black") ))
}









# Patches & landscape structure ----

#From Schneider et al., 2016 TE
Get_patches = function(landscape, state) {
  
  if (state=="+"){
    landscape[landscape %in% c(1,2)] ="+"
  }
  
  pattern = landscape
  pattern = pattern %in% state #keeping the state of interest
  map = rep(NA, times = length(landscape))
  old = rep(99, times = length(landscape)) #to compare
  
  while(!identical(old[pattern], map[pattern])) { 
    
    old = map
    count = as.integer(1)
    
    for(i in which(pattern)) { #for each cell of interest
      
      neighbors = map[x_with_border][x_to_evaluate[i]+interact] #get its neighbors
      if(all(is.na(neighbors)) ) { 
        map[i] = count #then no patch -> = 1
      } else {
        map[i] = min(neighbors, na.rm = TRUE) 
      }
      count = count +1
    }
    
  }
  
  map = as.factor(map)
  patchvec = as.vector(sapply(levels(map), function(i) length(which(map == i) ) )) 
  
  out = vector()
  if(length(patchvec) > 0) out = sort(patchvec) else out = NA
  return(out)
  
} 


mapping = function(width, height, boundary = "periodic", i_matrix = matrix(c(0,1,0,1,NA,1,0,1,0), ncol = 3, byrow = TRUE)) {
  
  X = matrix(as.integer(1:(width*height)), ncol = width, byrow =TRUE)
  X = cbind(X[,width], X, X[,1] )  
  X = rbind(X[height,], X, X[1,] ) 
  x_with_border = as.integer(t(X))
  
  assign("x_with_border", as.integer(t(X))  , envir = .GlobalEnv )
  
  assign("x_to_evaluate", sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]  )	, envir = .GlobalEnv )
  
  
  I = i_matrix	
  
  neighbours_in_I = which(is.finite(abs(I)/abs(I)), arr.in = TRUE)
  
  relrow = neighbours_in_I[,1]-which(is.na(I), arr.ind = TRUE)[1]
  relcol = neighbours_in_I[,2]-which(is.na(I), arr.ind = TRUE)[2]
  
  assign("interact", relrow * dim(X)[2] + relcol, envir = .GlobalEnv )
  
}

count  = function(x, neighbor) {
  
  neighbors = numeric(length = prod(x$dim))
  x_logical_with_border = (x$cells %in% neighbor)[x_with_border]
  for(k in interact) {
    neighbors = neighbors + x_logical_with_border[x_to_evaluate+k]
  }
  return(neighbors)  
}



Get_frequency_number_patches=function(landscape){
  
  if (is.matrix(landscape)){
    landscape=as.vector(landscape)
  }
  
  d_patch=tibble()
  for (i in c(1)){ 
    d_patch=rbind(d_patch,tibble(patch_size=as.numeric(Get_patches(landscape,i)))%>% add_column(., Species=i))
  }

  d_freq_patch=d_size=tibble()
  for (k in unique(d_patch$Species)){  #for each species and for total vegetation
    size_patches_k=filter(d_patch,Species==k)$patch_size
    cumbins=sort(unique(size_patches_k)) 
    d_freq_patch=rbind(d_freq_patch,tibble(Number=as.numeric(sapply(cumbins, function(k) length(which(size_patches_k >= k)) )),
                                           Frequency = as.numeric(sapply(cumbins, function(k) length(which(size_patches_k >= k))/length(size_patches_k) )),
                                           Size=cumbins,
                                           Species=k))
    
    d_size=rbind(d_size,tibble(Max_size=max(cumbins),
                               Species=k))
    
  }
  
  return(list(Patches_size=d_patch,Patches_frequency=d_freq_patch,Max_size=d_size))
  
}



# fitting Power Laws



fitpoly =  function(data , indices, modelout = FALSE) {
  model = lm(log(p) ~  - 1 + log(size) + I(log(size)^2), data = data[indices,] )
  if(modelout) {return(model)} else {return(coefficients(model)[2])} 
} 
fitlm =  function(data , indices, modelout = FALSE) {
  model =lm(log(p) ~  - 1 + log(size), data = data[indices,] )
  if(modelout) {return(model)} else {return(coefficients(model)[2])} 
} 























Get_params=function(delta,b,c,m,d,r,f,z){
  
  return(list(
    
    delta  =  delta ,
    b      =  b     ,
    c      =  c     ,
    m      =  m     ,
    d      =  d     ,
    r      =  r     ,
    f      =  f     ,
    z      =  z
  ))
}

Get_classical_param=function(delta= 0.1,b=0.6,c= 0.3,m=.1,d= 0.2,r=0.0001,f=0.9,z=4){
  return(Get_params(     delta= delta,b=b,c= c,m=m,d= d,r=r,f=f,z=z))
}

Get_initial_lattice=function(frac=c(.8,.1,.1),size=25){
  return(matrix(sample(1:3,replace = T,size=size*size,prob = frac),ncol=size,nrow=size))
}


fourneighbors = function(landscape, state = 1, bounds = 1) {
  
  neighborhood = matrix(
    c(0, 1, 0,
      1, 0, 1,
      0, 1, 0),
    nrow = 3)
  nb_neighbors = simecol::neighbors(
    x = landscape,
    state = state,
    wdist = neighborhood,
    bounds = bounds)
  
}


#for (k in 1:length(params))assign(names(params)[k],params[[k]])

CA_Kefi = function(init, params) {
  
  # Variables : 1 = vegetation, 2 = fertile, 3 = degraded   
  landscape = init
  rho_v = sum(landscape == 1) / length(landscape)
  
  # Neighbors :
  neigh_v = fourneighbors(landscape, state = 1, bounds = 1)
  
  
  with(params, {
    
    colonization = (delta * rho_v + (1 - delta) * neigh_v / z) *(b - c * rho_v )*dt
    
    # calculate regeneration, degradation & mortality rate
    death = m *dt
    regeneration = (r + f * neigh_v / z)*dt
    degradation = d*dt 
    
    # Apply rules
    rnum = runif(length(landscape)) # one random number between 0 and 1 for each cell
    landscape_update = landscape
    
    ## New vegetation
    landscape_update[which(landscape == 2 & rnum <= colonization)] = 1
    
    ## New fertile
    landscape_update[which(landscape == 1 & rnum <= death)] = 2
    landscape_update[which(landscape == 3 & rnum <= regeneration)] = 2
    
    ## New degraded 
    landscape_update[which(landscape == 2 & rnum > colonization & rnum <= colonization + degradation)] = 3
    
    
    #length(which(landscape == 1 & rnum <= death))
    #length(which(landscape == 2 & rnum > colonization & rnum <= colonization + degradation))
    #length(which(landscape == 3 & rnum <= regeneration))
    #length(which(landscape == 2 & rnum <= colonization))
    
    
    return(   list(Rho_v = sum(landscape_update == 1) / length(landscape_update),
                   Rho_f = sum(landscape_update == 2) / length(landscape_update),
                   Rho_D = 1-sum(landscape_update == 1) / length(landscape_update)-sum(landscape_update == 2) / length(landscape_update),
                   Landscape=landscape_update))
  })
  
}


Run_CA_kefi=function(time=seq(1,1000,1),params,ini,plot=F){
  
  
  d=tibble(Time=1,Rho_V=sum(ini == 1) / length(ini),Rho_F=sum(ini == 2) / length(ini),Rho_D=sum(ini == 3) / length(ini))
  state=list(Landscape=ini,Rho_v=d$Rho_V,Rho_f=d$Rho_F,Rho_D=d$Rho_D)
  
  for (k in 2:length(time)){
    
    params$dt=time[k]-time[k-1]
    
    state=CA_Kefi(state$Landscape,params = params)
    if (plot==T & k%%200==0){
      png(paste0("../Figures/State_",k,".png"))
      Plot_landscape(state$Landscape)
      dev.off()
    }
    
    d=rbind(d,tibble(Time=k,Rho_V=state$Rho_v,Rho_F=state$Rho_f,Rho_D=state$Rho_D))
  }
  state$Landscape[state$Landscape>2]=0
  
  return(list(d=d,State=state$Landscape))
  
}

Plot_dynamics=function(d,different_sim=F,simple=F){
  
  if (different_sim==T & simple==F){  
    d_melt=melt(d,id.vars=c("Time","Type"))
    p=ggplot(d_melt%>%mutate(.,variable=recode_factor(variable,"Rho_V"="Vegetation","Rho_F"="Fertile","Rho_D"="Degraded")))+
      geom_line(aes(x=Time,y=value,color=variable,linetype=Type),lwd=1)+scale_color_manual(values=c("#9DD696","#D2C48F","#777777"))+
      labs(x="Time steps",y="Fraction",color="")+the_theme
    return(p)
    
  } 
  
  if (different_sim==F & simple==F){  
    
    d_melt=melt(d,id.vars=c("Time"))
    p=ggplot(d_melt%>%mutate(.,variable=recode_factor(variable,"Rho_V"="Vegetation","Rho_F"="Fertile","Rho_D"="Degraded")))+
      geom_line(aes(x=Time,y=value,color=variable))+scale_color_manual(values=c("#9DD696","#D2C48F","#777777"))+
      labs(x="Time steps",y="Fraction",color="")+the_theme
    return(p)
    
  }
  
  if (simple == T){
    plot(d$Time,d$Rho_V,xlab="Time steps","l",ylab="Fraction",ylim=c(0,1),col="#9DD696")
    lines(d$Time,d$Rho_F,ylim=c(0,1),col="#D2C48F")
    lines(d$Time,d$Rho_D,ylim=c(0,1),col="#777777")
  }
  
}


Plot_empirical=function(id){
  d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
  image(as.matrix(read.table(paste0("../Data/Data_Biocom/landscapes/",d_biocom$File_ID[id],".txt"))),col=c("white","black"),axes=F)
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

Get_empirical_site=function(id){
  d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
  return(as.matrix(read.table(paste0("../Data/Data_Biocom/landscapes/",d_biocom$File_ID[id],".txt"))))
}

Get_path_empirical_raster=function(id){
  
  d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
  land=as.matrix(read.table(paste0("../Data/Data_Biocom/landscapes/",d_biocom$File_ID[id],".txt")))
  
  
  tiff(paste0("./Decreasing_resolution/Tiff_landscapes/Landscape_",id,".tif"), compression = "lzw")
  
  image(land,col=c("white","black"),xaxs="i",yaxs="i")

  dev.off()    
  
  
  return(paste0("./Decreasing_resolution/Tiff_landscapes/Landscape_",id,".tif"))
}


Get_sumstat=function(landscape){
  
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
  
  
  d=tibble(rho_p=cover,
         nb_neigh=mean_nb_neigh,clustering=mean_clustering,
         skewness=spatial_ews[2],variance=spatial_ews[1],moran_I=spatial_ews[3],
         Spectral_ratio=spectral_ratio,PLR=PLR,PL_expo=alpha_exp)
  
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






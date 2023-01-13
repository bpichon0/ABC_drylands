rm(list=ls())
library(simecol)
library(tidyverse);library(reshape2);library(latex2exp)
library(animation);library(magick)
library(deSolve);library(rootSolve)
library(FME);library(ggpubr);library(spatialwarnings)



the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90, face = "italic"),
                                strip.text.x = element_text(size = 10, face = "italic"),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))



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
  image(landscape,xaxt = "n",yaxt ="n",col=c("#9DD696","#D2C48F","#777777") )
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



fitPL = function(psd, p_spanning) {
  
  "
  psd = patch size distribution 
  p_spanning = lower limitnàfnthe up-bent power law
  
  "
  
  # code of fitted classes
  
  
  out = list()
  out$best = NA
  out$AIC = vector("numeric", length = 3)
  out$dAIC = vector("numeric", length = 3)
  
  # criteria for vegetated state & desert state
  
  ##### linear power law model for parameter estimation
  PLlm = lm(I(log(p)) ~  1 - I(log(size)) , data = psd) 
  
  ###########
  
  try( {out$TPLdown = nls(I(log(p)) ~ alpha * log(size) + Sx * (1 - size) , 
                          data = psd,
                          start = list(alpha =  PLlm$coefficients, Sx = 1/1000),
                          #algorithm = "port",
                          trace = FALSE
  )}, silent = TRUE
  )    
  
  if(!is.null(out$TPLdown) & !coefficients(out$TPLdown)["Sx"] <= 0) {
    out$AIC[1] = AIC(out$TPLdown) 
  } else {
    out$TPLdown = list(NA)
    out$AIC[1] = NA
  }
  
  #####
  
  try({out$PL = nls(I(log(p)) ~ alpha * log(size), 
                    data = psd,
                    start = list( alpha =  PLlm$coefficients ),
                    trace = FALSE,
                    nls.control(maxiter = 50)
  )}, silent = TRUE
  )
  
  if(!is.null(out$PL)) {
    out$AIC[2] = AIC(out$PL)
  } else {
    out$PL  = list(NA)
    out$AIC[2] = NA
  }
  
  ###########
  
  
  try({out$TPLup = nls(I(log(p)) ~  log(b) + log(1+(size^(alpha))/b ) , 
                       data = psd,
                       start = list( alpha =  PLlm$coefficients, b = p_spanning ) , 
                       nls.control(maxiter = 50)
  )}, silent = TRUE
  )
  
  
  if(!is.null(out$TPLup)) {
    out$AIC[3] = AIC(out$TPLup) 
  } else { 
    #result$fit$summary$TPLup  = list(NA)
    out$TPLup  = list(NA)
    out$AIC[3] = NA
  }
  
  ###########
  
  out$dAIC =   out$AIC -min(out$AIC, na.rm = TRUE)
  
  out$best = which.min(out$AIC)+1
  
  return(out)
} 




Get_summary_stat=function(landscape){
  
  landscape[landscape<0]=0
  landscape=landscape>0
  # vegetation cover
  Cover=sum(landscape==T)/(dim(landscape)[1]*dim(landscape)[2])
  
  # number of neighbors
  mean_nb_neigh = mean(simecol::neighbors(x = landscape, state = 1, wdist = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3), bounds = 1)[which(landscape==1)])
  
  mapping(dim(landscape)[1],dim(landscape)[2])
  # 
  
  psd=Get_frequency_number_patches(landscape)$Patches_frequency
  colnames(psd)=c("n","p","size","Sp")
  
  p_spawning=tail(psd$p,1)
  
  a=fitPL(psd = psd,p_spanning = p_spawning)
  class = c("DEG", "DOWN","PL", "UP", "COV")[a$best]
  if (class %in% c("DEG","COV")){alpha=NA}
  if (class %in% c("UP")){alpha=coefficients(a$TPLup)["alpha"]}
  if (class %in% c("DOWN")){alpha=coefficients(a$TPLdown)["alpha"]}
  if (class %in% c("PL")){alpha=coefficients(a$PL)["alpha"]}
  
  # generic_sews(landscape)

  return(tibble(Cover=Cover,Neigh=mean_nb_neigh,Slope_psd=alpha))
  
}









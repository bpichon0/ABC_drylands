library(FME);library(ggpubr);library(spatialwarnings)
library(reshape2);library(abc);library(igraph);library(cluster)
library(FactoMineR) ;library(factoextra);library(pls)
library(missMDA);library(GGally);library(scales);library(magick)
library(png);library(EBImage);library(imager)
library(lme4);library(car);library(diptest);library(raster);library(ape)
library(abctools);library(viridis)
library(Hmisc);library(ggtext)


x = c("tidyverse", "ggpubr", "latex2exp", "deSolve", "reshape2", "simecol",
      "abc", "spatialwarnings", "FME","phaseR","ggpattern",
      "ggquiver", "scales","boot","RColorBrewer","ggnewscale","cluster","pls",
      "factoextra","FactoMineR","missMDA","GGally","diptest","raster","ape","abctools","viridis")

#install pacakges if not installed already
install.packages(setdiff(x, rownames(installed.packages())))


x = c("tidyverse", "ggpubr", "latex2exp", "deSolve", "reshape2", "simecol",
      "abc", "spatialwarnings", "FME","phaseR","ggpattern",
      "ggquiver", "scales","boot","RColorBrewer","ggnewscale","cluster","pls",
      "factoextra","FactoMineR","missMDA","GGally","diptest","raster","ape","abctools","viridis")
lapply(x, require, character.only = TRUE)


# dir.create("./Data/Simulations",showWarnings = F)
# dir.create("./Data/System_size",showWarnings = F)
# dir.create("./Figures/Final_figs",showWarnings = F)
# dir.create("./Figures/Final_figs/SI",showWarnings = F)


d_biocom=read.table("../Data_new/biocom_data.csv",sep=";")


the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill="white",color="white"),
                                strip.text.y = element_text(size = 10, angle = -90, face = "italic"),
                                strip.text.x = element_text(size = 10, face = "italic"),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))

my_pal=function(n){
  if (n%%2 ==0){
    return(c(brewer_pal(palette = "BrBG")(n)))
  }else{
    return(c(brewer_pal(palette = "BrBG")(n)[1:round(n/2)],"gray",brewer_pal(palette = "BrBG")(n)[(round(n/2)+2):n]))
  }
}

`%!in%` = Negate(`%in%`)

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

Plot_empirical=function(id,true_landscape=F){
  
  par(mfrow=c(1,1))
  
  if (true_landscape & is.character(id)){
    
    img1=readPNG(paste0("../Data/Data_Biocom/png_img/",id))
    
    grid::grid.raster(img1)
    
  }else {
    if (is.numeric(id)){
      d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
      image(t(apply(as.matrix(read.table(paste0("../Data/Data_Biocom/landscapes/",d_biocom$File_ID[id],".txt"))),2,rev)),col=c("white","black"),axes=F)
      
    }else {
      d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
      image(t(apply(as.matrix(read.table(paste0("../Data/Data_Biocom/landscapes/",d_biocom$File_ID[which(d_biocom$File_ID==gsub(".png","",id))],".txt"))),2,rev)),col=c("white","black"),axes=F)
    }
    
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

Centroid_patches=function(mat){
  
  d=tibble()
  
  raster_mat=raster::raster(mat)
  
  mat_clump=clump(raster_mat) #Get an id for each patch
  mat_clump_veg=as.matrix(mat_clump) #transforming into a matrix of number, each number = patch
  mat_clump_veg[is.na(mat_clump_veg)]=0
  
  table_id_patch=table(as.numeric(mat_clump_veg))[-1] #we only take vegetation, no bare soil
  
  for (i in 2:length(table_id_patch)){ #for each patch, we find the centroid and its size
    
    d=rbind(d,tibble(Size=length(which(as.numeric(mat_clump_veg)==as.numeric(names(table_id_patch)[i]))),ID=i,
                     centroid_x=mean(which(mat_clump_veg==as.numeric(names(table_id_patch)[i]),arr.ind = T)[,2]),
                     centroid_y=mean(which(mat_clump_veg==as.numeric(names(table_id_patch)[i]),arr.ind = T)[,1]),
                     window=nrow(mat)))
  }
  
  return (d)
}

Analyse_image=function(img_id){
  
  d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
  
  if (any(paste0(ifelse(as.numeric(strsplit(d_biocom$File_ID[img_id],"-")[[1]][1])<100,
                        ifelse(as.numeric(strsplit(d_biocom$File_ID[img_id],"-")[[1]][1])<9,paste0("00",d_biocom$File_ID[img_id]),paste0("0",d_biocom$File_ID[img_id])),
                        d_biocom$File_ID[img_id]
                 )
                 ,".png") %in% list.files(paste0("../Data/Data_Biocom/png_img/")))){
    img1=readPNG(paste0("../Data/Data_Biocom/png_img/",
                        ifelse(as.numeric(strsplit(d_biocom$File_ID[img_id],"-")[[1]][1])<100,
                               ifelse(as.numeric(strsplit(d_biocom$File_ID[img_id],"-")[[1]][1])<9,paste0("00",d_biocom$File_ID[img_id]),paste0("0",d_biocom$File_ID[img_id])),
                               d_biocom$File_ID[img_id]
                               )
                        ,".png"))
    
    grid::grid.raster(img1)
    
    
    plot_i=readline("Do we need to divide the data ?") 
    
    if(plot_i=="plot"){
      Plot_empirical(img_id)
      plot_i=readline("Do we need to divide the data ?") 
      dev.off(dev.list()["RStudioGD"])
      
    }
    
    if (plot_i %in% c("no","non","No","Non","n")){
      
      plot_j=readline("Herbs or shrubs ?") 
      
      return(plot_j)
      
    }else if(plot_i=="stop"){
      break
    }else { #multiple functional groups
      
      
      img1_tab = data.frame(expand.grid(x = seq.int(nrow(img1)),
                                        y = seq.int(ncol(img1))),
                            as.data.frame(matrix(img1, ncol = 3)))
      names(img1_tab) = c('x', 'y', 'red', 'green', 'blue')
      
      
      
      df = gather(img1_tab, channel, value, red, green, blue)
      df[ ,'channel'] = factor(df[ ,'channel'], levels = c("red", "green", "blue"),
                               ordered = TRUE)
      
      
      
      
      
      #classifying vegetation (2 types) vs bare soil
      
      km = kmeans(na.omit(img1_tab[ ,c( 'green','blue',"red")]),
                  centers = 3) #3 groups
      
      img1_tab[ ,'clust'] = NA
      img1_tab[!is.na(img1_tab[ ,'red']), 'clust'] = km[['cluster']]
      
      img1_tab[ ,'clust'] = with(img1_tab, as.integer(reorder(as.factor(clust), -green)))
      
      
      print(ggplot(img1_tab) +
              geom_raster(aes(x = y, y = x,
                              fill = as.factor(clust))) +
              coord_fixed() +
              theme_minimal() +
              scale_fill_manual(name = "Vegetation",
                                values = c('#F4EAA4', '#0A8E0B',"lightgreen","red"))+
              scale_y_reverse())
      
      plot_k=readline("Which vegetation needs to be more clustered ? \n (show or again for seing empirical data") 
      
      if (plot_k %in% c("again","show")){
        grid::grid.raster(img1)
        
        print(ggplot(img1_tab) +
                geom_raster(aes(x = y, y = x,
                                fill = as.factor(clust))) +
                coord_fixed() +
                theme_minimal() +
                scale_fill_manual(name = "Vegetation",
                                  values = c('#F4EAA4', '#0A8E0B',"lightgreen","red"))+
                scale_y_reverse())
        
        plot_k=readline("Which vegetation needs to be more clustered ? \n (show or again for seing empirical data") 
        
      }
      
      
      if(plot_k=="stop"){
        break
      }else{
        
        
        clust_filt = with(img1_tab,
                          gblur(matrix(!is.na(clust) & clust == as.numeric(plot_k),
                                       nrow = max(x), ncol = max(y)),
                                sigma = .5)) > .1
        
        img1_tab[,'clust_filt'] = img1_tab[,"clust"]
        
        img1_tab[which(as.vector(clust_filt)) ,'clust_filt'] = as.numeric(plot_k)
        
        print(ggplot(img1_tab) +
                geom_raster(aes(x = y, y = x,
                                fill = as.factor(clust_filt))) +
                coord_fixed() +
                theme_minimal() +
                scale_fill_manual(name = "Vegetation",
                                  values = c('#F4EAA4', '#0A8E0B',"lightgreen","red"))+
                scale_y_reverse())
        
        plot_l=readline("Good enought ? (high or low or else)") 
        
        index_low=index_high=1
        while(plot_l %in% c("high","low","more","less")){
          
          sig=.5
          if (plot_l %in% c("high","less")){
            
            clust_filt = with(img1_tab,
                              gblur(matrix(!is.na(clust) & clust == as.numeric(plot_k),
                                           nrow = max(x), ncol = max(y)),
                                    sigma = sig)) > 0.1*index_high
            
            img1_tab[,'clust_filt'] = img1_tab[,"clust"]
            
            img1_tab[which(as.vector(clust_filt)) ,'clust_filt'] = as.numeric(plot_k)
            print(ggplot(img1_tab) +
                    geom_raster(aes(x = y, y = x,
                                    fill = as.factor(clust_filt))) +
                    coord_fixed() +
                    theme_minimal() +
                    scale_fill_manual(name = "Vegetation",
                                      values = c('#F4EAA4', '#0A8E0B',"lightgreen","red"))+
                    scale_y_reverse())
            
            index_high=index_high+1
            
            
          }else if(plot_l=="stop"){
            break
          } else{
            
            clust_filt = with(img1_tab,
                              gblur(matrix(!is.na(clust) & clust == as.numeric(plot_k),
                                           nrow = max(x), ncol = max(y)),
                                    sigma = sig)) > 0.1/index_low
            
            img1_tab[,'clust_filt'] = img1_tab[,"clust"]
            
            img1_tab[which(as.vector(clust_filt)) ,'clust_filt'] = as.numeric(plot_k)
            
            print(ggplot(img1_tab) +
                    geom_raster(aes(x = y, y = x,
                                    fill = as.factor(clust_filt))) +
                    coord_fixed() +
                    theme_minimal() +
                    scale_fill_manual(name = "Vegetation",
                                      values = c('#F4EAA4', '#0A8E0B',"lightgreen","red"))+
                    scale_y_reverse())
            
            
            index_low=index_low+1
            
          }
          
          if (index_low >2 || index_high>2){
            sig=sig+.5
            index_low=1
            index_high=1
          }
          
          plot_l=readline("Good enought ? (high or low or else)") 
          
          
          
        }
        
        matrix_2=matrix(as.numeric(img1_tab[,'clust_filt']==2),nrow=sqrt(nrow(img1_tab)),ncol=sqrt(nrow(img1_tab)))
        matrix_3=matrix(as.numeric(img1_tab[,'clust_filt']==3),nrow=sqrt(nrow(img1_tab)),ncol=sqrt(nrow(img1_tab)))
        
        write.table(matrix_2,paste0("../Data/Data_Biocom/type_landscape/landscape_type1_",d_biocom$File_ID[img_id],".csv"),sep=";")
        write.table(matrix_3,paste0("../Data/Data_Biocom/type_landscape/landscape_type2_",d_biocom$File_ID[img_id],".csv"),sep=";")
        
        return("Coexistence")
        
      }
      
    }
  }else {
    return(NA)
  }
    
}

ABC_empirical=function(id,model="Eby_feedback"){
  
  #Classic Eby model
  
  d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models.csv",sep=";"),
               Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
               PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
    filter(., rho_p>0.05,Model==model,Pooling=="1")
  
  if (model=="Schneider_aggregation"){
    d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models_with_Schn_aggre.csv",sep=";"),
                 Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
                 PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
      filter(., rho_p>0.05,Model==model,Pooling=="1")
    
  }
  
  d2=tibble()
  Plot_empirical(id)
  target=Get_sumstat(pooling(Get_empirical_site(id),1)>.5)
  
  target[is.na(target)]=0
  
  
  if (any(which(target==0))){
    id_0=which(target==0)
    target=target[-id_0]
    all_sim=d_sim[,c(1:9)[-id_0]]
  }else {
    all_sim=d_sim[,c(1:9)]
  }
  
  test=abc(target=target,
           param = matrix(data = NA,nrow = nrow(all_sim),2),sumstat = all_sim,
           tol = 50/nrow(all_sim),method="rejection")
  
  if (length(target)==8){
    d2=rbind(d2,  rbind(tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))),c(NA,NA)))
  }else {
    d2=rbind(d2,  tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))))
  }
  
  
  target=Get_sumstat(pooling(Get_empirical_site(id),2)>.5)
  target[is.na(target)]=0
  
  
  if (any(which(target==0))){
    id_0=which(target==0)
    target=target[-id_0]
    all_sim=d_sim[,c(1:9)[-id_0]]
  }else {
    all_sim=d_sim[,c(1:9)]
  }
  
  test=abc(target=target,
           param = matrix(data = NA,nrow = nrow(all_sim),2),sumstat = all_sim,
           tol = 50/nrow(all_sim),method="rejection")
  
  if (length(target)==8){
    d2=cbind(d2,  rbind(tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))),c(NA,NA)))
  }else {
    d2=cbind(d2,  tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))))
  }
  
  #Eby model with *2 scaling
  
  
  
  d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models.csv",sep=";"),
               Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
               PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
    filter(., rho_p>0.05,Model==model,Pooling=="1/2")
  
  
  if (model=="Schneider_aggregation"){
    d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models_with_Schn_aggre.csv",sep=";"),
                 Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
                 PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
      filter(., rho_p>0.05,Model==model,Pooling=="1/2")
    
  }
  
  
  target=Get_sumstat(pooling(Get_empirical_site(id),1)>.5)
  target[is.na(target)]=0
  
  if (any(which(target==0))){
    id_0=which(target==0)
    target=target[-id_0]
    all_sim=d_sim[,c(1:9)[-id_0]]
  }else {
    all_sim=d_sim[,c(1:9)]
  }
  
  test=abc(target=target,
           param = matrix(data = NA,nrow = nrow(all_sim),2),sumstat = all_sim,
           tol = 50/nrow(all_sim),method="rejection")
  
  if (length(target)==8){
    d2=cbind(d2,  rbind(tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))),c(NA,NA)))
  }else {
    d2=cbind(d2,  tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))))
  }
  
  
  target=Get_sumstat(pooling(Get_empirical_site(id),2)>.5)
  target[is.na(target)]=0
  
  if (any(which(target==0))){
    id_0=which(target==0)
    target=target[-id_0]
    all_sim=d_sim[,c(1:9)[-id_0]]
  }else {
    all_sim=d_sim[,c(1:9)]
  }
  
  test=abc(target=target,
           param = matrix(data = NA,nrow = nrow(all_sim),2),sumstat = all_sim,
           tol = 50/nrow(all_sim),method="rejection")
  

  if (length(target)==8){
    d2=cbind(d2,  rbind(tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))),c(NA,NA)))
  }else {
    d2=cbind(d2,  tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))))
  }
  
  
  #Eby model with *3 scaling
  
  
  d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models.csv",sep=";"),
               Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
               PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
    filter(., rho_p>0.05,Model==model,Pooling=="1/3")
  
  
  if (model=="Schneider_aggregation"){
    d_sim=mutate(read.table("../Data/Step9_Spatial_resolution/All_sims_models_with_Schn_aggre.csv",sep=";"),
                 Spectral_ratio=log(Spectral_ratio),clustering=log(clustering),
                 PLR=ifelse(is.nan(PLR),0,PLR),PL_expo=ifelse(is.nan(PL_expo),0,PL_expo))%>%
      filter(., rho_p>0.05,Model==model,Pooling=="1/3")
    
  }
  
  
  
  target=Get_sumstat(pooling(Get_empirical_site(id),1)>.5)
  target[is.na(target)]=0
  
  if (any(which(target==0))){
    id_0=which(target==0)
    target=target[-id_0]
    all_sim=d_sim[,c(1:9)[-id_0]]
  }else {
    all_sim=d_sim[,c(1:9)]
  }
  
  test=abc(target=target,
           param = matrix(data = NA,nrow = nrow(all_sim),2),sumstat = all_sim,
           tol = 50/nrow(all_sim),method="rejection")
  
  if (length(target)==8){
    d2=cbind(d2,  rbind(tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))),c(NA,NA)))
  }else {
    d2=cbind(d2,  tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))))
  }
  
  
  target=Get_sumstat(pooling(Get_empirical_site(id),2)>.5)
  target[is.na(target)]=0
  
  if (any(which(target==0))){
    id_0=which(target==0)
    target=target[-id_0]
    all_sim=d_sim[,c(1:9)[-id_0]]
  }else {
    all_sim=d_sim[,c(1:9)]
  }
  
  test=abc(target=target,
           param = matrix(data = NA,nrow = nrow(all_sim),2),sumstat = all_sim,
           tol = 50/nrow(all_sim),method="rejection")
  
  if (length(target)==8){
    d2=cbind(d2,  rbind(tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))),c(NA,NA)))
  }else {
    d2=cbind(d2,  tibble(obs=as.numeric(target),sim=as.numeric(colMeans(test$ss,na.rm = T))))
  }
  
  colnames(d2)=c("Obs_1_Eby_1_OBS","Obs_1_Eby_1_SIM","Obs_2_Eby_1_OBS","Obs_2_Eby_1_SIM","Obs_1_Eby_.5_OBS","Obs_1_Eby_.5_SIM",
              "Obs_2_Eby_.5_OBS","Obs_2_Eby_.5_SIM","Obs_1_Eby_.3_OBS","Obs_3_Eby_.3_SIM","Obs_2_Eby_.3_OBS","Obs_2_Eby_.3_SIM")
  
  na_df=as.data.frame(matrix(NA,3,12))
  colnames(na_df)=colnames(d2)
  d2=rbind(d2,na_df)
  
  
  return(d2)
  
}

Plot_png=function(mat_id){
  
  d_biocom=read.table("../Data/Data_Biocom/biocom_data.csv",sep=";")
  
  img_1=readPNG(paste0("../Data/Data_Biocom/png_img/",
                       ifelse(as.numeric(strsplit(d_biocom$File_ID[mat_id],"-")[[1]][1])<100,
                              ifelse(as.numeric(strsplit(d_biocom$File_ID[mat_id],"-")[[1]][1])<9,paste0("00",d_biocom$File_ID[mat_id]),paste0("0",d_biocom$File_ID[mat_id])),
                              d_biocom$File_ID[mat_id]
                       ),".png"))
  grid::grid.raster(img_1)
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

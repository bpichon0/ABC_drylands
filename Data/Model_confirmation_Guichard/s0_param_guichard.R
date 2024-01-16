
param=expand.grid(d=seq(.01,.3,length.out=10),a0=c(.1,.2,.3),a2=c(.6,.8))

write.table(param,"./Parameters_guichard.csv",row.names = F,col.names = F,sep=";")

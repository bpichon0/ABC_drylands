
param=expand.grid(r=c(.11),d=c(.1),m=.2,c=.1,
                  f=c(.5,.75,.9),b=seq(.4,.9,length.out=20),delta=c(.1))

write.table(param,"./Parameters_kefi.csv",row.names = F,col.names = F,sep=";")

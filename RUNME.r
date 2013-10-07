#!/usr/bin/Rscript --no-init-file

maxiversion='Maxiv1.6'

args=commandArgs(TRUE)
sami=args[1]
guide=args[2]
specst=args[3]
region=args[4]
Ntile=as.numeric(args[5])
saminame=strsplit(sami,'.',fixed=T)[[1]][1]
sami=read.table(paste('Data/',sami,sep=''),header=T)
guide=read.table(paste('Data/',guide,sep=''),header=T)
specst=read.table(paste('Data/',specst,sep=''),header=T)

cat('\n\n\n############       Aaron\'s Easy SAMI Tiling Software       ############\n\n\n')

info=read.table("Params/samiparam.par",header = T,stringsAsFactors=F)
for(i in 1:length(info[1,])){assign(colnames(info)[i],info[info[,'Region']==region,i])}
varnames=read.table("Params/samivarnames.par",header = F,stringsAsFactors=F,col.names=c('CodeVariable','ColName','CatName'))

source('sourceall.rscript')

temp=samitile(sami,guide,specst,region,Ntile,makefiles=T,clean=T,printparams=T,makeplots=T,passnames=T,filenames=c(maxiversion,args[1:3]))

save(temp,file=paste(Loc,'T',StartNtile,'-',StartNtile+length(temp$tileinfo[,1])-1,'/OutFiles/tilingout.r',sep=''))

sami[match(temp$finpri[,1],sami[, varnames[varnames[, 1] == "ids", 2]]),varnames[varnames[,1]=='priority',2]]=temp$finpri[,2]

write.table(sami,file=paste(Loc,'T',StartNtile,'-',StartNtile+length(temp$tileinfo[,1])-1,'/OutFiles/',saminame,'T',StartNtile,'-',StartNtile+length(temp$tileinfo[,1])-1,'.dat',sep=''),quote=T,row.names=F)

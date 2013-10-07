samitile=function(sami,guide,specst,region='g09',Ntile=10,clean=F,makefiles=F,makeplots=F,printparams=F,passnames=T,filenames){

library('fields')
library('sphereplot')

info=read.table("Params/samiparam.par",header = T,stringsAsFactors=F)
for(i in 1:length(info[1,])){assign(colnames(info)[i],info[info[,'Region']==region,i])}
varnames=read.table("Params/samivarnames.par",header = F,stringsAsFactors=F,col.names=c('CodeVariable','ColName','CatName'))

if(printparams){
print(t(info[info[,'Region']==region,]),quote=F)
print(varnames,quote=F)
}

tiles={}
guidetiles={}
tileinfo={}
remids={}

if(makeplots){oldsami=sami}

for(i in StartNtile:(StartNtile+Ntile-1)){
    temp=besttile(sami=sami,guide=guide,specst=specst,region=region)
    tiles=rbind(tiles,cbind(Ntile=i,temp$sami))
    guidetiles=rbind(guidetiles,cbind(Ntile=i,temp$guide))
    sami[match(temp$outpri[,1],sami[, varnames[varnames[, 1] == "ids", 2]]),varnames[varnames[,1]=='priority',2]]=temp$outpri[,2]
    tileinfo=rbind(tileinfo,rbind(c(Ntile=i,temp$Nrem,temp$comp,temp$bestcoord)))
    remids=rbind(remids,cbind(Ntile=i,temp$remids))
    if(temp$comp==1){break}
}

colnames(tileinfo)[2:3]=c(varnames[varnames[,1]=='remtileinfo',2],varnames[varnames[,1]=='comptileinfo',2])
colnames(remids)[2]=varnames[varnames[,1]=='outids',2]

if(clean){
    if(file.exists(paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,sep=''))){
        files=list.files(paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,sep=''),recursive=T)
        invisible(file.remove(paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/',files,sep=''),recursive=T))
    }
}
    
if(makefiles){

dir.create(paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/OutFiles',sep=''), showWarnings = F, recursive = T)
dir.create(paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/Plots',sep=''), showWarnings = F, recursive = T)
dir.create(paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/Tiles',sep=''), showWarnings = F, recursive = T)

    for(i in StartNtile:(StartNtile+length(tileinfo[,1])-1)){
#   make tiles files
        cat('# Targets',Loc,'Y',Year,'S',Sem,'R',Run,'T',formatC(i,width=3,flag='0'),'.fld','\n','# ',tileinfo[tileinfo[,'Ntile']==i,'bestra'],' ',tileinfo[tileinfo[,'Ntile']==i,'bestdec'],'\n','# ',date(),'\n',file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/Tiles/Targets',Loc,'Y',Year,'S',Sem,'R',Run,'T',formatC(i,width=3,flag='0'),'.fld',sep=''),sep='')
        if(passnames){
        cat('#',filenames,'\n',file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/Tiles/Targets',Loc,'Y',Year,'S',Sem,'R',Run,'T',formatC(i,width=3,flag='0'),'.fld',sep=''),sep=' ',append=T)
        }
        suppressWarnings(write.table(rbind(tiles[tiles[,'Ntile']==i,2:length(tiles[1,])]),file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/Tiles/Targets',Loc,'Y',Year,'S',Sem,'R',Run,'T',formatC(i,width=3,flag='0'),'.fld',sep=''),quote=F,row.names=F,append=T))
        
#   make guide files
        cat('# Guide',Loc,'Y',Year,'S',Sem,'R',Run,'T',formatC(i,width=3,flag='0'),'.fld','\n','# ',tileinfo[tileinfo[,'Ntile']==i,'bestra'],' ',tileinfo[tileinfo[,'Ntile']==i,'bestdec'],'\n','# ',date(),'\n',file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/Tiles/Guide',Loc,'Y',Year,'S',Sem,'R',Run,'T',formatC(i,width=3,flag='0'),'.fld',sep=''),sep='')
        if(passnames){
        cat('#',filenames,'\n',file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/Tiles/Guide',Loc,'Y',Year,'S',Sem,'R',Run,'T',formatC(i,width=3,flag='0'),'.fld',sep=''),sep=' ',append=T)
        }
        suppressWarnings(write.table(rbind(guidetiles[guidetiles[,'Ntile']==i,2:length(guidetiles[1,])]),file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/Tiles/Guide',Loc,'Y',Year,'S',Sem,'R',Run,'T',formatC(i,width=3,flag='0'),'.fld',sep=''),quote=F,row.names=F,append=T))
    }

cat('#',date(),'\n',file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/OutFiles/TileInfo',Loc,'Y',Year,'S',Sem,'R',Run,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'.dat',sep=''),sep=' ')
if(passnames){
        cat('#',filenames,'\n',file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/OutFiles/TileInfo',Loc,'Y',Year,'S',Sem,'R',Run,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'.dat',sep=''),sep=' ',append=T)
        }
 suppressWarnings(write.table(tileinfo,file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/OutFiles/TileInfo',Loc,'Y',Year,'S',Sem,'R',Run,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'.dat',sep=''),quote=F,row.names=F,append=T))

cat('#',date(),'\n',file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/OutFiles/RemIDs',Loc,'Y',Year,'S',Sem,'R',Run,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'.dat',sep=''),sep=' ')
if(passnames){
        cat('#',filenames,'\n',file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/OutFiles/RemIDs',Loc,'Y',Year,'S',Sem,'R',Run,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'.dat',sep=''),sep=' ',append=T)
        }
 suppressWarnings(write.table(remids,file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/OutFiles/RemIDs',Loc,'Y',Year,'S',Sem,'R',Run,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'.dat',sep=''),quote=F,row.names=F,append=T))

}

if(makeplots){

    for(i in StartNtile:(StartNtile+length(tileinfo[,1])-1)){
        png(file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/Plots/Tile',Loc,'Y',Year,'S',Sem,'R',Run,'T',formatC(i,width=3,flag='0'),'.png',sep=''),width=800,height=800)
        plottile(tiles=tiles,guidetiles=guidetiles,tileinfo=tileinfo,remids=remids,Ntile=i,region=region)
        dev.off()
    }

png(file=paste(Loc,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'/Plots/Survey',Loc,'Y',Year,'S',Sem,'R',Run,'T',StartNtile,'-',StartNtile+length(tileinfo[,1])-1,'.png',sep=''),width=800,height=800)
plotsurvey(oldsami,tileinfo,remids,plotden=F,stop=PlotStop,region=region)
dev.off()

}

return=list(tiles=tiles,guidetiles=guidetiles,tileinfo=tileinfo,remids=remids,finpri=temp$outpri)
}


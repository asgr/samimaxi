plottile=function(tiles,guidetiles,tileinfo,remids,Ntile=1,region='g09'){

library('magicaxis')
library('plotrix')

info=read.table("Params/samiparam.par",header = T,stringsAsFactors=F)
for(i in 1:length(info[1,])){assign(colnames(info)[i],info[info[,'Region']==region,i])}

lims=c(RAmin,RAmin+RArange,DECmin,DECmin+DECrange)

varnames=read.table("Params/samivarnames.par",header = F,stringsAsFactors=F,col.names=c('CodeVariable','ColName','CatName'))

plotinfo=read.table("Params/samiplot.par",header = F,stringsAsFactors=F)

samipars=1:which.max(plotinfo[,1])
guidepars=(which.max(plotinfo[,1])+1):length(plotinfo[,1])

bestra=tileinfo[tileinfo[,'Ntile']==Ntile,'bestra']
bestdec=tileinfo[tileinfo[,'Ntile']==Ntile,'bestdec']

magplot(bestra,bestdec,xlim=c(bestra+fovradout,bestra-fovradout),ylim=c(bestdec-fovradout,bestdec+fovradout),asp=1,xlab='RA / Deg',ylab='Dec / Deg')

draw.circle(bestra,bestdec,radius=centre2bundleproximity/3600,lty=1)
draw.circle(bestra,bestdec,radius=fovradin,lty=1)
draw.circle(bestra,bestdec,radius=fovradout,lty=2)

tiles=rbind(tiles[tiles[,'Ntile']==Ntile,])
guidetiles=rbind(guidetiles[guidetiles[,'Ntile']==Ntile,])
remids=remids[remids[,'Ntile']==Ntile,]

for(i in 1:length(remids)){
draw.circle(tiles[tiles[,varnames[varnames[,1]=='outids',2]]==remids[i],varnames[varnames[,1]=='outra',2]],tiles[tiles[,varnames[varnames[,1]=='outids',2]]==remids[i],varnames[varnames[,1]=='outdec',2]],radius=bundle2bundleproximity/3600,lty=2)
}

tilesloc=samipars[match(tiles[,varnames[varnames[,1]=='outpriority',2]],plotinfo[samipars,1])]
points(tiles[,c('ra','dec')],col=plotinfo[tilesloc,2],pch=plotinfo[tilesloc,3],cex=plotinfo[tilesloc,4])

guidetilesloc=guidepars[match(guidetiles[,varnames[varnames[,1]=='outpriorityguide',2]],plotinfo[guidepars,1])]
points(guidetiles[,c('ra','dec')],col=plotinfo[guidetilesloc,2],pch=plotinfo[guidetilesloc,3],cex=plotinfo[guidetilesloc,4])

legend('topleft',legend=paste('Tile',Loc,'Y',Year,'S',Sem,'R',Run,'T',formatC(Ntile,width=3,flag='0'),sep=''))
legend('bottomright',legend=c(rev(plotinfo[1:8,5]),rev(plotinfo[guidepars,5]),plotinfo[9,5]),col=c(rev(plotinfo[1:8,2]),rev(plotinfo[guidepars,2]),plotinfo[11,2]),pch=c(rev(plotinfo[1:8,3]),rev(plotinfo[guidepars,3]),plotinfo[9,3]),pt.cex=c(rev(plotinfo[1:8,4]),rev(plotinfo[guidepars,4]),plotinfo[9,4]),bg='white')
}


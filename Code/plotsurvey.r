plotsurvey=function(sami,tileinfo,remids,stop=0.95,region='g09',plotden=F){

library('magicaxis')
library('plotrix')

info=read.table("Params/samiparam.par",header = T,stringsAsFactors=F)
for(i in 1:length(info[1,])){assign(colnames(info)[i],info[info[,'Region']==region,i])}

lims=c(RAmin,RAmin+RArange,DECmin,DECmin+DECrange)

varnames=read.table("Params/samivarnames.par",header = F,stringsAsFactors=F,col.names=c('CodeVariable','ColName','CatName'))

#Make variables nicer to use
ids=sami[,varnames[varnames[,1]=='ids',2]]
ra=sami[,varnames[varnames[,1]=='ra',2]]
dec=sami[,varnames[varnames[,1]=='dec',2]]
survey=sami[,varnames[varnames[,1]=='survey',2]]
priority=sami[,varnames[varnames[,1]=='priority',2]]

#Define some useful logical subsets
allmain=survey >= mainsurvey #Full main target list for the survey- useful for completeness densities etc
allfill=survey >= fillsurvey & survey < mainsurvey #Full filler list for the survey- useful for completeness densities etc
freemain=allmain & priority >= mainpri #Untargeted high priority (main survey) targets
freefill=allfill & priority >= fillpri #Untargeted low priority (filler) targets
regionlogic= ra>=lims[1] & ra<=lims[2] & dec>=lims[3] & dec<=lims[4]

if(plotden){denall=SphHat(ra[allmain],dec[allmain],gridres=c(0.1,0.1),lims=lims,rad=fovradin)}

col=rev(rainbow(80, start = 0, end = 2/3))

layout(cbind(c(1,3,5),c(2,4,5)),width=c(0.8,0.2))

par(mar=c(3.1,3.1,1,1))
begin=SphHat(ra[freemain],dec[freemain],gridres=c(0.1,0.1),lims=lims,rad=fovradin)
if(plotden){begin$z=1-begin$z/denall$z;zlim=c(0,1)}else{begin$z=begin$z/Nmaininbump;zlim=c(0,4)}
image(begin,asp=1,xlim=lims[2:1],ylim=lims[3:4],zlim=zlim,col=col,xlab='',ylab='')
points(ra[freemain],dec[freemain])
magaxis(labels=F,xlab='RA / Deg',ylab='Dec / Deg')
plot.new()
par(mar=c(3.1,1.1,1,1))
if(plotden){gradby=0.2}else{gradby=0.5}
color.legend(0,0,0.5,1,legend=seq(zlim[1],zlim[2],by=gradby),rect.col=col,gradient='y')


par(mar=c(3.1,3.1,1,1))
stoptile=tileinfo[min(c(max(which(tileinfo[,'Comp']<stop))+1,length(tileinfo[,1]))),1]
tarIDs=remids[remids[,1]<=stoptile,2]
tiled=ids %in% tarIDs
if(length(which(! tiled & freemain))>0){
end=SphHat(ra[! tiled & freemain],dec[! tiled & freemain],gridres=c(0.1,0.1),lims=lims,rad=fovradin)
	if(plotden){end$z=1-end$z/denall$z}else{end$z=end$z/Nmaininbump}
	image(end,asp=1,xlim=lims[2:1],ylim=lims[3:4],zlim=zlim,col=col,xlab='',ylab='')
	points(ra[! tiled & freemain],dec[! tiled & freemain])
}else{magplot(ra[! tiled & freemain],dec[! tiled & freemain],xlim=lims[2:1],ylim=lims[3:4])}
for(i in which(tileinfo[,'Ntile']<=stoptile)){
draw.circle(tileinfo[i,'bestra'],tileinfo[i,'bestdec'],radius=fovradin)
text(tileinfo[i,'bestra'],tileinfo[i,'bestdec'],tileinfo[i,'Ntile'],col='brown',cex=2)
}
magaxis(labels=F,xlab='RA / Deg',ylab='Dec / Deg')
plot.new()
par(mar=c(3.1,1.1,1,1))
if(plotden){gradby=0.2}else{gradby=0.5}
color.legend(0,0,0.5,1,legend=seq(zlim[1],zlim[2],by=gradby),rect.col=col,gradient='y')

par(mar=c(3.1,3.1,1,1))
magplot(tileinfo[,c('Ntile',varnames[varnames[,1]=='comptileinfo',2])],ylim=c(0,1),type='l',labels=c(T,T,F,F),xlab='No. Tiles',ylab='Completeness')
TotalTargets=length(which(allmain & regionlogic))
TotalTiled=length(which(tiled & regionlogic))
perfecttiles=TotalTargets/Nmaininbump
abline(tileinfo[1,varnames[varnames[,1]=='comptileinfo',2]]-tileinfo[1,'Ntile']/perfecttiles,1/perfecttiles,lty=3)
abline(h=stop,col='red',lty=2)
if(tileinfo[max(which(tileinfo[,'Ntile']<=stoptile)),varnames[varnames[,1]=='comptileinfo',2]]>=stop){abline(v=stoptile,col='red',lty=2)}
legend('bottomright',legend=c(
paste('Tiles=',stoptile-tileinfo[1,'Ntile']+1),
paste('Comp=',round(length(which(tiled & regionlogic & allmain))/length(which(regionlogic & allmain)),2),'/',round(length(which(tiled & regionlogic & (allmain | allfill)))/length(which(regionlogic & (allmain | allfill))),2) ),
paste('Remaining=',length(which((! tiled & freemain) & regionlogic & allmain)),'/',length(which((! tiled & freefill) & regionlogic & (allmain | allfill)))),
paste('Targets=',length(which(tiled & regionlogic & allmain)),'/',length(which(tiled & regionlogic & (allmain | allfill)))),
paste('Efficiency=',round(length(which(tiled & regionlogic & allmain))/(Nmaininbump*max(which(tileinfo[,'Ntile']<=stoptile))),2),'/',round(length(which(tiled & regionlogic & (allmain | allfill)))/(Nmaininbump*max(which(tileinfo[,'Ntile']<=stoptile))),2))
),
bty='o',bg='white',lty=c(2,2,NA,NA,NA),col=c('red','red','black','black','black'))
legend('topleft',legend=paste('Survey',Loc,'Y',Year,'S',Sem,'R',Run,'T',StartNtile,'-',StartNtile+length(tileinfo[,'Ntile'])-1,sep=''))
}

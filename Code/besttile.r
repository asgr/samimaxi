besttile=function(sami,guide,specst,region='g09'){

#Read in variable
info=read.table("Params/samiparam.par",header = T,stringsAsFactors=F)
for(i in 1:length(info[1,])){assign(colnames(info)[i],info[info[,'Region']==region,i])}

lims=c(RAmin,RAmin+RArange,DECmin,DECmin+DECrange)

varnames=read.table("Params/samivarnames.par",header = F,stringsAsFactors=F,col.names=c('CodeVariable','ColName','CatName'))

#Read in sami data
ids=sami[,varnames[varnames[,1]=='ids',2]]
ra=sami[,varnames[varnames[,1]=='ra',2]]
dec=sami[,varnames[varnames[,1]=='dec',2]]
survey=sami[,varnames[varnames[,1]=='survey',2]]
priority=sami[,varnames[varnames[,1]=='priority',2]]
mag=sami[,varnames[varnames[,1]=='mag',2]]

#Define some useful logical subsets
allmain=survey >= mainsurvey #Full main target list for the survey- useful for completeness densities etc
allfill=survey >= fillsurvey & survey < mainsurvey #Full filler list for the survey- useful for completeness densities etc
freemain=allmain & priority >= mainpri #Untargeted high priority (main survey) targets
freefill=allfill & priority >= fillpri #Untargeted low priority (filler) targets

#Read in guide data
raguide=guide[,varnames[varnames[,1]=='raguide',2]]
decguide=guide[,varnames[varnames[,1]=='decguide',2]]
#Remove guides that are too close to untargeted main survey targets
guide=guide[noclash(raguide,decguide,1:length(guide[,1]),ra[freemain],dec[freemain],ids[freemain],proximity=guide2bundleproximity)$ids1,]
idsguide=guide[,varnames[varnames[,1]=='idsguide',2]]
raguide=guide[,varnames[varnames[,1]=='raguide',2]]
decguide=guide[,varnames[varnames[,1]=='decguide',2]]
magguide=guide[,varnames[varnames[,1]=='magguide',2]]

#Read in standards data
specst=specst[order(specst[,varnames[varnames[,1]=='priorityspecst',2]],decreasing=T),]
#raspecst=specst[,varnames[varnames[,1]=='raspecst',2]]
#decspecst=specst[,varnames[varnames[,1]=='decspecst',2]]
#Remove standards that are too close to untargeted main survey targets
#specst=specst[noclash(raspecst,decspecst,1:length(specst[,1]),ra[freemain],dec[freemain],ids[freemain],proximity=bundle2bundleproximity)$ids1,]
idsspecst=specst[,varnames[varnames[,1]=='idsspecst',2]]
raspecst=specst[,varnames[varnames[,1]=='raspecst',2]]
decspecst=specst[,varnames[varnames[,1]=='decspecst',2]]
priorityspecst=specst[,varnames[varnames[,1]=='priorityspecst',2]]
magspecst=specst[,varnames[varnames[,1]=='magspecst',2]]

xseq=seq(RAmin,RAmin+RArange,by=0.05)
yseq=seq(DECmin,DECmin+DECrange,by=0.05)
posgrid=expand.grid(xseq,yseq)

#Greedy tiling density
if(TilingType=='greedy'){
temp=SphHatGuide(x=ra[freemain],y=dec[freemain],gx=posgrid[,1],gy=posgrid[,2],fovrad=fovradin,lims=lims,buffer=Buffer,multiplier=Multiplier)
}

#Dengreedy tiling density
if(TilingType=='dengreedy'){
tempall=SphHatGuide(x=ra[allmain],y=dec[allmain],gx=posgrid[,1],gy=posgrid[,2],fovrad=fovradin,lims=lims,buffer=Buffer,multiplier=Multiplier)
templeft=SphHatGuide(x=ra[freemain],y=dec[freemain],gx=posgrid[,1],gy=posgrid[,2],fovrad=fovradin,lims=lims,buffer=Buffer,multiplier=Multiplier)
temp=tempall
temp$gz=templeft$gz/tempall$gz
}

#Find the best position for the tile
tempmax=max(as.numeric(temp$gz),na.rm=TRUE)
besttile=which(temp$gz==tempmax,arr.ind=F)
besttile=resample(besttile,1)
bestra=temp$gx[besttile];bestdec=temp$gy[besttile];bestid=idsguide[temp$gloc[besttile]]

#Remove targets too close to the centre of the plate
sami=sami[noclash(ra,dec,1:length(sami[,1]),bestra,bestdec,1,proximity=centre2bundleproximity)$ids1,]
ids=sami[,varnames[varnames[,1]=='ids',2]]
ra=sami[,varnames[varnames[,1]=='ra',2]]
dec=sami[,varnames[varnames[,1]=='dec',2]]
survey=sami[,varnames[varnames[,1]=='survey',2]]
priority=sami[,varnames[varnames[,1]=='priority',2]]
mag=sami[,varnames[varnames[,1]=='mag',2]]

#Redefine some useful logical subsets
allmain=survey >= mainsurvey #Full main target list for the survey- useful for completeness densities etc
allfill=survey >= fillsurvey & survey < mainsurvey #Full filler list for the survey- useful for completeness densities etc
freemain=allmain & priority >= mainpri #Untargeted high priority (main survey) targets
freefill=allfill & priority >= fillpri #Untargeted low priority (filler) targets

#Remove fillers that collide with high priority targets
if(length(which(freefill))>0){
temp=noclash(ra[freefill],dec[freefill],which(freefill),ra[freemain],dec[freemain],1:length(sami[freemain,1]),proximity=bundle2bundleproximity)$ids1
freefill=rep(FALSE,length(freefill))
freefill[temp]=TRUE
}

#Remove guides that are too close to the centre of the plate
guide=guide[noclash(raguide,decguide,1:length(guide[,1]),bestra,bestdec,1,proximity=centre2bundleproximity)$ids1,]
idsguide=guide[,varnames[varnames[,1]=='idsguide',2]]
raguide=guide[,varnames[varnames[,1]=='raguide',2]]
decguide=guide[,varnames[varnames[,1]=='decguide',2]]
magguide=guide[,varnames[varnames[,1]=='magguide',2]]

#Find targets near besttile
allangs=round(sph2car(ra,dec,1,deg=T) %*% t(sph2car(bestra,bestdec,radius=1,deg=T)),10) #Measure all target angles relative to best tile
inner=acos(allangs)*180/pi<=fovradin #Find all objects within the SAMI FoV of the best tile
outer=acos(allangs)*180/pi>fovradin & acos(allangs)*180/pi<=fovradout #Find all objects outside the SAMI FoV of the best tile

#Check we have some targets in the FoV
if(length(inner) == 0){
   stop("No targets in this field of view! Stopping...")
}

#Find standards near besttile
allspecstangs=round(sph2car(raspecst,decspecst,1,deg=T) %*% t(sph2car(bestra,bestdec,radius=1,deg=T)),10) #Measure all target angles relative to best tile
innerspecst=acos(allspecstangs)*180/pi<=fovradin #Find all objects within the SAMI FoV of the best tile
outerspecst=acos(allspecstangs)*180/pi>fovradin & acos(allspecstangs)*180/pi<=fovradout #Find all objects outside the SAMI FoV of the best tile

#Check we have some standard stars in the FoV
if(length(innerspecst) == 0){
   stop("No standard stars in this field of view! Stopping...")
}


#Find guides near besttile
allguideangs=round(sph2car(raguide,decguide,1,deg=T) %*% t(sph2car(bestra,bestdec,radius=1,deg=T)),10) #Measure all guide angles relative to best tile
innerguide=acos(allguideangs)*180/pi<=fovradin #Find all guides within the SAMI FoV of the best tile
outerguide=acos(allguideangs)*180/pi>fovradin & acos(allguideangs)*180/pi<=fovradout #Find all guides outside the SAMI FoV of the best tile

#Check we have some guide stars in this FoV
if(length(innerguide) == 0){
   stop("No guide stars in this field of view! Stopping...")
}

#unpick using all main survey targets inside the SAMI FoV (fovradin) and also the extended region (fovradout)
priorder=unpick(ra[freemain & inner],dec[freemain & inner],ids=ids[freemain & inner],proximity=bundle2bundleproximity,Nsel=Nmainin)
priorderfill=unpick(ra[freefill & inner],dec[freefill & inner],ids=ids[freefill & inner],proximity=bundle2bundleproximity,Nsel=Nfillin)

#Put the tiling target priorities together for the best tile
outtile={}
if(length(which(priorder$Ncollide>0))>0 & Nmainin>0){outtile=cbind(priorder$ids[priorder$Ncollide>0],Pclusmainin)} #Clustered main survey objects within SAMI FoV
if(length(which(priorder$Ncollide==0))>0 & Nmainin>0){outtile=rbind(outtile,cbind(priorder$ids[priorder$Ncollide==0],Punclusmainin))} #Unclustered main survey objects within SAMI FoV
if(length(which(outer & freemain))>0 & Nmainout>0){outtile=rbind(outtile,cbind(resample(ids[outer & freemain],Nmainout),Pmainout))} #Main survey objects outside SAMI FoV
if(length(which(priorderfill$Ncollide>0))>0 & length(which(inner & freefill))>0 & Nfillin>0){outtile=rbind(outtile,cbind(priorderfill$ids[priorderfill$Ncollide>0],Pfillin))} #Clustered fillers within SAMI FoV
if(length(which(priorderfill$Ncollide>0))==0 & length(which(inner & freefill))>0 & Nfillin>0){outtile=rbind(outtile,cbind(priorderfill$ids[priorderfill$Ncollide==0],Pfillin))} #Unclustered fillers within SAMI FoV
if(length(which(outer & freefill))>0 & Nfillout>0){outtile=rbind(outtile,cbind(resample(ids[outer & freefill],Nfillout),Pfillout))} #Fillers outside SAMI FoV
outtile=rbind(outtile[order(outtile[,2],decreasing=T),])
#Next we reduce the priority for main survey inner targets that are beyond Nmaininbump- basically this is so we have 12 ultra high priority targets that we assume in the outer samitile code are the actual targets- this should mostly be the case.
if(length(outtile[,2])>Nmaininbump){
outtile[(Nmaininbump+1):length(outtile[,2]),2]=outtile[(Nmaininbump+1):length(outtile[,2]),2]-1
}
outtile=cbind(outtile,ra[match(outtile[,1],ids)],dec[match(outtile[,1],ids)],mag[match(outtile[,1],ids)],1)
#All potential galaxy IDs output to the Targets list
remids=as.numeric(outtile[,1])
#Now add the nearest standards
if(any(innerspecst)){
outtile=rbind(outtile,cbind(idsspecst[innerspecst],priorityspecst[innerspecst]+10+max(priorityspecst)+1,raspecst[innerspecst],decspecst[innerspecst],magspecst[innerspecst],0))
}
if(any(outerspecst)){
outtile=rbind(outtile,cbind(idsspecst[outerspecst],priorityspecst[outerspecst]+10,raspecst[outerspecst],decspecst[outerspecst],magspecst[outerspecst],0))
}
colnames(outtile)=c(varnames[varnames[,1]=='outids',2],varnames[varnames[,1]=='outpriority',2],varnames[varnames[,1]=='outra',2],varnames[varnames[,1]=='outdec',2],varnames[varnames[,1]=='outmag',2],varnames[varnames[,1]=='outtype',2])

#Put the tiling guide priorities together for the best tile
temp=SphHatGuide(x=ra[ids %in% outtile[outtile[,2]>=Punclusmainin,1]],y=dec[ids %in% outtile[outtile[,2]>=Punclusmainin,1]],gx=raguide,gy=decguide,fovrad=fovradin,lims=lims,buffer=Buffer,multiplier=Multiplier)
goodlength=length(which(ids %in% outtile[outtile[,2]>=Punclusmainin,1]))
outguide={}
if(length(which(temp$gz>=goodlength))>0){outguide=rbind(outguide,cbind(idsguide[temp$gloc[temp$gz>=goodlength]],PGall,goodlength))}
if(length(which(temp$gz>0 & temp$gz<goodlength & innerguide[temp$gloc]))>0){outguide=rbind(outguide,cbind(idsguide[temp$gloc[temp$gz>0 & temp$gz<goodlength & innerguide[temp$gloc]]],PGin,temp$gz[temp$gz>0 & temp$gz<goodlength & innerguide[temp$gloc]]))}
if(length(which(temp$gz>0 & temp$gz<goodlength & outerguide[temp$gloc]))>0){outguide=rbind(outguide,cbind(idsguide[temp$gloc[temp$gz>0 & temp$gz<goodlength & outerguide[temp$gloc]]],PGout,temp$gz[temp$gz>0 & temp$gz<goodlength & outerguide[temp$gloc]]))}
#outguide=outguide[outguide[,1]!=bestid,]
#outguide=rbind(cbind(bestid,PGbest,goodlength),outguide)
outguide=cbind(outguide,raguide[match(outguide[,1],idsguide)],decguide[match(outguide[,1],idsguide)],magguide[match(outguide[,1],idsguide)])
colnames(outguide)=c(varnames[varnames[,1]=='outidsguide',2],varnames[varnames[,1]=='outpriorityguide',2],varnames[varnames[,1]=='outNhipriguide',2],varnames[varnames[,1]=='outraguide',2],varnames[varnames[,1]=='outdecguide',2],varnames[varnames[,1]=='outmagguide',2])

#Generate the observed targets
#The osberved targets will be the top ranked 12 objects that are within the inner fovradin
remids=remids[outtile[,varnames[varnames[,1]=='outpriority',2]] %in% c(Pclusmainin,Pclusmainin-1,Punclusmainin,Punclusmainin-1,Pfillin,Pfillin-1)]
#Here we select the 12 top objects that meet the inner tile criterion
remids=remids[1:min(c(Nmaininbump,length(remids)))]
Nrem=length(remids)
priority[ids %in% remids]=0
freemain=allmain & priority >= mainpri #Untargeted high priority (main survey) targets
freefill=allfill & priority >= fillpri #Untargeted low priority (filler) targets

comp=1-length(which(allmain & freemain & ra>=lims[1] & ra<=lims[2] & dec>=lims[3] & dec<=lims[4]))/length(which(allmain & ra>=lims[1] & ra<=lims[2] & dec>=lims[3] & dec<=lims[4]))

return=list(sami=outtile,guide=outguide,outpri=cbind(ids,priority),comp=comp,remids=remids,Nrem=Nrem,bestcoord=c(bestra=bestra,bestdec=bestdec))
}


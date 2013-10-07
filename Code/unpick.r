unpick=function(ra, dec, ids, proximity = 229, Nsel=17, mean.neighbor=10){
Ncollide={}
if(length(ra)>0){
carttemp = rbind(sph2car(ra,dec,1,deg=T))
matches = fields.rdist.near(carttemp, carttemp, delta = sin((pi/180)*proximity/3600), mean.neighbor = mean.neighbor)
matches$ind=rbind(matches$ind)

keepids={}
remids={}
temp=matches
while(length(keepids)<Nsel & length(keepids)<length(ids) & length(temp$ind)>0){
remids=unique(c(remids,temp$ind[temp$ind[,1] %in% keepids,2],keepids))
temp$ind=rbind(temp$ind[! temp$ind[,1] %in% remids,])
if(length(temp$ind)>0){
    temptab=table(temp$ind[,1])
    temptab=cbind(as.numeric(names(temptab)),as.numeric(temptab))
    keepids=c(keepids,temptab[which.max(temptab[,2]),1])
    Ncollide=c(Ncollide,max(temptab[,2]))
    }
}

ids=ids[keepids]

ids=ids[order(Ncollide,decreasing=T)]
#minus 1 because we include sef matches
Ncollide=Ncollide[order(Ncollide,decreasing=T)]-1
}
return = list(ids=ids,Ncollide=Ncollide)
}


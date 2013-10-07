noclash=function(ra1, dec1, ids1, ra2, dec2, ids2, proximity = 150)
{
require(fields)
require(sphereplot)
carttemp1 = rbind(sph2car(ra1,dec1,1,deg=T))
carttemp2 = rbind(sph2car(ra2,dec2,1,deg=T))
matches = fields.rdist.near(carttemp1, carttemp2, delta = sin((pi/180)*proximity/3600), mean.neighbor = 4)
matches$ind = rbind(matches$ind)

if(length(matches)>0){
ids1=ids1[! 1:length(ids1) %in% matches$ind[,1]]
ids2=ids2[! 1:length(ids2) %in% matches$ind[,2]]
}else{ids1=ids1;ids2=ids2}
return = list(ids1=ids1,ids2=ids2)
}


SphHatGuide=function (x, y, gx, gy, lims=c(129,141,-1,3), fovrad=0.5, buffer=0, density = FALSE, multiplier = 1, orderN=T) 
{
loc=1:length(x)
loc=loc[x>lims[1] & x<lims[2] & y>lims[3] & y<lims[4]]
x=x[loc]
y=y[loc]

gloc=1:length(gx)
gloc=gloc[gx>lims[1]+buffer & gx<lims[2]-buffer & gy>lims[3]+buffer & gy<lims[4]-buffer]
gx=gx[gloc]
gy=gy[gloc]

    require(fields)
    if (any(!is.finite(x)) || any(!is.finite(y))) 
        stop("missing or infinite values in the data are not allowed")
    fullgrid = cbind(gx, gy)
    temp=fields.rdist.near(rbind(sph2car(long=gx,lat=gy,radius=1,deg=T)),rbind(sph2car(long=x,lat=y,radius=1,deg=T)), mean.neighbor = multiplier * 20, delta = sin(fovrad*pi/180))
    tempind=rbind(temp$ind)
    temptab = table(tempind[, 1])
    temptab = cbind(as.numeric(names(temptab)),as.numeric(temptab))
    gz=rep(0,length(gx))
    gz[temptab[,1]]=temptab[,2]
    if(orderN){
    gloc=gloc[order(gz,decreasing=T)]
    gx=gx[order(gz,decreasing=T)]
    gy=gy[order(gz,decreasing=T)]
    gz=gz[order(gz,decreasing=T)]
    }    
    return = list(loc=loc, x=x, y=y, gloc=gloc, gx = gx, gy = gy, gz = gz)
}


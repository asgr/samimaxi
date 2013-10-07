SphHat=function (x, y, rad = 1, gridres = c((max(x) - min(x))/25, (max(y)-min(y))/25), lims = c(range(x), range(y)), density = FALSE, multiplier = 1) 
{
    require(fields)
    nx <- length(x)
    ny <- length(y)
    n = c(1 + (lims[2] - lims[1])/gridres[1], 1 + (lims[4] - 
        lims[3])/gridres[2])
    if (length(y) != nx) 
        stop("data vectors must be the same length")
    if (any(!is.finite(x)) || any(!is.finite(y))) 
        stop("missing or infinite values in the data are not allowed")
    if (any(!is.finite(lims))) 
        stop("only finite values are allowed in 'lims'")
    gx <- seq(lims[1], lims[2], by = gridres[1])
    gy <- seq(lims[3], lims[4], by = gridres[2])
    fullgrid = expand.grid(gx, gy)
    temp = fields.rdist.near(rbind(sph2car(as.matrix(fullgrid),radius=1,deg=T)),rbind(sph2car(long=x,lat=y,radius=1,deg=T)),mean.neighbor=multiplier*ceiling(length(x)*pi*rad^2/((lims[2]-lims[1])*(lims[4]-lims[3]))),delta=sin(rad*pi/180))
    tempind=rbind(temp$ind)
    temptab = table(tempind[, 1])
    pad = rep(0, length(gx) * length(gy))
    pad[as.numeric(names(temptab))] = as.numeric(temptab)
    z <- matrix(pad, length(gx), length(gy))
    if (density) {
        z = z/(nx * pi * rad^2)
    }
    return = list(x = gx, y = gy, z = z)
}


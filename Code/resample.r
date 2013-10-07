resample=function (x, size, ...){
if (length(x)<size){size=length(x)}
if (length(x) <= 1) {
    if (!missing(size) && size == 0) 
        x[FALSE]
    else x
} else sample(x, size, ...)
}


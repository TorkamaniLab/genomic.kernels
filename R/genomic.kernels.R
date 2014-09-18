library("kernlab")

### UTILITY FUNCTIONS ###

gcheck <- function(x){
  if(!is.null(x)){
    if(!is.vector(x)) stop("input must be a vector")
    if(!is.numeric(x)) stop("input must be numeric")
  }
}

loci.weights <- function(x){
  x = as.matrix(x)
  p = colSums(x) / (2*nrow(x))
  1/sqrt(p*(1-p))
}

preprocess <- function(x){
  impute <- function(x){
    if(any(is.na(x))) x[is.na(x)] = median(x, na.rm=T)
    x
  }
  flip <- function(x) if(sum(x) < n) x else abs(2-x)

  n = nrow(x)
  apply(as.matrix(x), 
        2, 
        function(x) flip(impute(x)))
}


### KERNELS ###
weighted.linear.dot <- function(w){
  
  rval <- function(x,y=NULL){
    gcheck(x)
    gcheck(y)

    if(is.null(y)) sum(w*x^2) else 
      sum(w*x*y)   
  }
  
  return(new("weighted.linear.kernel", .Data=rval, kpar=list(w=w)))
}

setClass("weighted.linear.kernel",
         prototype=structure(.Data=function(){},
                            kpar=list()),
         contains=c("kernel"))



state.dot <- function(w){
  rval <- function(x, y=NULL){
    gcheck(x)
    gcheck(y)
    
    if(is.null(y)) sum(w*2) else
      sum(w*(2-abs(x-y)))
  }
  return(new("state.kernel", .Data=rval, kpar=list(w=w)))
}

setClass("state.kernel",
         prototype=structure(.Data=function(){},
                            kpar=list()),
         contains=c("kernel"))

### MATRIX KERNELS ###

kernelMatrix.weighted.linear.kernel <- function(kernel, x, y=NULL){
  
  w = kpar(kernel)$w
  
  k = function(x, y) crossprod(y, w*t(x))
  
  if(is.null(y)) k(x,x) else 
    k(x,y)
}

setMethod("kernelMatrix",
          signature(kernel="weighted.linear.kernel"),
          kernelMatrix.weighted.linear.kernel)

kernelMatrix.state.kernel <- function(kernel, x, y=NULL){
  
  w = kpar(kernel)$w
  
  k = function(x, y){
    x.t = t(x)
    y.t = t(y)
    apply(x.t, 2, function(z) colSums(w*(2-abs(x.t-z))))
  }
  
  if(is.null(y)) k(x,x) else
    k(x,y)

}
setMethod("kernelMatrix",
          signature(kernel="weighted.linear.kernel"),
          kernelMatrix.weighted.linear.kernel)




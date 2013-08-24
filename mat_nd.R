##----------------------------------------##
## The Network Deconvolution algorithm
## LICENSE: GPL-2
##
## AUTHOR: Jitao David Zhang <davidvonpku@gmail.com>
##         2013-08-24
##
## Adapted from the matlab code by Soheil Feizi.
##
## REFERENCES:
##   For more details, see the following paper:
##    Network Deconvolution as a General Method to Distinguish
##    Direct Dependencies over Networks
##    By: Soheil Feizi, Daniel Marbach,  Muriel Médard and Manolis Kellis
##    Nature Biotechnology
##----------------------------------------##

mi <- function(matrix, nbins=10, degree=3) {
  ## TODO: mi not done yet
  warning("mi has not been implemented")
  return(matrix)
}

boundnorm <- function(x) {
  if(min(x) != max(x)) {
    xmin <- min(x, na.rm=TRUE)
    xmax <- max(x, na.rm=TRUE)
    return((x - xmin)/(xmax-xmin))
  } else {
    message("x is a constant matrix\n")
    return(x)
  }
}

netDeconv <- function(mat, beta=0.9, alpha=1, obsOnly=TRUE) {
  if(nrow(mat)!=ncol(mat))
    stop("Input matrix must be a relevance work. If not, use 'mi' to computer mutual informations")

  ## processing the input matrix
  ## linearly mapping the input matrix to be between 0 and 1
  mat <- boundnorm(mat)

  ## diagonal values are filterd
  diag(mat) <- 0

  ## thresholding
  y <- quantile(mat, 1-alpha)
  matTh <- mat * (mat>=y)

  ## making the matrix symmetric if not already
  matTh <- (matTh + t(matTh))/2;

  ## eigen decomposition
  matEig <- eigen(matTh, symmetric=TRUE)
  matEv <- matEig$values
  lamN <- abs(pmin(min(matEv), 0))
  lamP <- abs(pmax(max(matEv), 0))
  m1 <- lamP*(1-beta)/beta
  m2 <- lamN*(1+beta)/beta
  m <- pmax(m1, m2);

  ## network deconvolution
  matEvPrimer <- matEv/(m+matEv);
  matNew1 <- matEig$vectors %*% diag(matEvPrimer) %*% chol2inv(matEig$vectors)

  
  if(obsOnly) {
    indEdges <- matTh>0
    indNonEdges <- matTh==0
    m1 <- max(mat*indNonEdges)
    m2 <- min(matNew1)
    matNew2 <- (matNew1 + pmax(m1-m2, 0)) * indEdges + (mat*indNonEdges)
  } else {
    m2 <- min(matNew1)
    matNew2 <- matNew1 + pmax(-m2, 0);
  }

  ## linearly mapping the deconvolved matrix to be between 0 and 1
  matNd <- boundnorm(matNew2)
}


## NOT TESTED YET WITH THE MATLAB CODE. USE ON YOUR OWN RISK
testData <- matrix(rnorm(10000), nrow=10)
testCor <- cor(t(testData))
testDeconv <- netDeconv(testCor)

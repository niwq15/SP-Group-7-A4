## The 'newt' function has inputs, including theta, func, grad, hess


newt <- function(theta,
                 func,
                 grad,
                 hess=NULL,
                 ...,
                 tol=1e-8,
                 fscale=1,
                 maxit=100,
                 max.half=20,
                 eps=1e-6) {
  ## use the 'func' function and 'theta' to get the objective function
  f <- func(theta,...) 
  
  ## use the 'grad' function and 'theta' to get the gradient vector
  gradval <- grad(theta,...)  
  
  ## check whether the objective and derivatives are not finite
  ## ? f, gradval is single value or a vector
  if (is.finite(f)==FALSE) {
    break
    return("The objective function is not finite at the initial theta value.") 
  }
  
  if (any(is.finite(gradval)==FALSE)) {
    break
    return("The derivatives are not finite at the initial theta value.") 
  }
  
  ## use the 'theta' values as initial x0
  xval <- theta
  len <- length(xval)
  
  ## loop over 2 to max iterations
  for (i in 2:maxit) {
    gradval <- grad(xval)
    
    ## If the objective or derivatives are finite, do the following steps.
    ## If the Hessian matrix is not supplied, we use an approximation by finite
    ## differencing of the gradient vector
    if (is.null(hess)) {
      ## finite difference Hessian
      Hfd <- matrix(0, len,len)
      
      for (i in 1:len) {## loop over parameters
        ## increase th1[i] by eps
        th1 <- xval; th1[i] <- th1[i]+ eps 
        
        ## compute resulting derivative values
        grad1 <- grad(th1,...)
        
        ## approximate -dl/dth[i]
        Hfd[i,] <- (grad1 - gradval)/eps 
      }
      
      ##make sure the matrix is symmetric
      Hfd <- (t(Hfd)+Hfd)/2
      chess <- chol(Hfd)
    } else {
      chess <- chol(hess(xval))
    }

    ## cholesky decomposition to get hessian inverse
    Hi <- backsolve(chess,forwardsolve(t(chess),diag(rep(1,len))))
    
    ## use formula
    xvalcheck <- xval - Hi %*% gradval
    
    ## set k to zero to loop over
    k <- 0
    
    ## set outputs for checks 
    f2 <- func(xval)
    f3 <- func(xvalcheck)
    
    ## while the step fails to reduce the objective reduce the step
    while (f3 > f2) {
      
      ## half step
      xvalcheck <- xval - (1/2)**k * Hi %*% gradval
      
      ## iterate k
      k <- k + 1
      
      ## evaluate new value so the loop terminates at some point
      f3 <- func(xvalcheck)
    }
    
    ## If we have an answer converging within tolerance then stop
    if ( all( abs(gradval) < tol * (abs(f2)+fscale) )) {
      
      ## get final values
      theta <- xvalcheck
      f <- func(theta, ...)
      g <- grad(theta, ...)
      if (is.null(hess)){
        ## finite difference Hessian
        Hfd <- matrix(0, len,len)
        
        for (i in 1:len) {## loop over parameters
          ## increase th1[i] by eps
          th1 <- xvalcheck; th1[i] <- th1[i]+ eps 
          
          ## compute resulting derivative values
          grad1 <- grad(th1,...)
          
          ## approximate -dl/dth[i]
          Hfd[i,] <- (grad1 - gradval)/eps 
        }
        
        ##make sure the matrix is symmetric
        Hfd <- (t(Hfd)+Hfd)/2
        chess <- chol(Hfd)
      } else {
      chess <- chol(hess(theta))
      }
      Hi <- backsolve(chess,forwardsolve(chess,diag(rep(1,len))))
      iter <- i
      
      returnlist <- list(f,theta,iter,g,Hi)
      names(returnlist) <- c("Objective value","x*","Iterations",
                             "Gradient at x*", "Inverse of Hessian at x*")
      ## return values
      return(returnlist)
    } ## else!
    
    ## assign new xval to continue newton method iterations over i
    xval <- xvalcheck
    
  } ## end newton iterations
} ## end newt


# x-squared function, or x^T x for n dimensional x. Works for any n.
# should return 0 vector with a min of 0
xtx <- function(x){
  t(x) %*% x
}

xtxdiv <- function(x){
  2*x
}

xtxhess <- function(x){
  diag(2,length(x))
}


# Simon's example, should return (1,1) with a min of 0
rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}


### function testing ###
k <- 0
theta <- c(0,0)
f <- rb(theta)
gradval <- gb(theta)
chess <- chol(hb(theta))
Hi <- backsolve(chess,forwardsolve(t(chess),diag(rep(1,length(theta)))))
theta2 <- theta - Hi %*% gradval
f2 <- rb(theta)
f3 <- rb(theta2)
while (f3 > f2) {
  
  theta2 <- theta - (1/2)**k * Hi %*% gradval
  k <- k + 1
  f3 <- rb(theta2)
}
f3 - f2

### Hessian inverse test
# create hessian
hess <- array(c(3,0,1,0,3,0,1,0,3),dim=c(3,3))

# save eigenvalues
eigensH <- eigen(hess)$values

## check if hessian is positive definite
# for positive definite all eigenvalues greater than 0
if (any(eigensH <= 0)){ ## if matrix is not positive definite
  
  print("Hessian not positive definite")
  # change this print to a return in the function
  
} else { ## matrix is positive definite
  
  # assign cholesky upper triangle
  chess <- chol(hess)
  
  # assign the hessian inverse (from notes)
  Hi <- backsolve(chess,forwardsolve(t(chess),diag(c(rep(1,length(eigensH))))))
  
  print(Hi)
  # remove this print in the function, result can be checked with 
  # print(solve(hess))
}
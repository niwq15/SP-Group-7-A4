## The 'nullhess' function helps to approximate the Hessian matrix when it is not supplied.
## The 'nullhess' function has inputs, 
## create a function to approximate the hessian when no input is given
nullhess <- function(theta,grad,...,eps=1e-6){
  
  len <- length(theta) ## find length of theta vector
  Hfd <- matrix(0, len,len) ## intialise a zero matrix for finite difference Hessian
  gradval <- grad(theta) ## calculate gradient using the input 'theta'
  
  for (i in 1:len) {## loop over parameters
    th1 <- theta; th1[i] <- th1[i] + eps ## increase i-th theta value by eps
    grad1 <- grad(th1,...) ## compute resulting gradient values
    ## approximate the i-th row in the hessian matrix
    Hfd[i,] <- (grad1 - gradval)/eps 
  }
  
  ## make sure the matrix is symmetric
  Hfd <- (t(Hfd)+Hfd)/2
  return(Hfd)
}



## The 'newt' function has inputs, theta (a vector of initial values), func (the objective function to minimize), grad (the gradient function),
## hess (the Hessian matrix function, 'hess = Null' when hessian function is not given), ... (any arguments of func, grad and hess except theta),
## tol (the convergence tolerance), fscale (an estimate of the magnitude of func near the optimum, which will be used in convergence testing), 
## maxit (the maximum number of Newton iterations to try before giving up), max.half (the maximum number of times a step should be halved to find
## a lower function value before concluding that the step has failed to improve the objective), eps (the finite difference intervals that will be 
## used when Hessian function is not supplied). 

## The 'newt' function 

## The 'newt' function will return f (the value of the objective function at the minimum), theta (the value of the parameters at the minimum), 
## iter (the number of iterations taken to reach the minimum), g (the gradient vector at the minimum), Hi (the inverse of the Hessian matrix at the minimum).

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6) {
  
  f <- func(theta,...) ## use the 'func' function and 'theta' to get the objective function
  gradval <- grad(theta,...) ## use the 'grad' function and 'theta' to get the gradient vector 
  
  ## Check whether the objective and derivatives are not finite at the initial theta
  ## The 'stop' function is used to jump out of the function and print an error message if a problem is detected
  ## whether use 'warning' instead??
  if (is.finite(f)==FALSE) {
    stop("The objective function is not finite at the initial theta value.")
  }
  if (any(is.finite(gradval)==FALSE)) {##note 'gradval' is a vector and need to check each element is finite
    stop("The derivatives are not finite at the initial theta value.")
  }
  
  ## If both the objective and derivatives are finite at the initial theta, do the following steps
  ## Use the 'theta' values as initial x0 
  len <- length(theta)
  
  for (i in 2:maxit) { ## loop over 2 to max iterations
    gradval <- grad(theta,...) # need?
    
    ## If the Hessian matrix is not supplied, we use the 'nullhess' function to calculate an approximation by finite 
    ## differencing of the gradient vector
    if (is.null(hess)) { ## if the 'hess' is not given
      hess <- nullhess ## use the 'nullhess' function instead
    } 
    ## Do we need check whether the hessian matrix is positive definite during iterations??
    chess <- chol(hess(theta,...)) ## solve with cholesky 
    Hi <- backsolve(chess,forwardsolve(t(chess),diag(rep(1,len)))) ## cholesky decomposition to get hessian inverse
    theta2 <- theta - Hi %*% gradval ## use the Newton's formula
    f2 <- func(theta2,...) ## the objective value for theta2
    k <- 0 ## initial value for the number of step halvings

    ## Compare the values of the objective function for the updated theta and the initial theta
    ## Possibility 1: if the step fails to reduce the objective we need to half the step
    while (f2 >= f) { 
      k <- k + 1 ## half step
      theta3 <- theta - ( 0.5 ^(k) )*( Hi %*% gradval) ## the step is halved 
      f2 <- func(theta3,...) ## evaluate the objective function value of the new theta
      ## Stop the halving steps if a smaller value is found
      if (f2 < f) {
        break
        theta2 <- theta3 ## store the better theta value in 'theta'
        n_half <- k
        gradval <- grad(theta,...)
      }
      ## Stop if we have tried max.half step halvings, but still fail to reduce the objective
      if (k == max.half) {
        stop("The step fails to improve the objective after trying 100 halved steps.")
        ## should return the value corresponding to initial theta value
      }
    }
    
    ## Possibility 2: if the step reduces the objective, then we need to test the convergence
    ## The convergence criteria: all elements of the gradient vector have absolute value less than the product of 'tol' and the absolute value of the objective 
    ## function plus fscale
    ## !! important change!!! mentioned in one piazza post @105 
    bound <- tol*abs(f2)+fscale
    if (all(abs(gradval) < rep(bound, len)) ) { ## If we have an answer converging within tolerance then stop
      ## get final values
      theta <- theta2
      f <- func(theta, ...)
      g <- grad(theta, ...)
      H <- hess(theta,...) ## the Hessian matrix
      ## Check whether the Hessian matrix is positive definite at convergence
      eigensH <- eigen(H)$values ## caculate the eigenvalues of 'H'
      if (any(eigensH <= 0)){## if there exists some non-positive eigenvalues
        warning("The Hessian is not positive definite.")
      } ## otherwise the Hessian matrix is positive definite
      chess <- chol(H) ## solve with cholesky
      Hi <- backsolve(chess,forwardsolve(chess,diag(rep(1,len))))
      iter <- i
      returnlist <- list(f,theta,iter,g,Hi)
      names(returnlist) <- c("Objective value","theta*","Number of iterations",
                             "Gradient at theta*", "Inverse of Hessian at theta*")
      return(returnlist) ## return values
    } else {
      ## If the new theta reduces the objective but fails to pass the convergence test, then update the theta and continue the above steps
      theta <- theta2
    }
    
    ## If maxit is reached without convergence
    if (i == maxit) {
      stop("After trying 100 iterations, the convergence is not found.")
    }
    
  } ## end newton iterations
} ## end newt


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

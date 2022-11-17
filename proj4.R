#### Team Members:
# Fatima Kasenally (S2443602); Wenqi Ni (s1792412); Cameron Allan (S1748084)

#### address of github repo:
# 


#### Contributions
#


# Wenqi created the  and added the warnings 

#

#### Overview 
# This file holds the self-contained code to implement Newton’s method for minimization of functions. 
# Newton’s method...
# We firstly create a 'nullhess' function to approximate the Hessian matrix by finite differencing of the gradient vector, which will be used 
# when the Hessian matrix is not supplied. Then we show the 'newt' function. 


## The 'nullhess' function helps to approximate the Hessian matrix when it is not supplied.
## The 'nullhess' function has inputs, 
## create a function to approximate the hessian when no input is given
nullhess <- function(theta,grad,...,eps=1e-6){
  gradval <- grad(theta,...) # calculate gradient using the input 'theta'
  len <- length(gradval) # store gradient vector length
  Hfd <- matrix(0, len,len) ## intialise a zero matrix for finite difference Hessian

  for (i in 1:len) {## loop over parameters
    th1 <- theta #store a copy of theta to operate on 
    th1[i] <- th1[i] + eps  # increase th1[i] by eps
    grad1 <- grad(th1,...) # compute gradient value at shifted theta
    
    #finite difference algebra: 
    # (g(x_i + eps) - g(x_i))  /  ((x_i + eps ) - x_i)
    # calculate the grad w.r.t theta_i 
    Hfd[i,] <- (grad1 - gradval)/eps #ith row of the hessian 
  }
  Hfd <- (t(Hfd)+Hfd)/2 ## make sure the matrix is symmetric
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





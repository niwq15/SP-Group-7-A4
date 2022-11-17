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

## The 'newt' function will apply the Newton's method to minimize the objective function and use Cholesky decomposition to calculate the inverse
## of Hessian matrix. For each iteration, if the updated objective is no less than the initial objective, we will use step halvings and try 
## at most 'max.half' times; if the update theta reduces the objective, we will try a convergence test and if such theta fails to pass it, repeat 
## the steps until we find one with convergence or having tried 'maxit' times. 
## The function will issue warnings for the following cases: 1. if the objective or derivatives are not finite at the initial theta; 2. if the method
## fails to reduce the objective after trying 'max.half' step halvings; 3. if the step fails to find some theta with convergence after trying 'maxit' 
## times; 4. if the Hessian matrix is not positive definite at convergence and during the steps. 

## The 'newt' function will return f (the value of the objective function at the minimum), theta (the value of the parameters at the minimum), 
## iter (the number of iterations taken to reach the minimum), g (the gradient vector at the minimum), Hi (the inverse of the Hessian matrix at the minimum).

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6) {
  
  f <- func(theta,...) #Store the function, evaluated at the first step theta
  gradval <- grad(theta,...) # Store the gradient vector evaluated at the first step theta 
  
  #Check that the function and gradient are continuous (smooth)
  if (is.finite(f)==FALSE) { #check the function at the start point is real 
    #When f is not real the function is not analytic. 
    stop("The objective function is not finite at the starting point. \n 
          The function is not analytic.")
  }
  
  
  if ( any( is.finite( gradval ) == FALSE ) ) { #check the gradient at the start is real 
    stop("The gradient is not finite at the starting point. \n 
          The function is not analytic.")
  }
  #If the Hessian is not supplied, calculate it through finite difference methods 
  if (is.null(hess)) { #no hessian supplied 
    # calculate the hessian using the finite difference method 
    h <- nullhess(theta,grad,...,eps=1e-6)
  } else { 
    #if the hessian matrix is supplied then use it to calculate the matrix
    h <- hess(theta,...)
  }
 
  #Minimization requires the hessian be positive (semi) definite 
  options(show.error.messages = TRUE)
  #Cholesky decomposition works only if matrix is positive definite 
  #try Cholesky decomposition, if it fails console returns an error 
  R <- try(chol(h) , stop("Hessian is not positive semi-definite at the minimum", call. = FALSE))
  #Hessian Inverse hi - using Cholesky (only works if chol works)
  hi <- chol2inv(chol(h))

  #Implement Newton's Method: 
  # theta_(k+1) = theta_k - f''(theta_k)^(-1) * f'(theta_k)
  # f''(theta_k)^(-1) = hi @ k th step 
  # f'(theta_k) = gradval @ k th step 
  # f(theta_k) = f @ k th step 
  # theta_k = k th param val scanned
  # theta_(k+1) = next theta point to scan (should be closer to the min)
    
  n_iter <- 1 # count the number of iterations 
  while (n_iter <= maxit ){ # iterate until maxit limit exceeded 
   
    theta2 <-  theta - hi %*% gradval #find the next step 
    f2 <- func(theta2,...)                #evaluate the fn at next step 
    # if the step increases the fn value rather than decreases, 
    # half the step size in the same direction (overstepped the min)
    n_half <- 0  # counter for times steps halved 
    while ( f2 >= f ){ #iterate until fn val of next step is lower than current step 
      n_half <- n_half + 1 #counter goes up 
      
      #if we have hit the halving limit without finding the minimum
      if (n_half >= max.half ) { 
        #submit error warning to the console
        stop( paste( "The next step cannot be found. \n 
                     The step size has been halved ", max.half, " times.")
              , call. = FALSE)
      } #else continue halving the step size
      
      #half the step size (again)
      theta2 <-  theta - ( hi %*% gradval ) * ( 0.5 ^(n_half) )
      #evaluate the function at the next step 
      f2 <- func(theta2,...) 
    }#end while (step halving)
   
    #Update the variables using the latest valid theta 
    theta <-  theta2
    f <- f2
    
    gradval <- grad(theta2,...) #next step's gradient vector 
    
    #next step's hessian matrix:
    if (is.null(hess)) { #no analytic hessian fn supplied  
      # calculate the hessian of the next step using the finite difference method 
      h <- nullhess(theta2,grad,...,eps=1e-6)
    } else { 
      #if the hessian matrix is supplied use it to calculate the hessian matrix
      #of the next step 
      h <- hess(theta2, ...)
    }
  
    #test that the hessian is still positive definite 
    R <- try(chol(h) , stop(paste("Hessian is not positive semi-definite at step ", n_iter), call. = FALSE))
    #Calculate the Hessian Inverse 'hi' - using Cholesky (only works if chol works)
    hi <- chol2inv(chol(h)
    
    #Test for convergence
    #if the gradient values are all zero (according to our convergence condition)
    if (all( abs(gradval) <= tol * (abs(f) + fscale ) ) ) { 
      #function returns the following parameters
      output <- list ( f = f, #function value at the minimum
                       theta = theta, #location of minimum
                       iter = n_iter, #number of iterations required to reach min
                       g = gradval, #gradient vector at the minimum 
                       Hi = hi) #inverse hessian at the minimum 
      return(output)
      
    }   
    n_iter <- n_iter + 1 #iteration num count increases
  }
  #If the code gets to this part number of iterations has exceeded limit :(
  warning("The maximum number of iterations has been exceeded")
  #function returns the following parameters at the last theta position searched
  output <- list ( f = f, #function value at the minimum
                   theta = theta, #location of minimum
                   iter = n_iter, #number of iterations required to reach min
                   g = gradval, #gradient vector at the minimum 
                   Hi = hi) #inverse hessian at the minimum 
  return(output)
  
}#end function 

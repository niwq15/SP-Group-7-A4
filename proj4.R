#### Team Members:
# Fatima Kasenally (S2443602); Wenqi Ni (s1792412); Cameron Allan (S1748084)

#### address of github repo:
# https://github.com/niwq15/SP-Group-7-A4


#### Contributions
#Fatima: added the finite test for the function and grad
# created a function to calculate the hessian via finite difference methods 
# created a test for hessian positive definiteness 
# implemented newton's method 
# implemented iteration limit (with error warning)
# implemented half step search (with error warning)
# feasibility checks of next step (f(theta_2), g(theta_2), h(theta_2) validation)
# calculated the inverse hessian (when applicable)
# implemented the convergence test 
# testing and validation of the code 

# Cameron:

# Wenqi: 
# created the function to approximate the Hessian matrix using finite difference method
# implemented the ways to perturb the Hessian matrix if it is not positive definite
# implemented newton's method and half step
# added most of the warnings, including checking whether the objective function
# and gradient are finite, the maximum number of halving steps, the maximum number of
# iterations. 

###### Overview 
# This file holds the self-contained code to implement Newton’s method for 
# minimization of functions. 
# The newton's method is based on the minimizing successive quadratic approximation
# For an objective function, its gradient vector and Hessian matrix will be used 
# to implement the method. The Hessian matrix needs to be positive definite, and 
# if not, we need to perturb it to be so. 
# First we create a 'nullhess' function to approximate the Hessian matrix via 
# finite difference methods using the gradient vector, which will be used 
# when the Hessian matrix is not supplied. 
# Then we show the 'newt' function, which implements Newton's method for 
# optimization. In the 'newt' function, some convergence test will be adopted. 


## The 'nullhess' function approximates the Hessian matrix when its analytic
## form is not supplied. It does this using the finite difference method. 
## The 'nullhess' function has inputs: 
    # theta : the start point of our search 
    # grad : the gradient vector at the start point 
    # eps : the size of the pertubation in theta 
## The output is the hessian matrix defined at theta
nullhess <- function(theta,grad,...,eps=1e-6){ 
  gradval <- grad(theta,...) # calculate gradient vector at 'theta'
  len <- length(gradval) # store gradient vector length
  Hfd <- matrix(0, len,len) # intialise a zero matrix for finite difference Hessian
  
  for (i in 1:len) {## loop over parameters
    th1 <- theta #store a copy of theta to operate on 
    th1[i] <- th1[i] + eps  # increase th1[i] by eps
    grad1 <- grad(th1,...) # compute gradient value at shifted theta
    
    #finite difference algebra: 
    # (g(theta_i + eps) - g(theta_i))  /  ((theta_i + eps ) - theta_i)
    # calculates the gradient's derivative w.r.t theta_i 
    Hfd[i,] <- (grad1 - gradval)/eps #ith row of the hessian 
  }
  Hfd <- (t(Hfd)+Hfd)/2 ## make sure the matrix is symmetric
  return(Hfd)
}


## The 'newt' function has inputs:
    # theta (a vector of initial values),
    # func (the objective function to minimize), 
    # grad (the gradient function),
    # hess (the Hessian matrix function, 'hess = Null' by default when hessian 
    #   function is not given),
    # ... (any additional arguments of func, grad and hess except theta),
    # tol (controls the convergence tolerance),
    # fscale (an estimate of the magnitude of func near the optimum,
    #   which will be used in convergence testing), 
    # maxit (the maximum number of Newton iterations to try before giving up),
    # max.half (the maximum number of times a step should be halved to find
    #    a lower function value before concluding that the step has failed to 
    #    improve the objective), 
    # eps (the finite difference intervals that will be used when Hessian 
    #   function is not supplied). 

## The 'newt' function will apply the Newton's method to minimize the objective 
## function and use Cholesky decomposition to calculate the inverse of Hessian
## matrix. For each iteration, if the updated objective is greater than the 
## initial objective we will half the step size to seach for a smaller objective
## value, this will be done at most 'max.half' times; if the updated theta 
## reduces the objective, we will try a convergence test. If the convergence
## test is failed we repeat the procedure until we find a new theta that 
## satifies the convergence requirements or the procedure has iterated 'maxit' 
## times. 

##The function will issue warnings for the following cases:
      ## 1. if the objective or derivatives are not finite at the initial theta;
      ## 2. if the method fails to reduce the objective after trying 'max.half' 
      ##    step halvings; 
      ## 3. if the step fails to find some theta that satifies convergence after
      ##    iterating newton's method 'maxit' times; 
      ## 4. if the Hessian matrix is not positive definite at convergence and 
      ##    during the steps. 

## The 'newt' function will return :
    # f (the value of the objective function at the minimum), 
    # theta (the value of the parameters at the minimum), 
    # iter (the number of iterations taken to reach the minimum), 
    # g (the gradient vector at the minimum), 
    # Hi (the inverse of the Hessian matrix at the minimum).

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,
                 max.half=20,eps=1e-6) {
  
  f <- func(theta,...) # store the function, evaluated at the first step theta
  gradval <- grad(theta,...) # store the gradient vector evaluated at the first
                              # step theta
  
  # Check the objective and derivatives are real at the initial theta
  # and therefore that the function is continuous.
  
  # The 'stop' function ends the method and returns an error message if there 
  # is a problem detected
  if (is.finite(f)==FALSE) {# if the function is not finite at the start point
    # When f is not real the function is not analytic. 
    stop("The objective function is not finite at the starting point. \n 
          The function is not analytic.")
  }  
  if ( any( is.finite( gradval ) == FALSE ) ) { # check the gradient at the start
                                                # is real 
    stop("The gradient is not finite at the starting point. \n 
          The function is not analytic.")
  }
  
  # If both the objective and derivatives are finite at the initial theta, 
  # proceed: 
  
  # Hessian validation:
  # If the Hessian is not supplied, calculate it via finite difference methods 
  if (is.null(hess)) { #no hessian supplied 
    # calculate the hessian using the finite difference method 
    h <- nullhess(theta,grad,...,eps=1e-6)
  } else { 
    #if the hessian matrix is supplied then use it to calculate the matrix
    h <- hess(theta,...)
  }
  
  # Minimization requires the hessian be positive (semi) definite :
  options(show.error.messages = TRUE) #console displays error messages 
  # Cholesky decomposition works only if matrix is positive definite 
  # Try Cholesky decomposition, the Hessian is perturbed. 

  ## The following part deals with the case that the Hessian matrix is not 
  ## positive definite and in this case we perturb the Hessian matrix to be so
  ## The approach is to add a multiple of the identity matrix and its matrix norm
  ## to it, trying till we get a positive definite Hessian matrix (tested by 
  ## Cholesky or eigen decomposition)
  R <- try(chol(h), silent = TRUE)
  if (inherits(R, "try-error")) {
    i <- 1
    while (i < maxit) { ## also iterate at most 100 times
      h1 <- h + (10^(-7+i))*norm(h)*diag(length(gradval)) 
      eigensH <- eigen(h1)$values
      if (all(eigensH > 0)){
          h <- h1
        break
      }
      i <- i + 1
    }
  }
  
  ## Hessian Inverse hi - using Cholesky (only works if chol works)
  # if the positive definite hessian was not found in 100 pertubations stop 
  # searching 
  hi <- try(chol2inv(chol(h)), stop("A positive definite Hessian was not found"),
            call. = FALSE)
  
  
  ## Implement Newton's Method: 
  ## Algrbraic Formulation of the kth iteration: 
  ## theta_(k+1) = theta_k - f''(theta_k)^(-1) * f'(theta_k)
  ## f''(theta_k)^(-1) = hi @ k th step 
  ## f'(theta_k) = gradval @ k th step 
  ## f(theta_k) = f @ k th step 
  ## theta_k = k th param val scanned
  ## theta_(k+1) = next theta point to scan (should be closer to the min)
  
  n_iter <- 1 # counter for the number of iterations 
  while (n_iter <= maxit ){ # try at most 'maxit' iterations before giving up
    
    theta2 <-  theta - hi %*% gradval # find the next step 
    f2 <- func(theta2,...)  # evaluate the fn at next step 
    
    ## Compare the values of the objective function for the updated theta and 
    ## the initial theta.
    
    ## Case 1: if the step fails to reduce the objective we need to half the 
    ## step size in the same direction (overstepped the min)
    n_half <- 0  # counter for times steps halved 
    while ( f2 >= f ){# iterate until fn val of next step (f2) is lower than 
                      # that of the current step (f)
      n_half <- n_half + 1 # counter goes up 
      # Stop if we have tried max.half step halvings, and failed to reduce the 
      # objective
      if (n_half >= max.half ) { 
        # submit error warning to the console
        stop( paste( "The next step cannot be found. \n 
                     The step size has been halved ", max.half, " times.")
              , call. = FALSE)
      } #else continue halving the step size
      
      # Halving the step size (again)
      theta2 <-  theta - ( hi %*% gradval ) * ( 0.5 ^(n_half) )
      f2 <- func(theta2,...) # evaluate the objective at the new theta 
    }# end while (step halving)
    
    ## Check half stepping does not lead to a non-finite objective function value
    if (is.finite(f2)==FALSE) { #check if the function is real 
      #When f2 is not real the function is not analytic at theta2. 
      stop("The function value at the new step is not real. \n 
          The function here is not analytic.")
    }
    
    ## Case 2: if the step reduces the objective, then we test for convergence
    ## Update the variables using the latest valid theta 
    theta <-  theta2 # update theta 
    f <- f2 #update objective function
    gradval <- grad(theta2,...) #next step's gradient vector 
    
    ## Next step's hessian matrix:
    if (is.null(hess)) { # no analytic hessian fn supplied  
      # Calculate hessian of the next step using the finite difference method 
      h <- nullhess(theta2,grad,...,eps=1e-6)
    } else { 
      # If the hessian matrix function is supplied use it to calculate 
      # the hessian matrix of the next step 
      h <- hess(theta2, ...)
    }
    
    ## Test whether the hessian is still positive definite (same method as before)
    R <- try(chol(h), silent = TRUE)
    if (inherits(R, "try-error")) {
      i <- 1
      while (i < maxit) {
        h1 <- h + (10^(-7+i))*norm(h)*diag(length(gradval)) 
        eigensH <- eigen(h1)$values
        if (all(eigensH > 0)){
            h <- h1
          break
        }
        i <- i + 1
      }
    }
    # Calculate the Hessian Inverse 'hi' using Cholesky 
    # (only works if chol works)
    hi <- try(chol2inv(chol(h)), stop("A positive definite Hessian was not found"),
              call. = FALSE)
    
    ## Test for convergence : 
    # If the gradient values are all approximately zero 
    # (according to our convergence condition)
    if (all( abs(gradval) < tol * (abs(f) + fscale ) ) ) { #if there is 
      # convergence then the minimum of the function has been found. 
      # function returns the following parameters:
      output <- list ( f = f, # function value at the minimum
                       theta = theta, # location of minimum
                       iter = n_iter, # number of iterations to reach min
                       g = gradval, # gradient vector at the minimum 
                       Hi = hi) # inverse hessian at the minimum 
      return(output) #function ends
    }   
    
    ## If the new theta reduces the objective but fails to pass the convergence 
    ## test, repeat the procedure above

    n_iter <- n_iter + 1 # iteration counter increases
  }
  # If this part of the method is reached, the number of iterations
  # has exceeded the limit (maxit)
  
  warning("The maximum number of iterations has been exceeded.")
  
  # The function returns the following parameters at the last theta position
  # searched: 
  output <- list ( f = f, # function value at the last theta
                   theta = theta, # value of the parameters at the last theta
                   iter = n_iter, # number of iterations executed
                   g = gradval, # gradient vector at the last theta
                   Hi = hi) # inverse hessian at the last theta  
  return(output) # return the list 
  
}#end function 

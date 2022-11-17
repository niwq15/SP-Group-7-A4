#For convenience define a function to calculate the hessian via finite difference methods 
#used in multiple places in the code 

nullhess <- function(theta,grad,...,eps=1e-6){
  gradval <- grad(theta,...) # calculate gradient using the input 'theta'
  len <- length(gradval) # store gradient vector length
  Hfd <- matrix(0, len,len) # intialise zeros matrix the size of the Hessian
  
  for (i in 1:len) {# calculate deriv w.r.t each param 
  
    th1 <- theta #store a copy of theta to operate on 
    th1[i] <- th1[i] + eps  # increase th1[i] by eps
    grad1 <- grad(th1,...) # compute gradient value at shifted theta
    
    #finite difference algebra: 
    # (g(x_i + eps) - g(x_i))  /  ((x_i + eps ) - x_i)
    # calculate the grad w.r.t theta_i 
    Hfd[i,] <- (grad1 - gradval)/eps #ith row of the hessian 
  }
  return(Hfd) #output the finite diff hessian 
}



# The 'newt' function has inputs, including theta, func, grad, hess


newt <- function(theta, #this is our initial search point 
                 func, #this is the functional form of the function
                 grad, #this is the functional form of the gradient 
                 hess=NULL, #this is the functional form of the hessian, 
                            #if no input is given this defaults to NULL 
                 ..., #allows newt to absorb the additional parameters required by 
                      #func, grad and hess
                 tol=1e-8, #tolerance controls severity of convergence requirements
                 fscale=1, #fscale is the order of the function val around the min
                 maxit=100, #maxit is the maximum number of iterations of 
                            #newton's method executed
                 max.half=20, #the maximum number of times we half our step size
                              #while searching for a next step that returns a lower 
                              #function value
                 eps=1e-6) {  #epsilon is the step in theta_i used to calculate 
                              #the derivative w.r.t. theta_i using finite differences

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
  
  while (n_iter < maxit ){ # iterate until maxit limit reached 
    
    theta2 <-  theta - hi %*% gradval #find the next step 
    
    f2 <- func(theta2,...)                #evaluate the fn at next step 
    
    
    # if the step increases the fn value rather than decreases, 
    # half the step size in the same direction (overstepped the min)
    n_half <- 0  # counter for times steps halved 
    while ( f2 > f ){ #iterate until fn val of next step is lower than current step 
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
    hi <- chol2inv(chol(h))
    
    
    
    
    
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
  
  
  
  
  
  
  
  ### warning is max number of iteration exceeded 
  
  
  
  
  
  
}





#Paste at the function end: 

#check hessian is positive semi-definite at the minimum.
#Code: 
options(show.error.messages = TRUE, call. = FALSE)
R <- try(chol(h) , stop("Hessian is not positive semi-definite at the minimum"))

#Hessian Inverse hi - using cholesky (only works if chol works)
hi <- chol2inv(chol(h))











#########  TESTING ZONE ###########



#Using TRY and STOP function to run the +ive hessian condition test 

options(show.error.messages = FALSE)
try(log("a"))
print(.Last.value)
options(show.error.messages = TRUE)


myfunc <- function(i){
  options(show.error.messages = TRUE)
  r <- try(log("a"), stop("try failed", call. = FALSE))
  #if (inherits(r ,"try-error")){
  #  stop()
  #}
  i<-i+1
  cat(i)
}

a <- myfunc(3)





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

th <- c(0,0)

rb(th)
gb(th)
hb(th)


# First error works correctly 

if (is.finite(Inf)==FALSE) { #check the function at the start point is real 
  #When f is not real the function is not analytic. 
  stop("The objective function is not finite at the starting point. \n 
          The function is not analytic.")
}

if (is.finite(rb(th))==FALSE) { #check the function at the start point is real 
  #When f is not real the function is not analytic. 
  stop("The objective function is not finite at the starting point. \n 
          The function is not analytic.")
}


# Second Error testing 


if ( any( is.finite( gb(c(0,0)))) == FALSE ) { #check the gradient at the start is real 
  stop("The gradient is not finite at the starting point. \n 
          The function is not analytic.")
}else{cat ("hi")}


if ( any( is.finite( gb(c(Inf,0)))) == FALSE ) { #check the gradient at the start is real 
  stop("The gradient is not finite at the starting point. \n 
          The function is not analytic.")
}else{cat ("hi")}

#this config also works for the 2nd warning

if ( any(is.finite(c(0, 0))==FALSE) ) { #check the gradient at the start is real 
  stop("The gradient is not finite at the starting point. \n 
          The function is not analytic.")
}else{cat("ello")}

if ( any(is.finite(gb(c(Inf, 0)))==FALSE) ) { #check the gradient at the start is real 
  stop("The gradient is not finite at the starting point. \n 
          The function is not analytic.")
}


##### Testing Hessian Finite Difference Function 


th <- c(0,0)
th1 <-c(1e-6,0)


rb(th)
gb(th)
hb(th)

gb(th)
gb(th1)

nullhess(th, gb)

### from the results we can see the analytic and numerical outputs are similar 
### error in element 1,2 comes from numerical error 


### test max half step warning:

if (1 >= 0) { 
  stop( paste( "The next step cannot be found. \n  
               The step size has been halved ", 5, " times.")
        , call. = FALSE)
}




#### Testing tolerance condition 

th <- c(0,0)


f<-rb(th)
gtest <- gb(th)
hb(th)

abs(gtest)

any(abs(gtest) <= measure)

abs(gtest) <= measure

all(abs(gtest) <= measure)

sum(abs(gtest) <= measure)
prod(abs(gtest) <= measure)
tol=1e-8
fscale = 1

measure <- tol * (abs(f) + fscale ) 





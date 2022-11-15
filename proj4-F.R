## The 'newt' function has inputs, including theta, func, grad, hess


newt <- function(theta, #this is our initial search point 
                 func, #this is the functional form of the function
                 grad, #this is the functional form of the gradient 
                 hess=NULL, #this is the functional form of the hessian, 
                            #when no hessian function is given this variables becomes NULL 
                 ..., #allows newt to absorb the additional parameters required by 
                      #func, grad and hess
                 tol=1e-8, #tolerance is the smallest difference between our steps permitted
                 fscale=1, #fscale is the order of the function's values around the minimum
                 maxit=100, #maxit is the maximum number of iterations of newton's method
                            #carried out
                 max.half=20, # the maximum number of times we half our step size
                              #while searching for a next step that returns a lower 
                              #function value
                 eps=1e-6) {  #epsilon is the step in x_i used to calculate the derivative 
                              #w.r.t. x_i using finite differences

  f <- func(theta,...) #Store the function, evaluated at the first setp theta
  
  gradval <- grad(theta,...) # Store the gradient vector evaluated at the first step theta 
  
  #Check that the function and gradient are continuous (smooth)
  #This is true when their values are finite 
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
    
    len <- length(gradval) # length of the gradient vector 
    
    Hfd <- matrix(0, len,len) ## finite difference Hessian
    
    for (i in 1:length(theta)) {# loop over each variable basis 
      
      th1 <- theta #store each dimension of theta to operate on 
      
      th1[i] <- th1[i] + eps  # increase th1[i] by eps
      
      grad1 <- grad(th1,...) # compute new derivative value where a parameter i
                             # has been shifted by epsilon 
      
      #approximating the derivative of the function with the finite difference method
      # (g(x_i + eps) - g(x_i))  /  ((x_i + eps ) - x_i)
      # All other variables x_j (j<>i) are kept the same 
      Hfd[i,] <- (grad1 - gradval)/eps #ith row of the hessian 
      #row 1 = derivatives of g(X) w.r.t x_1 
    }
    h <- Hfd #store the hessian
  } else { #if the hessian matrix is supplied then store it locally 
    h <- hess
  }
  
  #find the inverse hessian (no matter the definiteness)
  
  # Carry out QR decomposition of h (the hessian)
  
  R <- qr.R(qr(h))  # get the triangle matrix R as the output 
  Q <- qr.Q(qr(h))  # get the orthogonal matrix Q as the output 
  
  #Determine R inverse: 
  #make an identity matrix with the same shape as R 
  Identity_for_R <- diag(nrow = dim(R)[1])
  # use backwards substitution to obtain R inverse 
  # by solving for: R R^-1 = I
  R_inv <- backsolve(R, Identity_for_R)
  #Q inverse = Q transpose 
  
  #QR = h => hi = R_inv Q_transpose = inverse hessian 
  hi <- R_inv %*% t(Q)
  
  
  #Implement Newton's Method: 
  # theta_(k+1) = theta_k - f''(theta_k)^(-1) * f'(theta_k)
  # f''(theta_k)^(-1) = hi @ k th step 
  # f'(theta_k) = gradval @ k th step 
  # f(theta_k) = f @ k th step 
  # theta_k = k th param val scanned
  # theta_(k+1) = next theta point to scan (should be closer to the min)
  
  
  n_iter <- 0 # count the number of iterations 
  
  while (n_iter < maxit ){ # iteration until maximum iteration reached 
    theta2 <-  theta - hi %*% gradval #find the next step 
    f2 <- func(theta2,...)                #evaluate the fn at next step 
    
    #if the fn increases rather than decreases, half the step size in the same direction
    n_half <- 0  # counter for times steps halved 
    while ( f2 > f ){ #iterate until fn val decreases 
      n_half <- n_half + 1 #counter goes up 
      
      theta2 <-  theta - ( hi %*% gradval ) * ( 0.5 ^(n_half) )
      f2 <- func(theta2,...) 
      
      if (n_half >= max.half ) { 
        stop()}
      
      
    }
    
    
    
    
    
    
    n_iter <- n_iter + 1 #iteration num count increases
  }
  
  
  
  
  
  
  
  
  
  
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



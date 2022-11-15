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
  
  gradval <- grad(theta,...) # Store the graident vector evaluated at the first step theta 
  
  #Check that the function and gradient are continuous (smooth)
  #This is true when their values are finite 
  if (is.finite(f)==FALSE) { #check the function at the start point is real 
    #When f is not real the function is not analytic. 
    stop("The objective function is not finite at the starting point. \n 
          The function is not analytic.")
  }
  if (is.finite(gradval)==FALSE) { #check the gradient at the start is real 
    stop("The gradient is not finite at the starting point. \n 
          The function is not analytic.")
  }
  
  #If the Hessian is not supplied, calculate it through finite difference methods 
  
  if (hess=NULL) { #no hessian supplied 
    
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
  #f(x_(k+1)) = f(x_k) - f''(x_k)^(-1) * f'(x_k)
  #f''(x_k)^(-1) = hi @k  
  #f'(x_k) = gradval @ k 
  #f(x_k) = f @ k 
  
  
  
  
  
  
  
  
  
  
  
  
}





#Paste at the function end: 

#check hessian is positive semi-definite at the minimum.
#Code: 
options(show.error.messages = TRUE, call. = FALSE)
R <- try(chol(h) , stop("Hessian is not positive semi-definite at the minimum"))

#Hessian Inverse hi - using cholesky (only works if chol works)
hi <- chol2inv(chol(A))














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

a<-myfunc(3)









# x-squared function, or x^T x for n dimensional x. Works for any n.
# should return 0 vector with a min of 0
# R^n -> R
xtx <- function(x){
  t(x) %*% x
}

xtxdiv <- function(x){
  2*x
}

xtxhess <- function(x){
  diag(2,length(x))
}

# difference of cubes. R^2 -> R
diff3 <- function(x){
  x[1]**3 - x[2]**3
}
diff3derv <- function(x){
  c(3*x[1]**2,3*x[2]**2)
}
diff3hess <- function(x){
  h <- matrix(0,2,2)
  h[1,1] <- 8*x[1]
  h[2,2] <- -6*x[2]
  h
}


# Simon's example, should return (1,1) with a min of 0
# R^2 -> R
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


# wenqi function 1. R^2 -> R
W1f <- function(x,k=10){
  2*(x[2]-x[1]^2)^2 + (1-x[1])^2+k
}

W1g <- function(x,k=10){
  c(-8*x[1]*x[2] + 8*x[1]^3 - 2 + 2*x[1], 4*(x[2] - x[1]^2))
}

W1h <- function(x,k=10){
  h <- matrix(0,2,2)
  h[1,1] <- -8*x[2] + 24*x[1]^2 + 2
  h[2,2] <- 4
  h[1,2] <- h[2,1] <- -8*x[1]
  h
}

# wenqi function 2. R^3 -> R
W2f <- function(x){
  x[1]*exp(x[2] + x[3]^2)
}

W2g <- function(x){
  c(exp(x[2] + x[3]^2), x[1]*exp(x[2] + x[3]^2), 2*x[1]*x[3]*exp(x[2] + x[3]^2))
}

W2h <- function(x){
  h <- matrix(0,3,3)
  h[1,1] <- 0
  h[2,2] <- x[1]*exp(x[2] + x[3]^2)
  h[3,3] <- 2*x[1]*exp(x[2] + x[3]^2)
  h[1,2] <- h[2,1] <- exp(x[2] + x[3]^2)
  h[1,3] <- h[3,1] <- 2*x[3]*exp(x[2] + x[3]^2)
  h[2,3] <- h[3,2] <- 2*x[1]*x[3]*exp(x[2] + x[3]^2)
  h
}

# step function R -> R
step <- function(x){
  if (x < 0){
    x^2
  } else {
    1
  }
}
stepder <- function(x){
  if(x < 0){
    2*x
  } 
  if (x > 0) {
    0
  } else {
    Inf
  }
}
stephess <- function(x){
  if(x<0){
    2
  } 
  if (x > 0) {
    0
  } else {
    NULL
  }
}

# many minima function. R^2 -> R
xsinx <- function(x){
  x[1]*sin(x[1]) + x[2]*sin(x[2])
}

nablaxsinx <- function(x){
  c(sin(x[1]) + x[1]*cos(x[1]),sin(x[2]) + x[2]*cos(x[2]))
}

hessxsinx <- function(x){
  h <- matrix(0,2,2)
  h[1,1] <- 2*cos(x[1]) - x[1]*sin(x[1])
  h[2,2] <- 2*cos(x[2]) - x[2]*sin(x[2])
  h[1,2] <- h[2,1] <- 0
  h
}
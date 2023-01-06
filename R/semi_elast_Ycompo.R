semi_elast_Ycompo <- function(E_Y, b, simplex = T, D = ncol(Y_simp)) {
  
  # initialization
  n <- nrow(E_Y)
  res <- matrix(0, n, D)
  if (simplex)
    b <- matrix(clr(b), nrow = 1, ncol = D)
  
  for (k in 1:n) {
    z <- E_Y[k, ]
    
    if (simplex) {
      w <- diag(D) - rep(1, D) %*% t(z)
    } else {
      w <- (diag(D) - rep(1, D) %*% t(z)) %*% ilrBase(D = D)
    }
    
    res[k, ] <- w %*% t(b)
    
  }
  return(res)
}
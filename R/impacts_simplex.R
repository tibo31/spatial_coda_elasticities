map_elas_ternary <- function(df_coords, 
                              X, b_star, index_b, W, RHO, delta = 10^(-8),
                              range, size = 1, length = 0.1, cex.legend = 0.5,
                              legend.ternary = c("", "", "")) {
  
  # initialisation 
  n <- nrow(df_coords)
  D <- ncol(b_star) + 1
  res <- array(0, dim = c(n, n, D))
  W_Mat <- as(W, "Matrix")
  A_w <- Diagonal(n * (D - 1)) - kronecker(RHO, W_Mat)
  A_w_inv <- solve(A_w)
  E_Y_ilr <- A_w_inv %*% as.vector(X %*% b_star)
  E_Y <- as(ilrInv(matrix(E_Y_ilr, ncol = D - 1)), "matrix")
  
  # compute the semi-elasticities
  elast_sp <- sp_semi_elast_Ycompo(X, b_star, index_b = index_b, type = "appro",
                                   delta = delta, W = W, 
                                   RHO = matrix(c(0.65, 0, 0.18, 0.63), 2, 2))
  
  # the prediction 
  pred_y_simp <- E_Y * (1 + delta * cbind(apply(elast_sp[ , , 1], 2, sum),
                                          apply(elast_sp[ , , 2], 2, sum),
                                          apply(elast_sp[ , , 3], 2, sum)))
  
  # the dimension of the triangle: s1 is the basis, s2 is the height
  s1 <- size * range / 10
  s2 <- s1 * sqrt(3) / 2
  
  for (k in 1:n) {
    # define the individual triangle 
    x <- df_coords[k, 1]
    y <- df_coords[k, 2]
    A <- c(x - s1 / 2, y - s2 / 3)
    B <- c(x + s1 / 2, y - s2 / 3)
    C <- c(x, y + 2 * s2 / 3)
    
    lines(c(A[1], B[1]), c(A[2], B[2]), col = rgb(0.5, 0.5, 0.5)) # "#E16A86") # r
    lines(c(A[1], C[1]), c(A[2], C[2]), col = rgb(0.5, 0.5, 0.5)) # "#50A315") # g
    lines(c(B[1], C[1]), c(B[2], C[2]), col = rgb(0.5, 0.5, 0.5)) # "#009ADE") # b
    
    # legend 
    text(c(A[1], B[1], C[1]), c(A[2], B[2], C[2]), 
         legend.ternary, pos = c(3, 3, 2), 
         cex = cex.legend)
    
    # the adjusted values in the simplex
    E_Y_simplex_x <- E_Y[k, 1] * A[1] + E_Y[k, 2] * B[1] + E_Y[k, 3] * C[1] 
    E_Y_simplex_y <- E_Y[k, 1] * A[2] + E_Y[k, 2] * B[2] + E_Y[k, 3] * C[2] 
    points(E_Y_simplex_x, E_Y_simplex_y, pch = 16, cex = 0.5)
    
    # the predicted values in the triangle
    Y_pred_x <- pred_y_simp[k, 1] * A[1] + pred_y_simp[k, 2] * B[1] + pred_y_simp[k, 3] * C[1] 
    Y_pred_y <- pred_y_simp[k, 1] * A[2] + pred_y_simp[k, 2] * B[2] + pred_y_simp[k, 3] * C[2] 
    points(Y_pred_x, Y_pred_y, pch = 16, cex = 0.5, col = "#E16A86")
    
    q4 <- colorspace::qualitative_hcl(3, palette = "Dark 3")
    
    arrows(E_Y_simplex_x, E_Y_simplex_y, Y_pred_x, Y_pred_y, 
           col = "#E16A86", length = length)
    
    Y_simplex_x_start <- E_Y_simplex_x
    Y_simplex_y_start <- E_Y_simplex_y
    
    for(l in 1:n) {
      
      pred_y_simp_ind <- c(
        elast_sp[l, k, 1] * (delta * E_Y[k, 1]) + E_Y[k, 1],
        elast_sp[l, k, 2] * (delta * E_Y[k, 2]) + E_Y[k, 2],
        elast_sp[l, k, 3] * (delta * E_Y[k, 3]) + E_Y[k, 3]
      )
      
      Y_simplex_x_end <- pred_y_simp_ind[1] * A[1] + pred_y_simp_ind[2] * B[1] +
        pred_y_simp_ind[3] * C[1]
      Y_simplex_y_end <- pred_y_simp_ind[1] * A[2] + pred_y_simp_ind[2] * B[2] +
        pred_y_simp_ind[3] * C[2]
      
      if(l != 1) {
        my_vec_x <- (Y_simplex_x_start - E_Y_simplex_x)
        my_vec_y <- (Y_simplex_y_start - E_Y_simplex_y)
        Y_simplex_x_end <- Y_simplex_x_end + my_vec_x
        Y_simplex_y_end <- Y_simplex_y_end + my_vec_y
      }
      
      arrows(Y_simplex_x_start, Y_simplex_y_start,
             Y_simplex_x_end, Y_simplex_y_end,
             col = ifelse(l == k, "#50A315", "#009ADE"),
             length = length)
      
      Y_simplex_x_start <- Y_simplex_x_end
      Y_simplex_y_start <- Y_simplex_y_end
      
    }
  }
  
}

my_rgb_coda <- function(my_data, max_intensity = 1.5) {
  
  # check that data are CODA
  # stopifnot(all(apply(my_data, 1, sum) == sum(my_data[1, ])))
  
  my_don_col <- character(nrow(my_data))
  # give color to each composition
  for(i in 1:nrow(my_data)) {
    if(max(my_data[i, ]) <= 1 / max_intensity) {
      my_don_col[i] <- rgb(my_data[i, 1] * max_intensity, 
                           my_data[i, 2] * max_intensity,
                           my_data[i, 3] * max_intensity)
    } else {
      my_rate <- 1 / max(my_data[i, ])
      my_rgb_local <- my_rate * my_data[i, ]
      my_don_col[i] <- rgb(my_rgb_local[1], 
                           my_rgb_local[2],
                           my_rgb_local[3])
    }
  }
  return(my_don_col)
}
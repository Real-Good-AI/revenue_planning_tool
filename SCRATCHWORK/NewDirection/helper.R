get_likelihood_sigma <- function(row, dist_mat, rho, data){
      # Compute Matérn Covariance Matrix
      mat_cov <- Matern(d = dist_mat, smoothness = row[1], range = rho, phi = 1) + diag(row[2], dim(dist_mat)[1])
      
      # Compute sigma^2
      sigma.squared <- (1/length(data)) * (data %*% chol2inv(chol(mat_cov)) %*% data)[1,1] 
      
      # Compute likelihood
      likelihood <- dmvnorm(data, 
                            mean = rep(0, length(data)), 
                            sigma = sigma.squared * mat_cov,
                            log = TRUE)
      return(likelihood)
}
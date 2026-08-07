library(readr)
library(dplyr)
library(tidyverse)
library(fields)
library(mvtnorm)
library(doParallel)

candidates_all <- readRDS("ORG_random_candidates50.rds")
df <- readRDS("df_eligible_orgs.rds") |> select(-DATA_COUNT, -keep)
all_orgs <- unique(df$EIN2)

# distance matrix of years
dist_mat_all <- as.matrix(dist(sort(unique(df$TAX_YEAR) - min(df$TAX_YEAR) + 1), diag=TRUE, upper=TRUE))

best_all <- readRDS("ORG_best_empirical_all.rds")

rho = 1

cluster <- makeCluster(6)
registerDoParallel(cluster)

# note: .inorder = FALSE so the outer list won't be in trial order... I think for our purposes ok
res <- foreach(trialID = 1:length(candidates_all), .packages = c("dplyr", "fields", "mvtnorm"), .inorder = FALSE) %do% { 
      reject_null <- list()
      # print(paste0("Trial: ", trialID))
      for (orgID in all_orgs){
            ID <- paste0(orgID, ".", trialID)
            # fileConn <- file(paste0("errors/Trial", trialID, "_errors.txt"))
            if (!(ID %in% names(best_all))){
                  # fileConn <- file(paste0("errors/Trial", trialID, "_errors.txt"))
                  # write(x = ID, file = fileConn, append = TRUE, sep = "\n")
                  next
            }
            # close(fileConn)
            # pull their time series
            my_data <- df |> filter(EIN2 == orgID) |>
                  mutate(TOT_REV = log1p(TOT_REV),
                         TOT_EXP = log1p(TOT_EXP),
                         TOT_ASSET = log1p(TOT_ASSET))
            all_years <- my_data$TAX_YEAR - min(df$TAX_YEAR) + 1 # all the years we have data for this county
            log_rev <- my_data |> merge(data.frame(TAX_YEAR = 1991:2019), all = TRUE) |> arrange(TAX_YEAR) |> pull(TOT_REV)
            
            candidate <- candidates_all[[trialID]][[orgID]] - min(df$TAX_YEAR) + 1
            
            # pull best parameters
            best.curr <- best_all[[ID]] # best.curr$left$nu
            
            # Determine the left and right time series corresponding to the current candidate change point
            left_years <- all_years[1:which(all_years == candidate)-1]
            right_years <- all_years[which(all_years == candidate):length(all_years)]
            
            ### Null Hypothesis Model ###
            curr_all_years <- c(left_years, right_years) # Redundant for single candidate
            centered_rev <- log_rev[curr_all_years] - mean(log_rev[curr_all_years]) # center the log revenue values so they are mean 0
            
            dist_mat <- dist_mat_all[curr_all_years, curr_all_years]
            dist_mat <- (dist_mat - min(dist_mat)) / (max(dist_mat) - min(dist_mat))
            Sigma_0 <- Matern(d = dist_mat,
                              smoothness = best.curr$null$nu,
                              range = rho,
                              phi = 1) + diag(best.curr$null$nugget, dim(dist_mat)[1])
            
            sigma.sqrd_0 <- (1/length(centered_rev)) * (centered_rev %*% chol2inv(chol(Sigma_0)) %*% centered_rev)[1,1]
            Sigma_0 <- sigma.sqrd_0 * Sigma_0
            
            ### Alternative Hypothesis Model ###
            # center the log revenue values so they are mean 0
            centered_rev.L <- log_rev[left_years] - mean(log_rev[left_years])
            centered_rev.R <- log_rev[right_years] - mean(log_rev[right_years])
            
            # compute and normalize distance matrices
            dist_mat.L <- dist_mat_all[left_years, left_years]
            dist_mat.L <- (dist_mat.L - min(dist_mat.L)) / (max(dist_mat.L) - min(dist_mat.L))
            
            dist_mat.R <- dist_mat_all[right_years, right_years]
            dist_mat.R <- (dist_mat.R - min(dist_mat.R)) / (max(dist_mat.R) - min(dist_mat.R))
            
            # Compute Matern covariance matrices
            Sigma_L <- Matern(d = dist_mat.L,
                              smoothness = best.curr$left$nu,
                              range = rho,
                              phi = 1) + diag(best.curr$left$nugget, dim(dist_mat.L)[1])
            sigma.sqrd_L <- (1/length(left_years)) * (centered_rev.L %*% chol2inv(chol(Sigma_L)) %*% centered_rev.L)[1,1]
            Sigma_L <- sigma.sqrd_L * Sigma_L
            
            Sigma_R <- Matern(d = dist_mat.R,
                              smoothness = best.curr$right$nu,
                              range = rho,
                              phi = 1) + diag(best.curr$right$nugget, dim(dist_mat.R)[1])
            sigma.sqrd_R <- (1/length(right_years)) * (centered_rev.R %*% chol2inv(chol(Sigma_R)) %*% centered_rev.R)[1,1]
            Sigma_R <- sigma.sqrd_R * Sigma_R
            
            ### Likelihood Ratio Test ###
            L_0 <- dmvnorm(centered_rev, 
                           mean = rep(0, length(centered_rev)), 
                           sigma = Sigma_0,
                           log = TRUE)
            
            L_left <- dmvnorm(centered_rev.L, 
                              mean = rep(0, length(centered_rev.L)), 
                              sigma = Sigma_L,
                              log = TRUE)
            
            L_right <- dmvnorm(centered_rev.R, 
                               mean = rep(0, length(centered_rev.R)), 
                               sigma = Sigma_R,
                               log = TRUE)
            
            L_alt <- L_left + L_right
            
            
            # likelihood ratio test statistic
            lrt <- -2*L_0 + 2*L_alt
            doF <- 4 #  (mu1, mu2, Sigma1, Sigma2, nu1, nu2, nugget1, nugget2) vs (mu, Sigma, nu, nugget)
            reject_null[[orgID]] <- (lrt > qchisq(p = 0.05, df = doF, lower.tail = FALSE))
            
      }
      # reject_null_all[[as.character(trialID)]] <- reject_null
      reject_null
}

res.final <- unlist(lapply(res, function(x){sum(unlist(x)) / length(x)}))
hist(res.final)
sd(res.final)

stopCluster(cl = cluster)

saveRDS(res, "org_empirical_lrt_res.rds")

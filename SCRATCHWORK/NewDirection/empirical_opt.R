library(data.table)
library(readr)
library(dplyr)
library(tidyverse)
library(fields)
library(mvtnorm)
library(doParallel)
source("helper.R")

candidates_all <- readRDS("random_candidates50.rds")
df <- readRDS("county_long.rds") |> 
      filter(TAX_YEAR < 2020) |> 
      select(all_of(c("TAX_YEAR", "geoid_2010", "TOT_REV", "TOT_EXP", "TOT_ASSET", "county_full", "state_name", "REGION", "DIVISION", "TotDmgADJ", "n_disasters", "TotPerCapADJ", "FEMA_threshADJ"))) |>
      filter(geoid_2010 %in% names(candidates_all$`1`)) |>
      arrange(geoid_2010, TAX_YEAR)

# distance matrix of years
dist_mat_all <- as.matrix(dist(sort(unique(df$TAX_YEAR) - min(df$TAX_YEAR) + 1), diag=TRUE, upper=TRUE))

init_conds <- readRDS("empirical_init_conds_50.rds")$centers

best_parameters <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("type", "geoid", "nu", "nugget", "trial"))

errs <- c()
combos <- expand.grid(unique(df$geoid_2010), as.character(1:length(candidates_all))) |> arrange(Var1, Var2)
midPt <- nrow(combos)/2

cluster <- makeCluster(6)
registerDoParallel(cluster)

start.time <- Sys.time() 
resBest <- foreach(i = (midPt+1):nrow(combos), .combine = 'rbind', .packages = c("dplyr", "fields", "mvtnorm"), .verbose = TRUE, .errorhandling = "remove") %dopar% {
      countyID <- combos$Var1[i]
      trialID <- combos$Var2[i]
      print(paste("County ID:", countyID, "Trial Number:", trialID))
      # pull their time series
      my_data <- df |> filter(geoid_2010 == countyID) |>
            mutate(TOT_REV = log1p(TOT_REV),
                   TOT_EXP = log1p(TOT_EXP),
                   TOT_ASSET = log1p(TOT_ASSET))
      all_years <- my_data$TAX_YEAR - min(df$TAX_YEAR) + 1 # all the years we have data for this county
      log_rev <- my_data |> merge(data.frame(TAX_YEAR = 1991:2019), all = TRUE) |> arrange(TAX_YEAR) |> pull(TOT_REV)
      
      candidate <- candidates_all[[trialID]][[countyID]] - min(df$TAX_YEAR) + 1
      left_years <- all_years[1:which(all_years == candidate)-1]
      right_years <- all_years[which(all_years == candidate):length(all_years)]
      
      ### Null Model ###
      curr_all_years <- c(left_years, right_years)
      centered_rev <- log_rev[curr_all_years] - mean(log_rev[curr_all_years])
      dist_mat <- dist_mat_all[curr_all_years, curr_all_years]
      dist_mat <- (dist_mat - min(dist_mat)) / (max(dist_mat) - min(dist_mat))
      
      likelihood_wrapper.null <- function(pars.vec){
            return(-1 * get_likelihood_sigma(row = pars.vec,
                                             dist_mat = dist_mat,
                                             rho = 1,
                                             data = centered_rev))
      }
      
      tst <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("nu.init", "nugg.init", "nu", "nugget", "likelihood", "trial"))
      i <- 2      
      for (i in 1:nrow(init_conds)){
            res <- constrOptim(theta = init_conds[i,], 
                               f = likelihood_wrapper.null,
                               grad = NULL,
                               ui = rbind(c(1,0),c(-1,0),c(0,1)), #rbind(c(1,0),c(-1,0),c(0,1))
                               ci = c(0.001,-3.5,0)) #c(0.001,-3.5, 0)
            tst[nrow(tst)+1,] <- c(c(init_conds[i,1], init_conds[i,2], res$par[1], res$par[2], -1*res$value, i))
      }
      res.null <- as.data.frame(list(type = "null",
                                     geoid = countyID,
                                     nu = tst[which.max(tst$likelihood),]$nu,
                                     nugget = tst[which.max(tst$likelihood),]$nugget,
                                     trial = trialID))
      
      ### Alternative Model ###
      
      # Left
      centered_rev.L <- log_rev[left_years] - mean(log_rev[left_years])
      if (var(centered_rev.L) == 0){
            cat(i, file = file.path("errors/", paste0(i, ".txt")))
      }
      dist_mat.L <- dist_mat_all[left_years, left_years]
      dist_mat.L <- (dist_mat.L - min(dist_mat.L)) / (max(dist_mat.L) - min(dist_mat.L))
      
      likelihood_wrapper.left <- function(pars.vec){
            return(-1 * get_likelihood_sigma(row = pars.vec,
                                             dist_mat = dist_mat.L,
                                             rho = 1,
                                             data = centered_rev.L))
      }
      
      tst <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("nu.init", "nugg.init", "nu", "nugget", "likelihood", "trial"))
      for (i in 1:nrow(init_conds)){
            res <- constrOptim(theta = init_conds[i,], 
                               f = likelihood_wrapper.left,
                               grad = NULL,
                               ui = rbind(c(1,0),c(-1,0),c(0,1)), #rbind(c(1,0),c(-1,0),c(0,1))
                               ci = c(0.001,-3.5,0)) #c(0.001,-3.5, 0)
            tst[nrow(tst)+1,] <- c(c(init_conds[i,1], init_conds[i,2], res$par[1], res$par[2], -1*res$value, i))
      }
      res.left <- as.data.frame(list(type = "left",
                                     geoid = countyID,
                                     nu = tst[which.max(tst$likelihood),]$nu,
                                     nugget = tst[which.max(tst$likelihood),]$nugget,
                                     trial = trialID))
      
      # Right
      centered_rev.R <- log_rev[right_years] - mean(log_rev[right_years])
      if (var(centered_rev.R) == 0){
            cat(i, file = file.path("errors/", paste0(i, ".txt")))
      }
      dist_mat.R <- dist_mat_all[right_years, right_years]
      dist_mat.R <- (dist_mat.R - min(dist_mat.R)) / (max(dist_mat.R) - min(dist_mat.R))
      
      likelihood_wrapper.right <- function(pars.vec){
            return(-1 * get_likelihood_sigma(row = pars.vec,
                                             dist_mat = dist_mat.R,
                                             rho = 1,
                                             data = centered_rev.R))
      }
      
      tst <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("nu.init", "nugg.init", "nu", "nugget", "likelihood", "trial"))
      for (i in 1:nrow(init_conds)){
            res <- constrOptim(theta = init_conds[i,], 
                               f = likelihood_wrapper.right,
                               grad = NULL,
                               ui = rbind(c(1,0),c(-1,0),c(0,1)), #rbind(c(1,0),c(-1,0),c(0,1))
                               ci = c(0.001,-3.5,0)) #c(0.001,-3.5, 0)
            tst[nrow(tst)+1,] <- c(c(init_conds[i,1], init_conds[i,2], res$par[1], res$par[2], -1*res$value, i))
      }
      rbind(res.null,
            res.left,
            as.data.frame(list(type = "right",
                               geoid = countyID,
                               nu = tst[which.max(tst$likelihood),]$nu,
                               nugget = tst[which.max(tst$likelihood),]$nugget,
                               trial = trialID)))
      
}
timeDiff <- Sys.time() - start.time
print(timeDiff)

stopCluster(cl = cluster)

saveRDS(resBest, "best_empirical_2ndHalf.rds")


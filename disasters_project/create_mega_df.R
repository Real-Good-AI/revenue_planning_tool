library(data.table)
library(readr)

pf <- readRDS("data/pf_processed.rds")
pz <- readRDS("data/pz_processed.rds")
pz.ce <- readRDS("data/pz_ce_processed.rds")

# Keep only assets, rev, and expenses variables for private foundations to match PZ
pf <- pf[, !c("PF_01_REV_GRO_PROFIT_BOOKS", "PF_02_LIAB_TOT_EOY_BV", "PF_01_EXP_CONTR_PAID_BOOKS", "PF_01_EXP_TOT_EXP_DISBMT_DISBMT"), with = FALSE]

# rename relevant columns
setnames(pf,
         old = c("PF_02_ASSET_TOT_EOY_BV", "PF_01_REV_TOT_BOOKS", "PF_01_EXP_TOT_EXP_DISBMT_BOOKS", "F990_TOTAL_ASSETS_RECENT"),
         new = c("TOT_ASSET", "TOT_REV", "TOT_EXP", "SIZE"))

pf[, SOURCE := "PF"]
pz[, SOURCE := "PZ"]
pz.ce[, SOURCE := "PZ.CE"]

dt <- rbindlist(list(pf, pz, pz.ce), use.names = TRUE)
dt <- unique(dt)

rm(pf, pz, pz.ce)

# exploring duplicates in the mega dt
source("../SCRIPTS/clean_helper.R")

key_var <- c("EIN2", "TAX_YEAR")
dup_groups <- dt[, .N, by = key_var][N > 1] #
dups <- dt[dup_groups, on = key_var]

dollar_cols <- c("TOT_ASSET", "TOT_REV", "TOT_EXP")
diagnostics <- dups[, compare_pair_dt(.SD, dollar_cols), by = .(EIN2, TAX_YEAR)] #information about the duplicates
probs <- seq(0,1,0.1)
quantile(diagnostics$max_abs_diff, probs = probs)

# There are 370 (EIN2, TAX_YEAR) pairs where the records are exact duplicates EXCEPT in the source column... how to deal with?

dt[, DATA_COUNT := .N, by = c("EIN2")] #recalculate because now may have from different sources... but will be overcount because of duplicates
quantile(dt$DATA_COUNT, probs = probs)

# saveRDS(dt, "data/mega.rds") # uncomment this line of code to save the file

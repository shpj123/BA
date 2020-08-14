library(MASS)
library(lavaan)
library(tidyverse)
library(doRNG)
library(doParallel)
#### data generation 
ncores <- -1 + detectCores()

create_model <- function(ntrials) {
  # lavaan Syntax for LDS
  result <- ""
  prefix <- "Y"
  result <- paste0(result,
                   "eta_comp =~ ",
                   paste0("1*", prefix, 1:ntrials, "_comp", collapse = " + "))
  result <- paste0(
    result,
    "\n",
    "eta_incomp =~ ",
    paste0("1*", prefix, 1:ntrials, "_incomp", collapse = " + ")
  )
  result <- paste0(result,
                   "\n",
                   paste0(prefix, 1:ntrials, "_comp ~ 0*1", collapse = "\n"))
  result <- paste0(result,
                   "\n",
                   paste0(prefix, 1:ntrials, "_incomp ~ 0*1", collapse = "\n"))
  result <- paste0(
    result,
    "\n",
    paste0(
      prefix,
      1:ntrials,
      "_comp ~~ eps*",
      prefix,
      1:ntrials,
      "_comp",
      collapse = "\n"
    )
  )
  result <- paste0(
    result,
    "\n",
    paste0(
      prefix,
      1:ntrials,
      "_incomp ~~ eps*",
      prefix,
      1:ntrials,
      "_incomp",
      collapse = "\n"
    )
  )
  result <- paste0(
    result,
    "\n",
    "eta_comp ~ 1\n",
    "diff =~ eta_incomp\n",
    "diff ~ mu_diff*1\n",
    "eta_incomp ~ 0*1 + 1*eta_comp\n",
    "eta_comp ~~ eta_comp + diff + 0*eta_incomp\n",
    "diff ~~ var_diff*diff + 0*eta_incomp\n",
    "eta_incomp ~~ 0*eta_incomp"
  )
  result
}


run_single_sim <-
  function(nsim, sample_size, effect_size, ntrials) {
    samp <-
      get_data(sample_size = sample_size,
               effect_size = effect_size,
               ntrials = ntrials)
    
    singleresult <- tryCatch({
      model <- create_model(ntrials = ntrials)
      fit <- sem(model, samp)
      
      mu_diff <-
        parTable(fit)$est[parTable(fit)$label == "mu_diff"]
      se_mu_diff <-
        parTable(fit)$se[parTable(fit)$label == "mu_diff"]
      var_diff <-
        parTable(fit)$est[parTable(fit)$label == "var_diff"]
      p_value <- (1 - pnorm(abs(mu_diff / se_mu_diff))) * 2
      cohens_d <- mu_diff / sqrt(var_diff)
      
      data.frame(
        mu_diff = mu_diff,
        se_mu_diff = se_mu_diff,
        var_diff = var_diff,
        p_value = p_value,
        cohens_d = cohens_d,
        sample_size = sample_size,
        effect_size = effect_size,
        ntrials = ntrials,
        method = "LDS",
        msg = NA,
        nsim = nsim
      )
    }, warning = function(w) {
      model <- create_model(ntrials = ntrials)
      fit <- sem(model, samp)
      
      mu_diff <-
        parTable(fit)$est[parTable(fit)$label == "mu_diff"]
      se_mu_diff <-
        parTable(fit)$se[parTable(fit)$label == "mu_diff"]
      var_diff <-
        parTable(fit)$est[parTable(fit)$label == "var_diff"]
      p_value <- (1 - pnorm(abs(mu_diff / se_mu_diff))) * 2
      cohens_d <- mu_diff / sqrt(var_diff)
      
      data.frame(
        mu_diff = mu_diff,
        se_mu_diff = se_mu_diff,
        var_diff = var_diff,
        p_value = p_value,
        cohens_d = cohens_d,
        sample_size = sample_size,
        effect_size = effect_size,
        ntrials = ntrials,
        method = "LDS",
        msg = w$message,
        nsim = nsim
      )
    }, error = function(e) {
      data.frame(
        mu_diff = NA,
        se_mu_diff = NA,
        var_diff = NA,
        p_value = NA,
        cohens_d = NA,
        sample_size = sample_size,
        effect_size = effect_size,
        ntrials = ntrials,
        method = "LDS",
        msg = e$message,
        nsim = nsim
      )
    })
    
    ttest <-
      t.test(samp$sum_compatible, samp$sum_incompatible, paired = T)
    
    mu_diff <- mean(samp$sum_incompatible - samp$sum_compatible)
    var_diff <- var(samp$sum_incompatible - samp$sum_compatible)
    p_value <- ttest$p.value
    cohens_d <- mu_diff / sqrt(var_diff)
    
    singleresult <- bind_rows(
      singleresult,
      data.frame(
        mu_diff = mu_diff,
        se_mu_diff = sqrt(var_diff / sample_size),
        var_diff = var_diff,
        p_value = p_value,
        cohens_d = cohens_d,
        sample_size = sample_size,
        effect_size = effect_size,
        ntrials = ntrials,
        method = "t-Test",
        nsim = nsim
      )
    )
    
    singleresult
  }


run_sim <- function(...) {
  dots <- list(...)
  results <-
    lapply(1:2000, function(nsim)
      do.call(run_single_sim, c(nsim, dots)))
  
  do.call(rbind, results)
}

#### ntrials = 2 | 5
get_data <- function(sample_size, effect_size, ntrials) {
  mu <-
    c(1, 1 + effect_size) # due to the specific variance and corvariance used in this simulation
  vcov <- matrix(0.5, nrow = 2, ncol = 2)
  diag(vcov) <- 1
  resid <- sqrt(2)
  
  samp <- mvrnorm(n = sample_size, mu = mu, Sigma = vcov)
  samp <- as.data.frame(samp)
  names(samp) <- c("true_compatible", "true_incompatible")
  for (ntrial in 1:ntrials) {
    samp <- cbind(samp,
                  samp$true_compatible + rnorm(sample_size, sd = resid))
    names(samp)[ncol(samp)] <- paste0("Y", ntrial, "_comp")
  }
  
  for (ntrial in 1:ntrials) {
    samp <- cbind(samp,
                  samp$true_incompatible + rnorm(sample_size, sd = resid))
    names(samp)[ncol(samp)] <- paste0("Y", ntrial, "_incomp")
  }
  
  samp$sum_compatible <- rowMeans(samp[, 3:(3 + ntrials - 1)])
  samp$sum_incompatible <-
    rowMeans(samp[, (3 + ntrials):(3 + 2 * ntrials - 1)])
  samp
}

conditions <-
  expand.grid(
    sample_size =  c(20, 40, 80),
    effect_size = c(0, 0.2, 0.5, 0.8),
    ntrials = c(2, 5)
  )


cl <- makeCluster(ncores)
registerDoParallel(cl)
set.seed(1016)
results_ntrials2_5 <- foreach(
  sample_size = conditions$sample_size,
  effect_size = conditions$effect_size,
  ntrials = conditions$ntrials,
  .export = names(.GlobalEnv),
  .packages = c("lavaan", "MASS", "tidyverse"),
  .combine = "bind_rows"
) %dorng% {
  run_sim(sample_size, effect_size, ntrials)
}
stopCluster(cl)

## for Mac: Multi cores
# results_multi <- do.call(mcmapply,
# c(
# FUN = run_sim,
# conditions,
# mc.cores = ncores,
# SIMPLIFY = FALSE
# ))
# allresults <- results_multi %>% bind_rows()


## SINGLE CORE
# results_single_ntrials2_5 <-
#     apply(conditions, 1, function(cond)
#         run_sim(
#             sample_size = cond[1],
#             effect_size = cond[2],
#             ntrials = cond[3]
#         ))
# results_single_ntrials2_5 <- results_single_ntrials2_5 %>% bind_rows()
##


#### ntrials = 18 divided averagely into 6 groups, thus in output 6
group <- 6
groupntrials <- 18 / group

get_data <- function(sample_size, effect_size, ntrials) {
  mu <-
    c(1, 1 + effect_size) # due to the specific variance and corvariance used in this simulation
  vcov <- matrix(0.5, nrow = 2, ncol = 2)
  diag(vcov) <- 1
  resid <- sqrt(2)
  
  samp <- mvrnorm(n = sample_size, mu = mu, Sigma = vcov)
  samp <- as.data.frame(samp)
  names(samp) <- c("true_compatible", "true_incompatible")
  
  for (ntrial in 1:ntrials) {
    groupsum <- 0
    for (groupntrial in 1:groupntrials) {
      groupsum <-
        groupsum + samp$true_compatible + rnorm(sample_size, sd = resid)
    }
    samp <- cbind(samp, groupsum / groupntrials)
    names(samp)[ncol(samp)] <- paste0("Y", ntrial, "_comp")
  }
  
  for (ntrial in 1:ntrials) {
    groupsum <- 0
    for (groupntrial in 1:groupntrials) {
      groupsum <-
        groupsum + samp$true_incompatible + rnorm(sample_size, sd = resid)
    }
    samp <- cbind(samp, groupsum / 3)
    names(samp)[ncol(samp)] <- paste0("Y", ntrial, "_incomp")
  }
  
  samp$sum_compatible <- rowMeans(samp[, 3:(3 + ntrials - 1)])
  samp$sum_incompatible <-
    rowMeans(samp[, (3 + ntrials):(3 + 2 * ntrials - 1)])
  samp
}

conditions <-
  expand.grid(
    sample_size =  c(20, 40, 80),
    effect_size = c(0, 0.2, 0.5, 0.8),
    ntrials = group
  )


cl <- makeCluster(ncores)
registerDoParallel(cl)
set.seed(108)
results_ntrials18 <- foreach(
  sample_size = conditions$sample_size,
  effect_size = conditions$effect_size,
  ntrials = conditions$ntrials,
  .export = names(.GlobalEnv),
  .packages = c("lavaan", "MASS", "tidyverse"),
  .combine = "bind_rows"
) %dorng% {
  run_sim(sample_size, effect_size, ntrials)
}
stopCluster(cl)


## for Mac: Multi cores
# results_multi <- do.call(mcmapply,
# c(
# FUN = run_sim,
# conditions,
# mc.cores = ncores,
# SIMPLIFY = FALSE
# ))
# results_multi <- results_multi %>% bind_rows()
# allresults <- bind_rows(allresults, results_multi)


## SINGLE CORE
# results_single_ntrials18 <-
#     apply(conditions, 1, function(cond)
#         run_sim(
#             sample_size = cond[1],
#             effect_size = cond[2],
#             ntrials = cond[3]
#         ))
# results_single_ntrials18 <- results_single_ntrials18 %>% bind_rows()
# ##

allresults <- bind_rows(results_ntrials2_5, results_ntrials18)

library(feather)
write_feather(allresults, "allresults.feather")


#### data analysis
error_improper_results <- allresults %>%
        group_by(sample_size, effect_size, ntrials, method) %>%
        summarise(
                error_improper_n = n(),
                error_rate = mean(grepl(
                        "covariance matrix of latent variables", msg
                )),
                improper_rate = mean(grepl("some estimated lv", msg)),
                d_missing_rate = mean(is.na(cohens_d)),
                rest_improper_rate = improper_rate - d_missing_rate
        )

## exclude samples with errors and improper solutions
wronglines <- !is.na(allresults$msg)
allcleaned <-
        anti_join(allresults,
                  allresults[wronglines,],
                  by = c("sample_size", "effect_size", "ntrials", "nsim"))

## exclude samples with errors
errorlines <-
        grepl("covariance matrix of latent variables", allresults$msg)
errorcleaned <-
        anti_join(allresults,
                  allresults[errorlines, ],
                  by = c("sample_size", "effect_size", "ntrials", "nsim"))

## exclude samples with improper solutions
improperlines <- grepl("some estimated lv", allresults$msg)
impropercleaned <-
        anti_join(allresults,
                  allresults[improperlines, ],
                  by = c("sample_size", "effect_size", "ntrials", "nsim"))

tidyup <- function(df) {
        df[df$ntrials == 2, ]$ntrials = 0.5
        df[df$ntrials == 5, ]$ntrials = 0.71
        df[df$ntrials == 6, ]$ntrials = 0.9
        for (spalte in c("sample_size", "effect_size", "ntrials", "method")) {
                df[[spalte]] <-
                        as.factor(df[[spalte]])
        }
        df
}

analyze <- function(df, type) {
        ## rate of p < .05
        p_results <- df %>%
                group_by(sample_size, effect_size, ntrials, method) %>%
                summarise(p_rate = mean(p_value < 0.05), p_n = n())
        
        ## exclude samples with missing est. cohen's d
        dmissinglines <- is.na(df$cohens_d)
        dcleaned <- anti_join(df,
                              df[dmissinglines,],
                              by = c("sample_size", "effect_size", "ntrials", "nsim"))
        
        ## deviation of median of est. cohen's d / variance of est. cohen's d / MAD / RMSE
        d_results <- dcleaned %>%
                group_by(sample_size, effect_size, ntrials, method) %>%
                summarise(
                        cohens_d_deviation = median(cohens_d - effect_size),
                        cohens_d_sd = sd(cohens_d),
                        MAD = mad(cohens_d - effect_size),
                        #RMSE = sqrt(mean((cohens_d - effect_size) ^ 2)),
                        d_n = n()
                )
        
        compareresults <-
                Reduce(
                        function(...)
                                full_join(
                                        ...,
                                        by = c("sample_size", "effect_size", "ntrials", "method")
                                ),
                        list(p_results,
                             d_results,
                             error_improper_results)
                ) %>%
                tidyup()
        
        write_feather(compareresults,
                      paste(type, "compareresults.feather"))
}

analyze(errorcleaned, "errorcleaned")
analyze(impropercleaned, "impropercleaned")
analyze(allcleaned, "allcleaned")
analyze(allresults, "all")


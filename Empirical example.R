library(tidyverse)
library(lavaan)
library(foreign)

dat_raw <-
  read.spss("WS19empra.sav", # please contact the me for the data set
            use.value.labels = T,
            to.data.frame = T)

dat <- dat_raw %>%
  group_by(vpn, avatar_rot, spos_ou, rpos) %>%
  mutate(com = (avatar_rot < 0 &
                  spos_ou != rpos) |
           (avatar_rot > 0 & spos_ou == rpos)) %>%
  mutate(com = factor(
    com,
    levels = c(T, F),
    labels = c("compatible", "incompatible")
  ))

set.seed(1016)
dat_numbered <- dat %>%
  ungroup() %>%
  dplyr::select(vpn, rt, com) %>%
  group_by(vpn, com) %>%
  mutate(trial = sample(1:n()))

## the least number of vaild trials for a single VP & condition: 120
# validtrials <- dat_numbered %>%
# group_by(vpn, com) %>%
# summarise(valid = max(trial))
# mintrials <- min(validtrials$valid)

dat_120_long <- dat_numbered %>%
  ungroup() %>%
  filter(trial <= 120)

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

# 120 trials divided averagely into 6 groups
group <- 6

dat_grouped <- dat_120_long %>%
  ungroup() %>%
  mutate(com = recode(com, compatible = "comp", incompatible = "incomp")) %>%
  group_by(vpn, com) %>%
  mutate(trial = (trial %% group) + 1) %>%
  group_by(vpn, com, trial) %>%
  summarise(rt = mean(rt)) %>%
  pivot_wider(
    names_from = c("trial", "com"),
    values_from = "rt",
    names_prefix = "Y",
    names_sep = "_"
  )

dat_grouped$sum_compatible <-
  rowMeans(dat_grouped[, 2:(2 + group - 1)])
dat_grouped$sum_incompatible <-
  rowMeans(dat_grouped[, (2 + group):(2 + group * 2 - 1)])

model <- create_model(group)

fit <- sem(model, data = dat_grouped)


ldsoutcome <- function(fit) {
  mu_diff =
    parTable(fit)$est[parTable(fit)$label == "mu_diff"]
  se_mu_diff =
    parTable(fit)$se[parTable(fit)$label == "mu_diff"]
  var_diff =
    parTable(fit)$est[parTable(fit)$label == "var_diff"]
  p_value = (1 - pnorm(abs(mu_diff / se_mu_diff))) * 2
  cohens_d = mu_diff / sqrt(var_diff)
  
  outcome <- data.frame(mu_diff,
                        se_mu_diff,
                        var_diff,
                        p_value,
                        cohens_d)
}

lds = ldsoutcome(fit)

## reliability of sum score
indicators <- dat_120_long %>%
  pivot_wider(
    id_cols = vpn,
    names_from = c("com", "trial"),
    values_from = "rt"
  ) %>%
  dplyr::select(-vpn)

library(Lambda4)
rel_sum_inc <- indicators %>%
  select(starts_with("inc")) %>%
  lambda3() %>%
  .$lambda3 %>%
  .$Unstandardized

rel_sum_com <- indicators %>%
  select(starts_with("com")) %>%
  lambda3() %>%
  .$lambda3 %>%
  .$Unstandardized

toutcome <- function(dataset) {
  ttest <-
    t.test(dataset$sum_compatible, dataset$sum_incompatible, paired = T)
  
  mu_diff = mean(dataset$sum_incompatible - dataset$sum_compatible)
  var_diff = var(dataset$sum_incompatible - dataset$sum_compatible)
  p_value = ttest$p.value
  cohens_d = mu_diff / sqrt(var_diff)
  t = ttest$statistic
  
  outcome <- tibble(mu_diff,
                    var_diff,
                    p_value,
                    cohens_d,
                    t)
}

outcome120 <- toutcome(dat_grouped)

outcomeall <- dat %>%
  ungroup() %>%
  group_by(vpn, com) %>%
  summarize(mRT = mean(rt)) %>%
  pivot_wider(names_from = com,
              names_prefix = "sum_",
              values_from = mRT) %>%
  toutcome()

for (obj in c("outcomeall",
              "outcome120",
              "lds",
              "fit",
              "rel_sum_inc",
              "rel_sum_com")) {
  cat('\n---', obj, '---\n')
  print(get(obj))
}

fitMeasures(fit)

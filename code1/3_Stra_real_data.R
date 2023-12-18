library(tidyverse)
library(sandwich)
library(haven)
set.seed(2021)
###############################################################################
##                        Data Wrangling                                     ##
## Y: outcome, X: pretreatment covariates B: stratumn indicator Z: treatment ##
###############################################################################
dt_raw <- read_dta("data/ok/opportunity knocks replication files/OKgradesUpdate_Feb5_2010 anonymized.dta")
dt_pre <- dt_raw %>%
  filter(s_first_year == 0) %>% 
  dplyr::select(
    # stratum id
    s_group_quart,
    # treatment
    `T`,
    # outcome
    avggradefall2008, 
    # stratum
    s_male, 
    s_quartile,
    # cont X
    s_hsgrade3,
    s_gpapreviousyear,
    s_age,
    s_mtongue_english,
    s_liveathome,
    s_highfundsconcern,
  ) %>%  mutate(
    B = structure(as.integer(as.factor(s_group_quart)), "label" = "stratum id"),
    Z = structure(as.integer(`T`), "label" = "treatment"),
    .keep = "unused",
    .before = everything()
  ) %>%  filter(!is.na(avggradefall2008))
X <- dt_pre %>% dplyr::select(s_hsgrade3:s_highfundsconcern) %>% mutate_all(~{replace_na(.,mean(.,na.rm = TRUE))}) %>% as.matrix()
Y <- dt_pre$avggradefall2008;Z <- dt_pre$Z;B <- dt_pre$B
w <- tibble(B,Z) %>% group_by(B) %>% mutate(nk=n(),n1k=sum(Z),
                  n0k = nk-n1k,
                  w=nk^2/(Z*n1k+(1-Z)*n0k)^2*(Z*n1k/(n1k-1)+(1-Z)*n0k/(n0k-1))) %>% pull(w) 
H <- B %>% unique() %>% length(); k <- ncol(X); n <- nrow(X)

###############################################################################
##                        Regression Analysis                                ##
###############################################################################

# compute the design matrix 
design_BZ <- map2_dfc(.x = rep(c(0,1),times= H), .y = rep(c(1:H),each= 2), .f= ~{
  ((B==.y)&(Z==.x)) %>% as.numeric()
}) %>% as.matrix()

  # regression with no interaction  ------------
  lm_tom <- lm(Y~X+design_BZ-1, weights = w)
  bhat <- lm_tom %>% coef() %>% .[1:k]
  mean_vec <- lm_tom %>% coef() %>% .[-(1:k)]
  pik <- table(B)/n
  est_coef <- rep(pik,each=2)
  est_coef <- ifelse(c(1:(2*H))%%2==1, -1,1)*est_coef
  htau_str_tom <- (est_coef*mean_vec) %>% sum()
  V_hw_tom <- est_coef%*% vcovHC(lm_tom,type = 'HC2')[-(1:k),-(1:k)]%*%est_coef
  # ToM regression adjusted estimator and endpoints of the confidence interval----------------
  htau_str_tom
  htau_str_tom+sqrt(V_hw_tom)*qnorm(0.025)
  htau_str_tom+sqrt(V_hw_tom)*qnorm(0.975)
  
  # Numerical results in Table 1--------------
  htaux_str <- X %>% as_tibble() %>% bind_cols(B=B,Z=Z) %>% group_by(B,Z) %>% summarise_all(mean) %>% 
    ungroup() %>% dplyr::select(-c(B,Z)) %>% summarise_all(~{sum(.*est_coef)}) %>% as.matrix()
  htaux_str
  bhat
  
  # naive estimator ----------------------
  lm_naive <- lm(Y~design_BZ-1)
  mean_vec <- lm_naive %>% coef() 
  htau_str_naive <- (est_coef*mean_vec) %>% sum()
  V_naive <- est_coef%*% vcovHC(lm_naive,type = 'HC2')%*%est_coef
  
  # unadjusted estimator and endpoints of the confidence intervals--------------
  htau_str_naive
  htau_str_naive+sqrt(V_naive)*qnorm(0.025)
  htau_str_naive+sqrt(V_naive)*qnorm(0.975)
  
 

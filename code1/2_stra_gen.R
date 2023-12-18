rm(list = ls())
library(tidyverse)
library(MASS)
library(doParallel)
library(sandwich)
registerDoParallel(48)

####################################################################################
##                               1.   Setting Global Parameter                    ##
####################################################################################
# strat='MS+FL' or 'MS' or 'FL', num_rep: number of replications
strata <- 'MS+FL'; delta <- 0.25; num_rep <- 1000


####################################################################################
##             2. Parallel Computing for Randomized  Experment                    ##
####################################################################################
sim_res <- foreach(seed=1:48,.combine='bind_rows',.packages = c('MASS','tidyverse','sandwich')) %dopar% {
  set.seed(seed)
  
  # s: number of covariates
  setups <- expand_grid(
    s     = seq(1,30,by=4),
    SNR1 = c(0.25,0.5,1,2),
    SNR0 = c(0.25,0.5,1,2)
  ) %>% mutate(case = row_number())
  
  
  #setups <- expand_grid(
  #  s     = 25,
  #  SNR1 = 0.25,
  #  SNR0 = 0.25
  #) %>% mutate(case = row_number())
  beta_vec <- rt(30,3); Delta_vec <- rt(30,3)
  map_dfr(1:nrow(setups), ~{
    
    par <- setups %>% slice(.x)
    attach(par)
    ## 1 generating finite sample-----------------------------------------------------------
    
    # many small blocks----------------------------------------------------
    # K: number of strata
    
    if(strata=='MS'){
      K <- 20
      n <- ceiling(runif(K,10,20))
    }
    
    # hybrid blocks----------------------------------------------------
    if(strata=='MS+FL'){
      K<- 12
      n <- c(ceiling(runif(2,140,160)),ceiling(runif(10,10,20)))
    }
    
    # a few large blocks---------------------------------------------
    if(strata=='FL'){
      K<- 2
      n <- ceiling(runif(K,140,160))
    }
    
    beta <- beta_vec[1:s]; Delta <- Delta_vec[1:s]
    beta1 <- beta+Delta*delta;  beta0 <- beta-Delta*delta
    
    # Generating dataset
    p <- rbeta(K,4,5); rho <- 0.4; N <- sum(n)
    X <- mvrnorm(N, rep(0,times=s),Sigma=rho + (1-rho)*diag(s))
    colnames(X) <- paste('X',1:s,sep = '')
    random_effect_beta1 <- map_dfr(1:K,~{
      mvrnorm(n=1,mu=rep(0,times = s),Sigma=diag(s)) %>% matrix(nrow=1) %>% .[rep(1,times=n[.x]),] %>% data.frame()
    })
    random_effect_mu1 <- rep(rnorm(K),times=n)
    random_effect_beta0 <- map_dfr(1:K,~{
      mvrnorm(n=1,mu=rep(0,times = s),Sigma=diag(s)) %>% matrix(nrow=1) %>% .[rep(1,times=n[.x]),] %>% data.frame()
    })
    random_effect_mu0 <- rep(rnorm(K),times=n)
    mu1 <- rep(rnorm(K, mean=3, sd=1),times=n)
    mu0 <- rep(rnorm(K, mean=2, sd=1),times=n) 
    Y1 <- rt(1,3)+random_effect_mu1+rowSums(X*random_effect_beta1) + X%*%beta1
    Y1 <- Y1+ rnorm(N)*sd(Y1)/sqrt(SNR1)
    Y0 <- rt(1,3)+random_effect_mu0 + rowSums(X*random_effect_beta0) + X%*%beta0
    Y0 <- Y0 + rnorm(N)*sd(Y0)/sqrt(SNR0)
    tau <- mean(Y1)-mean(Y0)
    # B: strata indicator
    B <- rep(1:K, times = n)
    n0 <- (n*(1-p)) %>% ceiling() %>% pmax(2) %>% pmin(n-2); n1 <- n-n0
    pik <- n/N
    
    sigma2_unadj <- map(1:K,~{
      Y1b <- Y1[(B==.x),]
      Y0b <- Y0[(B==.x),]
      n1k <- n1[.x]
      n0k <- n0[.x]
      nk <- n[.x]
      (var(Y1b)*nk/n1k + var(Y0b)*nk/n0k - var(Y0b-Y1b))*pik[.x]
    }) %>% Reduce('+',.)/N
    
    Sxx <- map(1:K,~{
      Xb <- X[(B==.x),]
      n1k <- n1[.x]
      n0k <- n0[.x]
      nk <- n[.x]
      if(s==1){(var(Xb)*nk/n1k + var(Xb)*nk/n0k)*pik[.x]}
      else{(cov(Xb)*nk/n1k + cov(Xb)*nk/n0k)*pik[.x]}
    }) %>% Reduce('+',.)
    
    # 2. simulation -----------------------------------------------
    res_mat<- map_dfr(1:num_rep,~{
      Z <- map2(n,n1,~{
        z <- rep(0,times = .x)
        z[sample(.x,.y,replace = FALSE)] <- 1
        z
      }) %>% unlist()
      w <- tibble(B,Z) %>% group_by(B) %>% mutate(nk=n(),n1k=sum(Z),
                                                  n0k = nk-n1k,
                                                  w=nk^2/(Z*n1k+(1-Z)*n0k)^2*(Z*n1k/(n1k-1)+(1-Z)*n0k/(n0k-1))) %>% pull(w) 
      Y <- Y1*Z+Y0*(1-Z)
      
      
      # 2.1 tyranny of the minority-------------------------------------------------
      
      ##----------------------------------------------------------------
      ##tom regression
      design_BZ <- map2_dfc(.x = rep(c(0,1),times= K), .y = rep(c(1:K),each= 2), .f= ~{
        ((B==.y)&(Z==.x)) %>% as.numeric()
      }) %>% as.matrix()
      lm.pool <- lm(Y~X+design_BZ-1, weights = w)
      e <- residuals(lm.pool)
      mean_vec <- lm.pool %>% coef() %>% .[-(1:s)]
      est_coef <- rep(pik,each=2)
      est_coef <- ifelse(c(1:(2*K))%%2==1, -1,1)*est_coef
      tau_tyr <- (est_coef*mean_vec) %>% sum()
      V_hw <- est_coef%*% vcovHC(lm.pool,type = 'HC2')[-(1:s),-(1:s)]%*%est_coef
      ##----------------------------------------------------------------
      
      ##----------------------------------------------------------------
      # equivalent formula of tom regression
      # design_BZ <- tibble(1,Z) %>%
      #   bind_cols(map_dfc(2:K,  ~{
      #     ((B==.x)-pik[.x])*Z
      #   })) %>%
      #   bind_cols(map_dfc(2:K, ~{
      #     (B==.x)-pik[.x]
      #   } )) %>% as.matrix()
      # lm.pool <- lm(Y~X+design_BZ-1, weights = w)
      # tau_tyr <- lm.pool %>% coef() %>% .[s+2]
      # 
      # model_x <- model.matrix(lm.pool); crossXinv <-solve(t(model_x)%*%diag(w)%*%model_x)
      # Omega <- cross2(0:1,1:K) %>% lapply(., function(x){
      #   cl_j <- (B==x[[2]])&(Z==x[[1]]); e_j <- e[cl_j]; w_j <- w[cl_j]
      #   X_j <- model_x[cl_j,,drop=FALSE]
      #   t(X_j)%*%diag(w_j)%*%tcrossprod(e_j)%*%diag(w_j)%*%X_j*(N-1)/(N-s-2*K)*(2*K)/(2*K-1)
      # }) %>% Reduce('+',.)
      # leverage <- model_x%*%solve(t(model_x)%*%diag(w)%*%model_x)%*%t(model_x)%*%diag(w) %>% diag()
      # Omega <- t(model_x)%*%diag(w)%*%diag(e^2/(1-leverage))%*%diag(w)%*%model_x
      # crossXinv[s+2,]
      # 
      # vcovHC(lm.pool,type = 'HC3')[s+2,s+2]
      ##----------------------------------------------------------------
      
      # 2.2 plugin estimator----------------------------------------------
      
      Sxt_hat <- map(1:K,~{
        Y1b <- Y[(B==.x)&(Z==1),]
        Y0b <- Y[(B==.x)&(Z==0),]
        X1b <- X[(B==.x)&(Z==1),]
        X0b <- X[(B==.x)&(Z==0),]
        n1k <- n1[.x]
        n0k <- n0[.x]
        nk <- n[.x]
        (cov(X1b,Y1b)*nk/n1k + cov(X0b,Y0b)*nk/n0k)*pik[.x]
      }) %>% Reduce('+',.)
      
      Stt_hat <- map(1:K,~{
        Y1b <- Y[(B==.x)&(Z==1),]
        Y0b <- Y[(B==.x)&(Z==0),]
        n1k <- n1[.x]
        n0k <- n0[.x]
        nk <- n[.x]
        (var(Y1b)*nk/n1k + var(Y0b)*nk/n0k)*pik[.x]
      }) %>% Reduce('+',.)
      
      # plugin estimator of variance (2 variations) -------------------------------------------
      R2_hat <- (t(Sxt_hat)%*%solve(Sxx,Sxt_hat)/Stt_hat) %>% min(1) %>% max(0)
      V_plg <- (1-R2_hat)*Stt_hat/N
      
      b_plg <- solve(Sxx,Sxt_hat)
      Y_adj_plg <- Y-X%*%b_plg
      temp <- tibble(Y_adj_plg=Y_adj_plg,Z,B)  %>% group_by(B,Z) %>% 
        summarise_all(mean)
      temp1 <- temp %>% filter(Z==1); temp0 <- temp %>% filter(Z==0)
      tau1_plg <- (temp1 %>% pull(Y_adj_plg)*pik) %>% sum()
      tau0_plg <- (temp0 %>% pull(Y_adj_plg)*pik) %>% sum()
      tau_plg <- tau1_plg-tau0_plg
      
      
      # 2.3 coverage and interval length-------------------------------
      tau_plg_V_plg <- as.numeric(abs(tau-tau_plg)<sqrt(V_plg)*qnorm(0.975))
      tau_tyr_V_hw <- as.numeric(abs(tau-tau_tyr)<sqrt(V_hw)*qnorm(0.975))
      len_neyman <-   2*sqrt(Stt_hat/N)*qnorm(0.975)
      len_V_plg <- 2*sqrt(V_plg)*qnorm(0.975)
      len_V_hw <- 2*sqrt(V_hw)*qnorm(0.975)
      
      tibble( tau_tyr, tau_plg, 
              len_V_hw = len_V_hw/len_neyman, len_V_plg = len_V_plg/len_neyman,
              tau_plg_V_plg, tau_tyr_V_hw)
    }) 
    detach(par)
    
    
    res_mat %>% dplyr::select(-c(1:2)) %>% summarise_all(mean) %>% 
      bind_cols(res_mat %>% summarise(across(.col=1:2, list(relative_rmse = ~{sqrt(mean((.x-tau)^2))/sqrt(sigma2_unadj)}, 
                                                            relative_bias = ~{abs(mean(.x-tau))/sqrt(sigma2_unadj)}))))
  })%>% 
    bind_cols(setups) %>% bind_cols(seed=seed) %>%  mutate(s = as.factor(s))
} 

save(sim_res,strata, file = paste0('stra_',strata ,Sys.time(),'.Rdata' ))  
    
    
  
 
####################################################################################
##                           3. Plotting the Figures                              ##
####################################################################################
library(latex2exp)

  sim_res  %>% dplyr::select(c('tau_plg_V_plg', 'tau_tyr_V_hw', 'SNR0', 'SNR1', 's')) %>% rename(k = s) %>%
    group_by(k,SNR1,SNR0) %>% summarise_all(median) %>% mutate(k = as.numeric(as.character(k))) %>%
    pivot_longer(cols = 4:5,names_to = 'Method',values_to='Coverage') %>%
    mutate(Method = fct_relevel(Method, 'tau_tyr_V_hw', 'tau_plg_V_plg' )) %>%
    ggplot(.,aes(x=k,y=Coverage,shape = Method, linetype = Method)) + geom_point(size=2) + geom_line() +
    scale_shape_manual(values=c(1,3), labels = unname(TeX(c("$(\\hat{\\tau}^{tom},\\hat{V}_{HC2,str}^{tom})",
                                                           "$(\\hat{\\tau}^{plg},\\hat{V}_{str}^{plg})")))) +
    scale_linetype_manual(values=c(1,4), labels = unname(TeX(c("$(\\hat{\\tau}^{tom},\\hat{V}_{HC2,str}^{tom})",
                                                            "$(\\hat{\\tau}^{plg},\\hat{V}_{str}^{plg})")))) +
    facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
    geom_hline(yintercept = 0.95, linetype = 'dashed') + theme(text = element_text(size = 25))+

    ylab('Coverage Probability')
  
  
  ggsave(filename = paste0('str_cp_',strata,'.pdf'),units='in',height = 8.5, width = 14)
  
  sim_res %>% mutate(s = as.numeric(as.character(s))) %>% rename(k = s)  %>% dplyr::select(c('tau_tyr_relative_rmse',  'tau_plg_relative_rmse',,
                                                                                            'SNR0', 'SNR1', 'k')) %>%
    group_by(k,SNR1,SNR0) %>% summarise_all(median) %>%
    pivot_longer(cols = 4:5,names_to = 'Method',values_to='RMSE') %>%
    mutate(Method = fct_relevel(Method,  'tau_tyr_relative_rmse',  'tau_plg_relative_rmse')) %>%
    ggplot(.,aes(x=k,y=RMSE,shape = Method, linetype = Method)) + geom_point(size=2) + geom_line() +
    facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
    geom_hline(yintercept = 1, linetype = 'dashed') + theme(text = element_text(size = 25))+
    scale_shape_manual(values=c(1,3), labels = unname(TeX(c("$\\hat{\\tau}^{tom}",
                                                         "$\\hat{\\tau}^{plg}"))))+
    scale_linetype_manual(values=c(1,4), labels = unname(TeX(c("$\\hat{\\tau}^{tom}",
                                                            "$\\hat{\\tau}^{plg}"))))+
    ylab('Relative RMSE')
  
  ggsave(filename = paste0('str_relative_rmse_',strata,'.pdf'),units='in',height = 8.5, width = 14)
  
  sim_res%>% mutate(s = as.numeric(as.character(s))) %>% rename(k = s)  %>% dplyr::select(c('len_V_hw',  'len_V_plg',
                                                                                                'SNR0', 'SNR1', 'k')) %>%
    group_by(k,SNR1,SNR0) %>% summarise_all(median) %>%
    pivot_longer(cols = 4:5,names_to = 'Method',values_to='Len') %>%
    mutate(Method = fct_relevel(Method,  'len_V_hw',  'len_V_plg',)) %>%
    ggplot(.,aes(x=k,y=Len,shape = Method, linetype = Method)) + geom_point(size=2) + geom_line() + 
    facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
    geom_hline(yintercept = 1, linetype = 'dashed') + theme(text = element_text(size = 25))+
    scale_shape_manual(values=c(1,3),labels = unname(TeX(c("$\\hat{V}_{HC2,str}^{tom}",
                                                        "$\\hat{V}_{str}^{plg}"))))+
    scale_linetype_manual(values=c(1,4),labels = unname(TeX(c("$\\hat{V}_{HC2,str}^{tom}",
                                                        "$\\hat{V}_{str}^{plg}"))))+
    ylab('Relative Length')
  
  ggsave(filename = paste0('str_cilen_',strata,'.pdf'),units='in',height = 8.5, width = 14)



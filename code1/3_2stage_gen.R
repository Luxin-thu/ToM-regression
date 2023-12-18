rm(list=ls())
library(tidyverse)
library(doParallel)
library(sandwich)
library(ggplot2)
registerDoParallel(48)

####################################################################################
##                               1.   Setting Global Parameter                    ##
####################################################################################
#f: sampling proportion in the first stage, num_rep: number of replication of the experiments 
f <- 0.01; num_rep <- 1000; delta <- 0.25


####################################################################################
##             2. Parallel Computing for Randomized  Experment                    ##
####################################################################################
sim_res <- foreach(seed=1:48,.combine='bind_rows', .packages = c('tidyverse','sandwich')) %dopar% {
  set.seed(seed)
  setups <- expand_grid(
    k     = seq(2,17,3),
    SNR1 = c(0.25,0.5,1,2),
    SNR0 = c(0.25,0.5,1,2),
  ) %>% mutate(case = row_number())
  #setups <- expand_grid(
  #  k     = 2,
  #  SNR1 = c(0.5,2),
  #  SNR0 = c(0.5,2)
  #) %>% mutate(case = row_number())
  
  p1 <- 0.3; N <- 10000; n <- f*N; p0 <-  1-p1; n1 <- n*p1; n0 <- n*p0
  Xmat <- matrix(rnorm(20*N), nrow = N, ncol = 20)
  Xmat <- scale(Xmat, scale = FALSE)
  beta_vec <- rt(20,3); Delta_vec <- rt(20,3)
  beta1_vec <- beta_vec+Delta_vec*delta;  beta0_vec <- beta_vec-Delta_vec*delta
  beta1_vec <- beta1_vec[order(abs(beta_vec), decreasing = TRUE)]
  beta0_vec <- beta0_vec[order(abs(beta_vec), decreasing = TRUE)]
  
  map_dfr(1:nrow(setups), ~{
    # Generating finite population  
    par <- setups %>% slice(.x)
    attach(par)
    X <- Xmat[,1:k]; W <- X[,1:2]; beta1 <- beta1_vec[1:k]; beta0 <- beta0_vec[1:k] 
    Y1 <- rt(1,3) + X%*%beta1  + rnorm(N)*sd(X%*%beta1)/sqrt(SNR1)
    Y0 <- rt(1,3) + X%*%beta0  + rnorm(N)*sd(X%*%beta0)/sqrt(SNR0)
    tau <- mean(Y1)-mean(Y0); sigma2_unadj <- (var(Y1)/p1+var(Y0)/(1-p1)-f*var(Y1-Y0))/n
    # Generating assignments
    res_mat <- map_dfr(1:num_rep, ~{
      S <- sample(1:N, n,replace= FALSE); Xs <- X[S,] %>% scale(scale = FALSE)
      Ws <- W %>% scale(scale = FALSE) %>% .[S,]
      Y1s <- Y1[S]; Y0s <- Y0[S]; Z <- rep(0,times = n)
      Z[sample(n, size = n*p1, replace=FALSE)] <- 1
      Ys <- Y1s*Z + Y0s*(1-Z); w <- Z/p1^2+(1-Z)/p0^2
      
      # ToM regression-----------
      Zmp0 <- Z-p0
      lm_tom <- lm(Ys~1+Z+ Xs +Zmp0:Ws, weights = w)
      
      V_hw <- vcovHC(lm_tom,type = 'HC3')[2,2]
      tau_tom <- lm_tom %>% coef() %>% .[2]
      cover_tom <- as.numeric(abs(tau-tau_tom)<sqrt(V_hw)*qnorm(0.975))
      
      
      # Yang's method-----------------
      h_gamma1 <- solve(cov(Ws[Z==1,]),cov(Ws[Z==1,],Ys[Z==1]))
      h_gamma0 <- solve(cov(Ws[Z==0,]),cov(Ws[Z==0,],Ys[Z==0]))
      h_beta1 <- solve(cov(Xs[Z==1,]),cov(Xs[Z==1,],Ys[Z==1]))
      h_beta0 <- solve(cov(Xs[Z==0,]),cov(Xs[Z==0,],Ys[Z==0]))
      h_beta <- p0*h_beta1+p1*h_beta0; h_gamma <- h_gamma1-h_gamma0
      Ys_beta_gamma <- Ys-Xs%*%h_beta-Z*p1*Ws%*%h_gamma+(1-Z)*p0*Ws%*%h_gamma
      tau_beta_gamma <- lm(Ys_beta_gamma~1+Z) %>% coef() %>% .[2]
      V_beta_gamma <- (var(Ys_beta_gamma[Z==1])/p1+var(Ys_beta_gamma[Z==0])/p0)/n
      cover_beta_gamma <- as.numeric(abs(tau-tau_beta_gamma)<sqrt(V_beta_gamma)*qnorm(0.975))
      
      # unadjust-------------------
      neyman_unadj <- var(Ys[Z==1])/n1+var(Ys[Z==0])/n0
      
      tibble(tau_tom, tau_beta_gamma, 
             len_hw = sqrt(V_hw/neyman_unadj),  len_plg =  sqrt(V_beta_gamma/neyman_unadj), 
             cover_tom,  cover_beta_gamma)
    })
    detach(par)
    res_mat %>% dplyr::select(-c(1:2)) %>% summarise_all(mean) %>% 
      bind_cols(res_mat %>% summarise(across(.col=1:2, list(relative_rmse = ~{sqrt(mean((.x-tau)^2))/sqrt(sigma2_unadj)}))))
  })%>% 
    bind_cols(setups) %>% bind_cols(seed=seed)
} 

save(sim_res, file=paste0('crs_',Sys.time(),'.Rdata'))



####################################################################################
##                           3. Plotting the Figures                              ##
####################################################################################
library(latex2exp)

sim_res  %>% dplyr::select(c('cover_tom', 'cover_beta_gamma', 
                                                          'SNR0', 'SNR1', 'k')) %>%
  group_by(k,SNR1,SNR0) %>% summarise_all(median) %>%
  pivot_longer(cols = 4:5,names_to = 'Method',values_to='Coverage') %>%
  mutate(Method = fct_relevel(Method,  'cover_tom', 'cover_beta_gamma')) %>%
  ggplot(.,aes(x=k,y=Coverage,shape = Method, linetype = Method)) + geom_point(size=2) + geom_line() + 
  scale_shape_manual(values=c(1,3), labels = unname(TeX(c("$(\\hat{\\tau}_{crs}^{tom},\\hat{V}^{tom}_{HC3,crs})",
                                                       "$(\\hat{\\tau}_{crs}^{plg},\\hat{V}_{crs}^{plg})")))) +
  scale_linetype_manual(values=c(1,3), labels = unname(TeX(c("$(\\hat{\\tau}_{crs}^{tom},\\hat{V}^{tom}_{HC3,crs})",
                                                          "$(\\hat{\\tau}_{crs}^{plg},\\hat{V}_{crs}^{plg})")))) +
  facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
  geom_hline(yintercept = 0.95, linetype = 'dashed') + theme(text = element_text(size = 25))+
  
  ylab('Coverage Probability')

ggsave(filename = paste0('crs_cp','.pdf'),units='in',height = 8.5, width = 14)

sim_res  %>% dplyr::select(c('tau_tom_relative_rmse', 'tau_beta_gamma_relative_rmse',
                                                            'SNR0', 'SNR1', 'k')) %>%
  group_by(k,SNR1,SNR0) %>% summarise_all(median) %>%
  pivot_longer(cols = 4:5,names_to = 'Method',values_to='RMSE') %>%
  mutate(Method = fct_relevel(Method,  'tau_tom_relative_rmse', 'tau_beta_gamma_relative_rmse')) %>%
  ggplot(.,aes(x=k,y=RMSE,shape = Method,linetype = Method)) + geom_point(size=2) + geom_line() + 
  facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
  geom_hline(yintercept = 1, linetype = 'dashed') + theme(text = element_text(size = 25))+
  scale_shape_manual(values=c(1,3), labels = unname(TeX(c("$\\hat{\\tau}^{tom}_{crs}",
                                                       "$\\hat{\\tau}^{plg}_{crs}"))))+
  scale_linetype_manual(values=c(1,4), labels = unname(TeX(c("$\\hat{\\tau}^{tom}_{crs}",
                                                          "$\\hat{\\tau}^{plg}_{crs}"))))+
  ylab('Relative RMSE')

ggsave(filename = paste0('crs_relative_rmse','.pdf'),units='in',height = 8.5, width = 14)

sim_res  %>% dplyr::select(c('len_hw', 'len_plg',  'SNR0', 'SNR1', 'k')) %>%
  group_by(k,SNR1,SNR0) %>% summarise_all(median) %>%
  pivot_longer(cols = 4:5,names_to = 'Method',values_to='Len') %>%
  mutate(Method = fct_relevel(Method,  'len_hw', 'len_plg')) %>%
  ggplot(.,aes(x=k,y=Len,shape = Method, linetype = Method)) + geom_point(size=2) + geom_line() + 
  facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
  geom_hline(yintercept = 1, linetype = 'dashed') + theme(text = element_text(size = 25))+
  scale_shape_manual(values=c(1,3),labels = unname(TeX(c("$\\hat{V}_{HC3,crs}^{tom}",
                                                      "$\\hat{V}_{crs}^{plg}"))))+
  scale_linetype_manual(values=c(1,4),labels = unname(TeX(c("$\\hat{V}_{HC3,crs}^{tom}",
                                                      "$\\hat{V}_{crs}^{plg}"))))+
  ylab('Relative Length')

ggsave(filename = paste0('crs_cilen','.pdf'),units='in',height = 8.5, width = 14)


rm(list=ls())
library(MASS)
library(tidyverse)
library(doParallel)
library(sandwich)
registerDoParallel(48)

####################################################################################
##                               1.   Setting Global Parameter                    ##
####################################################################################
# r: treated proportion of clusters, num_rep: number of iteration of cluster randomized experiments
# m: number of clusters
r<-0.3; num_rep <- 1000;   m <- 50; delta <- 0.25


####################################################################################
##             2. Parallel Computing for Randomized  Experment                    ##
####################################################################################
sim_res <- foreach(seed=1:48,.combine='bind_rows', .packages = c('MASS','tidyverse','sandwich')) %dopar% {
  set.seed(seed)
  # setups <- expand_grid(
  #   k     = 9,
  #   SNR1 = 2,
  #   SNR0 = 2,
  # ) %>% mutate(case = row_number())
  
  setups <- expand_grid(
    k     = seq(1,9,by=2),
    SNR1 = c(0.25,0.5,1,2),
    SNR0 = c(0.25,0.5,1,2),
  ) %>% mutate(case = row_number())
  beta_vec <- rt(9,3); Delta_vec <- rt(9,3)
  map_dfr(1:nrow(setups), ~{
    par <- setups %>% slice(.x)
    attach(par)
    # generating finite sample--------------------------
    
    ni  <- sample(4:10, m, replace = TRUE);rho <- 0.4
    n <- sum(ni); Sigma <- rho + (1-rho)*diag(k)
    X <- mvrnorm(n, rep(0,times = k),Sigma = Sigma)
    X <- scale(X, scale = FALSE)
    r1 <- r; r0 <- 1-r; m1 <- m*r1; m0 <- m*r0
    beta <- beta_vec[1:k]; Delta <- Delta_vec[1:k]
    beta1 <- beta+Delta*delta;  beta0 <- beta-Delta*delta
    random_effect_beta1 <- map_dfr(1:m,~{
      mvrnorm(n=1,mu=rep(0,times = k),Sigma=0.25*diag(k)) %>% matrix(nrow=1) %>% .[rep(1,times=ni[.x]),] %>% data.frame()
    })
    random_effect_intercept1 <- rep(rnorm(m),ni)*0.5
    random_effect_beta0 <- map_dfr(1:m,~{
      mvrnorm(n=1,mu=rep(0,times = k),Sigma=0.25*diag(k)) %>% matrix(nrow=1) %>% .[rep(1,times=ni[.x]),] %>% data.frame()
    })
    random_effect_intercept0 <- rep(rnorm(m),ni)*0.5
    Y1 <- random_effect_intercept1 + rowSums(X*random_effect_beta1) + X%*%beta1
    Y1 <- Y1 + rt(1,3) + rnorm(n)*sd(Y1)/sqrt(SNR1)
    Y0 <- random_effect_intercept0 + rowSums(X*random_effect_beta0) + X%*%beta0
    Y0 <- Y0 + rt(1,3) + rnorm(n)*sd(Y0)/sqrt(SNR0)
    tau <- mean(Y1-Y0)
    cluster <- rep(1:m,ni)
    Y1tld <- tibble(Y1=Y1,cluster) %>% group_by(cluster) %>% summarise_all(sum) %>% pull(Y1)*m/n
    Y0tld <- tibble(Y0=Y0,cluster) %>% group_by(cluster) %>% summarise_all(sum) %>% pull(Y0)*m/n
    Xtld <- tibble(data.frame(X),cluster) %>% group_by(cluster) %>% summarise_all(sum) %>% 
      dplyr::select(-c(cluster)) %>% bind_cols(ni=ni-mean(ni)) %>% as.matrix()
    
    e1 <- lm(formula = Y1tld~1+Xtld) %>% resid()
    e0 <- lm(formula = Y0tld~1+Xtld) %>% resid()
    sigma2_adj <- (var(e1)/r+var(e0)/(1-r)-var(e1-e0))/m
    sigma2_unadj <- (var(Y1tld)/r+var(Y0tld)/(1-r)-var(Y1tld-Y0tld))/m
    sigma2_adj/sigma2_unadj
    # simulation----------------------------------
    res_mat <- map_dfr(1:num_rep,~{
      Z <- rep(0, times=m)
      Z[sample(m,m1,replace = FALSE)] <- 1
      Ytld <- Z*Y1tld+(1-Z)*Y0tld
      
      m_lin <- lm(formula = Ytld~Z+Xtld+Z:Xtld)
      tau_lin <- m_lin %>% coef() %>% .[2]
      
      w <- ifelse(Z==1, (r1)^(-2), (r0)^(-2))
      m_tyr <- lm(formula = Ytld~Z+Xtld, weights = w)
      tau_tyr <- m_tyr %>%  coef() %>% .[2]
      
      tau_naive <- lm(formula = Ytld~Z) %>% coef() %>% .[2]
      
      hw_lin <- vcovHC(m_lin,type = 'HC0')[2,2]
      
      hw_tyr <- vcovHC(m_tyr, type = 'HC0')[2,2]
      
      neyman_unadj <- var(Y1tld[Z==1])/m1+var(Y0tld[Z==0])/m0
      
      tau_lin_hw_lin <- as.numeric(abs(tau-tau_lin)<sqrt(hw_lin)*qnorm(0.975))
      tau_tyr_hw_lin <- as.numeric(abs(tau-tau_tyr)<sqrt(hw_lin)*qnorm(0.975))
      tau_tyr_hw_tyr <- as.numeric(abs(tau-tau_tyr)<sqrt(hw_tyr)*qnorm(0.975))
      len_hw_lin <- 2*sqrt(hw_lin)*qnorm(0.975)
      len_hw_tyr <- 2*sqrt(hw_tyr)*qnorm(0.975)
      len_neyman <- 2*sqrt(neyman_unadj)*qnorm(0.975)
      
      out <- tibble(tau_lin_hw_lin, tau_tyr_hw_lin, 
                    tau_tyr_hw_tyr, 
                    len_hw_lin = len_hw_lin/len_neyman, 
                    len_hw_tyr = len_hw_tyr/len_neyman, 
                    tau_lin, tau_tyr)
    })
    detach(par)
    res_mat %>% dplyr::select(-c(6:7)) %>% summarise_all(mean) %>% 
      bind_cols(res_mat %>% summarise(across(.col=6:7, list(relative_rmse = ~{sqrt(mean((.x-tau)^2))/sqrt(sigma2_unadj)}))))
  })%>% 
    bind_cols(setups) %>% bind_cols(seed=seed)
}


# saving data and plots ----------------------------------------
save(sim_res, file=paste0('cluster',Sys.time(),'.Rdata'))

####################################################################################
##                           3. Plotting the Figures                              ##
####################################################################################
library(latex2exp)
sim_res %>% mutate(k = as.integer(as.character(k)))  %>% dplyr::select(c('tau_lin_hw_lin', 'tau_tyr_hw_lin', 'tau_tyr_hw_tyr', 
                                                           'SNR0', 'SNR1', 'k')) %>%
  group_by(k,SNR1,SNR0) %>% summarise_all(median) %>%
  pivot_longer(cols = 4:6,names_to = 'Method',values_to='Coverage') %>%
  mutate(Method = fct_relevel(Method,  'tau_tyr_hw_tyr', 'tau_lin_hw_lin', 'tau_tyr_hw_lin')) %>%
  ggplot(.,aes(x=k,y=Coverage,shape = Method, linetype = Method)) + geom_point(size=2) + geom_line() + 
  scale_x_continuous(labels = function(x) as.integer(x), breaks = seq(1,9,by=2))+
  scale_shape_manual(values=c(1,2,4), labels = unname(TeX(c("$(\\hat{\\tau}^{tom}_{cl},\\hat{V}_{HC0,cl}^{tom})",
                                                       "$(\\hat{\\tau}^{lin}_{cl},\\hat{V}_{HC0,cl}^{lin})",
                                                       "$(\\hat{\\tau}^{tom}_{cl},\\hat{V}_{HC0,cl}^{lin})")))) +
  scale_linetype_manual(values=c(1,3,5), labels = unname(TeX(c("$(\\hat{\\tau}^{tom}_{cl},\\hat{V}_{HC0,cl}^{tom})",
                                                            "$(\\hat{\\tau}^{lin}_{cl},\\hat{V}_{HC0,cl}^{lin})",
                                                            "$(\\hat{\\tau}^{tom}_{cl},\\hat{V}_{HC0,cl}^{lin})")))) +
  facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
  geom_hline(yintercept = 0.95, linetype = 'dashed') + theme(text = element_text(size = 25))+
  ylab('Coverage Probability')

ggsave(filename = paste0('cl_cp','.pdf'),units='in',height = 8.5, width = 14)


sim_res %>% mutate(k = as.integer(as.character(k)))  %>% dplyr::select(c('tau_lin_relative_rmse', 'tau_tyr_relative_rmse',
                                                           'SNR0', 'SNR1', 'k')) %>%
  group_by(k,SNR1,SNR0) %>% summarise_all(median) %>%
  pivot_longer(cols = 4:5,names_to = 'Method',values_to='RMSE') %>%
  mutate(Method = fct_relevel(Method,  'tau_tyr_relative_rmse',  'tau_lin_relative_rmse')) %>%
  ggplot(.,aes(x=k,y=RMSE,shape = Method, linetype = Method)) + geom_point(size=2) + geom_line() + 
  scale_x_continuous(labels = function(x) as.integer(x), breaks = seq(1,9,by=2))+
  facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
  geom_hline(yintercept = 1, linetype = 'dashed') + theme(text = element_text(size = 25))+
  scale_shape_manual(values=1:2, labels = unname(TeX(c("$\\hat{\\tau}_{cl}^{tom}","$\\hat{\\tau}_{cl}^{lin}"))))+
  scale_linetype_manual(values=c(1,3), labels = unname(TeX(c("$\\hat{\\tau}_{cl}^{tom}","$\\hat{\\tau}_{cl}^{lin}"))))+
  ylab('Relative RMSE')
ggsave(filename = paste0('cl_relative_rmse','.pdf'),units='in',height = 8.5, width = 14)


sim_res %>% mutate(k = as.integer(k))  %>% dplyr::select(c('len_hw_lin', 'len_hw_tyr',  'SNR0', 'SNR1', 'k')) %>%
  group_by(k,SNR1,SNR0) %>% summarise_all(median) %>%
  pivot_longer(cols = 4:5,names_to = 'Method',values_to='Len') %>%
  mutate(Method = fct_relevel(Method,  'len_hw_tyr', 'len_hw_lin' )) %>%
  ggplot(.,aes(x=k,y=Len,shape = Method, linetype = Method)) + geom_point(size=2) + geom_line() + 
  scale_x_continuous(labels = function(x) as.integer(x), breaks = seq(1,9,by=2))+
  facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
  geom_hline(yintercept = 1, linetype = 'dashed') + theme(text = element_text(size = 25))+
  scale_shape_manual(values=1:2,labels = unname(TeX(c("$\\hat{V}_{HC0,cl}^{tom}","$\\hat{V}_{HC0,cl}^{lin}"))))+
  scale_linetype_manual(values=c(1,3),labels = unname(TeX(c("$\\hat{V}_{HC0,cl}^{tom}","$\\hat{V}_{HC0,cl}^{lin}"))))+
  ylab('Relative Length')
ggsave(filename = paste0('cl_cilen','.pdf'),units='in',height = 8.5, width = 14)




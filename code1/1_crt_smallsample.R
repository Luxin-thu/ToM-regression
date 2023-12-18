rm(list=ls())
library(MASS)
library(tidyverse)
library(doParallel)
library(sandwich)
num_cl <- 48
registerDoParallel(num_cl)

####################################################################################
##                              1.  Setting Global Parameter                      ##
####################################################################################

# setting n=20 or n=40, n:sample size, num_rep: number of replications
n <- 20; num_rep<-1000

# 3 settings r=0.3, r=0.4, r=0.5--------------
r<-0.3; r1 <- r; r0 <- 1-r; n1 <- n*r1; n0 <- n*r0; delta <- 0.25

####################################################################################
##           2. Parallel Computing for Randomized  Experment                      ##
####################################################################################

sim_res <- foreach(seed=1:100,.combine='bind_rows', .packages = c('MASS','tidyverse','sandwich','glmnet')) %dopar% {
  set.seed(seed)
  # Simulation setups
  setups <- expand_grid(
    k     = c(1,2,3),
    SNR1 = c(0.25,0.5,1,2),
    SNR0 = c(0.25,0.5,1,2)
  ) %>% mutate(case = row_number())
  
  # setups <- expand_grid(
  #   k     = c(1,2,3,4,5,6),
  #   SNR1 = c(0.25,0.5,1,2),
  #   SNR0 = c(0.25,0.5,1,2)
  # ) %>% mutate(case = row_number())

    # setups <- expand_grid(
    #   k     = c(1,2),
    #   SNR1 = c(0.5),
    #   SNR0 = c(0.5)
    # ) %>% mutate(case = row_number())
  
  Xmat <- matrix(rnorm(20*n), nrow = n, ncol = 20)
  Xmat <- scale(Xmat, scale = FALSE)
  beta_vec <- rt(20,3); Delta_vec <- rt(20,3)
  map_dfr(1:nrow(setups), ~{
    par <- setups %>% slice(.x)
    attach(par)
    # Generating finite population--------------------------
    beta <- beta_vec[1:k]; Delta <- Delta_vec[1:k]
    beta1 <- beta+Delta*delta;  beta0 <- beta-Delta*delta
    X <- Xmat[,1:k]
    if(k>1){    
      Y1 <- (X%*%beta1)  + rnorm(n)*sd(X%*%beta1)/sqrt(SNR1)
      Y0 <- (X%*%beta0)  + rnorm(n)*sd(X%*%beta0)/sqrt(SNR0)}
    
    if(k==1){
      Y1 <- (X*beta1)  + rnorm(n)*sd(X*beta1)/sqrt(SNR1)
      Y0 <- (X*beta0)  + rnorm(n)*sd(X*beta0)/sqrt(SNR0)  
    }
    
    tau <- mean(Y1-Y0)
    
    e1 <- lm(formula = Y1~1+X) %>% resid()
    e0 <- lm(formula = Y0~1+X) %>% resid()
    sigma2_adj <- (var(e1)/r+var(e0)/(1-r)-var(e1-e0))/(n)
    sigma2_unadj <- (var(Y1)/r+var(Y0)/(1-r)-var(Y1-Y0))/n
    # Generating assignments----------------------------------
    res_mat <- map_dfr(1:num_rep,~{
      Z <- rep(0, times=n)
      Z[sample(n,n1,replace = FALSE)] <- 1
      Y <- Z*Y1+(1-Z)*Y0; w <- Z/r1^2+(1-Z)/r0^2
      
      m_tom <- lm(Y~Z+ X, weights = w)
      tau_tom <- m_tom %>% coef() %>% .[2]
      hw0_tom<- vcovHC(m_tom, type = 'HC0')[2,2]
      tau_tom_hw0_tom<- as.numeric(abs(tau-tau_tom)<sqrt(hw0_tom)*qnorm(0.975))
      
      m_lin <- lm(formula = Y~Z+X+Z:X)
      tau_lin <- m_lin %>% coef() %>% .[2]
      hw0_lin <- vcovHC(m_lin,type = 'HC0')[2,2]
      tau_lin_hw0_lin <- as.numeric(abs(tau-tau_lin)<sqrt(hw0_lin)*qnorm(0.975))
      
      tau_tom_hw0_lin <- as.numeric(abs(tau-tau_tom)<sqrt(hw0_lin)*qnorm(0.975))
      if(k>1){h_beta <- solve(cov(X),cov(X[Z==1,],Y[Z==1])*r0+cov(X[Z==0,],Y[Z==0])*r1)
      tau_plgin <- mean((Y-X%*%h_beta)[Z==1])-mean((Y-X%*%h_beta)[Z==0])
      V_plgin <- var((Y-X%*%h_beta)[Z==1])/n1+var((Y-X%*%h_beta)[Z==0])/n0
      }else{
        h_beta <- (cov(X[Z==1],Y[Z==1])*r0+cov(X[Z==0],Y[Z==0])*r1)/var(X)
        tau_plgin <- mean((Y-X*h_beta)[Z==1])-mean((Y-X*h_beta)[Z==0])
        V_plgin <- var((Y-X*h_beta)[Z==1])/n1+var((Y-X*h_beta)[Z==0])/n0
      }
      
      tau_plgin_v_plgin <- as.numeric(abs(tau-tau_plgin)<sqrt(V_plgin)*qnorm(0.975))
      # Confidence interval length----------------------  
      neyman_unadj <- var(Y1[Z==1])/n1+var(Y0[Z==0])/n0
      len_hw0_lin <- 2*sqrt(hw0_lin)*qnorm(0.975)
      len_hw0_tom<- 2*sqrt(hw0_tom)*qnorm(0.975)
      len_neyman <- 2*sqrt(neyman_unadj)*qnorm(0.975)
      len_plgin <- 2*sqrt(V_plgin)*qnorm(0.975)
      
      out <- tibble(tau_lin_hw0_lin,  tau_tom_hw0_lin, 
                    tau_tom_hw0_tom,  tau_plgin_v_plgin,
                    len_hw0_lin = len_hw0_lin/len_neyman, 
                    len_hw0_tom = len_hw0_tom/len_neyman, 
                    len_plgin = len_plgin/len_neyman,
                    tau_lin,  tau_tom, tau_plgin)
    })
    detach(par)
    res_mat %>% dplyr::select(-c(8:10)) %>% summarise_all(mean) %>% 
      
      bind_cols(res_mat %>% summarise(across(.col=8:10, list(relative_rmse = ~{sqrt(mean((.x-tau)^2))/sqrt(sigma2_unadj)}, 
                                                             relative_bias = ~{abs(mean(.x-tau))/sqrt(sigma2_adj)}))))
  })%>% 
    bind_cols(setups) %>% bind_cols(seed=seed)
} 

# saving data----------------------------------------
save(sim_res,r, file=paste0('crt_n_',n,Sys.time(),'.Rdata'))



####################################################################################
##                        3.  Plotting the Figures                                ##
####################################################################################
library(latex2exp)

sim_res %>% mutate(k = as.integer(k))  %>% dplyr::select(c('tau_tom_hw0_tom', 'tau_tom_hw0_lin', 'tau_lin_hw0_lin', 
                                                          'SNR0', 'SNR1', 'k')) %>%
  group_by(k,SNR1,SNR0) %>% summarise_all(median) %>%
  pivot_longer(cols = 4:6,names_to = 'Method',values_to='Coverage') %>%
  mutate(Method = fct_relevel(Method,  'tau_tom_hw0_tom', 'tau_lin_hw0_lin','tau_tom_hw0_lin')) %>%
  ggplot(.,aes(x=k,y=Coverage,shape = Method, linetype = Method)) + geom_point(size=2) + geom_line() + 
  # scale_x_continuous(labels = function(x) as.integer(x), breaks = c(1,2,3))+
  scale_shape_manual(values=c(1,2,4), labels = unname(TeX(c("$(\\hat{\\tau}^{tom},\\hat{V}_{HC0}^{tom})",
                                                       "$(\\hat{\\tau}^{lin},\\hat{V}_{HC0}^{lin})",
                                                       "$(\\hat{\\tau}^{tom},\\hat{V}_{HC0}^{lin})")))) +
  scale_linetype_manual(values=c(1,3,5), labels = unname(TeX(c("$(\\hat{\\tau}^{tom},\\hat{V}_{HC0}^{tom})",
                                                       "$(\\hat{\\tau}^{lin},\\hat{V}_{HC0}^{lin})",
                                                       "$(\\hat{\\tau}^{tom},\\hat{V}_{HC0}^{lin})")))) +
  facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
  geom_hline(yintercept = 0.95, linetype = 'dashed') + theme(text = element_text(size = 25))+
  
  ylab('Coverage Probability')
   
ggsave(filename = paste0('crt_cp_n_',n,'.pdf'),units='in',height = 8.5, width = 14)


sim_res %>% mutate(k = as.integer(k))  %>% dplyr::select(c('tau_lin_relative_rmse', 'tau_tom_relative_rmse',
                                                            'SNR0', 'SNR1', 'k')) %>%
  group_by(k,SNR1,SNR0) %>% summarise_all(median) %>%
  pivot_longer(cols = 4:5,names_to = 'Method',values_to='RMSE') %>%
  mutate(Method = fct_relevel(Method,  'tau_tom_relative_rmse', 'tau_lin_relative_rmse')) %>%
  ggplot(.,aes(x=k,y=RMSE,shape = Method, linetype = Method)) + geom_point(size=2) + geom_line() + 
  # scale_x_continuous(labels = function(x) as.integer(x), breaks = c(1,2,3))+
  facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
  geom_hline(yintercept = 1, linetype = 'dashed') + theme(text = element_text(size = 25))+
  scale_shape_manual(values=1:2, labels = unname(TeX(c("$\\hat{\\tau}^{tom}",
                                                      "$\\hat{\\tau}^{lin}"))))+
  scale_linetype_manual(values=c(1,3), labels = unname(TeX(c("$\\hat{\\tau}^{tom}",
                                                       "$\\hat{\\tau}^{lin}"))))+
  ylab('Relative RMSE')

ggsave(filename = paste0('crt_relative_rmse_n_',n,'.pdf'),units='in',height = 8.5, width = 14)


sim_res %>% mutate(k = as.integer(k))  %>% dplyr::select(c('len_hw0_lin', 'len_hw0_tom',  'SNR0', 'SNR1', 'k')) %>%
  group_by(k,SNR1,SNR0) %>% summarise_all(median) %>%
  pivot_longer(cols = 4:5,names_to = 'Method',values_to='Len') %>%
  mutate(Method = fct_relevel(Method,  'len_hw0_tom', 'len_hw0_lin')) %>%
  ggplot(.,aes(x=k,y=Len,shape = Method, linetype = Method)) + geom_point(size=2) + geom_line() + 
  # scale_x_continuous(labels = function(x) as.integer(x), breaks = c(1,2,3))+
  facet_grid(SNR0~SNR1,labeller = label_both)+ theme(legend.position="bottom")+
  geom_hline(yintercept = 1, linetype = 'dashed') + theme(text = element_text(size = 25))+
  scale_shape_manual(values=1:2,labels = unname(TeX(c("$\\hat{V}_{HC0}^{tom}",
                                                  "$\\hat{V}_{HC0}^{lin}"))))+
  scale_linetype_manual(values=c(1,3),labels = unname(TeX(c("$\\hat{V}_{HC0}^{tom}","$\\hat{V}_{HC0}^{lin}"
                                                      ))))+
  ylab('Relative Length')

ggsave(filename = paste0('crt_cilen_n_',n,'.pdf'),units='in',height = 8.5, width = 14)

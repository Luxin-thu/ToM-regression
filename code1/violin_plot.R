rm(list=ls())
library(tidyverse)
library(gridExtra)
library(MASS)
set.seed(1)

leverage <- function(X){
  return(diag(X%*%solve(t(X)%*%X/n,t(X))/n))
}

#################################################################
##                  Setting Global Parameters                  ##
#################################################################

SNR1<- 0.25; SNR0 <- 2; q <- 1; r=0.3; n <- 100
rho = 0.4; delta <- 0.25
Sigma <- rho + (1-rho)*diag(20)
Xmat <- mvrnorm(n, rep(0,times = 20),Sigma = Sigma)
Xmat <- scale(Xmat, scale = FALSE)
beta_vec <- rt(20,3); Delta_vec <- rt(20,3)

#################################################################
##                  Simulating the Experiment                  ##
#################################################################

df_plotting <- map(c(2,8,14,20),~{
  k=.x
  beta <- beta_vec[1:k]; Delta <- Delta_vec[1:k]
  beta1 <- beta+Delta*delta;  beta0 <- beta-Delta*delta
  X <- Xmat[,1:k]
  Y1 <- scale(X%*%beta1)  + rnorm(n)/sqrt(SNR1)
  Y0 <- scale(X%*%beta0)  + rnorm(n)/sqrt(SNR0)
  
  r1 <- r; r0 <- 1-r; n1 <- n*r1; n0 <- n*r0
  Z <- rep(0, times=n)
  Z[sample(n,n1,replace = FALSE)] <- 1
  Y <- Z*Y1+(1-Z)*Y0
  
  ## ToM ------------------------------------
  w <- ifelse(Z==1, (n/n1)^2, (n/n0)^2)
  m_tom<- lm(formula = Y~Z+X, weights = w)
  
  sxt <- cov(X[Z==1,]); sxc <- cov(X[Z==0,]); hbxt <- colMeans(X[Z==1,]); hbxc <- colMeans(X[Z==0,])
  htau_x <-colMeans(X[Z==1,])-colMeans(X[Z==0,])
  
  #c_tom calibrated weights, h_tom leverage scores, Fc_tom influence function
  c_tom <- (1/n1-t(htau_x)%*%solve((n0-1)*sxc/r0^2+(n1-1)*sxt/r1^2)%*%(t(X)-hbxt)/r1^2)*Z+
    (1/n0+t(htau_x)%*%solve((n0-1)*sxc/r0^2+(n1-1)*sxt/r1^2)%*%(t(X)-hbxc)/r0^2)*(1-Z)
  h_tom <- leverage(cbind(1,Z,X)*w^{1/2}); IF_tom = drop(n*(2*Z-1)*(c_tom*resid(m_tom)/(1-h_tom)))
  Fc_tom <- sum((c_tom[Z==1]*n1-1)^2)/2+ sum((c_tom[Z==0]*n0-1)^2)/2
  
  ## Lin --------------------------------
  
  m_lin <- lm(formula = Y~Z+X+Z:X)
  
  sxt <- cov(X[Z==1,]); sxc <- cov(X[Z==0,]); hbxt <- colMeans(X[Z==1,]); hbxc <- colMeans(X[Z==0,])
  htau_x <-colMeans(X[Z==1,])-colMeans(X[Z==0,])
  
  #c_lin calibrated weights, h_lin leverage scores, Fc_lin influence function
  c_lin <- (1/n1-t(htau_x)%*%solve((n1-1)*sxt)%*%(t(X)-hbxt)*r0)*Z+
    (1/n0+t(htau_x)%*%solve((n0-1)*sxc)%*%(t(X)-hbxc)*r1)*(1-Z)
  h_lin <- leverage(cbind(1,Z,X,X*Z)); IF_lin = drop(n*(2*Z-1)*(c_lin*resid(m_lin)/(1-h_lin)))
  Fc_lin <- sum((c_lin[Z==1]*n1-1)^2)/2+ sum((c_lin[Z==0]*n0-1)^2)/2
  tibble(c_lin = drop(c_lin), h_lin = drop(h_lin), IF_lin = drop(IF_lin),c_tom = drop(c_tom), h_tom = drop(h_tom), IF_tom = drop(IF_tom), k=as.factor(k))
}) %>% bind_rows()


##################################################################
##                     Plotting the Figures                     ##
##################################################################
p1 <- df_plotting %>% dplyr::select(c(IF_lin, IF_tom, k)) %>% pivot_longer(cols = 1:2,names_to = 'Method',values_to='IF') %>% 
  ggplot(aes(y=IF, color=Method, x = k))+geom_violin() + ylab('Sample influence curves') + theme(text = element_text(size = 15),legend.position = "none") + scale_color_manual(values=1:2, labels = c('Lin','ToM'))

p2 <- df_plotting %>% dplyr::select(c(h_lin, h_tom, k)) %>% pivot_longer(cols = 1:2,names_to = 'Method',values_to='hi') %>% 
  ggplot(aes(y=hi, color=Method, x = k))+geom_violin()+ ylab('Leverage scores') + theme(text = element_text(size = 15),legend.position = "none") + scale_color_manual(values=1:2, labels = c('Lin','ToM'))

p3 <- df_plotting %>% dplyr::select(c(c_lin, c_tom, k)) %>% pivot_longer(cols = 1:2,names_to = 'Method',values_to='weight') %>% 
  ggplot(aes(y=weight, color=Method, x = k))+geom_violin() + ylab('Calibrated weights') + theme(legend.position = c (0.02, 0.98), legend.justification = c (0, 1),text = element_text(size = 15))  + scale_color_manual(values=1:2, labels = c('Lin','ToM'))

plot <- grid.arrange(p3, p2,  p1, 
             ncol = 3)

ggsave(file='violin_plot_weights.pdf', plot,width = 8, height = 4.5)

  
  


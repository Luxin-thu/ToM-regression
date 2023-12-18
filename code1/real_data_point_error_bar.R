library('ggpubr')



data <- tibble(setting = as.factor(c('OK','OK','less-polarization vs control','less-polarization vs control', 
                               'more-polarization vs control', 'more-polarization vs control')),
                   est = c(0.305,0.157,0.0982,0.0740,-0.0545,-0.0561),
                   se = c(0.562, 0.727,0.0275,0.0297,0.0281,0.0294),
                   Method = c('ToM', 'unadj','ToM', 'unadj', 'ToM', 'unadj')
                   ) %>% mutate(setting=fct_relevel(setting, 'OK','less-polarization vs control', 'more-polarization vs control'))
data <- data %>% mutate(ci = se*qnorm(0.975),len=c('2.203','2.850','0.108','0.116','0.110','0.115'))



data %>% filter(!setting%in%c('OK')) %>%  ggplot(aes(x=Method,y=est))+geom_point()+ylab('Estimator')+geom_errorbar(aes(ymin=est-ci, ymax=est+ci),width=.1)+geom_text(aes(label = len, y = est+ci, x=Method), hjust = -1, , size=7)+facet_wrap(.~setting, scales='free')+theme(text = element_text(size = 20),
                                                                                                                                                                                                                                                            strip.text.x = element_text(size = 20),
                                                                                                                                                                                                                                                            axis.text.x =element_text(size = 20),
                                                                                                                                                                                                                                                              axis.title.y =element_text(size = 20))
ggsave(filename = paste0('est_error_bar','.pdf'),units='in',height = 6, width = 14)

library(tidyverse)
library(sandwich)
library(leaps)
library(MASS)
library(formula.tools)
####################################################################################
##Target population covariates from the the website of United States Census Bureau##
####################################################################################

## Target population population covariate: age and sex
age_sex_data <- read.csv('data/social trust data/nc-est2020-agesex-res.csv')
age_sex_data1 <- age_sex_data %>%  dplyr::filter(AGE!=999&SEX!=0&AGE>=18) %>% dplyr::select(c('SEX','AGE','POPESTIMATE2020')) %>% 
  mutate(weights = POPESTIMATE2020/sum(POPESTIMATE2020)) %>% mutate(SEX=as.numeric(SEX==1))
female_population <- (age_sex_data1$SEX*age_sex_data1$weights) %>% sum()
age_population <- (age_sex_data1$AGE*age_sex_data1$weights) %>% sum()

## Target population population covariate: ratio of race
race_data <- read.csv('data/social trust data/SC-EST2020-ALLDATA5.csv')
## Race1 white and non-hispanic; Race2  Black or African American, non-Hispanic; Race3 Hispanic; Race4  Other, non-Hispanic; Race5 Two or more races, non-Hispanic 
race_data1 <- race_data %>% dplyr::select(-c('NAME')) %>% dplyr::filter(ORIGIN!=0&SEX!=0&AGE>=18) %>% mutate(race1=as.numeric(ORIGIN==1&RACE==1),
                                                                                         race2=as.numeric(ORIGIN==1&RACE==2), race3=as.numeric(ORIGIN==2), 
                                                                                         race4=1-race1-race2-race3) %>% 
  dplyr::select(c('POPESTIMATE2020','race1','race2','race3','race3','race4')) %>% mutate(weights = POPESTIMATE2020 /sum(POPESTIMATE2020))
race1_population <- (race_data1$race1*race_data1$weights) %>% sum()
race2_population <- (race_data1$race2*race_data1$weights) %>% sum()
race3_population <- (race_data1$race3*race_data1$weights) %>% sum()
race4_population <- (race_data1$race4*race_data1$weights) %>% sum()


####################################################################################
##                          Load Experimental Data                                ##
####################################################################################

library(haven)
data <- read_dta("data/social trust data/PerPolarization_SocialTrust_Experiment.dta")
col_use <- c('condition', 'perpolindx', 'gstrustindx', 'ameritrustindx', "qage" , "female", "qincome", "qeducation", "marital",
             "race", "income","nocollege")
data1 <- data[,colnames(data)%in%col_use]
data1 <- data1 %>% mutate(race1=as.numeric(race==1), race2 = as.numeric(race==2), race3= as.numeric(race==3), race4=as.numeric(race==4)) 
nrow(data1); data1 %>% summary()
data1 %>% group_by(condition) %>% summarise_all(~{mean(.,na.rm = TRUE)}) %>% .[,-1]
data1 %>% summarise_all(~{mean(.,na.rm = TRUE)}) %>% .[,-1]
data1 <- data1 %>% mutate_all(~{replace_na(., mean(.,na.rm = TRUE))}) 


####################################################################################
## Set the global parameter treatment=2 (treatment=1) to recover the results      ##
## of less-polarization (more-polarization) vs control                            ##
####################################################################################

# Generalized trust index condition 1=less, 2=more, 3=control-----------
treatment <- 1

####################################################################################
##                          Data  Wrangling                                       ##
####################################################################################

##  relevant variable to use in the analysis
data2 <- data1 %>% dplyr::select(c('condition','gstrustindx','qage','qincome','female','qeducation','marital',
                                   'nocollege','race1','race2','race3')) %>% filter(condition%in%c(treatment,3)) %>% 
                  mutate(condition = as.numeric(condition==treatment))


# pretreatment-covariate matrix including the main terms and quadratic terms of relevant variables --------------------
p1 <- data2$condition %>% mean(); p0 <- 1-p1
data2 <-data2  %>% mutate(race1_v = (condition-p0)*(race1-race1_population),race2_v= (condition-p0)*(race2-race2_population),
                          race3_v = (condition-p0)*(race3-race3_population),age_v = (condition-p0)*(qage-age_population), 
                          female_v = (condition-p0)*(female-female_population)) %>% mutate(qage=as.double(qage%/%10)) 


data3 <-  data2 %>% mutate(across(everything(), ~ {attributes(.x) <- NULL; .x})) %>%
  mutate(
    TRT = structure(condition, "label" = "treatment"),
    gstrustindx = structure(gstrustindx, "label" = "outcome"),
    .keep = "unused",
    .before = everything()
  ) %>% 
  mutate(
    across(ends_with('_v'),~{structure(.x, "label" = "v")})
  ) %>% 
  mutate(
    across(where(
      ~ is.null(attributes(.x))
    ), ~ {structure(.x, "label" = "x")}
    )
  )%>% 
  mutate(
    across(where(
      ~ ifelse(attr(.x, "label") == "x",
               setequal(unique(.x), c(0, 1)),
               F)
    ), ~ structure(as.integer(.x), "label" = "x")
    )
  )


inter_name <- data3 %>% 
  dplyr::select(where(~ attr(.x, "label") == "x")) %>% 
  names %>% 
  combn(2)

dt_inter <- map2_dfc(inter_name[1,], inter_name[2,], ~ {
  tibble('{.x}:{.y}':=structure(data3[[.x]]*data3[[.y]], label='interaction'))
})
# add quadratic terms 
data4 <- data3 %>% 
mutate(across(where(
  ~ attr(.x, "label") == "predictor" & !is.integer(.x)
), ~ structure(.x^2, label = "quadratic"),
.names = "{.col}^2"
)) %>% 
# add interaction 
bind_cols(dt_inter) 

####################################################################################
##                          Regression Analyssis                                  ##
####################################################################################
# ToM regression ----------------------------
 
w <- ifelse(data2$condition==2,p1^(-2),p0^(-2))

full_model <- lm(data=data4,formula = gstrustindx ~., weights = w)
lm_tom <- stepAIC(full_model, direction = "both", k=log(nrow(data4)), trace = FALSE)
lm_tom <- update(lm_tom, .~.+TRT)

summary(lm_tom)

# Confidence interval

(htau_tom <- lm_tom %>% coef() %>% .['TRT'])
(V_hw_tom <- vcovHC(lm_tom,type = 'HC3')['TRT','TRT'])
htau_tom+sqrt(V_hw_tom)*qnorm(0.025)
htau_tom+sqrt(V_hw_tom)*qnorm(0.975)

if(treatment ==2){
  hbeta <- lm_tom %>% coef() %>% .[c('`qage:qeducation`', '`qage:nocollege`', '`qincome:race1`','`qeducation:nocollege`')]

  Z <- data2$condition
  
  htaux <-data4[c('qage:qeducation', 'qage:nocollege', 'qincome:race1','qeducation:nocollege')] %>% 
    summarise_all(~{sum(.*Z)/sum(Z)-sum(.*(1-Z))/sum(1-Z)}) %>% as.matrix()
  # Numerical results in Table S.1
  print(hbeta)
  print(htaux)
  print(hbeta*htaux)
  
  # Unadjust method ----------------------------------
  lm_naive <- lm(data=data4,formula = gstrustindx~TRT, weights = w)
  
  
  (htau_naive <- lm_naive %>% coef() %>% .['TRT'])
  print(htau_naive)
  (V_hw_naive <- vcovHC(lm_naive,type = 'HC3')['TRT','TRT'])
  print(htau_naive+sqrt(V_hw_naive)*qnorm(0.025))
  print(htau_naive+sqrt(V_hw_naive)*qnorm(0.975))
}

if(treatment==1){
  
 hbeta <- lm_tom %>% coef() %>% .[c('qage','`qage:qeducation`','`qage:marital`','`qage:race3`','`qincome:race1`','`qincome:race3`' , '`female:race2`' , '`qeducation:marital`', '`nocollege:race1`')]
  
  Z <- data2$condition

  htaux <-data4[c('qage','qage:qeducation','qage:marital','qage:race3','qincome:race1','qincome:race3' , 'female:race2' , 'qeducation:marital', 'nocollege:race1')] %>%
    summarise_all(~{sum(.*Z)/sum(Z)-sum(.*(1-Z))/sum(1-Z)}) %>% as.matrix()
  # Numerical results in Table S.1
  print(hbeta)
  print(htaux)
  print(hbeta*htaux)
  
  # Unadjust method ----------------------------------
  lm_naive <- lm(data=data4,formula = gstrustindx~TRT, weights = w)
  
  
  (htau_naive <- lm_naive %>% coef() %>% .['TRT'])
  print(htau_naive)
  (V_hw_naive <- vcovHC(lm_naive,type = 'HC3')['TRT','TRT'])
  print(htau_naive+sqrt(V_hw_naive)*qnorm(0.025))
  print(htau_naive+sqrt(V_hw_naive)*qnorm(0.975))
  
}


save(res_mat, file= 'res.Rdata')



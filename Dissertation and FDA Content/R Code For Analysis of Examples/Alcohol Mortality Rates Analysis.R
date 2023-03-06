library(fda)
library(tidyverse)

setwd("~/Dissertation/Dissertation and FDA Content/Datasets for Examples/Alcohol Example")

#### Load Data ####

mortality_rates_wide <- read_csv("Alcohol Mortality Rates Per Country.csv", skip = 6)
mortality_rates_long <- pivot_longer(mortality_rates_wide,values_to="MortalityRate",cols=2:5,names_to = "Country")


ggplot(mortality_rates_long,aes(x=Year,y=MortalityRate,col=Country))+
  geom_point()+
  ylab("Age-standardised death rates per 100,000 people")


Years <- unique(mortality_rates_long$Year)


#### Function Fitting ####

basis_mortality_rates <- create.bspline.basis(range(Years),
                                   breaks=Years,
                                   norder=4)


observation_matrix <- data.matrix(mortality_rates_wide[,c(-1)])


GCV_func <- function(log_lambda,basis,observations,time_points,penalty){
  
  lambda <- 10^log_lambda
  
  fd_par_obj <- fdPar(basis,penalty,lambda)
  
  smoothbasisobj <- smooth.basis(time_points,observations,fd_par_obj)
  
  return(sum(smoothbasisobj$gcv))
}

optimised_function <- optimise(GCV_func,lower=0,upper=10,basis=basis_mortality_rates,observations=observation_matrix,time_points=Years,penalty=2)

minimum_log_lambda <- optimised_function$minimum

minimum_lambda <- 10^minimum_log_lambda

fd_par_obj <- fdPar(basis_mortality_rates,2,minimum_lambda)

sample_of_functions <- smooth.basis(Years,observation_matrix,fd_par_obj)$fd

sample_of_functions$fdnames$time <- "Years"

sample_of_functions$fdnames$values <- "Mortality Rate"


year_mesh <- seq(2001,2020,0.01)



eval_df <- as.data.frame(eval.fd(year_mesh,sample_of_functions)) %>%
  mutate(Year=year_mesh)

eval_df_long <- pivot_longer(eval_df,names_to = "Country",values_to = "MortalityRate",cols=1:4)



# gives visualisation of sample of functions

ggplot(mortality_rates_long,aes(x=Year,y=MortalityRate,col=Country))+
  geom_point()+
  ylab("Age-standardised death rates per 100,000 people")

ggplot(eval_df_long,aes(x=Year,y=MortalityRate,col=Country))+
  geom_line()+
  geom_point(data=mortality_rates_long,aes(x=Year,y=MortalityRate,col=Country),inherit.aes = FALSE,alpha=0.2)+
  ylab("Age-standardised death rates per 100,000 people")


#### Checking Residuals #### 

shap_test <- function(vector){
  
  shapiro.test(vector)$p.value
  
}

residuals_matrix <- eval.fd(Years,sample_of_functions)-observation_matrix

apply(t(residuals_matrix),1,shap_test)


#### Summary Functions ####

#gives visualisation of the mean of this sample
mean_func <- mean.fd(sample_of_functions)

mean_plot_df <- data.frame(Year=year_mesh,mean=eval.fd(mean_func,year_mesh))

ggplot(eval_df_long,aes(x=Year,y=MortalityRate,col=Country))+
  geom_line()+
  geom_point(data=mortality_rates_long,aes(x=Year,y=MortalityRate,col=Country),inherit.aes = FALSE,alpha=0.2)+
  geom_line(data=mean_plot_df,aes(x=Year,y=mean),inherit.aes = FALSE) + 
  ylab("Age-standardised death rates per 100,000 people")



#gives standard deviation of this sample
sd_func <- sd.fd(sample_of_functions)

sd_plot_df <- data.frame(Year=year_mesh,stand_dev=eval.fd(sd_func,year_mesh))

ggplot(sd_plot_df,aes(x=Year,y=stand_dev)) +
  geom_line() + 
  ylab("Standard Deviation in Curves")


#### fPCA of Data ####

principle_components_of_sample <- pca.fd(sample_of_functions)

func_eval <- eval.fd(year_mesh,principle_components_of_sample$harmonics)


PC_DF <- data.frame(
  Year = year_mesh,
  PC_1 = func_eval[,1],
  PC_2 = func_eval[,2]
)

ggplot(PC_DF,aes(x=Year,y=PC_1))+
  geom_line() +
  ylab("Principal Component Function 1")+
  ggtitle(paste0("PC 1 - Explains ",100*principle_components_of_sample$varprop[1],"% Of Variation"))

ggplot(PC_DF,aes(x=Year,y=PC_2))+
  geom_line() +
  ylab("Principal Component Function 2")+
  ggtitle(paste0("PC 2 - Explains ",100*principle_components_of_sample$varprop[2],"% Of Variation"))



scores_df <- data.frame(
  PC_1 = principle_components_of_sample$scores[,1],
  PC_2 = principle_components_of_sample$scores[,2],
  Country=colnames(observation_matrix)
)

ggplot(scores_df,aes(x=PC_1,y=PC_2,label=Country))+
  geom_point()+
  geom_text(hjust = 0, nudge_x = 1)+
  xlim(c(-25,45)) +
  xlab("Score of Principle Component 1") + 
  ylab("Score of Principle Component 2")


#### PC of Derivatives ####

derivs <- deriv.fd(sample_of_functions,Lfdobj = 1) 

eval_df_derivs <- as.data.frame(eval.fd(year_mesh,derivs)) %>%
  mutate(Year=year_mesh)

eval_df_derivs_long <- pivot_longer(eval_df_derivs,names_to = "Country",values_to = "MortalityRate",cols=1:4)


# gives visualisation of sample of functions

ggplot(eval_df_derivs_long,aes(x=Year,y=MortalityRate,col=Country))+
  geom_line()+
  ylab("Change in Age-standardised death rates per 100,000 people")


principle_components_of_derivs <- pca.fd(derivs)

func_eval_derivs <- eval.fd(year_mesh,principle_components_of_derivs$harmonics)


PC_DF_derivs <- data.frame(
  Year = year_mesh,
  PC_1 = func_eval_derivs[,1],
  PC_2 = func_eval_derivs[,2]
)

ggplot(PC_DF_derivs,aes(x=Year,y=PC_1))+
  geom_line() +
  ylab("Principal Component Function 1")

ggplot(PC_DF_derivs,aes(x=Year,y=PC_2))+
  geom_line()


scores_df_derivs <- data.frame(
  PC_1 = principle_components_of_sample$scores[,1],
  PC_2 = principle_components_of_sample$scores[,2],
  Country=colnames(observation_matrix)
)

ggplot(scores_df_derivs,aes(x=PC_1,y=PC_2,label=Country))+
  geom_point()+
  geom_text(hjust = 0, nudge_x = 1)+
  xlim(c(-25,45)) +
  xlab("Score of Principle Component 1") + 
  ylab("Score of Principle Component 2")




#### Individually Fitting Functions #### 


fit_individ_func <- function(time_points,obs_matrix,basis){
  
  rep_dims <-  colnames(obs_matrix)
  
  for (i in 1:ncol(obs_matrix)){
    
    dim_i <- rep_dims[i]
    
    optimised_function <- optimise(GCV_func,lower=0,upper=10,basis=basis,observations=obs_matrix[,i],time_points=time_points,penalty=2)
    
    minimum_log_lambda <- optimised_function$minimum
    
    minimum_lambda <- 10^minimum_log_lambda
    
    fd_par_obj <- fdPar(basis,2,minimum_lambda)
    
    fd_obj <- smooth.basis(time_points,obs_matrix[,i],fd_par_obj)$fd
    
    coefs_i <- fd_obj$coefs
    
    colnames(coefs_i) <- dim_i
    
    if (i==1){
      coefs <- coefs_i
    } else{
      coefs <- cbind(coefs, coefs_i)
    }

  }
  
  merged_fd <- fd(coef = coefs, basisobj =basis)
  merged_fd$fdnames$reps <- rep_dims
  
  return(merged_fd)
  
}


ind_fit_funcs <- fit_individ_func(Years,observation_matrix,basis_mortality_rates)


eval_df_indv <- as.data.frame(eval.fd(year_mesh,ind_fit_funcs)) %>%
  mutate(Year=year_mesh)

eval_df_indv_long <- pivot_longer(eval_df_indv,names_to = "Country",values_to = "MortalityRate",cols=1:4)



# gives visualisation of sample of functions

ggplot(mortality_rates_long,aes(x=Year,y=MortalityRate,col=Country))+
  geom_point()+
  ylab("Age-standardised death rates per 100,000 people")

ggplot(eval_df_indv_long,aes(x=Year,y=MortalityRate,col=Country))+
  geom_line()+
  geom_point(data=mortality_rates_long,aes(x=Year,y=MortalityRate,col=Country),inherit.aes = FALSE,alpha=0.2)+
  ylab("Age-standardised death rates per 100,000 people")


residuals_matrix <- eval.fd(Years,ind_fit_funcs)-observation_matrix

apply(t(residuals_matrix),1,shap_test)



#### Individually Fit Functions PCA ####


principle_components_of_sample_ind <- pca.fd(ind_fit_funcs)

func_eval_ind <- eval.fd(year_mesh,principle_components_of_sample_ind$harmonics)


PC_DF_ind <- data.frame(
  Year = year_mesh,
  PC_1 = func_eval_ind[,1],
  PC_2 = func_eval_ind[,2]
)

ggplot(PC_DF_ind,aes(x=Year,y=PC_1))+
  geom_line() +
  ylab("Principal Component Function 1")+
  ggtitle(paste0("PC 1 - Explains ",100*principle_components_of_sample_ind$varprop[1],"% Of Variation"))

ggplot(PC_DF_ind,aes(x=Year,y=PC_2))+
  geom_line() +
  ylab("Principal Component Function 2")+
  ggtitle(paste0("PC 2 - Explains ",100*principle_components_of_sample_ind$varprop[2],"% Of Variation"))


scores_df_ind <- data.frame(
  PC_1 = principle_components_of_sample_ind$scores[,1],
  PC_2 = principle_components_of_sample_ind$scores[,2],
  Country=colnames(observation_matrix)
)

ggplot(scores_df_ind,aes(x=PC_1,y=PC_2,label=Country))+
  geom_point()+
  geom_text(hjust = 0, nudge_x = 1)+
  xlim(c(-25,45)) +
  xlab("Score of Principle Component 1") + 
  ylab("Score of Principle Component 2")




#### Load in Unemployment Data ####

Scot_UnE <- read.csv("Scotland Unemployment.csv",skip=7) %>%
              rename("Year"=Important.notes,"Unemployment_Rate"=X) %>%
              filter(Year %in% Years) %>% 
              mutate(Year=as.numeric(Year)) %>%
              mutate(Country="Scotland")

Eng_UnE <- read.csv("England Unemployment.csv",skip=7) %>%
  rename("Year"=Important.notes,"Unemployment_Rate"=X) %>%
  filter(Year %in% Years) %>%
  mutate(Year=as.numeric(Year)) %>%
  mutate(Country="England")

Wales_UnE <- read.csv("Wales Unemployment.csv",skip=7) %>%
  rename("Year"=Important.notes,"Unemployment_Rate"=X) %>%
  filter(Year %in% Years) %>%
  mutate(Year=as.numeric(Year)) %>%
  mutate(Country="Wales")

NI_UnE <- read.csv("Northern Ireland Unemployment.csv",skip=7) %>%
  rename("Year"=Important.notes,"Unemployment_Rate"=X) %>%
  filter(Year %in% Years) %>%
  mutate(Year=as.numeric(Year)) %>%
  mutate(Country="Northern Ireland")

All_Countries_UnE <- bind_rows(Scot_UnE,Eng_UnE,Wales_UnE,NI_UnE)

All_Countries_UnE_Wide <- pivot_wider(All_Countries_UnE,names_from = Country,values_from = Unemployment_Rate)



#### Function Fitting ####

Years_UnE <- unique(All_Countries_UnE$Year)

basis_UnE <- create.bspline.basis(range(Years_UnE),
                                              breaks=Years_UnE,
                                              norder=4)


observation_matrix_UnE <- data.matrix(All_Countries_UnE_Wide[,c(-1)])

optimised_function_UnE <- optimise(GCV_func,lower=0,upper=10,basis=basis_UnE,observations=observation_matrix_UnE,time_points=Years_UnE,penalty=2)

minimum_log_lambda_UnE <- optimised_function_UnE$minimum

minimum_lambda_UnE <- 10^minimum_log_lambda_UnE

fd_par_obj_UnE <- fdPar(basis_UnE,2,minimum_lambda_UnE)

sample_of_functions_UnE <- fit_individ_func(Years_UnE,observation_matrix_UnE,basis_UnE)

sample_of_functions_UnE$fdnames$time <- "Years"

sample_of_functions_UnE$fdnames$values <- "Unemployment Rate"

eval_df_UnE <- as.data.frame(eval.fd(year_mesh,sample_of_functions_UnE)) %>%
  mutate(Year=year_mesh)

eval_df_long_UnE <- pivot_longer(eval_df_UnE,names_to = "Country",values_to = "Unemployment_Rate",cols=1:4)


# gives visualisation of sample of functions

ggplot(All_Countries_UnE,aes(x=Year,y=Unemployment_Rate,col=Country))+
  geom_point()+
  ylab("Unemployment Rate %")

ggplot(eval_df_long_UnE,aes(x=Year,y=Unemployment_Rate,col=Country))+
  geom_line()+
  geom_point(data=All_Countries_UnE,aes(x=Year,y=Unemployment_Rate,col=Country),inherit.aes = FALSE,alpha=0.2)+
  ylab("Unemployment Rate %")

#### Checking Residuals ####

shap_test <- function(vector){
  
  shapiro.test(vector)$p.value
  
}

residuals_matrix <- eval.fd(Years_UnE,sample_of_functions_UnE)-observation_matrix_UnE

apply(t(residuals_matrix),1,shap_test)


#### PC Analysis of Unemployment ####


principle_components_of_sample_UnE <- pca.fd(sample_of_functions_UnE)

func_eval_UnE <- eval.fd(year_mesh,principle_components_of_sample_UnE$harmonics)


PC_DF_UnE <- data.frame(
  Year = year_mesh,
  PC_1 = func_eval_UnE[,1],
  PC_2 = func_eval_UnE[,2]
)

ggplot(PC_DF_UnE,aes(x=Year,y=PC_1))+
  geom_line() +
  ylab("Principal Component Function 1")+
  ggtitle(paste0("PC 1 - Explains ",100*principle_components_of_sample_UnE$varprop[1],"% Of Variation"))

ggplot(PC_DF_UnE,aes(x=Year,y=PC_2))+
  geom_line() +
  ylab("Principal Component Function 2")+
  ggtitle(paste0("PC 2 - Explains ",100*principle_components_of_sample_UnE$varprop[2],"% Of Variation"))


scores_df_UnE <- data.frame(
  PC_1 = principle_components_of_sample_UnE$scores[,1],
  PC_2 = principle_components_of_sample_UnE$scores[,2],
  Country=colnames(observation_matrix_UnE)
)

ggplot(scores_df_UnE,aes(x=PC_1,y=PC_2,label=Country))+
  geom_point()+
  geom_text(hjust = 0, nudge_x = 0.08)+
  xlim(c(-3,2))+
  xlab("Score of Principle Component 1") + 
  ylab("Score of Principle Component 2")


#### Functional Regression ####

Countries <- unique(mortality_rates_long$Country)



sample_of_functions
sample_of_functions_UnE



fda::fRegress()

list_covariate <- list(rep(1,length(Countries)),sample_of_functions_UnE)

beta_fdPar <- fdPar(basis_mortality_rates,2)

beta_list <- list(beta_fdPar,beta_fdPar)

functional_regression <- fRegress(sample_of_functions,list_covariate,beta_list)

plot.fd(functional_regression$yhatfdobj)

estimate_of_response <- functional_regression$yhatfdobj

intercept <- functional_regression$betaestlist[[1]]

beta <- functional_regression$betaestlist[[2]]$fd

plot(beta)


#### GCV Plot (For Interest)####

loglambda <- seq(-6,6,0.25)

list_of_gcvs <- c()

list_of_dfs <- c()

for (log_lam in loglambda){
  
  lambda <- 10^log_lam
  
  fd_par_obj <- fdPar(basis_mortality_rates,2,lambda)
  
  smoothbasisobj <- smooth.basis(Years,as.matrix(mortality_rates_wide[,c(-1)]),fd_par_obj)
  
  list_of_gcvs <- c(list_of_gcvs,sum(smoothbasisobj$gcv))
  
  list_of_dfs <- c(list_of_dfs,fd_par_obj$df)
  
}

log_lambda_df <- data.frame(
  loglambda=loglambda,
  lambda=10^loglambda,
  GCV=list_of_gcvs
)


ggplot(log_lambda_df,aes(x=loglambda,y=GCV)) +
  geom_point() +
  geom_line()+
  geom_vline(xintercept = minimum_log_lambda)



#### SALSA (No Longer In Use) #### 

colnames(observation_matrix)

observation_matrix[1,]


SALSA_Fit_fda <- function(x_axis,obs_matrix,start_knots,min_knots,max_knots,degree,maxIter,gaps){
  
  rep_dims <-  colnames(obs_matrix)
  
  for (i in 1:ncol(obs_matrix)){
    
    dim <- rep_dims[i]
    
    data_SALSA <- data.frame(
      
      response=obs_matrix[,i],
      x_axis=x_axis
      
    )
  
    
    initialModel <- glm(response ~ -1, data=data_SALSA,family="Gaussian")
    
    varList <- c("x_axis")
    
    SALSA1DList <- list(fitnessMeasure="BIC", 
                        minKnots_1d=min_knots, maxKnots_1d=max_knots, 
                        startKnots_1d=start_knots, degree=degree, 
                        maxIterations=maxIter, gaps=gaps)
    
    # Run SALSA
    SALSA <- MRSea::runSALSA1D(initialModel=initialModel, 
                               salsa1dlist=SALSA1DList, 
                               varlist=varList, 
                               factorlist=NULL,
                               datain=data_SALSA,
                               splineParams=NULL,
                               suppress.printout=TRUE)
    
    SALSA_Fit <- SALSA$bestModel
    
    internal_knots <- SALSA_Fit$splineParams[[2]]$knots
    
    knots <- append(x_axis[1],internal_knots)
    
    knots <- append(knots,x_axis[length(x_axis)])
    
    #smoothing_parameter <- 
    
    basis <- create.bspline.basis(range(x_axis),
                                  breaks=knots,
                                  norder=degree+1)
    
    #fd_par_obj <- fdPar(basis,2,smoothing_parameter)
    
    #smoothbasisobj <- smooth.basis(x_axis,obs_matrix,fd_par_obj)
    
    #### COMPLETE THIS ####
    
    #merged_coefs <- cbind(fd1$coefs, fd2$coefs)
    #merged_fd <- fd(coef = merged_coefs, basisobj = minutebasis)
    #is.fd(merged_fd)
    
  }
  
  
  
  return(initialModel)
  
}


edit(getAnywhere('runSALSA1DFit'), file='source_rfcv.r')

edit(MRSea::runSALSA1D)

observation_matrix


SALSA_Fit_fda(Years,
              observation_matrix,
              1,
              1,
              30,
              3,
              30,
              0)


data <- mortality_rates_long %>%
  rename("response"=MortalityRate) %>%
  filter(Country=="England")

initialModel <- glm(response ~ -1, data=data)

varList <- c("Year")

SALSA1DList <- list(fitnessMeasure="BIC", 
                    minKnots_1d=1, maxKnots_1d=length(Years), 
                    startKnots_1d=1, degree=3, 
                    maxIterations=50, gaps=0)

# Run SALSA
SALSA <- MRSea::runSALSA1D(initialModel=initialModel, 
                           salsa1dlist=SALSA1DList, 
                           varlist=varList, 
                           factorlist=NULL,
                           datain=data,
                           splineParams=NULL,
                           suppress.printout=TRUE)


SALSA_Fit <- SALSA$bestModel

summary(SALSA_Fit)

SALSA_Fit$coefficients


fit_df <- data.frame(x=year_mesh,Fitted=predict(object=SALSA_Fit,newdata=data.frame(Year=year_mesh)),Country=rep("Scotland",length(year_mesh)))

SALSA_knot_plot <- ggplot(fit_df,aes(x=year_mesh,y=Fitted,col=Country)) +
  geom_line() +
  geom_point(data=mortality_rates_long,aes(x=Year,y=MortalityRate,col=Country),inherit.aes = FALSE,alpha=0.2)+
  ylab("Age-standardised death rates per 100,000 people") +
  facet_wrap(~Country)


SALSA_knot_plot


SALSA_Fit$splineParams


SALSA$bestModel$splineParams




  
  
  

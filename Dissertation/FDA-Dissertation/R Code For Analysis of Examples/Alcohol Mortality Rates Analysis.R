library(fda)
library(tidyverse)

setwd("~/Dissertation/FDA-Dissertation/Datasets for Examples")


mortality_rates_wide <- read_csv("Alcohol Mortality Rates Per Country.csv", skip = 6)
mortality_rates_long <- pivot_longer(mortality_rates_wide,values_to="MortalityRate",cols=2:5,names_to = "Country")


ggplot(mortality_rates_long,aes(x=Year,y=MortalityRate,col=Country))+
  geom_line()+
  ylab("Age-standardised death rates per 100,000 people")

Years <- unique(mortality_rates_long$Year)


basis_mortality_rates <- create.bspline.basis(range(Years),
                                   breaks=Years,
                                   norder=4)


observation_matrix <- as.matrix(mortality_rates_wide[,c(-1)])





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

sample_of_functions$fdnames$time <- Years

sample_of_functions$fdnames$values <- "Mortality Rate"



par(mfrow=c(1,1))

# gives visualisation of sample of functions
plot(sample_of_functions)

#gives visualisation of the mean of this sample
plot(mean.fd(sample_of_functions),ylim=c(8,30),col="black",lwd=2)
lines(sample_of_functions)

#gives standard deviation of this sample
plot(sd.fd(sample_of_functions))


principle_components_of_sample <- pca.fd(sample_of_functions)


principle_components_of_sample$varprop



















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



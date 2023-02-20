library(openair)
library(tidyverse)
library(fda)
library(worldmet)
library(openairmaps)
library(lubridate)
library(zoo)

#### Data Importing ####

# Shows All Information About All Monitoring Stations in The AURN network
station_data <- importMeta(source="aurn", all=TRUE)

# Shows information on stations which monitor PM2.5 still
stations_pm2.5 <- station_data %>%
  filter(Parameter_name == "PM2.5 particulate matter (Hourly measured)") %>%
  dplyr::select(site,code,start_date,end_date) %>%
  mutate(start_Year=year(start_date)) %>%
  filter(start_Year < 2008 & end_date == "ongoing") %>%
  dplyr::select(code) %>%
  t() %>%
  as.vector()
 
#gets data for particular site 
site_data <- importMeta(source = "aurn") %>%
  filter(site == "Dundee Mains Loan")

#site codes for sites in glasgow
glasgow_sites <- c("GLA4","GLKP","GHSR","GGWR")

# imports data for glasgow station for 2022
raw_data <- importAURN(
  site=stations_pm2.5,
  year=2008:2022,
  data_type="daily"
)

# stores the names of sites 
sites <- unique(raw_data$site)

# filters data for only columns of interest
data_longer <- raw_data %>%
  dplyr::select(date,site,pm2.5) %>%

# adds column which gives day data was taken on
  mutate(Day=day(date),Month=month(date),Year=year(date))


# creates new dataset of daily averages
daily_avg <- data_longer %>%
  mutate(Day_Month=as.Date(paste(Day,Month,Year,sep="-"),format="%d-%m-%Y")) %>%
  dplyr::select(-c(Year,Month,Day,date))

#pivots data so there is a column for each station
data_wider <- pivot_wider(data = daily_avg,names_from = "site",values_from = "pm2.5") 

#removes rows from data which contain NAs
data_wider_no_NAs <- data_wider %>%
  na.omit()

data_no_NAs_longer <- pivot_longer(data_wider_no_NAs,names_to = "site", values_to = "pm2.5",cols = 2:6)



# creates matrix of observations from wider dataframe with NAs removed
obs_matrix_no_NAs <- as.matrix(data_wider_no_NAs %>% dplyr::select(-c(Day_Month)))

obs_matrix_no_NAs_log <- log(obs_matrix_no_NAs)


# extracts time points being used and scales them so initial observation is taken as day 0
time_points <- (as.numeric(data_wider$Day_Month)-min(as.numeric(data_wider$Day_Month)))

# repeates process of data with NAs removed
time_points_no_NAs <- (as.numeric(data_wider_no_NAs$Day_Month)-min(as.numeric(data_wider_no_NAs$Day_Month)))


#### Plot Summaries of Data ####

# Plots PM2.5 emissions over time for each station 
ggplot(daily_avg,aes(x=Day_Month,y=pm2.5,col=site))+
  geom_line()+
  facet_wrap(~site)+
  ylim(c(0,50))

# Plots PM2.5 emissions with NAs removed
ggplot(data_no_NAs_longer,aes(x=Day_Month,y=pm2.5,col=site))+
  geom_line()+
  facet_wrap(~site)+
  ylim(c(0,50))





#### GCV Function Definition #### 

# creates function which outputs the GCV of a functional data object 
# fitted using a specific lambda 

GCV_func <- function(log_lambda,basis,observations,time_points,penalty){
  
  lambda <- 10^log_lambda
  
  fd_par_obj <- fdPar(basis,penalty,lambda)
  
  smoothbasisobj <- smooth.basis(time_points,observations,fd_par_obj)
  
  return(sum(smoothbasisobj$gcv))
}

## Functions Fitted Omitting NAs 

#creates B-Spline basis of order 7 with 60 basis functions

basis_no_NAs <- create.bspline.basis(range(time_points_no_NAs),
                                     nbasis = 60,
                                     norder=7)

# Finds optimal lambda and extracts from output

optimised_function_no_NAs <- optimise(GCV_func,lower=0,upper=10,
                                      basis=basis_no_NAs,observations=obs_matrix_no_NAs,
                                      time_points=time_points_no_NAs,penalty=2)

minimum_log_lambda_no_NAs <- optimised_function_no_NAs$minimum

minimum_lambda_no_NAs <- 10^minimum_log_lambda_no_NAs


# FD Par object created which specifies the type of smoothing as well as the lambda to be used

fd_par_obj_no_NAs <- fdPar(basis_no_NAs,2,minimum_lambda_no_NAs)

# Creates sample of functions

sample_of_functions_no_NAs <- smooth.basis(time_points_no_NAs,obs_matrix_no_NAs,fd_par_obj_no_NAs)


#### Plots of Sample of Functions NAs Removed

estimated_data <- as.data.frame(eval.fd(time_points,sample_of_functions_no_NAs$fd)) %>%
  mutate(Day_Month=data_wider$Day_Month) %>%
  dplyr::select(c(6,1,2,3,4,5))

estimated_data_longer <- pivot_longer(estimated_data,names_to = "site",values_to = "pm2.5",cols=c(2,3,4,5,6)) 

par(mfrow=c(1,1))

plot(sample_of_functions_no_NAs)

ggplot(estimated_data_longer,aes(x=Day_Month,y=pm2.5,col=site))+
  geom_line() +
  geom_line(daily_avg,mapping=aes(x=Day_Month,y=pm2.5,col=site),inherit.aes = FALSE,alpha=0.2)

#### Extract Fitted Values and Back-Fill ####

data_NAs_filled <- daily_avg

for (i in 1:length(data_longer$date)){
  
  check <- daily_avg$pm2.5[i] 
  
  if (is.na(check) == TRUE){
    
    date_to_fill <- daily_avg$Day_Month[i]
    
    site_to_fill <- daily_avg$site[i]
    
    data_NAs_filled$pm2.5[i] <- estimated_data_longer %>%
      filter(Day_Month==date_to_fill & site==site_to_fill) %>%
      dplyr::select(pm2.5) %>%
      as.vector() %>%
      as.numeric()
    
  }
  
}

#pivots data so there is a column for each station
data_wider_NAs_filled <- pivot_wider(data = data_NAs_filled,names_from = "site",values_from = "pm2.5")


# creates matrix of observations from wider dataframe 
obs_matrix_NAs_filled <- as.matrix(data_wider_NAs_filled %>% dplyr::select(-c(Day_Month)))



#### Functions Fitted With NAs replaced

# creates basis

basis <- create.bspline.basis(range(time_points),
                              nbasis = 60,
                              norder=7)

# Finds optimal lambda and extracts from output

optimised_function <- optimise(GCV_func,lower=0,upper=10,
                               basis=basis, observations=obs_matrix_NAs_filled,
                               time_points=time_points,penalty=2)

minimum_log_lambda <- optimised_function$minimum

minimum_lambda <- 10^minimum_log_lambda

# FD Par object created which specifies the type of smoothing as well as the lambda to be used
fd_par_obj <- fdPar(basis,2,minimum_lambda)

sample_of_functions <- smooth.basis(time_points,obs_matrix_NAs_filled,fd_par_obj)



func_plot <- as.data.frame(eval.fd(time_points,sample_of_functions$fd)) %>%
  mutate(Day_Month=data_wider$Day_Month) %>%
  dplyr::select(c(6,1,2,3,4,5))

func_plot_longer <- pivot_longer(func_plot,names_to = "site",values_to = "pm2.5",cols=c(2,3,4,5,6)) 



plot(sample_of_functions)

ggplot(daily_avg,mapping=aes(x=Day_Month,y=pm2.5,col=site))+
  geom_line(alpha=0.4)


ggplot(func_plot_longer,aes(x=Day_Month,y=pm2.5,col=site))+
  geom_line() +
  geom_line(daily_avg,mapping=aes(x=Day_Month,y=pm2.5,col=site),inherit.aes = FALSE,alpha=0.2)




#### Summaries of Data ####


plot(mean.fd(sample_of_functions$fd))

plot(sd.fd(sample_of_functions$fd))






#### PCA of Data ####

PCA_of_func <- pca.fd(sample_of_functions$fd)

scores_df <- data.frame(
  PC_1 = PCA_of_func$scores[,1],
  PC_2 = PCA_of_func$scores[,2]
)

PCs <- PCA_of_func$harmonics

ggplot(scores_df,aes(x=PC_1,y=PC_2))+
  geom_point()

plot(PCA_of_func)



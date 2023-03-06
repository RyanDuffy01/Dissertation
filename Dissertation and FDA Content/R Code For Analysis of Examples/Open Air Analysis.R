library(openair)
library(tidyverse)
library(fda)
library(worldmet)
library(openairmaps)

#### Data Importing ####

# Shows All Information About All Monitoring Stations in The AURN network
station_data <- importMeta(source="aurn", all=TRUE)[1:40,]

#gets data for particular site 
site_data <- importMeta(source = "aurn") %>%
  filter(site == "Dundee Mains Loan")

#site codes for sites in glasgow
glasgow_sites <- c("GLA4","GLKP","GHSR","GGWR")

# imports data for glasgow station for 2022
raw_data <- importAURN(
  site=glasgow_sites,
  year=2018:2022
  
)

#### Data Cleaning / Extraction of Relevant Info ####


# stores the names of sites 
sites <- unique(raw_data$site)

# filters data for only columns of interest
data_longer <- raw_data %>%
  dplyr::select(date,site,nox)

#pivots data so there is a column for each station
data_wider <- pivot_wider(data = data_longer,names_from = "site",values_from = "nox")

#removes rows from data which contain NAs
data_wider_no_NAs <- data_wider %>%
  na.omit()

# creates matrix of observations from wider dataframe with NAs removed
obs_matrix_no_NAs <- as.matrix(data_wider_no_NAs %>% dplyr::select(-c(date)))

# extracts time points being used and scales them so initial observation is taken as 0
time_points <- (as.numeric(data_wider$date)-min(as.numeric(data_wider$date)))*3.80517183071e-7

# repeates process of data with NAs removed
time_points_no_NAs <- (as.numeric(data_wider_no_NAs$date)-min(as.numeric(data_wider_no_NAs$date)))*3.80517183071e-7



#### Plot Summaries of Data ####

# Plots NOx emissions over time for each station 
ggplot(data_longer,aes(x=date,y=nox,col=site))+
  geom_line()+
  facet_wrap(~site)

# Plots Wind Speed over time for eaxh station
ggplot(raw_data,aes(x=date,y=ws,col=site))+
  geom_line()+
  facet_wrap(~site)


#### GCV Function Definition #### 

# creates function which outputs the GCV of a functional data object 
# fitted using a specific lambda 

GCV_func <- function(log_lambda,basis,observations,time_points,penalty){
  
  lambda <- 10^log_lambda
  
  fd_par_obj <- fdPar(basis,penalty,lambda)
  
  smoothbasisobj <- smooth.basis(time_points,observations,fd_par_obj)
  
  return(sum(smoothbasisobj$gcv))
}


#### Functions Fitted Omitting NAs ####

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


#### Plots of Sample of Functions NAs Removed ####

estimated_data <- as.data.frame(eval.fd(time_points,sample_of_functions_no_NAs$fd)) %>%
  mutate(date=data_wider$date) %>%
  dplyr::select(c(5,1,2,3,4))

estimated_data_longer <- pivot_longer(estimated_data,names_to = "site",values_to = "nox",cols=c(2,3,4,5)) 

par(mfrow=c(1,1))

plot(sample_of_functions_no_NAs)

ggplot(estimated_data_longer,aes(x=date,y=nox,col=site))+
  geom_line() +
  geom_line(data_longer,mapping=aes(x=date,y=nox,col=site),inherit.aes = FALSE,alpha=0.2)


#### Extract Fitted Values and Back-Fill ####


data_NAs_filled <- data_longer
  
for (i in 1:length(data_longer$date)){
  
  check <- data_longer$nox[i] 
  
  if (is.na(check) == TRUE){
    
    date_to_fill <- data_longer$date[i]
    
    site_to_fill <- data_longer$site[i]
    
    data_NAs_filled$nox[i] <- estimated_data_longer %>%
      filter(date==date_to_fill & site==site_to_fill) %>%
      dplyr::select(nox) %>%
      as.vector() %>%
      as.numeric()
    
  }
  
}

#pivots data so there is a column for each station
data_wider_NAs_filled <- pivot_wider(data = data_NAs_filled,names_from = "site",values_from = "nox")


# creates matrix of observations from wider dataframe 
obs_matrix_NAs_filled <- as.matrix(data_wider_NAs_filled %>% dplyr::select(-c(date)))



#### Functions Fitted With NAs replaced ####

# creates basis

basis <- create.bspline.basis(range(time_points),
                              nbasis = 40,
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


#### Plots of Sample of Functions NAs Removed ####

estimated_data_NAs_filled <- as.data.frame(eval.fd(time_points,sample_of_functions_no_NAs$fd)) %>%
  mutate(date=data_wider_NAs_filled$date) %>%
  dplyr::select(c(5,1,2,3,4))

estimated_data_NAs_filled_longer <- pivot_longer(estimated_data_NAs_filled,names_to = "site",values_to = "nox",cols=c(2,3,4,5)) 

par(mfrow=c(1,1))

plot(sample_of_functions)

ggplot(estimated_data_NAs_filled_longer,aes(x=date,y=nox,col=site))+
  geom_line() +
  geom_line(data_longer,mapping=aes(x=date,y=nox,col=site),inherit.aes = FALSE,alpha=0.2)


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


#### Functional Regression Of Data ####












#### TEST ZONE ####

closest_site_code <- function(longitude,latitude){
  
  list_of_codes <- c()
  
  for (i in 1:length(longitude)){
    
    meta_data <- getMeta(lat = latitude[i],lon = longitude[i],n=1,returnMap = FALSE)
    
    list_of_codes <- c(list_of_codes,meta_data$code)
  }
  
  list_of_codes
  
}

closest_site_code(station_data$longitude,station_data$longitude)


site_met <- importMeta()
my1 <- filter(site_met, code == "MY1")


aq <- importAURN(site = "my1", pollutant = c("no2", "nox", "pm2.5", "pm10"), year = 2020)
met <- worldmet::importNOAA(code = "037720-99999", year = 2020)
all <- left_join(select(met, date, ws, wd), aq, by = "date")
all

aq
met

data$pm2.5


# imports all data from all stations for a year
data <- importAURN(
  data_type = "annual",
  year=2016
)

# Gets data about all MET weather stations
getMeta()


# Gets data about MET stations nearest to specified latitude and longitude
h <- getMeta(lat = station_data$latitude[1],lon = station_data$longitude[1],n=1)



getMeta()











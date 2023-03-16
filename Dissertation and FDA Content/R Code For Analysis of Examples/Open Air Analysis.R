library(openair)
library(tidyverse)
library(fda)
library(worldmet)
library(openairmaps)
library(PostcodesioR)

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
  dplyr::select(date,site,nox) %>%
  mutate(nox=replace(nox,nox<0,0))


#pivots data so there is a column for each station
data_wider <- pivot_wider(data = data_longer,names_from = "site",values_from = "nox")

#removes rows from data which contain NAs
data_wider_no_NAs <- data_wider %>%
  na.omit()

data_longer_no_NAs <- data_wider_no_NAs %>%
  pivot_longer(cols=c(2:length(colnames(data_wider_no_NAs))),values_to = 'nox',names_to = 'site')

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

# Plots NOx emissions over time for each station 
ggplot(data_longer_no_NAs,aes(x=date,y=nox,col=site))+
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

ggplot(estimated_data_longer,aes(x=date,y=nox,col=site))+
  geom_line() +
  geom_line(data_longer,mapping=aes(x=date,y=nox,col=site),inherit.aes = FALSE,alpha=0.05)


#### Extract Fitted Values and Back-Fill ####

data_NAs_filled <- data_wider

rows_with_NAs <- rowSums(is.na(data_NAs_filled)) > 0

for (i in 1:length(rows_with_NAs)){
  
  if (rows_with_NAs[i] == TRUE){
    
    row_to_be_filled <- data_NAs_filled[i,]
    
    date_to_fill <- row_to_be_filled$date
    
    cols_to_be_filled <- is.na(row_to_be_filled)
    
    sites_to_fill <- colnames(row_to_be_filled)[cols_to_be_filled]
    
    for (j in 1:length(sites_to_fill)){
      
      value <- estimated_data %>%
        filter(date==date_to_fill) %>%
        dplyr::select(sites_to_fill[j]) %>%
        as.numeric()
      
      row_to_be_filled[,colnames(row_to_be_filled) == sites_to_fill[j]] <- value 
      
      data_NAs_filled[i,] <- row_to_be_filled
      
    }
  }
}




#pivots data so there is a column for each station
data_longer_NAs_filled <- pivot_longer(data = data_NAs_filled,names_to = "site",values_to = "nox",cols=2:length(colnames(data_NAs_filled)))


# creates matrix of observations from wider dataframe 
obs_matrix_NAs_filled <- as.matrix(data_NAs_filled %>% dplyr::select(-c(date)))



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


#### Checking Residuals ####

shap_test <- function(vector){
  
  shapiro.test(vector)$p.value
  
}

residuals_matrix <- eval.fd(time_points,sample_of_functions$fd)-obs_matrix_NAs_filled

residuals_df <- data.frame(residuals_matrix) 

colnames(residuals_df) <- colnames(residuals_matrix)

residuals_long <- residuals_df %>%
                    pivot_longer(values_to = 'nox',names_to = 'site',cols=1:length(colnames(residuals_df)))
  
  
ggplot(residuals_long,aes(x=nox))+
  geom_histogram()+
  facet_wrap(~site)

sapply(MASS::fitdistr(residuals_long$nox, "normal")$estimate)


ggplot(residuals_long, aes(sample = nox)) +
  stat_qq(distribution = qnorm) +
  stat_qq_line(distribution = qnorm) +
  facet_wrap(~site)





#### Plots of Sample of Functions NAs Removed ####

estimated_data_NAs_filled <- as.data.frame(eval.fd(time_points,sample_of_functions_no_NAs$fd)) %>%
  mutate(date=data_NAs_filled$date) %>%
  dplyr::select(c(5,1,2,3,4))

estimated_data_NAs_filled_longer <- pivot_longer(estimated_data_NAs_filled,names_to = "site",values_to = "nox",cols=c(2,3,4,5)) 

par(mfrow=c(1,1))

plot(sample_of_functions)

ggplot(estimated_data_NAs_filled_longer,aes(x=date,y=nox,col=site))+
  geom_line() +
  geom_line(data_longer,mapping=aes(x=date,y=nox,col=site),inherit.aes = FALSE,alpha=0.2)


#### Summaries of Data ####

mean_func <- mean.fd(sample_of_functions$fd)

sd_func <- sd.fd(sample_of_functions$fd)



#### PCA of Data ####

PCA_of_func <- pca.fd(sample_of_functions$fd)

scores_df <- data.frame(
  PC_1 = PCA_of_func$scores[,1],
  PC_2 = PCA_of_func$scores[,2]
)

PCs <- PCA_of_func$harmonics

ggplot(scores_df,aes(x=PC_1,y=PC_2))+
  geom_point()


#### Preparing SIMD Data ####

setwd("~/Dissertation/Dissertation and FDA Content/Datasets for Examples/Open Air/Poverty Index Data")

pov_ind_data <- read.csv("postcode_2020_1_simd2020v2.csv") %>%
  mutate(pc7=gsub("  "," ",pc7))

head(pov_ind_data)

get_post_func <- function(longitude,latitude){
  
  postcodes <- c()
  
  for (i in 1:length(longitude)){
    
    postcodes <- c(postcodes,reverse_geocoding(longitude=longitude[i],latitude=latitude[i])[[1]]$postcode)
    
  }
  
  return(postcodes)
  
}

#gets data for particular site 
site_data <- importMeta(source = "aurn") %>%
  filter(site %in% sites) %>%
  mutate(Postcode = get_post_func(longitude,latitude)) 


pov_ind_data_for_stations <- pov_ind_data %>% 
                                filter(pc7 %in% site_data$Postcode) %>%
                                rename("Postcode"=pc7) %>%
                                dplyr::select(c(Postcode,simd2020v2_sc_decile))
  
site_data_full <- full_join(site_data,pov_ind_data_for_stations,by="Postcode")




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




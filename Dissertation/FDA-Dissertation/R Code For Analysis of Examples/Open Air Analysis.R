library(openair)
library(tidyverse)
library(fda)
library(worldmet)
library(openairmaps)

# Shows All Information About All Monitoring Stations in The AURN network
station_data <- importMeta(source="aurn", all=TRUE)[1:40,]

# imports data for station YK10 for the years 2015 to 2020
data <- importAURN(
  site=c("YK10"),
  year=2022
  
)

as.numeric(data$date)


ggplot(data,aes(x=date,y=pm10))+
  geom_line()


ggplot(data,aes(x=date,y=ws))+
  geom_line()




non_na_entries <- which(is.na(data$no))
observations <- data$no[-non_na_entries]
time_points <- as.numeric(data$date[-non_na_entries])


basis <- create.bspline.basis(range(time_points),
                              nbasis = 15,
                              norder=4)


GCV_pos_func <- function(log_lambda,basis,observations,time_points,penalty){
  
  lambda <- 10^log_lambda
  
  fd_par_obj <- fdPar(basis,penalty,lambda)
  
  smoothbasisobj <- smooth.pos(time_points,observations,fd_par_obj)
  
  return(sum(smoothbasisobj$gcv))
}

optimised_function <- optimise(GCV_pos_func,lower=0,upper=10,basis=basis,observations=observations,time_points=time_points,penalty=2)

minimum_log_lambda <- 2.45

minimum_lambda <- 10^minimum_log_lambda

fd_par_obj <- fdPar(basis,2,minimum_lambda)

sample_of_functions <- smooth.pos(time_points,observations,fd_par_obj)$Wfdobj

sample_of_functions$fdnames$time <- Year

sample_of_functions$fdnames$values <- "NO Levels"




plot(as.numeric(data$date),data$no,type="l",ylim=c(0,4))
plot(exp(eval.fd(seq(range(time_points)[1],range(time_points)[2],length.out=4000),sample_of_functions)),type="l")


plot(sample_of_functions)
lines(time_points,observations)

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












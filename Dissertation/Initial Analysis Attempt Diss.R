library(openair)
library(tidyverse)
library(fda)
library(lubridate)

unique_vector <- function(data){
  data %>%
    t() %>%
    as.vector() %>%
    unique()
}

to_vector <- function(data){
  data %>%
    t() %>%
    as.vector()
}


aurn_detailed <- importMeta(source = "aurn", all = TRUE)

nox_sites <- aurn_detailed %>%
  filter(variable == "NOx") %>%
  filter(end_date == "ongoing") %>%
  dplyr::select(code) %>%
  unique_vector()

data <- importAURN(site=nox_sites[1:2],year=2021,pollutant=c("nox"),meta=TRUE) %>%
        dplyr::rename(Emissions="nox")

data <- data %>%
  mutate(sec_since_1st_date = as.numeric(date)-as.numeric(min(date))) %>%
  mutate(year_since_1st_date = sec_since_1st_date/(60*60*24*365)) %>% 
  as.data.frame() 

ggplot(data,aes(x=year_since_1st_date,y=Emissions)) +
  geom_line() +
  facet_wrap(~site)


data_for_pivot <- data %>%
  dplyr::select(year_since_1st_date,site,Emissions) %>%
  dplyr::arrange(year_since_1st_date)

Emission_mat <- pivot_wider(data=data_for_pivot,
                            values_from=Emissions,
                            names_from=site,
                            names_prefix = "") 

stat_1 <- Emission_mat %>%
  dplyr::select(names_of_stations[1]) %>%
  to_vector()

stat_2 <- Emission_mat %>%
  dplyr::select(names_of_stations[2]) %>%
  to_vector()


dates <- Emission_mat %>%
  dplyr::select(year_since_1st_date) %>%
  to_vector()



Emission_mat <- Emission_mat %>%
  dplyr::select(-year_since_1st_date) %>%
  as.matrix()

rownames(Emission_mat) <- dates

max_year <- max(dates)
min_year <- min(dates)

year_breaks <- seq(min_year,max_year,(max_year-min_year)/100)

year_breaks <- c(year_breaks,max_year)


spline_basis <- create.bspline.basis(rangeval = c(min_year,max_year),length(year_breaks)+4-2,breaks=year_breaks)


names_of_stations <- data %>%
                      dplyr::select(site) %>%
                      unique_vector()

labels <- list("Years Since First Observation","Station"=names_of_stations,"Emissions of NOx")



date_fp_obj <- fdPar(spline_basis, 2, 0.01)


fd_obj <- smooth.basis(dates,Emission_mat,date_fp_obj)



loglam <- seq(-6,0,0.25)
gcvsave <- c()
dfsave <- c()

for (i in 1:length(loglam)){
  lambda_i <- 10^loglam[i]
  fdPar_i <- fdPar(spline_basis, 2, lambda_i)
  fdobj_i <- smooth.basis(dates,Emission_mat,fdPar_i)
  gcvsave <- c(gcvsave,sum(fdobj_i$gcv))
  dfsave <- c(dfsave,fdobj_i$df)
}


gcv_v_loglam <- data.frame(log_lam=loglam,gcv=gcvsave,df=dfsave)

ggplot(gcv_v_loglam,aes(x=log_lam,y=dfsave))+
  geom_point()+
  geom_line()





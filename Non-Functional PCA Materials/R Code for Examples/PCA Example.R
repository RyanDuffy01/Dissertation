library(tidyverse)
library(ggplot2)

data <- data.frame(
  geneone= c(10,11,8,3,2,1),
  genetwo= c(6,4,5,3,2.8,1)
)

n <- length(data$geneone)

mean_vector <- t(as.data.frame(colSums(data)/n))

rownames(mean_vector) <- c("Mean")

data_min_mean <- as.matrix(sweep(data, 2, mean_vector))

data_min_mean

cov_matrix <- t(data_min_mean) %*% data_min_mean

eigenvecs_and_vals <- eigen(cov_matrix)

eigvec_PC1 <- eigenvecs_and_vals$vectors[,1]

eigvec_PC2 <- eigenvecs_and_vals$vectors[,2]

linegrad_PC1 <- eigvec_PC1[2]/eigvec_PC1[1]

linegrad_PC2 <- eigvec_PC2[2]/eigvec_PC2[1]


ggplot(as.data.frame(data_min_mean),aes(x=geneone,y=genetwo))+
  geom_point()+
  geom_abline(slope=linegrad_PC1,intercept=0)+
  geom_abline(slope=linegrad_PC2,intercept=0)


prcomp(data)







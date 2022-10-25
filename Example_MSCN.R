
###########################################################################
###########################################################################
###                                                                     ###
###                               EXAMPLE                               ###
###                                                                     ###
###########################################################################
###########################################################################

##### Authors 
# Hung Tong, University of Alabama
# Cristina Tortora, San Jose State University

# 1 cluster with outliers in both principal components and 1 elliptical cluster
# n = 600, n1 = 420, n2 = 180
# 10% outliers

source("mscnm.R") ###Change path
load('data_1x.rdata')    # load objects data_1x and outliers_1x

####Visualize the data
plot(data_1x, col = outliers_1x + 1, pch = 16)

true_labels <- rep(1:2, times = c(420, 180))
plot(data_1x, col = true_labels, pch = 16)

library(mice)

### Adding missing values
set.seed(1234)
data_m_1x <- as.matrix(ampute(data_1x, prop = 0.1, mech = "MAR")$amp)

## Running MSCN with missing value
## It may require some time
mod <- mixture_incomplete_mscn(data_m_1x, G = 2, max_iter = 20)

##Results
plot(data_1x, col = mod$cluster)
table(true_labels, mod$cluster)

plot(data_1x, col = mod$outlier + 1)
table(outliers_1x, mod$outlier)

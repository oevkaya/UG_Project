
#For reproducibility fix the seed 
set.seed(2107)

# Some R packages  --------------------------------------------------------
# load packages
install.packages("VineCopula") 
install.packages("copula") 
library(copula)
library(VineCopula)

#Loading the data set
library(readr)
aksehirPET <- read_csv("aksehirPET.csv")

# Calculate water balance as a PREC - PET from loaded data set
aksehirPET[, 1] <- NULL
aksehirPET$WB <- aksehirPET$PREC - aksehirPET$PET

#Create the time series object for further analysis
DATA_RAWTS <- ts(aksehirPET, start = c(1950, 1), frequency = 12)
colnames(DATA_RAWTS)

# variable names, dimension, sample size
vnames <- colnames(DATA_RAWTS)
d <- length(vnames)
N <- dim(DATA_RAWTS)[1]
summary(DATA_RAWTS[, -c(1,2)])

# Renaming with station abbreviation
aksehirRAWTS <- DATA_RAWTS
aksehirRAWTS <- aksehirRAWTS[, -c(1,2)]
# visualize data (time series) using stl function
library(stats)
require(graphics)

# Some visualization regarding decomposition
plot(stl(aksehirRAWTS[, "PREC"], "per"), main="Decomposition of PREC data")
plot(stl(aksehirRAWTS[, "PET"], "per"), main="Decomposition of PET data")

#Visualize all time series data in one graph

# Test stationarity of standardized time series (%95) based on ADF, PP & KPSS test
library(tseries)
pvals_d0 <- matrix(NA, nrow=5, ncol=3);
colnames(pvals_d0) <- c("adf_p","pp_p","kpss_p")
for (j in 1:5) {
  adf <- adf.test(aksehirRAWTS[,j]);
  pp <- pp.test(aksehirRAWTS[,j]);
  kpss <- kpss.test(aksehirRAWTS[,j]);
  pvals_d0[j,] <- c(adf$p.val,pp$p.val,kpss$p.val)
}

# ACF/PACF of original time series for PREC and PET 
op <- par(mfrow=c(2,2),mar = c(6, 6, 2, 0.2))
acf(aksehirRAWTS[,"PREC"], ylab="ACF of PREC", main = "", las = 2)
pacf(aksehirRAWTS[,"PREC"], ylab="PACF of PREC", main = "", las = 2)
acf(aksehirRAWTS[,"PET"], ylab="ACF of PET", main = "", las = 2)
pacf(aksehirRAWTS[,"PET"], ylab="PACF of PET", main = "", las = 2)
par(op)

# Focus on the deseasonalised and pre-processed univariate time series based on use of function SUI

# All functions loading from SIndices
library(car)
path_to_Rfn <- c("SIndices/R/")             # set the path
source_files <- list.files(path_to_Rfn, "*.R$")  # locate all .R files
library(purrr)
map(paste0(path_to_Rfn, source_files), source)  # source all R scripts!

SI_out <- SUI(mts = aksehirRAWTS, tdg = 0L, ar = 1L, ma = 0L, sgn = c(1,-1, 1, -1, 1) , scale = c(3, 6, 9, 12),
              powertrans = TRUE, na.fill = FALSE)

# Hongyi notes: SUI(.)$mits gives the deseasonalised and temporally indep ts.
aksehir_DetDes_TS <- SI_out$mits
colnames(aksehir_DetDes_TS)

# Possible ARIMA - SARIMA model fit for climatic variables PREC, PET

# Box Pierce and Ljung-Box test results of residuals to check autocorrelation
BP.LB.tests <- matrix(0, nrow=5, ncol=2)

# Store Box-Pierce and Ljung-Box test p values
for (i in 1:5) {
  bp <- Box.test(aksehir_DetDes_TS[, i], lag = 1,type="Box-Pierce")
  lb <- Box.test(aksehir_DetDes_TS[, i], lag = 1, type = "Ljung-Box")
  BP.LB.tests[i,] <- c(bp$p.val,lb$p.val)
}

###########################################################################
########### COPULA ANALYSIS PART ####################################

# Store residuals from the above process, it is the input for COPULA analysis part
res_DATA = aksehir_DetDes_TS

# pseudo observations/copula data
U_DATA <- pobs(res_DATA)

pairs(as.copuladata(U_DATA), margins = "unif", pch = 1, gap = 0.3, cex.labels = 2)

# Bivariate independence test results for each pair
ind_test_PRCP_PET <- BiCopIndTest(U_DATA[,"PREC"], U_DATA[,"PET"])
ind_test_PRCP_PET

Bicop_PRCP_PET_AIC <- BiCopSelect(U_DATA[,"PREC"], U_DATA[,"PET"], 
                                  familyset = NA, selectioncrit = "AIC", indeptest = TRUE)
Bicop_PRCP_PET_AIC

Bicop_PRCP_PET_BIC <- BiCopSelect(U_DATA[,"PREC"], U_DATA[,"PET"], 
                                  familyset = NA, selectioncrit = "BIC", indeptest = TRUE)
Bicop_PRCP_PET_BIC

# use gaussian copula for PREC/PET data based on AIC/BIC value
nc <- normalCopula(dim = 2, dispstr = "un")
(ncf <- fitCopula(nc, cbind(U_DATA[,"PREC"], U_DATA[,"PET"])))

#Contour and surface plots of fitted normal copula based on above selection
fit_normal <- BiCop(family = Bicop_PRCP_PET_BIC$family, 
                    par = Bicop_PRCP_PET_BIC$par, 
                    par2 = Bicop_PRCP_PET_BIC$par2)

plot(fit_normal) # (marginal normal) contour plot
plot(fit_normal, margins = "unif") # contour plot of actual copula density
plot(fit_normal, type = "surf") # surface plot of actual copula density

nc@parameters <- ncf@estimate
nc@parameters # SAME as what we obtained above 

# Calculate CDF for pseudo observations
PREC_PET_po <- pCopula(cbind(U_DATA[,"PREC"], U_DATA[,"PET"]), nc)
hist(PREC_PET_po)

# transform to U(0,1)
PREC_PET_uni <- pobs(as.matrix(PREC_PET_po))
hist(PREC_PET_uni)

# transform to N(0,1)
st_PREC_PET <- qnorm(PREC_PET_uni)
hist(st_PREC_PET)

st_PREC_PET  <- ts(st_PREC_PET, start = c(1950, 1), frequency = 12, names = "PREC.PET")
plot(st_PREC_PET, type = "h")

# Generate Bivariate Drought Indices using different time scales on PREC and PET jointly modeled by Gaussian copula above
COPDI36912 <- par(mfrow=c(2,2))
# For 3 month
ts3 <- 3
st_PREC_PET_3M  <- c(rep(0, ts3 - 1), st_PREC_PET)

COPDI3 <- ts(rowSums(embed(st_PREC_PET_3M, ts3))/sqrt(ts3), start = c(1950, 1),
             frequency = 12, names = paste0("COPDI", ts3))

plot(COPDI3, type = "l", las = 2, col="red", ylab="", main = "COPDI3")

# For 6 month
ts6 <- 6
st_PREC_PET_6M  <- c(rep(0, ts6 - 1), st_PREC_PET)

COPDI6 <- ts(rowSums(embed(st_PREC_PET_6M, ts6))/sqrt(ts6), start = c(1950, 1),
             frequency = 12, names = paste0("COPDI", ts6))

plot(COPDI6, type = "l", las = 2, col="red", ylab="", main = "COPDI6")

# For 9 month
ts9 <- 9
st_PREC_PET_9M  <- c(rep(0, ts9 - 1), st_PREC_PET)

COPDI9 <- ts(rowSums(embed(st_PREC_PET_9M, ts9))/sqrt(ts9), start = c(1950, 1),
             frequency = 12, names = paste0("COPDI", ts9))

plot(COPDI9, type = "l", las = 2, col="red", ylab="", main = "COPDI9")

# For 12 month
ts12 <- 12
st_PREC_PET_12M  <- c(rep(0, ts12 - 1), st_PREC_PET)

COPDI12 <- ts(rowSums(embed(st_PREC_PET_12M, ts12))/sqrt(ts12), start = c(1950, 1),
             frequency = 12, names = paste0("COPDI", ts12))

plot(COPDI12, type = "l", las = 2, col="red", ylab="", main = "COPDI12")

par(COPDI36912)

#######################################################################
######### OTHER DROUGHT INDEX GENERATION ###########################

# For comparison Purposes, Calculation of SPI and SPEI using SPEI package in R
library(SPEI)
# Calculation of SPI and SPEI for 3,6,9,12 months
#For SPI calculations
SPI3 <- spi(aksehirRAWTS[, "PREC"], scale = ts3)
SPI3 <- ts(SPI3$fitted, start = c(1950, 1), frequency = 12, names = "SPI3")
plot(SPI3 , type = "l", las = 2, col="red", ylab="", main="SPI3")

# Compare with the output of SI_out 
# NOT EXACTLY THE SAME SINCE DATA PRE_PROCESSING PARTS ARE DIFFERENT
plot(SI_out$SI$PREC[, 1], type = "l", las = 2, col="red", ylab="", main="SI_out3")

# For 6 month
SPI6 <- spi(aksehirRAWTS[, "PREC"], scale = ts6)
SPI6 <- ts(SPI6$fitted, start = c(1950, 1), frequency = 12, names = "SPI6")
plot(SPI6 , type = "l", las = 2, col="red", ylab="", main="SPI6")

# For 9 month
SPI9 <- spi(aksehirRAWTS[, "PREC"], scale = ts9)
SPI9 <- ts(SPI9$fitted, start = c(1950, 1), frequency = 12, names = "SPI9")
plot(SPI9 , type = "l", las = 2, col="red", ylab="", main="SPI9")

# For 12 month
SPI12 <- spi(aksehirRAWTS[, "PREC"], scale = ts12)
SPI12 <- ts(SPI12$fitted, start = c(1950, 1), frequency = 12, names = "SPI12")

plot(SPI12 , type = "l", las = 2, col="red", ylab="", main="SPI12")


#For SPEI calculations
SPEI3 <- spei(aksehirRAWTS[, "WB"], scale = ts3)
SPEI3 <- ts(SPEI3$fitted, start = c(1950, 1), frequency = 12, names = "SPEI3")
plot(SPEI3 , type = "l", las = 2, col="red", ylab="", main="SPEI3")

# Compare with the output of SI_out 
# NOT EXACTLY THE SAME SINCE DATA PRE_PROCESSING PARTS ARE DIFFERENT
# plot(SI_out$SI$WB[, 1], type = "l", las = 2, col="red", ylab="", main="SI_out3")

# For 6 month
SPEI6 <- spei(aksehirRAWTS[, "WB"], scale = ts6)
SPEI6 <- ts(SPEI6$fitted, start = c(1950, 1), frequency = 12, names = "SPEI6")
plot(SPEI6 , type = "l", las = 2, col="red", ylab="", main="SPEI6")

# For 9 month
SPEI9 <- spei(aksehirRAWTS[, "WB"], scale = ts9)
SPEI9 <- ts(SPEI9$fitted, start = c(1950, 1), frequency = 12, names = "SPEI9")
plot(SPEI9 , type = "l", las = 2, col="red", ylab="", main="SPEI9")

# For 12 month
SPEI12 <- spei(aksehirRAWTS[, "WB"], scale = ts12)
SPEI12 <- ts(SPEI12$fitted, start = c(1950, 1), frequency = 12, names = "SPEI12")

plot(SPEI12 , type = "l", las = 2, col="red", ylab="", main="SPEI12")


# Compare our calculation with SMI function result on PREC and PET 
SM_PRPE <- SMI(mts = DATA_RAWTS[, c('PREC', 'PET')], ar = 1, ma = 0, sgn = c(1,-1), scale = c(3, 6, 9, 12),
              order = NULL, familyset = 1:6, indeptest = TRUE,
              powertrans = TRUE, na.fill = FALSE,
              method = "M", MVNtest = FALSE)
# Fitted Copula model seems roughly the same 
SM_PRPE$RVM
# Store the time series object of multiplicative standardized multivariate index
SM_PRPE_ts <- SM_PRPE$SIM

# Compare our calculation with SMI function result on PREC and RH 
SM_PRRH <- SMI(mts = DATA_RAWTS[, c('PREC', 'RH')], ar = 1, ma = 0, sgn = c(1,1), scale = c(3, 6, 9, 12),
              order = NULL, familyset = 1:6, indeptest = TRUE,
              powertrans = TRUE, na.fill = FALSE,
              method = "M", MVNtest = FALSE)
# Fitted Copula model seems roughly the same 
SM_PRRH$RVM
# Store the time series object of multiplicative standardized multivariate index
SM_PRRH_ts <- SM_PRRH$SIM

# Compare our calculation with SMI function result on PREC and TMean 
SM_PRTM <- SMI(mts = DATA_RAWTS[, c('PREC', 'TMean')], ar = 1, ma = 0, sgn = c(1,-1), scale = c(3, 6, 9, 12),
               order = NULL, familyset = 1:6, indeptest = TRUE,
               powertrans = TRUE, na.fill = FALSE,
               method = "M", MVNtest = FALSE)
# Fitted Copula model seems roughly the same 
SM_PRTM$RVM
# Store the time series object of multiplicative standardized multivariate index
SM_PRTM_ts <- SM_PRTM$SIM

# Compare our calculation with SMI function result on PREC, PET and RH
SM_PRPERH <- SMI(mts = DATA_RAWTS[, c('PREC', 'PET', 'RH')], ar = 1, ma = 0, sgn = c(1,-1, 1), scale = c(3, 6, 9, 12),
               order = NULL, familyset = 1:6, indeptest = TRUE,
               powertrans = TRUE, na.fill = FALSE,
               method = "M", MVNtest = FALSE)
# Fitted Copula model seems roughly the same 
SM_PRPERH$RVM
# Store the time series object of multiplicative standardized multivariate index
SM_PRPERH_ts <- SM_PRPERH$SIM

saved_DATA = aksehirPET[, c(1,2)]

saved_DATA <- cbind(saved_DATA, SPI3, SPEI3, SM_PRPE_ts[,1], SM_PRRH_ts[,1], SM_PRTM_ts[,1], SM_PRPERH_ts[,1],
                    SPI6, SPEI6, SM_PRPE_ts[,2], SM_PRRH_ts[,2], SM_PRTM_ts[,2], SM_PRPERH_ts[,2], 
                    SPI9, SPEI9, SM_PRPE_ts[,3], SM_PRRH_ts[,3], SM_PRTM_ts[,3], SM_PRPERH_ts[,3],
                    SPI12, SPEI12, SM_PRPE_ts[,4], SM_PRRH_ts[,4], SM_PRTM_ts[,4], SM_PRPERH_ts[,4])
colnames(saved_DATA) <- c("Month", "Year", "SPI3", "SPEI3", "SMI(PREC-PET3)", "SMI(PREC-RH3)", "SMI(PREC-TMean3)", "SMI(PREC-PET-RH3)", 
                          "SPI6", "SPEI6", "SMI(PREC-PET6)", "SMI(PREC-RH6)", "SMI(PREC-TMean6)", "SMI(PREC-PET-RH6)", 
                          "SPI9", "SPEI9", "SMI(PREC-PET9)", "SMI(PREC-RH9)", "SMI(PREC-TMean9)", "SMI(PREC-PET-RH9)", 
                          "SPI12", "SPEI12", "SMI(PREC-PET12)", "SMI(PREC-RH12)", "SMI(PREC-TMean12)", "SMI(PREC-PET-RH12)")

save(saved_DATA, file = "data.RData")
# load("data.RData")
write.csv(saved_DATA, "saved_DATA.csv")

# import the r script "Run_Theory.R" for the calculations of 
# Duration, Severity, Mean Intensity,Peak intensity and Interarrival time
install.packages("ggpubr") 
library(ggpubr)
source("Run_Theory.R")

# Change the working directory to '/SMI_CSV_4*6'
# in order to save all the csv file in a folder
setwd(paste(getwd(),'/DroughtChar .csv files', sep = ''))

# Initialise values
lst_SMI <- c("SPI", "SPEI", "SMI(PREC-PET)", "SMI(PREC-RH)", "SMI(PREC-TMean)", "SMI(PREC-PET-RH)")
lst_SMI <- rep(lst_SMI, 4)
lst_DSI <- c()
for (j in 1:24) {
  scale <- ((j+5)%/%6)*3 # scale = 3,6,9,12
  
  if (j%%6==1) {
    lst_SMI[j] <- paste(lst_SMI[j], scale, sep = '')
  } else if (j%%6==2) {
      lst_SMI[j] <- paste(lst_SMI[j], scale, sep = '')
  } else if (j%%6==3) {
        lst_SMI[j] <- paste(c(substr(lst_SMI[j], 1, 12), substr(lst_SMI[j], 13,14)), collapse=toString(scale))
  } else if (j%%6==4) {
          lst_SMI[j] <- paste(c(substr(lst_SMI[j], 1, 11), substr(lst_SMI[j], 12,13)), collapse=toString(scale))
  } else if (j%%6==5) {
            lst_SMI[j] <- paste(c(substr(lst_SMI[j], 1, 14), substr(lst_SMI[j], 15,16)), collapse=toString(scale))
  } else if (j%%6==0) {
              lst_SMI[j] <- paste(c(substr(lst_SMI[j], 1, 15), substr(lst_SMI[j], 16,17)), collapse=toString(scale))
  }
  
  SMI_val <- run_DSI(na.remove(saved_DATA[,lst_SMI[j]]))
  lst_DSI <- cbind(lst_DSI, SMI_val)
  colnames(lst_DSI)[j] <- lst_SMI[j]
  # Add a NA value in the last row of Interarrival when its number of values 
  # is small than other characteristics
  if (length(lst_DSI[,j][[7]]) != length(lst_DSI[,j][[6]])){
    lst_DSI[,j][[7]] <- append(lst_DSI[,j][[7]], NA)
  }
  # Replace with a NA value in the last row of Interarrival when it's not an NA
  lst_DSI[,j][[7]][length(lst_DSI[,j][[7]])] <- NA
  
  # Export the result in the form of CSV
  write.csv(lst_DSI[,j], paste(lst_SMI[j],'.csv', sep = ''))
}

# Change the working directory back to the original
setwd(substr(getwd(),1,nchar(getwd())-23))

# Plotting the bivariate index coming from the package function SMI
plot(SM_PRPE$SIM[, 1], type = "l", las = 2, col="red", ylab="")
lines(COPDI3, type = "l", las = 2, col="blue", ylab="")
# plot(SM_PRPE$SIM[, 1], type = "l", las = 2, col="red", ylab="")
# lines(COPDI6, type = "l", las = 2, col="blue", ylab="")
# lines(COPDI9, type = "l", las = 2, col="blue", ylab="")
# lines(COPDI12, type = "l", las = 2, col="blue", ylab="")


# fitdistr()

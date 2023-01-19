#' Compute drought duration and severity based on run theory
#' 
#' The input data is monthly drought indices.
#' Duration is defined as the length of consecutive time series 
#' when drought index is below the threshold value (e.g., -1). 
#' Severity is defined as the summation of drought index below the threshold.
#' This analysis based on run theory is also referred to as threshold level method.
#' Here the standardized drought index (SDI) is used as 
#' the example to compute the drought characteristics. 
#' Other univariate and multivariate drought indices can also be used.

#### applied the run theory to a time series ()
## SDI: a numeric vector for Standardized Drought Index with no NA values
## threshold: a numeric value in which the features (below) of run theory is measured
run_DSI <- function(SDI, threshold = -.5)
{
  
  dataBase <- data.frame(DroughIndex = SDI) %>%
    transform(masked = ifelse(DroughIndex >= threshold, 1, 0)) %>%
    transform(index = cumsum(masked), index_rev = cumsum(abs(masked-1)))
  
  ####  Duration, Severity, Intensity ####  
  
  dataBase[dataBase$masked == 0, ] %>%
    by(., .$index, function(z){
      
      data.frame(D = dim(z)[1],
                 S = abs(sum(z$DroughIndex)),
                 MI = abs(sum(z$DroughIndex))/dim(z)[1],
                 PI = abs(min(z$DroughIndex)),
                 date_ini = row.names(z)[1],
                 date_fin = row.names(z)[nrow(z)])
      
    }) %>% do.call(rbind, .) -> df1
  
  
  ####  Interarrival ####  
  
  dataBase[dataBase$masked != 0, ] %>%
    by(., .$index_rev, function(z){
      
      data.frame(Int = dim(z)[1])
      
    }) %>% unlist() -> Int
  
  # first condition 
  
  if ((dataBase$masked[1]) == 1){
    
    Int <- Int[-1] 
    
  } else {
    
    Int <- Int
    
  }
  
  # second contidition
  
  if (dataBase$masked[length(dataBase$masked)] == 1) {
    
    Int <- Int + df1$D
    
  } else { 
    
    n <- c(Int, 0 ) + df1$D
    Int <- n[-length(n)]
    
  }
  
  
  return(list(Duration = as.numeric(df1$D),
              Severity = as.numeric(df1$S),
              Mean_Intensity = as.numeric(df1$MI),
              Peak_Intensity = as.numeric(df1$PI),
              Date_Ini_Ev = as.character(df1$date_ini),
              Date_Fin_Ev = as.character(df1$date_fin),
              Interarrival = as.numeric(Int)))
  
}

drought_char <- run_DSI(SI_out$SI$PREC[, 1], threshold = 0)

# Plotting the Drought Characteristics

## plot
library(ggplot2)

# Function takes the drought characteristic properties
# Hongyi comments: need to run library(forcats) here
library(forcats)
plot_SDI <- function(Drought_char){
 df <-  lapply(Drought_char %>% names(), 
         function(x){
           data.frame(value = Drought_char[[x]], 
                      feature = x) }) %>% 
    .[c(1, 2, 3, 4)] %>% # deleting Dates output to plot
    do.call(rbind, .)
    
    ggplot(df) + 
    geom_histogram(aes(x = value), color="darkblue", fill="lightblue") +
    geom_density(aes(x = value), alpha=.2, fill="#FF6666") +
    facet_wrap(~ fct_relevel(feature, "Duration", "Severity", "Mean_Intensity", "Peak_Intensity"), "free", nrow = 2, ncol = 2) + 
    ggtitle("Drought Characteristics Summary") +
    xlab("Months") + ylab("Frequency") + theme_pubclean() + 
    theme(strip.text.x = element_text( size = 12, color = "black", face = "bold.italic") )
}

plot_SDI(Drought_char = drought_char)

# Collecting 
main_char <- cbind(drought_char$Duration, drought_char$Severity, drought_char$Mean_Intensity, drought_char$Peak_Intensity, drought_char$Interarrival)
colnames(main_char) <- c('Duration', 'Severity', 'Mean_Intensity', 'Peak_Intensity', 'Interarrival Time')
# Randomization on Duration variable
main_char[, 'Duration'] <- main_char[, 'Duration'] + runif(dim(main_char)[1], min = -0.5, max = 0.5 )

# First pair copula Visualization, pairs copula data 
udat_char <- pobs(main_char[, 1:4])
pairs.copuladata(udat_char, labels = colnames(main_char))

# Storing all the necessary drought characteristics 
DATA_Dchar <- main_char
write.csv(DATA_Dchar, "aksehir_DChar.csv")





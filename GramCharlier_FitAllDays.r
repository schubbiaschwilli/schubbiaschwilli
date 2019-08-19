load("...\\EurexOptionsDaxPlus.RData")

BSMOption<-function(PutCallFlag,S,T,K,r,sigma,Greek)
{
   if(PutCallFlag=="Call" && Greek=="Price" && T>0)
   {
      d1=(log(S/K)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
      d2=d1 - sigma * sqrt(T)
      return(S * pnorm(d1, mean=0, sd=1) - K * exp(-r*T) * pnorm(d2, mean=0, sd=1)) 
   }
   if(PutCallFlag=="Call" && Greek=="Price" && T==0)
   {return(pmax(S-K,0))} 
   if(PutCallFlag=="Put" && T>0)
   {return(-1 * BSMOption("Call",S,T,K,r,-1 * sigma,Greek))}
   if(PutCallFlag=="Put" && Greek=="Price" && T==0)
   {return(pmax(K-S,0))}
}

GramCharlierOption<-function(PutCallFlag,S,T,K,r,VolatilityPerYear,OnePeriodSkewness,OnePeriodKurtosis,Greek)
{
   Skewness <- OnePeriodSkewness / sqrt(T)
   Kurtosis <- OnePeriodKurtosis / T
   Volatility <- VolatilityPerYear * sqrt(T/12)
   OnePeriodr <- r /12
   
   if(PutCallFlag=="Call" && Greek=="Price" && T>0)
   {
      d <- (log(S / K) + T * OnePeriodr + Volatility^2 / 2) / Volatility
      return(S * pnorm(d, mean=0, sd=1) - K*exp(-OnePeriodr*T)*pnorm(d - Volatility, mean=0, sd=1) + S * dnorm(d, mean=0, sd=1) * Volatility * (Skewness/6 * (2*Volatility - d) - (Kurtosis/24) * (1 - d^2 + 3 * d * Volatility - 3 * Volatility^2)))
   }
   if(PutCallFlag=="Call" && Greek=="Price" && T==0)
      {return(pmax(S-K,0))} 
   if(PutCallFlag=="Put" && T>0)
      {return(GramCharlierOption("Call",S,T,K,r,VolatilityPerYear,OnePeriodSkewness,OnePeriodKurtosis,"Price") + K * exp(-OnePeriodr * T)-S)}
   if(PutCallFlag=="Put" && Greek=="Price" && T==0)
      {return(pmax(K-S,0))}
}

# Create df for results
GCModel <- unique(EurexOptionsDax["Handelstag"])
GCModel <- data.frame(GCModel)
GCModel$VolatilityPerYear <- -1
GCModel$OnePeriodSkewness <- -1
GCModel$OnePeriodKurtosis <- -1
GCModel$RMSE_S <- -1

### Fit
library(minpack.lm)
IterationControl<-nls.control(maxiter=2^10)
startvalue<-list(VolatilityPerYear = 0.8, OnePeriodSkewness=-0.6, OnePeriodKurtosis=0) 

for(i in 1:nrow(GCModel))
{
   data <- subset(EurexOptionsDax, OptionType == 'Call' & Handelstag == as.POSIXct(GCModel[i,]$Handelstag, tz="UTC") & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)
   
   # Fit via nlsLM
   Price<-data$TaeglicherAbrechnungspreis
   OptionType<-data[1,"OptionType"]
   S<-data[1,"SchlusspreisBasiswert"]
   r<-log(1+data[1,"EONIA"])
   T<-data$t_delta
   K<-data$StrikePrice

   result<-nlsLM(Price~GramCharlierOption(OptionType,S,T,K,r,VolatilityPerYear,OnePeriodSkewness,OnePeriodKurtosis,'Price'), start=startvalue, control=IterationControl)

   GCModel[i,"VolatilityPerYear"]<-as.list(coef(result))$VolatilityPerYear
   GCModel[i,"OnePeriodSkewness"]<-as.list(coef(result))$OnePeriodSkewness
   GCModel[i,"OnePeriodKurtosis"]<-as.list(coef(result))$OnePeriodKurtosis
   
   GCModel[i,"RMSE_S"] <- sqrt(sum((Price - GramCharlierOption(OptionType,S,T,K,r,GCModel[i,"VolatilityPerYear"],GCModel[i,"OnePeriodSkewness"],GCModel[i,"OnePeriodKurtosis"],'Price'))^2) / nrow(data)) / S
   
   # Calc predicted RMSE/S
   if(i>1)
   {
      GCModel[i,"RMSE_S_Previous_Vs_Current"] <- sqrt(sum((Price - GramCharlierOption(OptionType,S,T,K,r,GCModel[i-1,"VolatilityPerYear"],GCModel[i-1,"OnePeriodSkewness"],GCModel[i-1,"OnePeriodKurtosis"],'Price'))^2) / nrow(data)) / S
   }
}

# Create plots
library("reshape2")
library("ggplot2")

# Plot parameter timeseries
tmp <- melt(GCModel[,c(1,2,3,4)], id="Handelstag")
tmp$Handelstag <- as.Date(tmp$Handelstag)
ggplot(data=tmp, aes(x=Handelstag, y=value, color=factor(variable, labels = c("Volatility Per Year","One Period Skewness","One Period Kurtosis")))) + 
   geom_line(size = 1) + ggtitle("Gram-Charlier-Model-Parameter") + scale_x_date(date_minor_breaks = "1 month") + theme_bw(base_size=16) + theme(legend.position="bottom") + 
   theme(legend.title=element_blank())

# Model RMSE_S v. RMSE_S VolaFittedPrevious
tmp <- melt(GCModel[,c(1,5,6)], id="Handelstag")  # convert to long format
tmp$Handelstag <- as.Date(tmp$Handelstag)
ggplot(data=tmp, aes(x=Handelstag, y=value, color=factor(variable, labels = c("RMSE / S","RMSE / S with predicted parameters from previous market data")))) + geom_line(size = 1) + ggtitle("Gram-Charlier-Model") + scale_x_date(date_minor_breaks = "1 month") + theme_bw(base_size=16) + theme(legend.position="bottom") + theme(legend.title=element_blank())

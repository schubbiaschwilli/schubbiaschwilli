load("...\\EurexOptionsDaxPlus.RData")
load("...\\EurexOptionsDaxMJDParameterResults.RData")

library(minpack.lm)
library(scatterplot3d)
library("reshape2")
library("ggplot2")

BSMOption <- function(PutCallFlag,S,T,K,r,sigma,Greek)
{  # Vergleiche 'Four Things You Might Not Know About The Black-Scholes-Formula'
   if(PutCallFlag=="Call" && Greek=="Price")
   {
      # Price Call vor Fälligkeit
      d1<-(log(S/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
      d2<-d1 - sigma * sqrt(T)
      value<-S*pnorm(d1, mean=0, sd=1) - K*exp(-r*T)*pnorm(d2, mean=0, sd=1)
      return(value) 
   }
   if(PutCallFlag=="Put")
   {
      return(-1*BSMOption("Call",S,T,K,r,-1*sigma,Greek))
   }
}

MJDOption <- function(PutCallFlag,S,T,K,r,sigma,lambda,m,delta,Greek,N)
{
   value<-0

   for(i in 0:N) 
   {
      S_n<-S*exp(i*m+i*delta^2/2-lambda*exp(m+delta^2/2)*T+lambda*T);
      sigma_n<-sqrt(sigma^2+i*delta^2/T);
      value<-value+exp(-lambda*T)*(lambda*T)^i/gamma(i+1)*BSMOption(PutCallFlag,S_n,T,K,r,sigma_n,Greek);
   }
   return(value)
}

### Fit ###
#for(i in 1:nrow(MJDParameter))
for(i in 1:1)
{
   # Get Data
   MarketData <- subset(EurexOptionsDax, OptionType == 'Call' & Handelstag == as.POSIXct(MJDParameter$Handelstag[i], tz="UTC") & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)
   
   Price<-MarketData$TaeglicherAbrechnungspreis
   OptionType<-MarketData[i,"OptionType"]
   S<-MarketData[i,"SchlusspreisBasiswert"]
   r<-log(1+MarketData[i,"EONIA"])
   T<-MarketData$t_delta
   K<-MarketData$StrikePrice
   
   # Fit via nlsLM
   IterationControl<-nls.control(maxiter=2^10)
   startvalue<-list(sigma=MJDParameter$sigma_init[i], lambda=MJDParameter$lambda_init[i], m=MJDParameter$m_init[i], delta=MJDParameter$delta_init[i])
   
   result<-nlsLM(Price~MJDOption(OptionType,S,T,K,r,sigma,lambda,m,delta,'Price',170), start=startvalue, control=IterationControl)
   
   MJDParameter$sigma[i] <- summary(result)$parameters[1]
   MJDParameter$lambda[i] <- summary(result)$parameters[2]
   MJDParameter$m[i] <- summary(result)$parameters[3]
   MJDParameter$delta[i] <- summary(result)$parameters[4]
   MJDParameter$NumberOf[i] <- nrow(MarketData)
   MJDParameter$RSS[i] <- deviance(result)
}

### Plot Single Day ###
MJDPrice <- function(K, T) 
{MJDOption(PutCallFlag, S, T, K, r, sigma, lambda ,m , delta, Greek="Price", N=170)}

#for(i in 1:nrow(MJDParameter))
for(i in 1:1)
{
   # Get Data
   MarketData <- subset(EurexOptionsDax, OptionType == 'Call' & Handelstag == as.POSIXct(MJDParameter$Handelstag[i], tz="UTC") & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)
   
   OptionType<-MarketData[i,"OptionType"]
   S<-MarketData[i,"SchlusspreisBasiswert"]
   r<-log(1+MarketData[i,"EONIA"])
   PutCallFlag<- OptionType[1]
   
   sigma<-MJDParameter$sigma[i]; #Volatility
   lambda<-MJDParameter$lambda[i];  #Expected average number of events / Intensität
   m<-MJDParameter$m[i];	#Expected value / Erwartungswert
   delta<-MJDParameter$delta[i];	#Std.Dev. of jumps / Standardabweichung der Sprünge
   
   ### Def Plots
   layout(matrix(c(1,1,1,1,1,2,2,2,2),1))
   title<-paste("Handelstag: ", MJDParameter$Handelstag[i], sep="")
   subtitle<-paste("S:",S,"; r:", round(r, digits=4),"; sigma:", round(sigma, digits=4),"; lambda:", round(lambda, digits=4),"; m:", round(m, digits=4),"; delta:", round(delta, digits=4), "; RSS:", format(round(MJDParameter$RSS[i], digits=4), big.mark=","), sep="")
   ### Volatility Surface 
   # Data 
   Strike<-seq(min(MarketData$StrikePrice), max(MarketData$StrikePrice), length=50)
   t<-seq(min(MarketData$t_delta), max(MarketData$t_delta), length=50)
   Price <- outer(Strike, t, MJDPrice)
   # Plot Price: Scatter vs Surface/ persp(x, y, z,...)
   plot<-persp(Strike, t, Price, theta=45, phi=15, expand=0.5, col='gray', main = NULL, ticktype="detailed", sub=NULL)
   # Add real prices
   ObservedPrices <- trans3d(MarketData$StrikePrice, MarketData$t_delta, MarketData$TaeglicherAbrechnungspreis, plot)
   CalculatetPrices <- trans3d(MarketData$StrikePrice, MarketData$t_delta, MJDPrice(MarketData$StrikePrice, MarketData$t_delta), plot)
   points(ObservedPrices, col="red", pch=16)
   segments(ObservedPrices$x, ObservedPrices$y, CalculatetPrices$x, CalculatetPrices$y)
   ### Deltas 
   deltas<-MarketData$TaeglicherAbrechnungspreis - MJDOption(OptionType, S, MarketData$t_delta, MarketData$StrikePrice, MarketData$EONIA, sigma, lambda, m, delta, Greek="Price", N=170)
   scatterplot3d(MarketData$StrikePrice, MarketData$t_delta, deltas, pch=16, highlight.3d=TRUE, angle=45, type="h", box=FALSE, main=NULL, xlab="Strike", ylab="t", zlab="Delta MJD- vs Observed Price")
   
   mtext(title, side=3, line=-2, outer=TRUE, cex=1)
   mtext(subtitle, side=1, line=-2, outer=TRUE, cex=1)
}

# Plot Parameter over time
par(mfrow=c(4,1), mar=c(2, 2, 0.5, 0.1))
plot(y=MJDParameter$sigma, x=MJDParameter$Handelstag, type='l', ylab="sigma", xlab="", col='red')
plot(y=MJDParameter$lambda, x=MJDParameter$Handelstag, type='l', ylab="lambda", xlab="", col='black')
plot(y=MJDParameter$m, x=MJDParameter$Handelstag, type='l', ylab="m", xlab="", col='blue')
plot(y=MJDParameter$delta, x=MJDParameter$Handelstag, type='l', ylab="delta", xlab="", col='green')

# Plot RSS over time
plot(y=MJDParameter$RSS, x=MJDParameter$Handelstag, type='l', ylab="RSS", xlab="", col='black')

# Create Scatterplot
pairs(MJDParameter[, c("sigma", "lambda", "m", "delta")], pch = 21)

# Correlation table
cor(MJDParameter[c("sigma","lambda","m","delta")])

# Calc RMSE/S
#for(i in 1:nrow(MJDParameter))
for(i in 1:1)
{
   MarketData <- subset(EurexOptionsDax, OptionType == 'Call' & Handelstag == as.POSIXct(MJDParameter$Handelstag[i], tz="UTC") & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)
   
   OptionType<-MarketData[i,"OptionType"]
   S<-MarketData[i,"SchlusspreisBasiswert"]
   r<-log(1+MarketData[i,"EONIA"])

   # Calc predicted RMSE/S
   MJDParameter[i,"RMSE_S"] <- sqrt(sum((MarketData$TaeglicherAbrechnungspreis - MJDOption(OptionType, S, MarketData$t_delta, MarketData$StrikePrice, r, MJDParameter[i,]$sigma, MJDParameter[i,]$lambda, MJDParameter[i,]$m, MJDParameter[i,]$delta, Greek="Price", N=170))^2) / nrow(MarketData)) / S
   
   # Calc predicted RMSE/S
   if(i>1)
   {
      MJDParameter[i,"RMSE_S_Previous_Vs_Current"] <- sqrt(sum((MarketData$TaeglicherAbrechnungspreis - MJDOption(OptionType, S, MarketData$t_delta, MarketData$StrikePrice, r, MJDParameter[i-1,]$sigma, MJDParameter[i-1,]$lambda, MJDParameter[i-1,]$m, MJDParameter[i-1,]$delta, Greek="Price", N=170))^2) / nrow(MarketData)) / S
   }
}

# Model RMSE_S v. RMSE_S  Fitted Previous
tmp <- melt(MJDParameter[,c(1,13,14)], id="Handelstag")  # convert to long format
tmp$Handelstag <- as.Date(tmp$Handelstag)
ggplot(data=tmp, aes(x=Handelstag, y=value, color=factor(variable, labels = c("RMSE / S","RMSE / S with predicted parameters from previous market data")))) + geom_line(size = 1) + ggtitle("Merton-Jump-Diffusion-Model") + scale_x_date(date_minor_breaks = "1 month") + theme_bw(base_size=16) + theme(legend.position="bottom") + theme(legend.title=element_blank())

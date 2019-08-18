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

load("D:\\h_da\\_Blog\\Themen\\Options\\BSM_ImpliedVola_AllSingleOptions\\CreataDataFile\\EurexOptionsDaxPlus.RData")

# Create df for results
GCModel <- unique(EurexOptionsDax["Handelstag"])
GCModel <- data.frame(GCModel)
GCModel$VolatilityPerYear <- -1
GCModel$OnePeriodSkewness <- -1
GCModel$OnePeriodKurtosis <- -1
GCModel$RMSE_S <- -1

### Fit

print(paste("Start : " ,Sys.time(), sep = " : "))

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

print(paste("Stop : " ,Sys.time(), sep = " : "))

# Create plots
library("reshape2")
library("ggplot2")
library("latex2exp")
library("scatterplot3d")


###

#ggplot(GCModel[,c(1,2)], aes(x = Handelstag, y = VolatilityPerYear)) + theme_bw() + geom_line()
#ggplot(GCModel[,c(1,3)], aes(x = Handelstag, y = OnePeriodSkewness)) + theme_bw() + geom_line()
#ggplot(GCModel[,c(1,4)], aes(x = Handelstag, y = OnePeriodKurtosis)) + theme_bw() + geom_line()

png(filename="D:/h_da/_Blog/Themen/GramCharlier/GramCharlierModel_TS_Parameter.png", width=1200, height=800, pointsize=40)
par(mfrow=c(3,1)) 
tmp <- melt(GCModel[,c(1,2,3,4)], id="Handelstag")
tmp$Handelstag <- as.Date(tmp$Handelstag)
ggplot(data=tmp, aes(x=Handelstag, y=value, color=factor(variable, labels = c("Volatility Per Year","One Period Skewness","One Period Kurtosis")))) + 
   geom_line(size = 1) + ggtitle("Gram-Charlier-Model-Parameter") + scale_x_date(date_minor_breaks = "1 month") + theme_bw(base_size=16) + theme(legend.position="bottom") + 
   theme(legend.title=element_blank())
dev.off()

### Scatterplot & Correlation parameter
pairs(GCModel[,c(2,3,4)])
cor(GCModel[,c(2,3,4)], use="pairwise.complete.obs")

### Scatterplot & Correlation delta
delta <- cbind(diff(GCModel$VolatilityPerYear, lag = 1), diff(GCModel$OnePeriodSkewness, lag = 1), diff(GCModel$OnePeriodKurtosis, lag = 1))
colnames(delta) <- c('Volatility Per Year', 'One Period Skewness', 'One Period Kurtosis')

pairs(delta, pch = 21)

cor(delta)

# Model RMSE_S v. RMSE_S VolaFittedPrevious
png(filename="D:/h_da/_Blog/Themen/GramCharlier/GramCharlierModel_Results.png", width=1200, height=800, pointsize=40)
tmp <- melt(GCModel[,c(1,5,6)], id="Handelstag")  # convert to long format
tmp$Handelstag <- as.Date(tmp$Handelstag)
ggplot(data=tmp, aes(x=Handelstag, y=value, color=factor(variable, labels = c("RMSE / S","RMSE / S with predicted parameters from previous market data")))) + geom_line(size = 1) + ggtitle("Gram-Charlier-Model") + scale_x_date(date_minor_breaks = "1 month") + theme_bw(base_size=16) + theme(legend.position="bottom") + theme(legend.title=element_blank())
dev.off()

### Compare with DVF 0
load("D:\\h_da\\_Blog\\Themen\\Options\\DeterministicVolatilityFunction\\EurexOptionsDax_DVF_results.RData")
png(filename="D:/h_da/_Blog/Themen/GramCharlier/GramCharlier_vs_DVF_Model_0.png", width=1200, height=800, pointsize=40)
tmp <- rbind(melt(GCModel[,c(1,5,6)], id="Handelstag"), melt(DVF_results[,c(1:3)], id="Handelstag"))
tmp$Handelstag <- as.Date(tmp$Handelstag)
ggplot(data=tmp, aes(x=Handelstag, y=value, color=factor(variable, labels = c("Gram Charlier","Gram Charlier with predicted parameters from previous market data","DVF Model 0","DVF Model 0 with predicted parameters from previous market data")))) + 
   geom_line(size = 1) + ggtitle("RMSE / S") + scale_x_date(date_minor_breaks = "1 month") + theme_bw(base_size=16) + theme(legend.position="bottom") + 
   theme(legend.title=element_blank())
dev.off()

### Compare with DVF 1
load("D:\\h_da\\_Blog\\Themen\\Options\\DeterministicVolatilityFunction\\EurexOptionsDax_DVF_results.RData")
png(filename="D:/h_da/_Blog/Themen/GramCharlier/GramCharlier_vs_DVF_Model_1.png", width=1200, height=800, pointsize=40)
tmp <- rbind(melt(GCModel[,c(1,5,6)], id="Handelstag"), melt(DVF_results[,c(1,4,5)], id="Handelstag"))
tmp$Handelstag <- as.Date(tmp$Handelstag)
ggplot(data=tmp, aes(x=Handelstag, y=value, color=factor(variable, labels = c("Gram Charlier","Gram Charlier with predicted parameters from previous market data","DVF Model 1","DVF Model 1 with predicted parameters from previous market data")))) + 
   geom_line(size = 1) + ggtitle("RMSE / S") + scale_x_date(date_minor_breaks = "1 month") + theme_bw(base_size=16) + theme(legend.position="bottom") + 
   theme(legend.title=element_blank())
dev.off()

### Price ###
#day <- as.POSIXct("2016-12-16", tz="UTC")
day <- as.POSIXct("2016-5-25", tz="UTC")
data <- subset(EurexOptionsDax, OptionType == 'Call' & Handelstag == day & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)

# Fit via nlsLM
Price<-data$TaeglicherAbrechnungspreis
OptionType<-data[1,"OptionType"]
S<-data[1,"SchlusspreisBasiswert"]
r<-log(1+data[1,"EONIA"])
T<-data$t_delta
K<-data$StrikePrice

IterationControl<-nls.control(maxiter=2^10)
startvalue<-list(VolatilityPerYear = 0.3, OnePeriodSkewness=0, OnePeriodKurtosis=0) 
result<-nlsLM(Price~GramCharlierOption(OptionType,S,T,K,r,VolatilityPerYear,OnePeriodSkewness,OnePeriodKurtosis,'Price'), start=startvalue, control=IterationControl)

VolatilityPerYear <- as.list(coef(result))$VolatilityPerYear
OnePeriodSkewness <- as.list(coef(result))$OnePeriodSkewness
OnePeriodKurtosis <- as.list(coef(result))$OnePeriodKurtosis

### Price
GramCharlierPrice <- function(K, T)
{GramCharlierOption(OptionType,S,T,K,r,VolatilityPerYear[i],OnePeriodSkewness[i],OnePeriodKurtosis[i],Greek="Price")}

main <- paste("Gram-Charlier-Model; Handelstag ", format(day, "%d-%m-%Y"), sep='')

StrikePrice <- seq(min(data$StrikePrice), max(data$StrikePrice), length.out=51)
t_delta <- seq(min(data$t_delta), max(data$t_delta), length.out=51)
Price <- outer(StrikePrice, t_delta, "GramCharlierPrice")

p <- persp(StrikePrice, t_delta, Price , theta=30, phi=30, col="lightgrey", expand = 0.5, shade = 0.2, xlab="Strike", ylab="Maturity", zlab="Price", ticktype = "detailed", main=main)
obs <- trans3d(data$StrikePrice, data$t_delta, data$TaeglicherAbrechnungspreis, p)
points(obs, col="red",pch=16)
pred <- trans3d(data$StrikePrice, data$t_delta, GramCharlierPrice(data$StrikePrice, data$t_delta), p)
segments(obs$x, obs$y, pred$x, pred$y)


###
### Vola
fImpliedVola<-function(sigma,Price,OptionType,S,T,K,r)
   # the function for which the root is sought
{return(BSMOption(OptionType,S,T,K,r,sigma,"Price")-Price)}

ImpliedVola<-function(Price,OptionType,S,T,K,r)
{return(uniroot(fImpliedVola,c(-0.00001,1),tol=1e-10, Price,OptionType,S,T,K,r)$root)}

GramCharlierPriceIV <- function(Price, K, T)
{
   #if(BSMOption(OptionType,S,T,K,r,0.001,"Price") < Price)
   #{return(0)} 
   #else
   return(uniroot(fImpliedVola,c(-0.001,1),tol=1e-10,Price,OptionType,S,T,K,r)$root)
}

main <- paste("Gram-Charlier-Model - Handelstag ", format(day, "%d-%m-%Y"), sep='')

StrikePrice <- seq(min(data$StrikePrice), max(data$StrikePrice), length.out=51)
t_delta <- seq(min(data$t_delta), max(data$t_delta), length.out=51)
IV <- Price

#IV[4,2] - GramCharlierPrice(StrikePrice[4], t_delta[2])
for(i in 1:dim(Price)[1]){
   for(j in 1:dim(Price)[2])
   {IV[i,j] <- GramCharlierPriceIV(Price[i,j], StrikePrice[i], t_delta[j])}
}

for(i in 1:nrow(data)){
   data[i,"GramCharlierPriceIV"] <- GramCharlierPriceIV(GramCharlierPrice(data$StrikePrice[i], data$t_delta[i]), data$StrikePrice[i], data$t_delta[i])
}
#IV[] <- outer(StrikePrice, t_delta, "GramCharlierPriceIV")

p <- persp(StrikePrice, t_delta, IV , theta=-45, phi=30, col="lightgrey", expand = 0.5, shade = 0.2, xlab="Strike", ylab="Maturity", zlab="BSM Volatility", ticktype = "detailed", main=main)
obs <- trans3d(data$StrikePrice, data$t_delta, data$ImpliedVola, p)
points(obs, col="red",pch=16)
pred <- trans3d(data$StrikePrice, data$t_delta, data$GramCharlierPriceIV, p)
segments(obs$x, obs$y, pred$x, pred$y)



######### Plot All Single days - Prices #########

setwd("D:/h_da/_Blog/Themen/GramCharlier/Plots/Prices")

### Plot Single day

### Price
GramCharlierPrice <- function(K, T, VolatilityPerYear, OnePeriodSkewness, OnePeriodKurtosis)
{GramCharlierOption(OptionType,S,T,K,r,VolatilityPerYear,OnePeriodSkewness,OnePeriodKurtosis,Greek="Price")}

for(i in 1:nrow(GCModel)){
   
   data <- subset(EurexOptionsDax, OptionType == 'Call' & Handelstag == GCModel$Handelstag[i] & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)
   
   # Fit via nlsLM
   OptionType<-data[1,"OptionType"]
   S<-data[1,"SchlusspreisBasiswert"]
   r<-log(1+data[1,"EONIA"])

   filename <- paste("GramCharlier_", format(GCModel$Handelstag[i], "%Y%m%d"), ".png", sep="")
   
   png(filename, width=1200, height=800, pointsize=18)
   
   main <- paste("Gram-Charlier-Model; Handelstag ", format(GCModel$Handelstag[i], "%d-%m-%Y"), sep='')
   subtitle <- paste("Volatility Per Year: ", format(GCModel$VolatilityPerYear[i],digits=4), "; One-PeriodSkewness: ", format(GCModel$OnePeriodSkewness[i],digits=4), "; One-Period-Kurtosis: ", format(GCModel$OnePeriodKurtosis[i],digits=4), sep='')
   
   StrikePrice <- seq(min(data$StrikePrice), max(data$StrikePrice), length.out=51)
   t_delta <- seq(min(data$t_delta), max(data$t_delta), length.out=51)
   Price <- outer(StrikePrice, t_delta, "GramCharlierPrice", GCModel$VolatilityPerYear[i], GCModel$OnePeriodSkewness[i], GCModel$OnePeriodKurtosis[i])
   
   p <- persp(StrikePrice, t_delta, Price , theta=30, phi=30, col="lightgrey", expand = 0.5, shade = 0.2, xlab="Strike", ylab="Maturity", zlab="Price", ticktype = "detailed", main=main, sub=subtitle)
   obs <- trans3d(data$StrikePrice, data$t_delta, data$TaeglicherAbrechnungspreis, p)
   points(obs, col="red",pch=16)
   pred <- trans3d(data$StrikePrice, data$t_delta, GramCharlierPrice(data$StrikePrice, data$t_delta, GCModel$VolatilityPerYear[i], GCModel$OnePeriodSkewness[i], GCModel$OnePeriodKurtosis[i]), p)
   segments(obs$x, obs$y, pred$x, pred$y)
   
   dev.off()
}


######### Plot All Single days - Volas #########

fImpliedVola<-function(sigma,Price,OptionType,S,T,K,r)
   # the function for which the root is sought
{return(BSMOption(OptionType,S,T,K,r,sigma,"Price")-Price)}

ImpliedVola<-function(Price,OptionType,S,T,K,r)
{return(uniroot(fImpliedVola,c(-0.00001,1),tol=1e-10,Price,OptionType,S,T,K,r)$root)}

GramCharlierPriceIV <- function(Price,K,T)
{return(uniroot(fImpliedVola,c(-0.001,1),tol=1e-10,Price,OptionType,S,T,K,r)$root)}

setwd("D:/h_da/_Blog/Themen/GramCharlier/Plots/Volatility")
for(i in 1:nrow(GCModel)){
   
   data <- subset(EurexOptionsDax, OptionType == 'Call' & Handelstag == GCModel$Handelstag[i] & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)
   
   main <- paste("Gram-Charlier-Model - Handelstag ", format(GCModel$Handelstag[i], "%d-%m-%Y"), sep='')
   subtitle <- paste("Volatility Per Year: ", format(GCModel$VolatilityPerYear[i],digits=4), "; One-PeriodSkewness: ", format(GCModel$OnePeriodSkewness[i],digits=4), "; One-Period-Kurtosis: ", format(GCModel$OnePeriodKurtosis[i],digits=4), sep='')
   
   # Fit via nlsLM
   OptionType<-data[1,"OptionType"]
   S<-data[1,"SchlusspreisBasiswert"]
   r<-log(1+data[1,"EONIA"])
   
   filename <- paste("GramCharlier_", format(GCModel$Handelstag[i], "%Y%m%d"), ".png", sep="")
   
   StrikePrice <- seq(min(data$StrikePrice), max(data$StrikePrice), length.out=51)
   t_delta <- seq(min(data$t_delta), max(data$t_delta), length.out=51)
   IV <- outer(StrikePrice, t_delta, "GramCharlierPrice", GCModel$VolatilityPerYear[i], GCModel$OnePeriodSkewness[i], GCModel$OnePeriodKurtosis[i])
   
   #IV[4,2] - GramCharlierPrice(StrikePrice[4], t_delta[2])
   for(j in 1:dim(IV)[1]){
      for(k in 1:dim(IV)[2])
      {IV[j,k] <- GramCharlierPriceIV(IV[j,k], StrikePrice[j], t_delta[k])}
   }
   
   #IV_ <- IV
   
   #j <- 48
   #k <- 4
   
   #GramCharlierPriceIV(IV_[j,k], StrikePrice[j], t_delta[k])

   #GramCharlierPrice(StrikePrice[j], t_delta[k], GCModel$VolatilityPerYear[i], GCModel$OnePeriodSkewness[i], GCModel$OnePeriodKurtosis[i])
   
   #ImpliedVola(IV[j,k],OptionType,S,t_delta[k],StrikePrice[j],r)
   
   #GramCharlierPriceIV(IV[j,k], StrikePrice[j], t_delta[k])
   
   for(j in 1:nrow(data)){
      data[j,"GramCharlierPrice"] <- GramCharlierPrice(data$StrikePrice[j], data$t_delta[j], GCModel$VolatilityPerYear[i], GCModel$OnePeriodSkewness[i], GCModel$OnePeriodKurtosis[i])
   }
   

   for(j in 1:nrow(data)){
      data[j,"GramCharlierPriceIV"] <- GramCharlierPriceIV(GramCharlierPrice(data$StrikePrice[j], data$t_delta[j], GCModel$VolatilityPerYear[i], GCModel$OnePeriodSkewness[i], GCModel$OnePeriodKurtosis[i]), data$StrikePrice[j], data$t_delta[j])
   }
   
   #IV[] <- outer(StrikePrice, t_delta, "GramCharlierPriceIV")
   png(filename, width=1200, height=800, pointsize=18)
   
   par(mfrow = c(1,2))
   
   #mtext(main, side = 3, line = -2, outer = TRUE)
   #text(0.5,0.5,main,cex=2,font=2)
   p <- persp(StrikePrice, t_delta, IV , theta=-45, phi=30, col="lightgrey", expand = 0.5, shade = 0.2, xlab="Strike", ylab="Maturity", zlab="BSM Volatility", ticktype = "detailed")
   obs <- trans3d(data$StrikePrice, data$t_delta, data$ImpliedVola, p)
   points(obs, col="red",pch=16)
   pred <- trans3d(data$StrikePrice, data$t_delta, data$GramCharlierPriceIV, p)
   segments(obs$x, obs$y, pred$x, pred$y)
 
   p <- persp(StrikePrice, t_delta, IV , theta=45, phi=30, col="lightgrey", expand = 0.5, shade = 0.2, xlab="Strike", ylab="Maturity", zlab="BSM Volatility", ticktype = "detailed")
   obs <- trans3d(data$StrikePrice, data$t_delta, data$ImpliedVola, p)
   points(obs, col="red",pch=16)
   pred <- trans3d(data$StrikePrice, data$t_delta, data$GramCharlierPriceIV, p)
   segments(obs$x, obs$y, pred$x, pred$y)
   
   title(main, line = -3, outer=TRUE)
   mtext(side=1, subtitle, line = -4, outer=TRUE)
   
   dev.off()
}


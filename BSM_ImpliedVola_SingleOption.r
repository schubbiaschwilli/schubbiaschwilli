load(".../EurexOptionsDax.RData") 

day <- as.POSIXct("2016-01-04", tz="UTC")

# Select all calls & puts from one day
data <- subset(EurexOptionsDax, Handelstag == day)

Moneyness<-function(OptionType, SchlusspreisBasiswert, StrikePrice)
{
   if(OptionType == 'Call') Value <- SchlusspreisBasiswert/StrikePrice
   if(OptionType == 'Put') Value <- StrikePrice/SchlusspreisBasiswert
   return(Value)
}

# Add Moneyness
data$Moneyness <- mapply(Moneyness,data[,2],data[,4],data[,3])

# Select data with moneyness between 0.5 and 2
data <- subset(data, Moneyness >= 0.5 & Moneyness <= 2)

# Add Maturity and t-delta
Maturity <- function(month, year)
{
   d <- as.Date(ISOdate(year, month, 1, tz="UTC"))
   # wday: 0-6 day of the week, starting on Sunday.
   d <- d + (7 - as.POSIXlt(d)$wday - 2)
   if((as.POSIXlt(d, tz="UTC")$mon+1)!=month) d <- d + 7
   d <- d + 14
   return(d)
}

data$Maturity <- Sys.Date()

for(i in 1:nrow(data))
{data[i,]$Maturity <- Maturity(data[i,]$Verfalltermin_Monat, data[i,]$Verfalltermin_Jahr)}

data$t_delta <- as.numeric(difftime(data$Maturity, data$Handelstag, tz = "UTC", units = "days") / 365)

### Select t delta between 1/52 and 2
data <- subset(data, t_delta >= 1/52 & t_delta <= 2)

### Calc implied vola
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

data$ImpliedVola <- -1

fImpliedVola<-function(sigma,Price,OptionType,S,T,K,r)
# the function for which the root is sought
{return(BSMOption(OptionType,S,T,K,r,sigma,"Price")-Price)}

ImpliedVola<-function(Price,OptionType,S,T,K,r)
{return(uniroot(fImpliedVola,c(-0.01,1),tol=1e-10, Price,OptionType,S,T,K,r)$root)}

for(i in 1:nrow(data))
{
   Price <- data[i,]$TaeglicherAbrechnungspreis
   OptionType <- data[i,]$OptionType
   S <- data[i,]$SchlusspreisBasiswert
   r <- log(1+data[i,]$EONIA)
   T <- data[i,]$t_delta
   K <- data[i,]$StrikePrice

   data[i,]$ImpliedVola <- ImpliedVola(Price,OptionType,S,T,K,r)
}

### Plot smile
library(scatterplot3d);

# Calls
par(mfrow = c(1,2))
data_calls <- subset(data, OptionType == 'Call')
scatterplot3d(data_calls$Moneyness, data_calls$t_delta, data_calls$ImpliedVola, pch=16, highlight.3d=TRUE, angle = 30, xlab="Moneyness", ylab="t", zlab="Implied Vola", grid=T)
scatterplot3d(data_calls$Moneyness, data_calls$t_delta, data_calls$ImpliedVola, pch=16, highlight.3d=TRUE, angle = 315, xlab="Moneyness", ylab="t", zlab="Implied Vola", grid=T)
main <- paste("Implied Vola Calls (Handelstag : ", format(day, "%d-%m-%Y"), ")", sep='')
mtext(main, side = 3, line = -2, outer = TRUE)

# Puts
par(mfrow = c(1,2))
data_puts <- subset(data, OptionType == 'Put')
scatterplot3d(data_puts$Moneyness, data_puts$t_delta, data_puts$ImpliedVola, pch=16, highlight.3d=TRUE, angle = 30, xlab="Moneyness", ylab="t", zlab="Implied Vola", grid=T)
scatterplot3d(data_puts$Moneyness, data_puts$t_delta, data_puts$ImpliedVola, pch=16, highlight.3d=TRUE, angle = 315, xlab="Moneyness", ylab="t", zlab="Implied Vola", grid=T)
main <- paste("Implied Vola Puts (Handelstag : ", format(day, "%d-%m-%Y"), ")", sep='')
mtext(main, side = 3, line = -2, outer = TRUE)

# Boxplot 
boxplot(data$ImpliedVola ~ data$OptionType, xlab="Implied Vola", ylab="Option Type", horizontal=TRUE) 

### Plot Histogramms
library(ggplot2)
ggplot(data, aes(ImpliedVola, fill = OptionType)) + geom_histogram(binwidth = .025, position = "dodge")

load("...\\EurexOptionsDaxPlus.RData")

library("reshape2")
library("ggplot2")

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

# Create df for results
result <- unique(EurexOptionsDax["Handelstag"])
result <- data.frame(result)
colnames(result) <- c("Handelstag")
result$Fitted <- -1
result$ResidualSumOfSquares <- -1
result$DegreesOfFreedom <- -1

# Model Calibration
library(minpack.lm)

for(i in 1:nrow(result))
{
   data <- subset(EurexOptionsDax, OptionType == 'Call' & Handelstag == as.POSIXct(result[i,]$Handelstag, tz="UTC") & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)
   
   # Fit via nlsLM
   Price<-data$TäglicherAbrechnungspreis
   OptionType<-data[1,"OptionType"]
   S<-data[1,"SchlusspreisBasiswert"]
   r<-log(1+data[1,"EONIA"])
   T<-data$t_delta
   K<-data$StrikePrice
   
   startvalue<-list(sigma=0.2);
   IterationControl<-nls.control(maxiter=2^10);
   
   try(result_nlsLM<-nlsLM(Price~BSMOption(OptionType,S,T,K,r,sigma,'Price'), start=startvalue, control=IterationControl), silent=TRUE)
   
   result[i,"Fitted"]<-as.list(coef(result_nlsLM))$sigma;
   result[i,"ResidualSumOfSquares"]<-deviance(result_nlsLM);
   result[i,"DegreesOfFreedom"]<-df.residual(result_nlsLM);
}

# Create statistics
   data_statistcs <- subset(EurexOptionsDax, OptionType == 'Call' & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)
   mean <- aggregate(x=data_statistcs[,c(1,12)], by=list(data_statistcs$Handelstag), FUN=mean)[,c(1,3)]
   sd <- aggregate(x=data_statistcs[,c(1,12)], by=list(data_statistcs$Handelstag), FUN=sd)[,c(1,3)]

   colnames(sd) <- c("Handelstag", "StdDev")
   colnames(mean) <- c("Handelstag", "Mean")

   result <- merge(result, mean, by="Handelstag")
   result <- merge(result, sd, by="Handelstag")
   result$MSE <- result$ResidualSumOfSquares / result$Count

# Create plots
   # Model Calibration vs. Mean
   tmp <- melt(result[,c(1,2,5)], id="Handelstag")  # convert to long format
   tmp$Handelstag <- as.Date(tmp$Handelstag)
   ggplot(data=tmp, aes(x=Handelstag, y=value, colour=variable)) + geom_line(size = 1) + scale_x_date(date_minor_breaks = "1 month") + theme(legend.title=element_blank())

   # RSS vs StdDev
   plot(result[,1],result[,3], type='l', xlim=range(result[,1]), ylim=range(result[,3]), xlab="Handelstag", ylab="")
   par(new = TRUE)
   plot(result[,1],result[,6], type='l', col='red', axes=FALSE, xlab="", ylab="", cex = 1.5)
   axis(4, col.axis='red')
   legend("topright", c("RSS", "Std Dev"), col = c("black","red"), text.font=4, lty=1, cex = 0.75)
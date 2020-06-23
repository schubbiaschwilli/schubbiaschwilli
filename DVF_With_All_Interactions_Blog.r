load(".../EurexOptionsDaxPlus.RData")

library("reshape2")
library("ggplot2")
library("latex2exp")
library("scatterplot3d")

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

DVF_Model_sigma<-function(t_delta, StrikePrice, Model){
   # Intercept
   res <- unname(Model$coefficients[1])
   
   for(i in 2:length(Model$coefficients))
   {
      exp <- gsub(":", "*", names(Model$coefficients)[i])
      
      if(!is.na(unname(Model$coefficients[i]))){
         res <- res + unname(Model$coefficients[i]) * eval(parse(text=exp))
      }   
   }
   return(res)
}

Create_DVF_Model<-function(n){
   if(n>1){
      fStrikePrice <- c("StrikePrice", paste("I(StrikePrice^", 2:n, ")", sep=""))
      ft_delta <- c("t_delta", paste("I(t_delta^", 2:n, ")", sep=""))
   } else {
      fStrikePrice <- c("StrikePrice")
      ft_delta <- c("t_delta")
   }
   
   f <- expand.grid(fStrikePrice, ft_delta)
   formula <- ""
   
   for (i in 1:n){
      formula <- paste(formula, fStrikePrice[i], " + ", sep="")
      formula <- paste(formula, ft_delta[i], " + ", sep="")
   }
   
   for (i in 1:nrow(f)){
      formula <- paste(formula, f[i,1], " * ", f[i,2], sep="")
      if(i != nrow(f)){formula <- paste(formula, " + ", sep="")} 
   }
   
   formula <- paste("ImpliedVola ~ ", formula, sep="")
   return(formula)
}

DVF_Model <- unique(EurexOptionsDax["Handelstag"])
DVF_Model <- data.frame(DVF_Model)
DVF_Model$RSS <- -1
DVF_Model$RMSE_S <- -1

n <- 4
formula <- Create_DVF_Model(n)

for(i in 1:nrow(DVF_Model))
{
      data <- subset(EurexOptionsDax, OptionType == 'Call' & Handelstag == as.POSIXct(DVF_Model[i,]$Handelstag, tz="UTC") & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)
      
      Price<-data$TaeglicherAbrechnungspreis
      OptionType<-data[1,"OptionType"]
      S<-data[1,"SchlusspreisBasiswert"]
      r<-log(1+data[1,"EONIA"])
      T<-data$t_delta
      K<-data$StrikePrice
      
      Model <- lm(as.formula(formula), data=data)
      
      DVF_Model[i,"RSS"] <- sum((Price - BSMOption(OptionType,S,T,K,r,DVF_Model_sigma(T,K,Model),'Price'))^2)
      DVF_Model[i,"RMSE_S"] <- sqrt(sum((Price - BSMOption(OptionType,S,T,K,r,DVF_Model_sigma(T,K,Model),'Price'))^2) / nrow(data)) / S
      
      # Calc predicted RMSE/S 
      if(i>1)
      {
         DVF_Model[i,"RMSE_S_Previous_Vs_Current"] <- sqrt(sum((Price - BSMOption(OptionType,S,T,K,r,DVF_Model_sigma(T,K,previousModel),'Price'))^2) / nrow(data)) / S
      }
      previousModel<-Model
      
      # Plot Surface & Deltas for every single day
      if(i == 248){
         ### Def Plots
         layout(matrix(c(1,1,1,1,1,2,2,2,2),1))
         title<-paste("Handelstag: ", data$Handelstag[i], sep="")
         subtitle<-paste("RSS:", format(round(DVF_Model$RSS[i], digits=4), big.mark=","), sep="")
         
         ### Volatility Surface
         # Data 
         Strike<-seq(min(data$StrikePrice), max(data$StrikePrice), length=51)
         t<-seq(min(data$t_delta), max(data$t_delta), length=51)
         Price <- Strike %o% t
         
         for(j in 1:length(Strike)){
            for(k in 1:length(t)){
               Price[j,k] <- BSMOption(OptionType,S,t[k],Strike[j],r,DVF_Model_sigma(t[k],Strike[j],Model),'Price')
            }
         }
         
         # Plot Price: Scatter vs Surface/ persp(x, y, z,...)
         plot<-persp(Strike, t, Price, theta=45, phi=15, expand=0.5, col='gray', main = NULL, ticktype="detailed", sub=NULL);
         # Add real prices
         ObservedPrices <- trans3d(data$StrikePrice, data$t_delta, data$TaeglicherAbrechnungspreis, plot);
         points(ObservedPrices, col="red", pch=16);
         
         CalculatetPrices_ <- rep(0,length(data$StrikePrice))
         
         for(j in 1:length(CalculatetPrices_)){
            CalculatetPrices_[j] <- BSMOption(OptionType,S,data$t_delta[j],data$StrikePrice[j],r,DVF_Model_sigma(data$t_delta[j],data$StrikePrice[j],Model),'Price')
         }
         
         CalculatetPrices <- trans3d(data$StrikePrice, data$t_delta, CalculatetPrices_, plot);
         segments(ObservedPrices$x, ObservedPrices$y, CalculatetPrices$x, CalculatetPrices$y)
         
         ### Deltas 
         deltas<-data$TaeglicherAbrechnungspreis-CalculatetPrices_
         scatterplot3d(data$StrikePrice, data$t_delta, deltas, pch=16, highlight.3d=TRUE, angle=45, type="h", box=FALSE, main=NULL, xlab="Strike", ylab="t", zlab="Delta DVF- vs Observed Price")  
         mtext(title, side=3, line=-2, outer=TRUE, cex=1)
         mtext(subtitle, side=1, line=-2, outer=TRUE, cex=1)
      }
}

# Plot RMSE_S
tmp <- melt(DVF_Model[, c(1,3,4)], id="Handelstag")
tmp$Handelstag <- as.Date(tmp$Handelstag)
ggplot(data=tmp, aes(x=Handelstag, y=value, color=factor(variable, labels=c("RMSE/S","RMSE/S Steady State")))) + geom_line(size = 1) + ggtitle("RMSE/S") + scale_x_date(date_minor_breaks = "1 month") + theme_bw(base_size=16) + theme(legend.position="bottom") + theme(legend.title=element_blank())

# Plot RSS
tmp <- melt(DVF_Model[, c(1,2)], id="Handelstag")
tmp$Handelstag <- as.Date(tmp$Handelstag)
ggplot(data=tmp, aes(x=Handelstag, y=value, color=factor(variable, labels=c("RSS")))) + geom_line(size = 1) + ggtitle("RSS") + scale_x_date(date_minor_breaks = "1 month") + theme_bw(base_size=16) + theme(legend.position="bottom") + theme(legend.title=element_blank())


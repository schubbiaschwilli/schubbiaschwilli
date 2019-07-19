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

DVF_Model1_sigma<-function(K, a_0, a_1, a_2)
{  # Model 1 (a_0 + a_1 K + a_2 K^2)
   return(a_0 + a_1*K + a_2*(K^2))
}

DVF_Model2_sigma<-function(T, K, a_0, a_1, a_2, a_3, a_5)
{  # Model 2 (a_0 + a_1 K + a_2 K^2 + a_3 T + a_5 K*T)
   return(a_0 + a_1*K + a_2*(K^2) + a_3*T + a_5*K*T)
}

DVF_Model3_sigma<-function(T, K, a_0, a_1, a_2, a_3, a_4, a_5)
{  # Model 3 (a_0 + a_1 K + a_2 K^2 + a_3 T + a_4 T^2 + a_5 K*T)
   return(a_0 + a_1*K + a_2*(K^2) + a_3*T + a_4*(T^2) + a_5*K*T)
}

# Create df for results
DVF_Model0 <- unique(EurexOptionsDax["Handelstag"])
DVF_Model0 <- data.frame(DVF_Model0)
DVF_Model0$RMSE_S <- -1
DVF_Model0$a_0 <- 0

DVF_Model1 <- DVF_Model0
DVF_Model1$a_1 <- 0
DVF_Model1$a_2 <- 0

DVF_Model2 <-DVF_Model1
DVF_Model2$a_3 <- 0
DVF_Model2$a_5 <- 0

DVF_Model3 <- DVF_Model2
DVF_Model3$a_4 <- 0

### Calc DVF
for(i in 1:nrow(DVF_Model0))
{
   data <- subset(EurexOptionsDax, OptionType == 'Call' & Handelstag == as.POSIXct(DVF_Model3[i,]$Handelstag, tz="UTC") & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)
   
   Price<-data$TaeglicherAbrechnungspreis
   OptionType<-data[1,"OptionType"]
   S<-data[1,"SchlusspreisBasiswert"]
   r<-log(1+data[1,"EONIA"])
   T<-data$t_delta
   K<-data$StrikePrice
   
   # Model 0 (Constant)
   Model0 <- lm(ImpliedVola ~ 1, data = data)
   DVF_Model0[i,"a_0"]<-as.list(coef(Model0))[1]
   
   DVF_Model0[i,"RMSE_S"] <- sqrt(sum((Price - BSMOption(OptionType,S,T,K,r,DVF_Model0[i,"a_0"],'Price'))^2) / nrow(data)) / S
   
   # Model 1 (a_0 + a_1 K + a_2 K^2)
   Model1 <- lm(ImpliedVola  ~ StrikePrice + I(StrikePrice^2), data = data)
   DVF_Model1[i,"a_0"]<-as.list(coef(Model1))[1]
   DVF_Model1[i,"a_1"]<-as.list(coef(Model1))[2]
   DVF_Model1[i,"a_2"]<-as.list(coef(Model1))[3]
   
   DVF_Model1[i,"RMSE_S"] <- sqrt(sum((Price - BSMOption(OptionType,S,T,K,r,DVF_Model1_sigma(K,DVF_Model1[i,"a_0"],DVF_Model1[i,"a_1"],DVF_Model1[i,"a_2"]),'Price'))^2) / nrow(data)) / S
   
   # Model 2 (a_0 + a_1 K + a_2 K^2 + a_3 T + a_5 K*T)
   Model2 <- lm(ImpliedVola  ~ StrikePrice + I(StrikePrice^2) + t_delta + StrikePrice:t_delta, data = data)
   DVF_Model2[i,"a_0"]<-as.list(coef(Model2))[1]
   DVF_Model2[i,"a_1"]<-as.list(coef(Model2))[2]
   DVF_Model2[i,"a_2"]<-as.list(coef(Model2))[3]
   DVF_Model2[i,"a_3"]<-as.list(coef(Model2))[4]
   DVF_Model2[i,"a_5"]<-as.list(coef(Model2))[5]
   
   DVF_Model2[i,"RMSE_S"] <- sqrt(sum((Price - BSMOption(OptionType,S,T,K,r,DVF_Model2_sigma(T,K,DVF_Model2[i,"a_0"],DVF_Model2[i,"a_1"],DVF_Model2[i,"a_2"],DVF_Model2[i,"a_3"],DVF_Model2[i,"a_5"]),'Price'))^2) / nrow(data))  / S
   
   # Model 3 (a_0 + a_1 K + a_2 K^2 + a_3 T + a_4 T^2 + a_5 K*T)
   Model3 <- lm(ImpliedVola  ~ StrikePrice + I(StrikePrice^2) + t_delta + I(t_delta^2) + StrikePrice:t_delta, data = data)
   DVF_Model3[i,"a_0"]<-as.list(coef(Model3))[1]
   DVF_Model3[i,"a_1"]<-as.list(coef(Model3))[2]
   DVF_Model3[i,"a_2"]<-as.list(coef(Model3))[3]
   DVF_Model3[i,"a_3"]<-as.list(coef(Model3))[4]
   DVF_Model3[i,"a_4"]<-as.list(coef(Model3))[5]
   DVF_Model3[i,"a_5"]<-as.list(coef(Model3))[6]
   
   DVF_Model3[i,"RMSE_S"] <- sqrt(sum((Price - BSMOption(OptionType,S,T,K,r,DVF_Model3_sigma(T,K,DVF_Model3[i,"a_0"],DVF_Model3[i,"a_1"],DVF_Model3[i,"a_2"],DVF_Model3[i,"a_3"],DVF_Model3[i,"a_4"],DVF_Model3[i,"a_5"]),'Price'))^2) / nrow(data)) / S
   
   # Calc predicted RMSE/S 
   if(i>1)
   {
      DVF_Model0[i,"RMSE_S_Previous_Vs_Current"] <- sqrt(sum((Price - BSMOption(OptionType,S,T,K,r,DVF_Model0[i-1,"a_0"],'Price'))^2) / nrow(data)) / S
      DVF_Model1[i,"RMSE_S_Previous_Vs_Current"] <- sqrt(sum((Price - BSMOption(OptionType,S,T,K,r,DVF_Model1_sigma(K,DVF_Model1[i-1,"a_0"],DVF_Model1[i-1,"a_1"],DVF_Model1[i-1,"a_2"]),'Price'))^2) / nrow(data)) / S
      DVF_Model2[i,"RMSE_S_Previous_Vs_Current"] <- sqrt(sum((Price - BSMOption(OptionType,S,T,K,r,DVF_Model2_sigma(T,K,DVF_Model2[i-1,"a_0"],DVF_Model2[i-1,"a_1"],DVF_Model2[i-1,"a_2"],DVF_Model2[i-1,"a_3"],DVF_Model2[i-1,"a_5"]),'Price'))^2) / nrow(data)) / S
      DVF_Model3[i,"RMSE_S_Previous_Vs_Current"] <- sqrt(sum((Price - BSMOption(OptionType,S,T,K,r,DVF_Model3_sigma(T,K,DVF_Model3[i-1,"a_0"],DVF_Model3[i-1,"a_1"],DVF_Model3[i-1,"a_2"],DVF_Model3[i-1,"a_3"],DVF_Model3[i-1,"a_4"],DVF_Model3[i-1,"a_5"]),'Price'))^2) / nrow(data)) / S
   }
}

### Create plots
library("reshape2")
library("ggplot2")
library("latex2exp")
library("scatterplot3d")

# Model 0: RMSE_S v. RMSE_S VolaFittedPrevious
tmp <- melt(DVF_Model0[,c(1,2,4)], id="Handelstag")
tmp$Handelstag <- as.Date(tmp$Handelstag)

ggplot(data=tmp, aes(x=Handelstag, y=value, colour=factor(variable, labels = c("RMSE / S","RMSE / S with predicted parameters from previous market data")))) +
   geom_line(size = 1) + ggtitle(TeX("Model 0 : $\\sigma_{iv} = a_0$")) + scale_x_date(date_minor_breaks = "1 month") + 
   theme_bw(base_size=16) + theme(legend.position="bottom") + theme(legend.title=element_blank())

# Model 1: RMSE_S v. RMSE_S VolaFittedPrevious
tmp <- melt(DVF_Model1[,c(1,2,6)], id="Handelstag")
tmp$Handelstag <- as.Date(tmp$Handelstag)

ggplot(data=tmp, aes(x=Handelstag, y=value, colour=factor(variable, labels = c("RMSE / S","RMSE / S with predicted parameters from previous market data")))) +
   geom_line(size = 1) + ggtitle(TeX("Model 1 : $\\sigma_{iv} = a_0 + a_1 K + a_2 K^2$")) + scale_x_date(date_minor_breaks = "1 month") +
   theme_bw(base_size=16) + theme(legend.position="bottom") + theme(legend.title=element_blank())

# Model 2: RMSE_S v. RMSE_S VolaFittedPrevious
tmp <- melt(DVF_Model2[,c(1,2,8)], id="Handelstag")
tmp$Handelstag <- as.Date(tmp$Handelstag)

ggplot(data=tmp, aes(x=Handelstag, y=value, colour=factor(variable, labels = c("RMSE / S","RMSE / S with predicted parameters from previous market data")))) +
   geom_line(size = 1) + ggtitle(TeX("Model 2 : $\\sigma_{iv} = a_0 + a_1 K + a_2 K^2 + a_3 T + a_5 KT$")) + scale_x_date(date_minor_breaks = "1 month") +
   theme_bw(base_size=16) + theme(legend.position="bottom") + theme(legend.title=element_blank())

# Model 3 RMSE_S v. RMSE_S VolaFittedPrevious
tmp <- melt(DVF_Model3[,c(1,2,9)], id="Handelstag")
tmp$Handelstag <- as.Date(tmp$Handelstag)

ggplot(data=tmp, aes(x=Handelstag, y=value, colour=factor(variable, labels = c("RMSE / S","RMSE / S with predicted parameters from previous market data")))) +
   geom_line(size = 1) + ggtitle(TeX("Model 3 : $\\sigma_{iv} = a_0 + a_1 K + a_2 K^2 + a_3 T + a_4 T^2 + a_5 KT$")) + scale_x_date(date_minor_breaks = "1 month") +
   theme_bw(base_size=16) + theme(legend.position="bottom") + theme(legend.title=element_blank())

### Plot Single day
PlotRegression <- function(model, DVF, formula)
{
   main <- paste("DFV Model ", DVF, "; ", data$OptionType[1] , "s Handelstag ", format(day, "%d-%m-%Y"), sep='')
   
   StrikePrice <- seq(min(data$StrikePrice), max(data$StrikePrice), length.out=51)
   t_delta <- seq(min(data$t_delta), max(data$t_delta), length.out=51)
   IV <- outer(StrikePrice, t_delta, function(StrikePrice, t_delta) predict(model, data.frame(StrikePrice, t_delta)))
   
   p <- persp(StrikePrice, t_delta, IV, theta=30, phi=30, col="lightgrey", expand = 0.5, shade = 0.2, xlab="Strike", ylab="Maturity", zlab="Volatility", ticktype = "detailed", main=main, sub=TeX(formula))
   obs <- trans3d(data$StrikePrice, data$t_delta, data$ImpliedVola, p)
   pred <- trans3d(data$StrikePrice, data$t_delta, fitted(model), p)
   
   points(obs, col="red",pch=16)
   segments(obs$x, obs$y, pred$x, pred$y)
}

fff <- function(x){
   # Formula For Function
   if(x<0) {s="-"}
   else {s="+"}
   
   return(paste(s, "\\;", format(abs(x),digits=4), sep=""))
}

day <- as.POSIXct("2016-04-01", tz="UTC")
data <- subset(EurexOptionsDax, OptionType == 'Call' & Handelstag == day & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)

# Model 1 (a_0 + a_1 K + a_2 K^2)
model1 <- lm(ImpliedVola  ~ StrikePrice + I(StrikePrice^2), data = data)
formula<-paste("$\\sigma_{iv} \\; =\\; ", fff(model1$coefficients[1]), "\\;" ,fff(model1$coefficients[2]), "\\; K \\;", fff(model1$coefficients[3]), "\\; K^2$", sep="")
PlotRegression(model1, 1, formula)

# Model 2 (a_0 + a_1 K + a_2 K^2 + a_3 T + a_5 K*T)
model2 <- lm(ImpliedVola  ~ StrikePrice + I(StrikePrice^2) + t_delta + StrikePrice*t_delta, data = data)
formula<-paste("$\\sigma_{iv} \\; = \\; ", fff(model2$coefficients[1]), "\\;" , fff(model2$coefficients[2]), "\\; K \\;", fff(model2$coefficients[3]), "\\; K^2 \\; ", fff(model2$coefficients[4]), "\\; T \\;", fff(model2$coefficients[5]), "\\; KT$", sep="")
PlotRegression(model2, 2, formula)

# Model 3 (a_0 + a_1 K + a_2 K^2 + a_3 T + a_4 T^2 + a_5 K*T)
model3 <- lm(ImpliedVola  ~ StrikePrice + I(StrikePrice^2) + t_delta + I(t_delta^2) + StrikePrice*t_delta, data = data)
formula<-paste("$\\sigma_{iv}\\; =\\; ", fff(model3$coefficients[1]), "\\;", fff(model3$coefficients[2]), "\\; K \\;", fff(model3$coefficients[3]), "\\; K^2", fff(model3$coefficients[4]), "\\; T\\;", fff(model3$coefficients[5]), "\\; T^2 \\;", fff(model3$coefficients[6]), "\\; KT \\; $", sep="")
PlotRegression(model3, 3, formula)

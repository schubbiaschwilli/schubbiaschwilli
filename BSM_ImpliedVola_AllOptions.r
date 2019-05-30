load(".../EurexOptionsDaxPlus.RData")

   # Create statistics
      data_statistcs <- subset(EurexOptionsDax, OptionType == 'Call' & Moneyness >= 0.8 & Moneyness <= 1.2 & t_delta >= 1/12 & t_delta <= 1)
      
      max <- aggregate(x=data_statistcs[,c(1,12)], by=list(data_statistcs$Handelstag), FUN=max)[,c(1,3)]
      mean <- aggregate(x=data_statistcs[,c(1,12)], by=list(data_statistcs$Handelstag), FUN=mean)[,c(1,3)]
      min <- aggregate(x=data_statistcs[,c(1,12)], by=list(data_statistcs$Handelstag), FUN=min)[,c(1,3)]
      sd <- aggregate(x=data_statistcs[,c(1,12)], by=list(data_statistcs$Handelstag), FUN=sd)[,c(1,3)]
      
      colnames(max) <- c("Handelstag", "Maximum")
      colnames(sd) <- c("Handelstag", "StdDev")
      colnames(min) <- c("Handelstag", "Minimum")
      colnames(mean) <- c("Handelstag", "Mean")
      
      statistics <- data.frame(unique(data_statistcs$Handelstag))
      colnames(statistics) <- c("Handelstag")
      
      statistics <- merge(statistics, max, by="Handelstag")
      statistics <- merge(statistics, mean, by="Handelstag")
      statistics <- merge(statistics, min, by="Handelstag")
      statistics <- merge(statistics, sd, by="Handelstag")
      statistics$MeanPlusStddDev <- statistics$Mean + statistics$StdDev
      statistics$MeanMinusStddDev <- statistics$Mean - statistics$StdDev
      
      statistics <- statistics[,c(1,2,6,3,7,4)]
      colnames(statistics) = c("Handelstag", "Maximum", "Mean Plus StddDev", "Mean", "Mean Minus StddDev", "Minimum")
     
   # Plot TS
      library(ggplot2)   
      library(reshape)
      
      statistics.melted <- melt(statistics, variable_name="StatisticalSize", id = "Handelstag") 
      ggplot(data=statistics.melted, aes(x=Handelstag, y=value, color=StatisticalSize)) + geom_line(size=.75) + labs(x="Handelstag", y="", col="")
   
   # Plot Std.Dev.
      ggplot(data=sd, aes(x=Handelstag, y=StdDev)) + geom_line(color = "Blue", size = 1)

   # Get example data wito vola <> 0
      subset(EurexOptionsDax, OptionType == 'Call' & StrikePrice >= 8250 & StrikePrice <= 8550 & Verfalltermin_Monat == 4 & Verfalltermin_Jahr == 2016 & Handelstag == ISOdate(2016, 3, 17, 0, 0, tz = "UTC"))

   # Plot DAX
      DAX <- aggregate(x=EurexOptionsDax[,c(1,4)], by=list(EurexOptionsDax$Handelstag), FUN=max)[,-1]
      plot(DAX, type='l')

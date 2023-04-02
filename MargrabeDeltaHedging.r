S1 <- 100
S2 <- 100
mu1 <- 0.1
mu2 <- 0.1
sigma1 <- 0.25
sigma2 <- 0.25
rho <- -0.3
T <- 0.5

NumberOfPaths <- 10^4
NumberOfHedgepoints <- 52*T

dt <- T / NumberOfHedgepoints

MargrabeOption<-function(S1, S2, sigma1, sigma2, rho, T, Greek){
  sigma <- sqrt(sigma1^2 + sigma2^2 - 2*rho*sigma1*sigma2)
  d1 <- (log(S1/S2) + T * (sigma^2 / 2)) / (sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  if(Greek=="Price" && T>0) return(S1*pnorm(d1, mean=0, sd=1) - S2*pnorm(d2, mean=0, sd=1))
  if(Greek=="Price" && T==0) return(pmax(S2-S1, 0))
  if(Greek=="Delta1" && T>0) return(pnorm(d1, mean=0, sd=1))
  if(Greek=="Delta2" && T>0) return(-pnorm(d2, mean=0, sd=1))
  if((Greek=="Delta1" || Greek=="Delta2") && T==0) return(0)
}

results <- data.frame(S1=rep(NA, times=NumberOfPaths),
                      S2=rep(NA, times=NumberOfPaths),
                      Hedge=rep(NA, times=NumberOfPaths))

set.seed(sum(strtoi(charToRaw("MARGRABE"),16L)))

for(j in 1:NumberOfPaths){
  S1_temp <- S1
  S2_temp <- S2
  
  # new stock position S1 and S2
  alpha1 <- MargrabeOption(S1_temp, S2_temp, sigma1, sigma2, rho, T, "Delta1")
  alpha2 <- MargrabeOption(S1_temp, S2_temp, sigma1, sigma2, rho, T, "Delta2")
            
  # initial investment
  V_T <- MargrabeOption(S1_temp, S2_temp, sigma1, sigma2, rho, T, "Price")
  V_T <- V_T - (S1_temp*alpha1 + S2_temp*alpha2)
  
  for(i in 1:NumberOfHedgepoints){
    # Simulate N(0,1)
    Wt1 <- rnorm(1)	
    Wt2 <- Wt1 * rho + sqrt(1-(rho^2))*rnorm(1)
    
    S1_temp <- S1_temp*exp((mu1-sigma1^2/2)*dt+sigma1*sqrt(dt)*Wt1)
    S2_temp <- S2_temp*exp((mu2-sigma2^2/2)*dt+sigma2*sqrt(dt)*Wt2)
    
    # Vermögenszprozess - Kein Cashbond, also _OHNE_ Verzinsung !!!
    V_T <- (S1_temp*alpha1 + S2_temp*alpha2)
    
    # new stock position S1 and S2
    alpha1 <- MargrabeOption(S1_temp, S2_temp, sigma1, sigma2, rho, (T - i*dt), "Delta1")
    alpha2 <- MargrabeOption(S1_temp, S2_temp, sigma1, sigma2, rho, (T - i*dt), "Delta2")

    # New alpha - HedgeError re-adjustments
    V_T <- V_T - (S1_temp * alpha1 + S2_temp * alpha2)
    alpha1 <- alpha1 + V_T / 2 / S1_temp
    alpha2 <- alpha2 + V_T / 2 / S2_temp
  }
  
  V_T <- (S1_temp * alpha1 + S2_temp * alpha2)
  
  results$S1[j] <- S1_temp
  results$S2[j] <- S2_temp
  results$Hedge[j] <- V_T
}

results$PayOff <- pmax(results$S1-results$S2, 0)
results$Error <- results$PayOff - results$Hedge

### Plot Asset 1 vs 2
  exercised <- subset(results, S1>S2)
  notexercised <- subset(results, !S1>S2)

  plot(exercised$S1, exercised$S2, pch=20, col="green", xlab="Asset 1", ylab="Asset 2", xlim=range(results$S1), ylim=range(results$S2))
  points(notexercised$S1, notexercised$S2, pch=20, col="darkgrey")

### Plot Asset 1 vs 2 with Colour depending on Error
  rbPal <- colorRampPalette(c('red','green'))
  results$ErrorCol <- rbPal(10)[as.numeric(cut(results$Error, breaks=10))]
  plot(results$S1, results$S2, pch=20, cex=0.75, col=results$ErrorCol, xlab="Asset 1", ylab="Asset 2", xlim=range(results$S1), ylim=range(results$S2))
  
### Plot Errors
  par(mfrow=c(1,2))
  hist(results$Error, breaks=nclass.FD(results$Error), freq=F, col="darkgrey", border="black", main='', xlab='Error Delta Hedging', ylab='') 
  # Add NormalDistribution
  stdabw<-function(x) {n=length(x) ; sqrt(var(x) * n / (n-1))}
  x<-seq(floor(min(results$Error)), ceiling(max(results$Error)), by=0.05)
  points(x, dnorm(x, mean=mean(results$Error), sd=stdabw(results$Error)), type="l", lwd=2, col=2)
  
  ### QQ-Plot Error
  qqnorm(results$Error, main='')
  qqline(results$Error, lwd=2, col=2)
  
### Rotating plot/gif
  # https://rpubs.com/aagarwal29/179912
  PayOff <- function(S1, S2) return(pmax(S1-S2, 0))
  Asset1<-seq(min(results$S1), max(results$S1), length=51)
  Asset2<-seq(min(results$S2), max(results$S2), length=51)
  PayOff <- outer(Asset1, Asset2, PayOff)
  
  library(rgl)
  open3d(windowRect=c(200, 200, 800, 800))
  plot3d(results$S1, results$S2, results$Hedge, type="s", size=0.5, lit=TRUE, main="", sub="", xlab="Asset1", ylab="Asset2", zlab="Payoff/Hedge", col="red")
  surface3d(Asset1, Asset2, PayOff, alpha=0.5, front="lines", back="lines")
  spin3d(axis=c(0, 0, 1), rpm=30)
  # Create gif
  # movie3d(spin3d(axis=c(0, 0, 1), rpm=30), duration=2, movie="MargrabeSpin", clean=FALSE)
library(quantmod)
t = 1/252
From <- 2009
To <- 2018

PDF_MJD <- function(x, mu, sigma, lambda, m, delta){
   # Rama/Cont p. 111; gamma <-> mu, mu <-> m
   sum <- 0
   
   for(k in 0:170)
   {
      numerator <- (lambda*t)^k * exp(-((x - mu*t - k*m)^2) / (2*(sigma^2*t + k*delta^2)))
      denominator <- factorial(k) * sqrt(2*pi * (sigma^2*t + k*delta^2))
      sum = sum + numerator/denominator
   }
   return(sum*exp(-lambda*t))
}

CDF_MJD<-function(x, mu, sigma, lambda, m, delta)
{integrate(PDF_MJD, lower=-Inf, upper=x, mu=mu, sigma=sigma, lambda=lambda, m=m, delta=delta)$value}

Likelihood_PDF_MJD <- function(theta, x)
{ -sum(log(PDF_MJD(x, mu=theta[1], sigma=theta[2], lambda=theta[3], m=theta[4], delta=theta[5])))}

# Startvalues
mu <- 0.2689
sigma <- 0.0755
lambda <- 287.61
m <- -0.00066
delta <- 0.01107

# Get data
getSymbols("^GDAXI", index.class="POSIXct", from=as.Date(ISOdate(From-1, 12, 30)), to=as.Date(ISOdate(To, 12, 31)), src="yahoo")

   # Calc Log-Returns & na.omit(),
   log_returns <- log(na.omit(GDAXI$GDAXI.Close))
   log_returns <- diff(log_returns, lag=1)[-1]
   log_returns <- as.data.frame(log_returns$GDAXI.Close)
   log_returns <- as.matrix(log_returns$GDAXI.Close)

# Fit
results <- nlm(f=Likelihood_PDF_MJD, p=c(mu, sigma, lambda, m, delta), x=log_returns, iterlim=250)

#### Plot probability density function ###
stdabw <- function(x) {n=length(x) ; sqrt(var(x) * n / (n-1))}

main<-paste(From, " - ", To, sep="")
sub<-paste("mu: ", round(results$estimate[1],5) , "; sigma: ", round(results$estimate[2],5) , "; lambda: ", round(results$estimate[3],5) , "; m: ", round(results$estimate[4],5) , "; delta: ", round(results$estimate[5],5), sep="")
hist(log_returns, breaks = nclass.FD(log_returns), freq=FALSE, col="gray", xlab=sub, main=main) # main="Histogramm der stetigen Renditen vs Normalverteilung"
curve(dnorm(x, mean=mean(log_returns), sd=stdabw(log_returns)), col='blue', type="l", lwd=2, lty=2, add=TRUE)
curve(PDF_MJD(x, mu, sigma, lambda, m, delta), col='red', type="l", lwd=2, lty=1, add=TRUE)

#### Plot cumulative distribution function ###
# Data for CDF
NumberOfSupportPoints <- 1001
x <- seq(from=-0.075, to=.075, length.out=NumberOfSupportPoints)
y <- rep(0, NumberOfSupportPoints)

for(j in 1:NumberOfSupportPoints)
{y[j] <- CDF_MJD(x[j], mu, sigma, lambda, m, delta)}

# Plot
main<-paste(From, " - ", To, sep="")
sub<-paste("mu: ", round(results$estimate[1],5) , "; sigma: ", round(results$estimate[2],5), "; lambda: ", round(results$estimate[3],5) , "; m: ", round(results$estimate[4],5) , "; delta: ", round(results$estimate[5],5), sep="")

plot(ecdf(log_returns), lty=1, pch='+', lwd=1, ylab="CDF", col='blue', xlab=sub, main=main)
lines(x,y, col='red')
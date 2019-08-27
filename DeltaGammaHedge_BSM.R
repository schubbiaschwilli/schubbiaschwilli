BSMOption <- function(PutCallFlag,S,T,K,r,sigma,Greek)
# 'Four Things You Might Not Know About The Black-Scholes-Formula'
{
   if(PutCallFlag=="Call" && Greek=="Price" && T>0)
   {
      d1=(log(S/K)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
      d2=d1 - sigma * sqrt(T)
      return(S * pnorm(d1, mean=0, sd=1) - K * exp(-r*T) * pnorm(d2, mean=0, sd=1)) 
   }
   if(PutCallFlag=="Call" && Greek=="Delta" && T>0)
   {
      d1=(log(S/K)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
      return(pnorm(d1, mean=0, sd=1)) 
   }
   if((PutCallFlag=="Call" || PutCallFlag=="Put") && Greek=="Gamma" && T>0)
   {
      d1=(log(S/K)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
      return(dnorm(d1, mean=0, sd=1) / (S * sigma * sqrt(T))) 
   }
   if(PutCallFlag=="Call" && Greek=="Delta" && T==0)
      {return(0)}
   if(PutCallFlag=="Call" && Greek=="Price" && T==0)
      {return(pmax(S-K,0))} 
   if(PutCallFlag=="Put" && Greek=="Price" && T==0)
      {return(pmax(K-S,0))}
   if(PutCallFlag=="Put")
      {return(-1 * BSMOption("Call",S,T,K,r,-1 * sigma,Greek))}
}

varianz <- function(x) {n=length(x) ; var(x) * n / (n-1)}
stdabw <- function(x) {n=length(x) ; sqrt(var(x) * n / (n-1))}

# Parameter
PutCallFlag = "Call"; # Put or Call
S_0 = 100;
T = 0.5;
K = 100;
r = 0.05;
mu = 0.15;
sigma = 0.5;
Npaths = 1000 # 5 * 10^3;
Nhedgepoints = 52;

##################
set.seed(123456);
dt = T/Nhedgepoints;

# Data
S_t<-rep(0, times=Npaths);
Payoff<-rep(0, times=Npaths);

# DeltaHedge
Hedge_DeltaHedge<-rep(0, times=Npaths);
Error_DeltaHedge<-rep(0, times=Npaths);

# DeltaGammaHedge
Hedge_DeltaGammaHedge<-rep(0, times=Npaths);
Error_DeltaGammaHedge<-rep(0, times=Npaths);

# Vergleiche 'Four Things You Might Not Know About The Black-Scholes-Formula'

# Loop over j = 1 to Npaths
for(j in 1:Npaths) 
{
   S <- S_0;
   
   # DeltaHedge
   V_DeltaHedge <- BSMOption(PutCallFlag,S,T,K,r,sigma,"Price"); # initial investment
   a_DeltaHedge <- BSMOption(PutCallFlag,S,T,K,r,sigma,"Delta"); # stock position = delta
   b_DeltaHedge <- V_DeltaHedge-a_DeltaHedge*S; # rest in bank; self-fin. Cond.
   
   # DeltaGammaHedge
   K0 <- S
   # initial investment
   V_DeltaGammaHedge <- BSMOption(PutCallFlag,S,T,K,r,sigma,"Price") 
   #c=BScall(S,T-i*dt,K1,.,greek=gamma)/BScall(S,T-i*dt,K0,.,greek=gamma)
   c_DeltaGammaHedge <- BSMOption(PutCallFlag,S,T,K,r,sigma,"Gamma")/BSMOption(PutCallFlag,S,T,K0,r,sigma,"Gamma")
   #a=BScall(S,T-i*dt,K1,.,greek=delta)-c*BScall(S,T-i*dt,K0,.,greek=delta)
   a_DeltaGammaHedge <- BSMOption(PutCallFlag,S,T,K,r,sigma,"Delta")-c_DeltaGammaHedge*BSMOption(PutCallFlag,S,T,K0,r,sigma,"Delta")
   #b=V-a*S-c*BScall(S,T-i*dt,K0,., greek=price)
   b_DeltaGammaHedge <- V_DeltaGammaHedge-a_DeltaGammaHedge*S-c_DeltaGammaHedge * BSMOption(PutCallFlag,S,T,K0,r,sigma,"Price")
   
   # Loop over i=1 to Nhedgepoints
   for(i in 1:Nhedgepoints) 
   {
      eps <- rnorm(1); # eps=random.normal() simulate outcome of N(0,1)
      S <- S*exp((mu-0.5*sigma^2)*dt+sigma*sqrt(dt)*eps);
      
      # Delta-Hedge
      V_DeltaHedge<-a_DeltaHedge*S+b_DeltaHedge*exp(r*dt);
      a_DeltaHedge<-BSMOption(PutCallFlag,S,T-i*dt,K,r,sigma,"Delta");
      b_DeltaHedge<-V_DeltaHedge-a_DeltaHedge*S;
      
      # DeltaGammaHedging : K0 = K_liquid (S_T) ; K1 = illiquid (S0)
      
      # V=a*S+b*exp(r*dt)+c*BScall(S,T-i*dt, K0,.,greek=price)
      V_DeltaGammaHedge<-a_DeltaGammaHedge*S+b_DeltaGammaHedge*exp(r*dt)+c_DeltaGammaHedge*BSMOption(PutCallFlag,S,T-i*dt,K0,r,sigma,"Price")
      
      K0 <- S;
      
      #c=BScall(S,T-i*dt,K1,.,greek=gamma)/BScall(S,T-i*dt,K0,.,greek=gamma)
      if (i==Nhedgepoints)
         c_DeltaGammaHedge <- 0 else
            c_DeltaGammaHedge <- BSMOption(PutCallFlag,S,T-i*dt,K,r,sigma,"Gamma")/BSMOption(PutCallFlag,S,T-i*dt,K0,r,sigma,"Gamma")
      
      #a=BScall(S,T-i*dt,K1,.,greek=delta)-c*BScall(S,T-i*dt,K0,.,greek=delta)
      a_DeltaGammaHedge <- BSMOption(PutCallFlag,S,T-i*dt,K,r,sigma,"Delta")-c_DeltaGammaHedge*BSMOption(PutCallFlag,S,T-i*dt,K0,r,sigma,"Delta")
      
      #b=V-a*S-c* BScall(S,T-i*dt,K0,., greek=price)
      b_DeltaGammaHedge <- V_DeltaGammaHedge-a_DeltaGammaHedge*S-c_DeltaGammaHedge*BSMOption(PutCallFlag,S,T-i*dt,K0,r,sigma,"Price")
   }
   
   S_t[j] <- S;
   Hedge_DeltaHedge[j] <- V_DeltaHedge;
   Hedge_DeltaGammaHedge[j] <- V_DeltaGammaHedge;
   Payoff[j] <- BSMOption(PutCallFlag,S,0,K,r,sigma,"Price");
}

# Calc Error
Error_DeltaHedge<-Payoff - Hedge_DeltaHedge;
Error_DeltaGammaHedge<-Payoff - Hedge_DeltaGammaHedge;

### Plots

# Plot PayOff vs Hedge
   plot(S_t, Hedge_DeltaHedge, col = "red", type = "p", pch='.', cex=3, main = "Results Hedging" , ylab="Hedge / PayOff", xlab="")
   points(S_t, Hedge_DeltaGammaHedge, col = 'blue', type = "p", pch='.', cex=3)
   points(S_t, Payoff, col = 'green', type = "p", pch='.', cex=3)

   # Legende formatieren
   legend(if(PutCallFlag=="Call"){min(S_t)} else max(S_t)*3/4, max(Payoff), legend = c("Payoff","DeltaHedge","DeltaGammaHedge"), col = c("green", "red","blue"), pch=c(19,19,19), cex=1)

# Define Histogramms Errors and Boxplot Errors
   layout(rbind(c(1,2),c(1,2),c(3,4)))

	# Plot Hist from Error
		hist(Error_DeltaHedge, breaks = nclass.FD(Error_DeltaHedge), freq=FALSE, col = "gray", border="black", main = '' , xlab='Error DeltaHeding',ylab='relative frequencies') 
		# Add NormalDistribution
		x<-seq(floor(min(Error_DeltaHedge)), ceiling(max(Error_DeltaHedge)), by=0.05) 
		points(x, dnorm(x, mean=mean(Error_DeltaHedge), sd=stdabw(Error_DeltaHedge)), type = "l", lwd=2, col=2)

		hist(Error_DeltaGammaHedge, breaks = nclass.FD(Error_DeltaGammaHedge), freq=FALSE, col = "gray", border="black", main = '' , xlab='Error DeltaGammaHeding',ylab='relative frequencies') 
		# Add NormalDistribution
		x<-seq(floor(min(Error_DeltaGammaHedge)), ceiling(max(Error_DeltaGammaHedge)), by=0.05)
		points(x, dnorm(x, mean=mean(Error_DeltaGammaHedge), sd=stdabw(Error_DeltaGammaHedge)), type = "l", lwd=2, col=2)

   # Plot Boxplot
	   boxplot(Error_DeltaHedge, horizontal=TRUE, main='BoxPlot Error DeltaHedging')
	   boxplot(Error_DeltaGammaHedge, horizontal=TRUE, main='BoxPlot Error DeltaGammaHedging')

# Define QQ-Plots
  layout(rbind(c(1,2)))

	# QQ-Plot der Fehler
		qqnorm(Error_DeltaHedge, main='QQPlot Error DeltaHedging')
		qqline(Error_DeltaHedge, lwd=2 , col=2)

	# QQ-Plot der Fehler
		qqnorm(Error_DeltaGammaHedge, main='QQPlot Error DeltaGammaHedging')
		qqline(Error_DeltaGammaHedge, lwd=2 , col=2)

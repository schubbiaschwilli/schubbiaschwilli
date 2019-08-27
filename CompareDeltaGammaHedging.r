BSMOption<-function(PutCallFlag,S,T,K,r,sigma,Greek)
  # 'Four Things You Might Not Know About The Black-Scholes-Formula'
{
   if(PutCallFlag=="Call" && Greek=="Price" && T>0)
   {
      d1=(log(S/K)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
      d2=d1 - sigma * sqrt(T)
      return(S * pnorm(d1, mean=0, sd=1) - K * exp(-r*T) * pnorm(d2, mean=0, sd=1)) 
   }
   if(PutCallFlag=="Call" && Greek=="Price" && T==0)
      {return(pmax(S-K,0))} 
   if(PutCallFlag=="Call" && Greek=="Delta" && T>0)
   {
      d1=(log(S/K)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
      return(pnorm(d1, mean=0, sd=1))
   }
   if(PutCallFlag=="Call" && Greek=="Delta" && T==0)
      {return(0)}
   if(Greek=="Gamma" && T > 0)
   {      
      d1 <- (log(S/K)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
      value <- dnorm(d1, mean=0, sd=1) / (S * sigma * sqrt(T))
      return(value) 
   }
   else   
   if(PutCallFlag=="Put" && T>0)
      {return(-1 * BSMOption("Call",S,T,K,r,-1 * sigma,Greek))}
   if(PutCallFlag=="Put" && Greek=="Price" && T==0)
      {return(pmax(K-S,0))}
}

### Parameter ###
PutCallFlag = "Call"; # Put or Call
S_0 = 100;
T = 1;
K = S_0 * 1.15;
r = 0.05;
mu = 0.1;
sigma = 0.2
Npaths = 10000

Nhedgepoints_List<-c(seq(10, 100, 5), seq(100, 1000, 50))

# Create df for results
ResultsDeltaHedges <- data.frame(Nhedgepoints_List)
names(ResultsDeltaHedges) <- c("Nhedgepoints")
ResultsDeltaHedges$Stddev_Error <- 0
ResultsDeltaHedges$Variance_Error <- 0

ResultsDeltaGammaHedges <- ResultsDeltaHedges

varianz <- function(x) {n=length(x) ; var(x) * n / (n-1)}
stdabw<-function(x) {n=length(x) ; sqrt(var(x) * n / (n-1))}

### Simulation Engine ###
for(k in 1:length(Nhedgepoints_List)) 
{
  Nhedgepoints<-Nhedgepoints_List[k]
  set.seed(123456);
  dt = T/Nhedgepoints;
  
  # Data
  S_t<-rep(0, times=Npaths);
  Payoff<-rep(0, times=Npaths);
  Hedge_DeltaHedge<-rep(0, times=Npaths);
  Error_DeltaHedge<-rep(0, times=Npaths);

  Hedge_DeltaGammaHedge<-rep(0, times=Npaths);
  Error_DeltaGammaHedge<-rep(0, times=Npaths);
  
  # Loop over j = 1 to Npaths
  for(j in 1:Npaths)
  {
    S <- S_0;
    
    # DeltaHedge
    V_DeltaHedge<-BSMOption(PutCallFlag,S,T,K,r,sigma,"Price"); # initial investment
    a_DeltaHedge<-BSMOption(PutCallFlag,S,T,K,r,sigma,"Delta"); # stock position = delta
    b_DeltaHedge<-V_DeltaHedge-a_DeltaHedge*S; # rest in bank; self-fin. Cond.
    
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
      eps<-rnorm(1); # eps=random.normal() simulate outcome of N(0,1)
      S<-S*exp((mu-0.5*sigma^2)*dt+sigma*sqrt(dt)*eps);
      
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
    
    S_t[j]<-S;
    Hedge_DeltaHedge[j]<-V_DeltaHedge;
    Hedge_DeltaGammaHedge[j]<-V_DeltaGammaHedge;
    Payoff[j]<-BSMOption(PutCallFlag,S,0,K,r,sigma,"Price");
  }
  
  # Calc Error
  Error_DeltaHedge<-Payoff - Hedge_DeltaHedge;
  
  ResultsDeltaHedges$Stddev_Error[k] <- stdabw(Error_DeltaHedge)
  ResultsDeltaHedges$Variance_Error[k] <- varianz(Error_DeltaHedge)

  Error_DeltaGammaHedge<-Payoff - Hedge_DeltaGammaHedge;
  
  ResultsDeltaGammaHedges$Stddev_Error[k] <- stdabw(Error_DeltaGammaHedge)
  ResultsDeltaGammaHedges$Variance_Error[k] <- varianz(Error_DeltaGammaHedge)
}

library(minpack.lm)

### Plot Hedge-Error vs #Hedgepoints Std.Dev.
   ### Fit Delta-Hedge ###
   fitDeltaHedge <- nlsLM(Stddev_Error ~ a*Nhedgepoints^b, data=ResultsDeltaHedges, start=list(a=7, b=-0.5) , control=nls.control(maxiter=2^10))
   legendDeltaHedge <- paste("Delta-Hedging: y = ", round(as.list(coef(fitDeltaHedge))$a,4), " x^", round(as.list(coef(fitDeltaHedge))$b,4), sep="")

   plot(ResultsDeltaHedges$Nhedgepoints, ResultsDeltaHedges$Stddev_Error, ylab="Std.Dev.(Error)", xlab="# hedge points", ylim=c(0,6))
   lines((10:1000), as.list(coef(fitDeltaHedge))$a*(10:1000)^as.list(coef(fitDeltaHedge))$b, col = "red")

   ### Fit Delta-Gamma-Hedge ###
   fitDeltaGammaHedge <- nlsLM(Stddev_Error ~ a*Nhedgepoints^b, data=ResultsDeltaGammaHedges, start=list(a=4, b=-0.7) , control=nls.control(maxiter=2^10))
   legendDeltaGammaHedge <- paste("Delta-Gamma-Hedging: y = ", round(as.list(coef(fitDeltaGammaHedge))$a,4), " x^", round(as.list(coef(fitDeltaGammaHedge))$b,4), sep="")

   points(ResultsDeltaGammaHedges$Nhedgepoints, ResultsDeltaGammaHedges$Stddev_Error)
   lines((10:1000), as.list(coef(fitDeltaGammaHedge))$a*(10:1000)^as.list(coef(fitDeltaGammaHedge))$b, col = "blue")

   legend("topright", legend=c(legendDeltaHedge, legendDeltaGammaHedge), col = c("red", "blue"), lty=1, cex=1)
   
### Plot Hedge-Error vs #Hedgepoints Variance
   ### Fit Delta-Hedge ###
   fitDeltaHedge <- nlsLM(Variance_Error ~ a*Nhedgepoints^b, data=ResultsDeltaHedges, start=list(a=31, b=-1) , control=nls.control(maxiter=2^10))
   legendDeltaHedge <- paste("Delta-Hedging: y = ", round(as.list(coef(fitDeltaHedge))$a,4), " x^", round(as.list(coef(fitDeltaHedge))$b,4), sep="")
   
   plot(ResultsDeltaHedges$Nhedgepoints, ResultsDeltaHedges$Variance_Error, ylab="Variance(Error)", xlab="# hedge points", ylim=c(0,6))
   lines((10:1000), as.list(coef(fitDeltaHedge))$a*(10:1000)^as.list(coef(fitDeltaHedge))$b, col = "red")
   
   ### Fit Delta-Gamma-Hedge ###
   fitDeltaGammaHedge <- nlsLM(Variance_Error ~ a*Nhedgepoints^b, data=ResultsDeltaGammaHedges, start=list(a=4, b=-1) , control=nls.control(maxiter=2^10))
   legendDeltaGammaHedge <- paste("Delta-Gamma-Hedging: y = ", round(as.list(coef(fitDeltaGammaHedge))$a,4), " x^", round(as.list(coef(fitDeltaGammaHedge))$b,4), sep="")
   
   points(ResultsDeltaGammaHedges$Nhedgepoints, ResultsDeltaGammaHedges$Variance_Error)
   lines((10:1000), as.list(coef(fitDeltaGammaHedge))$a*(10:1000)^as.list(coef(fitDeltaGammaHedge))$b, col = "blue")
   
   legend("topright", legend=c(legendDeltaHedge, legendDeltaGammaHedge), col = c("red", "blue"), lty=1, cex=1)
   
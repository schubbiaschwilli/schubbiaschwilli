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
   if(PutCallFlag=="Call" && Greek=="Delta" && T==0)
      {return(0)}
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
Npaths = 1000

Nhedgepoints_List<-c(seq(10, 100, 10), seq(100, 1000, 50))

Error_DeltaHedges_Nhedgepoints<-rep(0, length(Nhedgepoints_List));
Error_DeltaHedges_stdabw<-rep(0, length(Nhedgepoints_List));
Error_DeltaHedges_varianz<-rep(0, length(Nhedgepoints_List));

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
   
   # Loop over j = 1 to Npaths
   for(j in 1:Npaths) 
   {
      S = S_0;
      
      # DeltaHedge
      V_DeltaHedge<-BSMOption(PutCallFlag,S,T,K,r,sigma,"Price"); # initial investment
      a_DeltaHedge<-BSMOption(PutCallFlag,S,T,K,r,sigma,"Delta"); # stock position = delta
      b_DeltaHedge<-V_DeltaHedge-a_DeltaHedge*S; # rest in bank; self-fin. Cond.
      
      K0 = S
      
      # Loop over i=1 to Nhedgepoints
      for(i in 1:Nhedgepoints) 
      {
         eps<-rnorm(1); # eps=random.normal() simulate outcome of N(0,1)
         S<-S*exp((mu-0.5*sigma^2)*dt+sigma*sqrt(dt)*eps);
         
         # Delta-Hedge
         V_DeltaHedge<-a_DeltaHedge*S+b_DeltaHedge*exp(r*dt);
         a_DeltaHedge<-BSMOption(PutCallFlag,S,T-i*dt,K,r,sigma,"Delta");
         b_DeltaHedge<-V_DeltaHedge-a_DeltaHedge*S;
   
         K0<-S;
      }
      
      S_t[j]<-S;
      Hedge_DeltaHedge[j]<-V_DeltaHedge;
      Payoff[j]<-BSMOption(PutCallFlag,S,0,K,r,sigma,"Price");
   }
   
   # Calc Error
   Error_DeltaHedge<-Payoff - Hedge_DeltaHedge;
   
   Error_DeltaHedges_Nhedgepoints[k]<-Nhedgepoints_List[k] 
   Error_DeltaHedges_stdabw[k]<-stdabw(Error_DeltaHedge)
   Error_DeltaHedges_varianz[k]<-varianz(Error_DeltaHedge)
}
### Plots ###
library(minpack.lm)

fit <- nlsLM(Error_DeltaHedges_stdabw ~ a*Error_DeltaHedges_Nhedgepoints^b, start=list(a=7.3787, b=-0.498) , control=nls.control(maxiter=2^10))

a<-as.list(coef(fit))$a
b<-as.list(coef(fit))$b             

plot(Error_DeltaHedges_Nhedgepoints, Error_DeltaHedges_stdabw, ylab="std. dev. Error Delta Hedges", xlab="# hedge points")
lines((10:1000), a*(10:1000)^b, col = "red")
legend("topright",legend=c(paste("y = ", round(a,4), " x^", round(b,4), sep="")), col = c("red"), lty=1, cex=1)

###
fit <- nlsLM(Error_DeltaHedges_varianz ~ a*Error_DeltaHedges_Nhedgepoints^b, start=list(a=53.54, b=-0.983) , control=nls.control(maxiter=2^10))

a<-as.list(coef(fit))$a
b<-as.list(coef(fit))$b             

plot(Error_DeltaHedges_Nhedgepoints, Error_DeltaHedges_varianz, ylab="Variance Error Delta Hedges", xlab="# hedge points")
lines((10:1000), a*(10:1000)^b, col = "red")
legend("topright",legend=c(paste("y = ", round(a,4), " x^", round(b,4), sep="")), col = c("red"), lty=1, cex=1)

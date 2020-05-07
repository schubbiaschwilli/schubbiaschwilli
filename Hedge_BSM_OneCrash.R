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

# Hedge
Nhedgepoints<-52
Npaths<-25000
dt = T/Nhedgepoints;

# Parameter
PutCallFlag <- "Call"; # Put or Call
S_0 <- 100
K <- 100
T <- 1
r <- 0.05
mu <- 0.1
sigma <- 0.2

# Jump
m <- -0.25

# Create Jumps & Dataframe
set.seed(234567)
Jumps <- data.frame(sample(1:Nhedgepoints, size=Npaths, replace=T), rep(NA, Npaths))
colnames(Jumps) <- c("Interval", "S")
   
# Data
S_t<-rep(0, times=Npaths);
Payoff<-rep(0, times=Npaths);
   
# DeltaHedge
Hedge_DeltaHedge<-rep(0, times=Npaths);
Error_DeltaHedge<-rep(0, times=Npaths);
   
# DeltaGammaHedge
Hedge_DeltaGammaHedge<-rep(0, times=Npaths);
Error_DeltaGammaHedge<-rep(0, times=Npaths);
   
# Init
set.seed(123456)
   
# Loop over j = 1 to Npaths
for(j in 1:Npaths) 
{
   S = S_0;
   
   # DeltaHedge
   V_DeltaHedge=BSMOption(PutCallFlag,S,T,K,r,sigma,"Price"); # initial investment
   a_DeltaHedge=BSMOption(PutCallFlag,S,T,K,r,sigma,"Delta"); # stock position = delta
   b_DeltaHedge=V_DeltaHedge-a_DeltaHedge*S; # rest in bank; self-fin. Cond.
   
   # DeltaGammaHedge
   K0 = S
   # initial investment
   V_DeltaGammaHedge=BSMOption(PutCallFlag,S,T,K,r,sigma,"Price") 
   #c=BScall(S,T-i*dt,K1,.,greek=gamma)/BScall(S,T-i*dt,K0,.,greek=gamma)
   c_DeltaGammaHedge=BSMOption(PutCallFlag,S,T,K,r,sigma,"Gamma")/BSMOption(PutCallFlag,S,T,K0,r,sigma,"Gamma")
   #a=BScall(S,T-i*dt,K1,.,greek=delta)-c*BScall(S,T-i*dt,K0,.,greek=delta)
   a_DeltaGammaHedge=BSMOption(PutCallFlag,S,T,K,r,sigma,"Delta")-c_DeltaGammaHedge*BSMOption(PutCallFlag,S,T,K0,r,sigma,"Delta")
   #b=V-a*S-c*BScall(S,T-i*dt,K0,., greek=price)
   b_DeltaGammaHedge=V_DeltaGammaHedge-a_DeltaGammaHedge*S-c_DeltaGammaHedge * BSMOption(PutCallFlag,S,T,K0,r,sigma,"Price")
   
   # Loop over i=1 to Nhedgepoints
   for(i in 1:Nhedgepoints){
      # Jumps (or Diffusion)
      if(Jumps[j, "Interval"]==i){
            Jumps[j, "S"] <- S
            S<-S*(1+m)
         } else {
            S<-S*exp((mu-0.5*sigma^2)*dt+sigma*sqrt(dt)*rnorm(1)) # Diffusion
      }
         
      # Delta-Hedge
      V_DeltaHedge=a_DeltaHedge*S+b_DeltaHedge*exp(r*dt);
      a_DeltaHedge=BSMOption(PutCallFlag,S,T-i*dt,K,r,sigma,"Delta");
      b_DeltaHedge=V_DeltaHedge-a_DeltaHedge*S;
      
      # DeltaGammaHedging : K0 = K_liquid (S_T) ; K1 = illiquid (S0)
      
      # V=a*S+b*exp(r*dt)+c*BScall(S,T-i*dt, K0,.,greek=price)
      V_DeltaGammaHedge=a_DeltaGammaHedge*S+b_DeltaGammaHedge*exp(r*dt)+c_DeltaGammaHedge*BSMOption(PutCallFlag,S,T-i*dt,K0,r,sigma,"Price")
      
      K0=S
      
      #c=BScall(S,T-i*dt,K1,.,greek=gamma)/BScall(S,T-i*dt,K0,.,greek=gamma)
      if (i==Nhedgepoints)
         c_DeltaGammaHedge=0 else
            c_DeltaGammaHedge=BSMOption(PutCallFlag,S,T-i*dt,K,r,sigma,"Gamma")/BSMOption(PutCallFlag,S,T-i*dt,K0,r,sigma,"Gamma")
      
      #a=BScall(S,T-i*dt,K1,.,greek=delta)-c*BScall(S,T-i*dt,K0,.,greek=delta)
      a_DeltaGammaHedge=BSMOption(PutCallFlag,S,T-i*dt,K,r,sigma,"Delta")-c_DeltaGammaHedge*BSMOption(PutCallFlag,S,T-i*dt,K0,r,sigma,"Delta")
      
      #b=V-a*S-c* BScall(S,T-i*dt,K0,., greek=price)
      b_DeltaGammaHedge=V_DeltaGammaHedge-a_DeltaGammaHedge*S-c_DeltaGammaHedge*BSMOption(PutCallFlag,S,T-i*dt,K0,r,sigma,"Price")
   }
   
   S_t[j] = S;
   Hedge_DeltaHedge[j] = V_DeltaHedge;
   Hedge_DeltaGammaHedge[j] = V_DeltaGammaHedge;
   Payoff[j] = BSMOption(PutCallFlag,S,0,K,r,sigma,"Price");
}

# Calc Error
Error_DeltaHedge<- Hedge_DeltaHedge - Payoff
Error_DeltaGammaHedge<- Hedge_DeltaGammaHedge - Payoff 
delta <- Hedge_DeltaGammaHedge - Hedge_DeltaHedge

# Plot PayOff vs Hedge
   par(mfrow = c(2, 2),oma = c(0,0,2,0))
   plot(S_t, Hedge_DeltaHedge, col = "red", type = "p", pch=20, cex=1, ylab="Hedge Error / PayOff", xlab="S(t)")
   points(S_t, Hedge_DeltaGammaHedge, col = 'blue', type = "p", pch=20, cex=1)
   points(S_t, Payoff, col = 'green', type = "p", pch=20, cex=1)
   
   # Histogramms
   hist(delta, breaks=nclass.FD(delta), freq=T, col="gray", main='' , xlab='Delta-Gamma- minus Delta-Hedge') 
   hist(Error_DeltaHedge, breaks=nclass.FD(Error_DeltaHedge), freq=T, col="gray", main='' , xlab='Error DeltaHedging') 
   hist(Error_DeltaGammaHedge, breaks=nclass.FD(Error_DeltaGammaHedge), freq=T, col="gray", main='' , xlab='Error DeltaGammaHedging') 
   mtext(paste("Strike ", K, sep=""), line=-2, outer=TRUE, cex=1.5)

   par(mfrow = c(2, 2),oma = c(0,0,2,0))
   plot(S_t, Hedge_DeltaHedge, col = "red", type = "p", pch=20, cex=1, ylab="Hedge Error / PayOff", xlab="S(t)")
   points(S_t, Hedge_DeltaGammaHedge, col = 'blue', type = "p", pch=20, cex=1)
   points(S_t, Payoff, col = 'green', type = "p", pch=20, cex=1)
      
   # Define Plots Histogramms Errors and Boxplot Errors
   hist(delta, breaks=nclass.FD(delta), freq=T, col="gray", main='' , xlab='Delta-Gamma- minus Delta-Hedge') 
   hist(Error_DeltaHedge, breaks=nclass.FD(Error_DeltaHedge), freq=T, col="gray", main='' , xlab='Error DeltaHedging') 
   hist(Error_DeltaGammaHedge, breaks=nclass.FD(Error_DeltaGammaHedge), freq=T, col="gray", main='' , xlab='Error DeltaGammaHedging') 
   
   mtext(paste("Strike ", K, sep=""), line=-2, outer=TRUE, cex=1.5)

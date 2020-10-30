pchisq_ding <- function(x, F, theta){
  chi2nc <- x
  
  if(x <= 0){
    return(NA)  
  } 

  lam <- theta / 2
  n <- 1
  u <- exp(-lam)
  v <- u
  X2 <- x / 2
  f2 <- F / 2
  
  tmp1 <- f2*log(X2) - X2 - gammln(f2 + 1)
  
  if(tmp1 > 700){
    return(chilarge(x, F, theta))
  }
  
  t <- exp(tmp1)
  
  if(t <= 1e-15){
    return(chilarge(x, F, theta)) 
  } 
  
  term <- v*t
  chi2nc <- term
                          
  repeat{
    if(F + 2*n - x > 0){
      bound <- t*x / (F + 2*n - x)
      if(bound < 1e-15){
        return(chi2nc)
      }
    }
    u <- u*lam / n
    v <- v + u
    t <- t*x / (F + 2*n)
    term <- v*t
    chi2nc <- chi2nc + term
    n <- n + 1    
  }
}

chilarge <- function(z2, v2, k2){
  z <- 0.5*z2
  v <- 0.5*v2
  k <- 0.5*k2
  
  h <- 1 - 2*(2*v + 2*k)*(2*v + 6*k) / ((2*v + 4*k)*(2*v + 4*k)*3)
  p <- (2*v + 4*k) / ((2*v + 2*k)*(2*v + 2*k))
  m <- (h - 1)*(1 - 3*h)
  A <- 1 - h*p*(1 - h + 0.5*(2 - h)*m*p) - (2*z / (2*v + 2*k)) ^ h
  b <- h*sqrt(2*p*(1 + m*p))
  rr <- A / b
  
  x <- pnorm(rr)
  return(1 - x)
}

gammln <- function(xx){
  # Define constants
  c1 <- 76.1800917294715
  c2 <- -86.5053203294168
  c3 <- 24.0140982408309
  c4 <- -1.23173957245015
  c5 <- 1.20865097386618E-03
  c6 <- -5.395239384953E-06
    
  if(xx <= 0){
    return(0)  
  }
  
  Y <- xx
  x <- xx
  tmp <- x + 5.5
  tmp <- tmp - (x + 0.5)*log(tmp)
  ser <- 1.00000000019001
  
  Y <- Y + 1
  ser <- ser + c1 / Y
  Y <- Y + 1
  ser <- ser + c2 / Y
  Y <- Y + 1
  ser <- ser + c3 / Y
  Y <- Y + 1
  ser <- ser + c4 / Y
  Y <- Y + 1
  ser <- ser + c5 / Y
  Y <- Y + 1
  ser <- ser + c6 / Y
  
  return(-tmp + log(2.506628274631*ser / x))
}

CEVOption<-function(PutCallFlag,S,T,K,r,sigma,alpha){

  if(alpha==1){
    return(BSMOption(PutCallFlag,S,T,K,r,sigma,Greek="Price"))
  }
  
  if(alpha!=1){
    q <- 0
    ny <- (sigma^2)/(2*(r-q)*(alpha-1))*(exp(2*(r-q)*(alpha-1)*T)-1)
    
    a <- (K*exp(-(r-q)*T))^(2*(1-alpha))/(((1-alpha)^2)*ny)
    b <- 1/(1-alpha)
    c <- (S^(2*(1-alpha)))/((1-alpha)^2*ny)
  }
  
  if(alpha>0 && alpha<1 && PutCallFlag=="Call"){   
    return(S*exp(-q*T)*(1-pchisq_ding(a, b+2, theta=c)) - K*exp(-r*T)*pchisq_ding(c, b, theta=a))
  }
  
  if(alpha>0 && alpha<1 && PutCallFlag=="Put"){
    return(K*exp(-r*T)*(1-pchisq_ding(c, b, theta=a)) - S*exp(-q*T)*pchisq_ding(a, b+2, theta=c))
  }
  
  if(alpha>1 && PutCallFlag=="Call"){
    return(S*exp(-q*T)*(1-pchisq_ding(c, -b, theta=a)) - K*exp(-r*T)*pchisq_ding(a, 2-b, theta=c))
  }
  
  if(alpha>1 && PutCallFlag=="Put"){
    return(K*exp(-r*T)*(1-pchisq_ding(a, 2-b, theta=c)) - S*exp(-q*T)*pchisq_ding(c, -b, theta=a))
  }
}

BSMOption<-function(PutCallFlag,S,T,K,r,sigma,Greek="Price"){
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

# Load data
data <- read.csv(".../CEVTestData.csv", stringsAsFactors=FALSE)

data[,9] <- as.POSIXct(data[,9])
data[,10] <- as.POSIXct(data[,10])

for(i in 1:nrow(data)){  
  if(data[i,8] == -1){
    data[i,9] = Sys.time()
    data[i,8] = CEVOption(PutCallFlag=data[i,1],S=data[i,3],T=data[i,2],K=data[i,4],r=data[i,5],sigma=data[i,7],alpha=data[i,6])
    data[i,10] = Sys.time()
  }
}

save(data, file=".../CEVTestData_R_With_Ding.rdata")

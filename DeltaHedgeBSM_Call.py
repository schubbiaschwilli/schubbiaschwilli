import math 
import numpy as np
import scipy.stats as si

# Vergleiche 'Four Things You Might Not Know About The Black-Scholes-Formula'
def BSMOption(PutCallFlag, S, T, K, r, sigma, Greek='Price'):
         
    if PutCallFlag == 'Call' and Greek == 'Price' and T > 0:
        d1 = (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
        d2 = (np.log(S / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
        result = (S * si.norm.cdf(d1, 0.0, 1.0) - K * np.exp(-r * T) * si.norm.cdf(d2, 0.0, 1.0))
    if PutCallFlag == 'Call' and Greek == 'Price' and T == 0:
        result = max(S-K,0)
    if PutCallFlag == 'Call' and Greek == 'Delta' and T > 0:
        result =  si.norm.cdf((np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T)), 0.0, 1.0) 
    if PutCallFlag == 'Call' and Greek == 'Delta' and T == 0:
        result = 0
    if PutCallFlag == 'Put':
        result = -BSMOption('Call', S, T, K, r, -sigma, Greek=Greek)
    if PutCallFlag == 'Put' and Greek == 'Price' and T == 0:
        result = max(K-S,0)
    return result

### Parameter ###
PutCallFlag = "Call" # Put or Call
S_0 = 100.0
T = 1.0
K = S_0 * 1.15
r = 0.05
mu = 0.1
sigma = 0.2
Nhedgepoints = 52
Npaths = 1000

dt = T/Nhedgepoints

# Results
S_t = [0.0] * Npaths
Payoff = [0.0] * Npaths
Hedge = [0.0] * Npaths
Error = [0.0] * Npaths

np.random.seed(12345)

### Simulation Engine ###
# Loop over Npaths
for j in range(0, Npaths):

   S = S_0
   
   # DeltaHedge
   V = BSMOption(PutCallFlag,S,T,K,r,sigma,"Price") # initial investment
   a = BSMOption(PutCallFlag,S,T,K,r,sigma,"Delta") # stock position<-delta
   b = V-a*S # rest in bank, self-fin. Cond.
   
   # Loop over Nhedgepoints
   for i in range(0, Nhedgepoints+1):
      eps = np.random.normal(0, 1) # eps=random.normal() simulate outcome of N(0,1)
      S = S*math.exp((mu-0.5*sigma ** 2)*dt+sigma*math.sqrt(dt)*eps)
      
      # Delta-Hedge
      V = a*S+b*math.exp(r*dt)
      a = BSMOption(PutCallFlag,S,T-i*dt,K,r,sigma,"Delta")
      b = V-a*S
   
   S_t[j] = S
   Hedge[j] = V
   Payoff[j] = BSMOption(PutCallFlag,S,0,K,r,sigma,"Price")

# Calc Error
Error = np.array(Payoff) - np.array(Hedge)

# Plot
import matplotlib.pyplot as plt
plt.scatter(S_t, Hedge, s=1, color='red', label='Hedge')
plt.scatter(S_t, Payoff, s=1, color='green', label='Pay Off')
lgnd = plt.legend(loc='upper left', fontsize = 'x-large')
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
plt.title('Results Delta Hedging')
plt.xlabel('S(t)')
plt.ylabel('Pay off / Hedge')
plt.show()

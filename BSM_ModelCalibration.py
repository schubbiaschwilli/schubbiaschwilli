import numpy as np
import scipy.stats as si
import pandas as pd
from scipy.optimize import least_squares

### BSM-Functions 
def BSMOption(PutCallFlag, S, T, K, r, sigma, Greek='Price'):
    if PutCallFlag == 'Call' and Greek == 'Price' and T > 0:
        d1 = (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
        d2= d1 - sigma * np.sqrt(T)
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

def BSMCostFunc(sigma, PutCallFlag, S, T, K, r, Price):
    sum = 0.0
    for i in range(len(Price)):
        # lambd, vbar, eta, rho, v0, r, q, tau, S0, K
        sum = sum + (BSMOption(PutCallFlag[i], S[i], T[i], K[i], r[i], sigma, Greek='Price') - Price[i])**2
    return sum

### Load data
rawdata = pd.read_hdf('...\EurexOptionsDaxPlus.h5')

### Create Dataframe for Results
result = pd.DataFrame(rawdata['Handelstag'].unique(), columns=['Handelstag'])
result['ImpliedVolatility'] = -1.0
result['RSS'] = -1.0

### Fit Implied Vola
for i in range(len(result)):   
    data = rawdata.loc[(rawdata['Handelstag'] == result['Handelstag'][i]) & (rawdata['OptionType'] == "Call") & (rawdata['Moneyness'] >= 0.8) & (rawdata['Moneyness'] <= 1.2) & (rawdata['t_delta'] >= 1/12) & (rawdata['t_delta'] <= 1)]
    data = data.reset_index()
    tmp = least_squares(fun=BSMCostFunc, x0=0.2, args=(data["OptionType"], data["SchlusspreisBasiswert"], data["t_delta"], data["StrikePrice"], np.log(1+data["EONIA"]), data["TaeglicherAbrechnungspreis"]), method='lm')
    result['ImpliedVolatility'][i] = tmp.x[0]
    result['RSS'][i] = tmp.fun[0]

### create statistic and merge data 
statistic_data = rawdata.loc[(rawdata['OptionType'] == "Call") & (rawdata['Moneyness'] >= 0.8) & (rawdata['Moneyness'] <= 1.2) & (rawdata['t_delta'] >= 1/12) & (rawdata['t_delta'] <= 1)]
statistics = statistic_data.groupby('Handelstag')['ImpliedVola'].agg(['max', 'mean', 'min', 'std', 'count'])
statistics['Handelstag'] = statistics.index
statistics['ImpliedVolaPlusStdDev'] = statistics['mean'] + statistics['std']
statistics['ImpliedVolaMinusStdDev'] = statistics['mean'] - statistics['std']

result = result.merge(statistics, left_on='Handelstag', right_on='Handelstag')
result['MSE'] = result['RSS'] / result['count'] 

### Plots
import matplotlib.pyplot as plt

### statistics
plt.plot('Handelstag', 'max', data=statistics, linewidth=0.5, label="Max")
plt.plot('Handelstag', 'ImpliedVolaPlusStdDev', data=statistics, linewidth=0.5, label="Mean + Std Dev")
plt.plot('Handelstag', 'mean', data=statistics, linewidth=0.5, label="Mean")
plt.plot('Handelstag', 'ImpliedVolaMinusStdDev', data=statistics, linewidth=0.5, label="Mean - Std Dev")
plt.plot('Handelstag', 'min', data=statistics, linewidth=0.5, label="Min")
plt.legend(loc=1, fontsize='xx-small')

### Implied volatilities
plt.plot('Handelstag', 'ImpliedVolatility', data=result, linewidth=0.5, label="Fit")
plt.plot('Handelstag', 'mean', data=result, linewidth=0.5, label="Mean")
plt.legend()

### RSS vs Std Dev
fig, ax1 = plt.subplots()
ax1.set_xlabel('Handelstag')
ax1.set_ylabel('RSS', color="blue")
ax1.plot('Handelstag', 'RSS', data=result, linewidth=0.5, label="RSS", color="blue")
ax2 = ax1.twinx()
ax2.set_ylabel('Std Dev', color="orange")
ax2.plot('Handelstag', 'std', data=result, linewidth=0.5, label="Std Dev", color="orange")
fig.tight_layout() 

### MSE vs Std Dev
fig, ax1 = plt.subplots()
ax1.set_xlabel('Handelstag')
ax1.set_ylabel('MSE', color="blue")
ax1.plot('Handelstag', 'MSE', data=result, linewidth=0.5, label="MSE", color="blue")
ax2 = ax1.twinx()
ax2.set_ylabel('Std Dev', color="orange")
ax2.plot('Handelstag', 'std', data=result, linewidth=0.5, label="Std Dev", color="orange")
fig.tight_layout()
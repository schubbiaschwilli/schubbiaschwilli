import pandas as pd
import numpy as np
import scipy.stats as si

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

vBSMOption = np.vectorize(BSMOption)

def DVF_Model2_sigma(T, K, a_0, a_1, a_2, a_3, a_5):
    # Model 2 (a_0 + a_1 K + a_2 K^2 + a_3 T + a_5 K*T)
    return a_0 + a_1*K + a_2*(K**2) + a_3*T + a_5*K*T

vDVF_Model2_sigma = np.vectorize(DVF_Model2_sigma)

def DVF_Model3_sigma(T, K, a_0, a_1, a_2, a_3, a_4, a_5):
    # Model 3 (a_0 + a_1 K + a_2 K^2 + a_3 T + a_4 T^2 + a_5 K*T)
    return a_0 + a_1*K + a_2*(K**2) + a_3*T + a_4*(T**2) + a_5*K*T

vDVF_Model3_sigma = np.vectorize(DVF_Model3_sigma)

# to open or create a HDF5 file
rawdata = pd.read_hdf('...\EurexOptionsDaxPlus.h5')

# Create dataframe for results
DVF_Model2_result = pd.DataFrame(rawdata['Handelstag'].unique(), columns=['Handelstag'])
DVF_Model2_result['a_0'] = 0.0
DVF_Model2_result['a_1'] = 0.0
DVF_Model2_result['a_2'] = 0.0
DVF_Model2_result['a_3'] = 0.0
DVF_Model2_result['a_5'] = 0.0
DVF_Model2_result['RMSE_S'] = 0.0
DVF_Model2_result['RMSE_S_Previous_Vs_Current'] = np.nan

DVF_Model3_result = pd.DataFrame(rawdata['Handelstag'].unique(), columns=['Handelstag'])
DVF_Model3_result['a_0'] = 0.0
DVF_Model3_result['a_1'] = 0.0
DVF_Model3_result['a_2'] = 0.0
DVF_Model3_result['a_3'] = 0.0
DVF_Model3_result['a_4'] = 0.0
DVF_Model3_result['a_5'] = 0.0
DVF_Model3_result['RMSE_S'] = 0.0
DVF_Model3_result['RMSE_S_Previous_Vs_Current'] = np.nan

import statsmodels.api as sm

for i in range(len(DVF_Model3_result)):   
    data = rawdata.loc[(rawdata['Handelstag'] == DVF_Model3_result['Handelstag'][i]) & (rawdata['OptionType'] == "Call") & (rawdata['Moneyness'] >= 0.8) & (rawdata['Moneyness'] <= 1.2) & (rawdata['t_delta'] >= 1/12) & (rawdata['t_delta'] <= 1)]
    data = data.reset_index()
    
    # Model 2 (a_0 + a_1 K + a_2 K^2 + a_3 T + a_5 K*T)
    X = np.column_stack((data['StrikePrice'], data['StrikePrice']**2, data['t_delta'], data['StrikePrice']*data['t_delta']))
    X = sm.add_constant(X)
    y = data['ImpliedVola']
    
    model = sm.OLS(y, X)
    results = model.fit()
        
    DVF_Model2_result['a_0'][i] = results.params[0]
    DVF_Model2_result['a_1'][i] = results.params[1]
    DVF_Model2_result['a_2'][i] = results.params[2]
    DVF_Model2_result['a_3'][i] = results.params[3]
    DVF_Model2_result['a_5'][i] = results.params[4]
        
    # Create RMSE
    DVF_Model2_result['RMSE_S'][i] = np.sqrt(np.mean((data['TaeglicherAbrechnungspreis'] - vBSMOption(data['OptionType'], data['SchlusspreisBasiswert'], data['t_delta'], data['StrikePrice'], data['EONIA'], vDVF_Model2_sigma(data['t_delta'], data['StrikePrice'], DVF_Model2_result['a_0'][i], DVF_Model2_result['a_1'][i], DVF_Model2_result['a_2'][i], DVF_Model2_result['a_3'][i], DVF_Model2_result['a_5'][i]), Greek='Price'))**2)) / data['SchlusspreisBasiswert'][0]
    
    if i > 0:
            DVF_Model2_result['RMSE_S_Previous_Vs_Current'][i] = np.sqrt(np.mean((data['TaeglicherAbrechnungspreis'] - vBSMOption(data['OptionType'], data['SchlusspreisBasiswert'], data['t_delta'], data['StrikePrice'], data['EONIA'], vDVF_Model2_sigma(data['t_delta'], data['StrikePrice'], DVF_Model2_result['a_0'][i-1], DVF_Model2_result['a_1'][i-1], DVF_Model2_result['a_2'][i-1], DVF_Model2_result['a_3'][i-1], DVF_Model2_result['a_5'][i-1]), Greek='Price'))**2)) / data['SchlusspreisBasiswert'][0]
    
    # Model 3 (a_0 + a_1 K + a_2 K^2 + a_3 T + a_4 T^2 + a_5 K*T)
    X = np.column_stack((data['StrikePrice'], data['StrikePrice']**2, data['t_delta'], data['t_delta']**2, data['StrikePrice']*data['t_delta']))
    X = sm.add_constant(X)
    y = data['ImpliedVola']
    
    model = sm.OLS(y, X)
    results = model.fit()
        
    DVF_Model3_result['a_0'][i] = results.params[0]
    DVF_Model3_result['a_1'][i] = results.params[1]
    DVF_Model3_result['a_2'][i] = results.params[2]
    DVF_Model3_result['a_3'][i] = results.params[3]
    DVF_Model3_result['a_4'][i] = results.params[4]
    DVF_Model3_result['a_5'][i] = results.params[5]
        
    # Create RMSE
    DVF_Model3_result['RMSE_S'][i] = np.sqrt(np.mean((data['TaeglicherAbrechnungspreis'] - vBSMOption(data['OptionType'], data['SchlusspreisBasiswert'], data['t_delta'], data['StrikePrice'], data['EONIA'], vDVF_Model3_sigma(data['t_delta'], data['StrikePrice'], DVF_Model3_result['a_0'][i], DVF_Model3_result['a_1'][i], DVF_Model3_result['a_2'][i], DVF_Model3_result['a_3'][i], DVF_Model3_result['a_4'][i], DVF_Model3_result['a_5'][i]), Greek='Price'))**2)) / data['SchlusspreisBasiswert'][0]

    if i > 0:
        DVF_Model3_result['RMSE_S_Previous_Vs_Current'][i] = np.sqrt(np.mean((data['TaeglicherAbrechnungspreis'] - vBSMOption(data['OptionType'], data['SchlusspreisBasiswert'], data['t_delta'], data['StrikePrice'], data['EONIA'], vDVF_Model3_sigma(data['t_delta'], data['StrikePrice'], DVF_Model3_result['a_0'][i-1], DVF_Model3_result['a_1'][i-1], DVF_Model3_result['a_2'][i-1], DVF_Model3_result['a_3'][i-1], DVF_Model3_result['a_4'][i-1], DVF_Model3_result['a_5'][i-1]), Greek='Price'))**2)) / data['SchlusspreisBasiswert'][0]

# Create Plots
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# Single day
import datetime 
Handelstag = datetime.datetime(2016, 5, 17)
data = rawdata.loc[(rawdata['Handelstag'] == Handelstag) & (rawdata['OptionType'] == "Call") & (rawdata['Moneyness'] >= 0.8) & (rawdata['Moneyness'] <= 1.2) & (rawdata['t_delta'] >= 1/12) & (rawdata['t_delta'] <= 1)]
data = data.reset_index()

# Fit data
X = np.column_stack((data['StrikePrice'], data['StrikePrice']**2, data['t_delta'], data['t_delta']**2, data['StrikePrice']*data['t_delta']))
X = sm.add_constant(X)
y = data['ImpliedVola']

model = sm.OLS(y, X)
results = model.fit()

# Create data for grid
strike = np.linspace(np.min(data['StrikePrice']), np.max(data['StrikePrice']), 51)
maturity = np.linspace(np.min(data['t_delta']), np.max(data['t_delta']), 51)
Strike, Maturity = np.meshgrid(strike, maturity)
ImpliedVola = DVF_Model3_sigma(Maturity, Strike, results.params[0], results.params[1], results.params[2], results.params[3], results.params[4], results.params[5])

# Plot
fig = plt.figure()
ax = plt.axes(projection='3d')
plt.title("DVF Model 3 - " + data['OptionType'][0] + " - " + Handelstag.strftime("%b %d %Y"))
ax.plot_wireframe(Strike, Maturity, ImpliedVola, color='grey')
ax.scatter(data['StrikePrice'], data['t_delta'], data['ImpliedVola'], c='red', s=3)
ax.set_xlabel('Strike')
ax.set_ylabel('Maturity')
ax.set_zlabel('Implied Volatility')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
ax.zaxis.set_major_locator(plt.MaxNLocator(7))

# Plot timeseries      
plt.plot('Handelstag', 'RMSE_S', data=DVF_Model2_result, linewidth=0.5, label="RMSE/S")
plt.plot('Handelstag', 'RMSE_S_Previous_Vs_Current', data=DVF_Model2_result, linewidth=0.5, label="RMSE/S Previous VsCurrent")
plt.legend(title="Model 2", loc=2, fontsize='small')
plt.show()

plt.plot('Handelstag', 'RMSE_S', data=DVF_Model3_result, linewidth=0.5, label="RMSE/S")
plt.plot('Handelstag', 'RMSE_S_Previous_Vs_Current', data=DVF_Model3_result, linewidth=0.5, label="RMSE/S Previous VsCurrent")
plt.legend(title="Model 3", loc=2, fontsize='small')
plt.show()

plt.plot('Handelstag', 'RMSE_S', data=DVF_Model2_result, linewidth=0.5, label="Model 2")
plt.plot('Handelstag', 'RMSE_S', data=DVF_Model3_result, linewidth=0.5, label="Model 3")
plt.legend(title="RMSE/S", loc=2, fontsize='small')
plt.show()

plt.plot('Handelstag', 'RMSE_S_Previous_Vs_Current', data=DVF_Model2_result, linewidth=0.5, label="Model 2")
plt.plot('Handelstag', 'RMSE_S_Previous_Vs_Current', data=DVF_Model3_result, linewidth=0.5, label="Model 3")
plt.legend(title="RMSE/S Steady State", loc=2, fontsize='small')
plt.show()

import numpy as np
import scipy.stats as si
import datetime 
import pandas as pd

print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

# to open or create a HDF5 file
rawdata = pd.read_hdf('D:\h_da\_Blog\Themen\Options\BSM_ImpliedVola_SingleOption_Python\EurexOptionsDax.h5')

# Select all calls & puts from one day
day = datetime.datetime(2016,1,4)
data = rawdata.loc[(rawdata['Handelstag'] == day)]

# Add Moneyness
def Moneyness(OptionType, SchlusspreisBasiswert, StrikePrice):
    if OptionType == 'Call':
        Value = SchlusspreisBasiswert/StrikePrice
    if OptionType == 'Put':
        Value = StrikePrice/SchlusspreisBasiswert
    return Value

data['Moneyness'] = data.apply(lambda data: Moneyness(data['OptionType'], data['SchlusspreisBasiswert'], data['StrikePrice']), axis=1)

# Add Maturity and t-delta
def Maturity(month, year):
    d = datetime.datetime(year,month,1)
    # https://docs.python.org/3/library/datetime.html#datetime.date.weekday
    # date.weekday(): Return the day of the week as an integer, where Monday is 0 and Sunday is 6. 
    d = d + datetime.timedelta(days=7) - datetime.timedelta(days=d.weekday()) - datetime.timedelta(days=3)
    if d.month != month:
        d = d + datetime.timedelta(days=7)
    d = d + datetime.timedelta(days=14)
    return d

def diff_dates(date1, date2):
    return abs(date2-date1).days

data['Maturity'] = data.apply(lambda data: Maturity(data['Verfalltermin_Monat'], data['Verfalltermin_Jahr']), axis=1)
data['t_Delta'] = data.apply(lambda data: diff_dates(data['Handelstag'], data['Maturity']), axis=1) / 365

data = data.loc[(data['Moneyness'] >= 0.5) & (data['Moneyness'] <= 2) & (data['t_Delta'] <= 2) & (data['t_Delta'] >= 1/52)]

### Calc implied vola
# Please note 'Four Things You Might Not Know About The Black-Scholes-Formula'
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

from scipy import optimize
data['ImpliedVola'] = -1
def BSMOptionFindImpliedVola(sigma, Price, PutCallFlag, S, T, K, r):
    return Price - BSMOption(PutCallFlag, S, T, K, r, sigma, Greek='Price')

for i in range(len(data)):
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.newton.html
    #scipy.optimize.bisect(f, a, b, args=(), xtol=2e-12, rtol=8.881784197001252e-16, maxiter=100, full_output=False, disp=True)[source]
    data['ImpliedVola'].iloc[i] = optimize.bisect(BSMOptionFindImpliedVola, a=-0.01, b=1, args=(data['TaeglicherAbrechnungspreis'].iloc[i], data['OptionType'].iloc[i], data['SchlusspreisBasiswert'].iloc[i], data['t_Delta'].iloc[i], data['StrikePrice'].iloc[i], data['EONIA'].iloc[i]), xtol=1e-8, maxiter=500)

calls = data.loc[(data['OptionType'] == "Call")]
puts = data.loc[(data['OptionType'] == "Put")]

### Plt 3D-Sccatterplot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# plot calls
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.title(day.strftime('Calls %Y-%m-%d'))
ax.scatter(calls['t_Delta'], calls['Moneyness'], calls['ImpliedVola'], c='r', s=5)
ax.set_xlabel('t')
ax.set_ylabel('Moneyness')
ax.set_zlabel('Implied Volatility')
ax.view_init(30, 120)
plt.savefig('D:/h_da/_Blog/Themen/Options/BSM_ImpliedVola_SingleOption_Python/Calls.png', dpi=180)
plt.show()
 
# plot puts
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.title(day.strftime('Puts %Y-%m-%d'))
ax.scatter(puts['t_Delta'], puts['Moneyness'], puts['ImpliedVola'], c='r', s=5)
ax.set_xlabel('t')
ax.set_ylabel('Moneyness')
ax.set_zlabel('Implied Volatility')
ax.view_init(30, 120)
plt.savefig('D:/h_da/_Blog/Themen/Options/BSM_ImpliedVola_SingleOption_Python/Puts.png', dpi=180)
plt.show()
 
# plot calls 1 with colormap
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.title(day.strftime('Calls %Y-%m-%d'))
cmhot = cm.get_cmap("hot")
c = np.abs(calls['ImpliedVola'])
pnt3d=ax.scatter(calls['t_Delta'], calls['Moneyness'], calls['ImpliedVola'], c=c, s=5,cmap=cmhot)
ax.set_xlabel('t')
ax.set_ylabel('Moneyness')
ax.set_zlabel('Implied Volatility')
cbar=plt.colorbar(pnt3d)
ax.view_init(30, 120)
plt.savefig('D:/h_da/_Blog/Themen/Options/BSM_ImpliedVola_SingleOption_Python/Calls1.png', dpi=180)
plt.show()
 
# plot calls 2 with colormap
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.title(day.strftime('Calls %Y-%m-%d'))
cmhot = cm.get_cmap("hot")
c = np.abs(calls['ImpliedVola'])
pnt3d=ax.scatter(calls['t_Delta'], calls['Moneyness'], calls['ImpliedVola'], c=c, s=5,cmap=cmhot)
ax.set_xlabel('t')
ax.set_ylabel('Moneyness')
ax.set_zlabel('Implied Volatility')
cbar=plt.colorbar(pnt3d)
ax.view_init(30, 315)
plt.savefig('D:/h_da/_Blog/Themen/Options/BSM_ImpliedVola_SingleOption_Python/Calls2.png', dpi=180)
plt.show()
 
# plot puts 1 with colormap
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.title(day.strftime('Puts %Y-%m-%d'))
cmhot = cm.get_cmap("hot")
c = np.abs(puts['ImpliedVola'])
pnt3d=ax.scatter(puts['t_Delta'], puts['Moneyness'], puts['ImpliedVola'], c=c, s=5,cmap=cmhot)
ax.set_xlabel('t')
ax.set_ylabel('Moneyness')
ax.set_zlabel('Implied Volatility')
cbar=plt.colorbar(pnt3d)
ax.view_init(30, 120)
plt.savefig('D:/h_da/_Blog/Themen/Options/BSM_ImpliedVola_SingleOption_Python/Puts1.png', dpi=180)
plt.show()

# plot puts 2 with colormap
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.title(day.strftime('Puts %Y-%m-%d'))
cmhot = cm.get_cmap("hot")
c = np.abs(puts['ImpliedVola'])
pnt3d=ax.scatter(puts['t_Delta'], puts['Moneyness'], puts['ImpliedVola'], c=c, s=5,cmap=cmhot)
ax.set_xlabel('t')
ax.set_ylabel('Moneyness')
ax.set_zlabel('Implied Volatility')
cbar=plt.colorbar(pnt3d)
ax.view_init(30, 315)
plt.savefig('D:/h_da/_Blog/Themen/Options/BSM_ImpliedVola_SingleOption_Python/Puts2.png', dpi=180)
plt.show()

# Hist
plt.style.use('bmh')
bins = np.linspace(0, 1, 39)

plt.hist([calls['ImpliedVola'], puts['ImpliedVola']], bins, label=['Calls', 'Puts'])
plt.legend(loc='upper right')
plt.savefig('D:/h_da/_Blog/Themen/Options/BSM_ImpliedVola_SingleOption_Python/Hist_Calls_Puts.png', dpi=180)
plt.show()

# Boxplot
plt.boxplot([calls[['ImpliedVola']], puts[['ImpliedVola']]])
plt.xticks([1, 2,], ['Calls', 'Puts'])
plt.savefig('D:/h_da/_Blog/Themen/Options/BSM_ImpliedVola_SingleOption_Python/Boxplot.png', dpi=180)
plt.show()

print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
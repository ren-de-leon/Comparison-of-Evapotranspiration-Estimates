
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sb
from scipy.stats.stats import pearsonr


data = pd.read_csv('min_temp.csv')
ds = pd.read_csv('max_temp.csv')
RH = pd.read_excel('RELATIVE HUMIDITY.xlsx')
ds1 = pd.read_excel('2013 (1).xlsx')
CD = pd.read_excel('2013 closed.xlsx')
dr = pd.date_range(start = '2013-01-01', end = '2013-12-31', freq= 'D' )
CD1 = pd.DataFrame(CD.values, index=dr)
CD1.columns=['Hnew','LEnew', 'Hm', 'LEm']


ds1['TIMESTAMP'] = pd.DatetimeIndex(ds1['TIMESTAMP'])
ds1 = ds1.set_index('TIMESTAMP')
ds1 = ds1.resample('H').mean()
ws = ds1['WS [m/s]']

rnett =  pd.read_csv("net_Radiation.csv") 
rnett.drop('Unnamed: 0', axis =1)

rnett['Date [annual]'] = pd.DatetimeIndex(rnett['Date [annual]'])
rnett = rnett.set_index('Date [annual]')
RH['Time[h]'] = pd.DatetimeIndex(RH['Time[h]'])
rnett = rnett.resample('H').mean()

new_data = pd.read_csv('new_data.csv')
new_data.drop('Unnamed: 0', axis =1)
new_data['TIME'] = pd.DatetimeIndex(new_data['TIME'])
new_data = new_data.set_index('TIME')
RH= RH.set_index(RH['Time[h]'])
new_data = new_data.resample('H').mean()
RH = RH.resample('H').mean()

WS =[]
for i , item in enumerate(ws.values):
    WS.append(item)

rh = []
for i , item in enumerate(RH.values):
    rh.append(item[0])
    
rnet = rnett.RNet
Rnet =[]
for i , item in enumerate(rnet.values):
    Rnet.append(item)

g = rnett.G
G =[]
for i , item in enumerate(g.values):
    G.append(item)
    
# T = data['T (Celcius)']
t = new_data['T (Celcius)']
T =[]
for i , item in enumerate(t.values):
    T.append(item)
    
Tx = ds.C
Tn = data['T (Celcius)']


# =============================================================================
#    HARGREAVES
# =============================================================================

lamb = []
for i, item in enumerate(T):
    l = (2.501 - (0.0236 * T[i]))
    lamb.append(l)
    
def Hargreaves(T, Tx, Tn, RNet, lamb):
    H=[]
    for i,item in enumerate(T):
        h = 0.0023*((T[i] + 17.8)*((Tx[i] - Tn[i])**0.5))* RNet[i] / (lamb[i])
        H.append(h)
    return H

H = Hargreaves(T, Tx, Tn, Rnet, lamb)
d = pd.DataFrame(H)    

# =============================================================================
#      EDDY COVARIANCE closed data 
# =============================================================================




l = CD1['LEnew']

l = l.astype(float)
le_med = np.nanmedian(l)
le = l.replace(np.nan, le_med)

LE =[]
for i , item in enumerate(le.values):
    LE.append(item)


new_data = pd.read_csv('new_data.csv')
new_data.drop('Unnamed: 0', axis =1)
new_data['TIME'] = pd.DatetimeIndex(new_data['TIME'])
new_data = new_data.set_index('TIME')
t = new_data['T (Celcius)']
T =[]
for i , item in enumerate(t.values):
    T.append(item)

lamb = []
for i, item in enumerate(T):
    l = (2.501 - (0.0236 * T[i]))
    lamb.append(l)

def formular(LE , lamb):
    pw = 997.77
    ECD = []
    sum1 = 0
    count = 0
    for i, item in enumerate(LE): 
        sum1 +=  (LE[i] / lamb[i])
        count +=1
        ECD.append(((1 / pw) * sum1) * 100)
        sum1 = 0
        count = 0
        # print(sum1, i)  
    return ECD 
ECD = formular(LE , lamb)

    
# =============================================================================
#   PRIESTLY TAYLOR 
# =============================================================================
def priestly_taylor (rnet, G):
    PT = []
    psyco = 0.054
    alpha = 1.3
    lv = 2453.0
    # slope = []
    for i ,item in enumerate(G):
        s = (4098 * (0.6108 * math.exp((17.27 * T[i])/(T[i] + 237.3))) )/ ((T[i] + 237.3)**2)
        et = (alpha * ((s)*(rnet[i] - G[i])/(lv*(s + psyco)))) * 1000
        et = et/10
        PT.append(et)
    return (PT)

PT = priestly_taylor(Rnet, G)

# def priestly_taylor2(rnet, G):
#     ET2=[]
#     alpha = 1.3
#     psyco = 0.054
#     lat_vap = 2.45
#     for i, item in enumerate(rnet):
#         s = (4098 * (0.6108 * math.exp((17.27 * T[i])/(T[i] + 237.3))) )/ ((T[i] + 237.3)**2)
#         et = alpha * (s * (rnet[i] - G[i]) / (lat_vap * (s + psyco))) 
#         ET2.append(et)
#     return (ET2)
    
# ET2 = priestly_taylor2(Rnet, G)
# d2 = pd.DataFrame(ET2)



# =============================================================================
#      EDDY COVARIANCE 
# =============================================================================


rnett =  pd.read_csv("net_Radiation.csv")
rnett.drop('Unnamed: 0', axis =1)
rnett['Date [annual]'] = pd.DatetimeIndex(rnett['Date [annual]'])
rnett = rnett.set_index('Date [annual]')

le = rnett['LE']
le_med = np.nanmedian(le)
le = le.replace(np.nan, le_med)

LE =[]
for i , item in enumerate(le.values):
    LE.append(item)


new_data = pd.read_csv('new_data.csv')
new_data.drop('Unnamed: 0', axis =1)
new_data['TIME'] = pd.DatetimeIndex(new_data['TIME'])
new_data = new_data.set_index('TIME')
t = new_data['T (Celcius)']
T =[]
for i , item in enumerate(t.values):
    T.append(item)

lamb = []
for i, item in enumerate(T):
    l = (2.501 - (0.0236 * T[i]))
    lamb.append(l)


def formular(LE , lamb):
    pw = 997.77
    EC = []
    sum1 = 0
    count = 0
    for i, item in enumerate(LE): 
        sum1 +=  (LE[i] / lamb[i])
        count +=1
        if count == 48:
            EC.append(((1 / pw) * sum1))
            sum1 = 0
            count = 0
            # print(sum1, i)  
    return EC
EC = formular(LE , lamb)


# =============================================================================
#      FAO 
# =============================================================================
ES =[]
for i, item in enumerate (T):
    es = 0.6108 * math.exp((17.27 * T[i]) / (T[i] + 273.3))
    ES.append(es)
    
    
EN=[]
for i, item in enumerate(rh):
    en = (rh[i] * ES[i]) / 100
    EN.append(en)


def FAO (rnet, G, T, windspeed , es , en):
    FAO = []
    psyco = 0.054
    for i, item in enumerate(rnet):
        s = (4098 * (0.6108 * math.exp((17.27 * T[i])/(T[i] + 237.3))) )/ ((T[i] + 237.3)**2)
        fao = (0.408 * s * (rnet[i] - G[i]) + psyco * (900/(T[i] +273)) * (es[i] - en[i]))  / (s + psyco * (1 + 0.34 * windspeed[i]))
        fao = fao / 4
        FAO.append(fao)
    return FAO
    
f = FAO(Rnet,G,T,WS,ES,EN)   





# =============================================================================
#           COMBINED DATA 
# =============================================================================

df = pd.date_range(start='2013-01-01 00:00:00' , end ='2013-12-31 23:00:00', periods = 8760)
df = pd.DataFrame(df, columns=['Date'])
df['Hargreaves'] = H
df['Priestly_Taylor'] = PT
df['FAO'] = f

df = df.set_index('Date')
# df1 = df.resample('H').mean()
df2 = df.resample('D').mean()
df2['Eddy_Cov'] = EC
#df2 = df.resample('D').mean()
df2['Eddy_Cov_closed'] = ECD

df3 = df2.resample('M').mean()

month = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# data_std_H = (df2.Hargreaves - (np.nanmean(df2.Hargreaves))) / np.nanstd(df2.Hargreaves)
# data_std_PT = (df2.Priestly_Taylor - (np.nanmean(df2.Priestly_Taylor))) / np.nanstd(df2.Priestly_Taylor)
# data_std_fao = (df2.FAO - (np.nanmean(df2.FAO))) / np.nanstd(df2.FAO)
# data_std_EC = (df2.Eddy_Cov - (np.nanmean(df2.Eddy_Cov))) / np.nanstd(df2.Eddy_Cov)



#==============================================================================
# #            PLOTTING FOR STD 
#  #=============================================================================
# fig, ax = plt.subplots(dpi = 300)
# ax.plot(df2.index, data_std_EC , color='b', label = 'Eddy_Cov')
# ax.plot(df2.index, data_std_H, color='k', label = 'Hargreaves')
# ax.plot(df2.index, data_std_fao,color= 'g', label = 'FAO')
# ax.plot(df2.index, data_std_PT,color= 'r', label = 'Priestly_Taylor')
# plt.legend(loc='lower right')


# =============================================================================
#      PLOTTING THE DATA FOR HOURLY 
# =============================================================================

# fig, ax = plt.subplots(dpi = 300)
# ax.plot(df1.index, df1.Hargreaves, color='k', label = 'Hargreaves')
# ax.plot(df1.index, df1.Priestly_Taylor, color='r', label = 'Priestly_Taylor')
# ax.plot(df1.index, df1.Eddy_Cov, color='b', label = 'Eddy_Cov')
# ax.plot(df1.index, df1.FAO, color='g', label = 'FAO')
# plt.legend(loc='lower right')


# =============================================================================
#      PLOTTING THE DATA 
# =============================================================================

fig, ax = plt.subplots(dpi = 300)
ax.plot(df2.index, df2.Hargreaves, color='k', label = 'Hargreaves')
ax.plot(df2.index, df2.Priestly_Taylor, color='r', label = 'Priestly_Taylor')
ax.plot(df2.index, df2.Eddy_Cov, color='b', label = 'Eddy_Cov')
ax.plot(df2.index, df2.Eddy_Cov_closed, color='indigo', label = 'Eddy_Cov_closed')
ax.plot(df2.index, df2.FAO, color='g', label = 'FAO')
# plt.legend(handles=['vertical_line'], loc='upper right')
plt.legend(bbox_to_anchor= (1, 1), loc='upper left')
plt.ylabel("Evapotranspiration ($mm day^{-1}$)")
# =============================================================================
#  monthly data
# =============================================================================

fig, ax = plt.subplots(dpi = 300)
ax.plot(month, df3.Hargreaves, color='k', label = 'Hargreaves')
ax.plot(month, df3.Priestly_Taylor, color='r', label = 'Priestly_Taylor')
ax.plot(month, df3.Eddy_Cov, color='b', label = 'Eddy_Cov')
ax.plot(month, df3.Eddy_Cov_closed, color='indigo', label = 'Eddy_Cov_closed')
ax.plot(month, df3.FAO, color='g', label = 'FAO')
# plt.legend(loc='upper left')
plt.legend(bbox_to_anchor= (1, 1), loc='upper left')

# #===============================================================================
# # fao and hargreaves (monthly)
# #==============================================================================

# # fig, ax = plt.subplots(dpi = 300)
# sb.regplot(df3.FAO, df3.Hargreaves, color='k')

# #==============================================================================
# # fao and Priestly_taylor (monthly)
# #==============================================================================
  
# sb.regplot(df3.FAO, df3.Priestly_Taylor, color='r' )

# #==============================================================================
# # fao and eddy_cov (monthly)
# #==============================================================================
  
# sb.regplot(df3.FAO, df3.Eddy_Cov, color='b' )




#  #===============================================================================
#  #    fao and hargreaves (daily)
# #==============================================================================

#     # fig, ax = plt.subplots(dpi = 300)
# sb.regplot(df2.FAO, df2.Hargreaves,color='k')

# #==============================================================================
# # fao and Priestly_taylor (daily)
# #==============================================================================
  
# sb.regplot(df2.FAO, df2.Priestly_Taylor, color='r' )



# #==============================================================================
# # fao and Eddy_Cov (daily)
# #==============================================================================
  
# sb.regplot(df2.FAO, df2.Eddy_Cov, color='b' )

#==============================================================================
# fao and eddy_cov (daily)
#==============================================================================
# fao = df2.FAO
# edd = df2.Eddy_Cov
# fig,ax = plt.subplot()
# sb.regplot(fao, edd, color='b',ax = ax )
# r,p = pearsonr(df2.FAO, df2.Eddy_Cov)
 
# ax.txt(0,4, f' r = {r}', horizontalalignment = 'center')

#==============================================================================
# error 
#=============================================================================
def isNan (num):
    return num != num

def arith_mean (x):
    count = 0
    sum = 0
    for item in x:
        if not isNan(item):  
            sum += item
            count += 1
    if count != 0:              
        return sum/count
    return 0

def sz_inv(x):
    return (np.size(x))**-1

def power(x,y):
    return x**y

def rms(x,y):
    return (np.sqrt(sz_inv(x)*np.nansum(power((x-y),2))))/arith_mean(y)

def bias(x,y):
    return np.nansum(x)/np.nansum(y)

#======================================================================
# l c for daily
#===================================================================
 
h = df2.Hargreaves
h_med = np.nanmedian(h)
h = h.replace(np.nan, h_med)

fao = df2.FAO
f_med = np.nanmedian(fao)
fao = fao.replace(np.nan, f_med)

ec = df2.Eddy_Cov
EC_med = np.nanmedian(ec)
ec = ec.replace(np.nan, EC_med)

ecd = df2.Eddy_Cov_closed
ECD_med = np.nanmedian(ecd)
ecd = ecd.replace(np.nan, ECD_med)

pt = df2.Priestly_Taylor
PT_med = np.nanmedian(pt)
pt = pt.replace(np.nan, PT_med)

r1,p = pearsonr(fao, h)
r2,p = pearsonr(fao, ec)
r3,p = pearsonr(fao, pt)
r4,p = pearsonr(fao, ecd)




#==============================================================================
#                bias, rmse & lc for daily values 
#==============================================================================

# fig, ax = plt.subplots()
# sb.regplot(df2.FAO, df2.Priestly_Taylor, color='r' , ax= ax)    # r-correlation, p-pvalue
# rmse = rms(df2.FAO, df2.Priestly_Taylor)        # root mean square error
# bias1 = bias(df2.FAO, df2.Priestly_Taylor)
# ax.text(0, 7,  f'bias = {round(bias1,2)}' , horizontalalignment='center', verticalalignment='center')
# ax.text(0, 6.5,  f'rmse = {round(rmse,2)}' , horizontalalignment='center', verticalalignment='center')



# fig, ax = plt.subplots()
# sb.regplot(df2.FAO, df2.Hargreaves, color='k', ax= ax)    # r-correlation, p-pvalue
# rmse = rms(df2.FAO, df2.Hargreaves)        # root mean square error
# bias1 = bias(df2.FAO, df2.Hargreaves)
# ax.text(0, 8,  f'bias = {round(bias1,2)}' , horizontalalignment='center', verticalalignment='center')
# ax.text(0, 7.5,  f'rmse = {round(rmse,2)}' , horizontalalignment='center', verticalalignment='center')

# fig, ax = plt.subplots()
# sb.regplot(df2.FAO, df2.Eddy_Cov, color='b', ax= ax)    # r-correlation, p-pvalue
# rmse = rms(df2.FAO, df2.Eddy_Cov)        # root mean square error
# bias1 = bias(df2.FAO, df2.Eddy_Cov)
# ax.text(0, 4,  f'bias = {round(bias1,2)}' , horizontalalignment='center', verticalalignment='center')
# ax.text(0, 3.5,  f'rmse = {round(rmse,2)}' , horizontalalignment='center', verticalalignment='center')



fig, ax = plt.subplots(figsize=(30,20), dpi=400)
# ax1 = fig.add_subplot(1,4,1)     # either use this one or the one below
ax1 = fig.add_subplot(2,2,1)        
sb.regplot(df2.FAO, df2.Eddy_Cov, color='b' )#, ax = ax1)    # r-correlation, p-pvalue
rmse = rms(df2.FAO, df2.Eddy_Cov)        # root mean square error
bias1 = bias(df2.FAO, df2.Eddy_Cov)
ax1.tick_params(axis ='x', labelsize=25)
ax1.tick_params(axis ='y', labelsize=25)
ax1.set_ylabel('Eddy Covariance', size =25)
ax1.set_xlabel('FAO', size =25)
ax1.text(0,  4,  f'bias = {round(bias1,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax1.text(0, 3.7,  f'rmse = {round(rmse,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax1.text(0, 3.4,  f'r = {round(r2,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)


# ax2 = fig.add_subplot(1,4,2)    # either use this one or the one below
ax2 = fig.add_subplot(2,2,2)
sb.regplot(df2.FAO, df2.Hargreaves, color='k' )#, ax = ax1)    # r-correlation, p-pvalue
rmse = rms(df2.FAO, df2.Hargreaves)        # root mean square error
bias1 = bias(df2.FAO, df2.Hargreaves)
ax2.tick_params(axis ='x', labelsize=25)
ax2.tick_params(axis ='y', labelsize=25)
ax2.set_ylabel('Hargreaves', size =25)
ax2.set_xlabel('FAO', size =25)
ax2.text(0, 8,  f'bias = {round(bias1,2)}' , horizontalalignment='center', verticalalignment='center',  size= 25)
ax2.text(0, 7.5,  f'rmse = {round(rmse,2)}' , horizontalalignment='center', verticalalignment='center',  size= 25)
ax2.text(0, 7.0,  f'r = {round(r1,2)}' , horizontalalignment='center', verticalalignment='center',  size= 25)


# ax3 = fig.add_subplot(1,4,3)     # either use this one or the one below
ax3 = fig.add_subplot(2,2,3)
sb.regplot(df2.FAO, df2.Priestly_Taylor, color='r')# , ax= ax)    # r-correlation, p-pvalue
rmse = rms(df2.FAO, df2.Priestly_Taylor)        # root mean square error
bias1 = bias(df2.FAO, df2.Priestly_Taylor)
ax3.tick_params(axis ='x', labelsize=25)
ax3.tick_params(axis ='y', labelsize=25)
ax3.set_ylabel('Priestly_Taylor', size =25)
ax3.set_xlabel('FAO', size =25)
ax3.text(0, 8,  f'bias = {round(bias1,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax3.text(0, 7.5,  f'rmse = {round(rmse,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax3.text(0, 7.0,  f'r = {round(r3,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)


# ax4 = fig.add_subplot(1,4,4)     # either use this one or the one below
ax4 = fig.add_subplot(2,2,4)
sb.regplot(df2.FAO, df2.Eddy_Cov_closed, color='indigo' )#, ax = ax1)    # r-correlation, p-pvalue
rmse = rms(df2.FAO, df2.Eddy_Cov_closed)        # root mean square error
bias1 = bias(df2.FAO, df2.Eddy_Cov_closed)
ax4.tick_params(axis ='x', labelsize=25)
ax4.tick_params(axis ='y', labelsize=25)
ax4.set_ylabel('Eddy Covariance closed', size =25)
ax4.set_xlabel('FAO', size =25)
ax4.text(0,  8,  f'bias = {round(bias1,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax4.text(0, 7,  f'rmse = {round(rmse,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax4.text(0, 6,  f'r = {round(r4,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)

#=========================================================================================================================
#     bias, rmse & lc for monthly values 
#=========================================================================================================================

h = df3.Hargreaves
h_med = np.nanmedian(h)
h = h.replace(np.nan, h_med)

fao = df3.FAO
f_med = np.nanmedian(fao)
fao = fao.replace(np.nan, f_med)

ec = df3.Eddy_Cov
EC_med = np.nanmedian(ec)
ec = ec.replace(np.nan, EC_med)

ecd = df3.Eddy_Cov_closed
ECD_med = np.nanmedian(ecd)
ecd = ecd.replace(np.nan, ECD_med)

pt = df3.Priestly_Taylor
PT_med = np.nanmedian(pt)
pt = pt.replace(np.nan, PT_med)

r1,p = pearsonr(fao, h)
r2,p = pearsonr(fao, ec)
r3,p = pearsonr(fao, pt)
r4,p = pearsonr(fao, ecd)

fig = plt.figure(figsize=(30,10), dpi=400)
ax1 = fig.add_subplot(1,3,1)
sb.regplot(df3.FAO, df3.Eddy_Cov, color='b' )#, ax = ax1)    # r-correlation, p-pvalue
rmse = rms(df3.FAO, df3.Eddy_Cov)        # root mean square error
bias1 = bias(df3.FAO, df3.Eddy_Cov)
ax1.tick_params(axis ='x', labelsize=25)
ax1.tick_params(axis ='y', labelsize=25)
ax1.set_ylabel('Eddy Covariance', size =25)
ax1.set_xlabel('FAO', size =25)
ax1.text(4,  2.1,  f'bias = {round(bias1,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax1.text(4, 2.0,  f'rmse = {round(rmse,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax1.text(4, 1.9,  f'r = {round(r2,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)

fig = plt.figure(figsize=(30,10), dpi=400)
ax1 = fig.add_subplot(1,3,1)
sb.regplot(df3.FAO, df3.Eddy_Cov_closed, color='indigo' )#, ax = ax1)    # r-correlation, p-pvalue
rmse = rms(df3.FAO, df3.Eddy_Cov_closed)        # root mean square error
bias1 = bias(df3.FAO, df3.Eddy_Cov_closed)
ax1.tick_params(axis ='x', labelsize=25)
ax1.tick_params(axis ='y', labelsize=25)
ax1.set_ylabel('Eddy Covariance closed', size =25)
ax1.set_xlabel('FAO', size =25)
ax1.text(4, 5,  f'bias = {round(bias1,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax1.text(4, 4.6,  f'rmse = {round(rmse,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax1.text(4, 4.2,  f'r = {round(r2,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)


ax2 = fig.add_subplot(1,3,2)
sb.regplot(df3.FAO, df3.Hargreaves, color='k' )#, ax = ax1)    # r-correlation, p-pvalue
rmse = rms(df3.FAO, df3.Hargreaves)        # root mean square error
bias1 = bias(df3.FAO, df3.Hargreaves)
ax2.tick_params(axis ='x', labelsize=25)
ax2.tick_params(axis ='y', labelsize=25)
ax2.set_ylabel('Hargreaves', size =25)
ax2.set_xlabel('FAO', size =25)
ax2.text(4, 5.8,  f'bias = {round(bias1,2)}' , horizontalalignment='center', verticalalignment='center',  size= 25)
ax2.text(4, 5.6,  f'rmse = {round(rmse,2)}' , horizontalalignment='center', verticalalignment='center',  size= 25)
ax2.text(4, 5.4,  f'r = {round(r1,2)}' , horizontalalignment='center', verticalalignment='center',  size= 25)


ax3 = fig.add_subplot(1,3,3)
sb.regplot(df3.FAO, df3.Priestly_Taylor, color='r')# , ax= ax)    # r-correlation, p-pvalue
rmse = rms(df3.FAO, df3.Priestly_Taylor)        # root mean square error
bias1 = bias(df3.FAO, df3.Priestly_Taylor)
ax3.tick_params(axis ='x', labelsize=25)
ax3.tick_params(axis ='y', labelsize=25)
ax3.set_ylabel('Priestly_Taylor', size =25)
ax3.set_xlabel('FAO', size =25)
ax3.text(4, 5.4,  f'bias = {round(bias1,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax3.text(4, 5.2,  f'rmse = {round(rmse,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)
ax3.text(4, 5.0,  f'r = {round(r3,2)}' , horizontalalignment='center', verticalalignment='center', size= 25)

#=============================================================================
# 
#=============================================================================

df3['m'] = [1,2,3,4,5,6,7,8,9,10,11,12]
df3.to_csv('df3.csv')

bdf3 = pd.read_csv('df3.csv')

width = 0.30
month = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

fig, ax3 = plt.subplots(dpi=300)
ax3.bar(bdf3.m - width, bdf3.FAO, width=0.30, color='g',align='center',label = 'fao')
ax3.bar(bdf3.m + width, bdf3.Eddy_Cov, width=0.30, color='b',align='center',label = 'Eddy_cov')
ax3.bar(bdf3.m , bdf3.Eddy_Cov_closed, width=0.30, color='indigo',align='center',label = 'Eddy_cov_closed')
plt.ylabel("Evapotranspiration ($mm day^{-1}$)")
plt.xlabel("Months of 2013")
plt.xticks(df3.m+width/3, month)
# plt.legend(loc ='best')
plt.legend(bbox_to_anchor= (1, 1), loc='upper left')

fig, ax3 = plt.subplots(dpi=300)
ax3.bar(bdf3.m - width, bdf3.FAO, width=0.30, color='g',align='center',label = 'fao')
ax3.bar(bdf3.m + width, bdf3.Hargreaves, width=0.30, color='k',align='center',label = 'Hargreaves')
ax3.bar(bdf3.m , bdf3.Priestly_Taylor, width=0.30, color='r',align='center',label = 'Priestly_Taylor')
plt.ylabel("Evapotranspiration ($mm day^{-1}$)")
plt.xlabel("Months of 2013")
plt.xticks(df3.m+width/3, month)
# plt.legend(loc ='best')
plt.legend(bbox_to_anchor= (1, 1), loc='upper left')



df3.to_csv('run')
al = pd.read_csv ('run')
np. mean(al)

df2.to_csv('all.csv')

pr =pd.read_csv('all.csv')
np.max(pr)

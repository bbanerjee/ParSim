
# coding: utf-8

# In[1]:


import numpy as np
from sklearn.svm import SVR
from sklearn import preprocessing
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV
import pandas as pd
import matplotlib.pyplot as plt


# In[2]:


data = pd.read_csv("./Sand_hydrostat.csv", header=0, skiprows=1)


# In[3]:


X = data.iloc[:,0]
y = data.iloc[:,1]
scaler = preprocessing.MinMaxScaler()


# In[4]:


clf = SVR(C=100.0, epsilon=0.2)
clf1000 = SVR(C=100000.0, epsilon=0.2)


# In[5]:


X_reshaped = X.values.reshape(-1,1)


# In[6]:


clf.fit(scaler.fit_transform(X_reshaped), y)
clf1000.fit(scaler.fit_transform(X_reshaped), y)


# In[7]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(X, y, '--', label='Hydrostat')
plt.plot(X, clf.predict(scaler.fit_transform(X_reshaped)), label='SVR fit (C=100)')
plt.plot(X, clf1000.predict(scaler.fit_transform(X_reshaped)), label='SVR fit (C=100000)')
plt.axis([0, 50, 0, 5000])
plt.xlabel('Engineering volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_Hydrostat.svg')


# In[8]:


unload_09 = pd.read_csv("./Sand_unload_09.csv", header=0, skiprows=1)
unload_18 = pd.read_csv("./Sand_unload_18.csv", header=0, skiprows=1)
unload_27 = pd.read_csv("./Sand_unload_27.csv", header=0, skiprows=1)
unload_36 = pd.read_csv("./Sand_unload_36.csv", header=0, skiprows=1)
unload_45 = pd.read_csv("./Sand_unload_45.csv", header=0, skiprows=1)


# In[9]:


eps_p_09 = unload_09.iloc[-1,0]
eps_e_09 = unload_09.iloc[:,0] - eps_p_09
eps_p_18 = unload_18.iloc[-1,0]
eps_e_18 = unload_18.iloc[:,0] - eps_p_18
eps_p_27 = unload_27.iloc[-1,0]
eps_e_27 = unload_27.iloc[:,0] - eps_p_27
eps_p_36 = unload_36.iloc[-1,0]
eps_e_36 = unload_36.iloc[:,0] - eps_p_36
eps_p_45 = unload_45.iloc[-1,0]
eps_e_45 = unload_45.iloc[:,0] - eps_p_45
p_09 = unload_09.iloc[:,1]
p_18 = unload_18.iloc[:,1]
p_27 = unload_27.iloc[:,1]
p_36 = unload_36.iloc[:,1]
p_45 = unload_45.iloc[:,1]


# In[10]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(eps_e_09, p_09, 'k-', label='Unload@9.04 strain')
plt.plot(eps_e_18, p_18, 'C0', label='Unload@18.08 strain')
plt.plot(eps_e_27, p_27, 'C3', label='Unload@27.12 strain')
plt.plot(eps_e_36, p_36, 'C1', label='Unload@36.16 strain')
plt.plot(eps_e_45, p_45, 'C4', label='Unload@45.20 strain')
plt.axis([0, 8, 0, 2000])
plt.xlabel('Elastic volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_ElasticStrain.svg')


# In[11]:


strains_09 = np.column_stack((eps_e_09, np.repeat(eps_p_09, eps_e_09.shape[0])))
strains_18 = np.column_stack((eps_e_18, np.repeat(eps_p_18, eps_e_18.shape[0])))
strains_27 = np.column_stack((eps_e_27, np.repeat(eps_p_27, eps_e_27.shape[0])))
strains_36 = np.column_stack((eps_e_36, np.repeat(eps_p_36, eps_e_36.shape[0])))
strains_45 = np.column_stack((eps_e_45, np.repeat(eps_p_45, eps_e_45.shape[0])))


# In[12]:


strains = np.concatenate((strains_09, strains_18, strains_27, strains_36, strains_45), axis=0)


# In[13]:


pressures = np.concatenate((p_09, p_18, p_27, p_36, p_45), axis=0)


# In[14]:


pressures


# In[15]:


curve_fitter = SVR(C=100.0, epsilon=0.2)


# In[16]:


curve_fitter.fit(scaler.fit_transform(strains), pressures)


# In[17]:


pressures_pred = curve_fitter.predict(scaler.fit_transform(strains_45))


# In[18]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(eps_e_09, p_09, 'k-', label='Unload@9.04 strain')
plt.plot(eps_e_18, p_18, 'C0', label='Unload@18.08 strain')
plt.plot(eps_e_27, p_27, 'C3', label='Unload@27.12 strain')
plt.plot(eps_e_36, p_36, 'C1', label='Unload@36.16 strain')
plt.plot(eps_e_45, p_45, 'C4', label='Unload@45.20 strain')
plt.plot(strains_45[:,0], pressures_pred, 'C6', label='SVR (C=100)')
plt.axis([0, 8, 0, 2000])
plt.xlabel('Elastic volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)


# In[19]:


data_09 = np.column_stack((strains_09, p_09))
data_18 = np.column_stack((strains_18, p_18))
data_27 = np.column_stack((strains_27, p_27))
data_36 = np.column_stack((strains_36, p_36))
data_45 = np.column_stack((strains_45, p_45))
data_all =  np.concatenate((data_09, data_18, data_27, data_36, data_45), axis=0)                       


# In[20]:


np.random.shuffle(data_all)


# In[21]:


data_all


# In[22]:


strains_shuffle = data_all[:,(0,1)]


# In[23]:


pressures_shuffle = data_all[:,2]


# In[52]:


fitter_shuffle = SVR(kernel='rbf', C=100000, epsilon=0.2)
fitter_shuffle.fit(scaler.fit_transform(strains_shuffle), pressures_shuffle)


# In[53]:


strains_scaled = scaler.fit_transform(strains)
start_09 = 0
end_09 = strains_09.shape[0]
start_18 = end_09
end_18 = start_18 + strains_18.shape[0]
start_27 = end_18
end_27 = start_27 + strains_27.shape[0]
start_36 = end_27
end_36 = start_36 + strains_36.shape[0]
start_45 = end_36
end_45 = start_45 + strains_45.shape[0]
pressures_pred_09 = fitter_shuffle.predict(strains_scaled[start_09:end_09,:])
pressures_pred_18 = fitter_shuffle.predict(strains_scaled[start_18:end_18,:])
pressures_pred_27 = fitter_shuffle.predict(strains_scaled[start_27:end_27,:])
pressures_pred_36 = fitter_shuffle.predict(strains_scaled[start_36:end_36,:])
pressures_pred_45 = fitter_shuffle.predict(strains_scaled[start_45:end_45,:])


# In[56]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.semilogy(eps_e_09, p_09, 'k--', label='Unload@9.04 strain')
plt.semilogy(eps_e_18, p_18, 'C0--', label='Unload@18.08 strain')
plt.semilogy(eps_e_27, p_27, 'C3--', label='Unload@27.12 strain')
plt.semilogy(eps_e_36, p_36, 'C1--', label='Unload@36.16 strain')
plt.semilogy(eps_e_45, p_45, 'C4--', label='Unload@45.20 strain')
plt.semilogy(strains_09[:,0], pressures_pred_09, 'k', linewidth=2, label='SVR (C=100000, eps=0.2)')
plt.semilogy(strains_18[:,0], pressures_pred_18, 'C0', linewidth=2)
plt.semilogy(strains_27[:,0], pressures_pred_27, 'C3', linewidth=2)
plt.semilogy(strains_36[:,0], pressures_pred_36, 'C1', linewidth=2)
plt.semilogy(strains_45[:,0], pressures_pred_45, 'C4', linewidth=2)
plt.axis([0, 8, 0, 2000])
plt.xlabel('Elastic volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_ElasticStrain_SVR_log.svg')


# In[27]:


#C_range = np.logspace(-2, 10, 13)
#epsilon_range = np.logspace(-9, 3, 13)
#param_grid = dict(epsilon=epsilon_range, C=C_range)


# In[28]:


#fits = []
#for C in C_range:
#    for epsilon in epsilon_range:
#        clf1 = SVR(C = C, epsilon=epsilon)
#        clf1.fit(scaler.fit_transform(strains_shuffle), pressures_shuffle)
#        fits.append((C, epsilon, clf1))    


# In[29]:


#scaler.fit_transform(strains_shuffle)


# In[30]:


#strains_shuffle


# In[57]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(eps_e_09, p_09, 'k--', label='Unload@9.04 strain')
plt.plot(eps_e_18, p_18, 'C0--', label='Unload@18.08 strain')
plt.plot(eps_e_27, p_27, 'C3--', label='Unload@27.12 strain')
plt.plot(eps_e_36, p_36, 'C1--', label='Unload@36.16 strain')
plt.plot(eps_e_45, p_45, 'C4--', label='Unload@45.20 strain')
plt.plot(strains_09[:,0], pressures_pred_09, 'k', linewidth=2, label='SVR (C=100000, eps=0.2)')
plt.plot(strains_18[:,0], pressures_pred_18, 'C0', linewidth=2)
plt.plot(strains_27[:,0], pressures_pred_27, 'C3', linewidth=2)
plt.plot(strains_36[:,0], pressures_pred_36, 'C1', linewidth=2)
plt.plot(strains_45[:,0], pressures_pred_45, 'C4', linewidth=2)
plt.axis([0, 8, 0, 2000])
plt.xlabel('Elastic volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_ElasticStrain_SVR.svg')


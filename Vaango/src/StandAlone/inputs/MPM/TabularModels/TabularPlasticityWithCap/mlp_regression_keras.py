
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
from keras.models import Sequential
from keras.layers import Dense
from keras.utils import plot_model
import matplotlib.pyplot as plt
from sklearn import preprocessing


# In[2]:


data = pd.read_csv("./Sand_hydrostat.csv", header=0, skiprows=1)


# In[3]:


train = data
X = train.iloc[:,0]
y = train.iloc[:,1]
scaler = preprocessing.MinMaxScaler()


# In[4]:


shuffled = np.column_stack((X,y))
print(shuffled.shape)
np.random.shuffle(shuffled)
shuffled


# In[5]:


def baseline_model():
    model = Sequential()
    model.add(Dense(32, input_dim=1, kernel_initializer='uniform', activation='relu'))
    model.add(Dense(64, kernel_initializer='normal'))
    model.add(Dense(1, kernel_initializer='uniform'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    return model


# In[6]:


model = baseline_model()


# In[7]:


X2D = shuffled[:,0].reshape(-1,1)
X2D.shape
y2D = shuffled[:,1]
model.fit(scaler.fit_transform(X2D), y2D, batch_size=192, epochs=5000, verbose=1, validation_split=0.2, shuffle=True)


# In[8]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(X, y, '--', label='Hydrostat', linewidth=3)
plt.plot(X, model.predict(scaler.fit_transform(X.values.reshape(-1,1))), label='MLP fit')
plt.axis([0, 50, 0, 5000])
plt.xlabel('Engineering volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_Hydrostat_MLP.svg')


# In[45]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(X, y, '--', label='Hydrostat', linewidth=3)
plt.plot(X, model.predict(scaler.fit_transform(X.values.reshape(-1,1))), label='MLP fit')
plt.axis([0, 30, 0, 500])
plt.xlabel('Engineering volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_Hydrostat_MLP_Zoom.svg')


# In[9]:


unload_09 = pd.read_csv("./Sand_unload_09.csv", header=0, skiprows=1)
unload_18 = pd.read_csv("./Sand_unload_18.csv", header=0, skiprows=1)
unload_27 = pd.read_csv("./Sand_unload_27.csv", header=0, skiprows=1)
unload_36 = pd.read_csv("./Sand_unload_36.csv", header=0, skiprows=1)
unload_45 = pd.read_csv("./Sand_unload_45.csv", header=0, skiprows=1)


# In[10]:


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
p_09_0 = unload_09.iloc[-1,1]
p_18_0 = unload_18.iloc[-1,1]
p_27_0 = unload_27.iloc[-1,1]
p_36_0 = unload_36.iloc[-1,1]
p_45_0 = unload_45.iloc[-1,1]
p_09 = unload_09.iloc[:,1] - p_09_0
p_18 = unload_18.iloc[:,1] - p_18_0
p_27 = unload_27.iloc[:,1] - p_27_0
p_36 = unload_36.iloc[:,1] - p_36_0
p_45 = unload_45.iloc[:,1] - p_45_0


# In[11]:


print(p_09_0, p_18_0, p_27_0, p_36_0, p_45_0)
p_45


# In[12]:


eps_p_30 = 0.5*(eps_p_27 + eps_p_36)
print(eps_p_27,eps_p_36, eps_p_30)


# In[13]:


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


# In[14]:


strains_09 = np.column_stack((eps_e_09, np.repeat(eps_p_09, eps_e_09.shape[0])))
strains_18 = np.column_stack((eps_e_18, np.repeat(eps_p_18, eps_e_18.shape[0])))
strains_27 = np.column_stack((eps_e_27, np.repeat(eps_p_27, eps_e_27.shape[0])))
strains_36 = np.column_stack((eps_e_36, np.repeat(eps_p_36, eps_e_36.shape[0])))
strains_45 = np.column_stack((eps_e_45, np.repeat(eps_p_45, eps_e_45.shape[0])))
strains_30 = np.column_stack((eps_e_45, np.repeat(eps_p_30, eps_e_45.shape[0])))


# In[15]:


print(strains_45)
print(p_45)


# In[16]:


strains = np.concatenate((strains_09, strains_18, strains_27, strains_36, strains_45), axis=0)
pressures = np.concatenate((p_09, p_18, p_27, p_36, p_45), axis=0)


# In[17]:


data_09 = np.column_stack((strains_09, p_09))
data_18 = np.column_stack((strains_18, p_18))
data_27 = np.column_stack((strains_27, p_27))
data_36 = np.column_stack((strains_36, p_36))
data_45 = np.column_stack((strains_45, p_45))
data_all =  np.concatenate((data_09, data_18, data_27, data_36, data_45), axis=0)    
np.random.shuffle(data_all)


# In[18]:


strains_shuffle = data_all[:,(0,1)]
pressures_shuffle = data_all[:,2]


# In[39]:


def baseline2D_model():
    model = Sequential()
    model.add(Dense(64, input_dim=2, kernel_initializer='normal', activation='sigmoid'))
    model.add(Dense(32, kernel_initializer='normal', activation='sigmoid'))
    model.add(Dense(32, kernel_initializer='normal', activation='relu'))
    model.add(Dense(1, kernel_initializer='normal'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    return model


# In[40]:


model2D = baseline2D_model()
print(strains_shuffle.shape)
print(pressures_shuffle.shape)


# In[41]:


model2D.fit(scaler.fit_transform(strains_shuffle), pressures_shuffle, batch_size=200, epochs=5000, verbose=1, validation_split=0.2, shuffle=True)


# In[42]:


strains = np.concatenate((strains, strains_30), axis=0)
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
start_30 = end_45
end_30 = start_30 + strains_30.shape[0]
pressures_pred_09 = model2D.predict(strains_scaled[start_09:end_09,:])
pressures_pred_18 = model2D.predict(strains_scaled[start_18:end_18,:])
pressures_pred_27 = model2D.predict(strains_scaled[start_27:end_27,:])
pressures_pred_36 = model2D.predict(strains_scaled[start_36:end_36,:])
pressures_pred_45 = model2D.predict(strains_scaled[start_45:end_45,:])
pressures_pred_30 = model2D.predict(strains_scaled[start_30:end_30,:])


# In[46]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(eps_e_09, p_09, 'k--', label='Unload@9.04 strain')
plt.plot(eps_e_18, p_18, 'C0--', label='Unload@18.08 strain')
plt.plot(eps_e_27, p_27, 'C3--', label='Unload@27.12 strain')
plt.plot(eps_e_36, p_36, 'C1--', label='Unload@36.16 strain')
plt.plot(eps_e_45, p_45, 'C4--', label='Unload@45.20 strain')
plt.plot(strains_09[:,0], pressures_pred_09, 'k', linewidth=2, label='MLP')
plt.plot(strains_18[:,0], pressures_pred_18, 'C0', linewidth=2)
plt.plot(strains_27[:,0], pressures_pred_27, 'C3', linewidth=2)
plt.plot(strains_36[:,0], pressures_pred_36, 'C1', linewidth=2)
plt.plot(strains_45[:,0], pressures_pred_45, 'C4', linewidth=2)
plt.plot(strains_30[:,0], pressures_pred_30, 'C7', linewidth=2, label="MLP (plastic strain = 0.3)")
plt.axis([0, 5, 0, 500])
plt.xlabel('Elastic volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_ElasticStrain_MLP.svg')


# In[24]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.semilogy(eps_e_09, p_09, 'k--', label='Unload@9.04 strain')
plt.semilogy(eps_e_18, p_18, 'C0--', label='Unload@18.08 strain')
plt.semilogy(eps_e_27, p_27, 'C3--', label='Unload@27.12 strain')
plt.semilogy(eps_e_36, p_36, 'C1--', label='Unload@36.16 strain')
plt.semilogy(eps_e_45, p_45, 'C4--', label='Unload@45.20 strain')
plt.semilogy(strains_09[:,0], pressures_pred_09, 'k', linewidth=2, label='MLP')
plt.semilogy(strains_18[:,0], pressures_pred_18, 'C0', linewidth=2)
plt.semilogy(strains_27[:,0], pressures_pred_27, 'C3', linewidth=2)
plt.semilogy(strains_36[:,0], pressures_pred_36, 'C1', linewidth=2)
plt.semilogy(strains_45[:,0], pressures_pred_45, 'C4', linewidth=2)
plt.semilogy(strains_30[:,0], pressures_pred_30, 'C7', linewidth=2)
plt.axis([0, 8, 0, 2000])
plt.xlabel('Elastic volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_ElasticStrain_MLP_log.svg')


# In[25]:


plot_model(model2D, show_shapes=True, to_file='MLP_model2D.png')


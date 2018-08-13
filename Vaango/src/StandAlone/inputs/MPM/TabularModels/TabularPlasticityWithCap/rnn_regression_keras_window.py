
# coding: utf-8

# In[1]:


import keras
import tensorflow as tf
import sys
import numpy as np
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pylab as plt
from sklearn import preprocessing
from scipy.interpolate import interp1d

print("python:{}, keras:{}, tensorflow: {}".format(sys.version, keras.__version__, tf.__version__))


# In[2]:


from keras.models import Sequential  
from keras.layers.core import Dense, Activation  
from keras.layers.recurrent import SimpleRNN


# In[3]:


data = pd.read_csv("./Sand_hydrostat.csv", header=0, skiprows=1)
x_scaler = preprocessing.MinMaxScaler()
y_scaler = preprocessing.MinMaxScaler()
X = data.iloc[:,0]
y = data.iloc[:,1]
X_scaled = x_scaler.fit_transform(X.reshape(-1,1))
X_scaled.reshape(X.shape[0],)
X_scaled.shape
y_scaled = y_scaler.fit_transform(y.reshape(-1,1))
y_scaled.reshape(X.shape[0],)
y_scaled.shape


# In[4]:


n_timesteps = 3
n_size = X.shape[0]
n_samples = int(n_size/n_timesteps)
indices = []
X_samples = []
y_samples = []
for sample in range(0, n_samples):
    #print("start=", sample*n_timesteps, "end = ", sample*n_timesteps+n_timesteps-1)
    indices.append(np.arange(sample*n_timesteps, sample*n_timesteps+n_timesteps))
    X_samples.append(X_scaled[indices[sample],0])
    y_samples.append(y_scaled[indices[sample],0])
X_samples = np.asarray(X_samples)
y_samples = np.asarray(y_samples)
print(X_samples.shape)
#print(indices)


# In[5]:


in_out_neurons = 1
hidden_neurons = 128


# In[6]:


X_rnn = X_samples.reshape(n_samples, n_timesteps, 1)
y_rnn = y_samples.reshape(n_samples, n_timesteps, 1)
print(X_rnn.shape)
print(y_rnn.shape)


# In[7]:


model = Sequential()
model.add(SimpleRNN(hidden_neurons, activation='relu', return_sequences=True, input_shape=(n_timesteps, in_out_neurons)))
model.add(Dense(in_out_neurons))
#model.add(Activation('linear'))
model.compile(loss="mean_squared_error", optimizer="adam")


# In[8]:


model.summary()


# In[9]:


model.fit(X_rnn, y_rnn, batch_size=5, epochs = 200)


# In[10]:


predicted = model.predict(X_samples[0].reshape(1, n_timesteps, 1))


# In[11]:


print(X_samples.shape)


# In[12]:


fig = plt.figure(figsize=[6,6])
ax = fig.add_subplot(111)
plt.plot(x_scaler.inverse_transform(X_samples[sample].reshape(n_timesteps,1)), y_scaler.inverse_transform(predicted.reshape(n_timesteps,1)), 'C0', label="RNN fit")
for sample in range(1, n_samples):
    pred = model.predict(X_samples[sample].reshape(1, n_timesteps, 1))
    plt.plot(x_scaler.inverse_transform(X_samples[sample].reshape(n_timesteps,1)), y_scaler.inverse_transform(pred.reshape(n_timesteps,1)), 'C0')
plt.plot(X, y, 'C1-', label='Hydrostat', linewidth=1)
plt.axis([0, 50, 0, 5000])
plt.xlabel('Engineering volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_Hydrostat_RNN.svg')


# In[33]:


fig = plt.figure(figsize=[6,6])
ax = fig.add_subplot(111)
plt.plot(x_scaler.inverse_transform(X_samples[sample].reshape(n_timesteps,1)), y_scaler.inverse_transform(predicted.reshape(n_timesteps,1)), 'C0', label="RNN fit")
for sample in range(1, n_samples):
    pred = model.predict(X_samples[sample].reshape(1, n_timesteps, 1))
    plt.plot(x_scaler.inverse_transform(X_samples[sample].reshape(n_timesteps,1)), y_scaler.inverse_transform(pred.reshape(n_timesteps,1)), 'C0')
plt.plot(X, y, 'C1-', label='Hydrostat', linewidth=1)
plt.axis([40, 45, 2000, 5000])
plt.xlabel('Engineering volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_Hydrostat_RNN_Zoom_Hi.svg')


# In[32]:


fig = plt.figure(figsize=[6,6])
ax = fig.add_subplot(111)
plt.plot(x_scaler.inverse_transform(X_samples[sample].reshape(n_timesteps,1)), y_scaler.inverse_transform(predicted.reshape(n_timesteps,1)), 'C0', label="RNN fit")
for sample in range(1, n_samples):
    pred = model.predict(X_samples[sample].reshape(1, n_timesteps, 1))
    plt.plot(x_scaler.inverse_transform(X_samples[sample].reshape(n_timesteps,1)), y_scaler.inverse_transform(pred.reshape(n_timesteps,1)), 'C0')
plt.plot(X, y, 'C1-', label='Hydrostat', linewidth=1)
plt.axis([0, 10, 0, 40])
plt.xlabel('Engineering volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_Hydrostat_RNN_Zoom_Lo.svg')


# In[13]:


unload_09 = pd.read_csv("./Sand_unload_09.csv", header=0, skiprows=1)
unload_18 = pd.read_csv("./Sand_unload_18.csv", header=0, skiprows=1)
unload_27 = pd.read_csv("./Sand_unload_27.csv", header=0, skiprows=1)
unload_36 = pd.read_csv("./Sand_unload_36.csv", header=0, skiprows=1)
unload_45 = pd.read_csv("./Sand_unload_45.csv", header=0, skiprows=1)
unload_09 = unload_09.sort_values(by=['(%)'],axis=0)
unload_18 = unload_18.sort_values(by=['(%)'],axis=0)
unload_27 = unload_27.sort_values(by=['(%)'],axis=0)
unload_36 = unload_36.sort_values(by=['(%)'],axis=0)
unload_45 = unload_45.sort_values(by=['(%)'],axis=0)


# In[14]:


eps_p_09 = unload_09.iloc[0,0]
eps_e_09 = unload_09.iloc[:,0] - eps_p_09
eps_p_18 = unload_18.iloc[0,0]
eps_e_18 = unload_18.iloc[:,0] - eps_p_18
eps_p_27 = unload_27.iloc[0,0]
eps_e_27 = unload_27.iloc[:,0] - eps_p_27
eps_p_36 = unload_36.iloc[0,0]
eps_e_36 = unload_36.iloc[:,0] - eps_p_36
eps_p_45 = unload_45.iloc[0,0]
eps_e_45 = unload_45.iloc[:,0] - eps_p_45
p_09 = unload_09.iloc[:,1]
p_18 = unload_18.iloc[:,1]
p_27 = unload_27.iloc[:,1]
p_36 = unload_36.iloc[:,1]
p_45 = unload_45.iloc[:,1]


# In[15]:


eps_09_scaler = preprocessing.MinMaxScaler()
eps_18_scaler = preprocessing.MinMaxScaler()
eps_27_scaler = preprocessing.MinMaxScaler()
eps_36_scaler = preprocessing.MinMaxScaler()
eps_45_scaler = preprocessing.MinMaxScaler()
p_09_scaler = preprocessing.MinMaxScaler()
p_18_scaler = preprocessing.MinMaxScaler()
p_27_scaler = preprocessing.MinMaxScaler()
p_36_scaler = preprocessing.MinMaxScaler()
p_45_scaler = preprocessing.MinMaxScaler()


# In[16]:


eps_e_09_scaled = eps_09_scaler.fit_transform(eps_e_09.values.reshape(-1,1)).reshape(eps_e_09.shape[0],)
eps_e_18_scaled = eps_18_scaler.fit_transform(eps_e_18.values.reshape(-1,1)).reshape(eps_e_18.shape[0],)
eps_e_27_scaled = eps_27_scaler.fit_transform(eps_e_27.values.reshape(-1,1)).reshape(eps_e_27.shape[0],)
eps_e_36_scaled = eps_36_scaler.fit_transform(eps_e_36.values.reshape(-1,1)).reshape(eps_e_36.shape[0],)
eps_e_45_scaled = eps_45_scaler.fit_transform(eps_e_45.values.reshape(-1,1)).reshape(eps_e_45.shape[0],)
eps_p_09_scaled = eps_09_scaler.fit_transform(eps_p_09)
eps_p_18_scaled = eps_18_scaler.fit_transform(eps_p_18)
eps_p_27_scaled = eps_27_scaler.fit_transform(eps_p_27)
eps_p_36_scaled = eps_36_scaler.fit_transform(eps_p_36)
eps_p_45_scaled = eps_45_scaler.fit_transform(eps_p_45)
p_09_scaled = p_09_scaler.fit_transform(p_09.values.reshape(-1,1)).reshape(p_09.shape[0],)
p_18_scaled = p_18_scaler.fit_transform(p_18.values.reshape(-1,1)).reshape(p_18.shape[0],)
p_27_scaled = p_27_scaler.fit_transform(p_27.values.reshape(-1,1)).reshape(p_27.shape[0],)
p_36_scaled = p_36_scaler.fit_transform(p_36.values.reshape(-1,1)).reshape(p_36.shape[0],)
p_45_scaled = p_45_scaler.fit_transform(p_45.values.reshape(-1,1)).reshape(p_45.shape[0],)


# In[17]:


n_eps_timesteps = 3
n_eps_size = p_45.shape[0]
n_eps_samples = int(n_eps_size/n_eps_timesteps)
eps_e_resampled = np.linspace(0, 1, n_eps_size-1)
p_09_resampled = np.interp(eps_e_resampled, eps_e_09_scaled, p_09_scaled)
p_18_resampled = np.interp(eps_e_resampled, eps_e_18_scaled, p_18_scaled)
p_27_resampled = np.interp(eps_e_resampled, eps_e_27_scaled, p_27_scaled)
p_36_resampled = np.interp(eps_e_resampled, eps_e_36_scaled, p_36_scaled)
p_45_resampled = np.interp(eps_e_resampled, eps_e_45_scaled, p_45_scaled)


# In[18]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(eps_e_resampled, p_09_resampled, 'k-', label='Unload@9.04 strain')
plt.plot(eps_e_resampled, p_18_resampled, 'C0', label='Unload@18.08 strain')
plt.plot(eps_e_resampled, p_27_resampled, 'C3', label='Unload@27.12 strain')
plt.plot(eps_e_resampled, p_36_resampled, 'C1', label='Unload@36.16 strain')
plt.plot(eps_e_resampled, p_45_resampled, 'C4', label='Unload@45.20 strain')
#plt.axis([0, 8, 0, 2000])
plt.xlabel('Elastic volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)


# In[19]:


strains_09 = np.column_stack((eps_e_resampled, np.repeat(eps_p_09, eps_e_resampled.shape[0])))
strains_18 = np.column_stack((eps_e_resampled, np.repeat(eps_p_18, eps_e_resampled.shape[0])))
strains_27 = np.column_stack((eps_e_resampled, np.repeat(eps_p_27, eps_e_resampled.shape[0])))
strains_36 = np.column_stack((eps_e_resampled, np.repeat(eps_p_36, eps_e_resampled.shape[0])))
strains_45 = np.column_stack((eps_e_resampled, np.repeat(eps_p_45, eps_e_resampled.shape[0])))


# In[20]:


print(strains_09.shape)
strains = np.concatenate((strains_09, strains_18, strains_27, strains_36, strains_45), axis=0)
pressures = np.concatenate((p_09_resampled, p_18_resampled, p_27_resampled, p_36_resampled, p_45_resampled), axis=0)


# In[21]:


eps_indices = []
eps_samples = []
p_samples = []
for sample in range(0, n_eps_samples*n_eps_timesteps):
    #print("start=", sample*n_eps_timesteps, "end = ", sample*n_eps_timesteps+n_eps_timesteps-1)
    eps_indices.append(np.arange(sample*n_eps_timesteps, sample*n_eps_timesteps+n_eps_timesteps))
    eps_samples.append(strains[eps_indices[sample]])
    p_samples.append(pressures[eps_indices[sample]])
eps_samples = np.asarray(eps_samples)
p_samples = np.asarray(p_samples)
p_samples = p_samples.reshape(p_samples.shape[0], p_samples.shape[1], 1)
print(eps_samples.shape)
print(p_samples.shape)


# In[22]:


eps_in_neurons = 2
eps_out_neurons = 1
eps_hidden_neurons = 128
eps_model = Sequential()
eps_model.add(SimpleRNN(hidden_neurons, activation='relu', return_sequences=True, input_shape=(n_eps_timesteps, eps_in_neurons)))
eps_model.add(Dense(eps_out_neurons))
#model.add(Activation('linear'))
eps_model.compile(loss="mean_squared_error", optimizer="adam")


# In[23]:


print(eps_samples)
#print(p_samples.shape)


# In[24]:


eps_model.summary()


# In[25]:


eps_model.fit(eps_samples, p_samples, batch_size=5, epochs = 200)


# In[26]:


print(n_eps_samples, strains_09.shape)
pred_09 = eps_model.predict(strains_09[0:n_eps_timesteps].reshape(1, n_eps_timesteps, 2))


# In[27]:


pred_09


# In[28]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(eps_e_resampled, p_09_resampled, 'k-', label='Unload@9.04 strain')
plt.plot(eps_e_resampled, p_18_resampled, 'C0', label='Unload@18.08 strain')
plt.plot(eps_e_resampled, p_27_resampled, 'C3', label='Unload@27.12 strain')
plt.plot(eps_e_resampled, p_36_resampled, 'C1', label='Unload@36.16 strain')
plt.plot(eps_e_resampled, p_45_resampled, 'C4', label='Unload@45.20 strain')
for sample in range(0, n_eps_samples):
    start = sample*n_eps_timesteps
    end = sample*n_eps_timesteps + n_eps_timesteps
    pred_09 = eps_model.predict(strains_09[start:end].reshape(1, n_eps_timesteps, 2))
    pred_18 = eps_model.predict(strains_18[start:end].reshape(1, n_eps_timesteps, 2))
    pred_27 = eps_model.predict(strains_27[start:end].reshape(1, n_eps_timesteps, 2))
    pred_36 = eps_model.predict(strains_36[start:end].reshape(1, n_eps_timesteps, 2))
    pred_45 = eps_model.predict(strains_45[start:end].reshape(1, n_eps_timesteps, 2))
    #print(pred)
    #print("start = ", start, "end = ", end, "strains_09 = ", strains_09[start:end].shape)
    plt.plot(strains_09[start:end,0], pred_09.reshape(n_eps_timesteps,1), 'k')
    plt.plot(strains_18[start:end,0], pred_18.reshape(n_eps_timesteps,1), 'C0')
    plt.plot(strains_27[start:end,0], pred_27.reshape(n_eps_timesteps,1), 'C3')
    plt.plot(strains_36[start:end,0], pred_36.reshape(n_eps_timesteps,1), 'C1')
    plt.plot(strains_45[start:end,0], pred_45.reshape(n_eps_timesteps,1), 'C4')
    
    
plt.xlabel('Normalized elastic volumetric strain', fontsize=16)
plt.ylabel('Normalized pressure', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_ElasticStrain_RNN.svg')


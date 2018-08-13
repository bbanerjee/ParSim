
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import PolylineIntersection as pl
from keras.models import Sequential
from keras.layers import Dense
from keras.utils import plot_model
import matplotlib.pyplot as plt
from sklearn import preprocessing
from keras.models import load_model
import h5py


# In[2]:


# Load the CSV data into Pandas dataframes
hydrostat = pd.read_csv("./DrySand_Hydrostat.csv", header=0, skiprows=7)
data_09 = pd.read_csv("./DrySand_LoadUnload_09.csv", header=0, skiprows=4)
data_18 = pd.read_csv("./DrySand_LoadUnload_18.csv", header=0, skiprows=4)
data_27 = pd.read_csv("./DrySand_LoadUnload_27.csv", header=0, skiprows=4)
data_36 = pd.read_csv("./DrySand_LoadUnload_36.csv", header=0, skiprows=4)
data_45 = pd.read_csv("./DrySand_LoadUnload_45.csv", header=0, skiprows=4)


# In[3]:


# Rename the columns of each dataframe
column_names = ["TotalStrainVol", "Pressure", "", "", "", ""]
hydrostat.columns = column_names
data_09.columns = column_names
data_18.columns = column_names
data_27.columns = column_names
data_36.columns = column_names
data_45.columns = column_names


# In[4]:


# Convert percent into strain, MPa to Pa
strain_fac = 0.01
pressure_fac = 1.0e6
hydrostat.TotalStrainVol *= strain_fac
hydrostat.Pressure *= pressure_fac
data_09.TotalStrainVol *= strain_fac
data_09.Pressure *= pressure_fac
data_18.TotalStrainVol *= strain_fac
data_18.Pressure *= pressure_fac
data_27.TotalStrainVol *= strain_fac
data_27.Pressure *= pressure_fac
data_36.TotalStrainVol *= strain_fac
data_36.Pressure *= pressure_fac
data_45.TotalStrainVol *= strain_fac
data_45.Pressure *= pressure_fac


# In[5]:


# Find the point at which unloading begins
p_max_09 = max(data_09.Pressure)
p_max_index_09 = data_09.Pressure.values.tolist().index(p_max_09)
p_max_18 = max(data_18.Pressure)
p_max_index_18 = data_18.Pressure.values.tolist().index(p_max_18)
p_max_27 = max(data_27.Pressure)
p_max_index_27 = data_27.Pressure.values.tolist().index(p_max_27)
p_max_36 = max(data_36.Pressure)
p_max_index_36 = data_36.Pressure.values.tolist().index(p_max_36)
p_max_45 = max(data_45.Pressure)
p_max_index_45 = data_45.Pressure.values.tolist().index(p_max_45)
p_max_index_00 = (np.abs(p_max_09 - hydrostat.Pressure.values)).argmin()
p_max_00 = hydrostat.Pressure.values[p_max_index_00]


# In[6]:


# Create separate dataframes for the unload data
data_09_unload = data_09[p_max_index_09:].copy()
data_18_unload = data_18[p_max_index_18:].copy()
data_27_unload = data_27[p_max_index_27:].copy()
data_36_unload = data_36[p_max_index_36:].copy()
data_45_unload = data_45[p_max_index_45:].copy()


# In[7]:


# Find plastic strains by intersecting the unload data with the pressure axis
pressure_axis = ((-1, 0),(1, 0))
poly_09_unload = list(data_09_unload[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
line_09_unload = (poly_09_unload[-1], poly_09_unload[-2])
plastic_strain_09 = pl.line_intersection(pressure_axis, line_09_unload)[0]
poly_18_unload = list(data_18_unload[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
line_18_unload = (poly_18_unload[-1], poly_18_unload[-2])
plastic_strain_18 = pl.line_intersection(pressure_axis, line_18_unload)[0]
poly_27_unload = list(data_27_unload[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
line_27_unload = (poly_27_unload[-1], poly_27_unload[-2])
plastic_strain_27 = pl.line_intersection(pressure_axis, line_27_unload)[0]
poly_36_unload = list(data_36_unload[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
line_36_unload = (poly_36_unload[-1], poly_36_unload[-2])
plastic_strain_36 = pl.line_intersection(pressure_axis, line_36_unload)[0]
poly_45_unload = list(data_45_unload[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
line_45_unload = (poly_45_unload[-1], poly_45_unload[-2])
plastic_strain_45 = pl.line_intersection(pressure_axis, line_45_unload)[0]
print(plastic_strain_09, plastic_strain_18, plastic_strain_27, plastic_strain_36, plastic_strain_45)


# In[8]:


# Reverse the order of the unload data to create elastic loading curves
data_00_load = hydrostat[:p_max_index_00]
data_09_load = data_09_unload.sort_index(ascending=False)
data_18_load = data_18_unload.sort_index(ascending=False)
data_27_load = data_27_unload.sort_index(ascending=False)
data_36_load = data_36_unload.sort_index(ascending=False)
data_45_load = data_45_unload.sort_index(ascending=False)


# In[9]:


# Simplify the loading dataframes and remove duplicates (if any)
data_00_all = data_00_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_09_all = data_09_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_18_all = data_18_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_27_all = data_27_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_36_all = data_36_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_45_all = data_45_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()


# In[10]:


# Compute max strain and pressure
total_strain_max = data_45_all['TotalStrainVol'].max()
pressure_max = data_45_all['Pressure'].max()


# In[11]:


# Rename strain variables
eps_p_09 = plastic_strain_09
eps_09 = data_09_all['TotalStrainVol'].values
eps_p_18 = plastic_strain_18
eps_18 = data_18_all['TotalStrainVol'].values
eps_p_27 = plastic_strain_27
eps_27 = data_27_all['TotalStrainVol'].values
eps_p_36 = plastic_strain_36
eps_36 = data_36_all['TotalStrainVol'].values
eps_p_45 = plastic_strain_45
eps_45 = data_45_all['TotalStrainVol'].values
eps_00 = np.linspace(0, total_strain_max)
eps_30 = np.linspace(0, total_strain_max)


# In[12]:


eps_p_00 = 0
eps_p_30 = 0.5*(eps_p_27 + eps_p_36)
print(eps_p_27,eps_p_36, eps_p_30)


# In[13]:


# Rename pressure variables
p_09 = data_09_all['Pressure'].values
p_18 = data_18_all['Pressure'].values
p_27 = data_27_all['Pressure'].values
p_36 = data_36_all['Pressure'].values
p_45 = data_45_all['Pressure'].values
print(p_09[1])


# In[14]:


# Scale the data
def scaled(x, min_x, max_x):
    return (x - min_x)/(max_x - min_x)

# Unscale the data
def unscaled(x, min_x, max_x):
    return min_x + x * (max_x - min_x)


# In[15]:


# Convert strains to scaled values
eps_09_scaled = scaled(eps_09, 0, total_strain_max)
eps_18_scaled = scaled(eps_18, 0, total_strain_max)
eps_27_scaled = scaled(eps_27, 0, total_strain_max)
eps_36_scaled = scaled(eps_36, 0, total_strain_max)
eps_45_scaled = scaled(eps_45, 0, total_strain_max)
eps_00_scaled = scaled(eps_00, 0, total_strain_max)
eps_30_scaled = scaled(eps_30, 0, total_strain_max)
eps_p_09_scaled = scaled(eps_p_09, 0, total_strain_max)
eps_p_18_scaled = scaled(eps_p_18, 0, total_strain_max)
eps_p_27_scaled = scaled(eps_p_27, 0, total_strain_max)
eps_p_36_scaled = scaled(eps_p_36, 0, total_strain_max)
eps_p_45_scaled = scaled(eps_p_45, 0, total_strain_max)
eps_p_00_scaled = scaled(eps_p_00, 0, total_strain_max)
eps_p_30_scaled = scaled(eps_p_30, 0, total_strain_max)


# In[16]:


# Convert pressures to scaled values
#p_09_scaled = scaled(p_09, 0, pressure_max)
#p_18_scaled = scaled(p_18, 0, pressure_max)
#p_27_scaled = scaled(p_27, 0, pressure_max)
#p_36_scaled = scaled(p_36, 0, pressure_max)
#p_45_scaled = scaled(p_45, 0, pressure_max)
p_09_scaled = scaled(p_09, 0, 1.0e6)
p_18_scaled = scaled(p_18, 0, 1.0e6)
p_27_scaled = scaled(p_27, 0, 1.0e6)
p_36_scaled = scaled(p_36, 0, 1.0e6)
p_45_scaled = scaled(p_45, 0, 1.0e6)


# In[17]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(eps_09, p_09, 'k-', label='Unload@9.04 strain')
plt.plot(eps_18, p_18, 'C0', label='Unload@18.08 strain')
plt.plot(eps_27, p_27, 'C3', label='Unload@27.12 strain')
plt.plot(eps_36, p_36, 'C1', label='Unload@36.16 strain')
plt.plot(eps_45, p_45, 'C4', label='Unload@45.20 strain')
#plt.axis([0, 40, 0, 2000])
plt.xlabel('Total volumetric strain', fontsize=16)
plt.ylabel('Pressure (Pa)', fontsize=16)
ax.legend(loc='best', fontsize=14)


# In[18]:


strains_09_scaled = np.column_stack((eps_09_scaled, np.repeat(eps_p_09_scaled, eps_09_scaled.shape[0])))
strains_18_scaled = np.column_stack((eps_18_scaled, np.repeat(eps_p_18_scaled, eps_18_scaled.shape[0])))
strains_27_scaled = np.column_stack((eps_27_scaled, np.repeat(eps_p_27_scaled, eps_27_scaled.shape[0])))
strains_36_scaled = np.column_stack((eps_36_scaled, np.repeat(eps_p_36_scaled, eps_36_scaled.shape[0])))
strains_45_scaled = np.column_stack((eps_45_scaled, np.repeat(eps_p_45_scaled, eps_45_scaled.shape[0])))
strains_30_scaled = np.column_stack((eps_30_scaled, np.repeat(eps_p_30_scaled, eps_30_scaled.shape[0])))
strains_00_scaled = np.column_stack((eps_00_scaled, np.repeat(eps_p_00_scaled, eps_00_scaled.shape[0])))


# In[19]:


strains_scaled = np.concatenate((strains_09_scaled, strains_18_scaled, strains_27_scaled, strains_36_scaled, strains_45_scaled), axis=0)
pressures_scaled = np.concatenate((p_09_scaled, p_18_scaled, p_27_scaled, p_36_scaled, p_45_scaled), axis=0)


# In[20]:


data_09_scaled = np.column_stack((strains_09_scaled, p_09_scaled))
data_18_scaled = np.column_stack((strains_18_scaled, p_18_scaled))
data_27_scaled = np.column_stack((strains_27_scaled, p_27_scaled))
data_36_scaled = np.column_stack((strains_36_scaled, p_36_scaled))
data_45_scaled = np.column_stack((strains_45_scaled, p_45_scaled))
data_all_scaled =  np.concatenate((data_09_scaled, data_18_scaled, data_27_scaled, data_36_scaled, data_45_scaled, 
                            data_09_scaled, data_09_scaled, data_09_scaled, data_09_scaled,
                            data_18_scaled, data_18_scaled, data_18_scaled, data_18_scaled,
                            data_27_scaled, data_27_scaled, data_27_scaled, data_27_scaled,
                            data_36_scaled, data_36_scaled, data_36_scaled, data_36_scaled,
                            data_45_scaled, data_45_scaled, data_45_scaled, data_45_scaled), axis=0)    
#np.random.seed(12345)
np.random.shuffle(data_all_scaled)


# In[21]:


strains_shuffle = data_all_scaled[:,(0,1)]
pressures_shuffle = data_all_scaled[:,2]


# In[22]:


def baseline2D_model():
    model = Sequential()
    model.add(Dense(64, input_dim=2, kernel_initializer='normal', activation='sigmoid'))
    model.add(Dense(32, kernel_initializer='normal', activation='sigmoid'))
    model.add(Dense(32, kernel_initializer='normal', activation='relu'))
    model.add(Dense(1, kernel_initializer='normal'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    return model


# In[23]:


model2D = baseline2D_model()
print(strains_shuffle.shape)
print(pressures_shuffle.shape)


# In[24]:


model2D.fit(strains_shuffle, pressures_shuffle, batch_size=192, epochs=3000, verbose=1, validation_split=0.2, shuffle=True)


# In[25]:


pressures_pred_09_scaled = model2D.predict(strains_09_scaled)
pressures_pred_18_scaled = model2D.predict(strains_18_scaled)
pressures_pred_27_scaled = model2D.predict(strains_27_scaled)
pressures_pred_36_scaled = model2D.predict(strains_36_scaled)
pressures_pred_45_scaled = model2D.predict(strains_45_scaled)
pressures_pred_30_scaled = model2D.predict(strains_30_scaled)
pressures_pred_00_scaled = model2D.predict(strains_00_scaled)


# In[26]:


strains_09 = unscaled(strains_09_scaled, 0, total_strain_max)
strains_18 = unscaled(strains_18_scaled, 0, total_strain_max)
strains_27 = unscaled(strains_27_scaled, 0, total_strain_max)
strains_36 = unscaled(strains_36_scaled, 0, total_strain_max)
strains_45 = unscaled(strains_45_scaled, 0, total_strain_max)
strains_30 = unscaled(strains_30_scaled, 0, total_strain_max)
strains_00 = unscaled(strains_00_scaled, 0, total_strain_max)


# In[27]:


pressure_max = 1.0e6
pressures_pred_09 = unscaled(pressures_pred_09_scaled, 0, pressure_max)
pressures_pred_18 = unscaled(pressures_pred_18_scaled, 0, pressure_max)
pressures_pred_27 = unscaled(pressures_pred_27_scaled, 0, pressure_max)
pressures_pred_36 = unscaled(pressures_pred_36_scaled, 0, pressure_max)
pressures_pred_45 = unscaled(pressures_pred_45_scaled, 0, pressure_max)
pressures_pred_30 = unscaled(pressures_pred_30_scaled, 0, pressure_max)
pressures_pred_00 = unscaled(pressures_pred_00_scaled, 0, pressure_max)


# In[28]:


lab_09 = str("Plastic strain = {0:.3f}".format(eps_p_09))
lab_18 = str("Plastic strain = {0:.3f}".format(eps_p_18))
lab_27 = str("Plastic strain = {0:.3f}".format(eps_p_27))
lab_36 = str("Plastic strain = {0:.3f}".format(eps_p_36))
lab_45 = str("Plastic strain = {0:.3f}".format(eps_p_45))
lab_00 = str("MLP: Plastic strain = {0:.3f}".format(eps_p_00))
lab_30 = str("MLP: Plastic strain = {0:.3f}".format(eps_p_30))
type(lab_09)


# In[29]:


def rescale(data, VARIABLE_TYPE):
    result = {
        "e" : lambda data : data*100.0,
        "p" : lambda data : data*1.0e-6
    }[VARIABLE_TYPE](data)
    return result


# In[30]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(rescale(eps_09, "e"), rescale(p_09, "p"), 'k--', label=lab_09)
plt.plot(rescale(eps_18, "e"), rescale(p_18, "p"), 'C0--', label=lab_18)
plt.plot(rescale(eps_27, "e"), rescale(p_27, "p"), 'C3--', label=lab_27)
plt.plot(rescale(eps_36, "e"), rescale(p_36, "p"), 'C1--', label=lab_36)
plt.plot(rescale(eps_45, "e"), rescale(p_45, "p"), 'C4--', label=lab_45)
plt.plot(rescale(strains_09[:,0], "e"), rescale(pressures_pred_09, "p"), 'k', linewidth=2, label='MLP')
plt.plot(rescale(strains_18[:,0], "e"), rescale(pressures_pred_18, "p"), 'C0', linewidth=2)
plt.plot(rescale(strains_27[:,0], "e"), rescale(pressures_pred_27, "p"), 'C3', linewidth=2)
plt.plot(rescale(strains_36[:,0], "e"), rescale(pressures_pred_36, "p"), 'C1', linewidth=2)
plt.plot(rescale(strains_45[:,0], "e"), rescale(pressures_pred_45, "p"), 'C4', linewidth=2)
plt.plot(rescale(strains_30[:,0], "e"), rescale(pressures_pred_30, "p"), 'C7', linewidth=2, label=lab_30)
plt.plot(rescale(strains_00[:,0], "e"), rescale(pressures_pred_00, "p"), 'C8', linewidth=2, label=lab_00)
plt.axis([0, 50, 0, 3000])
plt.xlabel('Total volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_ElasticStrain_MLP_Total.svg')


# In[31]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(rescale(eps_09, "e"), rescale(p_09, "p"), 'k--', label=lab_09)
plt.plot(rescale(strains_09[:,0], "e"), rescale(pressures_pred_09, "p"), 'k', linewidth=2, label='MLP fit')
plt.xlabel('Total volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
plt.grid(True)
fig.savefig('Fox_DrySand_ElasticStrain_MLP_Total_09.svg')


# In[32]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(rescale(eps_18, "e"), rescale(p_18, "p"), 'C0--', label=lab_18)
plt.plot(rescale(strains_18[:,0], "e"), rescale(pressures_pred_18, "p"), 'C0', linewidth=2, label='MLP fit')
plt.xlabel('Total volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
plt.grid(True)
fig.savefig('Fox_DrySand_ElasticStrain_MLP_Total_18.svg')


# In[33]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(rescale(eps_27, "e"), rescale(p_27, "p"), 'C3--', label=lab_27)
plt.plot(rescale(strains_27[:,0], "e"), rescale(pressures_pred_27, "p"), 'C3', linewidth=2, label='MLP fit')
plt.xlabel('Total volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
plt.grid(True)
fig.savefig('Fox_DrySand_ElasticStrain_MLP_Total_27.svg')


# In[34]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(rescale(eps_36, "e"), rescale(p_36, "p"), 'C1--', label=lab_36)
plt.plot(rescale(strains_36[:,0], "e"), rescale(pressures_pred_36, "p"), 'C1', linewidth=2, label='MLP fit')
plt.xlabel('Total volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
plt.grid(True)
fig.savefig('Fox_DrySand_ElasticStrain_MLP_Total_36.svg')


# In[35]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(rescale(eps_45, "e"), rescale(p_45, "p"), 'C4--', label=lab_45)
plt.plot(rescale(strains_45[:,0], "e"), rescale(pressures_pred_45, "p"), 'C4', linewidth=2, label='MLP fit')
plt.xlabel('Total volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
plt.grid(True)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_ElasticStrain_MLP_Total_45.svg')


# In[36]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(rescale(eps_27, "e"), rescale(p_27, "p"), 'C3--', label=lab_27)
plt.plot(rescale(eps_36, "e"), rescale(p_36, "p"), 'C1--', label=lab_36)
plt.plot(rescale(strains_30[:,0], "e"), rescale(pressures_pred_30, "p"), 'C7', linewidth=2, label=lab_30)
plt.axis([0, 40, -100, 1500])
plt.xlabel('Total volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
ax.legend(loc='best', fontsize=14)
plt.grid(True)
fig.savefig('Fox_DrySand_ElasticStrain_MLP_Total_30.svg')


# In[37]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(rescale(eps_09, "e"), rescale(p_09, "p"), 'k--', label=lab_09)
plt.plot(rescale(strains_00[:,0], "e"), rescale(pressures_pred_00, "p"), 'C8', linewidth=2, label=lab_00)
plt.axis([0, 20, -50, 500])
plt.xlabel('Total volumetric strain (%)', fontsize=16)
plt.ylabel('Pressure (MPa)', fontsize=16)
plt.grid(True)
ax.legend(loc='best', fontsize=14)
fig.savefig('Fox_DrySand_ElasticStrain_MLP_Total_00.svg')


# In[38]:


model2D.save('mlp_regression_keras_total_scaled.h5')


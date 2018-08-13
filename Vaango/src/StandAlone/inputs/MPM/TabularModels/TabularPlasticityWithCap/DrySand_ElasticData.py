
# coding: utf-8

# In[1]:


# %load DrySand_ElasticData.py
#  Compression is positive
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import PolylineIntersection as pl
import imp
pl = imp.reload(pl)


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
print(p_max_45, ",", p_max_index_45, ",", p_max_00, ",", p_max_index_00)


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


# Compute the volumetric elastic strains by subtracting plastic strains from total strains
data_09_unload.TotalStrainVol -= plastic_strain_09
data_18_unload.TotalStrainVol -= plastic_strain_18
data_27_unload.TotalStrainVol -= plastic_strain_27
data_36_unload.TotalStrainVol -= plastic_strain_36
data_45_unload.TotalStrainVol -= plastic_strain_45


# In[9]:


#print(data_09_unload)


# In[10]:


# Reverse the order of the unload data to create elastic loading curves
data_00_load = hydrostat[:p_max_index_00]
data_09_load = data_09_unload.sort_index(ascending=False)
data_18_load = data_18_unload.sort_index(ascending=False)
data_27_load = data_27_unload.sort_index(ascending=False)
data_36_load = data_36_unload.sort_index(ascending=False)
data_45_load = data_45_unload.sort_index(ascending=False)


# In[11]:


#print(data_09_load)


# In[12]:


# Create the unload data for 0 plastic strain from the loading data
#data_00_unload = data_00_load.sort_index(ascending=False)


# In[13]:


# Create tension data based on the unloading data
#data_00_tension = data_00_unload.copy()
#data_00_tension.TotalStrainVol *= (-1)
#data_00_tension.Pressure *= (-1)
#data_09_tension = data_09_unload.copy()
#data_09_tension.TotalStrainVol *= (-1)
#data_09_tension.Pressure *= (-1)
#data_18_tension = data_18_unload.copy()
#data_18_tension.TotalStrainVol *= (-1)
#data_18_tension.Pressure *= (-1)
#data_27_tension = data_27_unload.copy()
#data_27_tension.TotalStrainVol *= (-1)
#data_27_tension.Pressure *= (-1)
#data_36_tension = data_36_unload.copy()
#data_36_tension.TotalStrainVol *= (-1)
#data_36_tension.Pressure *= (-1)
#data_45_tension = data_45_unload.copy()
#data_45_tension.TotalStrainVol *= (-1)
#data_45_tension.Pressure *= (-1)


# In[14]:


#print(data_09_tension)


# In[15]:


# Add the tension to the compression data
#data_00_all = pd.concat([data_00_tension, data_00_load], ignore_index=True)
#data_09_all = pd.concat([data_09_tension, data_09_load], ignore_index=True)
#data_18_all = pd.concat([data_18_tension, data_18_load], ignore_index=True)
#data_27_all = pd.concat([data_27_tension, data_27_load], ignore_index=True)
#data_36_all = pd.concat([data_36_tension, data_36_load], ignore_index=True)
#data_45_all = pd.concat([data_45_tension, data_45_load], ignore_index=True)


# In[16]:


#print(data_09_all[400:450])


# In[17]:


# Simplify the loading dataframes and remove duplicates (if any)
data_00_all = data_00_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_09_all = data_09_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_18_all = data_18_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_27_all = data_27_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_36_all = data_36_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_45_all = data_45_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()


# In[18]:


# Compute the total volumetric strain
data_09_all.TotalStrainVol += plastic_strain_09
data_18_all.TotalStrainVol += plastic_strain_18
data_27_all.TotalStrainVol += plastic_strain_27
data_36_all.TotalStrainVol += plastic_strain_36
data_45_all.TotalStrainVol += plastic_strain_45


# In[19]:


#print(data_00_all)


# In[20]:


# Create an extra data set for tension interpolations
data_09_tension_all = data_00_all.copy()
data_09_tension_all.TotalStrainVol -= plastic_strain_09


# In[33]:


# Visual check of data set
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(data_09_tension_all.TotalStrainVol, data_09_tension_all.Pressure, 'C7')
plt.plot(data_00_all.TotalStrainVol, data_00_all.Pressure, 'C0')
plt.plot(data_09_all.TotalStrainVol, data_09_all.Pressure, 'C1')
plt.plot(data_18_all.TotalStrainVol, data_18_all.Pressure, 'C2')
plt.plot(data_27_all.TotalStrainVol, data_27_all.Pressure, 'C3')
plt.plot(data_36_all.TotalStrainVol, data_36_all.Pressure, 'C4')
plt.plot(data_45_all.TotalStrainVol, data_45_all.Pressure, 'C5')
plt.grid(True)
plt.axis([-0.1, 0.4, 0, 2e9])
plt.xlabel('Engineering volumetric strain', fontsize=16)
plt.ylabel('Pressure (Pa)', fontsize=16)
fig.savefig('Fox_DrySand_ElasticData.png')


# In[35]:


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(data_09_tension_all.TotalStrainVol, data_09_tension_all.Pressure, 'C7')
plt.plot(data_00_all.TotalStrainVol, data_00_all.Pressure, 'C0')
plt.plot(data_09_all.TotalStrainVol, data_09_all.Pressure, 'C1')
plt.plot(data_18_all.TotalStrainVol, data_18_all.Pressure, 'C2')
plt.plot(data_27_all.TotalStrainVol, data_27_all.Pressure, 'C3')
plt.plot(data_36_all.TotalStrainVol, data_36_all.Pressure, 'C4')
plt.plot(data_45_all.TotalStrainVol, data_45_all.Pressure, 'C5')
plt.grid(True)
plt.axis([-0.1, 0.2, 0, 1.0e8])
plt.xlabel('Engineering volumetric strain', fontsize=16)
plt.ylabel('Pressure (Pa)', fontsize=16)
fig.savefig('Fox_DrySand_ElasticData_Zoom.png')


# In[22]:


# Create data that can be written as JSON
elastic_data_json_09_t = {}
elastic_data_json_09_t["TotalStrainVol"] = data_09_tension_all.TotalStrainVol.values.tolist()
elastic_data_json_09_t["Pressure"] = data_09_tension_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_09_t)
elastic_data_json_00 = {}
elastic_data_json_00["TotalStrainVol"] = data_00_all.TotalStrainVol.values.tolist()
elastic_data_json_00["Pressure"] = data_00_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_00)
elastic_data_json_09 = {}
elastic_data_json_09["TotalStrainVol"] = data_09_all.TotalStrainVol.values.tolist()
elastic_data_json_09["Pressure"] = data_09_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_09)
elastic_data_json_18 = {}
elastic_data_json_18["TotalStrainVol"] = data_18_all.TotalStrainVol.values.tolist()
elastic_data_json_18["Pressure"] = data_18_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_18)
elastic_data_json_27 = {}
elastic_data_json_27["TotalStrainVol"] = data_27_all.TotalStrainVol.values.tolist()
elastic_data_json_27["Pressure"] = data_27_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_27)
elastic_data_json_36 = {}
elastic_data_json_36["TotalStrainVol"] = data_36_all.TotalStrainVol.values.tolist()
elastic_data_json_36["Pressure"] = data_36_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_36)
elastic_data_json_45 = {}
elastic_data_json_45["TotalStrainVol"] = data_45_all.TotalStrainVol.values.tolist()
elastic_data_json_45["Pressure"] = data_45_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_45)


# In[23]:


# Set up plastic strain list
plastic_strain_data = [-plastic_strain_09, 0, plastic_strain_09, plastic_strain_18, plastic_strain_27, plastic_strain_36, plastic_strain_45]
print(plastic_strain_data)


# In[24]:


# Set up elastic data list
vaango_elastic_data = {}
vaango_elastic_data["PlasticStrainVol"] = plastic_strain_data
vaango_elastic_data["Data"] = [elastic_data_json_09_t, elastic_data_json_00, elastic_data_json_09, elastic_data_json_18, elastic_data_json_27, elastic_data_json_36, elastic_data_json_45]
json.dumps(vaango_elastic_data)


# In[25]:


# Write JSON
meta_data_json = {}
meta_data_json["title"] = "Dry sand nonlinear elasticity data"

vaango_data = {}
vaango_data["Meta"] = meta_data_json
vaango_data["Data"] = vaango_elastic_data
json.dumps(vaango_data)

vaango_input_data = {}
vaango_input_data["Vaango_tabular_data"] = vaango_data
json.dumps(vaango_input_data, sort_keys=True)

with open('DrySand_ElasticData.json', 'w') as outputFile:
    json.dump(vaango_input_data, outputFile, sort_keys=True, indent=2)



# coding: utf-8

# In[35]:


# %load DrySand_CapData.py
#  Compression is positive
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import PolylineIntersection as pl
import imp
pl = imp.reload(pl)


# In[3]:


# Load the CSV data into Pandas dataframes
hydrostat = pd.read_csv("./DrySand_Hydrostat.csv", header=0, skiprows=7)
data_09 = pd.read_csv("./DrySand_LoadUnload_09.csv", header=0, skiprows=4)
data_18 = pd.read_csv("./DrySand_LoadUnload_18.csv", header=0, skiprows=4)
data_27 = pd.read_csv("./DrySand_LoadUnload_27.csv", header=0, skiprows=4)
data_36 = pd.read_csv("./DrySand_LoadUnload_36.csv", header=0, skiprows=4)
data_45 = pd.read_csv("./DrySand_LoadUnload_45.csv", header=0, skiprows=4)


# In[4]:


# Rename the columns of each dataframe
column_names = ["TotalStrainVol", "Pressure", "", "", "", ""]
hydrostat.columns = column_names
data_09.columns = column_names
data_18.columns = column_names
data_27.columns = column_names
data_36.columns = column_names
data_45.columns = column_names


# In[5]:


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


# In[21]:


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
p_min_00 = min(hydrostat.Pressure)
eps_min_00 = min(hydrostat.TotalStrainVol)

eps_max_09 = data_09.TotalStrainVol[p_max_index_09]
eps_max_18 = data_18.TotalStrainVol[p_max_index_18]
eps_max_27 = data_27.TotalStrainVol[p_max_index_27]
eps_max_36 = data_36.TotalStrainVol[p_max_index_36]
eps_max_45 = data_45.TotalStrainVol[p_max_index_45]
print(eps_max_45, ",", p_max_45, ",", p_max_index_45, ",", p_max_00, ",", p_max_index_00)
print(eps_min_00, p_min_00)


# In[41]:


# Find intersections of unloading start points with hydrostatic curve
hydrostat_poly = list(hydrostat[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
strain_axis_09 = ((eps_max_09, 0),(eps_max_09, 1.0e11))
[flag_09, pt_09, index_09] = pl.poly_intersection(strain_axis_09, hydrostat_poly)
strain_axis_18 = ((eps_max_18, 0),(eps_max_18, 1.0e11))
[flag_18, pt_18, index_18] = pl.poly_intersection(strain_axis_18, hydrostat_poly)
strain_axis_27 = ((eps_max_27, 0),(eps_max_27, 1.0e11))
[flag_27, pt_27, index_27] = pl.poly_intersection(strain_axis_27, hydrostat_poly)
strain_axis_36 = ((eps_max_36, 0),(eps_max_36, 1.0e11))
[flag_36, pt_36, index_36] = pl.poly_intersection(strain_axis_36, hydrostat_poly)
strain_axis_45 = ((eps_max_45, 0),(eps_max_45, 1.0e11))
[flag_45, pt_45, index_45] = pl.poly_intersection(strain_axis_45, hydrostat_poly)
pressure_09 = pt_09[1]
pressure_18 = pt_18[1]
pressure_27 = pt_27[1]
pressure_36 = pt_36[1]
pressure_45 = pt_45[1]
print(p_max_09, pressure_09, index_09, flag_09)
print(p_max_18, pressure_18, index_18, flag_18)
print(p_max_27, pressure_27, index_27, flag_27)
print(p_max_36, pressure_36, index_36, flag_36)
print(p_max_45, pressure_45, index_45, flag_45)


# In[40]:


# Visual check of data set
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(hydrostat.TotalStrainVol, hydrostat.Pressure, 'C7')
plt.plot(eps_max_09, p_max_09, 'C1x')
plt.plot(eps_max_18, p_max_18, 'C2x')
plt.plot(eps_max_27, p_max_27, 'C3x')
plt.plot(eps_max_36, p_max_36, 'C4x')
plt.plot(eps_max_45, p_max_45, 'C5x')
plt.plot(eps_max_09, pressure_09, 'C1o')
plt.plot(eps_max_18, pressure_18, 'C2o')
plt.plot(eps_max_27, pressure_27, 'C3o')
plt.plot(eps_max_36, pressure_36, 'C4o')
plt.plot(eps_max_45, pressure_45, 'C5o')
plt.grid(True)
#plt.axis([-0.1, 0.4, 0, 0.07e9])
plt.xlabel('Engineering volumetric strain', fontsize=16)
plt.ylabel('Pressure (Pa)', fontsize=16)
#fig.savefig('Fox_DrySand_ElasticData.png')


# In[7]:


# Create separate dataframes for the unload data
data_09_unload = data_09[p_max_index_09:].copy()
data_18_unload = data_18[p_max_index_18:].copy()
data_27_unload = data_27[p_max_index_27:].copy()
data_36_unload = data_36[p_max_index_36:].copy()
data_45_unload = data_45[p_max_index_45:].copy()


# In[8]:


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


# In[62]:


# Extract segments for interpolation
segment_00_09 = hydrostat[0:(index_09+1)]
segment_09_18 = hydrostat[index_09:(index_18+1)]
segment_18_27 = hydrostat[index_18:(index_27+1)]
segment_27_36 = hydrostat[index_27:(index_36+1)]
segment_36_45 = hydrostat[index_36:(index_45+1)]
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(hydrostat.TotalStrainVol, hydrostat.Pressure, 'C7')
plt.plot(segment_00_09.TotalStrainVol, segment_00_09.Pressure, 'C6')
plt.plot(segment_09_18.TotalStrainVol, segment_09_18.Pressure, 'C1')
plt.plot(segment_18_27.TotalStrainVol, segment_18_27.Pressure, 'C2')
plt.plot(segment_27_36.TotalStrainVol, segment_27_36.Pressure, 'C3')
plt.plot(segment_36_45.TotalStrainVol, segment_36_45.Pressure, 'C4')
plt.grid(True)
#plt.axis([-0.1, 0.4, 0, 0.07e9])
plt.xlabel('Engineering volumetric strain', fontsize=16)
plt.ylabel('Pressure (Pa)', fontsize=16)


# In[68]:


# Interpolate plastic strains to points on hydrostatic curve
def compute_t(ev, ev_min, ev_max):
    return (ev - ev_min)/(ev_max - ev_min)
def compute_ep(t, ep_min, ep_max):
    return (1 - t)*ep_min + t*ep_max;
def interp_ep(ev, ev_min, ev_max, ep_min, ep_max):
    t = compute_t(ev, ev_min, ev_max)
    return compute_ep(t, ep_min, ep_max)
def interp_seg(segment, ep_min, ep_max):
    ev_min = segment.TotalStrainVol.values[0]
    ev_max = segment.TotalStrainVol.values[-1]
    #print(ev_min, ev_max, ep_min, ep_max)
    ep = list(map(lambda ev: interp_ep(ev, ev_min, ev_max, ep_min, ep_max), segment.TotalStrainVol))
    return ep
def exterp_seg(segment_ext, segment, ep_min, ep_max):
    ev_min = segment.TotalStrainVol.values[0]
    ev_max = segment.TotalStrainVol.values[-1]
    #print(ev_min, ev_max, ep_min, ep_max)
    ep = list(map(lambda ev: interp_ep(ev, ev_min, ev_max, ep_min, ep_max), segment_ext.TotalStrainVol))
    return ep

ep_00_09 = exterp_seg(segment_00_09, segment_09_18, plastic_strain_09, plastic_strain_18)
ep_09_18 = interp_seg(segment_09_18, plastic_strain_09, plastic_strain_18)
ep_18_27 = interp_seg(segment_18_27, plastic_strain_18, plastic_strain_27)
ep_27_36 = interp_seg(segment_27_36, plastic_strain_27, plastic_strain_36)
ep_36_45 = interp_seg(segment_36_45, plastic_strain_36, plastic_strain_45)


# In[76]:


# Check interpolated data
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(ep_00_09, segment_00_09.Pressure, 'C6')
plt.plot(ep_09_18, segment_09_18.Pressure, 'C1')
plt.plot(ep_18_27, segment_18_27.Pressure, 'C2')
plt.plot(ep_27_36, segment_27_36.Pressure, 'C3')
plt.plot(ep_36_45, segment_36_45.Pressure, 'C4')
plt.grid(True)
#plt.axis([-0.1, 0.4, 0, 0.07e9])
#plt.axis([0.27, 0.35, 0.9e9, 2.0e9])
plt.xlabel('Plastic volumetric strain', fontsize=16)
plt.ylabel('Pressure (Pa)', fontsize=16)


# In[94]:


# Remove negative plastic strain data
seg_00_09 = pd.DataFrame({"PlasticStrainVol" : ep_00_09, "Pressure" : segment_00_09.Pressure})
seg_00_09_poly = list(seg_00_09.apply(tuple, axis=1))
ep_axis_00 = ((0, 0),(0, 1.0e11))
[flag_00, pt_00, index_00] = pl.poly_intersection(ep_axis_00, seg_00_09_poly)
#print(flag_00, pt_00, index_00)
#print(seg_00_09_poly[index_00])
seg_00_09.iloc[index_00] = dict(PlasticStrainVol=0, Pressure=pt_00[1])
seg_00_09_upd = seg_00_09[index_00:-1]


# In[96]:


# Join the segments
seg_09_18 = pd.DataFrame({"PlasticStrainVol" : ep_09_18, "Pressure" : segment_09_18.Pressure})
seg_18_27 = pd.DataFrame({"PlasticStrainVol" : ep_18_27, "Pressure" : segment_18_27.Pressure})
seg_27_36 = pd.DataFrame({"PlasticStrainVol" : ep_27_36, "Pressure" : segment_27_36.Pressure})
seg_36_45 = pd.DataFrame({"PlasticStrainVol" : ep_36_45, "Pressure" : segment_36_45.Pressure})
data_p_vs_ep = pd.concat([seg_00_09_upd, seg_09_18, seg_18_27, seg_27_36, seg_36_45], ignore_index=True)


# In[99]:


# Remove duplicates (if any)
data_all = data_p_vs_ep.copy().drop_duplicates()
print(data_p_vs_ep.size, data_all.size)
print(data_p_vs_ep)


# In[101]:


# Check data
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(data_all.PlasticStrainVol, data_all.Pressure, 'C6')
plt.grid(True)
#plt.axis([-0.1, 0.4, 0, 0.07e9])
#plt.axis([0.27, 0.35, 0.9e9, 2.0e9])
plt.xlabel('Plastic volumetric strain', fontsize=16)
plt.ylabel('Pressure (Pa)', fontsize=16)# In[18]:


# In[103]:


# Create data that can be written as JSON
data_json = {}
data_json["PlasticStrainVol"] = data_all.PlasticStrainVol.values.tolist()
data_json["Pressure"] = data_all.Pressure.values.tolist()
#json.dumps(data_json)


# In[104]:


# Write JSON
meta_data_json = {}
meta_data_json["title"] = "Dry sand hydrostatic strength cap data"

vaango_data = {}
vaango_data["Meta"] = meta_data_json
vaango_data["Data"] = data_json
#json.dumps(vaango_data)


# In[105]:


vaango_input_data = {}
vaango_input_data["Vaango_tabular_data"] = vaango_data
json.dumps(vaango_input_data, sort_keys=True)

with open('DrySand_HydrostaticCapData.json', 'w') as outputFile:
    json.dump(vaango_input_data, outputFile, sort_keys=True, indent=2)


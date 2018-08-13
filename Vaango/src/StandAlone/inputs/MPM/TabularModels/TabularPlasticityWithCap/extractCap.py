# %load DrySand_CapData.py
#  Compression is positive
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import PolylineIntersection as pl
from scipy import interpolate

def computeBulkModulus(data_load, eps_max):

  # Get the minimum strains
  eps_min = min(data_load['TotalStrainVol'])

  # Set up range of strains
  deps = data_load['TotalStrainVol'].values[-1] - \
         data_load['TotalStrainVol'].values[-2];
  neps = round((eps_max - eps_min)/(30*deps))
  eps_range = np.linspace(eps_min, eps_max, neps)

  # Create linear interpolation
  interp = interpolate.interp1d(data_load['TotalStrainVol'], 
                                data_load['Pressure'],
                                kind = 'cubic')

  # Compute interpolated values
  p_range = interp(eps_range)

  # Compute bulk modulus 
  K = np.diff(p_range)/np.diff(eps_range)

  # Create linear interpolation
  interp1 = interpolate.interp1d(eps_range[1:], K, kind='cubic',
                                 fill_value = 'extrapolate')

  # Compute interpolated bulk modulus 
  K_range = interp1(eps_range[1:])

  return eps_range, p_range, K_range


# Load the CSV data into Pandas dataframes
hydrostat = pd.read_csv("./DrySand_Hydrostat.csv", header=0, skiprows=7)
data_09 = pd.read_csv("./DrySand_LoadUnload_09.csv", header=0, skiprows=4)
data_18 = pd.read_csv("./DrySand_LoadUnload_18.csv", header=0, skiprows=4)
data_27 = pd.read_csv("./DrySand_LoadUnload_27.csv", header=0, skiprows=4)
data_36 = pd.read_csv("./DrySand_LoadUnload_36.csv", header=0, skiprows=4)
data_45 = pd.read_csv("./DrySand_LoadUnload_45.csv", header=0, skiprows=4)

# Rename the columns of each dataframe
column_names = ["TotalStrainVol", "Pressure", "", "", "", ""]
hydrostat.columns = column_names
data_09.columns = column_names
data_18.columns = column_names
data_27.columns = column_names
data_36.columns = column_names
data_45.columns = column_names

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

# Create separate dataframes for the unload data
data_09_unload = data_09[p_max_index_09:].copy().drop_duplicates()
data_18_unload = data_18[p_max_index_18:].copy().drop_duplicates()
data_27_unload = data_27[p_max_index_27:].copy().drop_duplicates()
data_36_unload = data_36[p_max_index_36:].copy().drop_duplicates()
data_45_unload = data_45[p_max_index_45:].copy().drop_duplicates()

# Reverse the order of the unload data to create elastic loading curves
data_09_load = data_09_unload.sort_index(ascending=False)
data_18_load = data_18_unload.sort_index(ascending=False)
data_27_load = data_27_unload.sort_index(ascending=False)
data_36_load = data_36_unload.sort_index(ascending=False)
data_45_load = data_45_unload.sort_index(ascending=False)

# Find intersections of unloading curves with hydrostatic curve
hydrostat_poly = list(hydrostat[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
data_09_poly = list(data_09_load[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
data_18_poly = list(data_18_load[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
data_27_poly = list(data_27_load[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
data_36_poly = list(data_36_load[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
data_45_poly = list(data_45_load[['TotalStrainVol', 'Pressure']].apply(tuple, axis=1))
[flag_09, pt_09, index_09] = pl.poly_intersection(data_09_poly, hydrostat_poly)
[flag_18, pt_18, index_18] = pl.poly_intersection(data_18_poly, hydrostat_poly)
[flag_27, pt_27, index_27] = pl.poly_intersection(data_27_poly, hydrostat_poly)
[flag_36, pt_36, index_36] = pl.poly_intersection(data_36_poly, hydrostat_poly)
[flag_45, pt_45, index_45] = pl.poly_intersection(data_45_poly, hydrostat_poly)
#print(pt_09, index_09, flag_09)
#print(pt_18, index_18, flag_18)
#print(pt_27, index_27, flag_27)
#print(pt_36, index_36, flag_36)
#print(pt_45, index_45, flag_45)

# Save the intersection points
eps_max_09 = pt_09[0]
eps_max_18 = pt_18[0]
eps_max_27 = pt_27[0]
eps_max_36 = pt_36[0]
eps_max_45 = pt_45[0]
pressure_09 = pt_09[1]
pressure_18 = pt_18[1]
pressure_27 = pt_27[1]
pressure_36 = pt_36[1]
pressure_45 = pt_45[1]
t_hydro_09 = pt_09[3]
t_hydro_18 = pt_18[3]
t_hydro_27 = pt_27[3]
t_hydro_36 = pt_36[3]
t_hydro_45 = pt_45[3]

# Compute initial bulk modulus
K_00 = (hydrostat['Pressure'].values[1] - hydrostat['Pressure'].values[0])/ \
       (hydrostat['TotalStrainVol'].values[1] - hydrostat['TotalStrainVol'].values[0])
eps_00 = [0, eps_max_45]
p_00 = [0, K_00*eps_max_45]
data_00_poly = [(eps_00[0], p_00[0]), (eps_00[1], p_00[1])]
hydrostat_poly0 = list(hydrostat[['TotalStrainVol', 'Pressure']][1:].apply(tuple, axis=1))
[flag_00, pt_00, index_00] = pl.poly_intersection(data_00_poly, hydrostat_poly0)
print(pt_00, index_00, flag_00)
index_00 = index_00 + 1
eps_max_00 = pt_00[0]
pressure_00 = pt_00[1]
t_hydro_00 = pt_00[3]
print(eps_max_00, pressure_00, t_hydro_00)
eps_00 = [0, eps_max_00]
p_00 = [0, K_00*eps_max_00]

# Compute the bulk modulus
eps_09, p_09, K_09 = computeBulkModulus(data_09_load, eps_max_09)
eps_18, p_18, K_18 = computeBulkModulus(data_18_load, eps_max_18)
eps_27, p_27, K_27 = computeBulkModulus(data_27_load, eps_max_27)
eps_36, p_36, K_36 = computeBulkModulus(data_36_load, eps_max_36)
eps_45, p_45, K_45 = computeBulkModulus(data_45_load, eps_max_45)

# Find plastic strains by intersecting the unload data with the pressure axis
pressure_axis = ((-1, 0),(1, 0))
line_09_load = (data_09_poly[0], data_09_poly[1])
eps_p_const_09 = pl.line_intersection(pressure_axis, line_09_load)[0]
line_18_load = (data_18_poly[0], data_18_poly[1])
eps_p_const_18 = pl.line_intersection(pressure_axis, line_18_load)[0]
line_27_load = (data_27_poly[0], data_27_poly[1])
eps_p_const_27 = pl.line_intersection(pressure_axis, line_27_load)[0]
line_36_load = (data_36_poly[0], data_36_poly[1])
eps_p_const_36 = pl.line_intersection(pressure_axis, line_36_load)[0]
line_45_load = (data_45_poly[0], data_45_poly[1])
eps_p_const_45 = pl.line_intersection(pressure_axis, line_45_load)[0]

# Compute eps_p = eps - pressure/bulk modulus
eps_p_00 = eps_00[1:] - p_00[1:]/K_00
eps_p_09 = eps_09[1:] - p_09[1:]/K_09
eps_p_18 = eps_18[1:] - p_18[1:]/K_18
eps_p_27 = eps_27[1:] - p_27[1:]/K_27
eps_p_36 = eps_36[1:] - p_36[1:]/K_36
eps_p_45 = eps_45[1:] - p_45[1:]/K_45

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
plt.plot(eps_00, p_00, 'C0x')
plt.plot(eps_09, p_09, 'C1')
plt.plot(eps_18, p_18, 'C2')
plt.plot(eps_27, p_27, 'C3')
plt.plot(eps_36, p_36, 'C4')
plt.plot(eps_45, p_45, 'C5')
plt.grid(True)
plt.xlabel('Engineering volumetric strain', fontsize=16)
plt.ylabel('Pressure (Pa)', fontsize=16)
#plt.show()

# Visual check of data set
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(eps_00[1:], K_00, 'C0x')
plt.plot(eps_09[1:], K_09, 'C1')
plt.plot(eps_18[1:], K_18, 'C2')
plt.plot(eps_27[1:], K_27, 'C3')
plt.plot(eps_36[1:], K_36, 'C4')
plt.plot(eps_45[1:], K_45, 'C5')
plt.grid(True)
plt.xlabel('Engineering volumetric strain', fontsize=16)
plt.ylabel('Bulk Modulus (Pa)', fontsize=16)
#plt.show()

# Visual check of data set
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(eps_00[1:], eps_p_00, 'C0x')
plt.plot(eps_09[1:], eps_p_09, 'C1')
plt.plot(eps_18[1:], eps_p_18, 'C2')
plt.plot(eps_27[1:], eps_p_27, 'C3')
plt.plot(eps_36[1:], eps_p_36, 'C4')
plt.plot(eps_45[1:], eps_p_45, 'C5')
plt.plot(0, 0, 'C0x')
plt.plot(eps_p_const_09, eps_p_const_09, 'C1x')
plt.plot(eps_p_const_18, eps_p_const_18, 'C2x')
plt.plot(eps_p_const_27, eps_p_const_27, 'C3x')
plt.plot(eps_p_const_36, eps_p_const_36, 'C4x')
plt.plot(eps_p_const_45, eps_p_const_45, 'C5x')
plt.grid(True)
plt.xlabel('Engineering volumetric strain', fontsize=16)
plt.ylabel('Plastic volumetric strain', fontsize=16)
#plt.show()

# Compute p and eps_v increments along hydrostat
p_diff = np.diff(hydrostat['Pressure'])
eps_diff = np.diff(hydrostat['TotalStrainVol'])

# Compute K values along hydrostat
K_00_max = K_00
K_09_max = max(K_09)
K_18_max = max(K_18)
K_27_max = max(K_27)
K_36_max = max(K_36)
K_45_max = max(K_45)

# Visual check of data set
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(0.0, K_00_max, 'C0x')
plt.plot(eps_max_09, K_09_max, 'C1x')
plt.plot(eps_max_18, K_18_max, 'C2x')
plt.plot(eps_max_27, K_27_max, 'C3x')
plt.plot(eps_max_36, K_36_max, 'C4x')
plt.plot(eps_max_45, K_45_max, 'C5x')
plt.grid(True)
plt.xlabel('Total volumetric strain', fontsize=16)
plt.ylabel('Bulk modulus max (Pa)', fontsize=16)

def interpData(segment, eps_min, eps_max):

  # Set up range of strains
  deps = segment['TotalStrainVol'].values[-1] - \
         segment['TotalStrainVol'].values[-2];
  neps = round((eps_max - eps_min)/deps)
  eps_range = np.linspace(eps_min, eps_max, neps)

  # Create linear interpolation
  interp = interpolate.interp1d(segment['TotalStrainVol'], 
                                segment['Pressure'],
                                kind = 'cubic')

  # Compute interpolated values
  p_range = interp(eps_range)

  return eps_range, p_range

# Extract hydrostat segments
seg_00_09 = hydrostat[index_00:index_09+2]
eps_00_09, p_00_09 = interpData(seg_00_09, eps_max_00, eps_max_09)
t_00_09 = (eps_00_09 - eps_max_00)/(eps_max_09 - eps_max_00)
K_00_09 = (1 - t_00_09)*K_00_max + t_00_09*K_09_max
plt.plot(eps_00_09, K_00_09, 'C1')

seg_09_18 = hydrostat[index_09-1:index_18+2]
eps_09_18, p_09_18 = interpData(seg_09_18, eps_max_09, eps_max_18)
t_09_18 = (eps_09_18 - eps_max_09)/(eps_max_18 - eps_max_09)
K_09_18 = (1 - t_09_18)*K_09_max + t_09_18*K_18_max
plt.plot(eps_09_18, K_09_18, 'C2')

seg_18_27 = hydrostat[index_18-1:index_27+2]
eps_18_27, p_18_27 = interpData(seg_18_27, eps_max_18, eps_max_27)
t_18_27 = (eps_18_27 - eps_max_18)/(eps_max_27 - eps_max_18)
K_18_27 = (1 - t_18_27)*K_18_max + t_18_27*K_27_max
plt.plot(eps_18_27, K_18_27, 'C3')

seg_27_36 = hydrostat[index_27-1:index_36+2]
eps_27_36, p_27_36 = interpData(seg_27_36, eps_max_27, eps_max_36)
t_27_36 = (eps_27_36 - eps_max_27)/(eps_max_36- eps_max_27)
K_27_36 = (1 - t_27_36)*K_27_max + t_27_36*K_36_max
plt.plot(eps_27_36, K_27_36, 'C4')

seg_36_45 = hydrostat[index_36-1:index_45+2]
eps_36_45, p_36_45 = interpData(seg_36_45, eps_max_36, eps_max_45)
t_36_45 = (eps_36_45 - eps_max_36)/(eps_max_45- eps_max_36)
K_36_45 = (1 - t_36_45)*K_36_max + t_36_45*K_45_max
plt.plot(eps_36_45, K_36_45, 'C5')

# Compute the plastic strains along the hydrostat segments
eps_p_00_09 = eps_00_09 - p_00_09/K_00_09
eps_p_09_18 = eps_09_18 - p_09_18/K_09_18
eps_p_18_27 = eps_18_27 - p_18_27/K_18_27
eps_p_27_36 = eps_27_36 - p_27_36/K_27_36
eps_p_36_45 = eps_36_45 - p_36_45/K_36_45

# Visual check of data set
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(eps_p_00_09, p_00_09, 'C1')
plt.plot(eps_p_09_18, p_09_18, 'C2')
plt.plot(eps_p_18_27, p_18_27, 'C3')
plt.plot(eps_p_27_36, p_27_36, 'C4')
plt.plot(eps_p_36_45, p_36_45, 'C5')
plt.grid(True)
plt.xlabel('Plastic volumetric strain', fontsize=16)
plt.ylabel('Pressure (Pa)', fontsize=16)

plt.show()

# Join the segments
df_00_09 = pd.DataFrame({"PlasticStrainVol" : eps_p_00_09, "Pressure" : p_00_09})
df_09_18 = pd.DataFrame({"PlasticStrainVol" : eps_p_09_18, "Pressure" : p_09_18})
df_18_27 = pd.DataFrame({"PlasticStrainVol" : eps_p_18_27, "Pressure" : p_18_27})
df_27_36 = pd.DataFrame({"PlasticStrainVol" : eps_p_27_36, "Pressure" : p_27_36})
df_36_45 = pd.DataFrame({"PlasticStrainVol" : eps_p_36_45, "Pressure" : p_36_45})
data_p_vs_ep = pd.concat([df_00_09, df_09_18, df_18_27, df_27_36, df_36_45], ignore_index=True)

# Remove duplicates (if any)
data_all = data_p_vs_ep.copy().drop_duplicates()
print(data_p_vs_ep.size, data_all.size)
print(data_p_vs_ep)

# Check data
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(data_all.PlasticStrainVol, data_all.Pressure, 'C6')
plt.grid(True)
#plt.axis([-0.1, 0.4, 0, 0.07e9])
#plt.axis([0.27, 0.35, 0.9e9, 2.0e9])
plt.xlabel('Plastic volumetric strain', fontsize=16)
plt.ylabel('Pressure (Pa)', fontsize=16)# In[18]:

# Create data that can be written as JSON
data_json = {}
data_json["PlasticStrainVol"] = data_all.PlasticStrainVol.values.tolist()
data_json["Pressure"] = data_all.Pressure.values.tolist()
#json.dumps(data_json)

# Write JSON
meta_data_json = {}
meta_data_json["title"] = "Dry sand hydrostatic strength cap data"

vaango_data = {}
vaango_data["Meta"] = meta_data_json
vaango_data["Data"] = data_json
#json.dumps(vaango_data)

vaango_input_data = {}
vaango_input_data["Vaango_tabular_data"] = vaango_data
json.dumps(vaango_input_data, sort_keys=True)

with open('DrySand_HydrostaticCapData.json', 'w') as outputFile:
    json.dump(vaango_input_data, outputFile, sort_keys=True, indent=2)


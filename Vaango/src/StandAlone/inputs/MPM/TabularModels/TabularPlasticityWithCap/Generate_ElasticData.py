# %load Generate_ElasticData.py
#  Compression is positive
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json

# function to integrate the bulk modulus
def computeK(eps_p_v):

  k0 = 5.0e9
  k1 = 10.0e9
  eps_p_ref = 1
  K = k0 + k1 * np.exp(eps_p_v/eps_p_ref)
  return K

def integrateK(eps_p_v):

  deps_v = 0.0001
  K = computeK(eps_p_v)
  eps_v_list = np.linspace(0, 1.0, 10)
  p_list = []
  for ii, eps_v in enumerate(eps_v_list):
     p_list.append(0.0)
     deps_v_list = [-deps_v] * int(eps_v/deps_v)
     for deps in deps_v_list:
       p_list[ii] += K*deps

  # Convert to compression +ve
  eps_bar_v_list = list(map(lambda ep : ep - eps_p_v, eps_v_list))
  p_bar_list = list(map(lambda p : -p, p_list))
  return eps_bar_v_list, p_bar_list

# Compute the p vs eps_v data at variout plastic strains
eps_v_bar, p_bar = integrateK(0.0)
data_00_load = pd.DataFrame(data = {'TotalStrainVol': eps_v_bar, 'Pressure': p_bar})
eps_v_bar, p_bar = integrateK(-0.2)
data_10_load = pd.DataFrame(data = {'TotalStrainVol': eps_v_bar, 'Pressure': p_bar})
eps_v_bar, p_bar = integrateK(-0.4)
data_20_load = pd.DataFrame(data = {'TotalStrainVol': eps_v_bar, 'Pressure': p_bar})
eps_v_bar, p_bar = integrateK(-0.6)
data_30_load = pd.DataFrame(data = {'TotalStrainVol': eps_v_bar, 'Pressure': p_bar})
eps_v_bar, p_bar = integrateK(-0.8)
data_40_load = pd.DataFrame(data = {'TotalStrainVol': eps_v_bar, 'Pressure': p_bar})
eps_v_bar, p_bar = integrateK(-1.0)
data_50_load = pd.DataFrame(data = {'TotalStrainVol': eps_v_bar, 'Pressure': p_bar})

# Simplify the loading dataframes and remove duplicates (if any)
data_00_all = data_00_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_10_all = data_10_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_20_all = data_20_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_30_all = data_30_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_40_all = data_40_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()
data_50_all = data_50_load[['TotalStrainVol', 'Pressure']].copy().drop_duplicates()

# Visual check of data set
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.plot(data_00_all.TotalStrainVol, data_00_all.Pressure, 'C0')
plt.plot(data_10_all.TotalStrainVol, data_10_all.Pressure, 'C1')
plt.plot(data_20_all.TotalStrainVol, data_20_all.Pressure, 'C2')
plt.plot(data_30_all.TotalStrainVol, data_30_all.Pressure, 'C3')
plt.plot(data_40_all.TotalStrainVol, data_40_all.Pressure, 'C4')
plt.plot(data_50_all.TotalStrainVol, data_50_all.Pressure, 'C5')
plt.grid(True)
#plt.axis([-0.1, 0.4, 0, 2e9])
plt.xlabel('Engineering volumetric strain', fontsize=16)
plt.ylabel('Pressure (Pa)', fontsize=16)
#fig.savefig('Fox_DrySand_ElasticData.png')
plt.show()


# Create data that can be written as JSON
elastic_data_json_00 = {}
elastic_data_json_00["TotalStrainVol"] = data_00_all.TotalStrainVol.values.tolist()
elastic_data_json_00["Pressure"] = data_00_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_00)
elastic_data_json_10 = {}
elastic_data_json_10["TotalStrainVol"] = data_10_all.TotalStrainVol.values.tolist()
elastic_data_json_10["Pressure"] = data_10_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_10)
elastic_data_json_20 = {}
elastic_data_json_20["TotalStrainVol"] = data_20_all.TotalStrainVol.values.tolist()
elastic_data_json_20["Pressure"] = data_20_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_18)
elastic_data_json_30 = {}
elastic_data_json_30["TotalStrainVol"] = data_30_all.TotalStrainVol.values.tolist()
elastic_data_json_30["Pressure"] = data_30_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_30)
elastic_data_json_40 = {}
elastic_data_json_40["TotalStrainVol"] = data_40_all.TotalStrainVol.values.tolist()
elastic_data_json_40["Pressure"] = data_40_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_40)
elastic_data_json_50 = {}
elastic_data_json_50["TotalStrainVol"] = data_50_all.TotalStrainVol.values.tolist()
elastic_data_json_50["Pressure"] = data_50_all.Pressure.values.tolist()
#json.dumps(elastic_data_json_50)


# Set up plastic strain list
plastic_strain_data = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
print(plastic_strain_data)


# Set up elastic data list
vaango_elastic_data = {}
vaango_elastic_data["PlasticStrainVol"] = plastic_strain_data
vaango_elastic_data["Data"] = [elastic_data_json_00, elastic_data_json_10, elastic_data_json_20, elastic_data_json_30, elastic_data_json_40, elastic_data_json_50]
json.dumps(vaango_elastic_data)

# Write JSON
meta_data_json = {}
meta_data_json["title"] = "Homel plastic coupled linear elasticity data"

vaango_data = {}
vaango_data["Meta"] = meta_data_json
vaango_data["Data"] = vaango_elastic_data
json.dumps(vaango_data)

vaango_input_data = {}
vaango_input_data["Vaango_tabular_data"] = vaango_data
json.dumps(vaango_input_data, sort_keys=True)

with open('TabularCapTest_11_Elastic.json', 'w') as outputFile:
    json.dump(vaango_input_data, outputFile, sort_keys=True, indent=2)


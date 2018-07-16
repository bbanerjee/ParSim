
# coding: utf-8

# In[22]:


# %load createDrySand_YieldData_JSON.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import PolylineIntersection as pl

data = pd.read_csv("./DrySand_YieldData.csv", header=0, skiprows=5)
data.columns = ["Pressure", "SqrtJ2"]


# In[25]:


# Convert MPa to Pa
data.Pressure *= 1.0e6
data.SqrtJ2 *= 1.0e6


# In[29]:


# Add tension cap point
pressure_axis = ((-1.0e6, 0),(1.0e6, 0))
polyline_data = list(data.apply(tuple, axis=1))
sqrtJ2_line = (polyline_data[0], polyline_data[1])
tension_cap_point = pl.line_intersection(pressure_axis, sqrtJ2_line)
print(tension_cap_point)

pressure_list = data.Pressure.values.tolist()
pressure_list.insert(0, tension_cap_point[0])
sqrtJ2_list = data.SqrtJ2.values.tolist()
sqrtJ2_list.insert(0, 0)


# In[30]:


yield_data_json = {}
yield_data_json["Pressure"] = pressure_list
yield_data_json["SqrtJ2"] = sqrtJ2_list
json.dumps(yield_data_json)


# In[31]:


meta_data_json = {}
meta_data_json["title"] = "Dry sand yield function data"


# In[32]:


vaango_data = {}
vaango_data["Meta"] = meta_data_json
vaango_data["Data"] = yield_data_json
json.dumps(vaango_data)


# In[33]:


vaango_input_data = {}
vaango_input_data["Vaango_tabular_data"] = vaango_data
json.dumps(vaango_input_data, sort_keys=True)


# In[34]:


with open('DrySand_YieldData.json', 'w') as outputFile:
    json.dump(vaango_input_data, outputFile, sort_keys=True, indent=2)


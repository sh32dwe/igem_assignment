#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra


# In[2]:


from cobra.io import load_model


# In[10]:


model = cobra.io.read_sbml_model("C:/Users/Dr Tauseef Ahmad/Downloads/e_coli_core.xml.gz")


# In[8]:


import numpy as np
import matplotlib.pyplot as plt


# In[11]:


initial_glucose_concentration = 10


# In[12]:


model.reactions.EX_glc__D_e.lower_bound = -initial_glucose_concentration


# In[21]:


total_time = 10  # Total simulation time in hours
time_step = 1/3600 
time_points = np.arange(0, total_time + time_step, time_step)


# In[22]:


growth_rates = []


# In[23]:


for t in time_points:
    solution = model.optimize()
    growth_rate = solution.objective_value
    growth_rates.append(growth_rate)


# In[24]:


plt.plot(time_points, growth_rates)
plt.xlabel('Time ()')
plt.ylabel('Growth rate (mmol/gDW/h)')
plt.title('Growth of E. coli over time')
plt.grid(True)
plt.show()


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra


# In[3]:


import numpy as np


# In[4]:


import matplotlib.pyplot as plt


# In[5]:


model = cobra.io.read_sbml_model("C:/Users/Dr Tauseef Ahmad/Downloads/e_coli_core.xml.gz")


# In[34]:


initial_glucose_concentration = 2  # Initial glucose concentration in millimoles per liter
total_time = 10  # Total simulation time in hours
time_step = 1/3600 


# In[35]:


time_points = np.arange(0, total_time + time_step, time_step)
glucose_concentrations = [initial_glucose_concentration]
growth_rates = []


# In[36]:


for t in range(1, len(time_points)):
    # Set the glucose uptake rate constraint based on the current glucose concentration
    model.reactions.EX_glc__D_e.lower_bound = -glucose_concentrations[t - 1]
    
    # Optimize the model to obtain flux distributions
    solution = model.optimize()
    
    # Extract the growth rate
    growth_rate = solution.objective_value
    
    # Decrement glucose concentration based on biomass production rate
    glucose_concentrations.append(glucose_concentrations[t-1] - growth_rate * time_step)

    # Append growth rate
    growth_rates.append(growth_rate)


# In[37]:


plt.plot(time_points[:-1], growth_rates, color='green')
plt.xlabel('Time (hours)')
plt.ylabel('Growth Rate')
plt.title('Growth Rate over Time')
plt.grid(True)
plt.show()


# In[39]:


plt.plot(time_points, glucose_concentrations, color='blue')
plt.xlabel('Time (hours)')
plt.ylabel('Glucose Concentration (millimoles per liter)')
plt.title('Glucose Concentration over Time')
plt.grid(True)
plt.show()


# In[43]:


plt.plot(glucose_concentrations[:-1], growth_rates, color='blue')
plt.ylabel('Growth rate')
plt.xlabel('Glucose Concentration (millimoles per liter)')
plt.title('Glucose Concentration vs Growth rate')
plt.grid(True)
plt.show()


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
from cobra.io import load_model


# In[2]:


#EColi=iJO1366
#Bacillus subtilis=iYO844
model_coli= load_model('iJO1366')


# In[3]:


print(model_coli.objective)


# In[11]:


biomass_reaction = model_coli.reactions.get_by_id('BIOMASS_Ec_iJO1366_core_53p95M')


# In[13]:


precursors = {}          #Initialising a dictionary to store precursor metabolites
for metabolite, coefficient in biomass_reaction.metabolites.items():
    if coefficient < 0:  # reactants have negative coefficients
        precursors[metabolite.id] = metabolite.name


# In[14]:


precursor_producing_reactions = {metabolite_id: [] for metabolite_id in precursors}


# In[17]:


for metabolite_id, (metabolite_name) in precursors.items():
    metabolite = model_coli.metabolites.get_by_id(metabolite_id)
    for reaction in model_coli.reactions:
        if metabolite in reaction.products:
            precursor_producing_reactions[metabolite_id].append(reaction.id)


# In[35]:


essential_genes = {metabolite_id: {} for metabolite_id in precursors}


# In[40]:


for metabolite_id, (metabolite_name) in precursors.items():
    metabolite = model_coli.metabolites.get_by_id(metabolite_id)
    for reaction_id in precursor_producing_reactions[metabolite_id]:
        reaction = model_coli.reactions.get_by_id(reaction_id)
        model_coli.objective = reaction
        original_flux = model_coli.slim_optimize()
        for gene in reaction.genes:
            with model_coli:
                gene.knock_out()
                solution = model_coli.optimize()
                new_flux = solution.objective_value
                if reaction_id not in essential_genes[metabolite_id]:
                       essential_genes[metabolite_id][reaction_id] = []
                if new_flux < 0.1* original_flux:
                        essential_genes[metabolite_id][reaction_id].append(gene.id)


# In[44]:


all_essential_genes = {}
for metabolite_id, (metabolite_name) in precursors.items():
    essential_gene_dict = essential_genes.get(metabolite_id, {})  # Get essential genes for the current metabolite
    all_essential_genes[metabolite_name] = essential_gene_dict

for metabolite_name, essential_gene_dict in all_essential_genes.items():
    print("Essential genes for reactions producing {}: {}".format(metabolite_name, essential_gene_dict))


# In[46]:


print(len(precursors.items()))


# In[47]:


print(len(model_coli.metabolites))


# In[ ]:





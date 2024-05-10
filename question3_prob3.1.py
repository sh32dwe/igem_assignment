#!/usr/bin/env python
# coding: utf-8

# In[3]:


import cobra
import pulp


# In[2]:


get_ipython().system('pip install pulp')


# In[4]:


model = cobra.io.read_sbml_model("C:/Users/Dr Tauseef Ahmad/Downloads/iMM904.xml")


# In[5]:


genes = [gene.id for gene in model.genes]


# In[6]:


prob = pulp.LpProblem("OptKnock", pulp.LpMaximize)


# In[35]:


# Create an empty list to store reactions producing 2-phenylethanol
reactions_producing_phenylethanol = []

# Iterate over all metabolites in the model
for metabolite in model.metabolites:
    # Check if the metabolite name contains "2 phenylethanol"
    if "2 phenylethanol" in metabolite.name:
        print("Metabolite ID:", metabolite.id)
        print("Metabolite Name:", metabolite.name)
        print("Metabolite Formula:", metabolite.formula)
        print("Reactions producing this metabolite:")
        # Iterate over all reactions producing this metabolite
        for reaction in metabolite.reactions:
            print("- Reaction ID:", reaction.id)
            print("  Reaction Name:", reaction.name)
            print("  Reaction Formula:", reaction.reaction)
            reactions_producing_phenylethanol.append(reaction)
# Check if any reactions producing 2-phenylethanol were found
if reactions_producing_phenylethanol:
    print("Total reactions producing 2-phenylethanol:", len(reactions_producing_phenylethanol))
else:
    print("No reactions producing 2-phenylethanol were found.")


# In[37]:


ethanol_production = []
for reaction in reactions_producing_phenylethanol:
    for product in reaction.products:
        if '2phetoh_c'== reaction.products or '2phetoh_e'== reaction.products or '2phetoh_m'== reaction.products:
            ethanol_production.append(reaction)
            


# In[38]:


ethanol_production_variable = pulp.LpVariable("Ethanol_production", lowBound=0)
prob += ethanol_production_variable


# In[39]:


ethanol_production_increments = {}


# In[40]:


for gene in genes:
    # Define binary variables for gene knockout states
    gene_knockout = pulp.LpVariable("Gene_" + gene, cat=pulp.LpBinary)

    # Add lower-level objective: ethanol production <= ethanol production when gene is knocked out
    prob += ethanol_production_variable <= model.slim_optimize(
        {reaction: 0 for reaction in ethanol_production}
    )


# In[41]:


prob.solve()


# In[42]:


for gene in genes:
    # Check if the variable for the gene knockout exists in the optimization problem
    if "Gene_" + gene in prob.variablesDict():
        gene_knockout_var = prob.variablesDict()["Gene_" + gene]
        if gene_knockout_var.varValue == 1:
            ethanol_production_with_gene = model.slim_optimize({gene: 0})
            ethanol_production_without_gene = pulp.value(ethanol_production_variable)
            ethanol_production_increment = (
                ethanol_production_without_gene - ethanol_production_with_gene
            )
            ethanol_production_increments[gene] = ethanol_production_increment


# In[43]:


if ethanol_production_increments:
    max_increment_gene = max(
        ethanol_production_increments.keys(),
        key=lambda gene: ethanol_production_increments[gene],
    )
    print(
        "Gene causing the highest increment in ethanol production:",
        max_increment_gene,
        "with an increment of",
        ethanol_production_increments[max_increment_gene],
    )
else:
    print("No gene knockout resulted in increased ethanol production.")


# In[ ]:





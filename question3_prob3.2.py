#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
import pulp
import itertools


# In[2]:


model = cobra.io.read_sbml_model("C:/Users/Dr Tauseef Ahmad/Downloads/iMM904.xml") 


# In[3]:


genes = [gene.id for gene in model.genes]


# In[4]:


prob = pulp.LpProblem("OptKnock", pulp.LpMaximize)


# In[5]:


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


# In[6]:


ethanol_production = []
for reaction in reactions_producing_phenylethanol:
    for product in reaction.products:
        if '2phetoh_c'== reaction.products or '2phetoh_e'== reaction.products or '2phetoh_m'== reaction.products:
            ethanol_production.append(reaction)


# In[7]:


ethanol_production_variable = pulp.LpVariable("Ethanol_production", lowBound=0)
prob += ethanol_production_variable


# In[8]:


phenylethanol_production_increments_pairs = {}


# In[10]:


# Loop through all possible combinations of two genes
for gene_pair in itertools.combinations(genes, 2):
    # Define binary variables for gene knockout states
    gene_knockout_vars = {
        gene: pulp.LpVariable("Gene_" + gene, cat=pulp.LpBinary)
        for gene in gene_pair
    }

    # Add lower-level objective: 2-phenylethanol production <= 2-phenylethanol production when genes are knocked out
    prob += ethanol_production_variable <= model.slim_optimize(
        {gene: 0 for gene in gene_pair}
    )


# In[11]:


prob.solve()


# In[12]:


for gene_pair in itertools.combinations(genes, 2):
    gene_knockout_vars = {
        gene: pulp.LpVariable("Gene_" + gene, cat=pulp.LpBinary)
        for gene in gene_pair
    }
    if all(gene_knockout_vars[gene].varValue == 1 for gene in gene_pair):
        phenylethanol_production_with_genes = model.slim_optimize({gene: 0 for gene in gene_pair})
        phenylethanol_production_without_genes = model.slim_optimize()
        phenylethanol_production_increment = (
            phenylethanol_production_with_genes - phenylethanol_production_without_genes
        )
        phenylethanol_production_increments_pairs[gene_pair] = phenylethanol_production_increment


# In[16]:


if phenylethanol_production_increments_pairs:
    max_increment_gene_pair = max(
        phenylethanol_production_increments_pairs.keys(),
        key=lambda pair: phenylethanol_production_increments_pairs[pair],
    )
    print(
        "Pair of genes causing the highest increment in 2-phenylethanol production:",
        max_increment_gene_pair,
        "with an increment of",
        phenylethanol_production_increments_pairs[max_increment_gene_pair],
    )
else:
    print("No pair of gene knockouts resulted in increased 2-phenylethanol production.")


# In[ ]:





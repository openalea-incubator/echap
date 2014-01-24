import pandas
import numpy
import matplotlib.pyplot as plt

from alinea.adel.WheatTillering import WheatTillering

import alinea.echap.architectural_data as archidb

plt.ion()

def test_axis_dynamics():
    pgen = archidb.Mercia_2010_plantgen()
    tdata = archidb.Tillering_data_2010()
    
    # fit by Mariem
    primary_proba = pgen['decide_child_axis_probabilities']
    plant_density = pgen['plants_number']
    ears_density = pgen['ears_density']
    ears_per_plant = float(ears_density) / plant_density
    nff = 12
    m = WheatTillering(primary_tiller_probabilities=primary_proba, ears_per_plant = ears_per_plant, nff=nff)
    
    #newfit
    obs = tdata['tillering']
    obs = obs[obs.Var=='Mercia']
    primary_proba={'T1': obs['MB'].values + obs['TC'].values}#handle damages to MB and simplify TC
    for w in ('T2','T3','T4','T5'):
        primary_proba[w] = obs[w].values
    ears_per_plant = obs['MB'].values + obs['TT'].values
    nff = float(obs['Nff'].values)
    new_m = WheatTillering(primary_tiller_probabilities=primary_proba, ears_per_plant = ears_per_plant, nff=nff)
    
    #test plantgen fit
    
    fit = m.axis_dynamics(plant_density = plant_density)
    new_fit = new_m.axis_dynamics(plant_density = float(tdata['ear_density_at_harvest']['Mercia']) / new_m.ears_per_plant)
    fit.plot('HS', ['total','primary','others'], style=['--r','--g','--b'])
    new_fit.plot('HS', ['total','primary','others'],style=['-r','-g','-b'])
    plt.plot([1,13],[tdata['plant_density_at_emergence']['Mercia'],tdata['ear_density_at_harvest']['Mercia']] ,'ob')
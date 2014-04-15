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
    primary_proba={'T2': obs['MB'].values + obs['TC'].values}#handle damages to MB and simplify TC
    for w in ('T1','T3','T4','T5'):
        primary_proba[w] = obs[w].values
    ears_per_plant = obs['MB'].values + obs['TT'].values
    nff = float(obs['Nff'].values)
    new_m = WheatTillering(primary_tiller_probabilities=primary_proba, ears_per_plant = ears_per_plant, nff=nff)
    
    # todo : add comparison to primary emission
    
    #test plantgen fit
    
    fit = m.axis_dynamics(plant_density = plant_density)
    new_fit = new_m.axis_dynamics(plant_density = float(tdata['ear_density_at_harvest']['Mercia']) / new_m.ears_per_plant)
    fit.plot('HS', ['total','primary','others'], style=['--r','--g','--b'])
    new_fit.plot('HS', ['total','primary','others'],style=['-r','-g','-b'])
    plt.plot([1,13],[tdata['plant_density_at_emergence']['Mercia'],tdata['ear_density_at_harvest']['Mercia']] ,'ob')
    
def test_adel():
    from alinea.echap.macros_annual_loop import setup_canopy
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data,plantgen2adel
    from alinea.adel.stand.stand import agronomicplot
    from alinea.adel.astk_interface import AdelWheat
    from alinea.adel.AdelR import devCsv
    from alinea.adel.postprocessing import axis_statistics, plot_statistics 
    
    pgen = archidb.Mercia_2010_plantgen()
    axeT_, dimT_, phenT_, _, _, _, _, _, _, _, _ = gen_adel_input_data(**pgen)
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
    devT = devCsv(axeT, dimT, phenT)

    nplants, positions, domain, domain_area, convUnit = agronomicplot(length=0.2, 
                                                            width=0.15, 
                                                            sowing_density=pgen['plants_density'], 
                                                            plant_density=pgen['plants_number'],
                                                            inter_row=0.125)
    adel = AdelWheat(nplants=nplants, positions = positions, nsect=1, devT=devT, seed= 1, sample='random')
    dd = range(0,1500,100)
    outs = [adel.get_exposed_areas(g, convert=True) for g in (adel.setup_canopy(age) for age in dd)]
    out = reduce(lambda x,y: pandas.concat([x,y],ignore_index=True), outs)
    axstat = axis_statistics(out, domain_area, convUnit)    
    res =  plot_statistics(axstat, nplants, domain_area)
    return res
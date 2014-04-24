import pandas
import numpy
import os
import matplotlib.pyplot as plt

from alinea.adel.WheatTillering import WheatTillering
import alinea.echap.architectural_data as archidb

from multiprocessing import Pool,freeze_support
import itertools

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
    
    #test plantgen fit&
    
    fit = m.axis_dynamics(plant_density = plant_density)
    
    new_fit = new_m.axis_dynamics(plant_density = float(tdata['ear_density_at_harvest']['Mercia']) / new_m.ears_per_plant)
    fit.plot('HS', ['total','primary','others'], style=['--r','--g','--b'])
    new_fit.plot('HS', ['total','primary','others'],style=['-r','-g','-b'])
    plt.plot([1,13],[tdata['plant_density_at_emergence']['Mercia'],tdata['ear_density_at_harvest']['Mercia']] ,'ob')

def setAdel(TT_stop, var):
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data,plantgen2adel
    from alinea.adel.stand.stand import agronomicplot
    from alinea.adel.astk_interface import AdelWheat
    from alinea.adel.AdelR import devCsv
    pgen = archidb.Mercia_2010_plantgen()
    #pgen = archidb.Rht3_2010_plantgen()
    #pgen = archidb.Tremie_2011_plantgen()
    tdata = archidb.Tillering_data_Mercia_Rht3_2010_2011()
    #tdata = archidb.Tillering_data_Tremie1_2011_2012()
    obs = tdata['tillering']
    obs = obs[obs.Var==var]
    #primary_proba={'T2': obs['MB'].values + obs['TC'].values}#handle damages to MB and simplify TC
    primary_proba={}
    for w in ('T1','T2','T3','T4'):
        primary_proba[w] = obs[w].values
    ears_per_plant = obs['MB'].values + obs['TT'].values
    nff = float(obs['Nff'].values)
    new_m = WheatTillering(primary_tiller_probabilities=primary_proba, ears_per_plant = ears_per_plant, nff=nff)
    #update avec les probas primares de Tillering data
    pgen['decide_child_axis_probabilities'] = primary_proba
    # replace plant_number par eardensity / ear_per_plant, avec earsperplant =tillering data['MB'].values + Tillering data['TT'].values
    pgen['plants_number'] = float(tdata['ear_density_at_harvest'][var]) / new_m.ears_per_plant
    #
    pgen['delais_TT_stop_del_axis'] -= TT_stop
    tx = pgen['dynT_user'].a_cohort[0]
    axeT_, dimT_, phenT_, _, _, _, _, _, _, _, _ = gen_adel_input_data(**pgen)
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
    devT = devCsv(axeT, dimT, phenT)

    nplants, positions, domain, domain_area, convUnit = agronomicplot(length=0.5, 
                                                            width=0.5, 
                                                            sowing_density=pgen['plants_density'], 
                                                            plant_density=pgen['plants_number'],
                                                            inter_row=0.125)
    adel = AdelWheat(nplants=nplants, positions = positions, nsect=10, devT=devT, seed= 1, sample='random')
    return adel, nplants, positions, domain, domain_area, convUnit
      
def sample_adel():
    adel, nplants, positions, domain, domain_area, convUnit = setAdel(0, 'Mercia')
    return adel, domain

def simLAI(TT_stop, var):
    from alinea.adel.postprocessing import axis_statistics, plot_statistics 
    
    adel, nplants, positions, domain, domain_area, convUnit = setAdel(TT_stop,var)
    dd = range(0,2000,100)
    outs = [adel.get_exposed_areas(g, convert=True) for g in (adel.setup_canopy(age) for age in dd)]
    new_outs = [df for df in outs if not df.empty]
    out = reduce(lambda x,y: pandas.concat([x,y],ignore_index=True), new_outs)
    axstat = axis_statistics(out, domain_area, convUnit)    
    res =  plot_statistics(axstat, nplants, domain_area)
    return res, tx
    #res['HS'] = res.ThermalTime*tx
    #print res.plot('HS','PAI_vert',color='b')

def compare_LAI():
    sim = test_adel()
    #Mercia
    #obs = pandas.DataFrame({'TT':[1137,1281,1796],'PAI':[3.27,4.03,4.05]})
    #Rht3
    obs = pandas.DataFrame({'TT':[1137,1281,1796],'PAI':[3.40,3.96,3.45]})
    #Tremie2012
    #obs = pandas.DataFrame({'TT':[923,1260,1550],'PAI':[1.1,4.31,4.43]})
    #Tremie2013
    #obs = pandas.DataFrame({'TT':[917.5,1018.3,1377,1612.25],'PAI':[1.73,2.7,2.78,1.95]})
    #sim.plot('ThermalTime','PAI_vert')
    #obs.plot('TT','PAI',style='o')
    sim['HS'] = sim.ThermalTime*tx
    sim.plot('HS','PAI_vert',color='g')
    obs['HS'] = obs.TT*tx
    obs.plot('HS','PAI',style='o',color='r')

def draft_TC(adel, domain,thermal_time):
    from alinea.adel.postprocessing import ground_cover
    
    echap_top_camera =  {'type':'perspective', 'distance':200., 'fov':50., 'azimuth':0, 'zenith':0.}
    g =adel.setup_canopy(thermal_time)
    return ground_cover(g,domain, camera=echap_top_camera, image_width = 4288, image_height = 2848)
    
    
def draft_light(adel, domain, z_level, thermal_time):
    from alinea.caribu.caribu_star import diffuse_source, run_caribu
    from alinea.caribu.label import Label
    
    g =adel.setup_canopy(thermal_time)
    scene = adel.scene(g)
    sources = diffuse_source(46)
    out = run_caribu(sources, scene, domain=domain, zsoil = z_level)
    labs = {k:Label(v) for k,v in out['label'].iteritems()}
    ei_soil = [out['Ei'][k] for k in out['Ei'] if labs[k].is_soil()]
    return numpy.mean(ei_soil)
    
def mat_ray(csv_file):
    from pandas import read_csv
    df = read_csv(csv_file,';')
    col1 = df.DATE             #Ajout de ma part dans le csv
    col2 = df.QR_PAR_Avg       #PAR a 2m du sol
    col3 = df.SE_PAR_Avg_1     #PAR surface du sol
    col7 = df.SE_PAR_Avg_5     #PAR ds le couvert (env 20 cm)
    col1 = col1.astype(int)
    col2 = col2.astype(int)
    col3 = col3.astype(int)
    col7 = col7.astype(int)
    plt.plot(col1, col2)
    plt.plot(col1, col3)
    plt.plot(col1, col7)
    
def func_star(a_b):
    return test_adel(*a_b)
def main():
    p = Pool()
    TT_stop=[0,600,400,200]
    #var = "Mercia"
    var = "Rht3"
    #var = "Tremie"
    p.map(func_star, itertools.izip(TT_stop, itertools.repeat(var)))
# if __name__=="__main__":
    # freeze_support()
    # main()


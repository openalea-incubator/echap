import pandas
import numpy
import os
import matplotlib.pyplot as plt

from alinea.adel.WheatTillering import WheatTillering
import alinea.echap.architectural_data as archidb

from alinea.echap.architectural_reconstructions import reconst_db

from multiprocessing import Pool,freeze_support
import itertools

plt.ion()

def get_reconstruction(name='Mercia', **args):
    fun = reconst_db[name]
    _, adel, domain, domain_area, convUnit, nplants = fun(**args)
    return adel, domain, domain_area, convUnit, nplants

def get_pgen(name='Mercia', original = False, dTT_stop = 0):
    fun = reconst_db[name]
    pgen, _, _, _, _, _ = fun(nplants=1,nsect=1,as_pgen=original, dTT_stop=dTT_stop)
    
    return pgen
    
    
# test new tillering parametrisation againts original one by Mariem and data
#
# usage (ipython):
# %pylab
# %run TestCalibration.py
# test_axis_dynamics('Mercia')
#
# Very strange simulation for Tremie
#    
    
def test_axis_dynamics(name='Mercia'):
# todo : add comparison to primary emission

    pgen = get_pgen(name, original=True)
    newpgen = get_pgen(name, original=False)
    if name is not 'Tremie':
        obs = archidb.Plot_data_Mercia_Rht3_2010_2011()[name]
    else:
        obs = archidb.Plot_data_Tremie_2011_2012()[name]
    
    def _getfit(pgen):
        primary_proba = pgen['decide_child_axis_probabilities']
        plant_density = pgen['plants_density']
        ears_density = pgen['ears_density']
        ears_per_plant = float(ears_density) / plant_density
        v=list(pgen['MS_leaves_number_probabilities'].values())
        k=list(pgen['MS_leaves_number_probabilities'].keys())
        nff = int(k[v.index(max(v))]) # key of maximal value
        m = WheatTillering(primary_tiller_probabilities=primary_proba, ears_per_plant = ears_per_plant, nff=nff)
        fit = m.axis_dynamics(plant_density = plant_density)
        return fit
    
    # fit by Mariem
    fit = _getfit(pgen)
    #newfit
    new_fit = _getfit(newpgen)

    fit.plot('HS', ['total','primary','others'], style=['--r','--g','--b'])
    new_fit.plot('HS', ['total','primary','others'],style=['-r','-g','-b'])
    plt.plot([1,13],[obs['plant_density_at_emergence'],obs['ear_density_at_harvest']] ,'ob')

    
def simLAI(adel, domain_area, convUnit, nplants):
    from alinea.adel.postprocessing import axis_statistics, plot_statistics 
     
    dd = range(0,2000,100)
    outs = [adel.get_exposed_areas(g, convert=True) for g in (adel.setup_canopy(age) for age in dd)]
    new_outs = [df for df in outs if not df.empty]
    out = reduce(lambda x,y: pandas.concat([x,y],ignore_index=True), new_outs)
    axstat = axis_statistics(out, domain_area, convUnit)    
    res =  plot_statistics(axstat, nplants, domain_area)
    return res
    #res['HS'] = res.ThermalTime*tx
    #print res.plot('HS','PAI_vert',color='b')

def compare_LAI(name='Mercia', dTT_stop=0, original='False', n=30):

    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, dTT_stop=dTT_stop, as_pgen=original)
    pgen = get_pgen(name, dTT_stop=dTT_stop, original=original)
    sim = simLAI(adel, domain_area, convUnit, nplants)
    obs = archidb.PAI_data()[name]
    #sim.plot('ThermalTime','PAI_vert')
    #obs.plot('TT','PAI',style='o')
    tx = pgen['dynT_user'].a_cohort[0]
    sim['HS'] = sim.ThermalTime*tx
    obs['HS'] = obs.TT*tx
    
    sim.plot('HS','PAI_vert',color='g')
    obs.plot('HS','PAI',style='or')

def draft_TC(g, adel, domain, zenith):
    from alinea.adel.postprocessing import ground_cover
    
    echap_top_camera =  {'type':'perspective', 'distance':200., 'fov':50., 'azimuth':0, 'zenith':zenith}
    #g =adel.setup_canopy(thermal_time)
    return ground_cover(g,domain, camera=echap_top_camera, image_width = 4288, image_height = 2848)
    
def comp_TC(name='Mercia', original='False', n=30, zenith=0, zen='0', dTT_stop=0): #zenith = 0 or 57

    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, dTT_stop=dTT_stop, as_pgen=original, nplants=n)    
    dd = range(400,2000,300)
    sim = [adel.setup_canopy(age) for age in dd]

    TC_sim= [draft_TC(g, adel, domain, zenith) for g in sim]
    
    obs = archidb.TC_data()[name+'_'+zen]
    
    n=0; tc=[]; tc_sen=[]
    while n<=5:
        tc.append(TC_sim[n]['green'])
        tc_sen.append(TC_sim[n]['senescent'])
        n = n + 1

    sim_green = pandas.DataFrame({'tt':dd, 'TCgreen':tc})
    sim_sen = pandas.DataFrame({'tt':dd, 'TCsem':tc_sen})
    sim_green.plot('tt','TCgreen',color='g')
    sim_sen.plot('tt','TCsem',color='y')
    obs.plot('TT','TC',style='or')
     
def draft_light(g, adel, domain, z_level):
    from alinea.caribu.caribu_star import diffuse_source, run_caribu
    from alinea.caribu.label import Label

    scene = adel.scene(g)
    sources = diffuse_source(46)
    out = run_caribu(sources, scene, domain=domain, zsoil = z_level)
    labs = {k:Label(v) for k,v in out['label'].iteritems() if not numpy.isnan(float(v))}
    ei_soil = [out['Ei'][k] for k in labs if labs[k].is_soil()]
    return numpy.mean(ei_soil)
   
def comp_light(name='Mercia', dTT_stop=0, original='False', n=30):

    adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, dTT_stop=dTT_stop, as_pgen=original, nplants=n)    
    dd = range(400,2000,300)
    sim = [adel.setup_canopy(age) for age in dd]

    light_sim_0 = [draft_light(g, adel, domain, 0) for g in sim]
    light_sim_20 = [draft_light(g, adel, domain, 20) for g in sim]
    
    obs = archidb.mat_data()[name]
    
    sim_df0 = pandas.DataFrame({'tt':dd, 'light0':light_sim_0})
    sim_df20 = pandas.DataFrame({'tt':dd, 'light20':light_sim_20})
    sim_df0.plot('tt','light0',color='b')
    sim_df20.plot('tt','light20',color='g')
    
    obs.plot('TT','light0',style='ob')
    obs.plot('TT','light20',style='gs')
    
 
    
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


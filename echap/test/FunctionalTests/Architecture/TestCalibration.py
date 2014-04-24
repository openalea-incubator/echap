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
# %run TestCalibbration.py
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

def mat_ray(data_file):
    header_row = ['TIMESTAMP','RECORD','QR_PAR_Avg','SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4','SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8']
    data = pandas.read_csv(data_file, parse_dates={'datetime':[0,0]}, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=1, decimal='.')
    
    #filtre valeur negative
    df = data
    df_header = data[:0]
    col = df_header.columns
    i=1
    while i<10:
        df = df[(df[col[i]]>0)]
        i=i+1
    #verif : any(dd.SE_PAR_Avg_1 <=0)
    
    #avg PAR niveau du sol ('SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4')
    dat0 = df.ix[:,'SE_PAR_Avg_1':'SE_PAR_Avg_4']
    dat0 = dat0.apply(numpy.mean,1)
    #avg PAR env 20 cm du sol ('SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8')
    dat20 = df.ix[:,'SE_PAR_Avg_5':'SE_PAR_Avg_8']
    dat20 = dat20.apply(numpy.mean,1)
    # tableau % dat0/dat200
    tab_prc0 = dat0/df.QR_PAR_Avg
    # tableau % dat20/dat200 
    tab_prc20 = dat20/df.QR_PAR_Avg
    return tab_prc0, tab_prc20
    
def meteo_20112012(dates): #date au format JJ-MM-AAAA (2012-03-09) = ['2012-04-11','2012-05-09']
    data_file = 'METEO_stationINRA_20112012.csv'
    tab_prc0, tab_prc20 = mat_ray(data_file)
    prc0 = []; prc20 = []
    TT = [1260, 1551]
    
    for date in dates :
        seq = pandas.date_range(start = date, periods=13, freq='H')
        t_deb = seq[12]
        args0 = tab_prc0.index
        i=0;j=0
        
        for arg0 in args0:
            if arg0==t_deb:
                val = tab_prc0[i]
                val2 = tab_prc20[i]
            else:
                i=i+1;  j=j+1
        prc0.append(val)
        prc20.append(val2)
        
    print prc0, prc20
    plt.plot(TT, prc0, color='r')
    plt.plot(TT, prc20, color='g')
    
    
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


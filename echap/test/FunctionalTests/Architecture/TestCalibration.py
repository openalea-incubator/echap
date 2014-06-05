import pandas
import numpy
import os
import matplotlib.pyplot as plt

from alinea.adel.WheatTillering import WheatTillering
import alinea.echap.architectural_data as archidb

from alinea.echap.architectural_reconstructions import reconst_db

from multiprocessing import Pool,freeze_support
import itertools

import openalea.plantgl.all as pgl
from alinea.adel.mtg_interpreter import *
from alinea.echap.weather_data import *
from alinea.caribu.caribu_star import diffuse_source, run_caribu

plt.ion()

def get_reconstruction(name='Mercia', **args):
    fun = reconst_db[name]
    pgen, adel, domain, domain_area, convUnit, nplants = fun(**args)
    return pgen, adel, domain, domain_area, convUnit, nplants

def get_pgen(name='Mercia', original = False, dTT_stop = 0):
    fun = reconst_db[name]
    pgen, _, _, _, _, _ = fun(nplants=1,nsect=1,as_pgen=original, dTT_stop=dTT_stop)
    return pgen
 
# ANGLES
def curve(name='Mercia'):
    from openalea.deploy.shared_data import shared_data
    import re
    
    if name is 'Mercia':
        data_file_xydb = shared_data(alinea.echap, 'xydb_GrignonMercia2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv') 
        
    header_row_xydb = ['variety','variety_code','harvest','plant','rank','ranktop','relative_ranktop','HS','inerv','x','y']
    header_row_srdb = ['rankclass','s','r']
    dfxy = pandas.read_csv(data_file_xydb, names=header_row_xydb, sep=',', index_col=0, skiprows=1, decimal='.')
    dfsr = pandas.read_csv(data_file_srdb, names=header_row_srdb, sep=',', index_col=0, skiprows=1, decimal='.')
    dfsr = dfsr.reset_index()
    # NIV1 
    # cle = rankclass
    dfxy['rankclass'] = 1
    dfxy['rankclass'][dfxy['ranktop'] <= 4] = 2
    # filtre sur la variete
    dfxy = dfxy.reset_index()
    dfxy = dfxy[dfxy['variety']=='Mercia']
    # creation de la colonne age
    dfxy['age'] = dfxy['HS'] - dfxy['rank'] + 1
    #cut intervalle de 0 a 1, etc.
    dfxy_cut = pandas.cut(dfxy.age, [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    dfxy['age'] = dfxy_cut
    # bidouille...
    dfxy['age'] = dfxy['age'].replace("(-1, 0]", 0)
    dfxy['age'] = dfxy['age'].replace("(0, 1]", 1)
    dfxy['age'] = dfxy['age'].replace("(1, 2]", 2)
    dfxy['age'] = dfxy['age'].replace("(2, 3]", 3)
    dfxy['age'] = dfxy['age'].replace("(3, 4]", 4)
    dfxy['age'] = dfxy['age'].replace("(4, 5]", 5)
    dfxy['age'] = dfxy['age'].replace("(5, 6]", 6)
    dfxy['age'] = dfxy['age'].replace("(6, 7]", 7)
    # etc...
    # dfxy['age'] = dfxy['age'].replace("(7, 8]", 8)

    # NIV2
    # creation du premier dico par rankclass
    groups = dfxy.groupby(["rankclass"])
    dxy={n:{} for n,g in groups}
    # creation du second dico par age
    for n,g in groups:
        gg = g.groupby('age')
        dxy[n] = {k:[] for k,ggg in gg}
    # remplissage x et y
    groups = dfxy.groupby(["inerv"])   
    for n,d in groups:
        dxy[int(d[['rankclass']].values[0])][int(d[['age']].values[0])].append(d.ix[:,['x','y']].to_dict('list'))
    
    # NIV3 : ajout des s et r
    # traitement de dfsr (creation de 2 listes) DE STRING
    dr = 0
    s1 = []; r1 = []
    s2 = []; r2 = []
    while dr<len(dfsr):
        if dfsr['rankclass'][dr] == 1:
            s1.append(dfsr['s'][dr])
            r1.append(dfsr['r'][dr])
            dr = dr + 1
        else :
            s2.append(dfsr['s'][dr])
            r2.append(dfsr['r'][dr])
            dr = dr + 1
    # ajout dans le dict de dict precedent
    # longueur rankclass=1 et =2
    rank1 = 0 ; list1 = 0 ; rank2 = 1 ; list2 = 0
    while rank1 < len(dxy[1]):
        while list1 < len(dxy[1][rank1]) :
            dxy[1][rank1][list1].update({'s':s1, 'r':r1})
            list1 = list1 + 1
        if list1 == len(dxy[1][rank1]) :
            list1 = 0
        rank1 = rank1 + 1
        
    while rank2 <= len(dxy[2]) :
        while list2 < len(dxy[2][rank2]) :
            dxy[2][rank2][list2].update({'s':s2, 'r':r2})
            list2 = list2 + 1
        if list2 == len(dxy[2][rank2]) :
            list2 = 0
        rank2 = rank2 + 1
    
    return dxy
#-------------------------------------------------------------------------------------
    
# test new tillering parametrisation against original one by Mariem and data
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
    
    if name is 'Mercia':
        obs = archidb.Plot_data_Mercia_Rht3_2010_2011()[name]
    elif name is 'Rht3':
        obs = archidb.Plot_data_Mercia_Rht3_2010_2011()[name]
    elif name is 'Tremie':
        obs = archidb.Plot_data_Tremie_2011_2012()[name]
    elif name is 'Tremie2':
        obs = archidb.Plot_data_Tremie_2012_2013()[name]
    
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
    
    #fit by Mariem
    fit = _getfit(pgen)
    #newfit
    new_fit = _getfit(newpgen)

    fit.plot('HS', ['total','primary','others'], style=['--r','--g','--b'])
    new_fit.plot('HS', ['total','primary','others'],style=['-r','-g','-b'])
    plt.plot([1,13],[obs['plant_density_at_emergence'],obs['ear_density_at_harvest']] ,'ob')    

#------------------------------------------------------------------------------------- 
    
def simLAI(adel, domain_area, convUnit, nplants):
    from alinea.adel.postprocessing import axis_statistics, plot_statistics 
     
    dd = range(0,3000,100)
    outs = [adel.get_exposed_areas(g, convert=True) for g in (adel.setup_canopy(age) for age in dd)]
    new_outs = [df for df in outs if not df.empty]
    out = reduce(lambda x,y: pandas.concat([x,y],ignore_index=True), new_outs)
    axstat = axis_statistics(out, domain_area, convUnit)    
    res =  plot_statistics(axstat, nplants, domain_area)
    return res
    #res['HS'] = res.ThermalTime*tx
    #print res.plot('HS','PAI_vert',color='b')

def compare_LAI(name='Mercia', dTT_stop=0, original=False, n=30):

    # ajout des angles
    dxy = curve(name)

    pgen, adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, dTT_stop=dTT_stop, as_pgen=original, dict=dxy)
    sim = simLAI(adel, domain_area, convUnit, nplants)
    obs = archidb.PAI_data()[name]
    tx = pgen['dynT_user'].a_cohort[0]
    
    #Graph LAI_vert et LAI_tot en fonction des TT pour Corinne
    '''sim.plot('ThermalTime','LAI_vert', color='g')
    sim.plot('ThermalTime','LAI_tot', color='r')
    plt.title('Mercia 2010/2011')
    plt.legend(("LAI_vert", "LAI_tot"), 'best')
    plt.show()'''
    
    sim['HS'] = sim.ThermalTime*tx
    obs['HS'] = obs.TT*tx
    
    sim.plot('HS','PAI_vert',color='y')
    obs.plot('HS','PAI',style='or')
    
#------------------------------------------------------------------------------------- 

def draft_TC(g, adel, domain, zenith, rep):
    from alinea.adel.postprocessing import ground_cover
    
    echap_top_camera =  {'type':'perspective', 'distance':200., 'fov':50., 'azimuth':0, 'zenith':zenith}
    gc, im, box = ground_cover(g,domain, camera=echap_top_camera, image_width = 4288, image_height = 2848, getImages=True, replicate=rep)
    
    return gc
    
def comp_TC(name='Mercia', original=False, n=30, zenith=0, dTT_stop=0): #zenith = 0 or 57

    if zenith==0:
        zen='0'; rep=1
    else:
        zen='57';rep=2
        
    # ajout des angles
    dxy = curve(name)

    _, adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, dTT_stop=dTT_stop, as_pgen=original, nplants=n, dict=dxy)  
    return adel
    dd = range(400,2300,300)
    sim = [adel.setup_canopy(age) for age in dd]

    TC_sim = [draft_TC(g, adel, domain, zenith, rep) for g in sim]
    
    obs = archidb.TC_data()[name+'_'+zen]
    
    n=0; tc=[]; tc_sen=[]
    while n<len(TC_sim):
        tc.append(TC_sim[n]['green'])
        tc_sen.append(TC_sim[n]['senescent'])
        n = n + 1

    sim_green = pandas.DataFrame({'tt':dd, 'TCgreen':tc, 'TC_sen':tc_sen})
    
    # Si on veut tracer les courbes 1-TC et TC_tot pr comparer avec graph rayonnement obs/sim
    '''
    sim_green['1-TC']=1-sim_green['TCgreen']
    sim_green['TCtot']=sim_green['TCgreen']+sim_green['TC_sen']
    sim_green['1-TCtot']=1-sim_green['TCtot']
    sim_green.plot('tt','1-TC',color='y',linestyle='--',linewidth=2)
    sim_green.plot('tt','1-TCtot',color='y',linewidth=2)
    obs.plot('TT','1-TC',style='oy')
    '''
    
    # Si on veut tracer TC green et senescence sur le même graph
    sim_sen = pandas.DataFrame({'tt':dd, 'TCsem':tc_sen})
    sim_green.plot('tt','TCgreen',color='g')
    sim_sen.plot('tt','TCsem',color='y')
    obs.plot('TT','TC',style='or')

#-------------------------------------------------------------------------------------
# GRAPH METEO
    
def mat_ray_obs(data_file):
    import re
    
    header_row = ['DATE','H','RECORD','QR_PAR_Avg','SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4','SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8']
    #data = pandas.read_csv(data_file, parse_dates={'datetime':[0,1]}, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=1, decimal='.')
    df = pandas.read_csv(data_file, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=1, decimal='.')
     
    #filtre seulement heure entre 9h et 19h pour chaque journee
    df = df[df['H']>9]; df = df[df['H']<19]
    
    #tableau non traité
    dfa=df; dfa = dfa.reset_index()
    dfa.columns = ['datetime','H','RECORD','QR_PAR_Avg','SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4','SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8']
    dfa = dfa.groupby(dfa['datetime'], axis=0).mean(); dfa = dfa.reset_index()
    da=0
    while da<len(dfa):
        dfa['datetime'][da] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})", "\\1-\\2-\\3", dfa['datetime'][da])
        da = da+1
    dfa['%SE1'] = dfa['SE_PAR_Avg_1']/dfa['QR_PAR_Avg']; dfa['%SE2'] = dfa['SE_PAR_Avg_2']/dfa['QR_PAR_Avg']
    dfa['%SE3'] = dfa['SE_PAR_Avg_3']/dfa['QR_PAR_Avg']; dfa['%SE4'] = dfa['SE_PAR_Avg_4']/dfa['QR_PAR_Avg']
    dfa['%SE5'] = dfa['SE_PAR_Avg_5']/dfa['QR_PAR_Avg']; dfa['%SE6'] = dfa['SE_PAR_Avg_6']/dfa['QR_PAR_Avg']
    dfa['%SE7'] = dfa['SE_PAR_Avg_7']/dfa['QR_PAR_Avg']; dfa['%SE8'] = dfa['SE_PAR_Avg_8']/dfa['QR_PAR_Avg']
    
    # traitement tableau de données obs
    #avg PAR niveau du sol ('SE_PAR_Avg_1','SE_PAR_Avg_2','SE_PAR_Avg_3','SE_PAR_Avg_4')
    dat0 = df.ix[:,'SE_PAR_Avg_1':'SE_PAR_Avg_4']; dat0 = dat0.apply(numpy.mean,1)
    #avg PAR env 20 cm du sol ('SE_PAR_Avg_5','SE_PAR_Avg_6','SE_PAR_Avg_7','SE_PAR_Avg_8')
    dat20 = df.ix[:,'SE_PAR_Avg_5':'SE_PAR_Avg_8']; dat20 = dat20.apply(numpy.mean,1)
    # tableau % dat0/dat200
    tab_prc0 = dat0/df.QR_PAR_Avg
    # tableau % dat20/dat200
    tab_prc20 = dat20/df.QR_PAR_Avg
    
    # chaine de traitement pour hauteur au niveau du sol
    tab0 = pandas.DataFrame(tab_prc0)
    tab0 = tab0.reset_index(); tab0.columns = ['datetime','%']
    tab0 = tab0.groupby(tab0['datetime'], axis=0).mean()
    tab0 = tab0.reset_index()
    # changement format date JJ/MM/AAAA en JJ-MM-AAAA
    dt=0
    while dt<len(tab0):
        tab0['datetime'][dt] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})", "\\1-\\2-\\3", tab0['datetime'][dt])
        dt = dt+1

    # chaine de traitement pour hauteur a 20 cm du sol
    tab20 = pandas.DataFrame(tab_prc20)
    tab20 = tab20.reset_index(); tab20.columns = ['datetime','%']
    tab20 = tab20.groupby(tab20['datetime'], axis=0).mean()
    tab20 = tab20.reset_index()
    # changement format date JJ/MM/AAAA en JJ-MM-AAAA
    dte=0
    while dte<len(tab20):
        tab20['datetime'][dte] = re.sub("([0-9]{2})/([0-9]{2})/([0-9]{4})", "\\1-\\2-\\3", tab20['datetime'][dte])
        dte = dte+1
    
    return dfa, tab0, tab20
    
def draft_light(g, adel, domain, z_level):
    from alinea.caribu.label import Label

    scene = adel.scene(g)
    #modelisation afin de voir si erreur
    #scene=plot3d(g)
    #pgl.Viewer.display(scene)

    sources = diffuse_source(46)
    out = run_caribu(sources, scene, domain=domain, zsoil = z_level)
    labs = {k:Label(v) for k,v in out['label'].iteritems() if not numpy.isnan(float(v))}
    ei_soil = [out['Ei'][k] for k in labs if labs[k].is_soil()]
    return float(numpy.mean(ei_soil))
    
#fonction a lancer pour obtenir graph 'rayonnement obs contre rayonnement sim'
def graph_meteo(name='Mercia', dTT_stop=0, original=False, n=30):
    from alinea.astk.TimeControl import thermal_time # Attention la fonction marche seulement avec freq='H' !!! (pas en jour)
    
    # PARTIE OBS ------------------------------------------------------------------------------------------------------------
    if name is 'Mercia' or 'Rht3': #Mercia & Rht3 2010/2011
        data_file = 'METEO_stationINRA_20102011.csv'
        # correspondance entre date et TT
        met = Boigneville_2010_2011()
        seq = pandas.date_range(start="2010-10-15", end="2011-06-20", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2010-10-15", end="2011-06-20")
        bid = bid[seq]
    if name is 'Tremie': #Tremie 2011/2012
        data_file = 'METEO_stationINRA_20112012.csv'
        # correspondance entre date et TT
        met = Boigneville_2011_2012()
        seq = pandas.date_range(start="2011-10-21", end="2012-07-18", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2011-10-21", end="2012-07-18")
        bid = bid[seq]
        
    dfa, tab0, tab20 = mat_ray_obs(data_file)
        
    # merger bid avec tab0 puis tab20 
        # mise en forme de bid
    bid = pandas.DataFrame(list(bid.values), index=bid.index)
    bid = bid.reset_index()
    bid.columns = ['datetime','TT']
    bid['datetime'] = bid['datetime'].map(lambda x: x.strftime('%d-%m-%Y'))
        # merge bid avec tab0 puis tab20, clasement par ordre croissant
    tab0 = tab0.merge(bid); tab20 = tab20.merge(bid)
    tab0 = tab0.sort(['TT']); tab20 = tab20.sort(['TT'])
    
    # PARTIE SIM ------------------------------------------------------------------------------------------------------------
    _, adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, dTT_stop=dTT_stop, as_pgen=original, nplants=n)    
    dd = range(900,2200,100)
    sim = [adel.setup_canopy(age) for age in dd]

    light_sim_0 = [draft_light(g, adel, domain, z_level=0) for g in sim]
    #light_sim_5 = [draft_light(g, adel, domain, z_level=5) for g in sim]
    light_sim_20 = [draft_light(g, adel, domain, z_level=20) for g in sim]
    #light_sim_25 = [draft_light(g, adel, domain, z_level=25) for g in sim]

    sim0 = pandas.DataFrame({'TT':dd, 'light0':light_sim_0})
    #sim5 = pandas.DataFrame({'TT':dd, 'light5':light_sim_5})
    sim20 = pandas.DataFrame({'TT':dd, 'light20':light_sim_20})
    #sim25 = pandas.DataFrame({'TT':dd, 'light25':light_sim_25})

    # GRAPHES ---------------------------------------------------------------------------------------------------------------
    #obs tous les points
    dfa = dfa.merge(bid); dfa = dfa.sort(['TT']); print dfa.head()
        #niveau du sol (point rouge)
    dfa.plot('TT','%SE1',style='om',markersize=4); dfa.plot('TT','%SE2',style='om',markersize=4); dfa.plot('TT','%SE3',style='om',markersize=4); dfa.plot('TT','%SE4',style='om',markersize=4)
        #20cm du sol (point bleu)
    dfa.plot('TT','%SE5',style='oc',markersize=4); dfa.plot('TT','%SE6',style='oc',markersize=4); dfa.plot('TT','%SE7',style='oc',markersize=4); dfa.plot('TT','%SE8',style='oc',markersize=4)
    #obs moyennes
    tab0.plot('TT','%',color='m')
    tab20.plot('TT','%',color='c')
    #sim
    sim0.plot('TT','light0',color='r',linewidth=2)  
    #sim5.plot('TT','light5',color='r',linestyle='--',linewidth=2)     
    sim20.plot('TT','light20',color='b',linewidth=2)
    #sim25.plot('TT','light25',color='b',linestyle='--',linewidth=2)
 
#-------------------------------------------------------------------------------------

# test parallelisation
'''
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
'''


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
    tx = pgen['dynT_user'].a_cohort[0]
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

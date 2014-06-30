import pandas
import matplotlib.pyplot as plt
from alinea.adel.postprocessing import axis_statistics, plot_statistics 
from alinea.echap.architectural_reconstructions import reconst_db
from alinea.echap.recorder import LeafElementRecorder
from alinea.astk.TimeControl import *
from alinea.adel.WheatTillering import WheatTillering
from alinea.astk.Weather import Weather
from openalea.deploy.shared_data import shared_data
from alinea.echap.weather_data import *

from alinea.caribu.caribu_star import rain_and_light_star
from macros_annual_loop import update_lesions
from alinea.echap.interfaces import record as do_record

plt.ion()

def get_reconstruction(name='Mercia', **args):
    fun = reconst_db[name]
    pgen, adel, domain, domain_area, convUnit, nplants = fun(**args)
    return pgen, adel, domain, domain_area, convUnit, nplants
    
    
def sim_LAI_methode1(adel, domain_area, convUnit, nplants):
    dd = range(0,2300,100)
    
    outs = [adel.get_exposed_areas(g, convert=True) for g in (adel.setup_canopy(age) for age in dd)]
    new_outs = [df for df in outs if not df.empty]
    out = reduce(lambda x,y: pandas.concat([x,y],ignore_index=True), new_outs)
    axstat = axis_statistics(out, domain_area, convUnit)    
    res =  plot_statistics(axstat, nplants, domain_area)
    return res

def sim_LAI_methode2(adel, domain_area, convUnit, nplants, adelpars):
    delta_dd = 100; dd_start = 0; dd_end = 2300
    dd = range(dd_start, dd_end, delta_dd)
    
    g = adel.setup_canopy(dd_start)
    gouts = [adel.get_exposed_areas(g, convert=True)]
    #gouts = adel.get_exposed_areas(g, convert=True)
    for i in range((dd_end - dd_start) / delta_dd):
        g = adel.grow_dd(g, delta_dd, adelpars)
        gouts.append(adel.get_exposed_areas(g, convert=True))
        #gouts = gouts.append(adel.get_exposed_areas(g, convert=True)) 
    new_gouts = [df for df in gouts if not df.empty]
    gout = reduce(lambda x,y: pandas.concat([x,y],ignore_index=True), new_gouts)
    axstat = axis_statistics(gout, domain_area, convUnit)    
    res =  plot_statistics(axstat, nplants, domain_area)
    return res
    
def recorder_f(weather, adel, domain, convUnit, seq, adelpars):
    Tbase = 0; delay = 100
    recorder = LeafElementRecorder()
    TTmodel = DegreeDayModel(Tbase)
    
    every_dd = thermal_time_filter(seq, weather, TTmodel, delay)
    canopy_timing = IterWithDelays(*time_control(seq, every_dd, weather.data))
    g = adel.setup_canopy(Tbase, adelpars)
    
    for i, controls in enumerate(canopy_timing):
        canopy_iter = controls
        if canopy_iter:
            print canopy_iter.value.index[-1]
            g = adel.grow_dd(g, delay, adelpars=adelpars)
            _=rain_and_light_star(g, light_sectors = '1', domain=domain, convUnit=convUnit)
            _=update_lesions(g, canopy_iter.dt, True)  
            _=do_record(g, canopy_iter.value, recorder)
    recorder.save_records('test_plot.csv')
    rec = recorder.get_records()
    return rec
    
def plot(dataframefile, what, t, by, axe, plant, style):
    from itertools import cycle, islice
    d = dataframefile
    d = d[d[what].notnull()]
    if plant == 'all':
        dm=d[d.axe==axe]
    else:
        dm=d[(d.axe==axe) & (d.plant==plant)]
    gr=dm.groupby(t)
    dmp=gr.agg(numpy.sum)
    dmp = dmp.reset_index()
    '''gr=dmp.groupby([by,t])
    dms=gr.agg(numpy.mean)
    dms=dms.reset_index()
    gr=dms.groupby(by)'''
    dmp.plot(t, what, style=style)
    
def plot_LAI(name='Mercia', dTT_stop=0, original=False, n=30):

    pgen, adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, dTT_stop=dTT_stop, as_pgen=original)
    tx = pgen['dynT_user'].a_cohort[0]
    
    adelpars = {'senescence_leaf_shrink' : 1, 'startLeaf' : -0.4, 'endLeaf' : 1.6, 'endLeaf1': 1.6, 'stemLeaf' : 1.2, 'epsillon' : 1e-6}
    
    # PLOT SIM LAI TESTCALIBRATION.PY ET AUTRE METHODE GROW_DD
    outs = sim_LAI_methode1(adel, domain_area, convUnit, nplants)
    gouts = sim_LAI_methode2(adel, domain_area, convUnit, nplants, adelpars)
    #print outs
    #print gouts
    
    # plot sim 1
    outs['HS'] = outs.ThermalTime*tx
    outs.plot('HS', 'LAI_vert', style='or')
    outs.plot('HS', 'LAI_tot', style='om')
    # plot sim 2
    gouts['HS'] = gouts.ThermalTime*tx
    gouts.plot('HS', 'LAI_vert', color='r')
    outs.plot('HS', 'LAI_tot', color='m')
    
    plt.show()
    
    # RECORDER
    weather = Boigneville_2010_2011()
    seq = pandas.date_range(start = "2010-10-15", periods = 6000, freq = 'H')  # 5000 pour le cycle complet 
    rec = recorder_f(weather, adel, domain, convUnit, seq, adelpars)
    # concordance date AAAA-MM-JJ en TT
    bid = thermal_time(seq, weather.data)
    bid = bid[seq]
    bid = pandas.DataFrame(list(bid.values), index=bid.index)
    bid = bid.reset_index()
    bid.columns = ['date','TT']
    bid['HS'] = bid['TT']*tx
    # ajout colonne TT dans rec
    rec['date'] = rec['date'].map(lambda x: x.strftime('%Y-%m-%d'))
    bid['date'] = bid['date'].map(lambda x: x.strftime('%Y-%m-%d'))
    rec_all = pandas.merge(bid,rec,on='date') 
    rec_all = rec_all.drop('TT', 1)
    # facteur echelle ordonnees 
    rec_all['green_area_unit'] = rec_all['green_area']/1000
    rec_all['tot_area'] = rec_all['green_area'] + rec_all['senesced_area']
    rec_all['tot_area_unit'] = rec_all['tot_area']/1000
    # plot recorder, colonne green_area
    plot(rec_all, what='green_area_unit', t = 'HS', by='ntop', axe = 'MS', plant='all', style='ob')
    plot(rec_all, what='tot_area_unit', t = 'HS', by='ntop', axe = 'MS', plant='all', style='og')

'''  
def test_setup(name='Mercia', dTT_stop=0, original=False, n=30):
    pgen, adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, dTT_stop=dTT_stop, as_pgen=original)
    ages = range(0, 2300, 100)
    recorder = LeafElementRecorder()
    for a in ages:
        g=adel.setup_canopy(a)
        header={'dd':a}
        for vid in g:
            if g.label(vid).startswith('LeafElement'):
                n = g.node(vid)
                header.update({'plant' : n.complex().complex().complex().complex().label,
                                'axe' : n.complex().complex().complex().label,
                                'metamer' : int(''.join(list(n.complex().complex().label)[7:])),
                                'organ' : n.complex().label,
                                'ntop' : n.complex().ntop,
                                'id' : n._vid})
 
                recorder.record(n,header)
        #adel.plot(g)
    recorder.save_records('test.csv')
    #return g, recorder

from macros_annual_loop import setup_canopy
from alinea.echap.recorder import LeafElementRecorder

g, adel, domain, domain_area, convUnit, nplants = setup_canopy(length=0.2, nsect=3)
    
    
def test_grow():
    g=adel.setup_canopy(100)
    recorder = LeafElementRecorder()
    for i in range(24):
        g=adel.grow_dd(g,100)
        header={'dd':adel.canopy_age}
        for vid in g:
            if g.label(vid).startswith('LeafElement'):
                n = g.node(vid)
                header.update({'plant' : n.complex().complex().complex().complex().label,
                                'axe' : n.complex().complex().complex().label,
                                'metamer' : int(''.join(list(n.complex().complex().label)[7:])),
                                'organ' : n.complex().label,
                                'ntop' : n.complex().ntop,
                                'id' : n._vid
                         })
 
                recorder.record(n,header)

        adel.plot(g)
    return g, recorder
        
'''  
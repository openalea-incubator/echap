# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:59:54 2013

@author: lepse
"""
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import mpl
from matplotlib import cm
from matplotlib.colors import Normalize, LogNorm
from pandas import *

from alinea.pearl.pearl import *
from alinea.pearl.pearl_leaf import *

from alinea.echap.interfaces import pesticide_surfacic_decay
from alinea.echap.milne_leaf import *
from alinea.echap.interfaces import pesticide_penetrated_decay
from alinea.echap.interception_leaf import *
from alinea.echap.interfaces import pesticide_interception
from alinea.echap.microclimate_leaf import *
from alinea.echap.interfaces import local_microclimate
from alinea.echap.interfaces import pesticide_efficacy
from alinea.echap.wheat_mtg import *

from alinea.echap.color_map import *

from alinea.pesticide_efficacy.pesticide_efficacy import *
from alinea.weather.global_weather import *

from datetime import datetime, timedelta
from openalea.deploy.shared_data import get_shared_data_path


############# update

def update_on_leaves(g, label = 'LeafElement'):
    """ Read weather data for a step of simulation and apply it to each leaf.
    """        
    vids = [vid for vid in g if g.label(vid).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.surfacic_doses={'Chlorothalonil':1,'Epoxiconazole':2}
        n.penetrated_doses={'Chlorothalonil':1,'Epoxiconazole':2}
        n.temp = 12
        n.rain_intensity = 0
        n.relative_humidity = 100 
        n.wetness = True
        n.microclimate = {}
    return g


def update_no_doses(g, label = 'LeafElement'):
    """ Read weather data for a step of simulation and apply it to each leaf.
    """        
    vids = [vid for vid in g if g.label(vid).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.temp = 12
        n.rain_intensity = 0
        n.relative_humidity = 100 
        n.wetness = True
        n.microclimate = {}
    return g


########################## initiate
#data
db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}, 
'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}}

compound_parameters = [{'Ae': 0.5, 'Ke': 7.0099999999999998, 'decay_rate': 0.069000000000000006,
 'dose_max_ha': 125, 'Ap': 0.70999999999999996, 'Kp': 6, 'type_code': 2, 
 'compound': 'Epoxiconazole'}, {'Ae': 1.0, 'Ke': 6.4900000000000002, 'decay_rate': 0.01, 
 'dose_max_ha': 1000, 'Ap': 0.0, 'Kp': 0, 'type_code': 3, 'compound': 'Chlorothalonil'}]

 # Loop
t = 0
dt = 1
nb_steps = 2

# Param
timestep=1
product_name='Opus'
dose=1.5
meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
label='LeafElement'
t_deb = "2000-10-01 01:00:00"

# Initialisation du mtg 
g = adel_mtg()
g = update_no_doses(g)
scene_geometry = g.property('geometry')

# models
interception_model = CaribuInterceptModel()
Pearl_decay_model = PearLeafDecayModel(db)
Milne_decay_model = PenetratedDecayModel()
climate_model = MicroclimateLeaf()
weather = Weather(data_file=meteo01_filepath)

    
########################## extraction d'un data frame des doses de pesticide sur le mtg

def get_df_out(time,g):
    sd = g.property('surfacic_doses')
    lab = g.property('label')
    sp = g.property('penetrated_doses')            
    recs = [(id,lab[id],comp,dose) for id,v in sd.iteritems() for comp, dose in v.iteritems() if lab[id] is 'LeafElement']
    ids,labs,compounds,surfdose = zip(*recs)
    dfs = DataFrame({'time':[time]*len(ids), 'id' : ids, 'label': labs,'compound' : compounds, 'surfacic_doses' : surfdose})
    if not 'penetrated_doses' in g.property_names():
        df=dfs        
        #dfp = DataFrame(columns=('time', 'id', 'label','compound', 'penetrated_dose'))
    else:
        recp = [(id,lab[id],comp,dose) for id,v in sp.iteritems() for comp, dose in v.iteritems() if lab[id] is 'LeafElement']
        idp,labp,compoundp,pendose = zip(*recp)
        dfp = DataFrame({'time':[time]*len(idp), 'id' : idp, 'label': labp,'compound' : compoundp, 'penetrated_doses' : pendose})
        df = merge(dfs, dfp, left_on=('compound', 'id', 'label', 'time'), right_on=('compound', 'id', 'label', 'time'), how='outer')    
    return df


########################## display

def plot_decay(out, leaf=12):
    #from pylab import *
    df = out[out['id']==leaf]
    plt.plot(df['time'], df['surfacic_doses'])
    plt.plot(df['time'], df['penetrated_doses'])
    plt.show()


def plot_pesticide(g, property_name='surfacic_doses', compound_name='Epoxiconazole', cmap=green_lightblue_blue):
    """ plot the plant with pesticide doses """
    if type(cmap) is str:
        try:
            _cmap = cm.get_cmap(cmap())
        except:
            raise Exception('This colormap does not exist')
    else:
        _cmap = cmap()

    green = (0,180,0)

    for v in g.vertices(scale=g.max_scale()): 
        n = g.node(v)
        if 'surfacic_doses' in n.properties():
            r,gg,b,s = _cmap(n.surfacic_doses[compound_name]*200)
            n.color = (int(r*255),int(gg*255),int(b*255))           
        else : 
            n.color = green
    scene = plot3d(g)
    Viewer.display(scene)
    return g


def dose_norm(dose, dose_max_ha):
    """ normalise doses(g.m-2) """
    dn = 0
    if dose_max_ha > 0 :
        dn = float(dose * 1e4) / dose_max_ha
    return dn

def plot_pesticide_norm(g, property_name='surfacic_doses', compound_name='Epoxiconazole', cmap=green_lightblue_blue, lognorm=False):
    """ plot the plant with pesticide doses """
    prop = g.property(property_name)
    keys = prop.keys()
    value = []
    for k, val in prop.iteritems():
        value.append(val[compound_name])
        val = np.array(value)

    if type(cmap) is str:
        try:
            _cmap = cm.get_cmap(cmap())
        except:
            raise Exception('This colormap does not exist')
    else:
        _cmap = cmap()

    green = (0,180,0)

    # norm = Normalize(vmin=0, vmax=max(v)) if not lognorm else LogNorm(vmin=0, vmax=max(v)) 
    # values = norm(v)
    for i in range(0,len(value)):
        val[i] = dose_norm(value[i], dose_max_ha)
    
    colors = (_cmap(val)[:,0:3])*255
    colors = np.array(colors,dtype=np.int).tolist()

    for vid in g.vertices(scale=g.max_scale()): 
        n = g.node(vid)
        if 'surfacic_doses' in n.properties():
            n.color = tuple(dict(zip(keys,colors))[vid])
        else : 
            n.color = green

    scene = plot3d(g)
    Viewer.display(scene)
    return g


########################## import .csv as dict

def products_from_csv(csvname, delimiter = ';'):
    """ 
    Read a csv of products parameters and import them in a dict.
    Expected columns are :
        - 'product' : commercial name of the product
        - 'compound' : name of the active compound of the product
        - 'dose' : dose of active compound in the product (g.l-1)    
    """
    tab = recfromcsv(csvname, delimiter = delimiter, case_sensitive = True)
    d = {}
    for i in range(0,len(tab['compound'])):
        d[tab[i][0]] = {tab[i][1]:tab[i][2]}
    return d


########################## tests decays

def test_surfacic():
    g = adel_mtg2()
    g = update_no_doses(g)
    db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}}
    decay_model = PearLeafDecayModel(db)
    g = pesticide_surfacic_decay(g, decay_model)
    return g

def test_penetrated():
    g = adel_mtg2()
    g = update_on_leaves(g)
    decay_model = PenetratedDecayModel()
    g = pesticide_penetrated_decay(g, decay_model)
    return g

def test_all_decay():
    g = adel_mtg()
    g = update_on_leaves(g)
    db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}}
    decay_model = PearLeafDecayModel(db)
    g = pesticide_surfacic_decay(g, decay_model)
    decay_model = PenetratedDecayModel()
    g = pesticide_penetrated_decay(g, decay_model)
    return g

def test_no_doses():
    g = adel_mtg()
    g = update_no_doses(g)
    db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}}
    decay_model = PearLeafDecayModel(db)
    g = pesticide_surfacic_decay(g, decay_model)
    decay_model = PenetratedDecayModel()
    g = pesticide_penetrated_decay(g, decay_model)
    return g

###################################### tests interception

def test_intercept():
    g = adel_mtg()
    g = update_on_leaves(g)
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    return g


def test_intercept_no_dose():
    g = adel_mtg()
    g = update_no_doses(g)
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    return g


def test_intercept_elevation():
    g = adel_mtg()
    g = update_no_doses(g)
    interception_model = CaribuInterceptModel(elevation=90.)
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    print g.property('surfacic_doses')
    interception_model = CaribuInterceptModel(elevation=45.)
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    print g.property('surfacic_doses')
    interception_model = CaribuInterceptModel(elevation=1.)
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    print g.property('surfacic_doses')
    return g

###################################### test microclimate

def test_microclimate():
    g = adel_mtg()
    g = update_no_doses(g)
    climate_model = MicroclimateLeaf()
    meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
    weather = Weather(data_file=meteo01_filepath)
    t_deb = "2000-10-01 01:00:00"
    g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
    print g.property('microclimate')
    return g


###################################### test efficacy

def test_efficacy():
    g = adel_mtg()
    g = update_no_doses(g)
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    efficacy_model = PesticideEfficacyModel()
    g = pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
    return g

def test_efficacy_nopest():
    g = adel_mtg()
    efficacy_model = PesticideEfficacyModel()
    g = pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
    return g

    
###################################### test lesions

def test_update():
    g = adel_mtg()
    g = update_no_doses(g)
    # Interception
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    # Microclimate
    climate_model = MicroclimateLeaf()
    meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
    weather = Weather(data_file=meteo01_filepath)
    g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
    # Efficacy
    efficacy_model = PesticideEfficacyModel()
    g = pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
    # Update lesion
    p = ParCycle()
    LesionsCycle = Lesions(p)
    g = update_lesions(g, LesionsCycle, label="LeafElement", timestep=1)
    return g


##################################### loop test

def test_decay_doses():
    db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}}
    t_deb = "2000-10-01 01:00:00"
    # Loop
    t = 0
    dt = 1
    nb_steps = 6
    # Initialisation du mtg 
    g = adel_mtg()
    g = update_no_doses(g)
    # models
    interception_model = CaribuInterceptModel()
    Pearl_decay_model = PearLeafDecayModel(db)
    Milne_decay_model = PenetratedDecayModel()
    climate_model = MicroclimateLeaf()
    meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
    weather = Weather(data_file=meteo01_filepath)
    # Interception
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
    t_deb = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[3]
    # sauvegarde etat initial
    out = get_df_out(0,g)
    # loop
    for i in range(nb_steps):
        t += dt        
        # Surfacic decay
        g = pesticide_surfacic_decay(g, Pearl_decay_model, timestep=dt)
        # Penetrated decay
        g = pesticide_penetrated_decay(g, Milne_decay_model, timestep=dt)
        df = get_df_out(t,g)
        out = out.append(df)
        plot_pesticide_norm(g, property_name='surfacic_doses', compound_name='Epoxiconazole', cmap=green_lightblue_blue, lognorm=False)
        #plot_pesticide(g, property_name='surfacic_doses', compound_name='Epoxiconazole', cmap=green_lightblue_blue)
        g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
        t_deb = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[3]
    return out


def test_local_meteo():
    # Loop
    t_deb = "2000-10-01 01:00:00"
    t = 0
    dt = 1
    nb_steps = 2
    # Initialisation du mtg 
    g = adel_mtg()
    g = update_no_doses(g)
    # models
    interception_model = CaribuInterceptModel()
    climate_model = MicroclimateLeaf()
    meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
    weather = Weather(data_file=meteo01_filepath)
    # Interception
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    print t_deb
    # Microclimate
    g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
    t_deb = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[3]
    print g.property('microclimate')
    # loop
    for i in range(nb_steps):
        t += dt    
        print t_deb
        # Microclimate
        g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
        t_deb = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[3]
        print g.property('microclimate')





   


    
    
    

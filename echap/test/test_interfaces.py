# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:59:54 2013

@author: lepse
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mpl
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
from alinea.pesticide_efficacy.pesticide_efficacy import *
from alinea.echap.interfaces import pesticide_efficacy
from alinea.echap.global_meteo import *

from alinea.adel.newmtg import *
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
import alinea.adel.fitting as fitting
from alinea.adel.AdelR import devCsv,setAdel,RunAdel,genGeoLeaf,genGeoAxe


############# Création d'un MTG

def adelR(nplants,dd):
    devT = devCsv('../../adel/example/data/axeTCa0N.csv','../../adel/example/data/dimTCa0N.csv','../../adel/example/data/phenTCa0N.csv','../../adel/example/data/earTCa0N.csv','../../adel/example/data/ssi2sen.csv')
    geoLeaf = genGeoLeaf()
    geoAxe = genGeoAxe()
    pars = setAdel(devT,geoLeaf,geoAxe,nplants)
    cantable = RunAdel(dd,pars)
    return pars,cantable

def leaves_db():
    import cPickle as Pickle
    fn = 'E:/openaleapkg/adel/adel/data/leaves_simple.db'
    f = open(fn)
    leaves = Pickle.load(f)
    f.close()
    leaves,discard = fitting.fit_leaves(leaves, 9)
    return leaves

def leaves_db_flow(fn):
    import cPickle as Pickle
    f = open(fn)
    leaves = Pickle.load(f)
    f.close()
    leaves,discard = fitting.fit_leaves(leaves, 9)
    return leaves

def adel_mtg():
    """ create a very simple adel mtg """
    d = {'plant':[1,1],'axe_id':['MS','T1'],'ms_insertion':[0,1],'numphy':[1,1], 
         'Laz': [0,90], 'Ll' :[3,3], 'Lv' :[3,3] , 'Lsen':[0,0], 'L_shape':[3,3], 'Lw_shape':[.3,.3], 'Linc':[0,0],
         'Einc':[0,45],'El':[1,1],'Ev':[1,1],'Esen':[0,0],'Ed': [0.1,0.1]}
    g=mtg_factory(d,adel_metamer,leaf_db=leaves_db(), leaf_sectors=1)
    g=mtg_interpreter(g)
    return g

def adel_mtg2(nb_sect=1):
    """ create a less simple adel mtg """
    p, d = adelR(3,1000)
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,leaf_db=leaves_db(),stand=[((0,0,0),0),((10,0,0),90), ((0,10,0), 0)])
    g=mtg_interpreter(g)
    return g

def adel_mtg3(nb_sect=1, leaf_db=None, d=None, p=None):
    """ create a less simple adel mtg """
    if p: # nb_plants
        size = int(ceil(sqrt(p)))
        stand = numpy.array([(i, j) for i in range(size) for j in range(size)])
        numpy.random.shuffle(stand)
        stand = [((i, j, 0),random.randint(0,90)) for i, j in stand[:p]]
    else:
        stand = [((0,0,0),0),((10,0,0),90), ((0,10,0), 0)]
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,leaf_db=leaf_db,stand=stand)
    g=mtg_interpreter(g)
    return g

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


def plot_pesticide(g, compound_name='Epoxiconazole', colmap='winter_r'):
    """ plot the plant with pesticide doses """
    cmap = mpl.cm.get_cmap(colmap)
    green = (0,180,0)
    for v in g.vertices(scale=g.max_scale()): 
        n = g.node(v)
        if 'surfacic_doses' in n.properties():
            r,gg,b,s= cmap(n.surfacic_doses[compound_name]*100)
            n.color = (int(r*255),int(gg*255),int(b*255))           
        else : 
            n.color = green
    scene = plot3d(g)
    Viewer.display(scene)
    return g


########################## tests decays

def test_surfacic():
    g = adel_mtg()
    g = update_no_doses(g)
    #db = {'Chlorothalonil':{}, 'Epoxiconazole':{}, 'Metconazole':{}}
    db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'FacTraDepRex':0, 
    'FacVolDepRex':0,'FacPenDepRex':0,'FacWasDepRex':0,'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'FacTraDepRex':0, 
    'FacVolDepRex':0,'FacPenDepRex':0,'FacWasDepRex':0,'FraDepRex':0}}
    decay_model = PearLeafDecayModel(db)
    g = pesticide_surfacic_decay(g, decay_model)
    return g

def test_penetrated():
    g = adel_mtg()
    g = update_on_leaves(g)
    decay_model = PenetratedDecayModel()
    g = pesticide_penetrated_decay(g, decay_model)
    return g

def test_all_decay():
    g = adel_mtg()
    g = update_on_leaves(g)
    db = {'Chlorothalonil':{}, 'Epoxiconazole':{}, 'Metconazole':{}}
    decay_model = PearLeafDecayModel(db)
    g = pesticide_surfacic_decay(g, decay_model)
    decay_model = PenetratedDecayModel()
    g = pesticide_penetrated_decay(g, decay_model)
    return g

def test_no_doses():
    g = adel_mtg()
    g = update_no_doses(g)
    db = {'Chlorothalonil':{}, 'Epoxiconazole':{}, 'Metconazole':{}}
    decay_model = PearLeafDecayModel(db)
    g = pesticide_surfacic_decay(g, decay_model)
    decay_model = PenetratedDecayModel()
    g = pesticide_penetrated_decay(g, decay_model)
    return g

###################################### tests interception

def test_intercept():
    g = adel_mtg()
    g = update_on_leaves(g)
    g = mtg_interpreter(g) 
    scene = plot3d(g) 
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, scene, interception_model, product_name='Opus new', dose=1.5)
    return g

def test_intercept_no_dose():
    g = adel_mtg()
    g = update_no_doses(g)
    g = mtg_interpreter(g)
    scene = plot3d(g) 
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, scene, interception_model, product_name='Opus new', dose=1.5)
    return g


###################################### test microclimate

def test_microclimate():
    g = adel_mtg()
    g = update_no_doses(g)
    g = mtg_interpreter(g)
    scene = plot3d(g) 
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, scene, interception_model, product_name='Opus new', dose=1.5)
    climate_model = CaribuMicroclimModel()
    g = local_microclimate(g, scene, climate_model, t_deb='2000-10-01 09:00:00', label='LeafElement', timestep=1)[0]
    print g.property('microclimate')
    return g


###################################### test efficacy

def test_efficacy():
    g = adel_mtg()
    g = update_no_doses(g)
    g = mtg_interpreter(g)
    scene = plot3d(g) 
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, scene, interception_model, product_name='Opus new', dose=1.5)
    efficacy_model = PesticideEfficacyModel()
    g = pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
    return g

def test_efficacy_nopest():
    g = adel_mtg()
    efficacy_model = PesticideEfficacyModel()
    g = pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
    return g


##################################### loop test

def test_decay_doses():
    db = {'Chlorothalonil':{}, 'Epoxiconazole':{}, 'Metconazole':{}}
    t_deb = '2000-10-01 01:00:00'
    # Loop
    t = 0
    dt = 1
    nb_steps = 4
    # Initialisation du mtg 
    g = adel_mtg()
    g = update_no_doses(g)
    g = mtg_interpreter(g)
    scene = plot3d(g) 
    # models
    interception_model = CaribuInterceptModel()
    Pearl_decay_model = PearLeafDecayModel(db)
    Milne_decay_model = PenetratedDecayModel()
    meteo_reader = Meteo()
    climate_model = CaribuMicroclimModel()
    # Interception
    g = pesticide_interception(g, scene, interception_model, product_name='Opus new', dose=1.5)
    g = local_microclimate(g, scene, climate_model, meteo_reader, t_deb=t_deb, label='LeafElement', timestep=dt)[0]
    t_deb = local_microclimate(g, scene, climate_model, meteo_reader, t_deb=t_deb, label='LeafElement', timestep=dt)[3]
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
        plot_pesticide(g, compound_name='Epoxiconazole', colmap='winter_r')
        g = local_microclimate(g, scene, climate_model, meteo_reader, t_deb=t_deb, label='LeafElement', timestep=dt)[0]
        t_deb = local_microclimate(g, scene, climate_model, meteo_reader, t_deb=t_deb, label='LeafElement', timestep=dt)[3]
    return out


def test_local_meteo():
    # Loop
    t_deb = '2000-10-01 02:00:00'
    t = 0
    dt = 1
    nb_steps = 2
    # Initialisation du mtg 
    g = adel_mtg()
    g = update_no_doses(g)
    g = mtg_interpreter(g)
    scene = plot3d(g) 
    # models
    interception_model = CaribuInterceptModel()
    climate_model = CaribuMicroclimModel()
    # Interception
    g = pesticide_interception(g, scene, interception_model, product_name='Opus new', dose=1.5)
    print t_deb
    # Microclimate
    g = local_microclimate(g, scene, climate_model, t_deb=t_deb, label='LeafElement', timestep=dt)[0]
    t_deb = local_microclimate(g, scene, climate_model, t_deb=t_deb, label='LeafElement', timestep=dt)[3]
    print g.property('microclimate')
    # loop
    for i in range(nb_steps):
        t += dt    
        print t_deb
        # Microclimate
        g = local_microclimate(g, scene, climate_model, t_deb=t_deb, label='LeafElement', timestep=dt)[0]
        t_deb = local_microclimate(g, scene, climate_model, t_deb=t_deb, label='LeafElement', timestep=dt)[3]
        print g.property('microclimate')





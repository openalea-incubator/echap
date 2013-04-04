# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:59:54 2013

@author: lepse
"""

from alinea.pearl.pearl import *
from alinea.pearl.pearl_leaf import *
from alinea.echap.interfaces import pesticide_surfacic_decay
from alinea.simcycle.milne_leaf import *
from alinea.echap.interfaces import pesticide_penetrated_decay
from alinea.echap.interception_leaf import *
from alinea.echap.interfaces import pesticide_interception

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
        n.surfacic_doses={'Chlorothalonil':1,'Epoxiconazole':2,'Metconazole':3}
        n.penetrated_doses={'Chlorothalonil':1,'Epoxiconazole':2,'Metconazole':3}
        n.temp = 12
        n.rain_intensity = 0
        n.relative_humidity = 100 
        n.wetness = True
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
    return g


########################## tests decays

def test_surfacic():
    g = adel_mtg()
    g = update_on_leaves(g)
    decay_model = PearLeafDecayModel()
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
    decay_model = PearLeafDecayModel()
    g = pesticide_surfacic_decay(g, decay_model)
    decay_model = PenetratedDecayModel()
    g = pesticide_penetrated_decay(g, decay_model)
    return g

def test_no_doses():
    g = adel_mtg()
    g = update_no_doses(g)
    decay_model = PearLeafDecayModel()
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
    g = pesticide_interception(g, interception_model)
    return g

def test_intercept_no_dose():
    g = adel_mtg()
    g = update_no_doses(g)
    g = mtg_interpreter(g)
    scene = plot3d(g) 
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, interception_model)
    return g


#product_name = 'Caramba'
#dose = 200000
#productsDB = {'Ignite':{'Epoxiconazole':83}, 'Nevo':{'Chlorothalonil':375}, 'Caramba':{'Metconazole':60}}
#elevation = 90
#azimuth = 0









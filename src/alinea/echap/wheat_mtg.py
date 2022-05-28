""" Few functions to call wheat MTGs. """
# Imports #########################################################################
from alinea.adel.newmtg import *
import alinea.adel.data_samples as adel_data
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
import alinea.adel.fitting as fitting
from alinea.adel.AdelR import setAdel,RunAdel,genGeoLeaf,genGeoAxe
from math import ceil, sqrt
import random
import numpy

# Plant ###########################################################################
def adelR(nplants,dd):
    devT = adel_data.devT()
    geoLeaf = genGeoLeaf()
    geoAxe = genGeoAxe()
    pars = setAdel(devT,geoLeaf,geoAxe,nplants)
    cantable = RunAdel(dd,pars)
    return pars,cantable

def adel_mtg():
    """ create a very simple adel mtg """
    return adel_data.adel_two_metamers()

def adel_one_leaf():
    """ create a very simple adel mtg """
    d = {'plant':[1],'axe_id':['MS'],'ms_insertion':[0],'numphy':[1], 
         'Laz': [0], 'Ll' :[3], 'Lv' :[3] , 'Lsen':[0], 'L_shape':[3], 'Lw_shape':[.3], 'Linc':[0],
         'Einc':[0],'El':[0],'Ev':[0],'Esen':[0],'Ed': [0.1],'HS_final':[2],'nff':[2]}
    g=mtg_factory(d,adel_metamer, leaf_sectors=1)
    g=mtg_interpreter(g)
    return g

def adel_mtg2(nb_sect=1):
    """ create a less simple adel mtg """
    p, d = adelR(3,1000)
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,stand=[((0,0,0),0),((10,0,0),90), ((0,10,0), 0)])
    g=mtg_interpreter(g)
    return g

def adel_mtg3(nb_sect=1, d=None, p=None):
    """ create a less simple adel mtg """
    p, d = adelR(3,1000)
    if p: # nb_plants
        size = int(ceil(sqrt(p)))
        stand = numpy.array([(i, j) for i in range(size) for j in range(size)])
        numpy.random.shuffle(stand)
        stand = [((int(i)-10*size/2., int(j)-10*size/2., 0),random.randint(0,90)) for i, j in 10*stand[:p]]
    else:
        stand = [((0,0,0),0),((10,0,0),90), ((0,10,0), 0)]
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,stand=stand)
    g=mtg_interpreter(g)
    return g
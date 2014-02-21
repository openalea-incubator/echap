''' Data and parameters for pesticide models
'''

import numpy
from copy import deepcopy

# finite value used to represent infinite dt50 (d)
infinite = 1000

# k (d-1) to dt50 converter
def dt50(k):
    return numpy.log(2) / k


    
# dosage (g.l-1) of active compounds in comercial products
#---------------------------------------------------------

products_data={'Opus': {'Epoxiconazole': 125},
               'Opus_new': {'Epoxiconazole': 83},
               'Banko 500': {'Chlorothalonil': 500}}
               
# Decay parameters
#-----------------

# Milne parameters (global  activity decay and dose-response curve paramaters). Star values, Paris.ppt from Alice Milne, last slide

milne_parameters = {'Epoxiconazole': {'Ae': 3.37 * 2.0 / 9,
                                      'Ke': 4.64,
                                      'decay_rate': 1. / 20.,#d-1
                                      'dose_max_ha': 125,
                                      'Ap': 3.42 * 2.0 / 9,
                                      'Kp': 6.61,
                                      'type_code': 2},
                    'Chlorothalonil':{'Ae': 1.25 * 2.0 / 9, 
                                      'Ke': 5.84, 
                                      'decay_rate': 1. / 40, 
                                      'dose_max_ha': 1000, 
                                      'Ap': 4.03 * 2.0 / 9, 
                                      'Kp': 4.13, 
                                      'type_code': 3}
                    }

#
# Pearl parameters

# * volatilisation (found in biblio): 

# MolarMass:            MolMas(g.mol-1)
# VapourPressure:       PreVapRef (Pa)
# TemVapPre:            TemRefVap (C)
# WatSolubility:        SlbWatRef (mg.L-1)
# TemWatSol:            TemRefSlb (C)
# ThicknessLay          ThiAirBouLay (m)
# TemDif:               TemRefDif (C)
# 4.3E-5                CofDifWatRef (m2.d-1), taken as constant
# CofDifAir             CofDifAirRef (m2.d-1)

# * Transformation / Penetration
# DT50Pen should be for penetration IN THE CUTICULE, currently it is for surface to cell (To be change with Nebila)
# DT50Pen:              DT50PenCrp (d) (from surface to cuticule)
# DT50Tra:              DT50TraCrp (d) 

# * washing : 
# FacWas:               FacWasCrp (m-1)

# by default, molecular parameters = those for pure active substance, DT50 = DT50 for formulated product (Nebila)

pearl_parameters = {'Epoxiconazole_pure': {'MolarMass':329.8, 
                                      'VapourPressure':0.00001, 
                                      'TemVapPre':20, 
                                      'WatSolubility':7.1, 
                                      'TemWatSol':20, 
                                      'ThicknessLay':0.0006,
                                      'DT50Pen':3.85,
                                      'DT50Tra':infinite,
                                      'CofDifAir':0.428,
                                      'TemDif':20,
                                      'FacWas':0.0001},
                    'Epoxiconazole': {'MolarMass':329.8, 
                                      'VapourPressure':0.00001, 
                                      'TemVapPre':20, 
                                      'WatSolubility':7.1, 
                                      'TemWatSol':20, 
                                      'ThicknessLay':0.0006,
                                      'DT50Pen':1.44,
                                      'DT50Tra':infinite,
                                      'CofDifAir':0.428,
                                      'TemDif':20,
                                      'FacWas':0.0001},
                    'Chlorothalonil': {'MolarMass':265.9,
                                       'VapourPressure':0.0000762,
                                       'TemVapPre':20,
                                       'WatSolubility':0.81,
                                       'TemWatSol':20, 
                                       'ThicknessLay':0.0006,
                                       'DT50Pen':3,
                                       'DT50Tra':2.1,
                                       'CofDifAir':0.428,
                                       'TemDif':20,
                                       'FacWas':0.0001}
                    }

# additional  or re-evaluated parameters for the echap decay model

# DT50Penc (d) for penetration from cuticule to cells
# DT50Inac (d) : dt50 for inactivation of prouct penetrated in the cell
echap_decay_parameters = deepcopy(pearl_parameters)
echap_decay_parameters['Epoxiconazole'].update({'DT50inac' : dt50(0.05)})
echap_decay_parameters['Chlorothalonil'].update({'DT50inac' : dt50(0.025)})                

# vitesse avance tracteur + volume to get ro : poitevin p17 pour tartrazine, quel est leulve pour opus ? cle meme que tartrazine ?
# espacement entre buse = 
               

# application date and  application doses (l.ha-1)
# Test2001
Test2001_treat1 = """date,dose, product_name
2001-04-25 01:00:00, 0, bug
2001-04-25 02:00:00, 1, Opus
2001-05-17 09:00:00, 0, bug
2001-05-17 10:00:00, 1, Opus
2001-05-17 11:00:00, 0, bug
"""
               


                   
#granulometry : 
# Buse : type of nozzle          
# Prouit : product (tartrazine for validation, and opus pour simulation)
# Vdep : nozzle displacement (covariable, mm.s-1) : used bothe speed (5/10 or only one)
# Vitesse : drop speed (m.s-1)
# Diam : drop diameter (micro.m)
# xc : distance to the nozzle center (mm)
# Vitesse          
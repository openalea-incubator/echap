''' Data and parameters for pesticide models
'''

# dosage (g.l-1) of active compounds in comercial products

products_data={'Opus': {'Epoxiconazole': 125},
               'Opus_new': {'Epoxiconazole': 83},
               'Banko 500': {'Chlorothalonil': 500}}

# application date and  application doses (l.ha-1)
# Test2001
Test2001_treat1 = """date,dose, product_name
2001-04-25 01:00:00, 0, bug
2001-04-25 02:00:00, 1, Opus
2001-05-17 09:00:00, 0, bug
2001-05-17 10:00:00, 1, Opus
2001-05-17 11:00:00, 0, bug
"""
               
# Milne parameters (global activity decay and dose-response curve paramaters)

milne_parameters = {'Epoxiconazole': {'Ae': 0.5,
                                      'Ke': 7.01,
                                      'decay_rate': 0.06,
                                      'dose_max_ha': 125,
                                      'Ap': 0.71,
                                      'Kp': 6,
                                      'type_code': 2},
                    'Chlorothalonil':{'Ae': 1.0, 
                                      'Ke': 6.49, 
                                      'decay_rate': 0.01, 
                                      'dose_max_ha': 1000, 
                                      'Ap': 0.0, 
                                      'Kp': 0, 
                                      'type_code': 3}
                    }

# Pearl parameters (surfacic decay)
#----------------------------------

# MolarMass:            MolMas(g.mol-1)
# VapourPressure:       PreVapRef (Pa)
# TemVapPre:            TemRefVap (C)
# WatSolubility:        SlbWatRef (mg.L-1)
# TemWatSol:            TemRefSlb (C)
# ThicknessLay          ThiAirBouLay (m)
# DT50Pen:              DT50PenCrp (d)
# DT50Tra:              DT50TraCrp (d)
# FacWas:               FacWasCrp (m-1)
# TemDif:               TemRefDif (C)
# 4.3E-5                CofDifWatRef (m2.d-1), taken as constant
# CofDifAir             CofDifAirRef (m2.d-1)

# by default, molecular parameters = pure active substance, DT50 = DT50 fromulated product (Nebila)

# vitesse avance tracteur + volume to get ro : poitevin p17 pour tartrazine, quel est leulve pour opus ? cle meme que tartrazine ?
# espacement entre buse = 

pearl_parameters = {'Epoxiconazole_pure': {'MolarMass':329.8, 
                                      'VapourPressure':0.00001, 
                                      'TemVapPre':20, 
                                      'WatSolubility':7.1, 
                                      'TemWatSol':20, 
                                      'ThicknessLay':0.0006,
                                      'DT50Pen':3.85,
                                      'DT50Tra':1000,
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
                                      'DT50Tra':1000,# inifite as k almost zero (Pierre)
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
                    
#granulometry : 
# Buse : type of nozzle          
# Prouit : product (tartrazine for validation, and opus pour simulation)
# Vdep : nozzle displacement (covariable, mm.s-1) : used bothe speed (5/10 or only one)
# Vitesse : drop speed (m.s-1)
# Diam : drop diameter (micro.m)
# xc : distance to the nozzle center (mm)
# Vitesse          
''' Data and parameters for pesticide models
'''

# dosage (g.l-1) of active compounds in comercial products

products_data={'Opus': {'Epoxiconazole': 125}, 
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

pearl_parameters = {'Epoxiconazole': {'MolarMass':329.8, 
                                      'VapourPressure':0.00001, 
                                      'TemVapPre':20, 
                                      'WatSolubility':7.1, 
                                      'TemWatSol':20, 
                                      'ThicknessLay':0.0006,
                                      'DT50Pen':0.33,
                                      'DT50Tra':0.433,
                                      'CofDifAir':0.428,
                                      'TemDif':20,
                                      'FacWas':0.0001}, 
                    'Chlorothalonil': {'MolarMass':265.9,
                                       'VapourPressure':0.0000762,
                                       'TemVapPre':20,
                                       'WatSolubility':0.81,
                                       'TemWatSol':20, 
                                       'ThicknessLay':0.0006,
                                       'DT50Pen':0.14,
                                       'DT50Tra':0.23,
                                       'CofDifAir':0.428,
                                       'TemDif':20,
                                       'FacWas':0.0001}
                    }
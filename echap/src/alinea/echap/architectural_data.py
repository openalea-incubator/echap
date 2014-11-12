import pandas
import numpy

from openalea.deploy.shared_data import shared_data
import alinea.echap
from alinea.adel.plantgen.plantgen_interface import read_plantgen_inputs

from math import sqrt

import matplotlib.pyplot as plt
plt.ion()

#----------------------------------------------------- dimension
def dimensions_data():
    dim = {}
    for var in ['Mercia','Rht3','Tremie12','Tremie13']:
        if var=='Tremie12' or var=='Tremie13':
            for nff in [12,13]:
                fn = shared_data(alinea.echap, var+'_dimT%d_user_base.csv'%(nff))
                dim[var,nff] = pandas.read_csv(fn)
        elif var=='Mercia' or var=='Rht3':
            for nff in [11,12]:
                fn = shared_data(alinea.echap, var+'_dimT%d_user_base.csv'%(nff))
                dim[var,nff] = pandas.read_csv(fn)
    return dim
    
# blade dimension data from Corinne/Tino scans in 2009/2010

def blade_dimensions_MerciaRht3_2009_2010():
    """ blade dimenstion from field experiment at Grignon in 2009-2010.
        Computed by R script (Christian) from scaned leaves database
    """
    fn = shared_data(alinea.echap, 'scaned_blade_dimensions_MerciaRht3_2010.csv')
    return pandas.read_csv(fn,na_values=('NA'),sep=' ')

#----------------------------------------------------- Fitted dimension

def Mercia_2011_fitted_dimensions():
    dim = {}
    for nff in [11,12,13]:
        fn = shared_data(alinea.echap, 'Mercia_dimT%d_user.csv'%(nff))
        dim[nff] = pandas.read_csv(fn)
    return dim
    
def Rht3_2011_fitted_dimensions():
    dim = {}
    for nff in [11,12]:
        fn = shared_data(alinea.echap, 'Rht3_dimT%d_user.csv'%(nff))
        dim[nff] = pandas.read_csv(fn)
    return dim
    
def Tremie12_fitted_dimensions():
    dim = {}
    for nff in [12,13]:
        fn = shared_data(alinea.echap, 'Tremie12_dimT%d_user.csv'%(nff))
        dim[nff] = pandas.read_csv(fn)
    return dim
    
def Tremie13_fitted_dimensions():
    dim = {}
    for nff in [12,13]:
        fn = shared_data(alinea.echap, 'Tremie13_dimT%d_user.csv'%(nff))
        dim[nff] = pandas.read_csv(fn)
        # on recupere L_sheath et L_internode, W_blade, W_sheath et W_internode de Tremie12
        fn_tremie12 = shared_data(alinea.echap, 'Tremie12_dimT%d_user.csv'%(nff))
        dim_tremie12 = pandas.read_csv(fn_tremie12)
        dim[nff]['L_sheath'] = dim_tremie12['L_sheath']
        dim[nff]['L_internode'] = dim_tremie12['L_internode']
        dim[nff]['W_blade'] = dim_tremie12['W_blade']
        dim[nff]['W_sheath'] = dim_tremie12['W_sheath']
        dim[nff]['W_internode'] = dim_tremie12['W_internode']
        dim[nff]['L_internode'] = dim_tremie12['L_internode']
        
    return dim

def dimension_fits():
    d = {'Mercia': Mercia_2011_fitted_dimensions(),
        'Rht3': Rht3_2011_fitted_dimensions(),
        'Tremie12': Tremie12_fitted_dimensions(),
        'Tremie13': Tremie13_fitted_dimensions()}
    return d
    

#----------------------------------------------------- Plantgen
def plantgen_as_dict(inputs, dynT, dimT):
    d={}
    d['dynT_user'], d['dimT_user'], d['plants_number'],d['plants_density'], d['decide_child_axis_probabilities'], d['MS_leaves_number_probabilities'], d['ears_density'], d['GL_number'], d['delais_TT_stop_del_axis'], d['TT_col_break'],d['inner_params'] =  read_plantgen_inputs(inputs, dynT, dimT)
    return d
    
def GL_number():
    GL = {'Mercia': pandas.DataFrame({'TT':[485.90,1159.15,1732.10,1818.05,1959.55], '11':[3.5,3.4,3.774,2.670,1.732], '12':[3.6,3.6,4.036,3.010,1.903], '13':[4.4,3.4,4.7,3.8,2.5], 'mediane':[3.6,3.4,4.036,3.01,1.903]}),
         'Rht3': pandas.DataFrame({'TT':[485.90,1159.15,1732.10,1818.05,1959.55], '11':[3.5,3.3,3.9,2.7,1.7], '12':[3.9,3.6,4.3,3,2.2], '10':[3.3,4,2.9,None,None], 'mediane':[3.5,3.6,3.9,2.85,1.95]}),
         'Tremie12': pandas.DataFrame({'TT':[928.9,1044.25,1268.7,1392.15,1483.95,1568.7,1659.1,1797.95,1905,2084.95], '12':[1.7,2.4,4.3,4.7,4.7,4.5,4.3,3.7,2.2,0.6], '13':[1.9,2.1,4.1,5.2,4.6,4.4,4.1,3.4,2.3,1.0], 'mediane':[1.8,2.25,4.2,4.95,4.65,4.45,4.2,3.55,2.25,0.8]}),
         'Tremie13': pandas.DataFrame({'TT':[565.8,915,1048.8,1284.4,1351.4,1414.7,1545.7,1639,1762.5,1883.7], '11':[4.2,2.6,3.5,3.6,3.6,3.2,2.6,2,1.1,0.2], '12':[4.4,2.9,3.1,4.3,4.2,3.9,3.1,2,1.3,0.4], 'mediane':[4.3,2.75,3.3,3.95,3.9,3.55,2.85,2,1.2,0.3]})}
    return GL
    
def SSI_number():
    SSI = {'Mercia': pandas.DataFrame({'TT':[485.90,1159.15,1732.10,1818.05,1959.55,2108.65], '11':[0.05,5.97,7.23,8.33,9.27,11], '12':[0.05,6.31,7.96,8.99,10.10,12.00], '13':[0,7.02,8.27,9.19,10.55,13], 'mediane':[0.05,6.31,7.96,8.99,10.1,12]}),
         'Rht3': pandas.DataFrame({'TT':[485.90,1159.15,1732.10,1818.05,1959.55,2108.65], '11':[0.04,5.94,7.12,8.32,9.30,11.00], '12':[0.09,6.25,7.69,8.97,9.75,12], '10':[0.05,5.25,7.13,None,None,10], 'mediane':[0.05,5.94,7.13,8.645,9.525,11]}),
         'Tremie12': pandas.DataFrame({'TT':[928.9,1044.25,1268.7,1392.15,1483.95,1568.7,1659.1,1797.95,1905,2084.95], '12':[5.7,6.1,6.6,7.3,8,8.2,8.3,9,10.4,12.1], '13':[5.6,6.4,6.9,7.3,8.1,8.2,8.8,9.6,10.6,11.7], 'mediane':[5.65,6.25,6.75,7.3,8.05,8.2,8.55,9.3,10.5,11.9]}),
         'Tremie13': pandas.DataFrame({'TT':[565.8,915,1048.8,1284.4,1351.4,1414.7,1545.7,1639,1762.5,1883.7], '11':[0.1,5.5,6.1,7.4,7.4,7.8,8.4,9,9.9,10.8], '12':[0.1,5.8,6.9,7.7,7.8,8.1,8.9,10,10.7,11.6], 'mediane':[0.1,5.65,6.5,7.55,7.6,7.95,8.65,9.5,10.3,11.2]})}
    return SSI
    
def HS_number():
    SSI = {'Mercia': pandas.DataFrame({'TT':[485.90,1159.15,1304.10,1732.10,1818.05,1959.55,2108.65], '11':[3.6,9.4,10.5,11,11,11,11], '12':[3.6,9.9,11.3,12.0,12.0,12.0,12.0], '13':[4.4,10.4,11.9,13,13,13,13], 'mediane':[3.6,9.9,11.3,12,12,12,12]}),
         'Rht3': pandas.DataFrame({'TT':[485.90,1159.15,1304.1,1732.10,1818.05,1959.55,2108.65], '11':[3.5,9.2,10.2,11.0,11.0,11.0,11.0], '12':[4.04,9.86,11.26,12,12,12,12], '10':[3.38,9.25,10.09,10,10,10,10], 'mediane':[3.5,9.25,10.2,11,11,11,11]}),
         'Tremie12': pandas.DataFrame({'TT':[928.9,1044.25,1190.95,1268.7,1392.15,1483.95,1568.7,1659.1,1797.95,1905,2084.95], '12':[7.6,8.5,10.1,11,12,12.6,12.6,12.7,12.7,12.7,12.6], '13':[7.5,8.5,10.2,11.0,12.5,12.7,12.6,12.9,12.9,12.9,12.6], 'mediane':[7.55,8.5,10.15,11,12.25,12.65,12.6,12.8,12.8,12.8,12.6]}),
         'Tremie13': pandas.DataFrame({'TT':[565.8,915,1048.8,1284.4,1351.4,1414.7,1545.7,1639,1762.5,1883.7], '11':[4.3,8.1,9.6,11,11,11,11,11,11,11], '12':[4.5,8.8,10,12,12,12,12,12,12,12], 'mediane':[4.4,8.45,9.8,11.5,11.5,11.5,11.5,11.5,11.5,11.5]})}
    return SSI


    
'''   
GL_Mercia = {'GL_number' : {1732.1: 4.03631578947369, 1818.05:3.01047368421053,1959.55:1.90263157894737, 2108.65:0.0},
            'TT_col_break' : 0.0}   
GL_Rht3 = {'GL_number' : {1732.1: 3.88352941176471, 1818.05:2.68005882352941,1959.55:1.69764705882353, 2108.65:0.0},
            'TT_col_break' : 0.0}   
GL_Tremie = {'GL_number' : {1483.95: 4.9025, 1568.7:4.73684210526316,
             1659.1:4.16814814814815, 1797.95:3.43777777777778,
             1905:2.35888888888889,2084.95:0.85578947368421,2150:0.00},
            'TT_col_break' : 0.0}
     
#---
     
def Mercia_2010_plantgen():
    dynT = shared_data(alinea.echap, 'Mercia_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Mercia_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
    
# composite ---
def Mercia_2010_nff11_plantgen():
    dynT = shared_data(alinea.echap, 'Mercia_dynT11_user.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT11_user.csv')
    inputs = shared_data(alinea.echap, 'Mercia_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
def Mercia_2010_nff12_plantgen():
    dynT = shared_data(alinea.echap, 'Mercia_dynT12_user.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT12_user.csv')
    inputs = shared_data(alinea.echap, 'Mercia_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
def Mercia_2010_nff13_plantgen():
    dynT = shared_data(alinea.echap, 'Mercia_dynT13_user.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT13_user.csv') #pas donnees excel => ajout de la ligne 13 a dimT12
    inputs = shared_data(alinea.echap, 'Mercia_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
    
def Rht3_2010_nff11_plantgen():
    dynT = shared_data(alinea.echap, 'Rht3_dynT11_user.csv')
    dimT = shared_data(alinea.echap, 'Rht3_dimT11_user.csv')
    inputs = shared_data(alinea.echap, 'Rht3_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
def Rht3_2010_nff12_plantgen():
    dynT = shared_data(alinea.echap, 'Rht3_dynT12_user.csv')
    dimT = shared_data(alinea.echap, 'Rht3_dimT12_user.csv')
    inputs = shared_data(alinea.echap, 'Rht3_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
#---
    
def Rht3_2010_plantgen():
    dynT = shared_data(alinea.echap, 'Rht3_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Rht3_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Rht3_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)

def Tremie_2011_plantgen():
    dynT = shared_data(alinea.echap, 'Tremie1_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Tremie1_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Tremie1_plantgen_inputs_MINnew.py')
    return plantgen_as_dict(inputs, dynT, dimT)

def Tremie_2012_plantgen():
    dynT = shared_data(alinea.echap, 'Tremie2_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Tremie2_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Tremie1_plantgen_inputs_MINnew.py')
    return plantgen_as_dict(inputs, dynT, dimT)

def HS_data():
    fn = shared_data(alinea.echap, 'HS_data_Mercia_Rht3_2010_2011.csv')
    #fn = shared_data(alinea.echap, 'HS_data_Tremie1_2011_2012.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    grouped = data.groupby(['variety'],as_index=False)
    return grouped.aggregate('mean')
'''  
#
# Plot data
#
def valSD(listeVal): #calcul moyenne/ecart type
    nombre = 0
    sommeValeurs = 0.0
    sommeCarres = 0.0
    for x in listeVal:
        nombre += 1
        sommeValeurs += x
        sommeCarres += x * x
    mean = sommeValeurs / nombre
    SD = sqrt(sommeCarres / nombre - mean * mean)
    return SD

def Plot_data_Mercia_Rht3_2010_2011():
    """
    Plot data for Boigneville 2010-2011 
    
    data include :
    - estimates of plant density at emergence (id Archi12). Emergence means 75% plants emerged
    - estimates of ear density at harvest (all axes, id Archi10)
    
    Notes :
    - Missing raw data for plant counts at emergence (Archi7, used Archi12 instead)
    """
    d={'Mercia': {
        'code_date': {'sowing':'2010-10-15', 'emergence':'2010-11-02', 'harvest':'2011-06-20'},
        'TT_date': {'2010-10-15':0, '2010-11-02':143, '2011-06-20':2109},
        'sowing_density': 220, 
        'plant_density_at_emergence' : 203, 
        'raw_ear_density_at_harvest':[507,440,427,357,430,410,427,480,363, 350, 497, 410, 493, 340, 407, 467, 490, 433, 547, 450, 527, 427, 483, 493],
        'inter_row': 0.15},
        'Rht3': {
        'code_date': {'sowing':'2010-10-15', 'emergence':'2010-11-02', 'harvest':'2011-06-20'},
        'TT_date': {'2010-10-15':0, '2010-11-02':143, '2011-06-20':2109},
        'sowing_density': 220, 
        'plant_density_at_emergence' : 211, 
        'raw_ear_density_at_harvest':[440,420,330,433,347,410,387,367,360,397,357,377,427,367,380,347,330,420,423,397,383,367,377,377],
        'inter_row': 0.15}}
    for g in ('Mercia','Rht3'):
        d[g]['ear_density_at_harvest'] = numpy.mean(d[g]['raw_ear_density_at_harvest'])
        d[g]['ear_density_at_harvest_SD'] = valSD(d[g]['raw_ear_density_at_harvest'])
    return d

def Plot_data_Tremie_2011_2012():
    """
    Plot data for Boigneville 2011-2012 
    
    Notes:
    - No data found for density of emergence as this trial was followed only after winter (winter trial was on appache)
    - plant density data at date 3 (11/04/2012) and 4 comes from axe counting/ LAI measurement data taken on  5 rank * 60 cm prelevement)
    """
    d = {
    'code_date':{'sowing': '2011-10-21','emergence': '2011-11-03', 'harvest':'2012-06-19', 'd3': '2012-04-11', 'd4': '2012-05-09', 'arvalis':'2012-03-20'}, 
    'TT_date': {'2011-10-21':0, '2011-11-03':141, '2012-06-19':2135, '2012-03-20':1006, '2012-04-11':1240, '2012-05-09':1515},
    'sowing_density': 280,
    'raw_ear_density_at_harvest':[479, 490, 598, 608, 538, 503, 493, 430, 458, 437, 486, 489, 465, 406],
    'inter_row':0.143,
    'plant_density':{'2012-03-20': [238, 304, 287, 237, 290, 301, 290, 287, 273],
                     '2012-04-11': [301, 273, 315], 
                     '2012-05-09':[257, 301, 263]},
    'axe_density': {'2012-04-11': [918, 956, 979],
                    '2012-05-09':[601, 585, 506]}, 
    'fertile_axis_density':{'2012-05-09':[545, 569, 443]} 
    }
    d['ear_density_at_harvest']=numpy.mean(d['raw_ear_density_at_harvest'])
    d['mean_plant_density'] = numpy.mean(reduce(lambda x,y:x+y,d['plant_density'].values()))
    return d
    
def Plot_data_Tremie_2012_2013():
    """
    Plot data for Boigneville 2012-2013 
    
    Notes
    - no date found for plant density counts at stage 2F and epis1cm (estimation epis 1cm : 9/04/2012)
    - no date found for countings of ear density
    - plant density and ear density were estimated on 2 ranks * 1m plots at stage 2F and epis1cm (0.3 m2), and on plots of 5 ranks * 0.6 m for LAI plots (0.45 m2) 
    - *** IMPORTANT ***Plant density data measured for LAI estimation in excel files have consider 4 ranks instead of 5 (confirmed by Benjamin) (ie density and LAI should be multiplied by 0.8)
    """
    d = {
   'code_date':{'sowing': '2012-10-29','emergence': '2012-11-19', '2F':'2013-01-03', 'epis1cm': '2013-04-09', 'LAI1':'2013-04-22', 'LAI2': '2013-05-13', 'harvest': '2013-06-09'},
   'TT_date': {'2012-10-29':0, '2012-11-19':140, '2013-04-22':941, '2013-05-13':1185, '2013-01-03':421, '2013-04-09':811, '2013-06-09': 1548},
   'sowing_density': 300,
   'plant_density':{'2013-01-03': [237, 287, 217, 237, 293, 220, 253, 213, 220, 253], #'2F'
                    '2013-04-09': [203, 203, 193, 207, 223, 203, 207, 207, 200, 197], #'epis1cm'
                    '2013-04-22':[328 * 0.8, 322 * 0.8, 356 * 0.8],
                    '2013-05-13':[272* 0.8, 361 * 0.8, 350 * 0.8]},
   'raw_ear_density_at_harvest':[643, 580, 693, 813, 663, 670, 693, 763, 567, 590, 707, 637, 617, 747, 760, 693, 653, 670],
    'inter_row':0.15
    }
    d['ear_density_at_harvest'] = numpy.mean(d['raw_ear_density_at_harvest'])
    d['ear_density_at_harvest_SD'] = valSD(d['raw_ear_density_at_harvest'])
    d['mean_plant_density'] = numpy.mean(reduce(lambda x,y:x+y,d['plant_density'].values()))
    return d

def Plot_data():
    d = Plot_data_Mercia_Rht3_2010_2011()
    d.update({'Tremie12':Plot_data_Tremie_2011_2012(),
              'Tremie13':Plot_data_Tremie_2012_2013()})
    return d
    
#
# LAI data
#
def PAI_photo_data():
    d = {'Mercia': pandas.DataFrame({'TT_date':[1137,1137,1137,1137,1137,1137,1137,1137,1137,1281,1281,1796,1796], 'PAI_vert_photo':[3.11,3.18,2.93,3.39,3.58,3.38,3.39,3.22,3.24,3.94,4.12,3.82,4.27]}),
    'Rht3': pandas.DataFrame({'TT_date':[1137,1137,1137,1137,1137,1137,1137,1137,1137,1137,1281,1281,1281,1796,1796,1796,1796,1796], 'PAI_vert_photo':[4.05,3.16,3.16,3.84,3.27,3.37,3.5,3.65,3,3.03,3.94,4.03,3.9,3.32,3.51,3.52,3.19,3.69]}),
    'Tremie12': pandas.DataFrame({'TT_date':[923,923,923,923,923,923,923,923,923,923,923,923,923,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1551,1551,1551,1551,1551,1551,1551,1551], 'PAI_vert_photo':[1.09,0.8,1.02,1.09,1.3,0.88,0.74,0.93,0.87,1.08,0.92,1.22,1.23,3.64,3.45,3.8,3.79,5.69,3.97,3.84,4.35,4.31,4.85,5.15,4.18,4.57,4.52,3.84,3.98,4.2,4.22,3.95,3.99,4.54,4.39]}),
    'Tremie13': pandas.DataFrame({'TT_date':[884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580], 'PAI_vert_photo':[1.84,2,1.91,1.43,1.76,2.07,1.82,1.52,1.42,1.96,1.73,1.65,1.61,1.67,1.48,1.52,1.75,1.7,2.07,1.76,1.82,1.69,1.64,2.35,2.15,2.41,2.4,2.13,2.55,2.28,2.76,2.69,2.41,3.07,3.01,3.24,3.32,3.1,2.25,2.38,2.26,2.3,2.17,2.91,2.78,2.97,2.89,2.69,3.13,3.3,3.25,3.22,2.68,2.78,2.51,2.83,2.65,2.79,2.79,2.99,2.95,3.12,2.84,2.79,3.08,2.62,2.99,2.64,2.76,2.97,2.56,2.7,2.62,2.67,2.82,2.72,2.65,2.06,1.87,1.96,1.84,1.96,2.27,2.14,2.17,1.94,2.04,2.07,2.02,1.88,1.68,1.84,1.84,1.83,1.8,1.74,1.74,1.85,2.05,2.06,2.05,2.04]})
    }
    return d

def LAI_biomasse_data():
    d = {'Tremie12': pandas.DataFrame({'TT_date':[923,923,923,1260,1260,1260,1551,1551,1551], 'LAI_vert_biomasse':[0.97,1.2,1.14,3.61,4.90,4.42,4.46,4.74,4.08]}),
    'Tremie13': pandas.DataFrame({'TT_date':[941,941,941,1060,1060,1060], 'LAI_vert_biomasse':[(4.45*0.8),(4.59*0.8),(5.27*0.8),(5.63*0.8),(5.29*0.8),(6.55*0.8)]})
    }
    return d
    
def _nmax(sub):
    sub['nmax'] = sub['id_Feuille'].max()
    return sub
    
def scan_dep_Tremie12():
    '''
    Leaf scan data for Tremie 2011-2012
    
    Found data are:
        - Manip 'archi/silhouette': Scan des feuilles aux dates 09/03/2012 (HS 7.4)
                                                     11/04/2012 (HS 10.5)
                                                     09/05/2013 (HS 13.06)
                                files : C:\Users\echap\Desktop\data_mariem\manipeEchap\Manipe_Boigneville20112012\1.Essai_Tremie\1.MesuresINRA\0.DatesBrute\D3_11.04.2012\Scanne_11042012\
                                C:\Users\echap\Desktop\data_mariem\manipeEchap\Manipe_Boigneville20112012\1.Essai_Tremie\1.MesuresINRA\0.DatesBrute\D4_09.05.2012\scannes_09.05.2012\
                                Avec fichier de notations associe
                                C:\Users\echap\Desktop\ECHAP_ARVALIS\Architecture\Archi pour publi\Tremie 2012\2. Longueur-largeur-surface 2N et DFE_ 2012.xlsx
    
            numerotation des feuilles coherente avec photos et vrai numero determine (plantes baguees)
            donnes dispo dans share/data raw_srdb_BoignevilleTremie12.csv. 

    notes:
    -target dates of interception are 13/04 (HS 10.67) and 11/05 (HS 13.4)
    - Pour avoir les donnnes utilisees dans le fichier interception, il faut exclure les feueuille stat=3 et appliquer fraction HS a la deniere feuille emergee

    '''
    data_file = shared_data(alinea.echap, 'raw_srdb_BoignevilleTremie12.csv')
    df = pandas.read_csv(data_file, decimal=',', sep=';')
    df = df.loc[:,:'stat']

    df = df[df['id_Axe']=='MB']
    df = df.sort(['HS','id_Feuille'])
    grouped = df.groupby(['HS','N plante'], as_index=False)
    df = grouped.apply(_nmax)
    df['ntop_cur'] = df['nmax'] - df['id_Feuille'] + 1
    
    # pas de donnees observees pour la feuille 11 a HS = 10.15
    # on cree donc une copie de id_feuille=10 pour ce HS que l'on change par id_Feuille=11 et les autres colonnes vides
    # ajout ligne id_feuille 11 pour HS 10.15
    df.ntop_cur[df.HS==10.15] = df.ntop_cur+1
    
    return df
    
def scan_dep_Tremie13():
    '''
    Scan des feuilles non depouillees pour le 22/04/2012 et le 03/05/2012
    C:\Users\echap\Desktop\ECHAP_ARVALIS\Architecture\Archi Anne_2012-2013\Analyse Scan LAI\
    '''
    data_file = shared_data(alinea.echap, 'Tremie13_scan.csv')
    df = pandas.read_csv(data_file, decimal=',', sep=';')
    
    df = df[df['id_Axe']=='MB']
    df = df.sort(['HS','id_Feuille'])
    grouped = df.groupby(['HS','N plante'], as_index=False)
    df = grouped.apply(_nmax)
    df['ntop_cur'] = df['nmax'] - df['id_Feuille'] + 1
    
    return df

def treatment_scan(name='Tremie13'): #Tremie13, pas de colonne stat => stat='all'

    if name is 'Tremie12':
        df = scan_dep_Tremie12()
    elif name is 'Tremie13':
        df = scan_dep_Tremie13()

    # on prend maintenant tjrs stat 1, 2 ou 3, plus besoin du parametre stat='all'
    '''if stat is 'all':
        print 'On prend l\'ensemble des donnees => statuts 1, 2 ou 3'
    else:
        df = df[df['stat']<3]
        print 'On prend seulement les donnees ayant les statuts 1 ou 2'
    '''    
    
    dater = []; moy_feuille = []; ntop = []; A_bl = []
    lst_date = df['HS'].unique()
    for HS in lst_date:
        dt = df[df['HS']==HS]
        
        # tableau avec notation depuis le haut (ntop_cur) et moyenne notation depuis le bas (moyenne id_Feuille)
        #data = dt.groupby(['ntop_cur'], as_index=False).mean()  
        # tableau avec notation depuis le bas (id_Feuille) et moyenne notation depuis le haut (moyenne ntop_cur)     
        data = dt.groupby(['id_Feuille'], as_index=False).mean()    

        nbr = data['ntop_cur'].count(); cpt = 0
        while cpt<nbr:
            dater.append(HS)
            moy_feuille.append(data['id_Feuille'][cpt])
            ntop.append(data['ntop_cur'][cpt])
            A_bl.append(data['A_bl'][cpt])
            cpt += 1
            
    surfSet = zip(dater,moy_feuille,ntop,A_bl)
    df_fin = pandas.DataFrame(data = surfSet, columns=['HS', 'moyenne id_Feuille', 'ntop_cur', 'Area A_bl'])
    return df_fin    
    
    # barplot chaque feuille par date
    '''for date in ['mercredi 9 mai 2012']:
        dt = df[df['date prelev']==date] 
        for leaf in dt['id_Feuille'].unique():  
            plt.figure()
            dy = dt[dt['id_Feuille']==leaf] 
            #dt.plot('id_Feuille', 'A_bl', label='A_bl', kind='bar')
            dy['A_bl'].hist()
            #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
            plt.suptitle("scan "+date+", num leaf : "+str(leaf))
            #plt.xlabel("num feuille")'''
        
    # plot moyenne feuille par date
    '''dfe = df.groupby(['date prelev','id_Feuille'], as_index=False).mean()
    for date in ['lundi 2 avril 2012', 'mercredi 11 avril 2012', 'mercredi 9 mai 2012', 'vendredi 9 mars 2012']: 
        plt.figure()
        dt = dfe[dfe['date prelev']==date]
        dt.plot('id_Feuille', 'A_bl', '-or', label = 'A_bl')
        dt.plot('id_Feuille', 'A_bl_green', '-og', label = 'A_bl_green')
        #titre, legende et titre axe x
        plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
        plt.suptitle("scan "+date)
        plt.xlabel("num feuille")'''

    
#
# TC data
#   
def TC_data():
    d = {'Mercia_0':pandas.DataFrame({'TT':[475,475,475,475,475,475,475,475,475,1137,1137,1137,1137,1137,1137,1137,1137,1137,1137,1281,1281,1281,1281,1281,1281,1281,1281,1281,1796,1796,1796,1796,1796,1796,1796,1796,1796],'TC':[0.17,0.16,0.13,0.16,0.12,0.13,0.13,0.14,0.10,0.66,0.66,0.66,0.59,0.58,0.58,0.66,0.65,0.65,0.64,0.69,0.72,0.73,0.72,0.67,0.69,0.81,0.80,0.71,0.52,0.61,0.51,0.53,0.62,0.56,0.68,0.67,0.77]}),
         'Rht3_0': pandas.DataFrame({'TT':[475,475,475,475,475,475,475,475,475,1137,1137,1137,1137,1137,1137,1137,1137,1137,1137,1137,1281,1281,1281,1281,1281,1281,1281,1281,1281,1796,1796,1796,1796,1796,1796,1796,1796,1796],'TC':[0.11,0.10,0.11,0.11,0.10,0.10,0.11,0.12,0.12,0.85,0.79,0.79,0.78,0.80,0.80,0.8,0.8,0.81,0.8,0.8,0.82,0.84,0.8,0.85,0.84,0.86,0.83,0.81,0.81,0.55,0.55,0.57,0.49,0.56,0.64,0.64,0.62,0.68]}),
         'Tremie12_0': pandas.DataFrame({'TT':[923,923,923,923,923,923,923,923,923,923,923,923,923,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1551,1551,1551,1551,1551,1551,1551,1551,1551,1551,1551,1551,1551,1551],'TC':[0.27,0.24,0.23,0.29,0.30,0.28,0.30,0.31,0.30,0.29,0.28,0.33,0.32,0.36,0.58,0.51,0.37,0.49,0.56,0.62,0.54,0.53,0.62,0.57,0.49,0.69,0.72,0.36,0.58,0.51,0.37,0.49,0.56,0.62,0.54,0.53,0.62,0.57,0.49,0.69,0.72]}),
         'Tremie13_0': pandas.DataFrame({'TT':[884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580],'TC':[0.54,0.53,0.51,0.48,0.55,0.51,0.52,0.53,0.52,0.51,0.55,0.55,0.56,0.59,0.51,0.55,0.54,0.52,0.47,0.50,0.53,0.59,0.55,0.51,0.47,0.55,0.51,0.52,0.51,0.5,0.49,0.47,0.44,0.45,0.63,0.57,0.54,0.46,0.56,0.58,0.62,0.62,0.62,0.6,0.62,0.63,0.67,0.63,0.76,0.71,0.7,0.65,0.73,0.67,0.75,0.69,0.74,0.76,0.65,0.62,0.63,0.62,0.6,0.73,0.69,0.66,0.7,0.73,0.62,0.61,0.75,0.85,0.59,0.51,0.55,0.54,0.58,0.46,0.7,0.49,0.47,0.56,0.44,0.57,0.51,0.5,0.63,0.67,0.42,0.5,0.54,0.56,0.41,0.58,0.63,0.66,0.62,0.52,0.57,0.45,0.55,0.57,0.5,0.46,0.48,0.5,0.51,0.54,0.44,0.55,0.67,0.64,0.56,0.58,0.55,0.55,0.64,0.60,0.62,0.58,0.57,0.56,0.64,0.57,0.5,0.52,0.57,0.56,0.55,0.57,0.51,0.55,0.52,0.52,0.49,0.56,0.52,0.4,0.55,0.58,0.51,0.53,0.53,0.51,0.46,0.53,0.53]}),
         'Mercia_57':pandas.DataFrame({'TT':[1137,1137,1137,1137,1137,1137,1137,1137,1137,1281,1281,1796,1796],'TC':[0.94,0.95,0.93,0.96,0.96,0.96,0.96,0.95,0.95,0.97,0.98,0.97,0.98]}),
         'Rht3_57': pandas.DataFrame({'TT':[1137,1137,1137,1137,1137,1137,1137,1137,1137,1137,1281,1281,1281,1796,1796,1796,1796,1796],'TC':[0.98,0.95,0.95,0.97,0.95,0.96,0.96,0.97,0.94,0.94,0.97,0.98,0.97,0.95,0.96,0.96,0.95,0.97]}),
         'Tremie12_57': pandas.DataFrame({'TT':[923,923,923,923,923,923,923,923,923,923,923,923,923,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1260,1551,1551,1551,1551,1551,1551,1551,1551],'TC':[0.64,0.53,0.61,0.64,0.70,0.56,0.50,0.58,0.55,0.64,0.58,0.68,0.68,0.97,0.96,0.97,0.97,0.99,0.97,0.97,0.98,0.98,0.99,0.99,0.98,0.99,0.99,0.97,0.98,0.98,0.98,0.97,0.98,0.99,0.98]}),
         'Tremie13_57': pandas.DataFrame({'TT':[884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,884,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,994,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1351,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580,1580],'TC':[0.82,0.85,0.83,0.73,0.81,0.85,0.82,0.76,0.73,0.84,0.8,0.78,0.78,0.79,0.75,0.76,0.80,0.79,0.85,0.81,0.82,0.79,0.78,0.89,0.86,0.89,0.89,0.86,0.91,0.88,0.92,0.92,0.89,0.94,0.94,0.95,0.95,0.94,0.88,0.89,0.88,0.88,0.87,0.93,0.92,0.94,0.93,0.92,0.95,0.95,0.95,0.95,0.92,0.92,0.90,0.93,0.91,0.93,0.93,0.94,0.94,0.95,0.93,0.93,0.94,0.91,0.94,0.91,0.92,0.94,0.91,0.92,0.91,0.92,0.93,0.92,0.92,0.85,0.82,0.84,0.82,0.84,0.88,0.86,0.87,0.84,0.85,0.85,0.85,0.83,0.79,0.82,0.82,0.82,0.81,0.8,0.8,0.82,0.85,0.85,0.85,0.85]})}
    return d
#
# mat data
#
def mat_data():
    d = {'Mercia':pandas.DataFrame({'TT':[1137,1281,1796],'light0':[],'light20':[]}),
         'Rht3': pandas.DataFrame({'TT':[1137,1281,1796],'light0':[],'light20':[]}),
         'Tremie': pandas.DataFrame({'TT':[1260,1550],'light0':[0.177,0.045],'light20':[0.633,0.141]})}
    return d
    
#
# leaf surface cm2
#
def surface_data():
    d = {'Tremie12': pandas.DataFrame({'TT':[905.5,1160.1,1239.4,1514.6],'surface':[9.27,21.72,28.03,30.55]}),
    'Tremie13': pandas.DataFrame({'TT':[941,1096.3],'surface':[9.35,13.15]})}
    return d

#GL/HS/SSI
def dynamique_data():
    dd = shared_data(alinea.echap, 'dynamique_data.csv')
    data = pandas.read_csv(dd,decimal='.',sep=';')
    return data
    
#-------------------------------------------------------------------------------  
# Utilities for processing tillering data

def _maxna(x):
    m = x.max()
    if any(x.isnull()) and m == 0:
        m = None
    return m

def emission_probabilities(df, last='T6'):
    grouped = df.groupby('N',as_index=False)
    em = grouped.agg(_maxna)
    # em = em.reset_index()
    s = em.ix[:,'TC':last].sum()
    n = em.ix[:,'TC':last].apply(lambda x: x.dropna().count())
    probas = s / n
    return probas.to_dict()
    
def plant_viability(df):
    grouped = df.groupby('Date', as_index=False)
    res = grouped.apply(lambda x: x['MB'].count() * 1.0 / len(x['MB']))
    return {'Date':res.index.tolist(), 'viability':res.tolist()}
  
def axis_dynamics(df):
    grouped = df.groupby('Date')
    s = grouped.agg('sum').ix[:,'TP':]
    n = grouped.agg(lambda x: x.apply(lambda x: x.dropna().count())).ix[:,'TP':]
    axis =  s / n
    axis = axis.replace(numpy.inf, numpy.nan)
    res = axis.to_dict('list')
    res.update({'Date':axis.index.tolist()})
    return res
    
def nff_probabilities(df):
    prop = 1. * numpy.bincount(map(int,df['Nff'].dropna())) / len(df['Nff'].dropna())
    d = dict(zip(range(len(prop)), prop))
    return {k:v for k,v in d.iteritems() if v > 0}
    
    
#-------------------------------------------------------------------------------  
# TO DO : remove ears_per_plant (should be part of reconstruction)


def Tillering_data_Mercia_Rht3_2010_2011():
    """Tillering data for Boigneville 2010-2011
    
    Data are from 36 tagged plant per cultivar, followed during the complete season.
    
    Data found in Archi11,Archi19, Archi33 were used to build a synthetic table containing :
        - presence/absence of living mainstem (column MS, NA means the plant is dead) and number of leaf emited (Nff)
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive)
        - total number of primary tillers (column TP) and/or secondary tillers (dead or alive, column TS), and/or total number of tiller present (primary or secondary, dead or alive), and/or total number of tiller with at least 2 ligulated leaves (column TT3F).
        - number of fertile tillers (column FT) at the end from the counting of tillers that are alive or that bears an ear 
        
    These data are aggregated to estimate primary emission probability, dynamics of mean number of axis present on plants, plant mortality and  number of ears per plant
        
    
    Notes :
    - No tillering data at date 4 and 5
    - when information is available at date 1 or 3, question marks at date 2 were replaced by confirmed tiller positions 
    - at date 3 , for the 2x 24 plants not measured in details, we invert column 'total primary' and 'total secondary' as it makes data much more consistent with date 2
    - at date 3 TT3F was estimated as total primary  + delta secondary
    - at date 3, on the 2 x 12 plants measured in details, it was possble to count presence/absence of living primary tillers (T1->T5). Living wes defined as Green leaves > 2 (minimal value of the GL model). 
"""
    
    fn = shared_data(alinea.echap, 'Tillering_data_Mercia_Rht3_2010_2011.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    date_code = {'d1':'2011-01-19', 'd2':'2011-04-18', 'd3':'2011-04-27', 'd4':'2011-05-26', 'd5':'2011-06-01', 'd6':'2011-06-09'}
    TT_code = {'d1':473, 'd2':1145, 'd3':1278, 'd4':1710, 'd5':1800, 'd6':1940}
    
    # infer emision at date 3 from Total primary present (TP)
    #add statistics for T6 as well as it may be present at date 3    
    def infer_d3(g):
        g.insert(13,'T6',0.0)
        g['T6'][g['MB']!=1] = numpy.nan # avoid infering T6 = 0 on dead plants
        TP = g['TP']
        date = g['Date']
        if all(TP[date==3].notnull()) and all(pandas.isnull(g.ix[date==3,'TC':'T6'])):
            infer = g.ix[date==2,'TC':'T6'] 
            if all(TP[date==3] > TP[date==2]):
                d = int(TP[date==3].values - TP[date==2].values)
                try:
                    lastT = int(max(2, max(numpy.where(infer > 0)[1])))
                except ValueError:
                    lastT = 2
                infer.ix[:,(lastT+1):(lastT + 1 + d)] = 1                
            g.ix[g['Date']==3,'TC':'T6'] = infer
        return g
        
    grouped = data.groupby(['Var', 'N'],as_index=False)
    newdata = grouped.apply(infer_d3)
    
    # compute emmission probability using notations before date 6
    edata = newdata[newdata['Date'] < 6]
    edata = edata.reset_index()
    grouped = edata.groupby('Var',as_index=False)
    emission = {k:emission_probabilities(v) for k,v in grouped}
    
    # compute nff probabilities
    grouped = data.groupby('Var',as_index=False)
    nff_prop = {k:nff_probabilities(v) for k,v in grouped}
    
    # compute ear_per_plante (including main stem) 
    eardata = newdata[newdata['Date'] >= 6]
    eardata = eardata.reset_index()
    grouped = eardata.groupby('Var',as_index=False)
    ears_per_plant = {k:  1  + (v['FT'].sum() / v['FT'].dropna().count()) for k,v in grouped}
    
    # compute plant viability
    grouped = data.groupby('Var',as_index=False)
    viability = {k:plant_viability(v) for k,v in grouped}
    
    #compute tillering dynamics
    axdyn = {k:axis_dynamics(v) for k,v in grouped}
    
    obs = {k:{'emission_probabilities': emission[k],
           'nff_probabilities': nff_prop[k],
           'ears_per_plant': ears_per_plant[k],
           'plant_survival': viability[k],
           'tillers_per_plant': axdyn[k], 
           'date_code' : date_code,
           'TT_code' : TT_code} for k in ('Mercia', 'Rht3')}
    return obs

def Tillering_data_Tremie12_2011_2012():
    """Tillering data for Boigneville 2011-2012
    
    Data come from different plants than those tagged at date 1 and date 6, from tagged plants at date 3 for presence / absence of T7 and at date 7 for counts of fertile tillers.
    Counting of axe per plant were also done at stage epis1cm on an independant set of plants, that was renamed here date 8.
    
    Data found in Archi2, Archi3, Archi14, Archi15 were :
        - presence/absence of living mainstem (column MS, NA means the plant is dead) and number of leaf emited (Nff)
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive)
        - total number of primary tillers (column TP) and/or secondary tillers (dead or alive, column TS), and/or total number of tiller present (primary or secondary, dead or alive), and/or total number of tiller with at least 2 ligulated leaves (column TT3F).  
        - number of fertile tillers (column FT) at the end from the counting of tillers that are alive or that bears an ear 
        
    These data are aggregated to estimate primary emission probability, dynamics of mean number of axis present on plants, plant mortality and  number of ears per plant
 
    Notes :
    - at Date 1, it was noted 'due to freezing of plants, it was difficult to determine main stem'
    - at date 2 there were damage due to to fly on Mb , T1, T2
    - At date 3, Mariem noted "T7 almost always present on scaned plants"
    - Emission of tiller 5 and 6 were never noted.
    - at date 7, values of fertile tiller  per plants seems buggy compared to other dates
    """
    
    fn = shared_data(alinea.echap, 'Tillering_data_Tremie1_2011_2012.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    date_code = {'d1':'2012-03-09', 'd2':'2012-04-02', 'd3':'2012-04-11', 'd4':'2012-05-09', 'd5':'2012-05-29', 'd6':'2012-06-12', 'd7':'2012-07-12', 'd8':'2012-04-04'}
    TT_code = {'d1':905, 'd2':1160, 'd3':1240, 'd4':1515, 'd5':1813, 'd6':2031, 'd7':2536, 'd8':1179}
    
    # compute emmission probability from data at date 1, and adding data at date 3 for tiller 7 only
    edata = data[data['Date'] == 1]
    edata['T7'] = data['T7'][data['Date'] == 3].values
    edata = edata.reset_index()
    emission = emission_probabilities(edata,last='T7')
    
    # compute nff probabilities
    nff_prop = nff_probabilities(data)

    # tiller dynamics
    axdyn = axis_dynamics(data)
    
    # compute ear_per_plant using data at date 6 and 7 and plot data of fertile tillers at date 4
    eardata = pandas.DataFrame(axdyn).ix[:,('Date', 'FT')].dropna()
    pdata = Plot_data_Tremie_2011_2012()
    ftd4 = 1.*numpy.array(pdata['fertile_axis_density']['2012-05-09'])  / numpy.array(pdata['plant_density']['2012-05-09']) - 1 #remove main stem to get FT count 
    # ear per plante at date 7 very strange (twice the value at other dates and non-existence of unfertile tillers): ears_per_plant taken as mean of counting at date 4 and 6
    ears_per_plant = 1 + numpy.array(ftd4.tolist() + eardata['FT'][eardata['Date'] == 6].tolist()).mean()

    
    obs = {'emission_probabilities': emission,
           'nff_probabilities': nff_prop,
           'ears_per_plant': ears_per_plant,
           'tillers_per_plant': axdyn, 
           'date_code' : date_code,
           'TT_code' : TT_code}
    return obs
    
def Tillering_data_Tremie13_2012_2013():
    """Tillering data for Boigneville 2012-2013
    
    Data for presence/absence of primary tillers come from tagged plants. Data at date 3 are from other other plants.
    
    Data found in were :
        - presence/absence of living mainstem (column MS, NA means the plant is dead) and number of leaf emited (Nff)
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive)
        - total number of primary tillers (column TP) and/or secondary tillers (dead or alive, column TS), and/or total number of tiller present (primary or secondary, dead or alive), and/or total number of tiller with at least 2 ligulated leaves (column TT3F). 
        
    These data are aggregated to estimate primary emission probability, dynamics of mean number of axis present on plants

    Notes :
    - No data found  for direct counting of fertile tillers per plant on tagged plants) 
    - data from date 3 are to be included
    """
    
    fn = shared_data(alinea.echap, 'Tillering_data_Tremie_2012_2013.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    date_code = {'d1':'2013-02-13', 'd2':'2013-03-29', 'd3': '2012-04-19'}
    TT_code = {'d1':566, 'd2':739, 'd3':915}
    
    # compute emmission probability from data on tagged plants at date 1 and 2
    edata = data[data['Date'] <= 2]
    edata = edata.reset_index()
    emission = emission_probabilities(edata,last='T5')
    
    # compute nff probabilities
    nff_prop = nff_probabilities(data)
    
    # tiller dynamics
    axdyn = axis_dynamics(data)

    obs = {'emission_probabilities': emission,
           'nff_probabilities': nff_prop,
           'tillers_per_plant': axdyn, 
           'date_code' : date_code,
           'TT_code' : TT_code}

    return obs
 

def Tillering_data():
    d = Tillering_data_Mercia_Rht3_2010_2011()
    d.update({'Tremie12': Tillering_data_Tremie12_2011_2012(),
             'Tremie13': Tillering_data_Tremie13_2012_2013()})
    return d
    
#-------------------------------------------------------------------------------  
# Angles / formes a plat  

def mean_nff():
    def _nff(ms_nff_probas):
        return sum([int(k)*v for k,v in ms_nff_probas.iteritems()])
    return {k:_nff(v['nff_probabilities']) for k,v in Tillering_data().iteritems()}
   
   
def leaf_curvature_data(name='Mercia'):

    def xy_reader(file):
        header_row_xydb = ['variety','variety_code','harvest','plant','rank','ranktop','relative_ranktop','HS','inerv','x','y']
        return pandas.read_csv(file, names=header_row_xydb, sep=',', index_col=False, skiprows=1, decimal='.')
        
    def sr_reader(file) :
        header_row_srdb = ['rankclass','s','r']
        return pandas.read_csv(file, names=header_row_srdb, sep=',', index_col=False, skiprows=1, decimal='.')
        
    if name is 'Mercia':
        data_file_xydb = shared_data(alinea.echap, 'xydb_GrignonMercia2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv') 
        # XY
        dfxy = xy_reader(data_file_xydb)
        # use Mercia + Rht3 
        #dfxy = dfxy[dfxy['variety'].isin(['Mercia','Rht3'])]
        #dfxy = dfxy.reset_index()
        # SR
        dfsr = sr_reader(data_file_srdb)

    if name is 'Rht3':
        data_file_xydb = shared_data(alinea.echap, 'xydb_GrignonMercia2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv') 
        # XY
        dfxy = xy_reader(data_file_xydb)
        #dfxy = dfxy[dfxy['variety'].isin(['Mercia','Rht3'])]
        #dfxy = dfxy.reset_index()
        # SR : use same as Mercia ???
        dfsr = sr_reader(data_file_srdb)
        
    if name is 'Tremie12':
        # fichier angle non dispo encore, on prend mercia en attendant
        data_file_xydb = shared_data(alinea.echap, 'xydb_GrignonMercia2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv')
        # XY
        dfxy = xy_reader(data_file_xydb)
        #dfxy = dfxy[dfxy['variety'].isin(['Mercia','Rht3'])]
        #dfxy = dfxy.reset_index()
        # SR
        dfsr = sr_reader(data_file_srdb)
        
    if name is 'Tremie13':
        # fichier angle non dispo encore, on prend mercia en attendant
        data_file_xydb = shared_data(alinea.echap, 'xydb_GrignonMercia2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv')
        # XY
        dfxy = xy_reader(data_file_xydb)
        #dfxy = dfxy[dfxy['variety'].isin(['Mercia','Rht3'])]
        #dfxy = dfxy.reset_index()
        # SR
        dfsr = sr_reader(data_file_srdb)
    
    return dfxy, dfsr
    

#
# Elaborated data
#    

# add mean dynamic point to HS, SSI, GL data using nff proba from tillering

def _meanDyn(name):
    df_GL = GL_number()[name]; df_HS = HS_number()[name]; df_SSI = SSI_number()[name]
    nff_proba = Tillering_data()[name]['nff_probabilities']
    df_GL['mean_pond'] = sum([df_GL[str(k)] * nff_proba[k] for k in nff_proba])
    df_HS['mean_pond'] = sum([df_HS[str(k)] * nff_proba[k] for k in nff_proba])
    df_SSI['mean_pond'] = sum([df_SSI[str(k)] * nff_proba[k] for k in nff_proba])
    return {'GL':df_GL, 'HS':df_HS, 'SSI':df_SSI}
    
def HS_GL_SSI_data():
    d={}
    for name in ('Mercia','Rht3', 'Tremie12', 'Tremie13'):
        d[name] = _meanDyn(name)
    return d
    
def _add_ghs(df, g):
    #hs_conv = HS_converter[g]
    df['Var'] = g
    #df['HS'] = hs_conv(df['TT'])
    df = df.sort('TT')
    return df

def PlantDensity():
    # Mercia /rht3
    ld = []
    for g in ('Mercia','Rht3'):
        pdata = Plot_data_Mercia_Rht3_2010_2011()[g]
        tdata = Tillering_data_Mercia_Rht3_2010_2011()[g]
        date,TT,density,SD = [],[],[],[]
        events = ['sowing', 'emergence', 'harvest']
        density = [pdata['sowing_density'], 
                   pdata['plant_density_at_emergence'],
                   pdata['ear_density_at_harvest'] / tdata['ears_per_plant']]
        ear_data = numpy.array(pdata['raw_ear_density_at_harvest']) / tdata['ears_per_plant']
        SD = [0,0,valSD(ear_data)]           
        for w in events:
            date.append(pdata['code_date'][w])
            TT.append(pdata['TT_date'][pdata['code_date'][w]])
        # add survival data
        for i,d in enumerate(tdata['plant_survival']['Date']):
            lab = 'd%d'%(d)
            events.append(lab)
            date.append(tdata['date_code'][lab])
            TT.append(tdata['TT_code'][lab])
            density.append(pdata['plant_density_at_emergence'] * tdata['plant_survival']['viability'][i])
            SD.append(0)
        df = pandas.DataFrame({'event': events, 'date':date, 'TT':TT, 'density':density})
        df['SD']=SD
        df = _add_ghs(df, g)
        ld.append(df)
    # Tremie12
    pdata = Plot_data_Tremie_2011_2012()
    tdata = Tillering_data_Tremie12_2011_2012()
    #hs_conv = HS_converter['Tremie12']
    date,TT,density,SD = [],[],[],[]
    events = ['sowing', 'harvest', 'd3', 'd4', 'arvalis']
    density = [pdata['sowing_density'], 
               pdata['ear_density_at_harvest'] / tdata['ears_per_plant'],
               numpy.mean(pdata['plant_density'][pdata['code_date']['d3']]),
               numpy.mean(pdata['plant_density'][pdata['code_date']['d4']]),
               numpy.mean(pdata['plant_density'][pdata['code_date']['arvalis']])
               ]
    ear_density_new=[]
    for i,item in enumerate(pdata['raw_ear_density_at_harvest']):
        lab = pdata['raw_ear_density_at_harvest'][i]/tdata['ears_per_plant']
        ear_density_new.append(lab)
    SD = [0, valSD(ear_density_new),
         valSD(pdata['plant_density'][pdata['code_date']['d3']]), 
         valSD(pdata['plant_density'][pdata['code_date']['d4']]), 
         valSD(pdata['plant_density'][pdata['code_date']['arvalis']])
         ]
    for w in events:
        date.append(pdata['code_date'][w])
        TT.append(pdata['TT_date'][pdata['code_date'][w]])
    df = pandas.DataFrame({'event': events, 'date' :date, 'TT':TT, 'density':density})
    df['SD']=SD
    df = _add_ghs(df, 'Tremie12')
    ld.append(df)
    # Tremie13
    pdata = Plot_data_Tremie_2012_2013()
    tdata = Tillering_data_Tremie13_2012_2013()
    #hs_conv = HS_converter['Tremie13']
    date,TT,density,SD = [],[],[],[]
    events = ['sowing', '2F', 'epis1cm', 'LAI1', 'LAI2']
    density = [pdata['sowing_density'], 
               numpy.mean(pdata['plant_density'][pdata['code_date']['2F']]),
               numpy.mean(pdata['plant_density'][pdata['code_date']['epis1cm']]),
               numpy.mean(pdata['plant_density'][pdata['code_date']['LAI1']]),
               numpy.mean(pdata['plant_density'][pdata['code_date']['LAI2']])
               ]
    SD = [0,
         valSD(pdata['plant_density'][pdata['code_date']['2F']]), 
         valSD(pdata['plant_density'][pdata['code_date']['epis1cm']]), 
         valSD(pdata['plant_density'][pdata['code_date']['LAI1']]),
         valSD(pdata['plant_density'][pdata['code_date']['LAI2']])
         ]
    for w in events:
        date.append(pdata['code_date'][w])
        TT.append(pdata['TT_date'][pdata['code_date'][w]])
    df = pandas.DataFrame({'event': events, 'date' :date, 'TT':TT, 'density':density})
    df['SD']=SD
    df = _add_ghs(df, 'Tremie13')
    ld.append(df)
    return reduce(lambda x,y : pandas.concat([x,y]), ld)

# converter axis_dyn output -> df
def _axdyn_df(tdb):
    date, TT = [], []
    for i,d in enumerate(tdb['tillers_per_plant']['Date']):
        lab = 'd%d'%(d)
        date.append(tdb['date_code'][lab])
        TT.append(tdb['TT_code'][lab])
    til = tdb['tillers_per_plant']
    return pandas.DataFrame({'date': date, 'TT': TT, 'TP' : til['TP'], 'TS':til['TS'], 'TPS': til['TT'], 'FT':til['FT'], 'TT3F':til['TT3F']})

def tillers_per_plant():
    """ Synthetic tables for tiller counting data
    """
    ld=[]
    for g in ('Mercia','Rht3'):
        tdb = Tillering_data_Mercia_Rht3_2010_2011()[g]
        df = _axdyn_df(tdb)
        df = _add_ghs(df, g)
        ld.append(df)
    #
    tdb =  Tillering_data_Tremie12_2011_2012()
    df = _axdyn_df(tdb)
    df = df[~df['TPS'].isnull()]#remove d3 where no data are available
    # add axe counts from plot data
    pdata = Plot_data_Tremie_2011_2012()
    axd = pdata['axe_density']
    faxd = pdata['fertile_axis_density']
    dfp= pandas.DataFrame({'date' : axd.keys(), 
                'TT': [pdata['TT_date'][k] for k in axd],
                'TPS': [numpy.mean(axd[k]) / numpy.mean(pdata['plant_density'][k]) - 1 for k in axd]})
    dfp['FT'] = numpy.nan
    d = faxd.keys()[0]
    dfp.ix[dfp['date']==d,'FT'] = numpy.mean(faxd[d])  / numpy.mean(pdata['plant_density'][d]) - 1
    df = pandas.concat([df, dfp])
    df = _add_ghs(df, 'Tremie12')
    ld.append(df)
    #
    tdb = Tillering_data_Tremie13_2012_2013()
    df = _axdyn_df(tdb)
    # add ears per plant, to get a point after regression
    pdata = Plot_data_Tremie_2012_2013()
    dfp= pandas.DataFrame({'date' :pdata['code_date']['harvest'],
                           'TT': pdata['TT_date'][pdata['code_date']['harvest']],
                           'FT': pdata['ear_density_at_harvest'] / pdata['mean_plant_density'] - 1}, index=[0])
    df = pandas.concat([df, dfp])
    df = _add_ghs(df, 'Tremie13')
    ld.append(df)    
    return reduce(lambda x,y : pandas.concat([x,y]), ld)
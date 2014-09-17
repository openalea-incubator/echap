import pandas
import numpy

from openalea.deploy.shared_data import shared_data
import alinea.echap
from alinea.adel.plantgen.plantgen_interface import read_plantgen_inputs



#----------------------------------------------------- Plantgen
def plantgen_as_dict(inputs, dynT, dimT):
    d={}
    d['dynT_user'], d['dimT_user'], d['plants_number'],d['plants_density'], d['decide_child_axis_probabilities'], d['MS_leaves_number_probabilities'], d['ears_density'], d['GL_number'], d['delais_TT_stop_del_axis'], d['TT_col_break'],d['inner_params'] =  read_plantgen_inputs(inputs, dynT, dimT)
    return d
    
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
'''
def Tremie_2012_plantgen():
    #dynT = shared_data(alinea.echap, 'Tremie2_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Tremie2_dimT_user.csv')
    #inputs = shared_data(alinea.echap, 'Tremie2_plantgen_inputs_MINnew.py')
    return plantgen_as_dict(inputs, dynT, dimT)
'''   
def HS_data():
    fn = shared_data(alinea.echap, 'HS_data_Mercia_Rht3_2010_2011.csv')
    #fn = shared_data(alinea.echap, 'HS_data_Tremie1_2011_2012.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    grouped = data.groupby(['variety'],as_index=False)
    return grouped.aggregate('mean')
#
# Plot data
#
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
        'ear_density_at_harvest' : 444,
        'raw_ear_density_at_harvest':[507,440,427,357,430,410,427,480,363, 350, 497, 410, 493, 340, 407, 467, 490, 433, 547, 450, 527, 427, 483, 493],
        'inter_row': 0.15},
        'Rht3': {
        'code_date': {'sowing':'2010-10-15', 'emergence':'2010-11-02', 'harvest':'2011-06-20'},
        'TT_date': {'2010-10-15':0, '2010-11-02':143, '2011-06-20':2109},
        'sowing_density': 220, 
        'plant_density_at_emergence' : 211, 
        'ear_density_at_harvest' : 384,
        'raw_ear_density_at_harvest':[440,420,330,433,347,410,387,367,360,397,357,377,427,367,380,347,330,420,423,397,383,367,377,377],
        'inter_row': 0.15}}
    return d

def Plot_data_Tremie_2011_2012():
    """
    Plot data for Boigneville 2011-2012 
    
    Notes:
    - No data found for density of emergence as this trial was followed only after winter (winter trial was on appache)
    - plant density data at date 3 (11/04/2012) and 4 comes from axe counting/ LAI measurement data taken on  5 rank * 60 cm prelevement)
    """
    d = {
    'code_date':{'sowing': '2011-10-21','emergence': '2011-11-03', 'harvest':'2012-06-19'}, 
    'TT_date': {'2011-10-21':0, '2011-11-03':141, '2012-06-19':2135, '2012-03-20':1006, '2012-04-11':1240, '2012-05-09':1515},
    'sowing_density': 280,
    'ear_density_at_harvest': 492,
    'raw_ear_density_at_harvest':[479, 490, 598, 608, 538, 503, 493, 430, 458, 437, 486, 489, 465, 406],
    'inter_row':0.143,
    'plant_density':{'2012-03-20': [238, 304, 287, 237, 290, 301, 290, 287, 273],
                     '2012-04-11': [301, 273, 315], 
                     '2012-05-09':[257, 301, 263]},
    'axe_density': {'2012-04-11': [918, 956, 979],
                    '2012-05-09':[601, 585, 506]}, 
    'fertile_axis_density':{'2012-05-09':[545, 569, 443]} 
    }
    return d
    
def Plot_data_Tremie_2012_2013():
    """
    Plot data for Boigneville 2012-2013 
    
    Notes
    - no date found for plant density counts at stage 2F and epis1cm (estimation epis 1cm : 9/04/2012)
    - no date found for countings of ear density
    - plant density and ear density were estimated on 2 ranks * 1m plots at stage 2F and epis1cm (0.3 m2), and on plots of 4 ranks * 0.6 m for LAI plots (0.36 m2) 
    - *** IMPORTANT ***Plant density data measured for LAI estimation seems bugy and more compatible with a 5 rank * 0.6 m plots dimension (ie density and LAI should be multiplied by 0.8)
    """
    d = {
   'code_date':{'sowing': '2012-10-29','emergence': '2012-11-19'},
   'TT_date': {'2012-10-29':0, '2012-11-19':140, '2013-04-22':941, '2013-05-13':1185, '2013-05-22':1284, '2013-04-09':811},
   'sowing_density': 300,
   'plant_density':{'2013-05-22': [237, 287, 217, 237, 293, 220, 253, 213, 220, 253], #'2F'
                    '2013-04-09': [203, 203, 193, 207, 223, 203, 207, 207, 200, 197], #'epis1cm'
                    '2013-04-22':[328 * 0.8, 322 * 0.8, 356 * 0.8],
                    '2013-05-13':[272* 0.8, 361 * 0.8, 350 * 0.8]},
   'ear_density_at_harvest' : 676,
   'raw_ear_density_at_harvest':[643, 580, 693, 813, 663, 670, 693, 763, 567, 590, 707, 637, 617, 747, 760, 693, 653, 670],
    'inter_row':0.15
    }
    return d
   
#
# LAI data
#
def PAI_data():
    '''
    d = {'Mercia':pandas.DataFrame({'TT':[1137,1281,1796],'PAI':[3.27,4.03,4.05]}),
         'Rht3': pandas.DataFrame({'TT':[1137,1281,1796],'PAI':[3.40,3.96,3.45]}),
         'Tremie': pandas.DataFrame({'TT':[923,1260,1550],'PAI':[1.1,4.31,4.43]})}
         #'Tremie': pandas.DataFrame({'TT':[1260,1550],'PAI':[3.8,4.7]})} #Publi Corinne
    '''
    
    d={'Mercia': {
        'code_date': {'2011-04-18', '2011-04-27', '2011-06-01'},
        'TT_date': {'2011-04-18':1137, '2011-04-27':1281, '2011-06-01':1796},
        'PAI_vert_moy_photo':{'2011-04-18':3.27, '2011-04-27':4.03, '2011-06-01':4.05},
        'PAI_vert_med_photo':{'2011-04-18':3.24, '2011-04-27':4.03, '2011-06-01':4.05},
        'PAI_vert_nbr_plt_photo':{'2011-04-18':9, '2011-04-27':2, '2011-06-01':2},
        'PAI_vert_ET_photo':{'2011-04-18':0.19, '2011-04-27':0.13, '2011-06-01':0.32} #ET=ecartype
        },
        'Rht3': {
        'code_date': {'2011-04-18', '2011-04-27', '2011-06-01'},
        'TT_date': {'2011-04-18':1137, '2011-04-27':1281, '2011-06-01':1796},
        'PAI_vert_moy_photo':{'2011-04-18':3.40, '2011-04-27':3.96, '2011-06-01':3.45},
        'PAI_vert_med_photo':{'2011-04-18':3.32, '2011-04-27':3.94, '2011-06-01':3.51},
        'PAI_vert_nbr_plt_photo':{'2011-04-18':10, '2011-04-27':3, '2011-06-01':5},
        'PAI_vert_ET_photo':{'2011-04-18':0.35, '2011-04-27':0.07, '2011-06-01':0.19}
        },
        'Tremie12': {
        'code_date': {'2012-03-09', '2012-04-11', '2012-05-09'},
        'TT_date': {'2012-03-09':923, '2012-04-11':1260, '2012-05-09':1551},
        'PAI_vert_moy_photo':{'2012-03-09':1.01, '2012-04-11':4.29, '2012-05-09':4.14},
        'PAI_vert_med_photo':{'2012-03-09':1.02, '2012-04-11':4.25, '2012-05-09':4.10},
        'PAI_vert_nbr_plt_photo':{'2012-03-09':13, '2012-04-11':14, '2012-05-09':8},
        'PAI_vert_ET_photo':{'2012-03-09':0.17, '2012-04-11':0.63, '2012-05-09':0.24},
        'LAI_vert_placette_biomasse':{'2012-03-09':0.43, '2012-04-11':0.43, '2012-05-09':0.43}, #surface placette en m2
        'LAI_vert_moy_biomasse':{'2012-03-09':1.10, '2012-04-11':4.31, '2012-05-09':4.43},
        'LAI_vert_med_biomasse':{'2012-03-09':1.14, '2012-04-11':4.66, '2012-05-09':4.41},
        'LAI_vert_nbr_plt_biomasse':{'2012-03-09':3, '2012-04-11':3, '2012-05-09':3},
        'LAI_vert_ET_biomasse':{'2012-03-09':0.1, '2012-04-11':0.6, '2012-05-09':0.3},
        'LAI_tot_moy_biomasse':{'2012-03-09':1.43}
        },
        'Tremie13': {
        'code_date': {'2013-04-17', '2013-04-26', '2013-05-29', '2013-06-13'},
        'TT_date': {'2013-04-17':918, '2013-04-26':1018, '2013-05-29':1377, '2013-06-13':1612},
        'PAI_vert_moy_photo':{'2013-04-17':1.7, '2013-04-26':2.75, '2013-05-29':2.73, '2013-06-13':1.88},
        'PAI_vert_med_photo':{'2013-04-17':1.69, '2013-04-26':2.84, '2013-05-29':2.69, '2013-06-13':1.84},
        'PAI_vert_nbr_plt_photo':{'2013-04-17':11, '2013-04-26':14, '2013-05-29':12, '2013-06-13':13},
        'PAI_vert_ET_photo':{'2013-04-17':0.16, '2013-04-26':0.41, '2013-05-29':0.14, '2013-06-13':0.13}
        }
      }
    return d

    
#
# TC data
#
def TC_data():
    d = {'Mercia_0':pandas.DataFrame({'TT':[475,1137,1281,1796],'TC':[0.14,0.63,0.73,0.61],'1-TC':[0.86,0.37,0.27,0.39]}),
         'Rht3_0': pandas.DataFrame({'TT':[475,1137,1281,1796],'TC':[0.11,0.80,0.83,0.59],'1-TC':[0.89,0.2,0.17,0.41]}),
         'Tremie_0': pandas.DataFrame({'TT':[923,1260,1550],'TC':[0.3,0.58,0.63],'1-TC':[0.7,0.42,0.37]}),
         'Mercia_57':pandas.DataFrame({'TT':[1137,1281,1796],'TC':[0.95,0.98,0.98]}),
         'Rht3_57': pandas.DataFrame({'TT':[1137,1281,1796],'TC':[0.956,0.975,0.96]}),
         'Tremie_57': pandas.DataFrame({'TT':[923,1260,1550],'TC':[0.59,0.97,0.985],'1-TC':[0.41,0.03,0.015]})}
    return d
      
#
# mat data
#
def mat_data():
    d = {'Mercia':pandas.DataFrame({'TT':[1137,1281,1796],'light0':[],'light20':[]}),
         'Rht3': pandas.DataFrame({'TT':[1137,1281,1796],'light0':[],'light20':[]}),
         'Tremie': pandas.DataFrame({'TT':[1260,1550],'light0':[0.177,0.045],'light20':[0.633,0.141]})}
    return d

#GL/HS/SSI
def dynamique_data():
    d={'Mercia': {
        'code_date': {'2011-01-19', '2011-04-18', '2011-04-27', '2011-05-26', '2011-06-01', '2011-06-09'},
        'TT_date': {'2011-01-19':485.9, '2011-04-18':1159.15, '2011-04-27':1304.10, '2011-05-26':1732.10, '2011-06-01':1818.05, '2011-06-09':1959.55},
        'HS_11_data':{'2011-01-19':3.6, '2011-04-18':9.4, '2011-04-27':10.5, '2011-05-26':11, '2011-06-01':11, '2011-06-09':11},
        'HS_12_data':{'2011-01-19':3.6, '2011-04-18':9.9, '2011-04-27':11.3, '2011-05-26':12, '2011-06-01':12, '2011-06-09':12},
        'HS_13_data':{'2011-01-19':4.4, '2011-04-18':10.4, '2011-04-27':11.9, '2011-05-26':13, '2011-06-01':13, '2011-06-09':13},
        'SSI_11_data':{'2011-01-19':0.05, '2011-04-18':5.97, '2011-05-26':7.23, '2011-06-01':8.33, '2011-06-09':9.27},
        'SSI_12_data':{'2011-01-19':0.05, '2011-04-18':6.31, '2011-05-26':7.96, '2011-06-01':8.99, '2011-06-09':10.10},
        'SSI_13_data':{'2011-01-19':0.00, '2011-04-18':7.02, '2011-05-26':8.27, '2011-06-01':9.19, '2011-06-09':10.55},
        'GL_11_data':{'2011-01-19':3.5, '2011-04-18':3.4, '2011-05-26':3.774, '2011-06-01':2.670, '2011-06-09':1.732},
        'GL_12_data':{'2011-01-19':3.6, '2011-04-18':3.6, '2011-05-26':4.036, '2011-06-01':3.010, '2011-06-09':1.903},
        'GL_13_data':{'2011-01-19':4.4, '2011-04-18':3.4, '2011-05-26':4.7, '2011-06-01':3.8, '2011-06-09':2.5}
        },
        'Rht3': {
        'code_date': {'2011-01-19', '2011-04-18', '2011-04-27', '2011-05-26', '2011-06-01', '2011-06-09'},
        'TT_date': {'2011-01-19':485.9, '2011-04-18':1159.15, '2011-04-27':1304.10, '2011-05-26':1732.10, '2011-06-01':1818.05, '2011-06-09':1959.55},
        'HS_11_data':{'2011-01-19':3.5, '2011-04-18':9.2, '2011-04-27':10.2, '2011-05-26':11, '2011-06-01':11, '2011-06-09':11},
        'HS_12_data':{'2011-01-19':4.04, '2011-04-18':9.86, '2011-04-27':11.26, '2011-05-26':12, '2011-06-01':12, '2011-06-09':12},
        'HS_13_data':{'2011-01-19':3.38, '2011-04-18':9.25, '2011-04-27':10.09, '2011-05-26':10, '2011-06-01':10, '2011-06-09':10},
        'SSI_11_data':{'2011-01-19':0.04, '2011-04-18':5.94, '2011-05-26':7.12, '2011-06-01':8.32, '2011-06-09':9.30},
        'SSI_12_data':{'2011-01-19':0.09, '2011-04-18':6.25, '2011-05-26':7.69, '2011-06-01':8.97, '2011-06-09':9.75},
        'SSI_13_data':{'2011-01-19':0.05, '2011-04-18':5.25, '2011-05-26':7.13},
        'GL_11_data':{'2011-01-19':3.5, '2011-04-18':3.3, '2011-05-26':3.9, '2011-06-01':2.7, '2011-06-09':1.7},
        'GL_12_data':{'2011-01-19':3.9, '2011-04-18':3.6, '2011-05-26':4.3, '2011-06-01':3.0, '2011-06-09':2.2},
        'GL_13_data':{'2011-01-19':3.3, '2011-04-18':4.0, '2011-05-26':2.9}
        },
        'Tremie12': {
        'code_date': {'2012-03-09', '2012-03-21', '2012-04-02', '2012-04-11', '2012-04-26', '2012-05-03', '2012-05-09', '2012-05-16', '2012-05-25', '2012-05-31', '2012-06-12'},
        'TT_date': {'2012-03-09':928.9, '2012-03-21':1044.25, '2012-04-02':1190.95, '2012-04-11':1268.7, '2012-04-26':1392.15, '2012-05-03':1483.95, '2012-05-09':1568.7, '2012-05-16':1659.1, '2012-05-25':1797.95, '2012-05-31':1905, '2012-06-12':2084.95},
        'HS_12_data':{'2012-03-09':7.6, '2012-03-21':8.5, '2012-04-02':10.1, '2012-04-11':11.0, '2012-04-26':12.0, '2012-05-03':12.6, '2012-05-09':12.6, '2012-05-16':12.7, '2012-05-25':12.7, '2012-05-31':12.7, '2012-06-12':12.6},
        'HS_13_data':{'2012-03-09':7.5, '2012-03-21':8.5, '2012-04-02':10.2, '2012-04-11':11.0, '2012-04-26':12.5, '2012-05-03':12.7, '2012-05-09':12.6, '2012-05-16':12.9, '2012-05-25':12.9, '2012-05-31':12.9, '2012-06-12':12.6},
        'SSI_12_data':{'2012-03-09':5.7, '2012-03-21':6.1, '2012-04-11':6.6, '2012-04-26':7.3, '2012-05-03':8.0, '2012-05-09':8.2, '2012-05-16':8.3, '2012-05-25':9.0, '2012-05-31':10.4, '2012-06-12':12.1},
        'SSI_13_data':{'2012-03-09':5.6, '2012-03-21':6.4, '2012-04-11':6.9, '2012-04-26':7.3, '2012-05-03':8.1, '2012-05-09':8.2, '2012-05-16':8.8, '2012-05-25':9.6, '2012-05-31':10.6, '2012-06-12':11.7},
        'GL_12_data':{'2012-03-09':1.7, '2012-03-21':2.4, '2012-04-11':4.3, '2012-04-26':4.7, '2012-05-03':4.7, '2012-05-09':4.5, '2012-05-16':4.3, '2012-05-25':3.7, '2012-05-31':2.2, '2012-06-12':0.6},
        'GL_13_data':{'2012-03-09':1.9, '2012-03-21':2.1, '2012-04-11':4.1, '2012-04-26':5.2, '2012-05-03':4.6, '2012-05-09':4.4, '2012-05-16':4.1, '2012-05-25':3.4, '2012-05-31':2.3, '2012-06-12':1.0}
        },
        'Tremie13': {
        'code_date': {},
        'TT_date': {},
        'HS_data':{},
        'SSI_data':{},
        'GL_data':{}
        }
      }
    return d
    
    
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
    - when information is available at date 1 or 3, question marks at date 2 were replaced by confirmed tiller positions """
    
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
        if all(TP[date==3].notnull()):
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

def Tillering_data_Tremie1_2011_2012():
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
    - Emission of tiller 5 and 6 were ever noted.
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
    
def Tillering_data_Tremie_2012_2013():
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
    
#-------------------------------------------------------------------------------  
# Angles / formes a plat  

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
        
    if name is 'Tremie':
        # fichier angle non dispo encore, on prend mercia en attendant
        data_file_xydb = shared_data(alinea.echap, 'xydb_GrignonMercia2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv')
        # XY
        #dfxy = xy_reader(data_file_xydb)
        #dfxy = dfxy[dfxy['variety'].isin(['Mercia','Rht3'])]
        dfxy = dfxy.reset_index()
        # SR
        dfsr = sr_reader(data_file_srdb)
        
    # add Lindex
    dfxy['Lindex'] = 1
    dfxy['Lindex'][dfxy['ranktop'] <= 4] = 2
    dfsr['Lindex'] = dfsr['rankclass']
    
    return dfxy, dfsr
    

       


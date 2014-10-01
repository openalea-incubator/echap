import pandas
import numpy

from openalea.deploy.shared_data import shared_data
import alinea.echap
from alinea.adel.plantgen.plantgen_interface import read_plantgen_inputs

#----------------------------------------------------- Fitted dimension

def Mercia_2010_fitted_dimensions():
    dim = {}
    for nff in [11,12,13]:
        fn = shared_data(alinea.echap, 'Mercia_dimT%d_user.csv'%(nff))
        dim[nff] = pandas.read_csv(fn)
    return dim
    
def Rht3_2010_fitted_dimensions():
    dim = {}
    for nff in [11,12]:
        fn = shared_data(alinea.echap, 'Rht3_dimT%d_user.csv'%(nff))
        dim[nff] = pandas.read_csv(fn)
    return dim
    
def Tremie12_fitted_dimensions():
    #dim = shared_data(alinea.echap, 'Tremie1_dimT_user.csv')
    dim = {}
    for nff in [12,13]:
        fn = shared_data(alinea.echap, 'Tremie12_dimT%d_user.csv'%(nff))
        dim[nff] = pandas.read_csv(fn)
    return dim

# ---------------------------------------------------- Fitted HS = f(TT)

class HaunStage(object):
    """ Handle HaunStage = f (ThermalTime) fits
    """
    
    def __init__(self, a_cohort = 1. / 110., TT_col_0 = 0):
        self.a_cohort = a_cohort
        self.TT_col_0 = TT_col_0
        
    def __call__(self, TT):
        return (numpy.array(TT) - self.TT_col_0) * self.a_cohort

HS_converter = {'Mercia': HaunStage(0.009380186, 101.4740799),
                'Rht3': HaunStage(0.008323522,65.00120875),
                'Tremie12': HaunStage(0.010736,236.2589295),
                'Tremie13': HaunStage(0.009380186, 101.4740799) # not fitted, use Mercia Fit
                }

#----------------------------------------------------- Plantgen
def plantgen_as_dict(inputs, dynT, dimT):
    d={}
    d['dynT_user'], d['dimT_user'], d['plants_number'],d['plants_density'], d['decide_child_axis_probabilities'], d['MS_leaves_number_probabilities'], d['ears_density'], d['GL_number'], d['delais_TT_stop_del_axis'], d['TT_col_break'],d['inner_params'] =  read_plantgen_inputs(inputs, dynT, dimT)
    return d
    
def GL_number():
    GL = {'Mercia': pandas.DataFrame({'TT':[485.90,1159.15,1732.10,1818.05,1959.55], '11':[3.5,3.4,3.774,2.670,1.732], '12':[3.6,3.6,4.036,3.010,1.903], '13':[4.4,3.4,4.7,3.8,2.5], 'mediane':[3.6,3.4,4.036,3.01,1.903]}),
         'Rht3': pandas.DataFrame({'TT':[485.90,1159.15,1732.10,1818.05,1959.55], '11':[3.5,3.3,3.9,2.7,1.7], '12':[3.9,3.6,4.3,3,2.2], '13':[3.3,4,2.9,None,None], 'mediane':[3.5,3.6,3.9,2.85,1.95]}),
         'Tremie12': pandas.DataFrame({'TT':[928.9,1044.25,1268.7,1392.15,1483.95,1568.7,1659.1,1797.95,1905,2084.95], '12':[1.7,2.4,4.3,4.7,4.7,4.5,4.3,3.7,2.2,0.6], '13':[1.9,2.1,4.1,5.2,4.6,4.4,4.1,3.4,2.3,1.0], 'mediane':[1.8,2.25,4.2,4.95,4.65,4.45,4.2,3.55,2.25,0.8]}),
         'Tremie13': pandas.DataFrame({'TT':[565.8,915,1048.8,1284.4,1351.4,1414.7,1545.7,1639,1762.5,1883.7], '11':[4.2,2.6,3.5,3.6,3.6,3.2,2.6,2,1.1,0.2], '12':[4.4,2.9,3.1,4.3,4.2,3.9,3.1,2,1.3,0.4], 'mediane':[4.3,2.75,3.3,3.95,3.9,3.55,2.85,2,1.2,0.3]})}
    return GL
'''   
GL_Mercia = {'GL_number' : {1732.1: 4.03631578947369, 1818.05:3.01047368421053,1959.55:1.90263157894737, 2108.65:0.0},
            'TT_col_break' : 0.0}   
GL_Rht3 = {'GL_number' : {1732.1: 3.88352941176471, 1818.05:2.68005882352941,1959.55:1.69764705882353, 2108.65:0.0},
            'TT_col_break' : 0.0}   
GL_Tremie = {'GL_number' : {1483.95: 4.9025, 1568.7:4.73684210526316,
             1659.1:4.16814814814815, 1797.95:3.43777777777778,
             1905:2.35888888888889,2084.95:0.85578947368421,2150:0.00},
            'TT_col_break' : 0.0} '''  
  
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
#
# Plot data
#
from math import sqrt
def valSD(listeVal):
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
    - *** IMPORTANT ***Plant density data measured for LAI estimation in excell fiiles have considere 4 ranks instead of 5 (confirmed by benjamin) (ie density and LAI should be multiplied by 0.8)
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
   
#
# LAI data
#
def PAI_data():
    d = {'Mercia': pandas.DataFrame({'code_date':['2011-04-18', '2011-04-27', '2011-06-01'],'TT_date':[1137,1281,1796],'PAI_vert_average_photo':[3.27,4.03,4.05],'PAI_vert_med_photo':[3.24,4.03,4.05],'PAI_vert_plt_photo':[9,2,2],'PAI_vert_SD_photo':[0.19,0.13,0.32]}), #SD=ecartype
         'Rht3': pandas.DataFrame({'code_date':['2011-04-18', '2011-04-27', '2011-06-01'],'TT_date':[1137,1281,1796],'PAI_vert_average_photo':[3.4,3.96,3.45],'PAI_vert_med_photo':[3.32,3.94,3.51],'PAI_vert_plt_photo':[10,3,5],'PAI_vert_SD_photo':[0.35,0.07,0.19]}),
         'Tremie12': pandas.DataFrame({'code_date':['2012-03-09', '2012-04-11', '2012-05-09'],'TT_date':[923,1260,1551],'PAI_vert_average_photo':[1.01,4.29,4.14],'PAI_vert_med_photo':[1.02,4.25,4.10],'PAI_vert_plt_photo':[13,14,8],'PAI_vert_SD_photo':[0.17,0.63,0.24],'LAI_vert_placette_biomasse':[0.43,0.43,0.43],'LAI_vert_average_biomasse':[1.10,4.31,4.43],'LAI_vert_med_biomasse':[1.14,4.66,4.41],'LAI_vert_plt_biomasse':[3,3,3],'LAI_vert_SD_biomasse':[0.1,0.6,0.3],'LAI_tot_average_biomasse':[1.43,None,None]}),
         'Tremie13': pandas.DataFrame({'code_date':['2013-04-17', '2013-04-26', '2013-05-29', '2013-06-13', '2013-04-22','2013-05-03'],'TT_date':[918,1018,1377,1612,941,1060],'PAI_vert_average_photo':[1.7,2.75,2.73,1.88,None,None],'PAI_vert_med_photo':[1.69,2.84,2.69,1.84,None,None],'PAI_vert_plt_photo':[11,14,12,13,None,None],'PAI_vert_SD_photo':[0.16,0.41,0.14,0.13,None,None],'LAI_vert_placette_biomasse':[None,None,None,None,0.36,0.36],'LAI_vert_average_biomasse':[None,None,None,None,(4.77*0.8),(5.82*0.8)],'LAI_vert_plt_biomasse':[None,None,None,None,15,15],'LAI_vert_SD_biomasse':[None,None,None,None,0.44,0.65]})}
         #'Tremie12': pandas.DataFrame({'TT':[1260,1550],'PAI':[3.8,4.7]})} #Publi Corinne
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
    
#-------------------------------------------------------------------------------  
# Angles / formes a plat  

# Processing of raw database to get data reabale or adel

# TO DO : automatic selection of tol_med = f(ntraj), based on histogram approx of normal distribution


def leaf_trajectories(dfxy, dfsr, bins = [-10, 0.5, 1, 2, 3, 4, 10], ntraj = 10, tol_med = 0.1):
    """
    Return a dynamic leaf database compatible with adel and appropriate dynamic key selecting functions ({Lindex:{num_traj:{age_class:xy_dict}}})
    
    - dfxy and srdb are dataframe containing the data
    - bins is None if leaf are static or define age clas otherwise
    - ntraj and tol_med control  the number of trajectory to sample among leaves of a given Lindex and of given age that have mean_angle +/- tol*med
    
    """
                   
    import random
          
    dfxy['age'] = dfxy['HS'] - dfxy['rank'] + 1
    #cut intervalle de 0 a 1, etc.
    dfxy_cut = pandas.cut(dfxy.age, bins)
    dfxy['age_class'] = dfxy_cut
    
    # use mean_angle to filter / group leaves
    mean_angle = dfxy.groupby('inerv').apply(lambda x: numpy.mean(abs(numpy.arctan2(numpy.diff(x['y'].values),numpy.diff(x['x'].values)))))

    #filter leaves that are above/below med*tol
    def filter_leaves(x):
        angles = mean_angle[set(x['inerv'])]
        med = angles.median()
        valid_angles = angles[(angles >= (1 - tol_med) * med) & (angles <= (1 + tol_med) * med)]
        return x[x['inerv'].isin(set(valid_angles.index))]
    
    validxy = dfxy.groupby(('Lindex','age_class'), group_keys=False).apply(filter_leaves)
    grouped = validxy.groupby(('Lindex','age_class'))
    
    # build trajectories
    trajectories = {k:[] for k in set(validxy['Lindex'])}
    for i in range(ntraj):
        for k in set(validxy['Lindex']):
            trajectories[k].append({})
            for t in set(validxy['age_class']):
                x = grouped.get_group((k,t))
                trajectories[k][i][t] = x.ix[x['inerv'] == random.sample(set(x['inerv']),1),['x','y']].to_dict('list')
    
    srdb = {k:v.ix[:,['s','r']].to_dict('list') for k, v in dfsr.groupby('Lindex')}
    
    return trajectories, srdb, bins
    


def leaf_curvature_data(name='Mercia', bins = [-10, 0.5, 1, 2, 3, 4, 10], ntraj = 10, tol_med = 0.1):

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
        dfxy = xy_reader(data_file_xydb)
        #dfxy = dfxy[dfxy['variety'].isin(['Mercia','Rht3'])]
        #dfxy = dfxy.reset_index()
        # SR
        dfsr = sr_reader(data_file_srdb)
        
    # add Lindex
    dfxy['Lindex'] = 1
    dfxy['Lindex'][dfxy['ranktop'] <= 4] = 2
    dfsr['Lindex'] = dfsr['rankclass']
    
    xy, sr, bins = leaf_trajectories(dfxy, dfsr, bins=bins, ntraj=ntraj, tol_med = tol_med)
    
    return xy, sr, bins
    

#
# Elaborated data
#    

def _add_ghs(df, g):
    hs_conv = HS_converter[g]
    df['Var'] = g
    df['HS'] = hs_conv(df['TT'])
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
    hs_conv = HS_converter['Tremie12']
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
    hs_conv = HS_converter['Tremie13']
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
    tdb =  Tillering_data_Tremie13_2012_2013()
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
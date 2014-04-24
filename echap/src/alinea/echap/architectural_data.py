import pandas

from openalea.deploy.shared_data import shared_data
import alinea.echap
from alinea.adel.plantgen.plantgen import read_plantgen_inputs

#-------------------------------------------------------------------------------  
def Mercia_2010_plantgen():
    dynT = shared_data(alinea.echap, 'Mercia_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Mercia_plantgen_inputs_MIN.py')
    d={}
    d['dynT_user'], d['dimT_user'], d['plants_number'],d['plants_density'], d['decide_child_axis_probabilities'], d['MS_leaves_number_probabilities'], d['ears_density'], d['GL_number'], d['delais_TT_stop_del_axis'], d['TT_col_break'],d['inner_params'] =  read_plantgen_inputs(inputs, dynT, dimT)
    return d
def Rht3_2010_plantgen():
    dynT = shared_data(alinea.echap, 'Rht3_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Rht3_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Rht3_plantgen_inputs_MIN.py')
    d={}
    d['dynT_user'], d['dimT_user'], d['plants_number'],d['plants_density'], d['decide_child_axis_probabilities'], d['MS_leaves_number_probabilities'], d['ears_density'], d['GL_number'], d['delais_TT_stop_del_axis'], d['TT_col_break'],d['inner_params'] =  read_plantgen_inputs(inputs, dynT, dimT)
    return d
def Tremie_2011_plantgen():
    dynT = shared_data(alinea.echap, 'Tremie1_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Tremie1_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Tremie1_plantgen_inputs_MINnew.py')
    d={}
    d['dynT_user'], d['dimT_user'], d['plants_number'],d['plants_density'], d['decide_child_axis_probabilities'], d['MS_leaves_number_probabilities'], d['ears_density'], d['GL_number'], d['delais_TT_stop_del_axis'], d['TT_col_break'],d['inner_params'] =  read_plantgen_inputs(inputs, dynT, dimT)
    return d
    
def HS_data():
    fn = shared_data(alinea.echap, 'HS_data_Mercia_Rht3_2010_2011.csv')
    #fn = shared_data(alinea.echap, 'HS_data_Tremie1_2011_2012.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    grouped = data.groupby(['variety'],as_index=False)
    return grouped.aggregate('mean')
    
#-------------------------------------------------------------------------------   
def _maxna(x):
    m = x.max()
    if any(x.isnull()) and m == 0:
        m = None
    return m
    
def diff(edata, s):
    e1 = edata.head(0)
    e2 = s.head(0) 
    e = list(set(e1) - set(e2))
    i = 0
    while i<(len(e)-1):
        s [e[i]] = None
        i = i + 1
    return s
    
def emis(data, d, n):
    res = d.ix[:,'Nff':] / n.ix[:,'Nff':] 
    for c in d.ix[:,'Nff':].columns.values:
        d.ix[:,c] = res.ix[:,c]
    # update TPE and compute MB and TT
    d['TPE'] = d.ix[:,'TC':'T5'].sum(axis=1)
    edata = data[data['Date'] >= 6]
    edata = edata.reset_index()
    grouped = edata.groupby(['Var','Date'],as_index=False)
    return grouped
    
#-------------------------------------------------------------------------------  
def Tillering_data_Mercia_Rht3_2010_2011():
    """Tillering data for Boigneville 2010-2011
    Found data are (id_data refers to column in data sythesis):
    - estimates of plant density at emergence (id Archi12)
    - estimates of ear density at harvest (all axes, id Archi10)
    - estimates of the number of elongated internodes on main stems
    - estimates of tiller dynamics, with three kind of data :
        - presence/absence of living mainstem (column MS) and number of leaf emited
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive), and/or  total number of primary tillers (column TPE) and secondary tillers emited (dead or alive, column TSE). These data are used for fitting emission probabilities
        - estimation of total living tillers (column TT) at the end from the counting of tillers that are alive or that bears an ear 
    Others data :
    - Number of secondary tiller to date 3 (Archi33)
    Notes :
    - No data to date 5
    - Missing raw data for plant counts at emergence (Archi7, used Archi12 instead)"""
    
    # summarise tillering data to mean number of tiller emited per plant (columns Tc->T5), mean number of leaves, and total tillers and MS present per plant at the end (column TT and MB)
    # this table was contrusted using Archi11,Archi19, Archi33
    fn = shared_data(alinea.echap, 'Tillering_data_Mercia_Rht3_2010_2011.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    # compute emmission probas of primary axis/ emission per plant of sencondary axis 
    edata = data[data['Date'] < 6]
    edata = edata.reset_index()
    grouped = edata.groupby(['Var', 'N'],as_index=False)
    d = grouped.agg(_maxna)
    d = d.reset_index()
    #---
    de = d.fillna(0)
    grouped1 = de.groupby(['Var','Date'],as_index=False)
    #---
    grouped = d.groupby(['Var','Date'],as_index=False)
    s = grouped.agg('sum')
    n = grouped1.agg(lambda x: float(x.dropna().count())) 
    diff(edata, s)
    d = s
    emis(data, d, n)
    s = grouped.agg('sum')
    n = grouped.agg(lambda x: x.dropna().count())
    d['MB'] = 1.0 * s['MB'] / n['MB']
    d['TT'] = 1.0 * s['TT'] / n['TT']
    obs = {'plant_density_at_emergence' : {'Mercia': 203, 'Rht3': 211}, 'ear_density_at_harvest' : {'Mercia': 444, 'Rht3' : 384}, 
            'tillering' : d}
    return obs

def Tillering_data_Tremie1_2011_2012():
    """Tillering data for Boigneville 2011-2012
    Expected data are :
    - estimates of plant density at emergence (id Archi1)
    - estimates of ear density at harvest (all axes) (id Archi16)
    - estimates of the number of elongated internodes on main stems
    - estimates of tiller dynamics, with three kind of data : 
        - presence/absence of living mainstem (column MS) and number of leaf emited
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive), and/or  total number of primary tillers (column TPE) and secondary tillers emited (dead or alive, column TSE). These data are used for fitting emission probabilities
        - estimation of total living tillers (column TT) at the end from the counting of tillers that are alive or that bears an ear
    Others data :
    - Number of secondary tiller to date 1 
    Notes :
    - No data to date 2, 4 & 5"""
    
    # summarise tillering data to mean number of tiller emited per plant (columns Tc->T7), mean number of leaves, and total tillers and MS present per plant at the end (column TT and MB)
    # this table was contrusted using Archi2, Archi3, Archi14, Archi15
    fn = shared_data(alinea.echap, 'Tillering_data_Tremie1_2011_2012.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    # compute emmission probas of primary axis/ emission per plant of sencondary axis 
    edata = data[data['Date'] < 6]
    edata = edata.reset_index()
    grouped = edata.groupby(['Var', 'N'],as_index=False)
    d = grouped.agg(_maxna)
    d = d.reset_index()
    #---
    de = d.fillna(0)
    grouped1 = de.groupby(['Var','Date'],as_index=False)
    #---
    grouped = d.groupby(['Var','Date'],as_index=False)
    s = grouped.agg('sum')
    n = grouped1.agg(lambda x: float(x.dropna().count())) 
    diff(edata, s)
    d = s
    emis(data, d, n)
    s = grouped.agg('sum')
    n = grouped.agg(lambda x: x.dropna().count())
    d['MB'] = 1.0 * s['MB'] / n['MB']
    d['TT'] = 1.0 * s['TT'] / n['TT']
    obs = {'plant_density_at_emergence' : {'Tremie': 278.55}, 'ear_density_at_harvest' : {'Tremie': 491}, 
            'tillering' : d}
    return obs
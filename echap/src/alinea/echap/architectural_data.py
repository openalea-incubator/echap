import pandas

from openalea.deploy.shared_data import shared_data
import alinea.echap
from alinea.adel.plantgen.plantgen import read_plantgen_inputs

def Mercia_2010_plantgen():
    dynT = shared_data(alinea.echap, 'Mercia_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Mercia_plantgen_inputs_MIN.py')
    d={}
    d['dynT_user'], d['dimT_user'], d['plants_number'],d['plants_density'], d['decide_child_axis_probabilities'], d['MS_leaves_number_probabilities'], d['ears_density'], d['GL_number'], d['delais_TT_stop_del_axis'], d['TT_col_break'],d['inner_params'] =  read_plantgen_inputs(inputs, dynT, dimT)
    return d

        
def HS_data_2010():
    fn = shared_data(alinea.echap, 'HS_data_Mercia_Rht3_2010_2011.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    grouped = data.groupby(['variety'],as_index=False)
    return grouped.aggregate('mean')
    
def Tillering_data_2010():
    """ Tillering data for Boigneville 2010-2011
    
    Expected data are :
    - estimates of plant density at emergence
    - estimates of ear density at harvest (all axes)
    - estimates of the number of elongated internodes on main stems
    - estimates of tiller dynamics, with three kind of data : 
        - presence/absence of living mainstem (column MS) and number of leaf emited
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive), and/or  total number of primary tillers (column TPE) and secondary tillers emited (dead or alive, column TSE). These data are used for fitting emission probabilities
        - estimation of total living tillers (column TT) at the end from the counting of tillers that are alive or that bears an ear """
    
    # summarise tillering data to mean number of tiller emited per plant (columns Tc->T5), mean number of leaves, and total tillers and MS present per plant at the end (column TT and MB)
    fn = shared_data(alinea.echap, 'Tillering_data_Mercia_Rht3_2010_2011.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    # compute emmission probas of primary axis/ emission per plant of sencondary axis 
    edata = data[data['Date'] < 6]
    edata = edata.reset_index()
    grouped = edata.groupby(['Var', 'N'],as_index=False)
    def _maxna(x):
        m = x.max()
        if any(x.isnull()) and m == 0:
            m = None
        return m
    d = grouped.agg(_maxna)
    d = d.reset_index()
    grouped = d.groupby(['Var','Date'],as_index=False)
    s = grouped.agg('sum')
    n = grouped.agg(lambda x: x.dropna().count())
    d = s
    d.ix[:,'Nff':] = 1.0 * d.ix[:,'Nff':] / n.ix[:,'Nff':]
    # update TPE and compute MB and TT
    d['TPE'] = d.ix[:,'TC':'T5'].sum(axis=1)
    edata = data[data['Date'] >= 6]
    edata = edata.reset_index()
    grouped = edata.groupby(['Var','Date'],as_index=False)
    s = grouped.agg('sum')
    n = grouped.agg(lambda x: x.dropna().count())
    d['MB'] = 1.0 * s['MB'] / n['MB']
    d['TT'] = 1.0 * s['TT'] / n['TT']
    obs = {'plant_density_at_emergence' : {'Mercia':203, 'Rht3':211},
            'ear_density_at_harvest' : {'Mercia': 444, 'Rht3' : 384}, 
            'tillering' : d}
    return obs
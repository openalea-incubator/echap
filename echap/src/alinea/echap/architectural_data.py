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
    fn = shared_data(alinea.echap, 'Tillering_data_Mercia_Rht3_2010_2011.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    grouped = data.groupby(['Date','Var'],as_index=False)
    d = grouped.agg('sum')
    n = grouped.agg(lambda x: x.dropna().count())
    d.ix[:,'Nff':] = 1.0 * d.ix[:,'Nff':] / n.ix[:,'Nff':]
    return d
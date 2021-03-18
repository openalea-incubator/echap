from openalea.deploy.shared_data import shared_data

from alinea.adel.plantgen.plantgen import *

def Mercia_MIN():
    import alinea.echap
    dynT = shared_data(alinea.echap, 'Mercia_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Mercia_plantgen_inputs_MIN.py')
    d={}
    d['dynT_user'], d['dimT_user'], d['plants_number'],d['plants_density'], d['decide_child_axis_probabilities'], d['MS_leaves_number_probabilities'], d['ears_density'], d['GL_number'], d['delais_TT_stop_del_axis'], d['TT_col_break'],d['inner_params'] =  read_plantgen_inputs(inputs, dynT, dimT)
    return d
    
# playing with plantgen
#d = Mercia_MIN()
#locals().update(d)
#from alinea.adel.plantgen import axeT, dimT, dynT, phenT, tools, params 
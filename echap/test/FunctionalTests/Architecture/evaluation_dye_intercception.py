""" Reconstruct a MTG wheat canopy at different treatment dates, save dye interception 
    and compare it with field data """

import pandas
from dye_interception import *
from openalea.deploy.shared_data import shared_data

def get_file_name(variety = 'Tremie12', nplants = 30,):
    return 'interception_'+variety.lower() + '_' + str(nplants) + 'pl.csv'

def get_file_path(variety = 'Tremie12', nplants = 30):
    filename = get_file_name(variety = variety, nplants = nplants)
    return str(shared_data(alinea.echap)/'interception_simulations'/filename)    

def save_dye_interception_single_variety(variety = 'Tremie12', nplants = 30,
                                         treatments = {'T1':1e4, 'T2':1e4},
                                         density=1, dimension=1):
    dfs = []
    for date, dose in treatments.iteritems():
        adel = None
        while adel is None:
            try:
                adel, g, df = repartition_at_application(name = variety, appdate = date, 
                                                         dose = dose, nplants = nplants, 
                                                         density = density, dimension = dimension)
            except:
                pass
        df = df.rename(columns = {'TT':'degree_days', 'plant':'num_plant', 
                                  'nff':'fnl', 'axe':'axis', 'ntop':'num_leaf_top'})
        df['date'] = df['date'].apply(lambda x: pandas.to_datetime(x, dayfirst=True))
        df['num_leaf_bottom'] = df['fnl'] - df['num_leaf_top'] + 1
    df = pandas.concat(dfs)
    df = df.reset_index(drop = True)
    df.to_csv(get_file_path(variety = variety, nplants = nplants), index = False)
    
def save_dye_interception(nplants = 30,treatments = {'T1':1e4, 'T2':1e4}, density=1, dimension=1):
    for variety in ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']:
        save_dye_interception_single_variety(variety = variety, nplants = nplants,
                                             treatments = treatments, density = density,
                                             dimension = dimension)
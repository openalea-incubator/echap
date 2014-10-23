''' Data for pesticide/dye intercption experiments
'''

import pandas
import alinea.echap
from openalea.deploy.shared_data import shared_data

def dye_interception():
    """ Dye interception data
    
    data are compiled mean/sd from dye interception.
    
    units : deposit per leaf (g per g.cm-2 applied)
    
    All volume mixed
    
    Raw data are in ECHAP_ARVALIS/Interception/Donnees brutes 3 annee
    
    original units : microg. per g.ha-1 (converted in g : g.cm-2 by multiplying by 100)
    """
    data_file = shared_data(alinea.echap, 'dye_interception.csv')
    df = pandas.read_csv(data_file, decimal=',', sep=';')
    return df
    
def dye_applications():
    """
    Date and measured Haunstage at date of application
    """
    apps = {'Mercia': {'T1': ('2011-04-19', 9.74),# mesure HS le 18/04, a corriger avec delta HS ?
                      'T2': ('2011-05-11', 12.8)},#HS not measured but estimated from HSconvert
            'Rht3': {'T1':('2011-04-19', 9.15), #mesure HS le 18
                     'T2':('2011-05-11', 12.48)},# Hsconvert estimation
            'Tremie12': {'T1':('2012-04-11', 10.98),
                     'T2':('2012-05-09', 12.63)},
            'Tremie13': {'T1':('2013-04-25', 8.7),#HSconvert estimation
                     'T2':('2013-05-17', 11.04) } #HS convert   etimation                
    }
    return apps
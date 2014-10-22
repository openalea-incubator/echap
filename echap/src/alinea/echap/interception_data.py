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
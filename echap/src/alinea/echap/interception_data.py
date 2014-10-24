''' Data for pesticide/dye interception experiments
'''

import pandas
import alinea.echap
from openalea.deploy.shared_data import shared_data

def dye_interception():
    """ Dye interception data
    
    data are  mean/sd deposit per leaf compiled from all sprayed volume experiment.
    To get comparable value, deposits are expressed relative to quantity sprayed
    
    units :  g per g.cm-2 applied
        
    Raw data are in ECHAP_ARVALIS/Interception/Donnees brutes 3 annee
    original units : microg per g.ha-1 (converted in g per g.cm-2 by multiplying original data by 100)
    
    Note:
        - at first date of application, for all cultivar, we hypothetise that the 1st leaf measured (F1) was in fact F2, as simulations indicate that F1 was small (5-10 cm2), hence probably not measured (confimed by Arvalis)
        - sd seems to be in original unit
    """
    data_file = shared_data(alinea.echap, 'dye_interception.csv')
    df = pandas.read_csv(data_file, decimal=',', sep=';')
    # shift leaf numbers at date 1
    decfeu = {'Mercia':'19/04/2011', 'Rht3':'19/04/2011', 'Tremie12':'11/04/2012', 'Tremie13':'25/04/2013'}
    for name,d in decfeu.iteritems():
        sel = (df['name'] == name) & (df['date'] == d) & (map(lambda x: x.startswith('F'),df['feuille']))
        df['feuille'][sel] = ('F2','F3','F4','F5')
        df['N feuille'][sel] = (2,3,4,5)
    return df
    
def dye_applications():
    """
    Date and measured Haunstage at date of application
    """
    apps = {'Mercia': {'T1': ('2011-04-19', 9.74),# mesure HS le 18/04, a corriger avec delta HS ?
                      'T2': ('2011-05-11', 12.8)},#HS not measured but estimated from HSconvert
            'Rht3': {'T1':('2011-04-19', 9.15), #mesure HS le 18
                     'T2':('2011-05-11', 12.48)},# Hsconvert estimation
            'Tremie12': {'date1':('2012-03-09', 7.55),#ajout date scan1
                     'date2':('2012-04-02', 10.15),#ajout date scan2
                     'T1':('2012-04-11', 10.98),
                     'T2':('2012-05-09', 12.63)},
            'Tremie13': {'date1':('2013-04-22', 8.36),#HSconvert estimation
                     'T1':('2013-04-25', 8.7),#HSconvert estimation
                     'date2':('2013-05-03', 9.7),#HS mesure le 02
                     'T2':('2013-05-17', 11.04) } #HS convert estimation                
    }
    return apps
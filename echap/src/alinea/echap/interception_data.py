''' Data for pesticide/dye interception experiments
'''

import pandas
import alinea.echap
from openalea.deploy.shared_data import shared_data
from alinea.echap.architectural_reconstructions import fit_HS
from alinea.echap.architectural_data import TT_lin

def dye_interception(todec = []):
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
    decfeu = {k:decfeu[k] for k in todec}
    for name,d in decfeu.iteritems():
        sel = (df['name'] == name) & (df['date'] == d) & (map(lambda x: x.startswith('F'),df['feuille']))
        df['feuille'][sel] = ('F2','F3','F4','F5')
        df['N feuille'][sel] = (2,3,4,5)
        # complete with missing leaf to allow barplot
        missing = df[sel].iloc[0]
        missing['feuille'] = 'F1'
        missing['N feuille'] = 1
        missing['mean'] = 0
        missing['sd'] = 0
        missing['IC'] = 0
        df = df.append(missing)
    
    df.rename(columns={'name':'var'}, inplace=True) 
    df.rename(columns={'N feuille':'ntop_cur'}, inplace=True)  
    df = df.drop('HS', 1) #HS not measured
    return df


#### To chek : retrieve true mesurement + date application + model hs on same graph

def TT_index():
    TT = TT_lin()
    return  TT.set_index(['label','Date'])
    
    
    
def tag_dates():
    tags = {'Mercia': {'T1': '2011-04-19', # date dye application
                       'T2': '2011-05-11'},
            'Rht3': {'T1':'2011-04-19',
                     'T2':'2011-05-11'},
            'Tremie12': {'scan_1' : '2012-03-09',#date scan (aka date1)
                         'scan_2' : '2012-04-02',
                         'scan_T1' : '2012-04-11',
                         'T1' : '2012-04-11',
                         'scan_T2': '2012-05-09',
                         'T2' : '2012-05-09'},
            'Tremie13': {'scan_T1' : '2013-04-22',
                         'T1' : '2013-04-25',
                         'scan_T2' : '2013-05-03',
                         'T2' : '2013-05-17'}
            }
    return tags                 
    
def dye_applications():
    """
    Date and measured Haunstage at date of application
    """
    apps = {
            'Mercia': {'T1-0.4': ('2011-04-19', 9.34),
                      'T1-0.2': ('2011-04-19', 9.54),
                      'T1': ('2011-04-19', 9.74),#HS mesure le 18/04, a corriger avec delta HS ?
                      'T1+0.2': ('2011-04-19', 9.94),
                      'T1+0.4': ('2011-04-19', 10.14),
                      'T2-2.5': ('2011-05-11', 10.3),
                      'T2-2': ('2011-05-11', 10.8),
                      'T2-1.5': ('2011-05-11', 11.3),
                      'T2-1': ('2011-05-11', 11.8),
                      'T2-0.5': ('2011-05-11', 12.3),
                      'T2-0.4': ('2011-05-11', 12.4),
                      'T2-0.2': ('2011-05-11', 12.6),
                      'T2': ('2011-05-11', 12.8),#HS not measured but estimated from HSconvert
                      'T2+0.2': ('2011-05-11', 13),
                      'T2+0.4': ('2011-05-11', 13.2),
                      'T2+0.5': ('2011-05-11', 13.3),
                      'T2+1': ('2011-05-11', 13.8),
                      'T2+1.5': ('2011-05-11', 14.3),
                      'T2+2': ('2011-05-11', 14.8),
                      'T2+2.5': ('2011-05-11', 15.3)},
            'Rht3': {'T1-0.4':('2011-04-19', 8.75),
                     'T1-0.2':('2011-04-19', 8.95),
                     'T1':('2011-04-19', 9.15), #HS mesure le 18
                     'T1+0.2':('2011-04-19', 9.35),
                     'T1+0.4':('2011-04-19', 9.55),
                     'T2-2.5': ('2011-05-11', 9.98),
                     'T2-2': ('2011-05-11', 10.48),
                     'T2-1.5': ('2011-05-11', 10.98),
                     'T2-1': ('2011-05-11', 11.48),
                     'T2-0.5': ('2011-05-11', 11.98),
                     'T2-0.4':('2011-05-11', 12.08),
                     'T2-0.2':('2011-05-11', 12.28),
                     'T2':('2011-05-11', 12.48),# HS convert estimation
                     'T2+0.2':('2011-05-11', 12.68),
                     'T2+0.4':('2011-05-11', 12.88),
                     'T2+0.5':('2011-05-11', 12.98),
                     'T2+1':('2011-05-11', 13.48),
                     'T2+1.5':('2011-05-11', 13.98),
                     'T2+2':('2011-05-11', 14.48),
                     'T2+2.5':('2011-05-11', 14.98)},
            'Tremie12': {'date1-0.4':('2012-03-09', 7.15),
                     'date1-0.2':('2012-03-09', 7.35),
                     'date1':('2012-03-09', 7.55),#ajout date scan1
                     'date1+0.2':('2012-03-09', 7.75),
                     'date1+0.4':('2012-03-09', 7.95),
                     'date2-0.4':('2012-04-02', 9.75),
                     'date2-0.2':('2012-04-02', 9.95),
                     'date2':('2012-04-02', 10.15),#ajout date scan2
                     'date2+0.2':('2012-04-02', 10.35),
                     'date2+0.4':('2012-04-02', 10.55),
                     'T1-0.4':('2012-04-11', 10.58),
                     'T1-0.2':('2012-04-11', 10.78),
                     'T1':('2012-04-11', 10.98),
                     'T1+0.2':('2012-04-11', 11.18),
                     'T1+0.4':('2012-04-11', 11.38),
                     'T2-2.5':('2012-05-09', 10.13),
                     'T2-2':('2012-05-09', 10.63),
                     'T2-1.5':('2012-05-09', 11.13),
                     'T2-1':('2012-05-09', 11.63),
                     'T2-0.5':('2012-05-09', 12.13),
                     'T2-0.4':('2012-05-09', 12.23),
                     'T2-0.2':('2012-05-09', 12.43),
                     'T2':('2012-05-09', 12.63),
                     'T2+0.2':('2012-05-09', 12.83),
                     'T2+0.4':('2012-05-09', 13.03),
                     'T2+0.5':('2012-05-09', 13.13),
                     'T2+1':('2012-05-09', 13.63),
                     'T2+1.5':('2012-05-09', 14.13),
                     'T2+2':('2012-05-09', 14.63),
                     'T2+2.5':('2012-05-09', 15.13)},
            'Tremie13': {'date1-0.4':('2013-04-22', 7.96),
                     'date1-0.2':('2013-04-22', 8.16),
                     'date1':('2013-04-22', 8.36),#HS convert estimation
                     'date1+0.2':('2013-04-22', 8.56),
                     'date1+0.4':('2013-04-22', 8.76),
                     'T1-0.4':('2013-04-25', 8.3),
                     'T1-0.2':('2013-04-25', 8.5),
                     'T1':('2013-04-25', 8.7),#HS convert estimation
                     'T1+0.2':('2013-04-25', 8.9),
                     'T1+0.4':('2013-04-25', 9.1),
                     'date2-0.4':('2013-05-03', 9.3),
                     'date2-0.2':('2013-05-03', 9.5),
                     'date2':('2013-05-03', 9.7),#HS mesure le 02
                     'date2+0.2':('2013-05-03', 9.9),
                     'date2+0.4':('2013-05-03', 10.1),
                     'T2-2.5':('2013-05-17', 8.54),
                     'T2-2':('2013-05-17', 9.04),
                     'T2-1.5':('2013-05-17', 9.54),
                     'T2-1':('2013-05-17', 10.04),
                     'T2-0.5':('2013-05-17', 11.54),
                     'T2-0.4':('2013-05-17', 10.64),
                     'T2-0.2':('2013-05-17', 10.84),
                     'T2':('2013-05-17', 11.04),#HS convert estimation 
                     'T2+0.2':('2013-05-17', 11.24),
                     'T2+0.4':('2013-05-17', 11.44),
                     'T2+0.5':('2013-05-17', 11.54),
                     'T2+1':('2013-05-17', 12.04),
                     'T2+1.5':('2013-05-17', 12.54),
                     'T2+2':('2013-05-17', 13.04),
                     'T2+2.5':('2013-05-17', 13.54)}                
    }
    return apps
    
def Petri_data(name='Tremie12'):
    HSconv = fit_HS()
    conv = HSconv[name]
    if name is 'Tremie12':
        data_file_T1 = shared_data(alinea.echap, 'petri_T1_20112012.csv')
        data_file_T2 = shared_data(alinea.echap, 'petri_T2_20112012.csv')
    if name is 'Tremie13':
        data_file_T1 = shared_data(alinea.echap, 'petri_T1_20122013.csv')
        data_file_T2 = shared_data(alinea.echap, 'petri_T2_20122013.csv')
    # HS T1 et T2
    dateT1, HS_T1 = dye_applications()[name]['T1']
    dateT2, HS_T2 = dye_applications()[name]['T2']
    # lecture des csv
    header_row=['PETRI','Volume','Niveau','Bloc','ABSORBANCE','DILUTION','concentration(mg/l)','ConcentrationArrondie(mg/l)','quantiteRetenue(mg)','quantite(g/ha)','rapportPoucentage(sol/emis)']
    df1 = pandas.read_csv(data_file_T1, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=0, decimal=',')
    df2 = pandas.read_csv(data_file_T2, dayfirst=True, names=header_row, sep=';', index_col=0, skiprows=0, decimal=',')
    # bon format pourcentage
    df1['rapportPoucentage(sol/emis)'] = df1['rapportPoucentage(sol/emis)'] / 100.
    df2['rapportPoucentage(sol/emis)'] = df2['rapportPoucentage(sol/emis)'] / 100.
    # ajout TT et conversion en HS
    df1['HS'] = HS_T1; df2['HS'] = HS_T2
    df1['TT'] = conv.TT(HS_T1); df2['TT'] = conv.TT(HS_T2)   
    # on filtre sur niveau = sol
    grouped1 = df1.groupby('Niveau'); grouped2 = df2.groupby('Niveau')
    petri1 = grouped1.get_group('sol'); petri2 = grouped2.get_group('sol')
    return petri1, petri2
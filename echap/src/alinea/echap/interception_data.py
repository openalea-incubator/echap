''' Data for pesticide/dye interception experiments
'''

import pandas
import numpy
import scipy.stats

import alinea.echap
from openalea.deploy.shared_data import shared_data
from alinea.echap.hs_tt import tt_hs_tag, derived_data_path
#from alinea.echap.architectural_reconstructions import fit_HS
#from alinea.echap.architectural_data import TT_lin, Pheno_data

        
def conf_int(lst, perc_conf=95):
    """
    Confidence interval - given a list of values compute the square root of
    the variance of the list (v) divided by the number of entries (n)
    multiplied by a constant factor of (c). This means that I can
    be confident of a result +/- this amount from the mean.
    The constant factor can be looked up from a table, for 95pcent confidence
    on a reasonable size sample (>=500) 1.96 is used.
    """
    
    n, v = len(lst), numpy.var(lst, ddof=1)
    c = scipy.stats.t.interval(perc_conf * 1.0 / 100, n-1)[1]
    
    return numpy.sqrt(v/n) * c

def TT_index():
    TT = TT_lin()
    return  TT.set_index(['label','Date'])
    
    
    
def tag_dates():
    tags = {'Mercia': {'T1-90':'2011-01-19',
                       'T1-1': '2011-04-18',
                       'T1': '2011-04-19',
                       'T2-14':'2011-04-27',
                       'T2': '2011-05-11',
                       'T2+15':'2011-05-26',
                       'T2+21':'2011-06-13'},
            'Rht3': {'T1-1': '2001-04-18',
                     'T1': '2011-04-19',
                     'T2-14':'2011-04-27',
                     'T2': '2011-05-11',
                     'T2+15':'2011-05-26'},
            'Tremie12': {'T1-33' : '2012-03-09',#date scan (aka date1)
                         'T1-9' : '2012-04-02',
                         'T1' : '2012-04-11',
                         'T2' : '2012-05-09',
                         'T2+34':'2012-06-12'},
            'Tremie13': {'T1-8':'2013-04-17',
                         'T1-6':'2013-04-19',
                         'T1-3' : '2013-04-22',
                         'T1' : '2013-04-25',
                         'T1+1' : '2013-04-26',
                         'T1+4': '2013-04-29',
                         'T2-15':'2013-05-02',
                         'T2-14' : '2013-05-03',
                         'T2' : '2013-05-17',
                         'T2+5' : '2013-05-22',
                         'T2+12' : '2013-05-29',
                         'T2+27' : '2013-06-13'}
            }
    return tags

def tag_treatments():
    tags= {'Mercia': {'application': ['T1','T2'],
                      'hs': ['T1-1','T2-14', 'T2+15']},
           'Rht3': {'application': ['T1','T2'],
                    'hs': ['T1-1','T2-14', 'T2+15']},
           'Tremie12': {'application': ['T1','T2'],
                        'hs': ['T1','T2'],
                        'scan': ['T1-33','T1-9','T1','T2'],
                        'silhouette': ['T1-33', 'T1', 'T2', 'T2+34']},
            'Tremie13': {'application': ['T1','T2'],
                        'hs': ['T1-3','T2-14','T2+5'],
                        'scan': ['T1-3','T2-14'],
                        'silhouette': ['T1+4']}}
    return tags    
    
def dye_interception(coef=0.045):
    """ Dye interception data
    
    coef : slope of the absobance = f(concentration) [mg.l-1] relationships 
           (absorbance = coef * concentration, ie conc = abs/coeff)
           In original files, in 2011 and 2012 a coef of 0.0448 was used whereas a coeff of 0.0453 was used in 2013
           here we use a mean value of 0.045
    
    units:
        volume : applied volumetric rate (l.ha-1)
        concentration: applied concentration (g.l-1)
        dilution: water volume used to dilute intercepted dye (cl)
        absorbance :optical read
        
        thus:
        quantity intercepted (mg) = abs / coef * dilution / 100
 
        deposit (g / g.cm-2 applied) = abs / coef * dilution / volume / concentration * 1000
       
        
    Raw data are in ECHAP_ARVALIS/Interception/Donnees brutes 3 annee
    original units : microg per g.ha-1 (converted in g per g.cm-2 by multiplying original data by 100)
    
    Note:
        - at first date of application, for all cultivar, we hypothetise that the 1st leaf measured (F1) was in fact F2, 
         as simulations indicate that F1 was small (5-10 cm2), hence probably not measured (confimed by Arvalis)
    """
    data_file = shared_data(alinea.echap, 'dye_interception.txt')
    df = pandas.read_csv(data_file, decimal=',', delim_whitespace=True)
    df['deposit'] = df['absorbance'] / coef * df['dilution'] / df['volume'] / df['concentration'] * 1000 
    #add aggregators
    df['axe'] = 'MS'
    df['ntop_cur'] = map(lambda x: int(x.split('F')[1]), df['leaf'])
    # hypothetic ntop_lig: one leaf is growing and has beeen measured, except if all leaves are ligulated
    df['ntop_lig'] = df['ntop_cur'] - 1
    # adjust
    #
    #  T1 : derniere feuille trop petite, pas mesuree a priori
    skiplast = {'T1':['Mercia','Rht3','Tremie12','Tremie13']}
    for t in skiplast:
        df.loc[(df['treatment'] == t) & (df['variety'].isin(skiplast[t])),'ntop_cur'] += 1
    #
    all_ligulated = {'T2':['Mercia','Rht3','Tremie12']}
    for t in all_ligulated:
        df.loc[(df['treatment'] == t) & (df['variety'].isin(all_ligulated[t])),'ntop_lig'] += 1
    #two_growing = {'T1':['Tremie13']}
    #for t in two_growing:
    #    df.loc[(df['treatment'] == t) & (df['variety'].isin(two_growing[t])),'ntop_lig'] -= 1    
    #index variety                                                           
    df=df.set_index('variety')
                                                                      
    return df


def petri_dye_interception(coef=0.045, diameter=8.5):
    """ Dye interception data in petri dishes placed above and below the canopy
    
    coef : slope of the absobance = f(concentration) [mg.l-1] relationships 
           (absorbance = coef * concentration, ie conc = abs/coeff)
           In original files, in 2011 and 2012 a coef of 0.0448 was used whereas
            a coeff of 0.0453 was used in 2013
           here we use a mean value of 0.045
    diameter: diameter of the petri dish (cm)
    
    units in the data columns:
        volume : applied volumetric rate of dye (l.ha-1)
        concentration: applied concentration (g.l-1)
        dilution: water volume used to dilute intercepted dye (cl)
        absorbance :optical read
        
    deposit (g / g.cm-2 applied) =
                        abs / coef * dilution / volume / concentration * 1000
    """

    data_file = shared_data(
        alinea.echap) / 'interception_data' / 'dye_interception_petri.csv'
    df = pandas.read_csv(data_file, decimal=',', sep=';')
    df['variety'] = numpy.where(df['year'] == 2012, 'Tremie12', 'Tremie13')
    df['deposit'] = df['absorbance'] / coef * df['dilution'] / df['volume'] / \
                    df['concentration'] * 1000
    df['po'] = df['deposit'] / (numpy.pi * diameter ** 2 / 4)
    treatment = df.pop('treatment')
    var = df['variety']
    df['daydate'] = numpy.where(var == 'Tremie12',
                                numpy.where(treatment == 'T1', '2012-04-11',
                                            '2012-05-09'),
                                numpy.where(treatment == 'T1', '2013-04-25',
                                            '2013-05-17'))
    return df


def petri_data(variety='Tremie12', tag='reference', level='soil'):
    df = None
    path = derived_data_path(None) / 'petri_dye_interception.csv'
    try:
        df = pandas.read_csv(path)
    except IOError:
        df = petri_dye_interception().loc[:,
             ('variety', 'daydate', 'level', 'po')]
        df.to_csv(path, index=False)

    if variety in df['variety'].values:
        df = df.loc[(df['variety'] == variety) & (df['level'] == level), :]
        tths = tt_hs_tag(variety, tag)
        df = df.merge(tths)
    return df


def gap_fraction():
    """ Gap fraction (non green fraction) estimated from vertical images
    """
    data = [pandas.read_csv(shared_data(alinea.echap,var + '_vertical_images.csv'),decimal=',', sep=';')
            for var in ['MerciaRht3','Tremie12','Tremie13']]
    df = pandas.concat(data)
    df['gap_fraction'] = (100 - df['pcent_veg']) / 100.
    tags = {'Mercia':{'19/01/2011':'T1-90','18/04/2011':'T1-1','27/04/2011':'T2-14','01/06/2011':'T2+21'},
            'Rht3':{'19/01/2011':'T1-90','18/04/2011':'T1-1','27/04/2011':'T2-14','01/06/2011':'T2+21'},
            'Tremie12':{'09/03/2012':'T1-33', '11/04/2012':'T1','09/05/2012':'T2'},
            'Tremie13':{'17/04/2013':'T1-8','26/04/2013':'T1+1', '29/05/2013':'T2+12','13/06/2013':'T2+27'}}
    df['treatment'] = map(lambda (var,d): tags[var][d], zip(df['variety'], df['date']))
    
    return df.groupby(['variety','treatment']).agg('mean').reset_index()
    
def scan_data():
    data_file = shared_data(alinea.echap, 'architectural_measurements/Compil_scan.csv')
    df = pandas.read_csv(data_file, decimal='.', sep=',')
    tags = {'Tremie12':{'09/03/2012':'T1-33', '02/04/2012':'T1-9', '11/04/2012':'T1','09/05/2012':'T2'},
            'Tremie13':{'22/04/2013':'T1-3', '03/05/2013':'T2-14'}}
    df['treatment'] = map(lambda (var,d): tags[var][d], zip(df['variety'], df['prelevement']))
    df = df[df['id_Axe'] == "MB"]
    #add aggregator
    df['axe'] = 'MS'
    df['metamer'] = df['rank']
    df['ntop_lig'] = df['Nflig'] - df['rank'] + 1
    #rem : for some date, nmax is not highest leaf, hence ntop_cur may be wrong
    def _ntopcur(x):
        nmax = x['metamer'].max()
        x['ntop_cur'] = nmax - x['metamer'] + 1
        return x
    df = df.groupby(['variety','treatment', 'N']).apply(_ntopcur).reset_index()
    #index variety                                                           
    df = df.set_index('variety')
    return df
 

def silhouettes():
    data_file_xydb = shared_data(alinea.echap, 'architectural_measurements/xydb_Boigneville_Tremie12_Tremie13.csv')
    df = pandas.read_csv(data_file_xydb)
    # Wrong data for plants 19, 20, 21 on harvest 2
    df = df[~((df['harvest']==2) & (df['plant'].isin([19, 20, 21])))]
    tags = {'Tremie12':{'1':'T1-33', '2': 'T1', '3': 'T2', '4':'T2+34'},
            'Tremie13': {'5': 'T1+4'}}
    df['treatment'] = map(lambda (var, h): tags[var][str(h)], zip(df['variety'], df['harvest']))
    # reduce to projection factor / hins
    def _compress(x):
        return pandas.Series({'HS':x['HS'].values[0], 'insertion_height':x['hins'].values[0], 
                              'h_projection': (x['x'].max() - x['x'].min()) / numpy.sum(numpy.sqrt(numpy.diff(x['x'])**2 + numpy.diff(x['y'])**2))})
                              
    data = df.groupby(['variety','treatment','plant','Nflig','rank','ranktop','relative_ranktop']).apply(_compress).reset_index()
    #add aggregator
    data['axe'] = 'MS'
    data = data.rename(columns={'rank':'metamer', 'ranktop': 'ntop', 'relative_ranktop': 'ntop_cur'})
    data['ntop'] = data['ntop'].round()
    data['ntop_cur'] = data['ntop_cur'].round()
    data['ntop_lig'] = data['Nflig'] - data['metamer'] + 1
    #index variety                                                           
    data = data.set_index('variety')
    return data
    
def haun_stages():
    tags = {'Mercia': {'2011-04-18':'T1-1', '2011-04-27':'T2-14','2011-05-26':'T2+15'},
            'Rht3': {'2011-04-18':'T1-1', '2011-04-27':'T2-14','2011-05-26':'T2+15'},
            'Tremie12':{'2012-03-09':'T1-33', '2012-04-02':'T1-9', '2012-04-11':'T1','2012-05-09':'T2','2012-06-12':'T2+34'},
            'Tremie13':{'2013-04-19':'T1-6','2013-04-22':'T1-3', '2013-04-29':'T1+4','2013-05-02':'T2-15', '2013-05-03':'T2-14','2013-05-22':'T2+5'}}
    
    def aggregate(data, column = 'HS', group = ['TT', 'label', 'nff'],
              func = numpy.mean, column_name = 'HS_mean'):
        df = data.groupby(group).agg(func).loc[:, column].to_frame()
        df = df.rename(columns={column : column_name})
        return df
    
    # Get and select data
    data = Pheno_data()
    tagged = data['archi_tagged']
    tagged = tagged.loc[:,('Date', 'label', 'TT', 'HS')].dropna()
    sampled = data['archi_sampled']
    sampled = sampled.loc[:,('Date', 'label','TT','HS')].dropna()
    data = pandas.concat([tagged, sampled])
      # Get mean, standard error and confidence interval for all NFFs for tagged plants and sampled plants
    df_HS = pandas.concat([aggregate(data, 'HS', ['Date', 'TT', 'label'], numpy.mean, 'HS_mean'),
                           aggregate(data, 'HS', ['Date', 'TT', 'label'], numpy.nanstd, 'HS_std'),
                           aggregate(data, 'HS', ['Date', 'TT', 'label'], conf_int, 'HS_conf')], axis=1)
                           
    df_HS = df_HS.reset_index()
    df_HS = df_HS[df_HS['Date'].isin(reduce(lambda x,y: x+y, [tags[v].keys() for v in tags],[]))]
    df_HS['treatment'] = map(lambda (var,d): tags[var][d.strftime('%Y-%m-%d')], zip(df_HS['label'], df_HS['Date']))
    df_HS = df_HS.reset_index().set_index('label')
    
    return df_HS
    
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
    

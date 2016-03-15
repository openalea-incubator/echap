import pandas
import numpy
try:
    import cPickle as pickle
except:
    import pickle
    
from openalea.deploy.shared_data import shared_data
import alinea.echap

def dateparse(x):
    return pandas.datetime.strptime(x, '%d/%m/%Y')
    
def TT_lin():
    filename = {'Mercia': 'MerciaRht3_TTlin_sowing.csv', 
               'Rht3': 'MerciaRht3_TTlin_sowing.csv', 
               'Tremie12': 'Tremie12_TTlin_sowing.csv', 
               'Tremie13': 'Tremie13_TTlin_sowing.csv'}
    def _read(filepath, label):
        df = pandas.read_csv(filepath, sep=';', decimal = ',')
        df['Date'] = df['Date'].apply(dateparse)
        df['label'] = label
        return df
    data = {k:_read(shared_data(alinea.echap)+ '/architectural_measurements/' + v, k) for k,v in filename.iteritems()}
    return pandas.concat(data.values())


#
# ____________________________________________________________________________Plot data
#

def Plot_data_Mercia_Rht3_2010_2011():
    """
    Plot data for Boigneville 2010-2011 
    
    data include :
    - estimates of plant density at emergence (id Archi12). Emergence means 75% plants emerged
    - estimates of ear density at harvest (all axes, id Archi10)
    
    Notes :
    - Missing raw data for plant counts at emergence (Archi7, used Archi12 instead)
    """
    d={'Mercia': {
        'code_date': {'sowing':'2010-10-15', 'emergence':'2010-11-02', 'harvest':'2011-06-20'},
        'TT_date': {'2010-10-15':0, '2010-11-02':143, '2011-06-20':2109},
        'sowing_density': 220, 
        'plant_density_at_emergence' : 203, 
        'raw_ear_density_at_harvest':[507,440,427,357,430,410,427,480,363, 350, 497, 410, 493, 340, 407, 467, 490, 433, 547, 450, 527, 427, 483, 493],
        'inter_row': 0.15},
        'Rht3': {
        'code_date': {'sowing':'2010-10-15', 'emergence':'2010-11-02', 'harvest':'2011-06-20'},
        'TT_date': {'2010-10-15':0, '2010-11-02':143, '2011-06-20':2109},
        'sowing_density': 220, 
        'plant_density_at_emergence' : 211, 
        'raw_ear_density_at_harvest':[440,420,330,433,347,410,387,367,360,397,357,377,427,367,380,347,330,420,423,397,383,367,377,377],
        'inter_row': 0.15}}
    for g in ('Mercia','Rht3'):
        d[g]['ear_density_at_harvest'] = numpy.mean(d[g]['raw_ear_density_at_harvest'])
        d[g]['ear_density_at_harvest_SD'] = numpy.nanstd(d[g]['raw_ear_density_at_harvest'])
    return d

def Plot_data_Tremie_2011_2012():
    """
    Plot data for Boigneville 2011-2012 
    
    Notes:
    - No data found for density of emergence as this trial was followed only after winter (winter trial was on appache)
    - plant density data at date 3 (11/04/2012) and 4 comes from axe counting/ LAI measurement data taken on  5 rank * 60 cm prelevement)
    """
    d = {
    'code_date':{'sowing': '2011-10-21','emergence': '2011-11-03', 'harvest':'2012-06-19', 'd3': '2012-04-11', 'd4': '2012-05-09', 'arvalis':'2012-03-20'}, 
    'TT_date': {'2011-10-21':0, '2011-11-03':141, '2012-06-19':2135, '2012-03-20':1006, '2012-04-11':1240, '2012-05-09':1515},
    'sowing_density': 280,
    'raw_ear_density_at_harvest':[479, 490, 598, 608, 538, 503, 493, 430, 458, 437, 486, 489, 465, 406],
    'inter_row':0.143,
    'plant_density':{'2012-03-20': [238, 304, 287, 237, 290, 301, 290, 287, 273],
                     '2012-04-11': [301, 273, 315], 
                     '2012-05-09':[257, 301, 263]},
    'axe_density': {'2012-04-11': [918, 956, 979],
                    '2012-05-09':[601, 585, 506]}, 
    'fertile_axis_density':{'2012-05-09':[545, 569, 443]} 
    }
    d['ear_density_at_harvest']=numpy.mean(d['raw_ear_density_at_harvest'])
    d['ear_density_at_harvest_SD'] = numpy.nanstd(d['raw_ear_density_at_harvest'])
    d['mean_plant_density'] = numpy.mean(reduce(lambda x,y:x+y,d['plant_density'].values()))
    return d
    
def Plot_data_Tremie_2012_2013():
    """
    Plot data for Boigneville 2012-2013 
    
    Notes
    - no date found for plant density counts at stage 2F and epis1cm (estimation epis 1cm : 9/04/2012)
    - no date found for countings of ear density
    - plant density and ear density were estimated on 2 ranks * 1m plots at stage 2F and epis1cm (0.3 m2), and on plots of 5 ranks * 0.6 m for LAI plots (0.45 m2) 
    - *** IMPORTANT ***Plant density data measured for LAI estimation in excel files have consider 4 ranks instead of 5 (confirmed by Benjamin) (ie density and LAI should be multiplied by 0.8)
    """
    d = {
   'code_date':{'sowing': '2012-10-29','emergence': '2012-11-19', '2F':'2013-01-03', 'epis1cm': '2013-04-09', 'LAI1':'2013-04-22', 'LAI2': '2013-05-13', 'harvest': '2013-06-09'},
   'TT_date': {'2012-10-29':0, '2012-11-19':140, '2013-04-22':941, '2013-05-13':1185, '2013-01-03':421, '2013-04-09':811, '2013-06-09': 1548},
   'sowing_density': 300,
   'plant_density':{'2013-01-03': [237, 287, 217, 237, 293, 220, 253, 213, 220, 253], #'2F'
                    '2013-04-09': [203, 203, 193, 207, 223, 203, 207, 207, 200, 197], #'epis1cm'
                    '2013-04-22':[328 * 0.8, 322 * 0.8, 356 * 0.8],
                    '2013-05-13':[272* 0.8, 361 * 0.8, 350 * 0.8]},
   'raw_ear_density_at_harvest':[643, 580, 693, 813, 663, 670, 693, 763, 567, 590, 707, 637, 617, 747, 760, 693, 653, 670],
    'inter_row':0.15
    }
    d['ear_density_at_harvest'] = numpy.mean(d['raw_ear_density_at_harvest'])
    d['ear_density_at_harvest_SD'] = numpy.nanstd(d['raw_ear_density_at_harvest'])
    d['mean_plant_density'] = numpy.mean(reduce(lambda x,y:x+y,d['plant_density'].values()))
    return d

def Plot_data():
    d = Plot_data_Mercia_Rht3_2010_2011()
    d.update({'Tremie12':Plot_data_Tremie_2011_2012(),
              'Tremie13':Plot_data_Tremie_2012_2013()})
    return d

#
# ____________________________________________________________________________Tillering data
#
# Utilities for processing tillering data
def _maxna(x):
    m = x.max()
    if any(x.isnull()) and m == 0:
        m = None
    return m

def emission_probabilities(df, last='T6'):
    grouped = df.groupby('N',as_index=False)
    em = grouped.agg(_maxna)
    # em = em.reset_index()
    s = em.ix[:,'TC':last].sum()
    n = em.ix[:,'TC':last].apply(lambda x: x.dropna().count())
    probas = s / n
    return probas.to_dict()
    
def plant_viability(df):
    grouped = df.groupby('Date', as_index=False)
    res = grouped.apply(lambda x: x['MB'].count() * 1.0 / len(x['MB']))
    return {'Date':res.index.tolist(), 'viability':res.tolist()}
  
def axis_dynamics(df):
    grouped = df.groupby('Date')
    s = grouped.agg('sum').ix[:,'TP':]
    n = grouped.agg(lambda x: x.apply(lambda x: x.dropna().count())).ix[:,'TP':]
    axis =  s / n
    axis = axis.replace(numpy.inf, numpy.nan)
    res = axis.to_dict('list')
    res.update({'Date':axis.index.tolist()})
    return res
    
def nff_probabilities(df):
    prop = 1. * numpy.bincount(map(int,df['Nff'].dropna())) / len(df['Nff'].dropna())
    d = dict(zip(range(len(prop)), prop))
    return {k:v for k,v in d.iteritems() if v > 0}
    

def Tillering_data_Mercia_Rht3_2010_2011():
    """Tillering data for Boigneville 2010-2011
    
    Data are from 36 tagged plant per cultivar, followed during the complete season.
    
    Data found in Archi11,Archi19, Archi33 were used to build a synthetic table containing :
        - presence/absence of living mainstem (column MS, NA means the plant is dead) and number of leaf emited (Nff)
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive)
        - total number of primary tillers (column TP) and/or secondary tillers (dead or alive, column TS), and/or total number of tiller present (primary or secondary, dead or alive), and/or total number of tiller with at least 2 ligulated leaves (column TT3F).
        - number of fertile tillers (column FT) at the end from the counting of tillers that are alive or that bears an ear 
        
    These data are aggregated to estimate primary emission probability, dynamics of mean number of axis present on plants, plant mortality and  number of ears per plant
        
    
    Notes :
    - No tillering data at date 4 and 5
    - when information is available at date 1 or 3, question marks at date 2 were replaced by confirmed tiller positions 
    - at date 3 , for the 2x 24 plants not measured in details, we invert column 'total primary' and 'total secondary' as it makes data much more consistent with date 2
    - at date 3 TT3F was estimated as total primary  + delta secondary
    - at date 3, on the 2 x 12 plants measured in details, it was possble to count presence/absence of living primary tillers (T1->T5). Living wes defined as Green leaves > 2 (minimal value of the GL model). 
"""
    
    fn = str(shared_data(alinea.echap)/'architectural_measurements'/'MerciaRht3_Tillering_data.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    date_code = {'d1':'2011-01-19', 'd2':'2011-04-18', 'd3':'2011-04-27', 'd4':'2011-05-26', 'd5':'2011-06-01', 'd6':'2011-06-09'}
    TT_code = {'d1':473, 'd2':1145, 'd3':1278, 'd4':1710, 'd5':1800, 'd6':1940}
    
    # infer emision at date 3 from Total primary present (TP)
    #add statistics for T6 as well as it may be present at date 3    
    def infer_d3(g):
        g.insert(13,'T6',0.0)
        g['T6'][g['MB']!=1] = numpy.nan # avoid infering T6 = 0 on dead plants
        TP = g['TP']
        date = g['Date']
        if all(TP[date==3].notnull()) and all(pandas.isnull(g.ix[date==3,'TC':'T6'])):
            infer = g.ix[date==2,'TC':'T6'] 
            if all(TP[date==3] > TP[date==2]):
                d = int(TP[date==3].values - TP[date==2].values)
                try:
                    lastT = int(max(2, max(numpy.where(infer > 0)[1])))
                except ValueError:
                    lastT = 2
                infer.ix[:,(lastT+1):(lastT + 1 + d)] = 1                
            g.ix[g['Date']==3,'TC':'T6'] = infer
        return g
        
    grouped = data.groupby(['Var', 'N'],as_index=False)
    newdata = grouped.apply(infer_d3)
    
    # compute emmission probability using notations before date 6
    edata = newdata[newdata['Date'] < 6]
    edata = edata.reset_index()
    grouped = edata.groupby('Var',as_index=False)
    emission = {k:emission_probabilities(v) for k,v in grouped}
    
    # compute nff probabilities
    grouped = data.groupby('Var',as_index=False)
    nff_prop = {k:nff_probabilities(v) for k,v in grouped}
    
    # compute ear_per_plante (including main stem) 
    eardata = newdata[newdata['Date'] >= 6]
    eardata = eardata.reset_index()
    grouped = eardata.groupby('Var',as_index=False)
    ears_per_plant = {k:  1  + (v['FT'].sum() / v['FT'].dropna().count()) for k,v in grouped}
    
    # compute plant viability
    grouped = data.groupby('Var',as_index=False)
    viability = {k:plant_viability(v) for k,v in grouped}
    
    #compute tillering dynamics
    axdyn = {k:axis_dynamics(v) for k,v in grouped}
    
    obs = {k:{'emission_probabilities': emission[k],
           'nff_probabilities': nff_prop[k],
           'ears_per_plant': ears_per_plant[k],
           'plant_survival': viability[k],
           'tillers_per_plant': axdyn[k], 
           'date_code' : date_code,
           'TT_code' : TT_code} for k in ('Mercia', 'Rht3')}
    return obs

def Tillering_data_Tremie12_2011_2012():
    """Tillering data for Boigneville 2011-2012
    
    Data come from sampled/scanned plants at date 1, from sampled + 'silhoueted' plants at date 6, from tagged plants at date 3 for presence / absence of T7 and at date 7 for counts of fertile tillers.
    Counting of axe per plant were also done at stage epis1cm on an independant set of plants, that was renamed here date 8.
    
    Data found in Archi2, Archi3, Archi14, Archi15 were :
        - presence/absence of living mainstem (column MS, NA means the plant is dead) and number of leaf emited (Nff)
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive)
        - total number of primary tillers (column TP) and/or secondary tillers (dead or alive, column TS), and/or total number of tiller present (primary or secondary, dead or alive), and/or total number of tiller with at least 2 ligulated leaves (column TT3F).  
        - number of fertile tillers (column FT) at the end from the counting of tillers that are alive or that bears an ear 
        
    These data are aggregated to estimate primary emission probability, dynamics of mean number of axis present on plants, plant mortality and  number of ears per plant
 
    Notes :
    - at Date 1, it was noted 'due to freezing of plants, it was difficult to determine main stem'
    - at date 2 there were damage due to to fly on Mb , T1, T2
    - At date 3, Mariem noted "T7 almost always present on scaned plants"
    - Emission of tiller 5 and 6 were never noted.
    - at date 7, values of fertile tiller  per plants seems buggy compared to other dates
    """
    
    fn = str(shared_data(alinea.echap)/'architectural_measurements'/'Tremie12_Tillering_data.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    date_code = {'d1':'2012-03-09', 'd2':'2012-04-02', 'd3':'2012-04-11', 'd4':'2012-05-09', 'd5':'2012-05-29', 'd6':'2012-06-12', 'd7':'2012-07-12', 'd8':'2012-04-04'}
    TT_code = {'d1':905, 'd2':1160, 'd3':1240, 'd4':1515, 'd5':1813, 'd6':2031, 'd7':2536, 'd8':1179}
    
    # compute emmission probability from data at date 1, and adding data at date 3 for tiller 7 only
    edata = data[data['Date'] == 1]
    edata.loc[:, 'T7'] = data['T7'][data['Date'] == 3].values
    edata = edata.reset_index()
    emission = emission_probabilities(edata,last='T7')
    
    # compute nff probabilities
    nff_prop = nff_probabilities(data)

    # tiller dynamics
    axdyn = axis_dynamics(data)
    
    # compute ear_per_plant using data at date 6 and 7 and plot data of fertile tillers at date 4
    eardata = pandas.DataFrame(axdyn).ix[:,('Date', 'FT')].dropna()
    pdata = Plot_data_Tremie_2011_2012()
    ftd4 = 1.*numpy.array(pdata['fertile_axis_density']['2012-05-09'])  / numpy.array(pdata['plant_density']['2012-05-09']) - 1 #remove main stem to get FT count 
    # ear per plante at date 7 very strange (twice the value at other dates and non-existence of unfertile tillers): ears_per_plant taken as mean of counting at date 4 and 6
    ears_per_plant = 1 + numpy.array(ftd4.tolist() + eardata['FT'][eardata['Date'] == 6].tolist()).mean()

    
    obs = {'emission_probabilities': emission,
           'nff_probabilities': nff_prop,
           'ears_per_plant': ears_per_plant,
           'tillers_per_plant': axdyn, 
           'date_code' : date_code,
           'TT_code' : TT_code}
    return obs
    
def Tillering_data_Tremie13_2012_2013():
    """Tillering data for Boigneville 2012-2013
    
    Data for presence/absence of primary tillers come from tagged plants. Data at date 3 are from other other plants.
    
    Data found in were :
        - presence/absence of living mainstem (column MS, NA means the plant is dead) and number of leaf emited (Nff)
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive)
        - total number of primary tillers (column TP) and/or secondary tillers (dead or alive, column TS), and/or total number of tiller present (primary or secondary, dead or alive), and/or total number of tiller with at least 2 ligulated leaves (column TT3F). 
        
    These data are aggregated to estimate primary emission probability, dynamics of mean number of axis present on plants

    Notes :
    - No data found  for direct counting of fertile tillers per plant on tagged plants) 
    - data from date 3 are to be included
    """
    
    fn = str(shared_data(alinea.echap)/'architectural_measurements'/'Tremie13_Tillering_data.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    date_code = {'d1':'2013-02-13', 'd2':'2013-03-29', 'd3': '2012-04-19'}
    TT_code = {'d1':566, 'd2':739, 'd3':915}
    
    # compute emmission probability from data on tagged plants at date 1 and 2
    edata = data[data['Date'] <= 2]
    edata = edata.reset_index()
    emission = emission_probabilities(edata,last='T5')
    
    # compute nff probabilities
    nff_prop = nff_probabilities(data)
    
    # tiller dynamics
    axdyn = axis_dynamics(data)
    
    #estimate ears per plant from plot data
    pdata = Plot_data_Tremie_2012_2013()
    ears_per_plant = pdata['ear_density_at_harvest'] / pdata['mean_plant_density']


    obs = {'emission_probabilities': emission,
           'nff_probabilities': nff_prop,
           'ears_per_plant': ears_per_plant,
           'tillers_per_plant': axdyn, 
           'date_code' : date_code,
           'TT_code' : TT_code}

    return obs
 

def Tillering_data():
    d = Tillering_data_Mercia_Rht3_2010_2011()
    d.update({'Tremie12': Tillering_data_Tremie12_2011_2012(),
             'Tremie13': Tillering_data_Tremie13_2012_2013()})
    return d

#
# ____________________________________________________________________________Pheno data (HS, SSI, GL)
#
#

def Pheno_data(pheno_dict = {}, sources = ['archi_tagged', 'archi_sampled', 'symptom_tagged'], 
                count = 0):
    """ Phenological data (HS, GL, SSI) found for 
        - Treated tagged plants followed for architectural notations
        - Treated tagged plant followed for symptom notation (SSI/GL)
        - some destructive samples (scan /silhouettes data
        
        Details in architectural_measurements R pre-processing scripts
    """

    src = sources[count]
    count += 1
    filename = 'Compil_Pheno_treated_' + src +'.csv'
    filepath = str(shared_data(alinea.echap)/'architectural_measurements'/filename)
    df = pandas.read_csv(filepath, sep = ',', decimal='.')
    df['Date'] = df['Date'].apply(dateparse)
    df_TT = TT_lin()
    df = pandas.merge(df, df_TT)
    pheno_dict[src] = df
    if count < len(sources):
        return Pheno_data(pheno_dict, sources, count)
    else:
        return pheno_dict
    
    
    # file_names = {'archi_tagged': 'Compil_Pheno_treated_archi_tagged.csv',
                  # 'archi_sampled':'Compil_Pheno_treated_archi_sampled.csv',
                  # 'symptom_tagged': 'Compil_Pheno_treated_symptom_tagged.csv'}
    # file_path = {k:str(shared_data(alinea.echap)/'architectural_measurements'/v) for k,v in file_names.iteritems()}
    # return {k:pandas.read_csv(v, sep = ',', decimal='.') for k,v in file_path.iteritems()}
 
#
# Reader for leaf by leaf data (G. Garin April 2015)
# 

def treated_archi_tagged_data(variety = 'Tremie12'):
    filename = variety + '_treated_archi_tagged.csv'
    file_path = shared_data(alinea.echap, filename)
    return pandas.read_csv(file_path, sep = ';')
    
def treated_symptom_tagged_data(variety = 'Tremie12'):
    filename = variety + '_treated_symptom_tagged.csv'
    file_path = shared_data(alinea.echap, filename)
    return pandas.read_csv(file_path, sep = ';')
    
def scan_dimensions_single_date(variety = 'Tremie12', date = '09/05/2012'):
    filename = variety+'_scan_sampled_plants_'+''.join([d[-2:] for d in date.split('/')])+'.txt'
    file_path = filepath = str(shared_data(alinea.echap)/'architectural_measurements'/filename)
    return pandas.read_csv(file_path, sep = "\t")


#----------------------------------------------------- dimension
    
def Dim_data():
    def read_dim_data(source = 'sampled'):
        filepath = shared_data(alinea.echap,
                    'architectural_measurements/Compil_Dim_treated_archi_'+source+'.csv')
        df = pandas.read_csv(filepath, na_values=('NA'), sep=',')
        df['Source'] = source
        return df
    sources = ['sampled', 'tagged']
    return pandas.concat(map(lambda source: read_dim_data(source), sources))
    
# blade dimension data from Corinne/Tino scans in 2009/2010

def blade_dimensions_MerciaRht3_2009_2010():
    """ blade dimension from field experiment at Grignon in 2009-2010.
        Computed by R script (Christian) from scaned leaves database
    """
    fn = shared_data(alinea.echap, 'scaned_blade_dimensions_MerciaRht3_2010.csv')
    return pandas.read_csv(fn,na_values=('NA'),sep=' ')


   
#
# LAI data
#
def PAI_photo_data():
    d = {'Mercia': pandas.DataFrame({'date':['2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-27','2011-04-27','2011-06-01','2011-06-01'], 'PAI_vert_photo':[3.11,3.18,2.93,3.39,3.58,3.38,3.39,3.22,3.24,3.94,4.12,3.82,4.27]}),
    'Rht3': pandas.DataFrame({'date':['2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-27','2011-04-27','2011-04-27','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01'], 'PAI_vert_photo':[4.05,3.16,3.16,3.84,3.27,3.37,3.5,3.65,3,3.03,3.94,4.03,3.9,3.32,3.51,3.52,3.19,3.69]}),
    'Tremie12': pandas.DataFrame({'date':['2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09'], 'PAI_vert_photo':[1.09,0.8,1.02,1.09,1.3,0.88,0.74,0.93,0.87,1.08,0.92,1.22,1.23,3.64,3.45,3.8,3.79,5.69,3.97,3.84,4.35,4.31,4.85,5.15,4.18,4.57,4.52,3.84,3.98,4.2,4.22,3.95,3.99,4.54,4.39]}),
    'Tremie13': pandas.DataFrame({'date':['2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13'], 'PAI_vert_photo':[1.84,2,1.91,1.43,1.76,2.07,1.82,1.52,1.42,1.96,1.73,1.65,1.61,1.67,1.48,1.52,1.75,1.7,2.07,1.76,1.82,1.69,1.64,2.35,2.15,2.41,2.4,2.13,2.55,2.28,2.76,2.69,2.41,3.07,3.01,3.24,3.32,3.1,2.25,2.38,2.26,2.3,2.17,2.91,2.78,2.97,2.89,2.69,3.13,3.3,3.25,3.22,2.68,2.78,2.51,2.83,2.65,2.79,2.79,2.99,2.95,3.12,2.84,2.79,3.08,2.62,2.99,2.64,2.76,2.97,2.56,2.7,2.62,2.67,2.82,2.72,2.65,2.06,1.87,1.96,1.84,1.96,2.27,2.14,2.17,1.94,2.04,2.07,2.02,1.88,1.68,1.84,1.84,1.83,1.8,1.74,1.74,1.85,2.05,2.06,2.05,2.04]})
    }
    #ajout convert date en TT 
    varieties=['Mercia','Rht3','Tremie12','Tremie13']   
    for var in varieties:
        d[var]['TT_date'] = 0.0
        #d[var].date=map(str,d[var].date)
        if var=='Mercia' or var=='Rht3':
            corresp = shared_data(alinea.echap, 'corresp_TTHS_mercia_rht3.csv')
            corresp = pandas.read_csv(corresp, decimal=',', sep=';')
        else:
            corresp = shared_data(alinea.echap, 'corresp_TTHS_'+var+'.csv')
            corresp = pandas.read_csv(corresp, decimal=',', sep=';')
        n=0; m=0
        while n<len(d[var]):
            if corresp['DATE'][m]==d[var]['date'][n]:
                d[var]['TT_date'][n] = corresp['TT'][m]
                n = n + 1
            else :
                m = m + 1
    return d

def LAI_biomasse_data():
    d = {'Tremie12': pandas.DataFrame({'TT_date':[923,923,923,1260,1260,1260,1551,1551,1551], 'LAI_vert_biomasse':[0.97,1.2,1.14,3.61,4.90,4.42,4.46,4.74,4.08]}),
    'Tremie13': pandas.DataFrame({'TT_date':[941,941,941,1060,1060,1060], 'LAI_vert_biomasse':[(4.45*0.8),(4.59*0.8),(5.27*0.8),(5.63*0.8),(5.29*0.8),(6.55*0.8)]})
    }
    return d
    
def _nmax(sub):
    sub['nmax'] = sub['id_Feuille'].max()
    return sub

def scan_dates(variety = 'Tremie12'):
    if variety == 'Tremie12':
        return ['09/03/2012', '02/04/2012', '11/04/2012', '09/05/2012']
    elif variety == 'Tremie13':
        return ['22/04/2013', '03/05/2013']
    else:
        raise ValueError("No scan available for this variety")
        
def scan_dep_Tremie12():
    '''
    Leaf scan data for Tremie 2011-2012
    
    Found data are:
        - Manip 'archi/silhouette': Scan des feuilles aux dates 09/03/2012 (HS 7.4)
                                                     11/04/2012 (HS 10.5)
                                                     09/05/2013 (HS 13.06)
                                files : C:\Users\echap\Desktop\data_mariem\manipeEchap\Manipe_Boigneville20112012\1.Essai_Tremie\1.MesuresINRA\0.DatesBrute\D3_11.04.2012\Scanne_11042012\
                                C:\Users\echap\Desktop\data_mariem\manipeEchap\Manipe_Boigneville20112012\1.Essai_Tremie\1.MesuresINRA\0.DatesBrute\D4_09.05.2012\scannes_09.05.2012\
                                Avec fichier de notations associe
                                C:\Users\echap\Desktop\ECHAP_ARVALIS\Architecture\Archi pour publi\Tremie 2012\2. Longueur-largeur-surface 2N et DFE_ 2012.xlsx
            stat = 1 : feuille entiere et en bon etat
            stat = 2 : longueur correcte mais surface affectee
            stat = 3 : feuille abimee OU en croissance
            numerotation des feuilles coherente avec photos et vrai numero determine (plantes baguees)
            donnes dispo dans share/data raw_srdb_BoignevilleTremie12.csv. 

    notes:
    -target dates of interception are 13/04 (HS 10.67) and 11/05 (HS 13.4)
    - Pour avoir les donnnes utilisees dans le fichier interception, il faut exclure les feueuille stat=3 et appliquer fraction HS a la deniere feuille emergee

    '''
    data_file = shared_data(alinea.echap, 'raw_srdb_BoignevilleTremie12.csv')
    df = pandas.read_csv(data_file, decimal=',', sep=';')
    df = df.loc[:,:'stat']

    df = df[df['id_Axe']=='MB']
    df = df.sort(['HS','id_Feuille'])
    grouped = df.groupby(['HS','N plante'], as_index=False)
    df = grouped.apply(_nmax)
    df['ntop_cur'] = df['nmax'] - df['id_Feuille'] + 1
    
    # pas de donnees observees pour la feuille 11 a HS = 10.15
    # on cree donc une copie de id_feuille=10 pour ce HS que l'on change par id_Feuille=11 et les autres colonnes vides
    # ajout ligne id_feuille 11 pour HS 10.15
    df.ntop_cur[df.HS==10.15] = df.ntop_cur+1
    
    return df
    
def scan_dep_Tremie13():
    '''
    Scan des feuilles non depouillees pour le 22/04/2012 et le 03/05/2012
    C:\Users\echap\Desktop\ECHAP_ARVALIS\Architecture\Archi Anne_2012-2013\Analyse Scan LAI\
    '''
    data_file = shared_data(alinea.echap, 'Tremie13_scan.csv')
    df = pandas.read_csv(data_file, decimal=',', sep=';')
    
    df = df[df['id_Axe']=='MB']
    df = df.sort(['HS','id_Feuille'])
    grouped = df.groupby(['HS','N plante'], as_index=False)
    df = grouped.apply(_nmax)
    df['ntop_cur'] = df['nmax'] - df['id_Feuille'] + 1
    
    return df

def treatment_scan(name='Tremie13'): #Tremie13, pas de colonne stat => stat='all'

    if name is 'Tremie12':
        df = scan_dep_Tremie12()
    elif name is 'Tremie13':
        df = scan_dep_Tremie13()

    # on prend maintenant tjrs stat 1, 2 ou 3, plus besoin du parametre stat='all'
    '''if stat is 'all':
        print 'On prend l\'ensemble des donnees => statuts 1, 2 ou 3'
    else:
        df = df[df['stat']<3]
        print 'On prend seulement les donnees ayant les statuts 1 ou 2'
    '''    
    
    dater = []; moy_feuille = []; ntop = []; A_bl = []
    lst_date = df['HS'].unique()
    for HS in lst_date:
        dt = df[df['HS']==HS]
        
        # tableau avec notation depuis le haut (ntop_cur) et moyenne notation depuis le bas (moyenne id_Feuille)
        data = dt.groupby(['ntop_cur'], as_index=False).mean()  
        # tableau avec notation depuis le bas (id_Feuille) et moyenne notation depuis le haut (moyenne ntop_cur)     
        #data = dt.groupby(['id_Feuille'], as_index=False).mean()    

        nbr = data['ntop_cur'].count(); cpt = 0
        while cpt<nbr:
            dater.append(HS)
            moy_feuille.append(data['id_Feuille'][cpt])
            ntop.append(data['ntop_cur'][cpt])
            A_bl.append(data['A_bl'][cpt])
            cpt += 1
            
    surfSet = zip(dater,moy_feuille,ntop,A_bl)
    df_fin = pandas.DataFrame(data = surfSet, columns=['HS', 'moyenne id_Feuille', 'ntop_cur', 'Area A_bl'])
    return df_fin    
    
    # barplot chaque feuille par date
    '''for date in ['mercredi 9 mai 2012']:
        dt = df[df['date prelev']==date] 
        for leaf in dt['id_Feuille'].unique():  
            plt.figure()
            dy = dt[dt['id_Feuille']==leaf] 
            #dt.plot('id_Feuille', 'A_bl', label='A_bl', kind='bar')
            dy['A_bl'].hist()
            #plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
            plt.suptitle("scan "+date+", num leaf : "+str(leaf))
            #plt.xlabel("num feuille")'''
        
    # plot moyenne feuille par date
    '''dfe = df.groupby(['date prelev','id_Feuille'], as_index=False).mean()
    for date in ['lundi 2 avril 2012', 'mercredi 11 avril 2012', 'mercredi 9 mai 2012', 'vendredi 9 mars 2012']: 
        plt.figure()
        dt = dfe[dfe['date prelev']==date]
        dt.plot('id_Feuille', 'A_bl', '-or', label = 'A_bl')
        dt.plot('id_Feuille', 'A_bl_green', '-og', label = 'A_bl_green')
        #titre, legende et titre axe x
        plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size':9})
        plt.suptitle("scan "+date)
        plt.xlabel("num feuille")'''

    
#
# TC data
#   
def TC_data():
    d = {'Mercia_0':pandas.DataFrame({'date':['2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01'],'TC':[0.17,0.16,0.13,0.16,0.12,0.13,0.13,0.14,0.10,0.66,0.66,0.66,0.59,0.58,0.58,0.66,0.65,0.65,0.64,0.69,0.72,0.73,0.72,0.67,0.69,0.81,0.80,0.71,0.52,0.61,0.51,0.53,0.62,0.56,0.68,0.67,0.77]}),
         'Rht3_0': pandas.DataFrame({'date':['2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-01-19','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-04-27','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01'],'TC':[0.11,0.10,0.11,0.11,0.10,0.10,0.11,0.12,0.12,0.85,0.79,0.79,0.78,0.80,0.80,0.8,0.8,0.81,0.8,0.8,0.82,0.84,0.8,0.85,0.84,0.86,0.83,0.81,0.81,0.55,0.55,0.57,0.49,0.56,0.64,0.64,0.62,0.68]}),
         'Tremie12_0': pandas.DataFrame({'date':['2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09'],'TC':[0.27,0.24,0.23,0.29,0.30,0.28,0.30,0.31,0.30,0.29,0.28,0.33,0.32,0.36,0.58,0.51,0.37,0.49,0.56,0.62,0.54,0.53,0.62,0.57,0.49,0.69,0.72,0.36,0.58,0.51,0.37,0.49,0.56,0.62,0.54,0.53,0.62,0.57,0.49,0.69,0.72]}),
         'Tremie13_0': pandas.DataFrame({'date':['2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13'],'TC':[0.54,0.53,0.51,0.48,0.55,0.51,0.52,0.53,0.52,0.51,0.55,0.55,0.56,0.59,0.51,0.55,0.54,0.52,0.47,0.50,0.53,0.59,0.55,0.51,0.47,0.55,0.51,0.52,0.51,0.5,0.49,0.47,0.44,0.45,0.63,0.57,0.54,0.46,0.56,0.58,0.62,0.62,0.62,0.6,0.62,0.63,0.67,0.63,0.76,0.71,0.7,0.65,0.73,0.67,0.75,0.69,0.74,0.76,0.65,0.62,0.63,0.62,0.6,0.73,0.69,0.66,0.7,0.73,0.62,0.61,0.75,0.85,0.59,0.51,0.55,0.54,0.58,0.46,0.7,0.49,0.47,0.56,0.44,0.57,0.51,0.5,0.63,0.67,0.42,0.5,0.54,0.56,0.41,0.58,0.63,0.66,0.62,0.52,0.57,0.45,0.55,0.57,0.5,0.46,0.48,0.5,0.51,0.54,0.44,0.55,0.67,0.64,0.56,0.58,0.55,0.55,0.64,0.60,0.62,0.58,0.57,0.56,0.64,0.57,0.5,0.52,0.57,0.56,0.55,0.57,0.51,0.55,0.52,0.52,0.49,0.56,0.52,0.4,0.55,0.58,0.51,0.53,0.53,0.51,0.46,0.53,0.53]}),
         'Mercia_57':pandas.DataFrame({'date':['2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-27','2011-04-27','2011-06-01','2011-06-01'],'TC':[0.94,0.95,0.93,0.96,0.96,0.96,0.96,0.95,0.95,0.97,0.98,0.97,0.98]}),
         'Rht3_57': pandas.DataFrame({'date':['2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-18','2011-04-27','2011-04-27','2011-04-27','2011-06-01','2011-06-01','2011-06-01','2011-06-01','2011-06-01'],'TC':[0.98,0.95,0.95,0.97,0.95,0.96,0.96,0.97,0.94,0.94,0.97,0.98,0.97,0.95,0.96,0.96,0.95,0.97]}),
         'Tremie12_57': pandas.DataFrame({'date':['2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-03-09','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-04-11','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09','2012-05-09'],'TC':[0.64,0.53,0.61,0.64,0.70,0.56,0.50,0.58,0.55,0.64,0.58,0.68,0.68,0.97,0.96,0.97,0.97,0.99,0.97,0.97,0.98,0.98,0.99,0.99,0.98,0.99,0.99,0.97,0.98,0.98,0.98,0.97,0.98,0.99,0.98]}),
         'Tremie13_57': pandas.DataFrame({'date':['2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-17','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-04-26','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-05-29','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13','2013-06-13'],'TC':[0.82,0.85,0.83,0.73,0.81,0.85,0.82,0.76,0.73,0.84,0.8,0.78,0.78,0.79,0.75,0.76,0.80,0.79,0.85,0.81,0.82,0.79,0.78,0.89,0.86,0.89,0.89,0.86,0.91,0.88,0.92,0.92,0.89,0.94,0.94,0.95,0.95,0.94,0.88,0.89,0.88,0.88,0.87,0.93,0.92,0.94,0.93,0.92,0.95,0.95,0.95,0.95,0.92,0.92,0.90,0.93,0.91,0.93,0.93,0.94,0.94,0.95,0.93,0.93,0.94,0.91,0.94,0.91,0.92,0.94,0.91,0.92,0.91,0.92,0.93,0.92,0.92,0.85,0.82,0.84,0.82,0.84,0.88,0.86,0.87,0.84,0.85,0.85,0.85,0.83,0.79,0.82,0.82,0.82,0.81,0.8,0.8,0.82,0.85,0.85,0.85,0.85]})}
    #ajout convert date en TT 
    varieties=['Mercia_0','Rht3_0','Tremie12_0','Tremie13_0','Mercia_57','Rht3_57','Tremie12_57','Tremie13_57']   
    for var in varieties:
        d[var]['TT'] = 0.0
        #d[var].date=map(str,d[var].date)
        if var=='Mercia_0' or var=='Mercia_57' or var=='Rht3_0' or var=='Rht3_57':
            corresp = shared_data(alinea.echap, 'corresp_TTHS_mercia_rht3.csv')
            corresp = pandas.read_csv(corresp, decimal=',', sep=';')
        elif var=='Tremie12_0' or var=='Tremie12_57':
            corresp = shared_data(alinea.echap, 'corresp_TTHS_Tremie12.csv')
            corresp = pandas.read_csv(corresp, decimal=',', sep=';')
        elif var=='Tremie13_0' or var=='Tremie13_57':
            corresp = shared_data(alinea.echap, 'corresp_TTHS_Tremie13.csv')
            corresp = pandas.read_csv(corresp, decimal=',', sep=';')
        n=0; m=0
        while n<len(d[var]):
            if corresp['DATE'][m]==d[var]['date'][n]:
                d[var]['TT'][n] = corresp['TT'][m]
                n = n + 1
            else :
                m = m + 1
    return d
#
# mat data
#
def mat_data():
    d = {'Mercia':pandas.DataFrame({'TT':[1137,1281,1796],'light0':[],'light20':[]}),
         'Rht3': pandas.DataFrame({'TT':[1137,1281,1796],'light0':[],'light20':[]}),
         'Tremie': pandas.DataFrame({'TT':[1260,1550],'light0':[0.177,0.045],'light20':[0.633,0.141]})}
    return d
    
#
# leaf surface cm2
#
def surface_data():
    d = {'Tremie12': pandas.DataFrame({'TT':[905.5,1160.1,1239.4,1514.6],'surface':[9.27,21.72,28.03,30.55]}),
    'Tremie13': pandas.DataFrame({'TT':[941,1096.3],'surface':[9.35,13.15]})}
    return d


    
#-------------------------------------------------------------------------------  
# Angles / formes a plat  

def mean_nff():
    def _nff(ms_nff_probas):
        return sum([int(k)*v for k,v in ms_nff_probas.iteritems()])
    return {k:_nff(v['nff_probabilities']) for k,v in Tillering_data().iteritems()}
   
def sr_data():
        
    def sr_reader(file) :
        header_row_srdb = ['rankclass','s','r']
        return pandas.read_csv(file, names=header_row_srdb, sep=',', index_col=False, skiprows=1, decimal='.')
    
    srdb = {'Mercia': shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv'),
            'Rht3': shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv'),
            'Tremie12': shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv'),
            'Tremie13': shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv')}

    return {k: sr_reader(srdb[k]) for k in srdb}

def leaf_curvature_data(name='Mercia'):

    def xy_reader(file):
        header_row_xydb = ['variety','variety_code','harvest','plant','rank','ranktop','relative_ranktop','HS','inerv','x','y','hins','side']
        return pandas.read_csv(file, names=header_row_xydb, sep=',', index_col=False, skiprows=1, decimal='.')
        
    def sr_reader(file) :
        header_row_srdb = ['rankclass','s','r']
        return pandas.read_csv(file, names=header_row_srdb, sep=',', index_col=False, skiprows=1, decimal='.')
        
    if name is 'Mercia':
        data_file_xydb = shared_data(alinea.echap, 'architectural_measurements/xydb_Grignon2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv') 
        # XY
        dfxy = xy_reader(data_file_xydb)
        # use Mercia + Rht3 
        #dfxy = dfxy[dfxy['variety'].isin(['Mercia','Rht3'])]
        #dfxy = dfxy.reset_index()
        # SR
        dfsr = sr_reader(data_file_srdb)

    if name is 'Rht3':
        data_file_xydb = shared_data(alinea.echap, 'architectural_measurements/xydb_Grignon2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv') 
        # XY
        dfxy = xy_reader(data_file_xydb)
        #dfxy = dfxy[dfxy['variety'].isin(['Mercia','Rht3'])]
        #dfxy = dfxy.reset_index()
        # SR : use same as Mercia ???
        dfsr = sr_reader(data_file_srdb)
        
    if name is 'Tremie12':
        # fichier angle non dispo encore, on prend mercia en attendant
        data_file_xydb = shared_data(alinea.echap, 'architectural_measurements/xydb_Grignon2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv')
        # XY
        dfxy = xy_reader(data_file_xydb)
        #dfxy = dfxy[dfxy['variety'].isin(['Mercia','Rht3'])]
        #dfxy = dfxy.reset_index()
        # SR
        dfsr = sr_reader(data_file_srdb)
        
    if name is 'Tremie13':
        # fichier angle non dispo encore, on prend mercia en attendant
        data_file_xydb = shared_data(alinea.echap, 'architectural_measurements/xydb_Grignon2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv')
        # XY
        dfxy = xy_reader(data_file_xydb)
        #dfxy = dfxy[dfxy['variety'].isin(['Mercia','Rht3'])]
        #dfxy = dfxy.reset_index()
        # SR
        dfsr = sr_reader(data_file_srdb)
    
    return dfxy, dfsr

# reder for nervs used for validation curvature
    
def xydb_reader(name = 'Mercia'):
    if name in ['Mercia', 'Rht3']:
        data_file_xydb = shared_data(alinea.echap, 'architectural_measurements/xydb_Grignon2010.csv')
        df = pandas.read_csv(data_file_xydb)
        df.loc[:, 'side'] = 1.
    elif name in ['Tremie12', 'Tremie13']:
        data_file_xydb = shared_data(alinea.echap, 'architectural_measurements/xydb_Boigneville_Tremie12_Tremie13.csv')
        df = pandas.read_csv(data_file_xydb)
        # Wrong data for plants 19, 20, 21 on harvest 2 (Redo)
        if name == 'Tremie12':
            df = df[~((df['harvest']==2) & (df['plant'].isin([19, 20, 21])))]
    return df[df['variety'] == name]

# interpolated median leaf from Grignon 2010

def read_trajectories(fn):     
    dat = pandas.read_csv(fn, names=['age','x','y','lindex'], sep=',', index_col=False, skiprows=1, decimal='.')
    ages = set(dat['age'])
    numage = numpy.array(sorted(list(ages)))
    agemed = (numage[1:]  + numage[:-1]) / 2.0
    bins = [numage[0] - 1.0] + agemed.tolist() + [numage[-1] + 1.0]
    dat['age_class'] = pandas.cut(dat['age'], bins,labels=False)
    grouped = dat.groupby(('lindex','age_class'))   
    trajs = {k:[{a:grouped.get_group((k,a)).ix[:,['x','y']] for a in set(dat['age_class'])}] for k in set(dat['lindex'])}
    return trajs, bins

def median_leaf_trajectories():
    """ interpolated xy for upper/lower leaves every 0.5 HS on Grignon 2009-2010 data
    """
    trajs = {'MerciaRht_byleafclass': 'MerciaRht_Grignon2010_byleafclass',
             'Tremie_byleafclass': 'Tremie_Boigneville2012_2013_byleafclass',
             'Soissons_byleafclass':'Soissons_Grignon2010_byleafclass',
             'MerciaRht': 'MerciaRht_Grignon2010',
             'Tremie': 'Tremie_Boigneville2012_2013',
             'Soissons':'Soissons_Grignon2010',
             'Tremie12':'Tremie12',
             'Tremie13':'Tremie13'}
    fn = {k:shared_data(alinea.echap, 'architectural_measurements/median_leaf_trajectories_' + trajs[k] + '.csv') for k in trajs}
    return {k:read_trajectories(fn[k]) for k in trajs}
    
#
# Elaborated data
#    


    
def _add_ghs(df, g):
    #hs_conv = HS_converter[g]
    df['Var'] = g
    #df['HS'] = hs_conv(df['TT'])
    df = df.sort('TT')
    return df

def PlantDensity():
    # Mercia /rht3
    ld = []
    pd = Plot_data_Mercia_Rht3_2010_2011()
    td = Tillering_data_Mercia_Rht3_2010_2011()
    for g in ('Mercia','Rht3'):
        pdata = pd[g]
        tdata = td[g]
        date,TT,density,SD = [],[],[],[]
        events = ['sowing', 'emergence', 'harvest']
        density = [pdata['sowing_density'], 
                   pdata['plant_density_at_emergence'],
                   pdata['ear_density_at_harvest'] / tdata['ears_per_plant']]
        ear_data = numpy.array(pdata['raw_ear_density_at_harvest']) / tdata['ears_per_plant']
        SD = [0,0,numpy.nanstd(ear_data)]           
        for w in events:
            date.append(pdata['code_date'][w])
            TT.append(pdata['TT_date'][pdata['code_date'][w]])
        # add survival data
        for i,d in enumerate(tdata['plant_survival']['Date']):
            lab = 'd%d'%(d+1)
            events.append(lab)
            date.append(tdata['date_code'][lab])
            TT.append(tdata['TT_code'][lab])
            density.append(pdata['plant_density_at_emergence'] * tdata['plant_survival']['viability'][i])
            SD.append(0)
        df = pandas.DataFrame({'event': events, 'date':date, 'TT':TT, 'density':density})
        df['SD']=SD
        df = _add_ghs(df, g)
        ld.append(df)
    # Tremie12
    pdata = Plot_data_Tremie_2011_2012()
    tdata = Tillering_data_Tremie12_2011_2012()
    #hs_conv = HS_converter['Tremie12']
    date,TT,density,SD = [],[],[],[]
    events = ['sowing', 'harvest', 'd3', 'd4', 'arvalis']
    density = [pdata['sowing_density'], 
               pdata['ear_density_at_harvest'] / tdata['ears_per_plant'],
               numpy.mean(pdata['plant_density'][pdata['code_date']['d3']]),
               numpy.mean(pdata['plant_density'][pdata['code_date']['d4']]),
               numpy.mean(pdata['plant_density'][pdata['code_date']['arvalis']])
               ]
    ear_density_new=[]
    for i,item in enumerate(pdata['raw_ear_density_at_harvest']):
        lab = pdata['raw_ear_density_at_harvest'][i]/tdata['ears_per_plant']
        ear_density_new.append(lab)
    SD = [0, numpy.nanstd(ear_density_new),
         numpy.nanstd(pdata['plant_density'][pdata['code_date']['d3']]), 
         numpy.nanstd(pdata['plant_density'][pdata['code_date']['d4']]), 
         numpy.nanstd(pdata['plant_density'][pdata['code_date']['arvalis']])
         ]
    for w in events:
        date.append(pdata['code_date'][w])
        TT.append(pdata['TT_date'][pdata['code_date'][w]])
    df = pandas.DataFrame({'event': events, 'date' :date, 'TT':TT, 'density':density})
    df['SD']=SD
    df = _add_ghs(df, 'Tremie12')
    ld.append(df)
    # Tremie13
    pdata = Plot_data_Tremie_2012_2013()
    tdata = Tillering_data_Tremie13_2012_2013()
    #hs_conv = HS_converter['Tremie13']
    date,TT,density,SD = [],[],[],[]
    events = ['sowing', '2F', 'epis1cm', 'LAI1', 'LAI2']
    density = [pdata['sowing_density'], 
               numpy.mean(pdata['plant_density'][pdata['code_date']['2F']]),
               numpy.mean(pdata['plant_density'][pdata['code_date']['epis1cm']]),
               numpy.mean(pdata['plant_density'][pdata['code_date']['LAI1']]),
               numpy.mean(pdata['plant_density'][pdata['code_date']['LAI2']])
               ]
    SD = [0,
         numpy.nanstd(pdata['plant_density'][pdata['code_date']['2F']]), 
         numpy.nanstd(pdata['plant_density'][pdata['code_date']['epis1cm']]), 
         numpy.nanstd(pdata['plant_density'][pdata['code_date']['LAI1']]),
         numpy.nanstd(pdata['plant_density'][pdata['code_date']['LAI2']])
         ]
    for w in events:
        date.append(pdata['code_date'][w])
        TT.append(pdata['TT_date'][pdata['code_date'][w]])
    df = pandas.DataFrame({'event': events, 'date' :date, 'TT':TT, 'density':density})
    df['SD']=SD
    df = _add_ghs(df, 'Tremie13')
    ld.append(df)
    return pandas.concat(ld)

# converter axis_dyn output -> df
def _axdyn_df(tdb):
    date, TT = [], []
    for i,d in enumerate(tdb['tillers_per_plant']['Date']):
        lab = 'd%d'%(d)
        date.append(tdb['date_code'][lab])
        TT.append(tdb['TT_code'][lab])
    til = tdb['tillers_per_plant']
    return pandas.DataFrame({'date': date, 'TT': TT, 'TP' : til['TP'], 'TS':til['TS'], 'TPS': til['TT'], 'FT':til['FT'], 'TT3F':til['TT3F']})

def tillers_per_plant():
    """ Synthetic tables for tiller counting data
    """
    ld=[]
    for g in ('Mercia','Rht3'):
        tdb = Tillering_data_Mercia_Rht3_2010_2011()[g]
        df = _axdyn_df(tdb)
        df = _add_ghs(df, g)
        ld.append(df)
    #
    tdb =  Tillering_data_Tremie12_2011_2012()
    df = _axdyn_df(tdb)
    df = df[~df['TPS'].isnull()]#remove d3 where no data are available
    # add axe counts from plot data
    pdata = Plot_data_Tremie_2011_2012()
    axd = pdata['axe_density']
    faxd = pdata['fertile_axis_density']
    dfp= pandas.DataFrame({'date' : axd.keys(), 
                'TT': [pdata['TT_date'][k] for k in axd],
                'TPS': [numpy.mean(axd[k]) / numpy.mean(pdata['plant_density'][k]) - 1 for k in axd]})
    dfp['FT'] = numpy.nan
    d = faxd.keys()[0]
    dfp.ix[dfp['date']==d,'FT'] = numpy.mean(faxd[d])  / numpy.mean(pdata['plant_density'][d]) - 1
    df = pandas.concat([df, dfp])
    df = _add_ghs(df, 'Tremie12')
    ld.append(df)
    #
    tdb = Tillering_data_Tremie13_2012_2013()
    df = _axdyn_df(tdb)
    # add ears per plant, to get a point after regression
    pdata = Plot_data_Tremie_2012_2013()
    dfp= pandas.DataFrame({'date' :pdata['code_date']['harvest'],
                           'TT': pdata['TT_date'][pdata['code_date']['harvest']],
                           'FT': pdata['ear_density_at_harvest'] / pdata['mean_plant_density'] - 1}, index=[0])
    df = pandas.concat([df, dfp])
    df = _add_ghs(df, 'Tremie13')
    ld.append(df)    
    return reduce(lambda x,y : pandas.concat([x,y]), ld)
    
def emission_probabilities_table():
    emdb = {}
    emdb['Mercia'] = Tillering_data_Mercia_Rht3_2010_2011()['Mercia']['emission_probabilities']
    emdb['Rht3'] = Tillering_data_Mercia_Rht3_2010_2011()['Rht3']['emission_probabilities']
    emdb['Tremie12'] = Tillering_data_Tremie12_2011_2012()['emission_probabilities']
    emdb['Tremie13'] = Tillering_data_Tremie13_2012_2013()['emission_probabilities']

    dfs = [pandas.DataFrame({'variety':var, 'tiller':probas.keys(), 'probability': probas.values()}) for  var,probas in emdb.iteritems()]
    return pandas.concat(dfs)

def variance(lst):
    """
    Uses standard variance formula (sum of each (data point - mean) squared)
    all divided by number of data points
    """
    mu = numpy.mean(lst)
    return 1.0/(len(lst)-1) * sum([(i-mu)**2 for i in lst]) if len(lst)>1 else 0.
        
def conf_int(lst, perc_conf=95):
    """
    Confidence interval - given a list of values compute the square root of
    the variance of the list (v) divided by the number of entries (n)
    multiplied by a constant factor of (c). This means that I can
    be confident of a result +/- this amount from the mean.
    The constant factor can be looked up from a table, for 95pcent confidence
    on a reasonable size sample (>=500) 1.96 is used.
    """
    from scipy.stats import t
    
    n, v = len(lst), variance(lst)
    c = t.interval(perc_conf * 1.0 / 100, n-1)[1]
    
    return numpy.sqrt(v/n) * c

def haun_stage_aggregated():
    """ Get mean, standard error and confidence interval for HS calibration data """
    def aggregate(data, group = ['TT', 'label', 'nff'],
                  func = numpy.mean, column_name = 'HS_mean'):
        df = data.groupby(group).agg(func).loc[:, 'HS'].to_frame()
        df = df.rename(columns={'HS':column_name})
        return df
    
    # Get and select data
    data = Pheno_data()
    tagged = data['archi_tagged']
    tagged = tagged.loc[tagged['HS'] < tagged['nff'],('label', 'nff', 'TT', 'HS')].dropna()
    sampled = data['archi_sampled']
    sampled = sampled.loc[(sampled['HS'] < sampled['nff']) | (numpy.isnan(sampled['nff'])),('label','TT','HS')].dropna()
    
    # Get mean, standard error and confidence interval for tagged plants by nff
    df_HS_tagged = pandas.concat([aggregate(tagged, ['TT', 'label', 'nff'], numpy.mean, 'HS_mean'),
                                  aggregate(tagged, ['TT', 'label', 'nff'], numpy.nanstd, 'HS_std'),
                                  aggregate(tagged, ['TT', 'label', 'nff'], conf_int, 'HS_conf')], axis=1)
    
    # Get mean, standard error and confidence interval for all NFFs for tagged plants and sampled plants
    (df_HS_tagged_global, df_HS_sampled_global) = map(lambda x: pandas.concat(
                                                  [aggregate(x, ['TT', 'label'], numpy.mean, 'HS_mean'),
                                                   aggregate(x, ['TT', 'label'], numpy.nanstd, 'HS_std'),
                                                   aggregate(x, ['TT', 'label'], conf_int, 'HS_conf')], axis=1),
                                                   (tagged, sampled))
    df_HS_tagged_global['source'] = 'tagged'
    df_HS_sampled_global['source'] = 'sampled'
    df_HS_global = pandas.concat([df_HS_tagged_global, df_HS_sampled_global])
    
    (df_HS_tagged, df_HS_global) = map(lambda x: x.reset_index(), (df_HS_tagged, df_HS_global))
    return df_HS_tagged, df_HS_global

def ssi_aggregated():
    """ Get mean, standard error and confidence interval for ssi calibration data """
    def aggregate(data, group = ['TT', 'label', 'nff'],
                  func = numpy.mean, column_name = 'SSI_mean'):
        df = data.groupby(group).agg(func).loc[:, 'SSI'].to_frame()
        df = df.rename(columns={'SSI':column_name})
        return df
    
    # Get and select data
    data = Pheno_data()
    tagged = data['archi_tagged']
    tagged = tagged.loc[tagged['SSI'] < tagged['nff'],('label', 'nff', 'TT', 'SSI')].dropna()
    sampled = data['archi_sampled']
    sampled = sampled.loc[(sampled['SSI'] < sampled['nff']) | (numpy.isnan(sampled['nff'])),('label','TT','SSI')].dropna()
    
    # Get mean, standard error and confidence interval for tagged plants by nff
    df_SSI_tagged = pandas.concat([aggregate(tagged, ['TT', 'label', 'nff'], numpy.mean, 'SSI_mean'),
                                  aggregate(tagged, ['TT', 'label', 'nff'], numpy.nanstd, 'SSI_std'),
                                  aggregate(tagged, ['TT', 'label', 'nff'], conf_int, 'SSI_conf')], axis=1)
    
    # Get mean, standard error and confidence interval for all NFFs for tagged plants and sampled plants
    (df_SSI_tagged_global, df_SSI_sampled_global) = map(lambda x: pandas.concat(
                                                  [aggregate(x, ['TT', 'label'], numpy.mean, 'SSI_mean'),
                                                   aggregate(x, ['TT', 'label'], numpy.nanstd, 'SSI_std'),
                                                   aggregate(x, ['TT', 'label'], conf_int, 'SSI_conf')], axis=1),
                                                   (tagged, sampled))
    df_SSI_tagged_global['source'] = 'tagged'
    df_SSI_sampled_global['source'] = 'sampled'
    df_SSI_global = pandas.concat([df_SSI_tagged_global, df_SSI_sampled_global])
    
    (df_SSI_tagged, df_SSI_global) = map(lambda x: x.reset_index(), (df_SSI_tagged, df_SSI_global))
    return df_SSI_tagged, df_SSI_global
    

  
def aggregateGL(data, group = ['TT', 'label', 'nff'],
              func = numpy.mean, column_name = 'GL_mean'):
    df = data.groupby(group).agg(func).loc[:, 'GL'].to_frame()
    df = df.rename(columns={'GL':column_name})
    return df    
 
def aggregate_GLsource(data, by=['TT', 'label', 'nff'], source='unknown'):
    aggregated = pandas.concat([aggregateGL(data, by, numpy.mean, 'GL_mean'),
                                aggregateGL(data, by, numpy.nanstd, 'GL_std'),
                                aggregateGL(data, by, conf_int, 'GL_conf')], axis=1)
    aggregated['source'] = source
    return aggregated 
    
def green_leaves_aggregated():
    """ Get mean, standard error and confidence interval for GL calibration data """

    # Get and select data
    data = Pheno_data()
    tagged = data['archi_tagged']
    tagged_obs = tagged.loc[~numpy.isnan(tagged['GL']), ('label', 'nff', 'TT', 'GL')]

    sampled = data['archi_sampled']
    sampled_obs_nff = sampled.loc[(~numpy.isnan(sampled['GL'])) &
                                  (~numpy.isnan(sampled['nff'])),
                                  ('label', 'nff', 'TT', 'GL')]
    sampled_obs_global = sampled.loc[~numpy.isnan(sampled['GL']), ('label', 'TT', 'GL')]

    symptom = data['symptom_tagged']
    symptom_obs = symptom.loc[~numpy.isnan(symptom['GL']), ('label', 'nff', 'TT', 'GL')]
    symptom_obs_ap = symptom.loc[(~numpy.isnan(symptom['GLap'])) &
                             (~numpy.isnan(symptom['SSIap'])), ('label', 'TT', 'GLap', 'SSIap')]
    symptom_obs_ap = symptom_obs_ap.rename(columns = {'GLap':'GL', 'SSIap':'SSI'})
    
    # Get mean, standard error and confidence interval of observed GL for plants with known nff
    sources = {'tagged':tagged_obs, 'sampled':sampled_obs_nff, 'symptom':symptom_obs}
    df_GL_obs_nff = pandas.concat([aggregate_GLsource(sources[k], ['TT', 'label', 'nff'], k) for k in sources])
    df_GL_obs_nff = df_GL_obs_nff.reset_index()
  
    # Get mean, standard error and confidence interval of observed GL for all NFFs for 
    # tagged plants, sampled plants and symptom plants
    sources = {'tagged':tagged_obs, 'sampled':sampled_obs_global, 'symptom':symptom_obs, 'symptom_ap':symptom_obs_ap}
    df_GL_obs_global = pandas.concat([aggregate_GLsource(sources[k], ['TT', 'label'], k) for k in sources])
    df_GL_obs_global = df_GL_obs_global.reset_index()
    df_GL_obs_global['nff'] = numpy.nan
    
    return df_GL_obs_nff, df_GL_obs_global
    
def green_leaves_estimated(HS_fit):

    def estimate_HS(label, nff, TT):
        return min(nff, HS_fit[label].HS(TT, nff))
    
    def estimate_GL(data_est):
        fun_HS = numpy.frompyfunc(estimate_HS, 3, 1)
        if 'nff' not in data_est.columns:
            data_est['nff'] = [numpy.nan for i in range(len(data_est))]
        data_est['GL'] = fun_HS(data_est['label'], data_est['nff'], data_est['TT']) - data_est['SSI']
        data_est['GL'] = data_est['GL'].astype(float)
        return data_est
 
    # Get and select data
    data = Pheno_data()
    tagged = data['archi_tagged']
    tagged_est_nff = tagged.loc[(numpy.isnan(tagged['GL'])) & 
                                (~numpy.isnan(tagged['nff'])) & 
                                (~numpy.isnan(tagged['SSI'])), 
                                ('label', 'nff', 'TT', 'GL', 'SSI')]
    tagged_est_global = tagged.loc[(numpy.isnan(tagged['GL'])) & 
                                    (numpy.isnan(tagged['nff'])) &
                                    (~numpy.isnan(tagged['SSI'])), 
                                    ('label', 'TT', 'GL', 'SSI')]
    sampled = data['archi_sampled']
    sampled_est_global = sampled.loc[(numpy.isnan(sampled['GL'])) & 
                                    (numpy.isnan(sampled['nff'])) &
                                    (~numpy.isnan(sampled['SSI'])), 
                                    ('label', 'TT', 'GL', 'SSI')]
    symptom = data['symptom_tagged']
    symptom_est = symptom.loc[(numpy.isnan(symptom['GL'])) &
                             (~numpy.isnan(symptom['SSI'])), ('label', 'nff', 'TT', 'GL', 'SSI')]
    symptom_est_ap = symptom.loc[(numpy.isnan(symptom['GLap'])) &
                                (~numpy.isnan(symptom['SSIap'])), ('label', 'TT', 'GLap', 'SSIap')]
    symptom_est_ap = symptom_est_ap.rename(columns = {'GLap':'GL', 'SSIap':'SSI'})
    
    # Get mean, standard error and confidence interval of estimated GL for plants with known nff
    (tagged_est_nff, symptom_est) = map(lambda x: estimate_GL(x), (tagged_est_nff, symptom_est))
    sources = {'tagged':tagged_est_nff, 'symptom':symptom_est}
    df_GL_est_nff = pandas.concat([aggregate_GLsource(sources[k], ['TT', 'label', 'nff'], k) for k in sources])
    df_GL_est_nff = df_GL_est_nff.reset_index()
    # Get mean, standard error and confidence interval of estimated GL for all NFFs for 
    # tagged plants, sampled plants and symptom plants
    (tagged_est_global, sampled_est_global, 
    symptom_est, symptom_est_ap) = map(lambda x: estimate_GL(x), 
                                        (tagged_est_global, sampled_est_global,
                                        symptom_est, symptom_est_ap))
    sources = {'tagged':tagged_est_global, 'sampled':sampled_est_global, 'symptom':symptom_est, 'symptom_ap':symptom_est_ap}
    df_GL_est_global = pandas.concat([aggregate_GLsource(sources[k], ['TT', 'label'], k) for k in sources])
    df_GL_est_global = df_GL_est_global.reset_index()
    df_GL_est_global['nff'] = numpy.nan
                                     
    return df_GL_est_nff, df_GL_est_global
    

 
def add_HS_green_leaves(data, HS_fit):
    """  
    """
       
    def get_HS(label, nff, TT):
        return HS_fit[label].HS(TT, nff)
        
    def get_HS_flag(label, nff):
        return HS_fit[label].HSflag(nff)
        
    def get_TT_flag(label, nff):
        return HS_fit[label].TTflag(nff)
        
    fun_HS = numpy.frompyfunc(get_HS, 3, 1)
    fun_HSflag = numpy.frompyfunc(get_HS_flag, 2, 1)
    fun_TT = numpy.frompyfunc(get_TT_flag, 2, 1)

    data['HS'] = fun_HS(data['label'], data['nff'], data['TT'])
    data['HSflag'] = fun_HSflag(data['label'], data['nff'])
    data['TTflag'] = fun_TT(data['label'], data['nff'])
    return data

    
def dimensions_aggregated(df_dim = Dim_data(), along = 'rank'):
    def _aggregate(data, group = ['label', 'Source', along], func = numpy.mean, column_name = 'mean'):
        df = data.groupby(group).agg(func).loc[:, ['L_blade', 'W_blade', 'A_blade',
                                                    'L_sheath', 'W_sheath', 'L_internode', 'H_col']]
        df = df.rename(columns={col:col+'_'+column_name for col in df.columns})
        return df.reset_index()

    def get_aggregate(data, group = ['label', 'Source', along]):
        funcs = {'mean':numpy.mean, 'std':numpy.std, 'conf':conf_int}
        df_ag = None
        for name, func in funcs.iteritems():
            if df_ag is None:
                df_ag = _aggregate(data, group=group, func=func, column_name=name)
            else:
                df_ag = pandas.merge(df_ag, _aggregate(data, group=group, 
                                                        func=func, column_name=name),
                                       on=group)
        return df_ag

    df_ag = get_aggregate(df_dim, group = ['label', 'Source', along])
    df_dim = df_dim[~numpy.isnan(df_dim['nff'])]
    if len(df_dim)>0:
        df_ag_nff = get_aggregate(df_dim, group = ['label', 'Source', 'nff', along])
    else:
        df_ag_nff = None
    return df_ag, df_ag_nff
    
class ReconstructionData(object):

    def __init__(self):
        self.Plot_data = Plot_data()
        self.Tillering_data = Tillering_data()
        self.Pheno_data = Pheno_data()
        self.Dimension_data = Dim_data()
        self.xy_data = median_leaf_trajectories()
        self.sr_data = sr_data()
        
    def save(self, filename):
        with open(filename, 'w') as output:
            pickle.dump(self, output)
 
def reconstruction_data(reset=False):
    filename = str(shared_data(alinea.echap)/'architectural_ReconstructionData.pckl')
    if not reset:
        try:
            with open(filename) as input:
                return pickle.load(input)
        except:
            pass
    Data = ReconstructionData()
    Data.save(filename)
    return Data
 
class ValidationData(object):
    
    def __init__(self):
        self.PlantDensity = PlantDensity()
        self.tillers_per_plant = tillers_per_plant()
        self.emission_probabilities = emission_probabilities_table()    
        self.haun_stage = haun_stage_aggregated()
        self.ssi = ssi_aggregated()
        self.dimensions = dimensions_aggregated()
        self.green_leaves = green_leaves_aggregated()
        
    def save(self, filename):
        with open(filename, 'w') as output:
            pickle.dump(self, output)
            
    def update_HS(self, HS_fit):
        self.green_leaves_estimated = green_leaves_estimated(HS_fit)
        self.green_leaves = [add_HS_green_leaves(df, HS_fit) for df in self.green_leaves]
        self.green_leaves_estimated = [add_HS_green_leaves(df, HS_fit) for df in self.green_leaves_estimated]
        
    
def validation_data(reset=False, HS_fit=None):
    filename = str(shared_data(alinea.echap)/'architectural_ValidationData.pckl')
    Data = None
    if not reset:
        try:
            with open(filename) as input:
                Data =  pickle.load(input)
        except:
            pass
    if Data is None:
        Data = ValidationData()
    if HS_fit is not None:
        Data.update_HS(HS_fit)
    Data.save(filename)
    return Data     

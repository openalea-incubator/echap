import pandas
import numpy

from openalea.deploy.shared_data import shared_data
import alinea.echap
from alinea.adel.plantgen.plantgen_interface import read_plantgen_inputs



#----------------------------------------------------- Plantgen
def plantgen_as_dict(inputs, dynT, dimT):
    d={}
    d['dynT_user'], d['dimT_user'], d['plants_number'],d['plants_density'], d['decide_child_axis_probabilities'], d['MS_leaves_number_probabilities'], d['ears_density'], d['GL_number'], d['delais_TT_stop_del_axis'], d['TT_col_break'],d['inner_params'] =  read_plantgen_inputs(inputs, dynT, dimT)
    return d
    
def Mercia_2010_plantgen():
    dynT = shared_data(alinea.echap, 'Mercia_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Mercia_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
    
# composite ---
def Mercia_2010_nff11_plantgen():
    dynT = shared_data(alinea.echap, 'Mercia_dynT11_user.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT11_user.csv')
    inputs = shared_data(alinea.echap, 'Mercia_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
def Mercia_2010_nff12_plantgen():
    dynT = shared_data(alinea.echap, 'Mercia_dynT12_user.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT12_user.csv')
    inputs = shared_data(alinea.echap, 'Mercia_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
def Mercia_2010_nff13_plantgen():
    dynT = shared_data(alinea.echap, 'Mercia_dynT13_user.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT13_user.csv') #pas donnees excel => ajout de la ligne 13 a dimT12
    inputs = shared_data(alinea.echap, 'Mercia_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
    
def Rht3_2010_nff11_plantgen():
    dynT = shared_data(alinea.echap, 'Rht3_dynT11_user.csv')
    dimT = shared_data(alinea.echap, 'Rht3_dimT11_user.csv')
    inputs = shared_data(alinea.echap, 'Rht3_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
def Rht3_2010_nff12_plantgen():
    dynT = shared_data(alinea.echap, 'Rht3_dynT12_user.csv')
    dimT = shared_data(alinea.echap, 'Rht3_dimT12_user.csv')
    inputs = shared_data(alinea.echap, 'Rht3_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)
#---
    
def Rht3_2010_plantgen():
    dynT = shared_data(alinea.echap, 'Rht3_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Rht3_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Rht3_plantgen_inputs_MIN.py')
    return plantgen_as_dict(inputs, dynT, dimT)

def Tremie_2011_plantgen():
    dynT = shared_data(alinea.echap, 'Tremie1_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Tremie1_dimT_user.csv')
    inputs = shared_data(alinea.echap, 'Tremie1_plantgen_inputs_MINnew.py')
    return plantgen_as_dict(inputs, dynT, dimT)
'''
def Tremie_2012_plantgen():
    #dynT = shared_data(alinea.echap, 'Tremie2_dynT_user.csv')
    dimT = shared_data(alinea.echap, 'Tremie2_dimT_user.csv')
    #inputs = shared_data(alinea.echap, 'Tremie2_plantgen_inputs_MINnew.py')
    return plantgen_as_dict(inputs, dynT, dimT)
'''   
def HS_data():
    fn = shared_data(alinea.echap, 'HS_data_Mercia_Rht3_2010_2011.csv')
    #fn = shared_data(alinea.echap, 'HS_data_Tremie1_2011_2012.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    grouped = data.groupby(['variety'],as_index=False)
    return grouped.aggregate('mean')
#
# Plot data
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
        'sowing_density': 220, 
        'plant_density_at_emergence' : 203, 
        'ear_density_at_harvest' : 444,
        'raw_ear_density_at_harvest':[507,440,427,357,430,410,427,480,363, 350, 497, 410, 493, 340, 407, 467, 490, 433, 547, 450, 527, 427, 483, 493],
        'inter_row': 0.15},
        'Rht3': {
        'code_date': {'sowing':'2010-10-15', 'emergence':'2010-11-02', 'harvest':'2011-06-20'},
        'sowing_density': 220, 
        'plant_density_at_emergence' : 211, 
        'ear_density_at_harvest' : 384,
        'raw_ear_density_at_harvest':[440,420,330,433,347,410,387,367,360,397,357,377,427,367,380,347,330,420,423,397,383,367,377,377],
        'inter_row': 0.15}}
    return d

def Plot_data_Tremie_2011_2012():
    """
    Plot data for Boigneville 2011-2012 
    
    Notes:
    - No data found for density of emergence as this trial was followed only after winter (winter trial was on appache)
    - plant density data at date 3 (11/04/2012) and 4 comes from axe counting/ LAI measurement data taken on  5 rank * 60 cm prelevement)
    """
    d = {
    'code_date':{'sowing': '2011-10-21','emergence': '2011-11-03', 'harvest':'2012-06-19'}, 
    'sowing_density': 280,
    'ear_density_at_harvest': 492,
    'raw_ear_density_at_harvest':[479, 490, 598, 608, 538, 503, 493, 430, 458, 437, 486, 489, 465, 406],
    'inter_row':0.143,
    'plant_density':{'2012-03-20': [238, 304, 287, 237, 290, 301, 290, 287, 273],
                     '2012-04-11': [301, 273, 315], 
                     '2012-05-09':[257, 301, 263]},
    'axe_density': {'2012-04-11': [918, 956, 979],
                    '2012-05-09':[601, 585, 506]}, 
    'fertile_axis_density':{'2012-05-09':[545, 569, 443]} 
    }
    return d
    
def Plot_data_Tremie_2012_2013():
    """
    Plot data for Boigneville 2012-2013 
    
    Notes
    - no date found for plant density counts at stage 2F and epis1cm (estimation epis 1cm : 9/04/2012)
    - no date found for countings of ear density
    - plant density and ear density were estimated on 2 ranks * 1m plots at stage 2F and epis1cm (0.3 m2), and on plots of 4 ranks * 0.6 m for LAI plots (0.36 m2) 
    - Plant density data measured for LAI estimation seems bugy and more compatible with a 5 rank * 0.6 m plots dimension (ie density and LAI should be multiplied by 0.8)
    """
    d = {
   'code_date':{'sowing': '2012-10-29','emergence': '2012-11-19'},
   'sowing_density': 300,
   'plant_density':{'2F': [237, 287, 217, 237, 293, 220, 253, 213, 220, 253],
                    'epis1cm': [203, 203, 193, 207, 223, 203, 207, 207, 200, 197], 
                    '2012-04-22':[328, 322, 356],
                    '2012-05-13':[272, 361, 350]},
   'ear_density_at_harvest' : 676,
   'raw_ear_density_at_harvest':[643, 580, 693, 813, 663, 670, 693, 763, 567, 590, 707, 637, 617, 747, 760, 693, 653, 670],
    'inter_row':0.15
    }
    return d
   
#
# LAI data
#
def PAI_data():
    d = {'Mercia':pandas.DataFrame({'TT':[1137,1281,1796],'PAI':[3.27,4.03,4.05]}),
         'Rht3': pandas.DataFrame({'TT':[1137,1281,1796],'PAI':[3.40,3.96,3.45]}),
         'Tremie': pandas.DataFrame({'TT':[923,1260,1550],'PAI':[1.1,4.31,4.43]})}
         #'Tremie': pandas.DataFrame({'TT':[1260,1550],'PAI':[3.8,4.7]})} #Publi Corinne
    return d
    
#
# TC data
#
def TC_data():
    d = {'Mercia_0':pandas.DataFrame({'TT':[475,1137,1281,1796],'TC':[0.14,0.63,0.73,0.61],'1-TC':[0.86,0.37,0.27,0.39]}),
         'Rht3_0': pandas.DataFrame({'TT':[475,1137,1281,1796],'TC':[0.11,0.80,0.83,0.59],'1-TC':[0.89,0.2,0.17,0.41]}),
         'Tremie_0': pandas.DataFrame({'TT':[923,1260,1550],'TC':[0.3,0.58,0.63],'1-TC':[0.7,0.42,0.37]}),
         'Mercia_57':pandas.DataFrame({'TT':[1137,1281,1796],'TC':[0.95,0.98,0.98]}),
         'Rht3_57': pandas.DataFrame({'TT':[1137,1281,1796],'TC':[0.956,0.975,0.96]}),
         'Tremie_57': pandas.DataFrame({'TT':[923,1260,1550],'TC':[0.59,0.97,0.985],'1-TC':[0.41,0.03,0.015]})}
    return d
      
#
# mat data
#
def mat_data():
    d = {'Mercia':pandas.DataFrame({'TT':[1137,1281,1796],'light0':[],'light20':[]}),
         'Rht3': pandas.DataFrame({'TT':[1137,1281,1796],'light0':[],'light20':[]}),
         'Tremie': pandas.DataFrame({'TT':[1260,1550],'light0':[0.177,0.045],'light20':[0.633,0.141]})}
    return d
    
#-------------------------------------------------------------------------------  
# Angles      
def leaf_curvature_data(name='Mercia'):
    import re
    
    if name is 'Mercia':
        data_file_xydb = shared_data(alinea.echap, 'xydb_GrignonMercia2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv') 
        # XY
        header_row_xydb = ['variety','variety_code','harvest','plant','rank','ranktop','relative_ranktop','HS','inerv','x','y']
        dfxy = pandas.read_csv(data_file_xydb, names=header_row_xydb, sep=',', index_col=0, skiprows=1, decimal='.')
        # filtre sur la variete pour Mercia/Rht3
        dfxy = dfxy.reset_index()
        dfxy1 = dfxy[dfxy['variety']=='Mercia'] 
        dfxy2 = dfxy[dfxy['variety']=='Rht3'] 
        dfxy = dfxy1.append(dfxy2)
        #dfxy = dfxy[dfxy['variety']=='Mercia']
        dfxy['rankclass'] = 1
        dfxy['rankclass'][dfxy['ranktop'] <= 4] = 2
       
    if name is 'Rht3':
        data_file_xydb = shared_data(alinea.echap, 'xydb_GrignonMercia2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv') 
        # XY
        header_row_xydb = ['variety','variety_code','harvest','plant','rank','ranktop','relative_ranktop','HS','inerv','x','y']
        dfxy = pandas.read_csv(data_file_xydb, names=header_row_xydb, sep=',', index_col=0, skiprows=1, decimal='.')
        # filtre sur la variete Mercia seulement
        dfxy = dfxy.reset_index()
        dfxy = dfxy[dfxy['variety']=='Mercia'] 
        # feuille du haut = feuille du bas pr coller a l'archi de rht3
        dfxy['rankclass'] = 2
            
    if name is 'Tremie':
        # fichier angle non dispo encore, on prend mercia en attendant
        data_file_xydb = shared_data(alinea.echap, 'xydb_GrignonMercia2010.csv') 
        data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonMercia2010.csv')
        header_row_xydb = ['variety','variety_code','harvest','plant','rank','ranktop','relative_ranktop','HS','inerv','x','y']
        #data_file_xydb = shared_data(alinea.echap, 'xydb_GrignonTremie2011.csv') 
        #data_file_srdb = shared_data(alinea.echap, 'srdb_GrignonTremie2011.csv') 
        #header_row_xydb = ['plant','rank','ranktop','HS','inerv','x','y']
        dfxy = pandas.read_csv(data_file_xydb, names=header_row_xydb, sep=',', index_col=0, skiprows=1, decimal='.')
        dfxy = dfxy.reset_index()
        dfxy = dfxy[dfxy['variety']=='Mercia']
        dfxy['rankclass'] = 1
        dfxy['rankclass'][dfxy['ranktop'] <= 4] = 2
       
    # SR
    header_row_srdb = ['rankclass','s','r']
    dfsr = pandas.read_csv(data_file_srdb, names=header_row_srdb, sep=',', index_col=0, skiprows=1, decimal='.')
    dfsr = dfsr.reset_index()
    
    # NIV1 
    # cle = rankclass
    
    # test seulement rankclass 2 = (feuille du haut = feuille du bas)
    #dfxy['rankclass'] = 2
    # si tout =
    #dfxy['rankclass'] = 1
    #dfxy['rankclass'][dfxy['ranktop'] <= 4] = 2
        
    # creation de la colonne age
    dfxy['age'] = dfxy['HS'] - dfxy['rank'] + 1
    #cut intervalle de 0 a 1, etc.
    dfxy_cut = pandas.cut(dfxy.age, [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    dfxy['age'] = dfxy_cut
    # remplacement intervalle par un int (ex : (x, y)=y seulement)
    # si marche pas : bidouille = 
    dfxy['age'] = dfxy['age'].replace("(-1, 0]", 0)
    dfxy['age'] = dfxy['age'].replace("(0, 1]", 1)
    dfxy['age'] = dfxy['age'].replace("(1, 2]", 2)
    dfxy['age'] = dfxy['age'].replace("(2, 3]", 3)
    dfxy['age'] = dfxy['age'].replace("(3, 4]", 4)
    dfxy['age'] = dfxy['age'].replace("(4, 5]", 5)
    dfxy['age'] = dfxy['age'].replace("(5, 6]", 6)
    dfxy['age'] = dfxy['age'].replace("(6, 7]", 7)
    dfxy['age'] = dfxy['age'].replace("(7, 8]", 8)
    dfxy['age'] = dfxy['age'].replace("(8, 9]", 9)
    dfxy['age'] = dfxy['age'].replace("(9, 10]", 10)
    dfxy['age'] = dfxy['age'].replace("(10, 11]", 11)
    dfxy['age'] = dfxy['age'].replace("(11, 12]", 12)
    '''
    # methode plus propre - ne marche pas
    dfxy = dfxy.reset_index()
    da = 0
    while da < len(dfxy):
        dfxy['age'][da] = re.sub("\(-*[0-9]{1,2}, ([0-9]{1,2})]", "\\1", dfxy['age'][da])
        da = da + 1'''

    # NIV2
    # creation du premier dico par rankclass
    groups = dfxy.groupby(["rankclass"])
    dxy={n:{} for n,g in groups}
    # creation du second dico par age
    for n,g in groups:
        gg = g.groupby('age')
        dxy[n] = {k:[] for k,ggg in gg}
    # remplissage x et y
    groups = dfxy.groupby(["inerv"]) 
    
    #dfxy['age'] = dfxy['age'].astype(int)
    for n,d in groups:
        dxy[int(d[['rankclass']].values[0])][int(d[['age']].values[0])].append(d.ix[:,['x','y']].to_dict('list'))
    
    # NIV3 : ajout des s et r
    # traitement de dfsr (creation de 2 listes) DE STRING
    dr = 0
    s1 = []; r1 = []
    s2 = []; r2 = []
    while dr<len(dfsr):
        if dfsr['rankclass'][dr] == 1:
            s1.append(dfsr['s'][dr])
            r1.append(dfsr['r'][dr])
            dr = dr + 1
        else :
            s2.append(dfsr['s'][dr])
            r2.append(dfsr['r'][dr])
            dr = dr + 1
            
    # ajout dans le dict de dict precedent
    # remplissage rankclass (rankclass=2 seulement pour rht3)
    if name is 'Rht3':
        rank2keys = dxy[2].keys()
        cpt2 = 0; rank2 = rank2keys[cpt2]; list2 = 0  
        while cpt2 < len(rank2keys):
            rank2 = rank2keys[cpt2]
            while list2 < len(dxy[2][rank2]) :
                dxy[2][rank2][list2].update({'s':s2, 'r':r2})
                list2 = list2 + 1
            if list2 == len(dxy[2][rank2]) :
                list2 = 0
            rank2 = rank2 + 1
            cpt2 = cpt2 + 1
    else:
        rank1keys = dxy[1].keys()
        cpt1 = 0; rank1 = rank1keys[cpt1]; list1 = 0
        while cpt1 < len(rank1keys):
            rank1 = rank1keys[cpt1]
            while list1 < len(dxy[1][rank1]) :
                dxy[1][rank1][list1].update({'s':s1, 'r':r1})
                list1 = list1 + 1
            if list1 == len(dxy[1][rank1]) :
                list1 = 0
            rank1 = rank1 + 1
            cpt1 = cpt1 + 1
        rank2keys = dxy[2].keys()
        cpt2 = 0; rank2 = rank2keys[cpt2]; list2 = 0  
        while cpt2 < len(rank2keys):
            rank2 = rank2keys[cpt2]
            while list2 < len(dxy[2][rank2]) :
                dxy[2][rank2][list2].update({'s':s2, 'r':r2})
                list2 = list2 + 1
            if list2 == len(dxy[2][rank2]) :
                list2 = 0
            rank2 = rank2 + 1
            cpt2 = cpt2 + 1
        
    
    return dxy
    
#-------------------------------------------------------------------------------  
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
    
    
 # deprecated   
    
def diff(edata, s):
    e1 = edata.head(0)
    e2 = s.head(0) 
    e = list(set(e1) - set(e2))
    i = 0
    while i<(len(e)-1):
        s [e[i]] = None
        i = i + 1
    return s
    
def emis(data, d, n):
    res = d.ix[:,'Nff':] / n.ix[:,'Nff':] 
    for c in d.ix[:,'Nff':].columns.values:
        d.ix[:,c] = res.ix[:,c]
    # update TPE and compute MB and TT
    d['TPE'] = d.ix[:,'TC':'T5'].sum(axis=1)
    edata = data[data['Date'] >= 6]
    edata = edata.reset_index()
    grouped = edata.groupby(['Var','Date'],as_index=False)
    return grouped
    
#-------------------------------------------------------------------------------  
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
    - when information is available at date 1 or 3, question marks at date 2 were replaced by confirmed tiller positions """
    
    fn = shared_data(alinea.echap, 'Tillering_data_Mercia_Rht3_2010_2011.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    date_code = {'d1':'2011-01-19','d2':'2011-04-18','d3':'2011-04-27', 'd4':'2011-05-26','d5':'2011-06-01','d6':'2011-06-09'}
    
    # infer emision at date 3 from Total primary present (TP)
    #add statistics for T6 as well as it may be present at date 3    
    def infer_d3(g):
        g.insert(13,'T6',0.0)
        g['T6'][g['MB']!=1] = numpy.nan # avoid infering T6 = 0 on dead plants
        TP = g['TP']
        date = g['Date']
        if all(TP[date==3].notnull()):
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
           'ears_per_plant': ears_per_plant[k],
           'plant_survival': viability[k],
           'tillers_per_plant': axdyn[k], 
           'date_code' : date_code} for k in ('Mercia', 'Rht3')}
    return obs

def Tillering_data_Tremie1_2011_2012():
    """Tillering data for Boigneville 2011-2012
    
    Data come from different plants than those tagged at date 1 and date 6, from tagged plants at date 3 for presence / absence of T7 and at date 7 for counts of fertile tillers.
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
    - Emission of tiller 5 and 6 were ever noted.
    - at date 7, values of fertile tiller  per plants seems buggy compared to other dates
    """
    
    fn = shared_data(alinea.echap, 'Tillering_data_Tremie1_2011_2012.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    date_code = {'d1':'2012-03-09','d2':'2012-04-02','d3':'2012-04-11', 'd4':'2012-05-09','d5':'2012-05-29','d6':'2012-06-12', 'd7':'2012-07-12', 'd8':'2012-04-04'}
    
    # compute emmission probability from data at date 1, and adding data at date 3 for tiller 7 only
    edata = data[data['Date'] == 1]
    edata['T7'] = data['T7'][data['Date'] == 3].values
    edata = edata.reset_index()
    emission = emission_probabilities(edata,last='T7')

    # tiller dynamics
    axdyn = axis_dynamics(data)
    
    # compute ear_per_plant using data at date 6 and 7 and plot data of fertile tillers at date 4
    eardata = pandas.DataFrame(axdyn).ix[:,('Date', 'FT')].dropna()
    pdata = Plot_data_Tremie_2011_2012()
    ftd4 = 1.*numpy.array(pdata['fertile_axis_density']['2012-05-09'])  / numpy.array(pdata['plant_density']['2012-05-09']) - 1 #remove main stem to get FT count 
    # ear per plante at date 7 very strange (twice the value at other dates and non-existence of unfertile tillers): ears_per_plant taken as mean of counting at date 4 and 6
    ears_per_plant = 1 + numpy.array(ftd4.tolist() + eardata['FT'][eardata['Date'] == 6].tolist()).mean()

    
    obs = {'emission_probabilities': emission,
           'ears_per_plant': ears_per_plant,
           'tillers_per_plant': axdyn, 
           'date_code' : date_code}
    return obs
    
def Tillering_data_Tremie_2012_2013():
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
    
    fn = shared_data(alinea.echap, 'Tillering_data_Tremie_2012_2013.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    date_code = {'d1':'2013-02-13','d2':'2013-03-29', 'd3': '2012-04-19'}
    
    # compute emmission probas of primary axis/ emission per plant of secondary axis 
    edata = data[data['Date'] < 6]
    edata = edata.reset_index()
    grouped = edata.groupby(['Var', 'N'],as_index=False)
    d = grouped.agg(_maxna)
    d = d.reset_index()
    #---
    de = d.fillna(0)
    grouped1 = de.groupby(['Var','Date'],as_index=False)
    #---
    grouped = d.groupby(['Var','Date'],as_index=False)
    s = grouped.agg('sum')
    n = grouped1.agg(lambda x: float(x.dropna().count())) 
    diff(edata, s)
    d = s
    emis(data, d, n)
    s = grouped.agg('sum')
    n = grouped.agg(lambda x: x.dropna().count())
    d['MB'] = 1.0 * s['MB'] / n['MB']
    d['TT'] = 1.0 * s['TT'] / n['TT']
    obs = {'plant_density_at_emergence' : {'Tremie': 204}, 'ear_density_at_harvest' : {'Tremie': 676}, 
            'tillering' : d}
    return obs
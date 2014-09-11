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
    - estimates of plant density at emergence (id Archi12)
    - estimates of ear density at harvest (all axes, id Archi10)
    
    Notes :
    - Missing raw data for plant counts at emergence (Archi7, used Archi12 instead)
    """
    d={'Mercia': {
        'sowing_density': 220, 
        'plant_density_at_emergence' : 203, 
        'ear_density_at_harvest' : 444,
        'inter_row': 0.15},
        'Rht3': {
        'sowing_density': 220, 
        'plant_density_at_emergence' : 211, 
        'ear_density_at_harvest' : 384,
        'inter_row': 0.15}}
    return d

def Plot_data_Tremie_2011_2012():
    d = {}
    d['Tremie'] = {'plant_density_at_emergence' : 279, 'ear_density_at_harvest' : 491,
    'inter_row':0.143
    }
    return d
    
def Plot_data_Tremie_2012_2013():
    d = {}
    d['Tremie2'] = {'plant_density_at_emergence' : 204, 'ear_density_at_harvest' : 676,
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
    em = em.reset_index()
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
    s = grouped.agg('sum').ix[:,'TP':'TT']
    n = grouped.agg(lambda x: x.apply(lambda x: x.dropna().count())).ix[:,'TP':'TT']
    axis =  s / n
    # add main stem to primary tillers and total tillers
    axis.ix[:,('TP','TT')] = axis.ix[:,('TP','TT')] + 1
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
    
    Data found in Archi11,Archi19, Archi33 were used to build a synthetic table containing :
        - presence/absence of living mainstem (column MS, NA means the plant is dead) and number of leaf emited (Nff)
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive)
        - total number of primary tillers (column TP) and/or secondary tillers (dead or alive, column TS), and/or total number of tiller present (primary or secondary, dead or alive). 
        - number of fertile tillers (column FT) at the end from the counting of tillers that are alive or that bears an ear 
        
    These data are aggregated to estimate primary emission probability, dynamics of mean number of axis present on plants, plant mortality and  number of ears per plant
        
    Others data available in raw files, but not included yet:
    - Number of secondary tiller with HS >=2 (Nt3f)  at date 3 (Archi33)
    
    Notes :
    - No data at date 5
    - when information is available at date 1 or 3, question marks at date 2 were replaced by confirmed tiller positions """
    
    fn = shared_data(alinea.echap, 'Tillering_data_Mercia_Rht3_2010_2011.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
    
    # infer emision at date 3 from Total primary present (TP)
    #add statistics for T6 as well as it may be present at date 3    
    def infer_d3(g):
        g.insert(13,'T6',0.0)
        g['T6'][g['MB']!=1] = numpy.nan # avoid infering T6 = 0 on dead plants
        TP = g['TP']
        date = g['Date']
        if TP[date==3].notnull():
            infer = g.ix[date==2,'TC':'T6'] 
            if TP[date==3] > TP[date==2]:
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
    emission = grouped.apply(emission_probabilities).to_dict()
    
    # compute ear_per_plante (including main stem) 
    eardata = newdata[newdata['Date'] >= 6]
    eardata = eardata.reset_index()
    grouped = eardata.groupby('Var',as_index=False)    
    ears_per_plant = grouped.apply(lambda x:  1  + (x['FT'].sum() / x['FT'].dropna().count())).to_dict()
    
    # compute plant viability
    grouped = data.groupby('Var',as_index=False)
    viability = grouped.apply(plant_viability).to_dict()
    
    #compute axis dynamics (with main stem added)
    axdyn = grouped.apply(axis_dynamics).to_dict()
    
    obs = {'emission_probabilities': emission,
           'ears_per_plant': ears_per_plant,
           'plant_viability': viability,
           'axis_per_plant': axdyn}
    return obs

def Tillering_data_Tremie1_2011_2012():
    """Tillering data for Boigneville 2011-2012
    Expected data are :
    - estimates of plant density at emergence (id Archi1)
    - estimates of ear density at harvest (all axes) (id Archi16)
    - estimates of the number of elongated internodes on main stems
    - estimates of tiller dynamics, with three kind of data : 
        - presence/absence of living mainstem (column MS) and number of leaf emited
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive), and/or  total number of primary tillers (column TPE) and secondary tillers emited (dead or alive, column TSE). These data are used for fitting emission probabilities
        - estimation of total living tillers (column TT) at the end from the counting of tillers that are alive or that bears an ear
    Others data :
    - Number of secondary tiller to date 1 
    Notes :
    - No data at date 2, 4 & 5"""
    
    # summarise tillering data to mean number of tiller emited per plant (columns Tc->T7), mean number of leaves, and total tillers and MS present per plant at the end (column TT and MB)
    # this table was contrusted using Archi2, Archi3, Archi14, Archi15
    fn = shared_data(alinea.echap, 'Tillering_data_Tremie1_2011_2012.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
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
    obs = {'plant_density_at_emergence' : {'Tremie': 278.55}, 'ear_density_at_harvest' : {'Tremie': 491}, 
            'tillering' : d}
    return obs
    
def Tillering_data_Tremie_2012_2013():
    """Tillering data for Boigneville 2012-2013
    Expected data are :
    - estimates of plant density at emergence (id Archi-)
    - estimates of ear density at harvest (all axes) (id Archi-)
    - estimates of the number of elongated internodes on main stems
    - estimates of tiller dynamics, with three kind of data : 
        - presence/absence of living mainstem (column MS) and number of leaf emited
        - presence/absence of primary tiller (column Tc-> Tn) whatever the state (dead or alive), and/or  total number of primary tillers (column TPE) and secondary tillers emited (dead or alive, column TSE). These data are used for fitting emission probabilities
        - estimation of total living tillers (column TT) at the end from the counting of tillers that are alive or that bears an ear
    Notes :
    - No data to date 3 and after this date"""
    
    # summarise tillering data to mean number of tiller emited per plant (columns Tc->T5), mean number of leaves, and total tillers and MS present per plant at the end (column TT and MB)
    fn = shared_data(alinea.echap, 'Tillering_data_Tremie_2012_2013.csv')
    data = pandas.read_csv(fn,decimal=',',sep='\t')
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
""" Reconstruct a MTG wheat canopy at different treatment dates, save dye interception 
    and compare it with field data """

import pandas
import numpy

from openalea.deploy.shared_data import shared_data
import alinea.echap

from alinea.echap.architectural_reconstructions import echap_reconstructions
import alinea.echap.interception_data as idata
from alinea.echap.interception_leaf import InterceptModel, pesticide_applications
from alinea.echap.interfaces import pesticide_interception
from alinea.echap.recorder import LeafElementRecorder
from alinea.echap.interfaces import record as do_record

# Intervalle de confiance
def mean(lst):
    return sum(lst) / float(len(lst))

def variance(lst):
    """
    Uses standard variance formula (sum of each (data point - mean) squared)
    all divided by number of data points
    """
    mu = mean(lst)
    return 1.0/(len(lst)-1) * sum([(i-mu)**2 for i in lst])

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
    
    return math.sqrt(v/n) * c

def df_interception_path(variety = 'Tremie12', nplants = 30):
    filename = 'interception_'+variety.lower() + '_' + str(nplants) + 'pl.csv'
    return str(shared_data(alinea.echap)/'interception_simulations'/filename)

def df_dates_emergence_path(variety='Tremie12', nplants=30):
    filename = 'dates_emergence_'+variety.lower() + '_' + str(nplants) + 'pl.csv'
    return str(shared_data(alinea.echap)/'interception_simulations'/filename)
    
def repartition_at_application(name='Tremie12', appdate='T1', dose=1e4, nplants=30, density=1, dimension=1,
                                reset=False, reset_data=False, numero_sim=1):
    """ 10000 l ha-1 is for 1 l/m2 """
    HS_applications = idata.dye_applications()
    # If adjustment needed : import echap parameters, modify parameters and create new reconstruction method
    # nlim_Mercia=3, nlim_Rht3=2, nlim_Tremie12=4, nlim_Tremie13=3
    #Reconstructions_appli = reconstructions.EchapReconstructions(nlim_factor = {'Mercia':nlim_Mercia, 'Rht3':nlim_Rht3, 'Tremie12':nlim_Tremie12, 'Tremie13':nlim_Tremie13} )
    reconst = echap_reconstructions(reset=reset, reset_data=reset_data)
    conv = reconst.HS_fit[name]
    adel = reconst.get_reconstruction(name=name, nplants=nplants, stand_density_factor = {name:density} , dimension=dimension)
    
    appdate_ref = ['T1-0.4', 'T1-0.2', 'T1', 'T1+0.2', 'T1+0.4',
            'T2-2.5', 'T2-2', 'T2-1.5', 'T2-1', 'T2-0.5', 'T2-0.4', 'T2-0.2', 'T2', 'T2+0.2', 'T2+0.4', 'T2+0.5', 'T2+1', 'T2+1.5', 'T2+2', 'T2+2.5']
    if appdate in appdate_ref:
        date, hs = HS_applications[name][appdate]
    else :
        hs = float(appdate); date = ''
    age = float(conv.TT(hs))
    g = adel.setup_canopy(age)
    
    recorder = LeafElementRecorder()
    applications= 'date,dose, product_name\n%s 10:00:00, %f, Tartrazine'%(date, dose)
    application_data = pesticide_applications(applications)
    interceptor = InterceptModel({'Tartrazine':{'Tartrazine': 1}}) # consider a 1g.l-1 tartrazine solution
    
    
    ########################### ATTENTION : si pas de domain, on est pas en couvert infini!!!!!!!!!!!!!!!!!!!!!!!!
    g,_ = pesticide_interception(g, interceptor, application_data, domain=adel.domain)
    do_record(g, application_data, recorder, header={'TT':age, 'HS':hs})
    df =  recorder.get_records()
    
    #add columns
    df['deposit_Tartrazine'] = df['surfacic_doses_Tartrazine'] * df['area']
    df['var'] = name
    df['treatment'] = appdate
    df['nb_plantes_sim'] = nplants
    df['numero_sim'] = numero_sim
        
    # compute ntop_cur :************************************BUG ADEL: some leaves have age < 0: filter them before ntop cur estimation ?
    gr = df.groupby(['HS','plant','axe'], group_keys=False)
    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        return sub
    df = gr.apply(_fun)
    df['ntop_cur'] = df['n_max'] - df['metamer'] + 1
    
    # compute leaf emergence
    df['leaf_emergence'] = df['TT'] - df['age'] 
    
    return adel, g, df
 
def dates_emergence(adel):
    df_phenT = adel.phenT()
    nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
    df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
    df_phenT['ntop_cur'] = df_phenT['nff'] - df_phenT['n']
    dates = [numpy.mean(df_phenT[df_phenT['ntop_cur']==lf]['tip']) for lf in range(1, int(max(df_phenT['n'])+2))]
    df_dates = pandas.DataFrame(dates, index = range(1, len(dates)+1), columns = ['leaf_emergence'])
    df_dates.index.name = 'ntop_cur'
    df_dates = df_dates.reset_index()
    return df_dates
 
def aggregate_by_leaf(df):
    """
    Aggregate interceptioin data by leaf and add colmun 'deposits' (= area * surfacic doses)
    """
    def first_val(x):
        return x[x.first_valid_index()]
        

    df = df[(df['hasEar'] == 1)]# do not consider aborting axes
    
    df = df.convert_objects()
    gr = df.groupby(['HS', 'plant','axe','metamer'], as_index=False)
    df_agg = gr.agg({'var': first_val,
                   'treatment': first_val,
                    'nb_plantes_sim':first_val,
                    'numero_sim':first_val,
                    'age':first_val,
                    'leaf_emergence': first_val,
                    'TT': first_val, 
                   'area': numpy.sum,
                   'date': first_val,
                   'green_area': numpy.sum,
                   'senesced_area': numpy.sum,
                   'id':lambda(x): '_'.join(map(str,x)),
                   'length':numpy.sum,
                   'ntop':first_val,
                   'ntop_cur':first_val,
                   'organ':first_val,
                   'penetrated_doses_Tartrazine': numpy.mean,
                   'surfacic_doses_Tartrazine':numpy.mean,
                   'deposit_Tartrazine': numpy.sum,
                   'lifetime': first_val,
                   'exposition': first_val,
                   })
    # strange transforms ?
    #df = df.rename(columns = {'TT':'degree_days', 'plant':'num_plant', 
    #                              'nff':'fnl', 'axe':'axis', 'ntop':'num_leaf_top'})
    #df['date'] = df['date'].apply(lambda x: pandas.to_datetime(x, dayfirst=True))
    #df['num_leaf_bottom'] = df['fnl'] - df['num_leaf_top'] + 1                   
    return df_agg

def interception_statistics(df_leaf, axis='MS'):
    data = df_leaf
    
    if axis == 'MS':
        data = data[data['axe'] == 'MS']
        
    sub = data.ix[:,('var','treatment','HS','TT','axe','metamer','ntop','ntop_cur','age','length','area','green_area','deposit_Tartrazine','exposition','lifetime')]
    dfmoy = sub.groupby(['var','treatment','HS','axe','ntop_cur']).mean().reset_index()
    dfsd = sub.groupby(['var','treatment','HS','axe','ntop_cur']).std().reset_index()   

    return dfmoy, dfsd
            

 
def run_dye_interception(variety = 'Tremie12', nplants = 30,
                                         treatments = {'T1':1e4, 'T2':1e4},
                                         nrep=1,
                                         density=1, dimension=1):
    
    dfint = []  
    for i in range(nrep):
        for date, dose in treatments.iteritems():
            adel, g, df = repartition_at_application(name=variety, appdate=date, 
                                                     dose=dose, nplants=nplants, 
                                                     density=density, dimension=dimension, numero_sim=i + 1)
            df_i = aggregate_by_leaf(df)
            dfint.append(df_i)

            
    df_interception = pandas.concat(dfint).reset_index(drop = True)
    df_interception.to_csv(df_interception_path(variety = variety, nplants = nplants), index = False)
    
    df_meanMS, df_sdMS = interception_statistics(df_interception)
    return df_meanMS, df_sdMS
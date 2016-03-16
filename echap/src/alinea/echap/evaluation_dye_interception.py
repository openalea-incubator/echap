""" Reconstruct a MTG wheat canopy at different treatment dates, save dye interception 
    and compare it with field data """

import pandas
import numpy
import os
from copy import deepcopy
from scipy.stats import t

from openalea.deploy.shared_data import shared_data
import alinea.echap

from alinea.echap.architectural_reconstructions import EchapReconstructions, echap_reconstructions, reconstruction_parameters
import alinea.echap.interception_data as idata
from alinea.echap.interception_leaf import InterceptModel, pesticide_applications
from alinea.echap.interfaces import pesticide_interception
from alinea.echap.recorder import LeafElementRecorder
from alinea.echap.interfaces import record as do_record

# Intervalle de confiance


def conf_int(lst, perc_conf=95):
    """
    Confidence interval - given a list of values compute the square root of
    the variance of the list (v) divided by the number of entries (n)
    multiplied by a constant factor of (c). This means that I can
    be confident of a result +/- this amount from the mean.
    The constant factor can be looked up from a table, for 95pcent confidence
    on a reasonable size sample (>=500) 1.96 is used.
    """  
    n, v = len(lst), numpy.var(lst, ddof = 1)
    c = t.interval(perc_conf * 1.0 / 100, n-1)[1]
    
    return numpy.sqrt(v/n) * c

def df_interception_path(variety = 'Tremie12', nplants = 30, simulation = 'reference'):
    filename = 'interception_'+ simulation + '_' + variety.lower() + '_' + str(nplants) + 'pl.csv'
    return str(shared_data(alinea.echap)/'interception_simulations'/filename)

def get_reconstruction(variety='Tremie13', nplants=30, reconstruction_pars=None, reset=False, reset_data=False):
    if reconstruction_pars is None: 
        reconst = echap_reconstructions(reset=reset, reset_data=reset_data)
    else:
        pars = reconstruction_parameters()
        pars.update(reconstruction_pars)
        reconst = EchapReconstructions(pars = pars, reset_data=reset_data)
    hsfit = reconst.HS_fit[variety]
    adel = reconst.get_reconstruction(name=variety, nplants=nplants)
    return adel, hsfit
    
def canopy_age(variety='Tremie13', date='T1'):
    TT = idata.TT_index()
    return round(TT.ix[(variety,date),0], 2)
    
def repartition_at_application(adel, hsfit, variety, date='T1', dose=1e4, num_sim=1):
    """ 10000 l ha-1 is for 1 l/m2 """
    tags = idata.tag_dates()
    if date in tags[variety]:
        tag = date
        date = tags[variety][date]
    else:
        tag=''
        
    recorder = LeafElementRecorder()
    applications= 'date,dose, product_name\n%s 10:00:00, %f, Tartrazine'%(date, dose)
    application_data = pesticide_applications(applications)
    interceptor = InterceptModel({'Tartrazine':{'Tartrazine': 1}}) # consider a 1g.l-1 tartrazine solution
    
    age = canopy_age(variety, date)
    g = adel.setup_canopy(age)   
    ########################### ATTENTION : si pas de domain, on est pas en couvert infini!!!!!!!!!!!!!!!!!!!!!!!!
    g,_ = pesticide_interception(g, interceptor, application_data, domain=adel.domain)
    do_record(g, application_data, recorder, header={'TT':age, 'HS':float(hsfit(age))})
    df =  recorder.get_records()
    #add midribs data
    midribs = adel.midrib_statistics(g)
    #add columns
    df['deposit_Tartrazine'] = df['surfacic_doses_Tartrazine'] * df['area']
    df['variety'] = variety
    df['treatment'] = tag
    df['nb_plantes_sim'] = adel.nplants
    df['domain_area'] = adel.domain_area * 10000 #m2 -> cm2
    df['numero_sim'] = num_sim
    # compute n_max, ntop_cur
    gr = df.groupby(['HS','plant','axe'], group_keys=False)
    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        sub['nflig'] = sub['metamer'][sub['is_ligulated'] > 0].max()
        return sub
    df = gr.apply(_fun)
    df['ntop_cur'] = df['n_max'] - df['metamer'] + 1
    df['ntop_lig'] = df['nflig'] - df['metamer'] + 1
    # compute leaf emergence
    df['leaf_emergence'] = df['TT'] - df['age'] 
    
    return g, df, midribs
 

def aggregation_types(what=None):
    """ aggregating functions for one simulation
    """
    def first_val(x):
        return x[x.first_valid_index()]
        
    types = {'variety': first_val,
                   'treatment': first_val,
                    'nb_plantes_sim':first_val,
                    'domain_area':first_val,
                    'numero_sim':first_val,
                    'age':first_val,
                    'leaf_emergence': first_val,
                    'TT': first_val, 
                    'nff':first_val,
                   'area': numpy.sum,
                   'date': first_val,
                   'green_area': numpy.sum,
                   'senesced_area': numpy.sum,
                   'id':lambda(x): '_'.join(map(str,x)),
                   'length':numpy.sum,
                   'ntop':first_val,
                   'ntop_cur':first_val,
                   'ntop_lig':first_val,
                   'organ':first_val,
                   'mature_length': first_val,
                   'penetrated_doses_Tartrazine': numpy.mean,
                   'surfacic_doses_Tartrazine':numpy.mean,
                   'deposit_Tartrazine': numpy.sum,
                   'lifetime': first_val,
                   'exposition': first_val
                   }
    if what is None:
        return types
    else:
        return {w:types[w] for w in what}
    
 
def aggregate_by_leaf(df):
    """
    Aggregate interceptioin data by leaf and add colmun 'deposits' (= area * surfacic doses)
    """

    df = df[(df['hasEar'] == 1)]# do not consider aborting axes
    
    df = df.convert_objects()
    gr = df.groupby(['HS', 'plant','axe','metamer'], as_index=False)
    df_leaf = gr.agg(aggregation_types())
    # strange transforms ?
    #df = df.rename(columns = {'TT':'degree_days', 'plant':'num_plant', 
    #                              'nff':'fnl', 'axe':'axis', 'ntop':'num_leaf_top'})
    #df['date'] = df['date'].apply(lambda x: pandas.to_datetime(x, dayfirst=True))
    #df['num_leaf_bottom'] = df['fnl'] - df['num_leaf_top'] + 1                   
    return df_leaf

def aggregate_by_axe(df_leaf):
    """
    Aggregate  leaf-aggregated interceptioin data by axe
    """

    res=[]
    what = ['variety','treatment', 'nb_plantes_sim','TT','nff',
                   'area', 'date', 'green_area','senesced_area', 'organ',
                   'surfacic_doses_Tartrazine','deposit_Tartrazine']
                   
    for name,gr in df_leaf.groupby(['HS', 'plant','axe','numero_sim'], as_index=False):
        df_agg = gr.groupby(['HS', 'plant','axe','numero_sim'], as_index=False).agg(aggregation_types(what))
        gr = gr.sort('metamer')
        frac = gr['length'] / gr['mature_length']
        ilig  = numpy.max(numpy.where(frac >= 1))
        df_agg['haun_stage'] = gr['metamer'].values[ilig] + frac[(ilig+1):].sum()
        df_agg['leaves'] = '_'.join(map(str, gr['metamer'].drop_duplicates()))
        res.append(df_agg)
    return pandas.concat(res)

def aggregate_by_plant(df_leaf):
    """
    Aggregate  leaf-aggregated interceptioin data by axe
    """

    what = ['variety','treatment', 'nb_plantes_sim','TT','nff',
                   'area', 'date', 'green_area','senesced_area', 'organ', 'domain_area',
                   'surfacic_doses_Tartrazine','deposit_Tartrazine']
    agg = df_leaf.groupby(['HS', 'plant','numero_sim'], as_index=False).agg(aggregation_types(what))
    plant_domain = agg['domain_area'] / agg['nb_plantes_sim']
    agg['fraction_intercepted'] = agg['deposit_Tartrazine'] / plant_domain
    agg['lai'] = agg['area'] / plant_domain
    agg['glai'] = agg['green_area'] / plant_domain
    return agg
    
def axis_statistics(df_sim, what='deposit_Tartrazine', err=conf_int, axis='MS'):
    data = aggregate_by_axe(df_sim)
    
    if axis == 'MS':
        data = data[data['axe'] == 'MS']
     
    sub = data.ix[:,['treatment',what]]     
    agg = sub.groupby('treatment').mean().reset_index()
    errag = sub.groupby('treatment').agg(err).reset_index()
    agg['ybar'] = agg[what]
    agg['yerr'] = errag[what]
    agg=agg.set_index('treatment')
    return agg

def plant_statistics(df_sim, what='deposit_Tartrazine', err=conf_int):
    data = aggregate_by_plant(df_sim)
         
    sub = data.ix[:,['treatment',what]]     
    agg = sub.groupby('treatment').mean().reset_index()
    errag = sub.groupby('treatment').agg(err).reset_index()
    agg['ybar'] = agg[what]
    agg['yerr'] = errag[what]
    agg=agg.set_index('treatment')
    return agg
    
def leaf_statistics(df_sim, what='deposit_Tartrazine', err=conf_int, by='ntop_cur', axis='MS'):

    data = df_sim
    if axis == 'MS':
        data = data[data['axe'] == 'MS']
    if not isinstance(what, list):
        what = [what]
    sub = data.ix[:,['treatment', by] + what]
    agg = sub.groupby(['treatment',by]).mean().reset_index()
    errag = sub.groupby(['treatment',by]).agg(err).reset_index()
    agg['ybar'] = agg[what[0]]
    agg['yerr'] = errag[what[0]]
    if by == 'ntop_cur':    
        agg['xbar'] = ['F' + str(int(i)) for i in agg[by]]
    elif by == 'ntop_lig':
        agg['xbar'] = ['Fl' + str(int(i)) for i in agg[by]]
    else:
        agg['xbar'] = ['L' + str(int(leaf)) for leaf in agg['metamer']]
    agg=agg.set_index('treatment')
    return agg
    

def pdict(value):
    """ create a parameter dict for all echap cultivar with value
    """
    return {k:deepcopy(value) for k in ('Mercia', 'Rht3','Tremie12', 'Tremie13')}    

def simulation_tags():
    tags = {'reference': {'dose' : pdict(1e4),
                          'reconstruction_pars' : None},
           }
    for shape in ('MerciaRht','Tremie','Soissons','Tremie12','Tremie13'):
        tags['shape_' + shape + '_byleafclass'] = {'dose': pdict(1e4),
                                  'reconstruction_pars': {'xy_data': pdict(shape + '_byleafclass'),
                                                          'top_leaves':pdict(4)
                                                          }
                                 }
        tags['shape_' + shape] = {'dose': pdict(1e4),
                                  'reconstruction_pars': {'xy_data': pdict(shape),
                                                          'top_leaves':pdict(0)
                                                          }
                                 }                                 
    return tags

def run_sim(adel, hsfit, var, date, dose, num_sim):
    g, df, midribs = repartition_at_application(adel, hsfit, var, date=date, dose=dose, num_sim=num_sim)
    df_i = aggregate_by_leaf(df)
    midribs = midribs.rename(columns={'leaf':'metamer'})
    midribs['plant'] = ['plant'+str(p) for p in midribs['plant']]
    # add domain area
    return df_i.merge(midribs)

    
def dye_interception(variety = 'Tremie12', nplants = 30, nrep = 1, simulation = 'reference', treatments=None, reset=False, reset_data=False):
                         
    path = df_interception_path(variety = variety, nplants = nplants, simulation = simulation)                      
    repetitions = range(1, nrep + 1)
    if treatments is None:
        treatments = idata.tag_treatments()[variety]['application']
    sim = simulation_tags()[simulation]
    new_sim = False
    dfint = []
    
    if os.path.exists(path):
        df_old = pandas.read_csv(path)
        done_reps = list(set(df_old['numero_sim']))
        missing = [rep for rep in repetitions if not rep in done_reps]
        #check all treatments are there
        tocheck = df_old.groupby('numero_sim')
        for i_sim, df in tocheck:
            done = df['treatment'].drop_duplicates().values
            dfint.append(df)
            to_do = [t for t in treatments if not t in done]
            if len(to_do) > 0:
                new_sim = True
                adel, hsfit = get_reconstruction(variety=variety, nplants=nplants, reconstruction_pars=sim['reconstruction_pars'], reset=reset, reset_data=reset_data)
                for t in to_do:
                     df_t = run_sim(adel, hsfit, variety, date=t, dose=sim['dose'][variety], num_sim=i_sim)
                     dfint.append(df_t)
    else:
        missing = repetitions

    if len(missing) > 0:
        new_sim = True
        for i in missing:
            adel, hsfit = get_reconstruction(variety=variety, nplants=nplants, reconstruction_pars=sim['reconstruction_pars'], reset=reset, reset_data=reset_data)
            for t in treatments:
                df_t = run_sim(adel, hsfit, variety, date=t, dose=sim['dose'][variety], num_sim=i)
                dfint.append(df_t)       
                
    df_interception = pandas.concat(dfint).reset_index(drop = True)
    if new_sim:
        df_interception.to_csv(path, index = False)
    
    df_interception = df_interception[df_interception['numero_sim'].isin(repetitions)]
    df_interception = df_interception[df_interception['treatment'].isin(treatments)]
    return df_interception
""" Reconstruct a MTG wheat canopy at different treatment dates, save dye interception 
    and compare it with field data """

import pandas
from openalea.deploy.shared_data import shared_data
from alinea.echap.architectural_reconstructions import echap_reconstructions
import alinea.echap.interception_data as idata
from alinea.echap.interception_leaf import InterceptModel, pesticide_applications
from alinea.echap.interfaces import pesticide_interception
from alinea.echap.recorder import LeafElementRecorder
from alinea.echap.interfaces import record as do_record

def get_file_name(variety = 'Tremie12', nplants = 30,):
    return 'interception_'+variety.lower() + '_' + str(nplants) + 'pl.csv'

def get_file_path(variety = 'Tremie12', nplants = 30):
    filename = get_file_name(variety = variety, nplants = nplants)
    return str(shared_data(alinea.echap)/'interception_simulations'/filename)    

def repartition_at_application(name, appdate='T1', dose=1e4, nplants=30, density=1, dimension=1,
                                reset=False, reset_data=False):
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
    age = conv.TT(hs)
    g = adel.setup_canopy(age)
    
    recorder = LeafElementRecorder()
    applications= 'date,dose, product_name\n%s 10:00:00, %f, Tartrazine'%(date, dose)
    application_data = pesticide_applications(applications)
    interceptor = InterceptModel({'Tartrazine':{'Tartrazine': 1}}) # consider a 1g.l-1 tartrazine solution
    
    
    ########################### ATTENTION : si pas de domain, on est pas en couvert infini!!!!!!!!!!!!!!!!!!!!!!!!
    g,_ = pesticide_interception(g, interceptor, application_data)
    do_record(g, application_data, recorder, header={'TT':age, 'HS':hs})
    df =  recorder.get_records()
    print 'repartition_at_application df.columns before ', df.columns
    return adel, g, df
    
def save_dye_interception_single_variety(variety = 'Tremie12', nplants = 30,
                                         treatments = {'T1':1e4, 'T2':1e4},
                                         density=1, dimension=1):
    dfs = []
    for date, dose in treatments.iteritems():
        adel = None
        while adel is None:
            try:
                adel, g, df = repartition_at_application(name = variety, appdate = date, 
                                                         dose = dose, nplants = nplants, 
                                                         density = density, dimension = dimension)
            except:
                pass
        df = df.rename(columns = {'TT':'degree_days', 'plant':'num_plant', 
                                  'nff':'fnl', 'axe':'axis', 'ntop':'num_leaf_top'})
        df['date'] = df['date'].apply(lambda x: pandas.to_datetime(x, dayfirst=True))
        df['num_leaf_bottom'] = df['fnl'] - df['num_leaf_top'] + 1
    df = pandas.concat(dfs)
    df = df.reset_index(drop = True)
    df.to_csv(get_file_path(variety = variety, nplants = nplants), index = False)
    
def save_dye_interception(nplants = 30,treatments = {'T1':1e4, 'T2':1e4}, density=1, dimension=1):
    for variety in ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']:
        save_dye_interception_single_variety(variety = variety, nplants = nplants,
                                             treatments = treatments, density = density,
                                             dimension = dimension)
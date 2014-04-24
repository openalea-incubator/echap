""" Reconstruction of Wheat for Boigneville data

Current reconstructions use fit from Mariem (pgen data) + new modification consigned here
"""
import alinea.echap.architectural_data as archidb
from alinea.adel.stand.stand import agronomicplot
from alinea.adel.astk_interface import AdelWheat

#generic function to be moved to  adel

def plantgen_to_devT(pgen):
    """ Creates devT tables from plantgen dict
    """
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data,plantgen2adel
    from alinea.adel.AdelR import devCsv
    
    axeT_, dimT_, phenT_, _, _, _, _, _, _, _, _ = gen_adel_input_data(**pgen)
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
    devT = devCsv(axeT, dimT, phenT)
    return devT

def plot_dimension(nplants, density, inter_row):
    """ return plot width and length to match nplants created by agronomic_plot
    """
    from math import sqrt
    
    inter_plant = 1. / inter_row / density
    nrow = max(1, int(sqrt(nplants) * 1. * inter_plant / inter_row))
    plant_per_row = max(1, int(nplants * 1. / nrow)) 
    plot_length = inter_plant * plant_per_row
    plot_width = inter_row * nrow
    return plot_length, plot_width
   

def setAdel(nplants, nsect, devT, sowing_density, plant_density, inter_row, seed, sample):    
    
    
    length, width = plot_dimension(nplants, sowing_density, inter_row)
    nplants, positions, domain, domain_area, convUnit = agronomicplot(length=length, 
                                                            width=width, 
                                                            sowing_density=sowing_density, 
                                                            plant_density=plant_density,
                                                            inter_row=inter_row)
                                                            
    adel = AdelWheat(nplants=nplants, positions = positions, nsect=nsect, devT=devT, seed= seed, sample=sample)
    return adel, domain, domain_area, convUnit, nplants
    
    
#---------- reconstructions    

reconst_db={}

def Mercia_2010(nplants=30, nsect=3, seed=1, sample='random', as_pgen=False, dTT_stop=0):

    pgen = archidb.Mercia_2010_plantgen()
    tdata = archidb.Tillering_data_Mercia_Rht3_2010_2011()['tillering']
    tdata = tdata[tdata.Var=='Mercia']
    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Mercia']
    #adapt reconstruction
    # TO DO get nff probailitiesfor main stem from tillerng data ?
    if not as_pgen:
        #handle fly damages to MB  by reducing proba of appearance of T2 (originally equal to 1) and simplify TC
        # BUG : MB is equal to 1 => should return  less than 1
        primary_proba={k:tdata[k].values for k in ('T1','T3','T4')}
        primary_proba['T2']= tdata['MB'].values + tdata['TC'].values
        ears_per_plant = tdata['MB'].values + tdata['TT'].values
        plant_density_at_harvest = float(pdata['ear_density_at_harvest']) / ears_per_plant
        nff = float(tdata['Nff'].values)
        # pgen adapt
        pgen['decide_child_axis_probabilities'] = primary_proba
        pgen['plants_density'] = plant_density_at_harvest #Mariem used plant_density_at_emergence
        pgen['plants_number'] = pgen['plants_density']#to be optimised from a precision wanted on tillering probas and used to compute plot dimension
        pgen['ears_density'] = pdata['ear_density_at_harvest']
        #hack for removing tillers
        pgen['delais_TT_stop_del_axis'] -= dTT_stop
        
    #generate reconstruction
    devT = plantgen_to_devT(pgen)
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, 
                                                           devT=devT, 
                                                           sowing_density=pdata['plant_density_at_emergence'], plant_density=plant_density_at_harvest, 
                                                           inter_row=pdata['inter_row'], 
                                                           seed=seed, sample=sample)
    
    return pgen, adel, domain, domain_area, convUnit, nplants
 
reconst_db['Mercia'] = Mercia_2010
 
age_at_application = {'Mercia_2010': {'2011-04-19': 1166,
                                      '2011-05-11' : 1500}}
                                
                                   
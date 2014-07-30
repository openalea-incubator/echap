""" Reconstruction of Wheat for Boigneville data

Current reconstructions use fit from Mariem (pgen data) + new modification consigned here
"""
import alinea.echap.architectural_data as archidb
from alinea.adel.stand.stand import agronomicplot
from alinea.adel.astk_interface import AdelWheat
from copy import copy
import pandas as pd

#generic function to be moved to adel

def plantgen_to_devT(pgen):
    """ Creates devT tables from plantgen dict
    """
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data, plantgen2adel
    from alinea.adel.AdelR import devCsv
    
    axeT_, dimT_, phenT_, _, _, _, _, _, _, _, _ = gen_adel_input_data(**pgen)
    
    # verif des sorties (vu avec Camille)
    #axeT_.to_csv('C:/Users/Administrateur/openaleapkg/echap/test/axeT_test.csv')
    #dimT_.to_csv('C:/Users/Administrateur/openaleapkg/echap/test/dimT_test.csv')
    #phenT_.to_csv('C:/Users/Administrateur/openaleapkg/echap/test/phenT_test.csv')
    
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
    devT = devCsv(axeT, dimT, phenT)
    return devT
    
def plantgen_to_devT_comp(pgen):
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data,plantgen2adel
    
    axeT_, dimT_, phenT_, _, _, _, _, _, _, _, _ = gen_adel_input_data(**pgen)
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
    return axeT, dimT, phenT

def plot_dimension(nplants, sowing_density, plant_density, inter_row):
    """ return plot width and length to match nplants created by agronomic_plot
    """
    from math import sqrt
    
    side = sqrt(1. / sowing_density * nplants)
    inter_plant = 1. / inter_row / sowing_density
    nrow = max(1, round(side / inter_row))
    plant_per_row = max(1, round(side / inter_plant)) 
    plot_length = inter_plant * plant_per_row
    plot_width = inter_row * nrow
    return plot_length, plot_width
   
def setAdel(nplants, nsect, devT, sowing_density, plant_density, inter_row, seed, sample, dep=7, dynamic_leaf_db=False, leaf_db=None, geoLeaf=None, incT=60, dinT=5,run_adel_pars = {'senescence_leaf_shrink' : 0.5,'startLeaf' : -0.4, 'endLeaf' : 1.6, 'endLeaf1': 1.6, 'stemLeaf' : 1.2,'epsillon' : 1e-6, 'HSstart_inclination_tiller': 1, 'rate_inclination_tiller': 30}, leaf_twist=180):    
     
    length, width = plot_dimension(nplants, sowing_density, plant_density, inter_row)
    nplants, positions, domain, domain_area, convUnit = agronomicplot(length=length, 
                                                            width=width, 
                                                            sowing_density=sowing_density, 
                                                            plant_density=plant_density,
                                                            inter_row=inter_row)
                                                            
    adel = AdelWheat(nplants=nplants, positions = positions, nsect=nsect, devT=devT, seed= seed, sample=sample, dep=dep, dynamic_leaf_db=dynamic_leaf_db, leaf_db=leaf_db, geoLeaf=geoLeaf, incT=incT, dinT=dinT, run_adel_pars = run_adel_pars, leaf_twist=leaf_twist)
    return adel, domain, domain_area, convUnit, nplants
    
    
#---------- reconstructions 

def geoLeaf(nlim=4,dazt=60,dazb=10):
    """ generate geoLeaf function for Adel """
    rcode = """
    geoLeaf <- list(
     Azim = function(a,n,nf) {{
            ntop = nf - n
            ifelse(ntop <= {ntoplim:d},
            180 + {dazTop:.2f} * (runif(1) - .5),
            180 + {dazBase:.2f} * (runif(1) - .5))
            }},
     Lindex = function(a,n,nf) {{
              ntop = nf - n
              ifelse(ntop <= {ntoplim:d}, 2, 1)
              }}
              )
    """
    return rcode.format(ntoplim = nlim, dazTop = dazt, dazBase = dazb)

# macros for general modification   
def new_pgen(pgen, nplants, primary_proba, tdata, pdata, dTT_stop):
        pgen = copy(pgen)
        ears_per_plant = tdata['MB'].values + tdata['TT'].values
        plant_density_at_harvest = float(pdata['ear_density_at_harvest']) / ears_per_plant
        nff = float(tdata['Nff'].values)
        # pgen adapt
        pgen['decide_child_axis_probabilities'] = primary_proba
        pgen['plants_density'] = plant_density_at_harvest #Mariem used plant_density_at_emergence
        pgen['plants_number'] = nplants
        pgen['ears_density'] = pdata['ear_density_at_harvest']
        #hack for removing tillers
        pgen['delais_TT_stop_del_axis'] -= dTT_stop
        return pgen
        
reconst_db={}

def fit_leaves(dxy, disc_level=9):
    import alinea.adel.fitting as fitting
    for k in dxy:
        leaves, discard = fitting.fit_leaves(dxy[k], disc_level)
        dxy[k] = leaves
    return dxy

def Mercia_2010(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=20, **kwds):

    pgen = archidb.Mercia_2010_plantgen()
    tdata = archidb.Tillering_data_Mercia_Rht3_2010_2011()['tillering']
    tdata = tdata[tdata.Var=='Mercia']
    shapes = archidb.leaf_curvature_data('Mercia')
    #pte bidouille pour modifier MB=1 qui ne prend pas en compte les effets de la mouche
    #tdata['MB']=0.9 
    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Mercia']
    #adapt reconstruction
    # TO DO get nff probabilities for main stem from tillering data ?
    if not as_pgen:
        #handle fly damages to MB  by reducing proba of appearance of T2 (originally equal to 1) and simplify TC
        # BUG : MB is equal to 1 => should return  less than 1
        primary_proba={k:tdata[k].values for k in ('T1','T3','T4')}
        primary_proba['T2']= tdata['MB'].values + tdata['TC'].values
        pgen = new_pgen(pgen, nplants, primary_proba, tdata, pdata, dTT_stop)    
    #generate reconstruction
    devT = plantgen_to_devT(pgen)
    # sans angles
    # adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample)
    # avec angles
    leaves = fit_leaves(shapes, disc_level)
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample, dynamic_leaf_db=True, leaf_db=leaves, geoLeaf=geoLeaf(), **kwds)
    
    return pgen, adel, domain, domain_area, convUnit, nplants
    
reconst_db['Mercia'] = Mercia_2010
    
def Mercia_composite_2010(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0):
    from alinea.adel.AdelR import devCsv
    
    pgen = archidb.Mercia_2010_plantgen()
    tdata = archidb.Tillering_data_Mercia_Rht3_2010_2011()['tillering']
    tdata = tdata[tdata.Var=='Mercia']
    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Mercia']
    if not as_pgen:
        primary_proba={k:tdata[k].values for k in ('T1','T3','T4')}
        primary_proba['T2']= tdata['MB'].values + tdata['TC'].values
        # pgen = new_pgen(pgen, nplants, primary_proba, tdata, pdata, dTT_stop) 
        # plantes a 11 talles
        pgen_1 = new_pgen(pgen, 11, primary_proba, tdata, pdata, dTT_stop) 
        pgen_1['MS_leaves_number_probabilities'] = {'11':1.0}
        pgen_1['ears_density'] = (444*11)/31
        pgen_1['plants_density'] = (158.8*11)/31
        # plantes a 12 talles
        pgen_2 = new_pgen(pgen, 19, primary_proba, tdata, pdata, dTT_stop)
        pgen_2['MS_leaves_number_probabilities'] = {'12':1.0}
        pgen_2['ears_density'] = (444*19)/31
        pgen_2['plants_density'] = (158.8*19)/31
        # plantes a 13 talles
        pgen_3 = new_pgen(pgen, 1, primary_proba, tdata, pdata, dTT_stop)
        pgen_3['MS_leaves_number_probabilities'] = {'13':1.0}
        pgen_3['ears_density'] = (444*1)/31
        pgen_3['plants_density'] = (158.8*1)/31
    #-----------------------------------------------------------------------------------------
    #generate reconstruction = MODIFICATION (3 sim plantgen)
    axeT, dimT, phenT = plantgen_to_devT_comp(pgen_1)
    sim = 2
    while sim <= 3 :
        if sim == 2:
            name=pgen_2
        else:
            name=pgen_3
        axeTs, dimTs, phenTs = plantgen_to_devT_comp(name)
        #modification des colonnes id_dim (commune a axeT et dimT) puis id_phen (commune a axeT et phenT)
        tx = sim * 10000
        axeTs['id_plt'] = axeTs['id_plt']+tx;
        axeTs['id_dim'] = axeTs['id_dim']+tx; dimTs['id_dim'] = dimTs['id_dim']+tx
        axeTs['id_phen'] = axeTs['id_phen']+tx; phenTs['id_phen'] = phenTs['id_phen']+tx
        axeT = axeT.append(axeTs, ignore_index=True)
        dimT = dimT.append(dimTs, ignore_index=True)
        phenT = phenT.append(phenTs, ignore_index=True)
        sim = sim + 1
    #-----------------------------------------------------------------------------------------
    devT = devCsv(axeT, dimT, phenT)
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample)
    return pgen, adel, domain, domain_area, convUnit, nplants
    
reconst_db['Mercia_comp'] = Mercia_composite_2010

def Rht3_2010(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=20):

    pgen = archidb.Rht3_2010_plantgen()
    tdata = archidb.Tillering_data_Mercia_Rht3_2010_2011()['tillering']
    tdata = tdata[tdata.Var=='Rht3']
    shapes = archidb.leaf_curvature_data('Rht3')
    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Rht3']
    #adapt reconstruction
    # TO DO get nff probabilities for main stem from tillering data ?
    if not as_pgen:
        primary_proba={k:tdata[k].values for k in ('T1','T2','T3','T4')}
        pgen = new_pgen(pgen, nplants, primary_proba, tdata, pdata, dTT_stop)   
        # plantes a 12 talles
        #pgen['MS_leaves_number_probabilities'] = {'12':1.0}        
    #generate reconstruction
    devT = plantgen_to_devT(pgen)
    # sans angles
    #adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample)
    # avec angles
    leaves = fit_leaves(shapes, disc_level)
    #adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample, dynamic_leaf_db=True, leaf_db=leaves, geoLeaf=geoLeaf())
    # CORRECTION DE LA DENSITE
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=147, plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample, dynamic_leaf_db=True, leaf_db=leaves, geoLeaf=geoLeaf(), incT=22, dinT=22)
   
    return pgen, adel, domain, domain_area, convUnit, nplants
    
reconst_db['Rht3'] = Rht3_2010

def Tremie_2011(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=20):
#BUG: nplants not generated as expected for this genotype !
    pgen = archidb.Tremie_2011_plantgen()
    tdata = archidb.Tillering_data_Tremie1_2011_2012()['tillering']
    tdata = tdata[tdata.Var=='Tremie']
    pdata = archidb.Plot_data_Tremie_2011_2012()['Tremie']
    shapes = archidb.leaf_curvature_data('Tremie')
    tdata['T5']=0; tdata['T6']=0
    #adapt reconstruction
    # TO DO get nff probabilities for main stem from tillering data ?
    if not as_pgen:
        primary_proba={k:tdata[k].values for k in ('T1','T2','T3','T4')} 
        pgen = new_pgen(pgen, nplants, primary_proba, tdata, pdata, dTT_stop) 
    #generate reconstruction
    devT = plantgen_to_devT(pgen)
    #sans angles
    #adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample)
    # avec angles
    leaves = fit_leaves(shapes, disc_level)
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=200, plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample, dynamic_leaf_db=True, leaf_db=leaves, geoLeaf=geoLeaf(), incT=22, dinT=22)
    
    return pgen, adel, domain, domain_area, convUnit, nplants
    
reconst_db['Tremie'] = Tremie_2011

'''
def Tremie_2012(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0):
    pgen = archidb.Tremie_2012_plantgen()
    tdata = archidb.Tillering_data_Tremie_2012_2013()['tillering']
    tdata = tdata[tdata.Var=='Tremie']
    pdata = archidb.Plot_data_Tremie_2012_2013()['Tremie']
    #adapt reconstruction
    # TO DO get nff probabilities for main stem from tillering data ?
    if not as_pgen:
        primary_proba={k:tdata[k].values for k in ('T1','T2','T3','T4','T5')} 
        pgen = new_pgen(pgen, nplants, primary_proba, tdata, pdata, dTT_stop) 
    #generate reconstruction
    devT = plantgen_to_devT(pgen)
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample)
    
    return pgen, adel, domain, domain_area, convUnit, nplants
    
reconst_db['Tremie2'] = Tremie_2012
'''
 
age_at_application = {'Mercia_2010': {'2011-04-19': 1166,
                                      '2011-05-11' : 1500}}
                                
                                   
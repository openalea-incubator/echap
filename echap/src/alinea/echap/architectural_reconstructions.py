""" Reconstruction of Wheat for Boigneville data

Current reconstructions use fit from Mariem (pgen data) + new modification consigned here
"""
import pandas
import numpy
from copy import deepcopy

from alinea.adel.astk_interface import AdelWheat
from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand

import alinea.echap.architectural_data as archidb

#generic function to be moved to adel

def plantgen_to_devT(pgen):
    """ Creates devT tables from plantgen dict
    """
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data, plantgen2adel
    from alinea.adel.AdelR import devCsv
    
    axeT_, dimT_, phenT_, phenT_abs, dimT_abs, dynT_, phenT_first, HS_GL_SSI_T, tilleringT, cardinalityT, config = gen_adel_input_data(**pgen)
    
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
    devT = devCsv(axeT, dimT, phenT)
    return devT
    
def plantgen_to_devT_comp(pgen_pars):
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data,plantgen2adel
    
    axeT_, dimT_, phenT_, _, _, _, _, _, _, _, _ = gen_adel_input_data(**pgen_pars)
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
    return axeT, dimT, phenT
    
    
def check_nff(devT):
    """ Count the probailitiy of occurence of MS with given nff
    """
    nffs = devT['axeT']['N_phytomer'][devT['axeT']['id_axis']=='MS']
    counts = {n:nffs.tolist().count(n) for n in set(nffs)}
    probas = {n:nffs.tolist().count(n) * 1.0 / len(nffs) for n in set(nffs)}
    return counts, probas
    
    
def check_primary_tillers(devT):
    """ Count/estimate probabilitie of occurence of primary tillers
    """
    import re 
    axis = devT['axeT']['id_axis']
    ms = [e for e in axis if re.match('MS',e)]
    tillers = [e for e in axis if re.match('T.$',e)]
    counts = {n:tillers.count(n) for n in set(tillers)}
    probas = {n:tillers.count(n) * 1.0 / len(ms) for n in set(tillers)}
    return counts, probas
    
#---
def sen(pgen):
    """ Creates devT tables from plantgen dict
    """
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data, plantgen2adel
    from alinea.adel.AdelR import devCsv
    
    _, _, _, _, _, _, _, HS_GL_SSI_T, _, _, _ = gen_adel_input_data(**pgen)
    
    # verif de la sortie pour obtenir le graph de SSI
    # HS_GL_SSI_T.to_csv('C:/Users/Administrateur/openaleapkg/echap/test/HS_GL_SSI_T_test.csv')
    
    return HS_GL_SSI_T
#---


   
def setAdel(**kwds):    
     
    adel = AdelWheat(**kwds)    
    return adel, adel.domain, adel.domain_area, adel.convUnit, adel.nplants
    
    
#---------- reconstructions 


# leaf geometry

def geoLeaf(nlim=4,dazt=60,dazb=10, Lindex_base = 1, Lindex_top = 2):
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
              ifelse(ntop <= {ntoplim:d}, {top_class:d}, {base_class:d})
              }}
              )
    """
    return rcode.format(ntoplim = nlim, dazTop = dazt, dazBase = dazb, top_class=Lindex_top, base_class= Lindex_base)


# Attention pour rht3 on veut geoleaf qui retourne 1 (= feuile du bas) quelque soit le rang !
    # creation de la colonne age et du selecteur d'age 
    
    
    
# macros for general modification   
def new_pgen(pgen, nplants, primary_proba, tdata, pdata, dTT_stop):
        pgenc = deepcopy(pgen)
        ears_per_plant = tdata['ears_per_plant']
        plant_density_at_harvest = float(pdata['ear_density_at_harvest']) / ears_per_plant
        #nff = float(tdata['emission_probabilities']['Nff'].values)
        # pgen adapt
        pgenc['decide_child_axis_probabilities'] = primary_proba
        pgenc['plants_density'] = plant_density_at_harvest #Mariem used plant_density_at_emergence
        pgenc['plants_number'] = nplants
        pgenc['ears_density'] = pdata['ear_density_at_harvest']
        #hack for removing tillers
        pgenc['delais_TT_stop_del_axis'] -= dTT_stop
        return pgenc
        
        
# Standard echap reconstruction protocol for plant development and axis dynamic       
#def echap_development(nplants, pgen, primary_proba, tdata, pdata, dTT_stop)
    #ici include mortality
        
reconst_db={}


# NORMAL
def Mercia_2010(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=7, **kwds):

    pgen = archidb.Mercia_2010_plantgen()
    tdb = archidb.Tillering_data_Mercia_Rht3_2010_2011()['Mercia']
    xy, sr, bins = archidb.leaf_curvature_data('Mercia')
    leaves = Leaves(xy, sr, geoLeaf=geoLeaf(), dynamic_bins = bins, discretisation_level = disc_level)
    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Mercia']
    
    stand = AgronomicStand(sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'],inter_row=pdata['inter_row'])
    #pte bidouille pour modifier MB=1 qui ne prend pas en compte les effets de la mouche
    #tdata['MB']=0.9 

    #adapt reconstruction
    # TO DO get nff probabilities for main stem from tillering data ?
    if not as_pgen:
        tdata = tdb['emission_probabilities']
        #handle fly damages to MB  by reducing proba of appearance of T2 (originally equal to 1) and simplify TC
        # BUG : MB is equal to 1 => should return  less than 1
        primary_proba={k:tdata[k] for k in ('T1','T3','T4')}
        #primary_proba['T2']= tdata['MB'].values + tdata['TC'].values
        pgen = new_pgen(pgen, nplants, primary_proba, tdb, pdata, dTT_stop)   
    #generate reconstruction
    devT = plantgen_to_devT(pgen)

    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, stand = stand , seed=seed, sample=sample, leaves = leaves, **kwds)
    
    return pgen, adel, domain, domain_area, convUnit, nplants
    
reconst_db['Mercia'] = Mercia_2010

def Rht3_2010(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=7, face_up=False, classic=False):

    pgen = archidb.Rht3_2010_plantgen()
    tdb = archidb.Tillering_data_Mercia_Rht3_2010_2011()['Rht3']
    xy, sr, bins = archidb.leaf_curvature_data('Rht3')
    leaves = Leaves(xy, sr, geoLeaf=geoLeaf(), dynamic_bins = bins, discretisation_level = disc_level)
    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Rht3']
    stand = AgronomicStand(sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'],inter_row=pdata['inter_row'])
    #adapt reconstruction
    # TO DO get nff probabilities for main stem from tillering data ?
    if not as_pgen:
        tdata = tdb['emission_probabilities']
        primary_proba={k:tdata[k] for k in ('T1','T2','T3','T4')}
        pgen = new_pgen(pgen, nplants, primary_proba, tdb, pdata, dTT_stop)   
        # plantes a 12 talles
        # pgen['MS_leaves_number_probabilities'] = {'12':1.0}        
    #generate reconstruction
    devT = plantgen_to_devT(pgen)
    
    #adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample, leaves = leaves, incT=22, dinT=22, face_up=face_up)
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, stand = stand , seed=seed, sample=sample, leaves = leaves, incT=22, dinT=22, face_up=face_up)
   
    return pgen, adel, domain, domain_area, convUnit, nplants
    
reconst_db['Rht3'] = Rht3_2010

# COMPOSITE
def Mercia_composite_2010(nplants=31, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=7, aborting_tiller_reduction = 1.0, **kwds):
    from alinea.adel.AdelR import devCsv

    pgen1 = archidb.Mercia_2010_nff11_plantgen()
    pgen2 = archidb.Mercia_2010_nff12_plantgen()
    pgen3 = archidb.Mercia_2010_nff13_plantgen()
    
    tdb = archidb.Tillering_data_Mercia_Rht3_2010_2011()['Mercia']
    xy, sr, bins = archidb.leaf_curvature_data('Mercia')
    leaves = Leaves(xy, sr, geoLeaf=geoLeaf(), dynamic_bins = bins, discretisation_level = disc_level)
    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Mercia']
    
    if not as_pgen:
        tdata = tdb['emission_probabilities']
        #handle fly damages to MB  by reducing proba of appearance of T2 (originally equal to 1) and simplify TC
        # BUG : MB is equal to 1 => should return  less than 1
        #primary_proba={k:tdata[k].values for k in ('T1','T3','T4')} #old
        primary_proba={k:tdata[k] for k in ('T1','T3','T4')}
        #primary_proba['T2']= tdata['MB'].values + tdata['TC'].values #old
        #primary_proba['T2']= tdata['T2'].values + tdata['TC'].values
        
        #pgen = new_pgen(pgen, nplants, primary_proba, tdata, pdata, dTT_stop)  
        #pgen['MS_leaves_number_probabilities'] = {'11':0.4,'12':0.6}    
        pgen_1 = new_pgen(pgen1, 11, primary_proba, tdb, pdata, dTT_stop)
        pgen_1['MS_leaves_number_probabilities'] = {'11':1.0}
        pgen_2 = new_pgen(pgen2, 19, primary_proba, tdb, pdata, dTT_stop)
        pgen_2['MS_leaves_number_probabilities'] = {'12':1.0}
        pgen_3 = new_pgen(pgen3, 1, primary_proba, tdb, pdata, dTT_stop)
        pgen_3['MS_leaves_number_probabilities'] = {'13':1.0}
           
    #generate reconstruction
    axeT, dimT, phenT = plantgen_to_devT_comp(pgen_1)

    # renumerotation des colonnes communes
    sim = 2
    while sim <= 3 :
        if sim == 2:
            name=pgen_2; plt=11
        else:
            name=pgen_3; plt=30
        axeTs, dimTs, phenTs = plantgen_to_devT_comp(name)
        tx = sim * 10000
        axeTs['id_plt'] = axeTs['id_plt']+plt
        #modification des colonnes id_dim (commune a axeT et dimT) puis id_phen (commune a axeT et phenT)
        axeTs['id_dim'] = axeTs['id_dim']+tx; dimTs['id_dim'] = dimTs['id_dim']+tx
        axeTs['id_phen'] = axeTs['id_phen']+tx; phenTs['id_phen'] = phenTs['id_phen']+tx
        axeT = axeT.append(axeTs, ignore_index=True)
        dimT = dimT.append(dimTs, ignore_index=True)
        phenT = phenT.append(phenTs, ignore_index=True) 
        sim = sim + 1
    devT = devCsv(axeT, dimT, phenT)

    #verif composition devT
    #print check_nff(devT)
    #print check_primary_tillers(devT)
    
    plant_density = pgen_2['plants_density']
    stand = AgronomicStand(sowing_density=pdata['plant_density_at_emergence'], plant_density=plant_density, inter_row=pdata['inter_row'])
    
    #adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=plant_density, inter_row=pdata['inter_row'], seed=seed, sample=sample, leaves=leaves, aborting_tiller_reduction=aborting_tiller_reduction, **kwds)
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, stand=stand, seed=seed, sample=sample, leaves=leaves, aborting_tiller_reduction=aborting_tiller_reduction, **kwds)
    
    return pgen_2, adel, domain, domain_area, convUnit, nplants
     
reconst_db['Mercia_comp'] = Mercia_composite_2010

#---Maquette interception
'''
def Mercia_composite_2_2010(nplants=31, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=20, aborting_tiller_reduction = 0.6, **kwds):
    from alinea.adel.AdelR import devCsv
    
    pgen1 = archidb.Mercia_2010_nff11_plantgen()
    pgen2 = archidb.Mercia_2010_nff12_plantgen()
    pgen3 = archidb.Mercia_2010_nff13_plantgen()
    
    #tdata = archidb.Tillering_data_Mercia_Rht3_2010_2011()['tillering'] #old
    #tdata = tdata[tdata.Var=='Mercia'] #old
    tdb = archidb.Tillering_data_Mercia_Rht3_2010_2011()['Mercia']
    
    shapes = archidb.leaf_curvature_data('Mercia')
    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Mercia']

    if not as_pgen:
        tdata = tdb['emission_probabilities']
        #primary_proba={k:tdata[k].values for k in ('T1','T3','T4')} #old
        primary_proba={k:tdata[k] for k in ('T1','T3','T4')}
        #primary_proba['T2']= tdata['MB'].values + tdata['TC'].values #old
        #primary_proba['T2']= tdata['T2'].values + tdata['TC'].values
           
        pgen_1 = new_pgen(pgen1, 11, primary_proba, tdb, pdata, dTT_stop)
        pgen_1['MS_leaves_number_probabilities'] = {'11':1.0}
        pgen_2 = new_pgen(pgen2, 19, primary_proba, tdb, pdata, dTT_stop)
        pgen_2['MS_leaves_number_probabilities'] = {'12':1.0}
        pgen_3 = new_pgen(pgen3, 1, primary_proba, tdb, pdata, dTT_stop)
        pgen_3['MS_leaves_number_probabilities'] = {'13':1.0}
           
    #generate reconstruction
    axeT, dimT, phenT = plantgen_to_devT_comp(pgen_1)

    # renumerotation des colonnes communes
    sim = 2
    while sim <= 3 :
        if sim == 2:
            name=pgen_2; plt=11
        else:
            name=pgen_3; plt=30
        axeTs, dimTs, phenTs = plantgen_to_devT_comp(name)
        tx = sim * 10000
        axeTs['id_plt'] = axeTs['id_plt']+plt
        #modification des colonnes id_dim (commune a axeT et dimT) puis id_phen (commune a axeT et phenT)
        axeTs['id_dim'] = axeTs['id_dim']+tx; dimTs['id_dim'] = dimTs['id_dim']+tx
        axeTs['id_phen'] = axeTs['id_phen']+tx; phenTs['id_phen'] = phenTs['id_phen']+tx
        axeT = axeT.append(axeTs, ignore_index=True)
        dimT = dimT.append(dimTs, ignore_index=True)
        phenT = phenT.append(phenTs, ignore_index=True) 
        sim = sim + 1
    devT = devCsv(axeT, dimT, phenT)

    #verif composition devT
    #print check_nff(devT)
    #print check_primary_tillers(devT)
    
    # avec angles
    leaves = fit_leaves(shapes, disc_level)
    
    plant_density = pgen_2['plants_density']
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=plant_density, inter_row=pdata['inter_row'], seed=seed, sample=sample, dynamic_leaf_db=True, leaf_db=leaves, geoLeaf=geoLeaf(), aborting_tiller_reduction=aborting_tiller_reduction, **kwds)
    
    return pgen_2, adel, domain, domain_area, convUnit, nplants
    
def Mercia_composite_3_2010(nplants=31, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=500, disc_level=20, aborting_tiller_reduction = 0.6, **kwds):
    from alinea.adel.AdelR import devCsv
    
    pgen1 = archidb.Mercia_2010_nff11_plantgen()
    pgen2 = archidb.Mercia_2010_nff12_plantgen()
    pgen3 = archidb.Mercia_2010_nff13_plantgen()
    
    #tdata = archidb.Tillering_data_Mercia_Rht3_2010_2011()['tillering'] #old
    #tdata = tdata[tdata.Var=='Mercia'] #old
    tdb = archidb.Tillering_data_Mercia_Rht3_2010_2011()['Mercia']
    
    shapes = archidb.leaf_curvature_data('Mercia')
    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Mercia']

    if not as_pgen:
        tdata = tdb['emission_probabilities']
        #primary_proba={k:tdata[k].values for k in ('T1','T3','T4')} #old
        primary_proba={k:tdata[k] for k in ('T1','T3','T4')}
        #primary_proba['T2']= tdata['MB'].values + tdata['TC'].values #old
        #primary_proba['T2']= tdata['T2'].values + tdata['TC'].values
            
        pgen_1 = new_pgen(pgen1, 11, primary_proba, tdb, pdata, dTT_stop)
        pgen_1['MS_leaves_number_probabilities'] = {'11':1.0}
        pgen_2 = new_pgen(pgen2, 19, primary_proba, tdb, pdata, dTT_stop)
        pgen_2['MS_leaves_number_probabilities'] = {'12':1.0}
        pgen_3 = new_pgen(pgen3, 1, primary_proba, tdb, pdata, dTT_stop)
        pgen_3['MS_leaves_number_probabilities'] = {'13':1.0}
           
    #generate reconstruction
    axeT, dimT, phenT = plantgen_to_devT_comp(pgen_1)

    # renumerotation des colonnes communes
    sim = 2
    while sim <= 3 :
        if sim == 2:
            name=pgen_2; plt=11
        else:
            name=pgen_3; plt=30
        axeTs, dimTs, phenTs = plantgen_to_devT_comp(name)
        tx = sim * 10000
        axeTs['id_plt'] = axeTs['id_plt']+plt
        #modification des colonnes id_dim (commune a axeT et dimT) puis id_phen (commune a axeT et phenT)
        axeTs['id_dim'] = axeTs['id_dim']+tx; dimTs['id_dim'] = dimTs['id_dim']+tx
        axeTs['id_phen'] = axeTs['id_phen']+tx; phenTs['id_phen'] = phenTs['id_phen']+tx
        axeT = axeT.append(axeTs, ignore_index=True)
        dimT = dimT.append(dimTs, ignore_index=True)
        phenT = phenT.append(phenTs, ignore_index=True) 
        sim = sim + 1
    devT = devCsv(axeT, dimT, phenT)

    #verif composition devT
    #print check_nff(devT)
    #print check_primary_tillers(devT)
    
    # avec angles
    leaves = fit_leaves(shapes, disc_level)
    
    plant_density = pgen_2['plants_density']
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=plant_density, inter_row=pdata['inter_row'], seed=seed, sample=sample, dynamic_leaf_db=True, leaf_db=leaves, geoLeaf=geoLeaf(), aborting_tiller_reduction=aborting_tiller_reduction, **kwds)
    
    return pgen_2, adel, domain, domain_area, convUnit, nplants
''' 
def Rht3_composite_1_2010(nplants=24, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=7, face_up=False, aborting_tiller_reduction = 1.0, classic=False, **kwds):
    from alinea.adel.AdelR import devCsv
    
    pgen1 = archidb.Rht3_2010_nff11_plantgen()
    pgen2 = archidb.Rht3_2010_nff12_plantgen()
    
    tdb = archidb.Tillering_data_Mercia_Rht3_2010_2011()['Rht3']
    xy, sr, bins = archidb.leaf_curvature_data('Rht3')
    leaves = Leaves(xy, sr, geoLeaf=geoLeaf(), dynamic_bins = bins, discretisation_level = disc_level)
    shapes = archidb.leaf_curvature_data('Rht3')
    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Rht3']

    if not as_pgen:
        tdata = tdb['emission_probabilities']
        primary_proba={k:tdata[k] for k in ('T1','T2','T3','T4')}
          
        pgen_1 = new_pgen(pgen1, 17, primary_proba, tdb, pdata, dTT_stop)
        pgen_1['MS_leaves_number_probabilities'] = {'11':1.0}
        pgen_2 = new_pgen(pgen2, 7, primary_proba, tdb, pdata, dTT_stop)
        pgen_2['MS_leaves_number_probabilities'] = {'12':1.0}
           
    #generate reconstruction
    axeT, dimT, phenT = plantgen_to_devT_comp(pgen_1)

    # renumerotation des colonnes communes
    sim = 2
    while sim <= 2 :
        if sim == 2:
            name=pgen_2; plt=11
        else:
            name=pgen_3; plt=30
        axeTs, dimTs, phenTs = plantgen_to_devT_comp(name)
        tx = sim * 10000
        axeTs['id_plt'] = axeTs['id_plt']+plt
        #modification des colonnes id_dim (commune a axeT et dimT) puis id_phen (commune a axeT et phenT)
        axeTs['id_dim'] = axeTs['id_dim']+tx; dimTs['id_dim'] = dimTs['id_dim']+tx
        axeTs['id_phen'] = axeTs['id_phen']+tx; phenTs['id_phen'] = phenTs['id_phen']+tx
        axeT = axeT.append(axeTs, ignore_index=True)
        dimT = dimT.append(dimTs, ignore_index=True)
        phenT = phenT.append(phenTs, ignore_index=True) 
        sim = sim + 1
    devT = devCsv(axeT, dimT, phenT)

    #verif composition devT
    #print check_nff(devT)
    #print check_primary_tillers(devT)
    
    plant_density = pgen_1['plants_density']
    stand = AgronomicStand(sowing_density=pdata['plant_density_at_emergence'], plant_density=plant_density, inter_row=pdata['inter_row'])
    
    #adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=plant_density, inter_row=pdata['inter_row'], seed=seed, sample=sample, leaves=leaves, aborting_tiller_reduction=aborting_tiller_reduction, incT=22, dinT=22, face_up=face_up, **kwds)
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, stand=stand, seed=seed, sample=sample, leaves=leaves, aborting_tiller_reduction=aborting_tiller_reduction, incT=22, dinT=22, face_up=face_up,**kwds)
    
    return pgen_1, adel, domain, domain_area, convUnit, nplants
''' 
def Rht3_composite_2_2010(nplants=24, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=20, face_up=True, aborting_tiller_reduction = 1.0, **kwds):
    from alinea.adel.AdelR import devCsv
    
    pgen1 = archidb.Rht3_2010_nff11_plantgen()
    pgen2 = archidb.Rht3_2010_nff12_plantgen()
    
    #tdata = archidb.Tillering_data_Mercia_Rht3_2010_2011()['tillering'] #old
    #tdata = tdata[tdata.Var=='Rht3'] #old
    tdb = archidb.Tillering_data_Mercia_Rht3_2010_2011()['Rht3']
    
    shapes = archidb.leaf_curvature_data('Rht3')
    pdata = archidb.Plot_data_Mercia_Rht3_2010_2011()['Rht3']

    if not as_pgen:
        tdata = tdb['emission_probabilities']
        
        #primary_proba={k:tdata[k].values for k in ('T1','T2','T3','T4')} 
        primary_proba={k:tdata[k] for k in ('T1','T2','T3','T4')}
          
        pgen_1 = new_pgen(pgen1, 17, primary_proba, tdb, pdata, dTT_stop)
        pgen_1['MS_leaves_number_probabilities'] = {'11':1.0}
        pgen_2 = new_pgen(pgen2, 7, primary_proba, tdb, pdata, dTT_stop)
        pgen_2['MS_leaves_number_probabilities'] = {'12':1.0}
           
    #generate reconstruction
    axeT, dimT, phenT = plantgen_to_devT_comp(pgen_1)

    # renumerotation des colonnes communes
    sim = 2
    while sim <= 2 :
        if sim == 2:
            name=pgen_2; plt=11
        else:
            name=pgen_3; plt=30
        axeTs, dimTs, phenTs = plantgen_to_devT_comp(name)
        tx = sim * 10000
        axeTs['id_plt'] = axeTs['id_plt']+plt
        #modification des colonnes id_dim (commune a axeT et dimT) puis id_phen (commune a axeT et phenT)
        axeTs['id_dim'] = axeTs['id_dim']+tx; dimTs['id_dim'] = dimTs['id_dim']+tx
        axeTs['id_phen'] = axeTs['id_phen']+tx; phenTs['id_phen'] = phenTs['id_phen']+tx
        axeT = axeT.append(axeTs, ignore_index=True)
        dimT = dimT.append(dimTs, ignore_index=True)
        phenT = phenT.append(phenTs, ignore_index=True) 
        sim = sim + 1
    devT = devCsv(axeT, dimT, phenT)

    #verif composition devT
    print check_nff(devT)
    print check_primary_tillers(devT)
    
    # avec angles
    leaves = fit_leaves(shapes, disc_level)
    
    plant_density = pdata['plant_density_at_emergence']
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=plant_density, inter_row=pdata['inter_row'], seed=seed, sample=sample, dynamic_leaf_db=True, leaf_db=leaves, geoLeaf=geoLeaf(), aborting_tiller_reduction=aborting_tiller_reduction, incT=22, dinT=22, face_up=face_up, **kwds)
    
    return pgen_1, adel, domain, domain_area, convUnit, nplants
'''

# lancement simulations maquette interception Nathalie
reconst_db['Mercia_maq1'] = Mercia_composite_2010
#reconst_db['Mercia_maq2'] = Mercia_composite_2_2010
#reconst_db['Mercia_maq3'] = Mercia_composite_3_2010
reconst_db['Rht3_maq1'] = Rht3_composite_1_2010
#reconst_db['Rht3_maq2'] = Rht3_composite_2_2010

#---

def Rht3_nomouche_2010(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=20, face_up=True):

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
        # pgen['MS_leaves_number_probabilities'] = {'12':1.0}        
    #generate reconstruction
    devT = plantgen_to_devT(pgen)
    # sans angles
    # adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample)
    # avec angles
    leaves = fit_leaves(shapes, disc_level)
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pdata['plant_density_at_emergence'], inter_row=pdata['inter_row'], seed=seed, sample=sample, dynamic_leaf_db=True, leaf_db=leaves, geoLeaf=geoLeaf(), incT=22, dinT=22, face_up=face_up)
   
    return pgen, adel, domain, domain_area, convUnit, nplants
    
reconst_db['Rht3_nomouche'] = Rht3_nomouche_2010

def Tremie12(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=7):
#BUG: nplants not generated as expected for this genotype !
    pgen = archidb.Tremie_2011_plantgen()
    tdb = archidb.Tillering_data_Tremie12_2011_2012()
    xy, sr, bins = archidb.leaf_curvature_data('Tremie')
    leaves = Leaves(xy, sr, geoLeaf=geoLeaf(), dynamic_bins = bins, discretisation_level = disc_level)
    pdata = archidb.Plot_data_Tremie_2011_2012()
    
    #stand = AgronomicStand(sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'],inter_row=pdata['inter_row']) # => pas de pdata['plant_density_at_emergence']
    stand = AgronomicStand(sowing_density=pdata['sowing_density'], plant_density=pgen['plants_density'],inter_row=pdata['inter_row'])

    #tdata['T5']=0; tdata['T6']=0
    #adapt reconstruction
    # TO DO get nff probabilities for main stem from tillering data ?
    if not as_pgen:
        #primary_proba={k:tdata[k].values for k in ('T1','T2','T3','T4')} #old
        tdata = tdb['emission_probabilities']
        primary_proba={k:tdata[k] for k in ('T1','T2','T3','T4')}
        pgen = new_pgen(pgen, nplants, primary_proba, tdb, pdata, dTT_stop)
    #generate reconstruction
    devT = plantgen_to_devT(pgen)

    #adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['sowing_density'], plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample, leaves=leaves, incT=22, dinT=22)
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, stand = stand , seed=seed, sample=sample, leaves = leaves, incT=22, dinT=22)
    
    return pgen, adel, domain, domain_area, convUnit, nplants
    
reconst_db['Tremie12'] = Tremie12
'''
def Tremie12_nogel(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=20):
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
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pdata['plant_density_at_emergence'], inter_row=pdata['inter_row'], seed=seed, sample=sample, dynamic_leaf_db=True, leaf_db=leaves, geoLeaf=geoLeaf(), incT=22, dinT=22)
    
    return pgen, adel, domain, domain_area, convUnit, nplants
    
reconst_db['Tremie12_nogel'] = Tremie12_nogel
'''

def Tremie13(nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, disc_level=7):
    # on prend le pgen de tremie1
    pgen = archidb.Tremie_2012_plantgen()
    
    tdb = archidb.Tillering_data_Tremie13_2012_2013()
    xy, sr, bins = archidb.leaf_curvature_data('Tremie')
    leaves = Leaves(xy, sr, geoLeaf=geoLeaf(), dynamic_bins = bins, discretisation_level = disc_level)
    pdata = archidb.Plot_data_Tremie_2012_2013()
    
    stand = AgronomicStand(sowing_density=pdata['sowing_density'], plant_density=pgen['plants_density'],inter_row=pdata['inter_row'])
    
    #adapt reconstruction
    # TO DO get nff probabilities for main stem from tillering data ?
    if not as_pgen:
        tdata = tdb['emission_probabilities']
        tdb['ears_per_plant'] = 1.*pdata['ear_density_at_harvest' ]/pdata['sowing_density']
        primary_proba={k:tdata[k] for k in ('T1','T2','T3','T4','T5')}
        pgen = new_pgen(pgen, nplants, primary_proba, tdb, pdata, dTT_stop) 
        
    pgen['MS_leaves_number_probabilities'] = {'11':0.71,'12':0.29}
    
    #generate reconstruction
    devT = plantgen_to_devT(pgen)

    #adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, sowing_density=pdata['plant_density_at_emergence'], plant_density=pgen['plants_density'], inter_row=pdata['inter_row'], seed=seed, sample=sample, dynamic_leaf_db=True, leaf_db=leaves, geoLeaf=geoLeaf(), incT=22, dinT=22)
    adel, domain, domain_area, convUnit, nplants = setAdel(nplants = nplants, nsect=nsect, devT=devT, stand = stand , seed=seed, sample=sample, leaves = leaves)
    
    return pgen, adel, domain, domain_area, convUnit, nplants
    
reconst_db['Tremie13'] = Tremie13

 
age_at_application = {'Mercia_2010': {'2011-04-19': 1166,
                                      '2011-05-11' : 1500}}
                                
                                   
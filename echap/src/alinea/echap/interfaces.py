"""Defines interfaces between g and the different models of echap project
"""

def pesticide_interception(g, scene, interception_model, product_name, dose, label='LeafElement'):
    """ 
    Interface between g and the interception model

    :Parameters:
    ----------
    - 'g' : MTG representing the canopy (and the soil)
    - 'scene'
    - 'interception_model' : A class embending the interception model and provide the following methods:    
        - 'interception_model.intercept(scene, product_name, dose)' : Return the dictionnary of scene_id: compound name of the product and surfacic doses (g.m-2)
        See :func:`~alinea.echap.interception_leaf.CaribuInterceptModel`
    - 'product_name' : Commercial name of the product 
    - 'dose' : Dose of product use in field (l.ha)
    
    :Returns:  
    --------
    - 'g' : Updated MTG representing the canopy (and the soil). 'surfacic_doses' property is added to g or updated if present.
    
    :Example:
    -------
      >>> g = MTG()
      >>> scene = plot3d(g)  
      >>> interception_model = CaribuInterceptModel()
      >>> pesticide_interception(g, scene, interception_model, product_name, dose)
      >>> return g
    """
    surf_dose = interception_model.intercept(product_name, dose, scene)
    if not 'penetrated_doses' in g.properties():
        vi = [vid for vid in surf_dose if g.label(vid).startswith(label)]   
        g.add_property('penetrated_doses')
        g.property('penetrated_doses').update(dict((v, {'Chlorothalonil': 0, 'Epoxiconazole': 0}) for v in vi))
    if not 'surfacic_doses' in g.properties():
        vi = [vid for vid in surf_dose if g.label(vid).startswith(label)]   
        g.add_property('surfacic_doses')
        g.property('surfacic_doses').update(dict((v, {'Chlorothalonil': 0, 'Epoxiconazole': 0}) for v in vi))
    surfacic_doses = g.property('surfacic_doses')
    for vid, nd in surf_dose.iteritems():
        if g.label(vid).startswith(label):
            for name, dose in nd.iteritems():
                if name in surfacic_doses[vid]:
                    surfacic_doses[vid][name] += dose
                else:
                    g.add_property('surfacic_doses')
                    g.property('surfacic_doses')[vid].update(nd)
    return g


def local_microclimate(g, scene, climate_model, mean_globalclimate, label='LeafElement'):
    """ 
    Interface between g and the microclimate model

    :Parameters:
    ----------
    - 'g' : MTG representing the canopy (and the soil)
    - 'scene'
    - 'interception_model' : A class embending the interception model and provide the following methods:    
        - 'climate_model.microclim(scene, product_name, dose)' : Return the dictionnary of scene_id: compound name of the product and surfacic doses (g.m-2)
        See :func:`~alinea.echap.interception_leaf.CaribuInterceptModel`
    - 'climate_model' : A class embending the microclimate model and provide the following methods:    
        - 'climate_model.microclim(microclimate, rain, scene)' : Return the dictionnary of scene_id: radiation and rain
        See :func:`~alinea.echap.interception_leaf.CaribuMicroclimModel`
    - 'energy' : 
    
    :Returns:  
    --------
    - 'g' : Updated MTG representing the canopy (and the soil). 'microclimate' property is added to g or updated if present.
    
    :Example:
    -------
      >>> g = MTG()
      >>> scene = plot3d(g)  
      >>> climate_model = CaribuMicroclimModel()
      >>> microclimate(g, scene, climate_model, rain)
      >>> return g
    """
    local_meteo = climate_model.microclim(mean_globalclimate, scene)
    g.add_property('microclimate')
    g.property('microclimate').update(local_meteo)
    return g


def pesticide_surfacic_decay(g, decay_model, label='LeafElement', timestep=24):
    """ 
    Interface between g and the decay model of Pearl

    :Parameters:
    ----------
    - 'g' : MTG representing the canopy (and the soil), doses are stored in the MTG as a property
    - 'decay_model' : A class embending the decay model and provide the following methods:    
        - 'decay_model.decay_and_penetrate(compound_name, compound_dose, microclimate, timestep)': Return for one compound the decayed surfacic dose (g.m-2), the penetrated amount and the loss to the environment.
        See :func:`~alinea.pearl.pearl_leaf.PearLeafDecayModel`

    :Returns:    
    ----------
    - 'g' : Updated MTG representing the canopy (and the soil). 'surfacic_doses' and 'penetrated_doses' properties are added to g or updated if present.

    :Example:
    -------
      >>> g = MTG()
      >>> scene = plot3d(g)  
      >>> db = {'Chlorothalonil':{}, 'Epoxiconazole':{}, 'Metconazole':{}}
      >>> interception_model = CaribuInterceptModel()
      >>> pesticide_interception(g, scene, interception_model, product_name, dose)
      >>> decay_model = PearLeafDecayModel(db)'Tair':temperature[vid]    temperature = g.property('temp')
      >>> g = pesticide_surfacic_decay(g, decay_model)
      >>> return g
    """
    surfacic_doses = g.property('surfacic_doses') 
    microclimate = g.property('microclimate') 
    if not 'penetrated_doses' in g.property_names():
        g.add_property('penetrated_doses')
    penetrated_doses = g.property('penetrated_doses')
    for vid, d in surfacic_doses.iteritems():
        if g.label(vid).startswith(label):
            for compound_name,compound_dose in d.iteritems():
                new_dose,penetrated_amount,loss = decay_model.decay_and_penetrate(compound_name,compound_dose,microclimate[vid],timestep)
                surfacic_doses[vid][compound_name] = new_dose
                if vid in penetrated_doses:
                    if compound_name in penetrated_doses:
                        penetrated_doses[vid][compound_name] += penetrated_amount
                    else:
                        penetrated_doses[vid].update({compound_name:penetrated_amount})
                else:
                    penetrated_doses[vid] = {compound_name : penetrated_amount}
    return g


def pesticide_penetrated_decay(g, decay_model, label='LeafElement', timestep=1):
    """ 
    Interface between g and the decay model of penetrated doses of pesticide

    :Parameters:
    ----------
    - 'g' : MTG representing the canopy (and the soil), doses are stored in the MTG as a property
    - 'decay_model' : : A class embending the penetrated doses decay model and provide the following methods:    
        - 'decay_model.decay(d,timestep)' : Model of penetrated pesticide decay. Return for one compound the dictionnary of compound name of the productof and the decayed penetrated dose (g.m-2).
        See :func:`~alinea.echap.milne_leaf.PenetratedDecayModel`, :func:`~alinea.simcycle.pesticide`

    :Returns:
    -------
    - 'g' : Updated MTG representing the canopy (and the soil). 'penetrated_doses' property is updated to g.

    :Example:
    -------   
      >>> g = MTG()
      >>> scene = plot3d(g)  
      >>> db = {'Chlorothalonil':{}, 'Epoxiconazole':{}, 'Metconazole':{}}
      >>> interception_model = CaribuInterceptModel()
      >>> pesticide_interception(g, scene, interception_model, product_name, dose)
      >>> decay_model = PearLeafDecayModel(db)
      >>> g = pesticide_surfacic_decay(g, decay_model)
      >>> decay_model = PenetratedDecayModel()
      >>> g = pesticide_penetrated_decay(g, decay_model)
      >>> return g      
    """
    penetrated_doses = g.property('penetrated_doses')
    for vid, d in penetrated_doses.iteritems():
        if g.label(vid).startswith(label):
            penetrated_doses = decay_model.decay(d,timestep)
    return g


def pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1):
    """ 
    Interface between g and the efficacy model of pesticide compounds

    :Parameters:
    ----------
    - 'g' : MTG representing the canopy (and the soil), doses are stored in the MTG as a property
    - 'efficacy_model' : : A class embending the efficacy model of surfacic and penetrated doses of pesticide and provide the following methods:    
        - 'efficacy_model.efficacy(surfacic_doses, penetrated_doses, timestep)' : Model of pesticide efficacy. Return the dictionnary of protectant and eradicant efficacy (value between 0 and 1) respectively for surfacic doses and penetrated doses of pesticide.
        See :func:`~alinea.pesticide_efficacy.pesticide_efficacy.PesticideEfficacyModel`, :func:`~alinea.simcycle.pesticide`

    :Returns:
    -------
    - 'g' : Updated MTG representing the canopy (and the soil). 'global_efficacy' property is updated to g.

    :Example:
    -------   
      >>> g = MTG()
      >>> scene = plot3d(g)  
      >>> interception_model = CaribuInterceptModel()
      >>> pesticide_interception(g, scene, interception_model, product_name, dose)
      >>> efficacy_model = PesticideEfficacyModel()
      >>> pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
      >>> return g      
    """
    if not 'surfacic_doses' in g.properties(): 
        vids = [vid for vid in g if g.label(vid).startswith(label)]
        for v in vids : 
            n = g.node(v)
            n.surfacic_doses = {'Chlorothalonil':0,'Epoxiconazole':0}
            n.penetrated_doses = {'Chlorothalonil':0,'Epoxiconazole':0}
    penetrated_doses = g.property('penetrated_doses')
    surfacic_doses = g.property('surfacic_doses') 
    if not 'global_efficacy' in g.properties():
        g.add_property('global_efficacy')
    vids = [vid for vid in surfacic_doses if g.label(vid).startswith(label)]
    for v in vids :     
        g.property('global_efficacy').update({v: efficacy_model.efficacy(surfacic_doses[v], penetrated_doses[v], timestep)})
    return g


#    for vid, surf_doses in surfacic_doses.iteritems():
#        if g.label(vid).startswith(label):
#            g.property('global_efficacy')[vid] = efficacy_model.efficacy(surfacic_doses[vid], penetrated_doses[vid], timestep)
#    return g




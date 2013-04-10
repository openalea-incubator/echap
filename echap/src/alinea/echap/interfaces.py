"""defines interfaces between g and the different models of echap project
"""

def pesticide_interception(g, scene, interception_model, label='LeafElement'):
    """ Interface between g and the interception model of Caribu
    :Parameters:
    ----------
    - `g` : MTG representing the canopy (and the soil) doses are stored in the MTG as a property
    - `interception_model` : Caribu model

    :Returns:  
    --------
    - `g` : Updated MTG representing the canopy (and the soil)
    
    Example
    -------
      >>> g = MTG()
      >>> scene = plot3d(g)  
      >>> interception_model = CaribuInterceptModel()
      >>> pesticide_interception(g, scene, interception_model)
      >>> return g
    """
    surf_dose = interception_model.intercept(scene)
    if not 'surfacic_doses' in g.properties():
        g.add_property('surfacic_doses')
        g.property('surfacic_doses').update(surf_dose)
    else:
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
    
# todo: hauteur de crop à récuperer pour faire tourner pearl

def pesticide_surfacic_decay(g, decay_model, label='LeafElement', timestep=24):
    """ Interface between g and the decay model of Pearl
    :Parameters:
    ----------
    - `g` : MTG representing the canopy (and the soil) doses are stored in the MTG as a property
    - `decay_model` :  model that compute for one compound the decayed surfacic dose (g.m2-1), the penetrated amount and the loss to the environment.
    -  doses is a dict compound_name: ammount    
    :Returns:    
    ----------
    g : Updated MTG representing the canopy (and the soil)
    
    Example
    -------
      >>> g = MTG()
      >>> dose = 1
      >>> coumpound = 
      >>> initiate(g, dose, compound)
      >>> dt = 1
      >>> nb_steps = 1000
      >>> for i in range(nb_steps):
      >>>     update_climate(g)
      >>>     infect(g, dt)
      >>>     update(g,dt)
      >>> return g
    
    """
    surfacic_doses = g.property('surfacic_doses') 
    temperature = g.property('temp')
    if not 'penetrated_doses' in g.property_names():
        g.add_property('penetrated_doses')
    penetrated_doses = g.property('penetrated_doses')
    for vid, d in surfacic_doses.iteritems():
        if g.label(vid).startswith(label):
            for compound_name,compound_dose in d.iteritems():
                microclimate = {'Tair':temperature[vid]}
                new_dose,penetrated_amount,loss = decay_model.decay_and_penetrate(compound_name,compound_dose,microclimate,timestep)
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
    """ Interface between g and the decay model of penetrated doses of pesticide
    :Parameters:
    ----------
    - 'g' : MTG representing the canopy (and the soil) doses are stored in the MTG as a property
    - `decay_model` : Model of penetrated pesticide decay (see simcycle.pesticide)
    
    :Returns:
    -------
    g : Updated MTG representing the canopy (and the soil)
    
    Example
    -------    
    """
    penetrated_doses = g.property('penetrated_doses')
    for vid, d in penetrated_doses.iteritems():
        if g.label(vid).startswith(label):
            penetrated_doses = decay_model.decay(d,timestep)
    return g

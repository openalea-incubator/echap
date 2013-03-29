"""defines interfaces between g and the different models of echap project
"""

def pesticide_surfacic_decay(g, decay_model, label='LeafElement'):
    """ Interface between g and the decay model of Pearl
    Parameters:
    ----------
    - `g` : MTG representing the canopy (and the soil) doses are stored in the MTG as a property
    - `decay_model` : Pearl model that allows compute the fate of pesticides on plant leaves and loss to the environment.
    -  doses is a dict compound_name: ammount    
    Returns    
    ----------
    g : Updated MTG representing the canopy (and the soil)
    """
    surfacic_doses = g.property('surfacic_doses') 
    penetrated_doses = g.property('penetrated_doses') 
    for vid, d in surfacic_doses.iteritems():
        if g.label(vid).startswith(label):
            for compound_name,compound_dose in d.iteritems():
                leaf = g.node(vid)
                microclimate = {'Tair':leaf.temp}
                surfacic_doses[vid][compound_name] = decay_model.decay(compound_name,compound_dose,microclimate)[0]
    for vid, d in penetrated_doses.iteritems():
        if g.label(vid).startswith(label):
            for compound_name,compound_dose in d.iteritems():
                leaf = g.node(vid)
                microclimate = {'Tair':leaf.temp}
                penetrated_doses[vid][compound_name] = decay_model.decay(compound_name,compound_dose,microclimate)[3]
    return g

    
def pesticide_penetrated_decay(g, decay_model, label='LeafElement'):
    """ Interface between g and the decay model of penetrated doses of pesticide
    Parameters
    ----------
    - 'g' : MTG representing the canopy (and the soil) doses are stored in the MTG as a property
    - `decay_model` : Model of penetrated pesticide decay (see simcycle.pesticide)
    Returns
    -------
    g : Updated MTG representing the canopy (and the soil)
    """
    vids = [vid for vid in g if g.label(vid).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.penetrated_doses_deg = n.penetrated_doses.copy()
        n.penetrated_active_doses = decay_model.decay(g.property('penetrated_doses')[v])
    return g


def update_surf(g, dt, label="LeafElement"):
    """ Update the surfacic doses of pesticide on the MTG.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy (and the soil)
        'doses' are stored in the MTG as a property
    dt: int
        Time step of the simulation

    Returns
    -------
    g: MTG
        Updated MTG representing the canopy (and the soil)
    
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
    doses = g.property('surfacic_doses')
    labels = g.property('label') 
    
    for vid, d in doses.iteritems():
        d = [dos for dos in d if dos.active]
        doses[vid] = d
        for doses in d:
            # proposition 1 on fait ici une correspondance nom attendus dans g <-> noms caracterisant un environnement de lesion (classe a faire ??)
            leaf=g.node(vid)
            #proposition 2 : les conventions de nommage noms sont definies ailleurs (ou???) et on passe juste une leaf qui repond a cette convention
            doses.update(dt, leaf)
          
    return g,

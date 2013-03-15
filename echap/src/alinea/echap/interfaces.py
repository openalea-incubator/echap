"""defines interfaces between g and the different models of echap project
"""



def pesticide_surfacic_decay(g, 
             decay_model,
             label="LeafElement"):
    """ Disperse spores of the lesions of fungus identified by fungus_name.

    :Parameters:
      - `g` : MTG representing the canopy (and the soil).
      - `decay_model` : model that allows positioning each DU in stock on g.

    :Example:
      >>> 
    doses is a dict compound_name:ammount
    """
    doses = g.property('surfacic_doses') 
    for vid, d in doses.iteritems():
        if g.label(vid).startswith(label):
            for compound_name,compound_dose in d.iteritems():
                leaf = g.node(vid)
                microclimate = {'Tair':leaf.temp}
                doses[vid][compound_name] = decay_model.decay(compound_name,compound_dose,microclimate)
    return g

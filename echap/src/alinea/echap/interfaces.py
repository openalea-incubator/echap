"""defines interfaces between g and the different models of echap project
"""



def pesticide_surfacic_decay(g, 
             decay_model,
             label="LeafElement"):
    """ Interface between g and the decay model of Pearl

    Parameters:
      - `g` : MTG representing the canopy (and the soil).
      - `decay_model` : Pearl model that allows compute the fate of pesticides on plant leaves and loss to the environment.

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

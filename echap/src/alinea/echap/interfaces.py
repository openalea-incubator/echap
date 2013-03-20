"""defines interfaces between g and the different models of echap project
"""
from alinea.pearl.pearl_leaf import *


def pesticide_surfacic_decay(g, decay_model, label="LeafElement"):
    """ Interface between g and the decay model of Pearl

    Parameters:
      - `g` : MTG representing the canopy (and the soil).
      - `decay_model` : Pearl model that allows compute the fate of pesticides on plant leaves and loss to the environment.

    doses is a dict compound_name:ammount
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

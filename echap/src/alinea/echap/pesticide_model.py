''' Models for pesticide, interception, decay and afficacy
'''

from alinea.pearl.pearl_leaf import PearLeafDecayModel
from alinea.echap.milne_leaf import PenetratedDecayModel
from alinea.pesticide_efficacy.pesticide_efficacy import PesticideEfficacyModel
from alinea.echap.interception_leaf import InterceptModel, pesticide_applications

import alinea.echap.pesticide_data as pest_db

from alinea.echap.interfaces import pesticide_surfacic_decay, pesticide_penetrated_decay, pesticide_efficacy, pesticide_interception

surfacic_decay_model = PearLeafDecayModel(pest_db.pearl_parameters)
#penetrated_decay_model = PenetratedDecayModel(pest_db.milne_parameters)
penetrated_decay_model = PenetratedDecayModel()
efficacy_model = PesticideEfficacyModel()

pesticide_interception_model = InterceptModel(pest_db.products_data)

def update_pesticides(g, weather_data):
    g = pesticide_surfacic_decay(g, surfacic_decay_model, weather_data)
    g = pesticide_penetrated_decay(g, penetrated_decay_model, weather_data)
    g = pesticide_efficacy(g, efficacy_model, weather_data)
    return g
    
def pesticide_intercept(g, application_data, label='LeafElement'):
    return pesticide_interception(g, pesticide_interception_model, application_data, label='LeafElement')
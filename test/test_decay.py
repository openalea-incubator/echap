""" Test decay function calls on a simplistic static structure"""

from alinea.echap.interfaces import pesticide_surfacic_decay
has_pearl = True
try:
    from alinea.pearl.pearl_leaf import PearLeafDecayModel
except ImportError:
    has_pearl = False
    

from alinea.astk.Weather import sample_weather
import alinea.adel.data_samples as adel_samples
from alinea.echap.add_property import add_surfacic_doses, add_microclimate


def test_global_climate(dt = 12):
    seq, weather = sample_weather(dt)
    g = adel_samples.adel_two_metamers()
    add_surfacic_doses(g)
    global_climate = weather.get_weather(seq)

    if has_pearl:
        model = PearLeafDecayModel()
        pesticide_surfacic_decay(g,model, global_climate)
    return g
    
def test_local_climate(dt = 12):
    seq, weather = sample_weather(dt)
    g = adel_samples.adel_two_metamers()
    add_surfacic_doses(g)
    global_climate = weather.get_weather(seq)
    add_microclimate(g,global_climate)
    if has_pearl:
        model = PearLeafDecayModel()

        pesticide_surfacic_decay(g,model, global_climate)
    return g
 
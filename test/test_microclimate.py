import openalea #add opeanalea for conda install package with develop package. remove when bugfix
from alinea.adel.astk_interface import initialise_stand
from alinea.astk.Weather import sample_weather
import alinea.adel.data_samples as adel_samples

from alinea.echap.microclimate_leaf import microclimate_leaf


def test_microclimate_leaf(dt=12):
    seq, weather = sample_weather(dt)
    wdata = weather.get_weather(seq)
    g, domain_area, domain, convUnit = adel_samples.adel_two_metamers_stand()
    microclimate_leaf(g,wdata, domain=domain, convUnit = convUnit)
    return g

def test_microclimate_leaf_realistic(dt=12):
    seq, weather = sample_weather(dt)
    wdata = weather.get_weather(seq)
    g, wheat, domain_area, domain, convUnit = initialise_stand(age=200)
    microclimate_leaf(g,wdata, domain=domain, convUnit = convUnit)
    return g



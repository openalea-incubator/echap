from openalea.deploy.shared_data import shared_data
import alinea.echap.architectural_data as archidb
from alinea.echap.architectural_reconstructions import reconst_db
#from alinea.adel.AdelR import devCsv
#import alinea.echap

'''
def devT_Mercia():
    axeT = shared_data(alinea.echap, 'Mercia_axeT.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT.csv')
    earT = shared_data(alinea.echap, 'Mercia_earT.csv')
    phenT = shared_data(alinea.echap, 'Mercia_phenT.csv')
    ssisenT = shared_data(alinea.echap, 'ssi2sen.csv')
    return devCsv(axeT,dimT,phenT,earT,ssisenT)
'''
    
from alinea.adel.stand.stand import agronomicplot
from alinea.adel.astk_interface import AdelWheat

def get_reconstruction(name='Mercia', **args):
    fun = reconst_db[name]
    pgen, adel, domain, domain_area, convUnit, nplants = fun(**args)
    return pgen, adel, domain, domain_area, convUnit, nplants
    
def reconst(name='Mercia', dTT_stop=0, original=False, n=30):
    pgen, adel, domain, domain_area, convUnit, nplants = get_reconstruction(name, nplants = n, dTT_stop=dTT_stop, as_pgen=original)
    return pgen, adel, domain, domain_area, convUnit, nplants

'''
#def setup_canopy(age = 300, nsect=3, length = 0.2, width=0.15, sowing_density=150, plant_density=150, inter_row=0.125, seed=1, sample='random'):
def setup_canopy(age = 1166, nsect=3, length = 0.4, width=0.30, sowing_density = 150, plant_density= 150, inter_row=0.125, seed=1, sample='random'):    
    #devT = devT_Mercia()
    
    pgen, _, _, _, _, _ = get_reconstruction(name, nplants = n, dTT_stop=dTT_stop, as_pgen=original)
    nplants, positions, domain, domain_area, convUnit = agronomicplot(length=length, 
                                                            width=width, 
                                                            sowing_density=sowing_density, 
                                                            plant_density=plant_density,
                                                            inter_row=inter_row)
    adel = AdelWheat(nplants=nplants, positions = positions, nsect=nsect, devT=devT, seed= seed, sample=sample)
    g = adel.setup_canopy(age)
    return g, adel, domain, domain_area, convUnit, nplants
'''

from alinea.alep.protocol import update
from alinea.septo3d.cycle.alep_objects import Septo3DFungus
fungus = Septo3DFungus()

def update_lesions(g,dt,activate=True):
    return update(g,dt)
    
from alinea.echap.interfaces import pesticide_surfacic_decay, pesticide_penetrated_decay, pesticide_efficacy, pesticide_interception

from alinea.pearl.pearl_leaf import PearLeafDecayModel
db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}}
PearlLeaf = PearLeafDecayModel(db)

from alinea.echap.milne_leaf import PenetratedDecayModel
Milne_leaf = PenetratedDecayModel()

from alinea.pesticide_efficacy.pesticide_efficacy import PesticideEfficacyModel
Milne_efficacy = PesticideEfficacyModel()

def update_pesticides(g, weather_data):
    g = pesticide_surfacic_decay(g, PearlLeaf, weather_data)
    g = pesticide_penetrated_decay(g, Milne_leaf, weather_data)
    g = pesticide_efficacy(g, Milne_efficacy, weather_data)
    return g
  
  
from alinea.alep.protocol import disperse
from alinea.popdrops.alep_interface import PopDropsEmission, PopDropsTransport

from alinea.septo3d.dispersion.alep_interfaces import Septo3DEmission, Septo3DTransport

def dispersion(g, weather_data, domain, domain_area, convUnit):
    emission = Septo3DEmission(domain=domain, convUnit = convUnit)
    transport = Septo3DTransport(domain = domain, domain_area =domain_area, convUnit=convUnit)
    g,dus = disperse(g,emission, transport, fungus.name, weather_data=weather_data)
    return g
 
 
from alinea.alep.protocol import external_contamination
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum, Septo3DSoilContamination

def contamination(g, weather_data, level, domain, domain_area, convUnit):
    inoc = SoilInoculum(fungus, sporulating_fraction =level, domain_area = domain_area, convUnit=convUnit)
    contaminator = Septo3DSoilContamination(domain = domain, domain_area =domain_area, convUnit=convUnit)
    g = external_contamination(g, inoc, contaminator, weather_data)
    return g
  
  
from alinea.echap.interception_leaf import InterceptModel
productsDB={'Opus': {'Epoxiconazole': 125}, 'Banko 500': {'Chlorothalonil': 500}}
interception = InterceptModel(productsDB)

def pesticide_intercept(g, application_data, label='LeafElement'):
    return pesticide_interception(g, interception, application_data, label='LeafElement')
from alinea.adel.AdelR import devCsv
from openalea.deploy.shared_data import shared_data
import alinea.echap

def devT_Mercia():
    axeT = shared_data(alinea.echap, 'Mercia_axeT.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT.csv')
    earT = shared_data(alinea.echap, 'Mercia_earT.csv')
    phenT = shared_data(alinea.echap, 'Mercia_phenT.csv')
    ssisenT = shared_data(alinea.echap, 'ssi2sen.csv')
    
    return devCsv(axeT,dimT,phenT,earT,ssisenT)
    
from alinea.adel.stand.stand import agronomicplot
from alinea.adel.astk_interface import AdelWheat

def setup_canopy(age = 300, nsect=3, length = 0.2, width=0.15, sowing_density=150, plant_density=150, inter_row=0.125, seed=1, sample='random'):
    devT = devT_Mercia()
    nplants, positions, domain, domain_area, convUnit = agronomicplot(length=length, 
                                                            width=width, 
                                                            sowing_density=sowing_density, 
                                                            plant_density=plant_density,
                                                            inter_row=inter_row)
    adel = AdelWheat(nplants=nplants, positions = positions, nsect=nsect, devT=devT, seed= seed, sample=sample)
    g = adel.setup_canopy(age)
    return g, adel, domain, domain_area, convUnit, nplants

from alinea.alep.protocol import update
from alinea.septo3d.cycle.alep_objects import GrowthControlModel
growth_control = GrowthControlModel()

def update_lesions(g,dt,activate=True):
    return update(g,dt,growth_control)
    
from alinea.echap.interfaces import pesticide_surfacic_decay, pesticide_penetrated_decay, pesticide_efficacy

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
    return pesticide_efficacy(g, Milne_efficacy, weather_data)
    
from alinea.alep.protocol import disperse
from alinea.septo3d.dispersion.alep_interfaces import Septo3DEmission, Septo3DTransport

def dispersion(g, weather_data, domain, domain_area, convUnit):
    emission = Septo3DEmission(domain=domain, convUnit = convUnit)
    transport = Septo3DTransport(domain = domain, domain_area =domain_area, convUnit=convUnit)
    g,dus = disperse(g,emission, transport, 'septo3d', weather_data=weather_data)
    return g
 
from alinea.alep.protocol import external_contamination
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum, Septo3DSoilContamination
from alinea.septo3d.cycle.alep_objects import Septo3D_DU
generator = Septo3D_DU()
 
def contamination(g, weather_data, level, domain, domain_area, convUnit):
    inoc = SoilInoculum(DU_generator=generator, sporulating_fraction =level, domain_area = domain_area, convUnit=convUnit)
    contaminator = Septo3DSoilContamination(domain = domain, domain_area =domain_area, convUnit=convUnit)
    return external_contamination(g, inoc, contaminator, weather_data)
    
    
from recorder import LeafElementRecorder
from alinea.echap.interfaces import record as do_record #to avoid confusion with numpy record
echap_recorder = LeafElementRecorder()

from architectural_reconstructions import Mercia_2010, age_at_application

from alinea.alep.protocol import update
from alinea.septo3d.cycle.alep_objects import GrowthControlModel
growth_control = GrowthControlModel()

def update_lesions(g,dt,activate=True):
    return update(g,dt,growth_control)
    
from alinea.echap.interfaces import pesticide_surfacic_decay, pesticide_penetrated_decay, pesticide_efficacy, pesticide_interception

from alinea.pearl.pearl_leaf import PearLeafDecayModel
from pesticide_data import pearl_parameters
PearlLeaf = PearLeafDecayModel(pearl_parameters)

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
from alinea.septo3d.dispersion.alep_interfaces import Septo3DEmission, Septo3DTransport
from alinea.septo3d.dispersion.alep_interfaces import SimpleTransport

def dispersion(g, weather_data, domain, domain_area, convUnit):
    emission = Septo3DEmission(domain=domain, convUnit = convUnit)
    transport = Septo3DTransport(domain = domain, domain_area =domain_area, convUnit=convUnit)
    g,dus = disperse(g,emission, transport, 'septo3d', weather_data=weather_data)
    return g
 
def simple_dispersion(g, weather_data, domain, domain_area, convUnit):
    emission = Septo3DEmission(domain=domain, convUnit = convUnit)
    transport = SimpleTransport(domain = domain, domain_area =domain_area, convUnit=convUnit)
    g,dus = disperse(g,emission, transport, 'septo3d', weather_data=weather_data)
    return g
 
from alinea.alep.protocol import external_contamination
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum, Septo3DSoilContamination
from alinea.septo3d.dispersion.alep_interfaces import SimpleSoilInoculum, SimpleContamination
from alinea.septo3d.cycle.alep_objects import Septo3D_DU
generator = Septo3D_DU()
 
def contamination(g, weather_data, level, domain, domain_area, convUnit):
    inoc = SoilInoculum(DU_generator=generator, sporulating_fraction =level, domain_area = domain_area, convUnit=convUnit)
    contaminator = Septo3DSoilContamination(domain = domain, domain_area =domain_area, convUnit=convUnit)
    g = external_contamination(g, inoc, contaminator, weather_data)
    return g
  
def simple_contamination(g, weather_data, level, domain, domain_area, convUnit):
    inoc = SimpleSoilInoculum(DU_generator=generator, sporulating_fraction =level, domain_area = domain_area, convUnit=convUnit)
    contaminator = SimpleContamination(domain = domain, domain_area =domain_area, convUnit=convUnit)
    g = external_contamination(g, inoc, contaminator, weather_data)
    return g  


from alinea.echap.interception_leaf import InterceptModel, pesticide_applications
from pesticide_data import products_data
interception = InterceptModel(products_data)

def pesticide_intercept(g, application_data, label='LeafElement'):
    return pesticide_interception(g, interception, application_data, label='LeafElement')
    
    
def application_data(product='Opus', appdate = '2011-04-19', dose = 0.5):
    applications= 'date,dose, product_name\n%s 10:00:00, %f, %s'%(appdate, dose, product)
    return pesticide_applications(applications)

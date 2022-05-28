# Imports
from alinea.echap.imports_echap import *

# Initiate
wheat = AdelWheat()
g,_ = new_canopy(wheat, age=100)
#grow_canopy(g, wheat, time_control)
# Début de la simul
t_deb = "2000-10-01 01:00:00"
# db pearl
db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}}
# Pesticide interception
pest_calendar = {'datetime':['2000-10-01 05:00:00', '2000-10-01 08:00:00', '2000-10-02 01:00:00'], 'dose':[1.5, 2, 1.2], 'product_name':['Opus', 'Banko 500', 'Opus']}

# Meteo microclimate
meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')

# Pesticide interception
pesticide_interception_model = CaribuInterceptModel(pest_calendar=pest_calendar)
# Efficacy
efficacy_model = PesticideEfficacyModel()
# Microclimate
climate_model = MicroclimateLeaf()
weather = Weather(data_file=meteo01_filepath)
# Rain interception
rain_interception_model = RapillyInterceptionModel()
# Septo3d
stock = create_stock(N=100,par=None)
inoculator = RandomInoculation()
dispersor = RandomDispersal()
# Pearl
Pearl_decay_model = PearLeafDecayModel(db)
# Milne
Milne_decay_model = PenetratedDecayModel()
# Fake controler in septo3d.cycle.alep_objects
controler = GrowthControlModel()

# Timer
nbsteps = 4
microclimate_timing = TimeControl(steps = nbsteps, weather = weather, model = climate_model, start_date = t_deb)
pest_timing = TimeControl(steps = nbsteps, model = pesticide_interception_model, start_date = t_deb)
meteo_timing = TimeControl(delay = 1, steps = nbsteps)
rain_timing = TimeControl(steps = nbsteps, weather = weather, model = rain_interception_model, start_date = t_deb)
septo_timing = TimeControl(delay = 1, steps = nbsteps)

timer = TimeControler(microclim = microclimate_timing, pest = pest_timing, meteo = meteo_timing, rain = rain_timing, septo = septo_timing)

# Initialisation de la maladie
g = initiate(g, stock, inoculator)
# healthy_surface
g = set_initial_properties_g(g, surface_leaf_element=5.)

for tc in timer:
    # Pesticide interception (timer)
    g,_ = pesticide_interception(g, pesticide_interception_model, tc['pest'], label='LeafElement')
    # Pesticide efficacy
    g = pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
    # Microclimate (timer)
    g,_ = local_microclimate(g, climate_model, tc['microclim'], label='LeafElement')
    print(g.property('microclimate'))
    # Rain interception (timer)
    g = rain_interception(g, rain_interception_model, tc['rain'], label='LeafElement', geometry = 'geometry')
    # Infect (timer)
    g = infect(g, tc['septo'].dt)
    # Update
    g = update(g, tc['septo'].dt, controler)
    # Disperse
    if tc['rain'].dt > 0 :
        g = disperse(g, dispersor, fungus_name='septo3d', label="LeafElement")
    # Surfacic decay
    g = pesticide_surfacic_decay(g, Pearl_decay_model, timestep=1)
    # Penetrated decay
    g = pesticide_penetrated_decay(g, Milne_decay_model, timestep=1)




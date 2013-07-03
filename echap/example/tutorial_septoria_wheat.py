# Imports
from alinea.echap.imports_echap import *

# Initiate
g = adel_mtg()
# Début de la simul
t_deb = "2000-10-01 01:00:00"

# Meteo microclimate
meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
# Pesticide calendar
treatment01_filepath = get_shared_data_path(['alinea/echap'], 'treatment01.csv')

# Pesticide interception
pesticide_interception_model = CaribuInterceptModel(pest_calendar=treatment01_filepath)
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
# Fake controler in septo3d.cycle.alep_objects
controler = GrowthControlModel()

# Timer
nbsteps = 2
pest_timing = TimeControl(steps = nbsteps, model = pesticide_interception_model, start_date = t_deb)
meteo_timing = TimeControl(delay = 1, steps = nbsteps)
rain_timing = TimeControl(steps = nbsteps, weather = weather, model = rain_interception_model, start_date = t_deb)
septo_timing = TimeControl(delay = 1, steps = nbsteps)

timer = TimeControler(pest = pest_timing, meteo = meteo_timing, rain = rain_timing, septo = septo_timing)

# Initialisation
g = initiate(g, stock, inoculator)
# healthy_surface
g = set_initial_properties_g(g, surface_leaf_element=5.)

for tc in timer:
    # Pesticide interception (timer)
    g = pesticide_interception(g, pesticide_interception_model, tc['pest'], label='LeafElement')
    # Pesticide efficacy
    g = pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
    # Microclimate (timer)
    g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
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




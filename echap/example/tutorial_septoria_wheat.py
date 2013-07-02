# Imports
from alinea.echap.imports_echap import *


# Initiate
g0 = adel_mtg()

# Septo3d
stock = create_stock(N=100,par=None)
inoculator = RandomInoculation()
g0 = initiate(g0, stock, inoculator)
# healthy_surface
g0 = set_initial_properties_g(g0, surface_leaf_element=5.)

# Meteo microclimate
meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
t_deb = "2000-10-01 01:00:00"
# Pesticide calendar
treatment01_filepath = get_shared_data_path(['alinea/echap'], 'treatment01.csv')
# Models
calendar = Calendar(data_file=treatment01_filepath)
pesticide_interception_model = CaribuInterceptModel()
climate_model = MicroclimateLeaf()
weather = Weather(data_file=meteo01_filepath)
rain_interception_model = RapillyInterceptionModel()

g0 = local_microclimate(g0, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]

# Fake controler in septo3d.cycle.alep_objects
controler = GrowthControlModel()

# Dispersal
dispersor = RandomDispersal()

g = copy.copy(g0)

# Timer
nbsteps = 4
pest_timing = TimeControl(steps = nbsteps, treat=treat, model = pesticide_interception_model, start_date = t_deb)
meteo_timing = TimeControl(delay = 1, steps = nbsteps)
wheat_timing = TimeControl(delay = 1, steps = nbsteps)
septo_timing = TimeControl(delay = 1, steps = nbsteps)
rain_timing = TimeControl(steps = nbsteps, weather=weather, model = rain_interception_model, start_date = t_deb)
plot_timing = TimeControl(delay=1, steps = nbsteps)

timer = TimeControler(rain = rain_timing, septo = septo_timing, ploting = plot_timing, wheat = wheat_timing, meteo = meteo_timing)

for tc in timer:
    g = infect(g, tc['septo'].dt)
    print count_lesions_by_leaf(g)
    print g.property('lesions')
    g = update(g, tc['septo'].dt, controler)
    g = rain_interception(g, rain_interception_model, tc['rain'], label='LeafElement', geometry = 'geometry')
    if tc['rain'].dt > 0 :
        g = disperse(g, dispersor, fungus_name='septo3d', label="LeafElement")




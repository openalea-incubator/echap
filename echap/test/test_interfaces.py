################### Imports
from alinea.echap.imports_echap import *

########################## tests decays

def test_surfacic():
    g = adel_mtg2()
    g = update_no_doses(g)
    db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}}
    decay_model = PearLeafDecayModel(db)
    g = pesticide_surfacic_decay(g, decay_model)
    return g

def test_penetrated():
    g = adel_mtg2()
    g = update_on_leaves(g)
    decay_model = PenetratedDecayModel()
    g = pesticide_penetrated_decay(g, decay_model)
    return g

def test_all_decay():
    g = adel_mtg()
    g = update_on_leaves(g)
    db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}}
    decay_model = PearLeafDecayModel(db)
    g = pesticide_surfacic_decay(g, decay_model)
    decay_model = PenetratedDecayModel()
    g = pesticide_penetrated_decay(g, decay_model)
    return g

def test_no_doses():
    g = adel_mtg()
    g = update_no_doses(g)
    db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}}
    decay_model = PearLeafDecayModel(db)
    g = pesticide_surfacic_decay(g, decay_model)
    decay_model = PenetratedDecayModel()
    g = pesticide_penetrated_decay(g, decay_model)
    return g

###################################### tests interception

def test_intercept():
    g = adel_mtg()
    g = update_on_leaves(g)
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    return g


def test_intercept_no_dose():
    g = adel_mtg()
    g = update_no_doses(g)
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    return g


def test_intercept_elevation():
    g = adel_mtg()
    g = update_no_doses(g)
    interception_model = CaribuInterceptModel(elevation=90.)
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    print g.property('surfacic_doses')
    interception_model = CaribuInterceptModel(elevation=45.)
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    print g.property('surfacic_doses')
    interception_model = CaribuInterceptModel(elevation=1.)
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    print g.property('surfacic_doses')
    return g

###################################### test microclimate

def test_microclimate():
    g = adel_mtg()
    g = update_no_doses(g)
    climate_model = MicroclimateLeaf()
    meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
    weather = Weather(data_file=meteo01_filepath)
    t_deb = "2000-10-01 01:00:00"
    g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
    print g.property('microclimate')
    return g


###################################### test efficacy

def test_efficacy():
    g = adel_mtg()
    g = update_no_doses(g)
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    efficacy_model = PesticideEfficacyModel()
    g = pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
    return g

def test_efficacy_nopest():
    g = adel_mtg()
    efficacy_model = PesticideEfficacyModel()
    g = pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
    return g

    
###################################### test lesions _ DU

def test_initiate_DU():
    """ Check if 'initiate' from 'alep.protocol.py' deposits DU on the MTG.
    """
    g = adel_mtg()
    stock = create_stock(N=100,par=None)
    inoculator = RandomInoculation()
    g = initiate(g, stock, inoculator)
    plot_DU(g)
    return g

def test_infect():
    """ Check if 'infect' from 'alep.protocol.py' leads to infection by dispersal units on the MTG.
    """
    g = adel_mtg()
    # Initiate
    stock = create_stock(N=100,par=None)
    inoculator = RandomInoculation()
    g = initiate(g, stock, inoculator)
    # Infect
    g = infect(g, dt=1)
    plot_lesions(g)
    return g

def test_update():
    """ Check if 'update' from 'alep.protocol.py' provokes the growth of a lesion instantiated on the MTG.
    """
    # Initiate and infect
    g = adel_mtg()
    g = update_no_doses(g)
    stock = create_stock(N=100,par=None)
    inoculator = RandomInoculation()
    g = initiate(g, stock, inoculator)
    g = infect(g, dt=1)
    # microclimate pour la methode update
    meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
    t_deb = "2000-10-01 01:00:00"
    climate_model = MicroclimateLeaf()
    weather = Weather(data_file=meteo01_filepath)
    g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
    # Efficacy pour la methode update
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    efficacy_model = PesticideEfficacyModel()
    g = pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
    # healthy_surface
    g = set_initial_properties_g(g, surface_leaf_element=5.)
    # Update Lesions
    growth_control_model = GrowthControlModel()
    g = update(g, dt=1, growth_control_model = growth_control_model)
    return g


def test_disperse():
    """ Check if 'disperse' from 'protocol.py' disperse new dispersal units on the MTG.
    """
    from alinea.echap.imports_echap import *
    # Initiate and infect
    g = adel_mtg()
    stock = create_stock(N=100,par=None)
    inoculator = RandomInoculation()
    g = initiate(g, stock, inoculator)
    g = infect(g, dt=1)
    # microclimate pour la methode update
    meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
    t_deb = "2000-10-01 01:00:00"
    climate_model = MicroclimateLeaf()
    weather = Weather(data_file=meteo01_filepath)
    g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
    global_climate = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[2]
    # Efficacy pour la methode update
    interception_model = CaribuInterceptModel()
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    efficacy_model = PesticideEfficacyModel()
    g = pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
    # healthy_surface
    g = set_initial_properties_g(g, surface_leaf_element=5.)
    # Update Lesions
    growth_control_model = GrowthControlModel()
    g = update(g, dt=1, growth_control_model = growth_control_model)
    # Rain interception model
    rain_interception_model = RapillyInterceptionModel()
    rain_timing = TimeControl(delay = 24, steps = 2, model = rain_interception_model, weather = weather)
    timer = TimeControler(rain = rain_timing)
    t = timer.next()
    g = rain_interception(g, rain_interception_model, t['rain'], label='LeafElement', geometry = 'geometry')
    # Disperse
    dispersor = RandomDispersal()
    g = disperse(g, dispersor, fungus_name='septo3d', label="LeafElement", activate=True)
    return g

les=g.property('lesions')[12][0]
les.AfficheLesion()
les.devLes(Svert=5, dt=1, T=[28], PPFD=[200], Rh=[90], efficacy={'protectant':0, 'eradicant' : 0})
les.AfficheLesion()


##################################### loop test

def test_decay_doses():
    db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'TemDif':20, 'FacWas':0.0001, 'FraDepRex':0}}
    t_deb = "2000-10-01 01:00:00"
    # Loop
    t = 0
    dt = 1
    nb_steps = 6
    # Initialisation du mtg 
    g = adel_mtg()
    g = update_no_doses(g)
    # models
    interception_model = CaribuInterceptModel()
    Pearl_decay_model = PearLeafDecayModel(db)
    Milne_decay_model = PenetratedDecayModel()
    climate_model = MicroclimateLeaf()
    meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
    weather = Weather(data_file=meteo01_filepath)
    # Interception
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
    t_deb = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[3]
    # sauvegarde etat initial
    out = get_df_out(0,g)
    # loop
    for i in range(nb_steps):
        t += dt        
        # Surfacic decay
        g = pesticide_surfacic_decay(g, Pearl_decay_model, timestep=dt)
        # Penetrated decay
        g = pesticide_penetrated_decay(g, Milne_decay_model, timestep=dt)
        df = get_df_out(t,g)
        out = out.append(df)
        plot_pesticide_norm(g, property_name='surfacic_doses', compound_name='Epoxiconazole', cmap=green_lightblue_blue, lognorm=False)
        #plot_pesticide(g, property_name='surfacic_doses', compound_name='Epoxiconazole', cmap=green_lightblue_blue)
        g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
        t_deb = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[3]
    return out


def test_local_meteo():
    # Loop
    t_deb = "2000-10-01 01:00:00"
    t = 0
    dt = 1
    nb_steps = 2
    # Initialisation du mtg 
    g = adel_mtg()
    g = update_no_doses(g)
    # models
    interception_model = CaribuInterceptModel()
    climate_model = MicroclimateLeaf()
    meteo01_filepath = get_shared_data_path(['alinea/echap'], 'meteo01.csv')
    weather = Weather(data_file=meteo01_filepath)
    # Interception
    g = pesticide_interception(g, interception_model, product_name='Opus', dose=1.5)
    print t_deb
    # Microclimate
    g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
    t_deb = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[3]
    print g.property('microclimate')
    # loop
    for i in range(nb_steps):
        t += dt    
        print t_deb
        # Microclimate
        g = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
        t_deb = local_microclimate(g, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[3]
        print g.property('microclimate')





   


    
    
    

"""Defines interfaces between g and the different models of echap project
"""

def setup_canopy(plant_model, age = 0):
    g = plant_model.setup_canopy(age)
    return g

def grow_canopy(g,plant_model,time_control):
    plant_model.grow(g,time_control)
    return g


def pesticide_interception(g, interception_model, product_name, dose, label='LeafElement'):
    """ 
    Interface between g and the interception model

    :Parameters:
    ----------
    - 'g' : MTG representing the canopy (and the soil)
    - 'interception_model' : A class embending the interception model and provide the following methods:    
        - 'interception_model.intercept(product_name, dose)' : Return the dictionnary of scene_id: compound name of the product and surfacic doses (g.m-2)
        See :class:`~alinea.echap.interception_leaf.CaribuInterceptModel`
    - 'product_name' (str)
        Commercial name of the product (e.g. "Opus", "Banko 500") 
    - 'dose' (float)
        Dose of product use in field (l.ha)
    - label (str)
        default "LeafElement"

    :Returns:  
    --------
    - 'g' : Updated MTG representing the canopy (and the soil). 'surfacic_doses' property is added to g or updated if present.
    
    :Example:
    -------
      >>> g = MTG()
      >>> interception_model = CaribuInterceptModel()
      >>> pesticide_interception(g, interception_model, product_name, dose)
      >>> return g
    """
    scene_geometry = g.property('geometry')
    surf_dose = interception_model.intercept(product_name, dose, scene_geometry)
    if not 'penetrated_doses' in g.properties():
        vi = [vid for vid in surf_dose if g.label(vid).startswith(label)]   
        g.add_property('penetrated_doses')
        g.property('penetrated_doses').update(dict((v, {'Chlorothalonil': 0, 'Epoxiconazole': 0}) for v in vi))
    if not 'surfacic_doses' in g.properties():
        vi = [vid for vid in surf_dose if g.label(vid).startswith(label)]   
        g.add_property('surfacic_doses')
        g.property('surfacic_doses').update(dict((v, {'Chlorothalonil': 0, 'Epoxiconazole': 0}) for v in vi))
    surfacic_doses = g.property('surfacic_doses')
    for vid, nd in surf_dose.iteritems():
        if g.label(vid).startswith(label):
            for name, dose in nd.iteritems():
                if name in surfacic_doses[vid]:
                    surfacic_doses[vid][name] += dose
                else:
                    g.add_property('surfacic_doses')
                    g.property('surfacic_doses')[vid].update(nd)
    return g


def local_microclimate(g, weather, climate_model, t_deb="2000-10-01 01:00:00", label='LeafElement', timestep=1):
    """ 
    Interface between g and the microclimate model

    :Parameters:
    -----------
    - 'g' : MTG representing the canopy (and the soil)
    - 'weather' : A class embending the weather reader and provide the following methods:    
        - weather.get_weather(timestep, t_deb)
        - weather.str_to_datetime(t_deb)
        - weather.add_global_radiation(globalclimate)
        - weather.add_vapor_pressure(globalclimate)
        - weather.next_date(timestep, t_deb)
        - weather.PPFD_to_global(PAR)
        - weather.Psat(T)
        - weather.humidity_to_vapor_pressure(humidity, Tair)
            See :class:`alinea.weather.global_weather.Weather`
    - 'climate_model' : A class embending the microclimate model and provide the following methods:    
        - 'climate_model.microclim(mean_globalclimate, scene)' : Return the dictionnary of scene_id: radiation and rain
        See :class:`~alinea.echap.microclimate_leaf.MicroclimateLeaf`
    - 't_deb' (str) format = "%Y-%m-%d %H:%M:%S"
        The start date to run the simulation. Default ""2000-10-01 01:00:00""
    - label (str)
        default "LeafElement"
    - timestep (int)
        The timestep of the simulation. Default 1

    :Returns:  
    --------
    - 'g' : Updated MTG representing the canopy (and the soil). 'microclimate' property (radiation and rain) is added to g or updated if present.
    - mean_globalclimate (Pandas dataframe)
        Climate variables averaged over the time step of the simulation
    - globalclimate (Pandas dataframe)
        Climate variables from t_deb to the end of the timestep
    - t_deb (str) format = "%Y-%m-%d %H:%M:%S"
        The new strat date (updated with the timestep) for the next simulation step

    :Example:
    ---------
      >>> g = MTG() 
      >>> weather = Weather()
      >>> climate_model = MicroclimateLeaf()
      >>> local_microclimate(g, weather, climate_model, t_deb, label='LeafElement', timestep)
      >>> return g, mean_globalclimate, globalclimate, t_deb
    """    scene_geometry = g.property('geometry')

    # Convert str into datetime
    t_deb = weather.str_to_datetime(t_deb)

    # Temporary : for tests
    mean_globalclimate, globalclimate = weather.get_weather(timestep, t_deb)
    local_meteo = climate_model.microclim(weather, scene_geometry, timestep, t_deb)
    g.add_property('microclimate')
    g.property('microclimate').update(local_meteo)

    # Update t_deb
    t_deb = weather.next_date(timestep, t_deb)

    return g, mean_globalclimate, globalclimate, t_deb


def pesticide_surfacic_decay(g, decay_model, label='LeafElement', timestep=1):
    """ 
    Interface between g and the decay model of Pearl

    :Parameters:
    -----------
    - 'g' : MTG representing the canopy (and the soil), doses are stored in the MTG as a property
    - 'decay_model' : A class embending the decay model and provide the following methods:    
        - 'decay_model.decay_and_penetrate(compound_name, compound_dose, microclimate, timestep)': Return for one compound the decayed surfacic dose (g.m-2), the penetrated amount and the loss to the environment.
        See :class:`alinea.pearl.pearl_leaf.PearLeafDecayModel`
    - label (str)
        default "LeafElement"
    - timestep (int)
        The timestep of the simulation. Default 1

    :Returns:    
    ----------
    - 'g' : Updated MTG representing the canopy (and the soil). 'surfacic_doses' and 'penetrated_doses' properties are added to g or updated if present.

    :Example:
    -------
      >>> g = MTG()
      >>> db = {'Chlorothalonil':{}, 'Epoxiconazole':{}}
      >>> interception_model = CaribuInterceptModel()
      >>> pesticide_interception(g, interception_model, product_name, dose)
      >>> decay_model = PearLeafDecayModel(db)
      >>> g = pesticide_surfacic_decay(g, decay_model)
      >>> return g
    """
    surfacic_doses = g.property('surfacic_doses') 
    microclimate = g.property('microclimate') 
    if not 'penetrated_doses' in g.property_names():
        g.add_property('penetrated_doses')
    penetrated_doses = g.property('penetrated_doses')
    for vid, d in surfacic_doses.iteritems():
        if g.label(vid).startswith(label):
            for compound_name,compound_dose in d.iteritems():
                new_dose,penetrated_amount,loss = decay_model.decay_and_penetrate(compound_name,compound_dose,microclimate[vid],timestep)
                surfacic_doses[vid][compound_name] = new_dose
                if vid in penetrated_doses:
                    if compound_name in penetrated_doses:
                        penetrated_doses[vid][compound_name] += penetrated_amount
                    else:
                        penetrated_doses[vid].update({compound_name:penetrated_amount})
                else:
                    penetrated_doses[vid] = {compound_name : penetrated_amount}
    return g


def pesticide_penetrated_decay(g, decay_model, label='LeafElement', timestep=1):
    """ 
    Interface between g and the decay model of penetrated doses of pesticide

    :Parameters:
    ----------
    - 'g' : MTG representing the canopy (and the soil), doses are stored in the MTG as a property
    - 'decay_model' : : A class embending the penetrated doses decay model and provide the following methods:    
        - 'decay_model.decay(d,timestep)' : Model of penetrated pesticide decay. Return for one compound the dictionnary of compound name of the productof and the decayed penetrated dose (g.m-2).
        See :class:`~alinea.echap.milne_leaf.PenetratedDecayModel`
        See :func:`alinea.simcycle.pesticide`
    - label (str)
        default "LeafElement"
    - timestep (int)
        The timestep of the simulation. Default 1

    :Returns:
    -------
    - 'g' : Updated MTG representing the canopy (and the soil). 'penetrated_doses' property is updated to g.

    :Example:
    -------   
      >>> g = MTG()
      >>> db = {'Chlorothalonil':{}, 'Epoxiconazole':{}, 'Metconazole':{}}
      >>> interception_model = CaribuInterceptModel()
      >>> pesticide_interception(g, interception_model, product_name, dose)
      >>> decay_model = PearLeafDecayModel(db)
      >>> g = pesticide_surfacic_decay(g, decay_model)
      >>> decay_model = PenetratedDecayModel()
      >>> g = pesticide_penetrated_decay(g, decay_model)
      >>> return g      
    """
    penetrated_doses = g.property('penetrated_doses')
    for vid, d in penetrated_doses.iteritems():
        if g.label(vid).startswith(label):
            penetrated_doses = decay_model.decay(d,timestep)
    return g


def pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1):
    """ 
    Interface between g and the efficacy model of pesticide compounds

    :Parameters:
    ----------
    - 'g' : MTG representing the canopy (and the soil), doses are stored in the MTG as a property
    - 'efficacy_model' : : A class embending the efficacy model of surfacic and penetrated doses of pesticide and provide the following methods:    
        - 'efficacy_model.efficacy(surfacic_doses, penetrated_doses, timestep)' : Model of pesticide efficacy. Return the dictionnary of protectant and eradicant efficacy (value between 0 and 1) respectively for surfacic doses and penetrated doses of pesticide.
        See :class:`alinea.pesticide_efficacy.pesticide_efficacy.PesticideEfficacyModel` 
        See :func:`alinea.simcycle.pesticide`
    - label (str)
        default "LeafElement"
    - timestep (int)
        The timestep of the simulation. Default 1

    :Returns:
    ---------
    - 'g' : Updated MTG representing the canopy (and the soil). 'global_efficacy' property is updated to g.

    :Example:
    ---------   
      >>> g = MTG()
      >>> interception_model = CaribuInterceptModel()
      >>> pesticide_interception(g, interception_model, product_name, dose)
      >>> efficacy_model = PesticideEfficacyModel()
      >>> pesticide_efficacy(g, efficacy_model, label='LeafElement', timestep=1)
      >>> return g      
    """
    if not 'surfacic_doses' in g.properties(): 
        vids = [vid for vid in g if g.label(vid).startswith(label)]
        for v in vids : 
            n = g.node(v)
            n.surfacic_doses = {'Chlorothalonil':0,'Epoxiconazole':0}
            n.penetrated_doses = {'Chlorothalonil':0,'Epoxiconazole':0}
    penetrated_doses = g.property('penetrated_doses')
    surfacic_doses = g.property('surfacic_doses') 
    if not 'global_efficacy' in g.properties():
        g.add_property('global_efficacy')
    vids = [vid for vid in surfacic_doses if g.label(vid).startswith(label)]
    for v in vids :     
        g.property('global_efficacy').update({v: efficacy_model.efficacy(surfacic_doses[v], penetrated_doses[v], timestep)})
    return g


def rain_interception(g, rain_interception_model, time_control, label='LeafElement', geometry = 'geometry'):
    """ 
    Interface between g and the rain interception model

    :Parameters:
    ------------
    - 'g': MTG representing the canopy (and the soil)
    - 'rain_interception_model' : : A class embending the rain interception model and provide the following methods:    
        - 'rain_interception_model.intercept(scene_geometry, time_control)' : Model of rain interception. return : fraction_runoff = f(normales, D_pdf, v_pdf), impacted_surface = f(pdf,intensities, duration, surfaces), splashed_droplets_per_square_meter = f(v_pdf, intensities)
        See :class:`alinea.popdrops.rain.RainInterceptionModel` 
    - time_control
    - label (str)
        default "LeafElement"
    - geometry (str)
        default "geometry"

    :Returns:
    ---------
    returns the fraction of intercepted rain
    
    :Example:
    ---------   
    """
    

    if time_control.dt > 0:
        scene_geometry = g.property('geometry')
        if not 'rain' in g.properties():
            g.add_property('rain')
        rain_leaf, rain_fate = rain_interception_model.intercept(scene_geometry, time_control)
        rain = dict([(k,{'water':v}) for k,v in rain_leaf.iteritems()])
        g.property('rain').update(rain_fate)
        vids = [vid for vid in g.property('rain')]
        for v in vids : 
            g.property('rain')[v].update(rain[v])
    return g








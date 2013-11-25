"""Defines interfaces between g and the different models of echap project
"""

class EchapInterfacesError(Exception): pass


def pesticide_interception(g, interception_model, application_data, label='LeafElement'):
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
    - application data : a panda dataframe with (date), product_name and dose (l.ha-1) in columns and time-indexed
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
    dose = application_data[['dose']].values[0][0]
    if dose > 0:
        scene_geometry = g.property('geometry')
        product_name = application_data[['product_name']].values[0][0]
        
        surf_dose = interception_model.intercept(scene_geometry, product_name, dose)
        if not 'penetrated_doses' in g.properties():
            g.add_property('penetrated_doses')
        if not 'surfacic_doses' in g.properties():
            g.add_property('surfacic_doses')
        surfacic_doses = g.property('surfacic_doses')
        penetrated_doses = g.property('penetrated_doses')
        for vid, nd in surf_dose.iteritems():
            if g.label(vid).startswith(label):
                for name, dose in nd.iteritems():
                    if vid in surfacic_doses:
                        if name in surfacic_doses[vid]:
                            surfacic_doses[vid][name] += dose
                        else:
                            surfacic_doses[vid][name] = nd[name]
                            penetrated_doses[vid][name] = 0
                    else:
                        surfacic_doses[vid] = nd
                        penetrated_doses[vid] = dict([(n,0) for n in nd])
    return g, interception_model


def local_microclimate(g, climate_model, weather_data, label='LeafElement'):
    """ 
    Interface between g and the microclimate model

    :Parameters:
    -----------
    - 'g' : MTG representing the canopy (and the soil)
    - 'climate_model' : A class embending the microclimate model and provide the following methods:    
        - 'climate_model.microclim(mean_globalclimate, scene)' : Return the dictionnary of scene_id: radiation and rain
        See :class:`~alinea.echap.microclimate_leaf.MicroclimateLeaf`
    - 'weather_data' : A panda dataframe with climatic data. 
    - label (str)
        default "LeafElement"


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
    """
     scene_geometry = g.property('geometry')
    local_meteo = climate_model.microclim(scene_geometry, weather_data)
    g.add_property('microclimate')
    g.property('microclimate').update(local_meteo)

    return g, climate_model


def pesticide_surfacic_decay(g, decay_model, weather_data, label='LeafElement'):
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
    if 'surfacic_doses' in g.property_names():
        surfacic_doses = g.property('surfacic_doses')
        if 'microclimate' not in g.property_names():
            print 'could not compute decay : TODO : use weather data instead'
        else:
            timestep = len(weather_data)
            microclimate = g.property('microclimate') 
            if not 'penetrated_doses' in g.property_names():
                g.add_property('penetrated_doses')
            penetrated_doses = g.property('penetrated_doses')
            for vid, d in surfacic_doses.iteritems():
                if g.label(vid).startswith(label):
                    for compound_name,compound_dose in d.iteritems():
                        if vid in microclimate:
                            new_dose,penetrated_amount,loss = decay_model.decay_and_penetrate(compound_name,compound_dose,microclimate[vid],timestep)
                        else:
                            new_dose,penetrated_amount,loss = 0, 0, 0
                            #print 'EchapInterfacesError : KeyErreur', vid
                        surfacic_doses[vid][compound_name] = new_dose
                        if vid in penetrated_doses:
                            if compound_name in penetrated_doses:
                                penetrated_doses[vid][compound_name] += penetrated_amount
                            else:
                                penetrated_doses[vid].update({compound_name:penetrated_amount})
                        else:
                            penetrated_doses[vid] = {compound_name : penetrated_amount}
    return g


def pesticide_penetrated_decay(g, decay_model, weather_data, label='LeafElement'):
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
    if 'penetrated_doses' in g.property_names():
        penetrated_doses = g.property('penetrated_doses')
        for vid, d in penetrated_doses.iteritems():
            if g.label(vid).startswith(label):
                penetrated_doses = decay_model.decay(d,weather_data)
    return g


def pesticide_efficacy(g, efficacy_model, weather_data, label='LeafElement'):
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
    
    if 'surfacic_doses' in g.property_names():
        surfacic_doses = g.property('surfacic_doses') 
    else: 
        surfacic_doses = {}
    if 'penetrated_doses' in g.property_names():
        penetrated_doses = g.property('penetrated_doses')
    else:
        penetrated_doses = {}
        
    if (len(surfacic_doses) + len(penetrated_doses)) > 0:
        if not 'global_efficacy' in g.properties():
            g.add_property('global_efficacy')
        vids = set().union(surfacic_doses, penetrated_doses)  
        for v in vids :
            if g.label(v).startswith(label):
                g.property('global_efficacy').update({v: efficacy_model.efficacy(surfacic_doses.get(v,{}), penetrated_doses.get(v,{}))})
    else:
        if 'global_efficacy' in g.properties():
            g.remove_property('global_efficacy')
    return g


def rain_interception(g, rain_interception_model, weather_data, label='LeafElement'):
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


    :Returns:
    ---------
    returns the fraction of intercepted rain

    :Example:
    ---------   
    """

    rain_data = weather_data[['rain']]
    if rain_data.sum() > 0:
        intensity = rain_data.mean()
        scene_geometry = g.property('geometry')
        if not 'rain' in g.properties():
            g.add_property('rain')
        rain_leaf, rain_fate = rain_interception_model.intercept(scene_geometry, intensity)
        rain = dict([(k,{'water':v}) for k,v in rain_leaf.iteritems()])
        g.property('rain').update(rain_fate)
        #vids = [vid for vid in g.property('rain')]
        #for v in vids : 
        #    g.property('rain')[v].update(rain[v])
    return g


def record(g, weather_data, recorder, label = 'LeafElement'):
    """
    tentative protocol for recording data during a simulation
    """
    date = weather_data.index[0].to_datetime()
    for vid in g:
        if g.label(vid).startswith(label):
            n = g.node(vid)
            header = {'date' : date,
                      'plant' : n.complex().complex().complex().complex().label,
                      'axe' : n.complex().complex().complex().label,
                      'metamer' : n.complex().complex().label,
                      'organ' : n.complex().label,
                      'id' : n._vid
                         }
            recorder.record(n, header)
            
    return g

def save_records(recorder, path):
    recorder.save_records(path)
    return path

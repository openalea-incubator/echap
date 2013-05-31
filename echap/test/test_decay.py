from alinea.echap.wheat_mtg import *

########################## extraction d'un data frame des surfacic_doses au cours du temps de pesticide sur le mtg

def get_df_out(time,g):
    sd = g.property('surfacic_doses')
    lab = g.property('label')
    sp = g.property('penetrated_doses')            
    recs = [(id,lab[id],comp,dose) for id,v in sd.iteritems() for comp, dose in v.iteritems() if lab[id] is 'LeafElement']
    ids,labs,compounds,surfdose = zip(*recs)
    dfs = DataFrame({'time':[time]*len(ids), 'id' : ids, 'label': labs,'compound' : compounds, 'surfacic_doses' : surfdose})
    if not 'penetrated_doses' in g.property_names():
        df=dfs        
        #dfp = DataFrame(columns=('time', 'id', 'label','compound', 'penetrated_dose'))
    else:
        recp = [(id,lab[id],comp,dose) for id,v in sp.iteritems() for comp, dose in v.iteritems() if lab[id] is 'LeafElement']
        idp,labp,compoundp,pendose = zip(*recp)
        dfp = DataFrame({'time':[time]*len(idp), 'id' : idp, 'label': labp,'compound' : compoundp, 'penetrated_doses' : pendose})
        df = merge(dfs, dfp, left_on=('compound', 'id', 'label', 'time'), right_on=('compound', 'id', 'label', 'time'), how='outer')    
    return df


########################## display

def plot_decay(out, leaf=12):
    #from pylab import *
    df = out[out['id']==leaf]
    plt.plot(df['time'], df['surfacic_doses'])
    plt.plot(df['time'], df['penetrated_doses'])
    plt.show()


def plot_pesticide(g, property_name='surfacic_doses', compound_name='Epoxiconazole', cmap=green_lightblue_blue):
    """ plot the plant with pesticide doses """
    if type(cmap) is str:
        try:
            _cmap = cm.get_cmap(cmap())
        except:
            raise Exception('This colormap does not exist')
    else:
        _cmap = cmap()

    green = (0,180,0)

    for v in g.vertices(scale=g.max_scale()): 
        n = g.node(v)
        if 'surfacic_doses' in n.properties():
            r,gg,b,s = _cmap(n.surfacic_doses[compound_name]*200)
            n.color = (int(r*255),int(gg*255),int(b*255))           
        else : 
            n.color = green
    scene = plot3d(g)
    Viewer.display(scene)
    return g


##################################### loop test interception

def test_decay_doses():
    db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'FacTraDepRex':0, 
    'FacVolDepRex':0,'FacPenDepRex':0,'FacWasDepRex':0,'FraDepRex':0}, 
    'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
    'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'FacTraDepRex':0, 
    'FacVolDepRex':0,'FacPenDepRex':0,'FacWasDepRex':0,'FraDepRex':0}}
    t_deb=datetime(2000, 10, 01, 9, 00, 00)
    # Loop
    t = 0
    dt = 1
    nb_steps = 1
    # Initialisation des mtg 
    g90 = adel_mtg()
    g45 = adel_mtg()   
    g1 = adel_mtg()   
    # models
    interception_model_90 = CaribuInterceptModel(elevation=90.)
    interception_model_45 = CaribuInterceptModel(elevation=45.)
    Pearl_decay_model = PearLeafDecayModel(db)
    Milne_decay_model = PenetratedDecayModel()
    climate_model = MicroclimateLeaf()
    weather = Weather()
    # Interception
    g90 = pesticide_interception(g.copy, scene, interception_model_90, product_name='Opus new', dose=1.5)
    print g90.property('surfacic_doses')
    g90 = local_microclimate(g90, scene, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]

    g45 = pesticide_interception(g.copy, scene, interception_model_45, product_name='Opus new', dose=1.5)
    print g45.property('surfacic_doses')
    g45 = local_microclimate(g45, scene, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]

    t_deb = local_microclimate(g90, scene, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[3]

    # sauvegarde etat initial
    out90 = get_df_out(0,g90)
    print out90
    out45 = get_df_out(0,g45)
    print out45

    # loop
    for i in range(nb_steps):
        t += dt        

        # Surfacic decay
        g90 = pesticide_surfacic_decay(g90, Pearl_decay_model, timestep=dt)
        g45 = pesticide_surfacic_decay(g45, Pearl_decay_model, timestep=dt)

        # Penetrated decay
        g90 = pesticide_penetrated_decay(g90, Milne_decay_model, timestep=dt)
        df90 = get_df_out(t,g90)
        out90 = out90.append(df90)
        print out90

        g45 = pesticide_penetrated_decay(g45, Milne_decay_model, timestep=dt)
        df45 = get_df_out(t,g45)
        out45 = out45.append(df45)
        print out45

        g90 = local_microclimate(g90, scene, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]
        g45 = local_microclimate(g45, scene, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]

        t_deb = local_microclimate(g90, scene, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[3]

    return out90, out45



out90, out45 = test_decay_doses()
out45[(out45['compound']=='Epoxiconazole') & (out45['id']==12)]
out90[(out90['compound']=='Epoxiconazole') & (out90['id']==12)]





import copy

db = {'Epoxiconazole':{'MolarMass':329.8, 'VapourPressure':0.00001, 'TemVapPre':20, 'WatSolubility':7.1, 'TemWatSol':20, 
'ThicknessLay':0.0006,'DT50Pen':0.33,'DT50Tra':0.433, 'CofDifAir':0.428, 'FacTraDepRex':0, 
'FacVolDepRex':0,'FacPenDepRex':0,'FacWasDepRex':0,'FraDepRex':0}, 
'Chlorothalonil':{'MolarMass':265.9, 'VapourPressure':0.0000762, 'TemVapPre':20, 'WatSolubility':0.81, 'TemWatSol':20, 
'ThicknessLay':0.0006,'DT50Pen':0.14,'DT50Tra':0.23, 'CofDifAir':0.428, 'FacTraDepRex':0, 
'FacVolDepRex':0,'FacPenDepRex':0,'FacWasDepRex':0,'FraDepRex':0}}
t_deb=datetime(2000, 10, 01, 9, 00, 00)
# Loop
t = 0
dt = 1
nb_steps = 1

# Initialisation
g = adel_mtg()
product_name = 'Opus'
dose = 1.5
label='LeafElement'
scene_geometry = g.property('geometry')    

# models
interception_model_90 = CaribuInterceptModel(elevation=90.)
interception_model_45 = CaribuInterceptModel(elevation=45.)
Pearl_decay_model = PearLeafDecayModel(db)
Milne_decay_model = PenetratedDecayModel()
climate_model = MicroclimateLeaf()
weather = Weather()
# Interception
g90 = copy.copy(g)
scene90 = plot3d(g90)
g45 = copy.copy(g)
scene45 = plot3d(g45)
g90.property = copy.copy(g.property)
print g90.property('surfacic_doses')
print g45.property('surfacic_doses')
print g.property('surfacic_doses')

interception_model = CaribuInterceptModel()
doses = interception_model.intercept(product_name, dose, scene_geometry)

g90 = pesticide_interception(g90, scene90, interception_model_90, product_name='Opus new', dose=1.5)
print g90.property('surfacic_doses')
print g45.property('surfacic_doses')
print g.property('surfacic_doses')

id(g90.property('surfacic_doses'))

id(g45.property('surfacic_doses'))

id(g.property('surfacic_doses'))

g90 = local_microclimate(g90, scene90, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]

g45 = pesticide_interception(g45, scene45, interception_model_45, product_name='Opus new', dose=1.5)
print g90.property('surfacic_doses')
print g45.property('surfacic_doses')
print g.property('surfacic_doses')


g45 = local_microclimate(g45, scene45, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[0]

t_deb = local_microclimate(g90, scene90, weather, climate_model, t_deb=t_deb, label='LeafElement', timestep=1)[3]


# sauvegarde etat initial
out90 = get_df_out(0,g90)
print out90
out45 = get_df_out(0,g45)
print out45
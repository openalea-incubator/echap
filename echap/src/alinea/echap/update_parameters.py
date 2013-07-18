from alinea.echap.wheat_mtg import *
from pandas import *


#################### Update

def update_on_leaves(g, label = 'LeafElement'):
    """ Read weather data for a step of simulation and apply it to each leaf.
    """        
    vids = [vid for vid in g if g.label(vid).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.surfacic_doses={'Chlorothalonil':1,'Epoxiconazole':2}
        n.penetrated_doses={'Chlorothalonil':1,'Epoxiconazole':2}
        n.temp = 12
        n.rain_intensity = 0
        n.relative_humidity = 100 
        n.wetness = True
        n.microclimate = {}
    return g


########################## extraction d'un data frame des doses de pesticide sur le mtg

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


########################## import .csv as dict

def products_from_csv(csvname, delimiter = ';'):
    """ 
    Read a csv of products parameters and import them in a dict.
    Expected columns are :
        - 'product' : commercial name of the product
        - 'compound' : name of the active compound of the product
        - 'dose' : dose of active compound in the product (g.l-1)    
    """
    tab = recfromcsv(csvname, delimiter = delimiter, case_sensitive = True)
    d = {}
    for i in range(0,len(tab['compound'])):
        d[tab[i][0]] = {tab[i][1]:tab[i][2]}
    return d
    
    
############## Add Lesions on g

def lesions_infect():
    g = adel_mtg()
    stock = create_stock(N=10,par=None)
    inoculator = RandomInoculation()
    g = initiate(g, stock, inoculator)
    g = infect(g, dt=1)
    return g


def update_plot(g):
    # Count lesions by id & add it as MTG property ####################################
    #nb_lesions_by_leaf = count_lesions_by_leaf(g, label = 'lf')
    #set_property_on_each_id(g, 'nb_lesions', nb_lesions_by_leaf, label = 'lf')
                       
    # Visualization ###################################################################
    g = alep_colormap(g, 'nb_lesions', cmap=green_white(levels=10), lognorm=False)
    trunk_ids = [n for n in g if g.label(n).startswith('tronc')]
    brown = (100,70,30)
    for id in trunk_ids:
        trunk = g.node(id)
        trunk.color = brown
    scene = plot3d(g)
    Viewer.display(scene)
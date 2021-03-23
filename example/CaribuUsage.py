from alinea.adel.newmtg import *
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
import alinea.adel.fitting as fitting
from alinea.adel.AdelR import devCsv,setAdel,RunAdel,genGeoLeaf,genGeoAxe

from alinea.caribu.CaribuScene import CaribuScene

############# Initialise le mtg
def adelR(nplants,dd):
    devT = devCsv('../../adel/example/data/axeTCa0N.csv','../../adel/example/data/dimTCa0N.csv','../../adel/example/data/phenTCa0N.csv','../../adel/example/data/earTCa0N.csv','../../adel/example/data/ssi2sen.csv')
    geoLeaf = genGeoLeaf()
    geoAxe = genGeoAxe()
    pars = setAdel(devT,geoLeaf,geoAxe,nplants)
    cantable = RunAdel(dd,pars)
    return pars,cantable

def leaves_db():
    import pickle as Pickle
    fn = 'E:/openaleapkg/adel/adel/data/leaves_simple.db'
    f = open(fn)
    leaves = Pickle.load(f)
    f.close()
    leaves,discard = fitting.fit_leaves(leaves, 9)
    return leaves

def leaves_db_flow(fn):
    import pickle as Pickle
    f = open(fn)
    leaves = Pickle.load(f)
    f.close()
    leaves,discard = fitting.fit_leaves(leaves, 9)
    return leaves

def adel_mtg():
    """ create a very simple adel mtg """
    d = {'plant':[1,1],'axe_id':['MS','T1'],'ms_insertion':[0,1],'numphy':[1,1], 
         'Laz': [0,90], 'Ll' :[3,3], 'Lv' :[3,3] , 'Lsen':[0,0], 'L_shape':[3,3], 'Lw_shape':[.3,.3], 'Linc':[0,0],
         'Einc':[0,45],'El':[1,1],'Ev':[1,1],'Esen':[0,0],'Ed': [0.1,0.1]}
    g=mtg_factory(d,adel_metamer,leaf_db=leaves_db(), leaf_sectors=1)
    g=mtg_interpreter(g)
    return g

def adel_mtg2(nb_sect=1):
    """ create a less simple adel mtg """
    p, d = adelR(3,1000)
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,leaf_db=leaves_db(),stand=[((0,0,0),0),((10,0,0),90), ((0,10,0), 0)])
    g=mtg_interpreter(g)
    return g

def adel_mtg3(nb_sect=1, leaf_db=None, d=None, p=None):
    """ create a less simple adel mtg """
    if p: # nb_plants
        size = int(ceil(sqrt(p)))
        stand = numpy.array([(i, j) for i in range(size) for j in range(size)])
        numpy.random.shuffle(stand)
        stand = [((i, j, 0),random.randint(0,90)) for i, j in stand[:p]]
    else:
        stand = [((0,0,0),0),((10,0,0),90), ((0,10,0), 0)]
    g=mtg_factory(d,adel_metamer, leaf_sectors=nb_sect,leaf_db=leaf_db,stand=stand)
    g=mtg_interpreter(g)
    return g

def update_on_leaves(g, label = 'LeafElement'):
    """ Read weather data for a step of simulation and apply it to each leaf.
    """        
    vids = [vid for vid in g if g.label(vid).startswith(label)]
    for v in vids : 
        n = g.node(v)
        n.surfacic_doses={'Chlorothalonil':1,'Epoxiconazole':2,'Metconazole':3}
        n.penetrated_doses={'Chlorothalonil':1,'Epoxiconazole':2,'Metconazole':3}
        n.temp = 12
        n.rain_intensity = 0
        n.relative_humidity = 100 
        n.wetness = True
    return g
#######################

# Definition du g
g = adel_mtg()
g = update_on_leaves(g)


#scene from a g
g=mtg_interpreter(g) #mtginterpreter calcule la geometrie
scene = plot3d(g) # on produit la scene

########################################################### Caribu interception
# source emission = lumiere intensite 1, oriente selon vecteur (0,0 -1)
source = (1,(0,0,-1))  
c_scene = CaribuScene()    
idmap = c_scene.add_Shapes(scene)    
c_scene.addSources(source)
output = c_scene.runCaribu(infinity=False)
resultat = c_scene.output_by_id(output, idmap)['Einc']


################################################################## Caribu light
# source emission = lumiere intensite 1, oriente selon vecteur (0,0 -1)
import alinea.caribu.sky_tools.turtle as turtle
# Pour rayonnement diffu
source = turtle
#sectors='16'
energy, emission, direction, elevation, azimuth = turtle.turtle(sectors='46', energy=1) 
# Tout est fait avec 1 en energie de départ
sources = list(zip(energy, direction))
# Caribu
c_scene = CaribuScene()    
idmap = c_scene.add_Shapes(scene)    
c_scene.addSources(sources)
output = c_scene.runCaribu(infinity=False)

# pour la pluie -> einc, et pour un rayonnement -> einc / area
# résultat à mettre sur g en écrasant l'ancienne valeur (dans local_meteo)

microclimate_A = g.property('microclimate_A')
Einc = c_scene.output_by_id(output, idmap)['Einc']
Area = c_scene.output_by_id(output, idmap)['Area']
for eid, e in Einc.items():
    for aid, a in Area.items():
        if eid == aid:
            microclimate_A[eid] = {'radiation': e / a} 

microclimate_E = g.property('microclimate_E')
EiInf = c_scene.output_by_id(output, idmap)['EiInf']
EiSup = c_scene.output_by_id(output, idmap)['EiSup']
for Infid, e in EiInf.items():
    for Supid, a in EiSup.items():
        if Infid == Supid:
            microclimate_E[Infid] = {'radiation': e + a} 







# General useful imports
import numpy as np
import pandas
import matplotlib.pyplot as plt
plt.ion()

# Imports for wheat
import alinea.adel.data_samples as adel_data
from alinea.alep.architecture import set_properties

# Imports alep
from alinea.alep.protocol import *
from alinea.alep.disease_outputs import SeptoRecorder


# Septo 3D
from alinea.septo3d.cycle.alep_objects import Septo3DFungus
s3d_fungus = Septo3DFungus()

#pesticide
from alinea.echap.milne_leaf import PenetratedDecayModel
from alinea.pesticide_efficacy.pesticide_efficacy import PesticideEfficacyModel, as_doses
from alinea.echap.interfaces import pesticide_surfacic_decay, pesticide_penetrated_decay, pesticide_efficacy

#doses en gMA/ha
doses = {'Te':0,'A':13.75,'B':27.5,'C':55,'D':110,'E':220}

senescence = {'J':np.array([15,17,21,24,31]),
              'psen':np.array([0,0.01,0.106,0.709,1])}
      

def as_doses(doses, name='Chlorothalonil'):
    return [{name:d/100} for d in doses]

def example_one_leaf(fungus=s3d_fungus, growth_controler = None, infection_controler = None, nb_dus=10, nb_steps=1000, leaf_area=7.5*.5, initial_dose = 'A', mult_dose= 1, delta=7):

    wdata = {'rain':[0.], 'temperature_air':[20.],'relative_humidity':[100.], 'PPFD':[0.]}
    # Initiate wheat MTG with one leaf and good conditions for disease
    g = adel_data.adel_one_leaf()
    if leaf_area is not None:
        set_properties(g, label='LeafElement', area=leaf_area, green_area=leaf_area, senesced_area = 0)
    if fungus.name == 'septo3d':
        set_properties(g, label='LeafElement', microclimate = wdata)
    else:
        set_properties(g, label='LeafElement', wetness=True, temp=22., rain_intensity=0.)
        
    #dose en g/m2
    dose = doses[initial_dose] / 1e4
    set_properties(g, label='LeafElement', penetrated_doses = {'Chlorothalonil':dose*mult_dose}, surfacic_doses = {'Chlorothalonil':dose*mult_dose})


    surfacic_decay_model = PenetratedDecayModel()
#penetrated_decay_model = PenetratedDecayModel(pest_db.milne_parameters)
    penetrated_decay_model = PenetratedDecayModel()
    efficacy_model = PesticideEfficacyModel()  

#senescence
    #hypothetic days after senescence when 'killing' senescence ocurs
    sene = np.interp(range(nb_steps),(senescence['J']+delta)*24,senescence['psen']*leaf_area)

    
    # Initiate disease111
    leaf = g.node(10)
    dispersal_units=[fungus.dispersal_unit() for i in range(nb_dus)]
    for du in dispersal_units:
        du.set_position([random.random() * leaf.length, 0])
    leaf.dispersal_units = dispersal_units

    # Initiate output recorder
    recorder = SeptoRecorder(vids=[10], group_dus=False)
    date_sequence=pandas.date_range("2014-01-01 01:00:00", periods=nb_steps+1, freq='H')
    

    
    # Simulation loop
    for i in range(nb_steps):
        # Update wheat healthy area
        #update_healthy_area(g, label = 'LeafElement')
        
        leaf.green_area = leaf_area - sene[i]
        leaf.senesced_area = sene[i]
        # Update dispersal units and lesions
        
        pesticide_efficacy(g, efficacy_model, wdata)
        infect(g, 1, infection_controler, label='LeafElement')
        update(g, 1, growth_controler, senescence_model=None, label='LeafElement')
        
        # Get outputs
        recorder.record(g, date_sequence[i])
    recorder.get_complete_dataframe()
   
    return g, recorder
    
    
def plot_lesion_states(rec, ax=None, **kwds):
    outputs = pandas.DataFrame({#'disease area': rec.data.leaf_disease_area / rec.data.leaf_area,
                                #'surface incubating': rec.data.surface_inc / rec.data.leaf_area,
                                'surface chlorotic': rec.data.surface_chlo / rec.data.leaf_area*100,
                                'surface necrotic': (rec.data.leaf_senesced_area+rec.data.surface_spo) / rec.data.leaf_area*100,
                                'surface sporulating': rec.data.surface_spo / rec.data.leaf_area*100})
    ax = outputs.plot(ax=ax, **kwds)
    ax.set_xlabel('Days')
    ax.set_ylabel('%Surface')
    ax.set_ylim([0,100])
    return ax
    

def run_sim(dus = 10, dose='A', delta = 7, mult_dose=1):  
    g,res=example_one_leaf(nb_dus=dus, nb_steps=1000, initial_dose = dose, delta=delta, mult_dose=mult_dose)
    orange = (1,102./255,0)
    maron=(153./255,51./255,0)
    bleu=(79./255,129./255,189./255)
    ax=plot_lesion_states(res, style='-', color=[orange, maron,bleu], linewidth=2)
    return g,res, ax
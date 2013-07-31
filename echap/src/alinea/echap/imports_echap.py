############################ Imports
# Utils
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize, LogNorm
from pandas import *
from datetime import datetime, timedelta
from openalea.deploy.shared_data import get_shared_data_path

# Utils echap - alep
#from alinea.echap.update_parameters import *
from alinea.echap.tests_nodes import *
from alinea.alep.disease_outputs import *
from alinea.astk.plantgl_utils import *

# Pearl
from alinea.pearl.pearl import *
from alinea.pearl.pearl_leaf import *
from alinea.echap.interfaces import pesticide_surfacic_decay
# Milne
from alinea.echap.milne_leaf import *
from alinea.echap.interfaces import pesticide_penetrated_decay
# Interception
from alinea.echap.interception_leaf import *
from alinea.echap.interfaces import pesticide_interception
# Microclimate
from alinea.weather.global_weather import *
from alinea.echap.microclimate_leaf import *
from alinea.echap.interfaces import local_microclimate
# Efficacy
from alinea.pesticide_efficacy.pesticide_efficacy import *
from alinea.echap.interfaces import pesticide_efficacy
# Rain interception
from alinea.septo3d.Rapilly import *
from alinea.echap.interfaces import rain_interception
# Septo3d
from alinea.septo3d.cycle.alep_objects import *
# Popdrops
from alinea.popdrops.PesticideInterception import *
from alinea.popdrops.RainInterception import *
from alinea.popdrops.Rain import *
from alinea.popdrops.Drops import *

# Alep protocol
from alinea.alep.inoculation import RandomInoculation
from alinea.alep.dispersal import RandomDispersal
from alinea.alep.protocol import initiate
from alinea.alep.protocol import infect
from alinea.alep.protocol import update
from alinea.alep.protocol import disperse

# Wheat
#from alinea.echap.wheat_mtg import *
from alinea.adel.astk_interface import AdelWheat
from alinea.astk.plant_interface import *

# Time control
from alinea.astk.TimeControl import *
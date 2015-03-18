""" Generate files of simulated canopy properties for varied wheat varieties
    with multiprocessing """

from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from evaluation_canopy_properties import *

def run_variety(variety = 'Tremie12'):
    run_and_save_multi_simu(variety = variety, nplants = 30, age_range = [400, 2600],
                            time_steps = [20, 100], nrep = 1)

# TODO : Set cmdlines that manage several jobs on One open Povray
# if __name__ == '__main__':
    # varieties = ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']
    # nb_cpu = cpu_count()
    # pymap(run_variety, varieties, nb_cpu)
    
for variety in ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']:
    run_variety(variety)
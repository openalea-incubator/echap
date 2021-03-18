# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:16:51 2013

@author: lepse
"""

from numpy import exp
import pandas

from alinea.echap.pesticide_data import milne_parameters

def _dose_decay(decay_rate, initial_dose, days):
    active_dose = initial_dose * exp(-decay_rate * days) 
    return active_dose


class SimcyclePesticide(Exception): pass       


def milne_leaf(initial_dose, compound_parameters, hours):
    """ Milne decay model applied to penetrated doses

        :Parameters:
        ------------
        - `initial_dose` (float)
            Mass of product present per unit surface inside the leaf (g.m-2)
        - `compound_parameters` (dict)
            A dict of compound parameters read in a .csv file 
        - `hours` (int)
            Timestep (hours)

        :Returns:
        ---------
        - 'initial_dose' (float)
            Remaining mass of product per unit surface inside the leaf after dt (g.m-2) 
    """
    days = float(hours) / 24
    ptable = compound_parameters
    for name, dose in initial_dose.iteritems():
        try:
            initial_dose[name] = _dose_decay(ptable[name]['decay_rate'], initial_dose[name], days)
        except KeyError:
            raise SimcyclePesticide('product %s not found in parameter dict'%name)
    return initial_dose

    
class MilneModel(object):
    """ A simple user-interface for manipulating Milne model for a single compound    
    """
    def __init__(self, compound_parameters):
        self.compound_parameters = compound_parameters
        
    def decay(self, dose, time_sequence):
        hours = range(len(time_sequence))
        doses = map(lambda x: _dose_decay(self.compound_parameters['decay_rate'], dose, float(x) / 24), hours)
        res = pandas.DataFrame({'datetime' : time_sequence, 'dose' :doses})
        res.index = res['datetime']
        return res

class PenetratedDecayModel(object):
    """ Adaptor for Milne decay model and compliying echap pesticide_penetrated_decay model protocol
    """
    def __init__(self,compound_parameters = milne_parameters):
        self.compound_parameters = compound_parameters
        
    def decay(self, initial_dose, weather_data):
        """ Return the penetrated active doses 

        :Parameters:
        ----------
        - 'initial_dose' (float) 
            Total amount of penetrated dose in the leaf element
        - 'dt' (int) 
            Timestep of the simulation (day)

        :Returns:
        -------
        - 'active_dose' (float)
            Total amount of penetrated active dose in the leaf element after decay
        """    
        hours = len(weather_data)
        active_dose = milne_leaf(initial_dose, self.compound_parameters, hours)
        return active_dose
    def decay_and_penetrate(self, name, dose, microclimate, dt):
        """ make all the product penetrates to simulate fully the model of Milne"""
        newdose = 0
        loss = 0
        penetrated = dose
        return newdose,penetrated,loss


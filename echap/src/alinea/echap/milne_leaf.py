# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:16:51 2013

@author: lepse
"""
from numpy import recfromcsv
from numpy import exp


def compounds_from_csv(csvname, delimiter = ';'):
    """ 
    Read a csv of compounds parameters and import them in a dict.
    Expected columns are :
        - 'compound'
        - 'dose_max_ha'
        - 'type_code'
        - 'Ap'
        - 'Kp'
        - 'Ae'
        - 'Ke'
        - 'decay_rate'    
    """
    tab = recfromcsv(csvname, delimiter = delimiter, case_sensitive = True)
    d = [dict(zip(tab.dtype.names, data)) for data in tab]
    return d


def _dose_decay(decay_rate, initial_dose, days):
    active_dose = initial_dose * exp(-decay_rate * days) 
    return active_dose


class SimcyclePesticide(Exception): pass       


def milne_leaf(initial_dose, compound_parameters, days):
    """ Milne decay model applied to penetrated doses

        :Parameters:
        ----------

        - `initial_dose` (float) - Mass of product present per unit surface inside the leaf (g.m-2)
        - `compound_parameters` : A dict of compound parameters read in a .csv file 
        - `days` : time step (day)

        :Returns:
        ----------

        - 'initial_dose' (float) - Remaining mass of product per unit surface inside the leaf after dt (g.m-2) 
    """
    ptable = dict([(p['compound'],p) for p in compound_parameters])
    for name, dose in initial_dose.iteritems():
        try:
            initial_dose[name] = _dose_decay(ptable[name]['decay_rate'], initial_dose[name], days)
        except KeyError:
            raise SimcyclePesticide('product %s not found in parameter dict'%name)
    return initial_dose


class PenetratedDecayModel(object):
    """ Adaptor for Milne decay model and compliying echap pesticide_penetrated_decay model protocol
    """
    def __init__(self,products_parameters = compounds_from_csv('E:/openaleapkg/echap/src/alinea/echap/products_parameters.csv')):
        self.compound_parameters = products_parameters
    def decay(self, initial_dose, dt):
        """ Return the penetrated active doses 

        :Parameters:
        ----------
        - initial_dose: Total amount of penetrated dose in the leaf element
        - dt: Time step of the simulation (day)

        :Returns:
        -------
        - active_dose: Total amount of penetrated active dose in the leaf element after decay
        """    
        active_dose = milne_leaf(initial_dose, self.compound_parameters, dt)
        return active_dose
    def decay_and_penetrate(self, name, dose, microclimate, dt):
        """ make all the product penetrates to simulate fully the model of Milne"""
        newdose = 0
        loss = 0
        penetrated = dose
        return newdose,penetrated,loss


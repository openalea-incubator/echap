# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:16:51 2013

@author: lepse
"""
from numpy import recfromcsv
from numpy import exp


def compounds_from_csv(csvname, delimiter = ';') :
    """ 
    Read a csv of compounds parameters and import them in a dict.
    """
    tab = recfromcsv(csvname, delimiter = delimiter, case_sensitive = True)
    d = [dict(zip(tab.dtype.names, data)) for data in tab]
    return d


products_parameters = compounds_from_csv('E:/openaleapkg/echap/src/alinea/echap/products_parameters.csv')


def _dose_decay(decay_rate, initial_dose, days):
    active_dose = initial_dose * exp(-decay_rate * days) 
    return active_dose


class SimcyclePesticide(Exception): pass       


def milne_leaf(initial_dose, compound_parameters, days):
    ptable = dict([(p['compound'],p) for p in compound_parameters])
    for name, dose in initial_dose.iteritems():
        try:
            initial_dose[name] = _dose_decay(ptable[name]['decay_rate'], initial_dose[name], days)
        except KeyError:
            raise SimcyclePesticide('product %s not found in parameter dict'%name)
    return initial_dose


class PenetratedDecayModel(object):
    """ Adaptor for Milne decay model
    """
    def __init__(self,compound_parameters={}):
        self.compound_parameters = products_parameters
    def decay(self, initial_dose, dt):
        active_dose = milne_leaf(initial_dose, self.compound_parameters, dt)
        return active_dose
    def decay_and_penetrate(self, name, dose, microclimate, dt):
        """ make all the product penetrates to simulate fully the model of Milne"""
        newdose = 0
        loss = 0
        penetrated = dose
        return newdose,penetrated,loss


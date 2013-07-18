# -*- coding: utf-8 -*-
"""
Created on Thu Apr 04 11:25:45 2013

@author: lepse
"""
from pandas import *
from datetime import datetime, timedelta
from numpy import exp
from alinea.astk.caribu_interface import *
from alinea.astk.TimeControl import * 


def product_dose(product_name, dose, productsDB):
    """ 
    :Parameters:
    ------------
    - product_name (str)
        Commercial name of the product
    - dose (float)
        Application dose in l.ha-1
    - productsDB (dict)
        Dict of products name, active compounds and concentration of active compounds (in g.l-1)

    :Returns:
    ---------
    - active_dose (float)
        Dose of active compound in g.m-2
    """
    for prod, sub in productsDB.iteritems():
        if prod == product_name:
            for compound_name, compound_dose in sub.iteritems():
                active_dose = (compound_dose * dose) * 10**-4
    return compound_name, active_dose


def interception_dose(product_name, dose, scene_geometry, productsDB, elevation, azimuth):
    """ Implement pesticide interception_model using Caribu model """ 
# TODO appeller la méthode 'convert_units' de caribu interface pour convertir les m en cm
    compound_name, active_dose = product_dose(product_name, dose, productsDB)
    received_dose = emission_inv(elevation, active_dose)
    sources = (received_dose,(vecteur_direction(elevation,azimuth)))
    out_moy = run_caribu(sources, scene_geometry)
    Einc = out_moy['Einc']
    return compound_name, Einc 


class CaribuInterceptModel(object):
    """ Adaptor for Caribu model compliying echap pesticide_interception model protocol 
    productsDB (dict)
        - product: Comercial name of the product
        - compound: Active compound of the product
        - dose: Concentration of active compound in g.l-1
    elevation
    azimuth
    pest_calendar (dict)
        - datetime: List of the dates and times of the treatments (str format = "%Y-%m-%d %H:%M:%S" or datetime format)
        - dose: List of the doses of each treatment in l.ha-1 (float)
        - product_name: List of the comercial names of the products used for the treatments (str)
    """
    def __init__(self, productsDB={'Opus': {'Epoxiconazole': 125}, 'Banko 500': {'Chlorothalonil': 500}}, elevation=90, azimuth=0, pest_calendar={}): 
        self.productsDB = productsDB
        self.elevation = elevation
        self.azimuth = azimuth
        self.pest_calendar = pest_calendar

    def timing(self, start_date = "", steps = 1, delay=1, weather=None):
        """ compute timing and time_control_sets for a simulation between start and stop. return False when there is no treatment
        """
        pest_data = self.pest_calendar
        data = DataFrame(pest_data)

        def str_to_datetime(t_deb):
            format = "%Y-%m-%d %H:%M:%S"
            if isinstance(t_deb, str):
                t_deb = datetime.strptime(t_deb,format)
            return t_deb

        istart = str_to_datetime(start_date)
        stop = istart + timedelta(hours=steps-1)
        i = istart
        step = []
        while i <= stop:
            step.append(i)
            i+= timedelta(hours=1)
        event = []
        id = 0
        for i in step:
            if i == str_to_datetime(data['datetime'][id]):#.to_datetime()
                event.append([data['dose'][id], data['product_name'][id]])
                id += 1
            else:
                event.append(False)
        
        return (TimeControlSet(dose = x[0], product = x[1], dt = len(x)) if x else TimeControlSet(dose=None, product=None, dt=0) for x in event)


    def intercept(self, time_control, scene_geometry):
        """ Return the surfacic doses intercept on each leaf and stem element

        :Parameters:
        ------------
        - time_control
        - scene_geometry

        :Returns:
        ---------
        - doses: (float)
            Dict of doses (g.m-2) calculated with the interception model for each leaf and stem element
        """       
        compound_name, Einc = interception_dose(time_control.product, time_control.dose, scene_geometry, self.productsDB, self.elevation, self.azimuth)
        doses = dict([(k,{compound_name:v}) for k,v in Einc.iteritems()])
        return doses







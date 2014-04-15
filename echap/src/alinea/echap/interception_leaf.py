# -*- coding: utf-8 -*-
"""
Created on Thu Apr 04 11:25:45 2013

@author: lepse
"""
import pandas
import os
import StringIO
from numpy import exp
from alinea.astk.caribu_interface import *
from alinea.caribu.caribu_star import rain_and_light_star, caribu_rain_star



def pesticide_applications(data, sep=','):
    """ Construct a pesticide application agenda from user data
    
    - data a file of string with columns: 
        * date, with 'yyyy-mm-dd hh:mm:ss' format
        * dose (unit??), 
        * product_name
    
    return a time indexed panda dataframe
    """
    if os.path.isfile(str(data)):
       path_or_buffer = data
    else:
       path_or_buffer = StringIO.StringIO(data)
    calendar = pandas.read_csv(path_or_buffer, sep=sep, skipinitialspace = True)
    calendar.index = pandas.to_datetime(calendar['date'])
    calendar = calendar.sort_index()
     
    return calendar 



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
    compound_name = None
    active_dose = 0
    for prod, sub in productsDB.iteritems():
        if prod == product_name:
            for compound_name, compound_dose in sub.iteritems():
                active_dose = (compound_dose * dose) * 10**-4
    return compound_name, active_dose


def interception_dose(product_name, dose, scene_geometry, productsDB, elevation, azimuth):
    """ Implement pesticide interception_model using Caribu model """ 
#TODO appeller la méthode 'convert_units' de caribu interface pour convertir les m en cm
    compound_name, active_dose = product_dose(product_name, dose, productsDB)
    received_dose = emission_inv(elevation, active_dose)
    sources = (received_dose,(vecteur_direction(elevation,azimuth)))
 #   out_moy = turtle_interception(sectors='1', scene_geometry, energy=received_dose, output_by_triangle = False, convUnit = 0.01):
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



    def intercept(self, scene_geometry, product_name, dose):
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
        compound_name, Einc = interception_dose(product_name, dose, scene_geometry, self.productsDB, self.elevation, self.azimuth)
        doses = dict([(k,{compound_name:v}) for k,v in Einc.iteritems()])
        return doses



class InterceptModel(object):
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
    def __init__(self, productsDB={'Opus': {'Epoxiconazole': 125}, 'Banko 500': {'Chlorothalonil': 500}}): 
        self.productsDB = productsDB

    def intercept(self, g, product_name, dose, domain = None, convUnit = 0.01, label='LeafElement'):
        compound_name, active_dose = product_dose(product_name, dose, self.productsDB)
        doses={}
        if not 'rain_exposed_area' in g.properties():
            g = caribu_rain_star(g, output_by_triangle = False, domain = domain, convUnit = convUnit, dt = 1)
        rain_exposed_area = g.property('rain_exposed_area')
        for vid in g:
            if g.label(vid).startswith(label):
                if vid in rain_exposed_area:
                    doses[vid] = {compound_name: rain_exposed_area[vid] * active_dose}
        return doses








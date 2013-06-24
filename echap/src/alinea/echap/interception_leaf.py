# -*- coding: utf-8 -*-
"""
Created on Thu Apr 04 11:25:45 2013

@author: lepse
"""
from numpy import exp
from alinea.astk.caribu_interface import *



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

    """
    def __init__(self, productsDB={'Opus': {'Epoxiconazole': 125}, 'Banko 500': {'Chlorothalonil': 500}}, elevation=90, azimuth=0): 
        self.productsDB = productsDB
        self.elevation = elevation
        self.azimuth = azimuth

    def intercept(self, product_name, dose, scene_geometry):
        """ Return the surfacic doses intercept on each leaf and stem element

        :Parameters:
        ------------
        - product_name: (str)
            Commercial name of the product
        - dose: (float)
            Application dose of product in l.ha-1
        - scene_geometry

        :Returns:
        ---------
        - doses: (float)
            Dict of doses (g.m-2) calculated with the interception model for each leaf and stem element
        """       
        compound_name, Einc = interception_dose(product_name, dose, scene_geometry, self.productsDB, self.elevation, self.azimuth)
        doses = dict([(k,{compound_name:v}) for k,v in Einc.iteritems()])
        return doses

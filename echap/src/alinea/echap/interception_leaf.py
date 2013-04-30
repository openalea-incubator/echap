# -*- coding: utf-8 -*-
"""
Created on Thu Apr 04 11:25:45 2013

@author: lepse
"""

from alinea.caribu.CaribuScene import CaribuScene
from numpy import recfromcsv
from numpy import exp
from math import radians, degrees, sin , cos

def products_from_csv(csvname, delimiter = ';'):
    """ 
    Read a csv of products parameters and import them in a dict.
    Expected columns are :
        - 'product' : commercial name of the product
        - 'compound' : name of the active compound of the product
        - 'dose' : dose of active compound in the product (g.l-1)    
    """
    tab = recfromcsv(csvname, delimiter = delimiter, case_sensitive = True)
    d = {}
    for i in range(0,len(tab['compound'])):
        d[tab[i][0]] = {tab[i][1]:tab[i][2]}
    return d


def product_dose(product_name, dose, productsDB):
    """ 
    :Parameters:
    ----------

    - product_name: Commercial name of the product
    - dose: Application dose in l.ha-1
    - productsDB: dict of products name, active compounds and concentration of active compounds

    :Returns:
    ----------

    - active_dose: Dose of active compound in g.m-2
    """
    for prod, sub in productsDB.iteritems():
        if prod == product_name:
            for compound_name, compound_dose in sub.iteritems():
                active_dose = (compound_dose * dose) * 10**-4
    return compound_name, active_dose


def _vecteur_direction(elevation,azimuth):
    theta = radians(90 - elevation)
    phi = radians(azimuth)
    return sin(theta) * cos(phi),sin(theta) * sin(phi),  -cos(theta)


def _emission_inv(elevation, product_name, dose, productsDB):
    """ return energy of emmision for a source of a given direction and of a given energy received on a horizontal surface """
    theta = radians(90 - elevation)
    compound_name, active_dose = product_dose(product_name, dose, productsDB) 
    receive_dose = active_dose * abs(cos(theta))
    return compound_name, receive_dose


def interception_dose(product_name, dose, productsDB, elevation, azimuth, scene):
    """ Implement pesticide interception_model using Caribu model """
    compound_name, receive_dose = _emission_inv(elevation, product_name, dose, productsDB)
    source = (receive_dose,(_vecteur_direction(elevation,azimuth))) # Quantité reçue
    c_scene = CaribuScene()    
    idmap = c_scene.add_Shapes(scene)      
    c_scene.addSources(source)
    c_scene.getIncidentEnergy()
    output = c_scene.runCaribu(infinity=False)
    Einc = c_scene.output_by_id(output, idmap)['Einc']
    return compound_name, Einc 


class CaribuInterceptModel(object):
    """ Adaptor for Caribu model compliying echap pesticide_interception model protocol """
    def __init__(self, productsDB=products_from_csv('E:/openaleapkg/echap/src/alinea/echap/products_names.csv'), elevation=90, azimuth=0): 
        self.productsDB = productsDB
        self.elevation = elevation
        self.azimuth = azimuth
    def intercept(self, product_name, dose, scene):
        """ Return the surfacic doses intercept on each leaf and stem element

        :Parameters:
        ----------

        - product_name: Commercial name of the product
        - dose: Application dose in l.ha-1
        - scene: Scene containing the simulated system

        :Returns:
        -------

        - doses: Dict of doses calculated with the interception model for each leaf and stem element
        """       
        compound_name, Einc = interception_dose(product_name, dose, self.productsDB, self.elevation, self.azimuth, scene)
        doses = dict([(k,{compound_name:v}) for k,v in Einc.iteritems()])
        return doses

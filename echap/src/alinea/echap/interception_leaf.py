# -*- coding: utf-8 -*-
"""
Created on Thu Apr 04 11:25:45 2013

@author: lepse
"""

from alinea.caribu.CaribuScene import CaribuScene
from numpy import recfromcsv
from numpy import exp
from math import radians, degrees, sin , cos

# def products_from_csv(csvname, delimiter = ';') :
    # """ 
    # Read a csv of products names and import them in a dict.
    # """
    # tab = recfromcsv(csvname, delimiter = delimiter, case_sensitive = True)
    # d = [dict(zip(tab.dtype.names, data)) for data in tab]
    # b = dict([(k,{d:v}) for k,v in b.iteritems()])
    # return b

# productsDB = products_from_csv('E:/openaleapkg/echap/src/alinea/echap/products_names.csv')


def product_dose(product_name, dose, productsDB):
    ''' Parameters:
            product_name: Commercial name of the product
            dose: Application dose in l/ha
            productsDB: Dict of products name and concentration of active compounds
        Outputs:
            active_dose: Dose of active compound in g/cm²
    '''
    for prod, sub in dict(productsDB).iteritems():
        if prod == product_name:
            for compound_name, compound_dose in sub.iteritems():
                active_dose = (compound_dose * dose) * 10**-8
    return compound_name, active_dose

""" Implement pesticide interception_model using Caribu model
"""

def _vecteur_direction(elevation,azimuth):
    theta = radians(90 - elevation)
    phi = radians(azimuth)
    return sin(theta) * cos(phi),sin(theta) * sin(phi),  -cos(theta)

def _emission_inv(elevation, product_name, dose, productsDB):
    """ return energy of emmision for a source of a given direction and of a given energy received on a horizontal surface"""
    theta = radians(90 - elevation)
    name, dose = product_dose(product_name, dose, productsDB) 
    receive_dose = dose * abs(cos(theta))
    return name, receive_dose

def interception_dose(product_name, dose, productsDB, elevation, azimuth, scene):
    name, receive_dose = _emission_inv(elevation, product_name, dose, productsDB)
    source = (receive_dose,(_vecteur_direction(elevation,azimuth))) # Quantité reçue
    c_scene = CaribuScene()    
    idmap = c_scene.add_Shapes(scene)      
    c_scene.addSources(source)
    c_scene.getIncidentEnergy()
    output = c_scene.runCaribu(infinity=False)
    Einc = c_scene.output_by_id(output, idmap)['Einc']
    return name, Einc 

class CaribuInterceptModel(object):
    """ Adaptor for Caribu model compliying echap pesticide interception_model protocol 
    """
    def __init__(self, product_name='Ignite', dose=200, productsDB={'Ignite':{'Epoxiconazole':83}, 'Nevo':{'Chlorothalonil':375}, 'Caramba':{'Metconazole':60}}, elevation=90, azimuth=0):
        self.dose = dose
        self.product_name = product_name     
        self.productsDB = productsDB
        self.elevation = elevation
        self.azimuth = azimuth
    def intercept(self, scene):
        compound_name, Einc = interception_dose(self.product_name, self.dose, self.productsDB, self.elevation, self.azimuth, scene)
        doses = dict([(k,{compound_name:v}) for k,v in Einc.iteritems()])
        return doses

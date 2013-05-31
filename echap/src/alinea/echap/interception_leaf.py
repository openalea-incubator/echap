# -*- coding: utf-8 -*-
"""
Created on Thu Apr 04 11:25:45 2013

@author: lepse
"""

from alinea.caribu.CaribuScene import CaribuScene
from numpy import exp
from math import radians, degrees, sin , cos
from openalea.plantgl import all as pgl


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


def geom2shape(vid, mesh):
    shape = pgl.Shape(mesh)
    shape.id = vid
    return shape


def run_caribu(sources, scene_geometry):
    c_scene = CaribuScene()
    shapes=[geom2shape(k,v) for k,v in scene_geometry.iteritems()]
    idmap = c_scene.add_Shapes(shapes)    
    c_scene.addSources(sources)
    output = c_scene.runCaribu(infinity=False)
    out_moy = c_scene.output_by_id(output, idmap)
    out_tri = c_scene.output_by_id(output, idmap, aggregate = False)
    return out_moy, out_tri


def interception_dose(product_name, dose, scene_geometry, productsDB, elevation, azimuth):
    """ Implement pesticide interception_model using Caribu model """
    compound_name, receive_dose = _emission_inv(elevation, product_name, dose, productsDB)
    sources = (receive_dose,(_vecteur_direction(elevation,azimuth))) # Quantité reçue
    out_moy,_ = run_caribu(sources, scene_geometry)
    Einc = out_moy['Einc']
    return compound_name, Einc 


class CaribuInterceptModel(object):
    """ Adaptor for Caribu model compliying echap pesticide_interception model protocol 
    
    productsDB: Dict
        - product: Comercial name of the product
        - compound: Active compound of the product
        - dose: Concentration of active compound in g.l-1

    """
    def __init__(self, productsDB={'Opus': {'Epoxiconazole': 83}, 'Banko 500': {'Chlorothalonil': 500}}, elevation=90, azimuth=0): 
        self.productsDB = productsDB
        self.elevation = elevation
        self.azimuth = azimuth

    def intercept(self, product_name, dose, scene_geometry):
        """ Return the surfacic doses intercept on each leaf and stem element

        :Parameters:
        ----------
        - product_name: Str
            Commercial name of the product
        - dose: Float
            Application dose of product in l.ha-1
        - scene_geometry

        :Returns:
        -------
        - doses: fLOAT
            Dict of doses (g.m-2) calculated with the interception model for each leaf and stem element
        """       
        compound_name, Einc = interception_dose(product_name, dose, scene_geometry, self.productsDB, self.elevation, self.azimuth)
        doses = dict([(k,{compound_name:v}) for k,v in Einc.iteritems()])
        return doses

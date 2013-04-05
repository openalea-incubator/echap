
# This file has been generated at Fri Apr 05 14:36:30 2013

from openalea.core import *


__name__ = 'Alinea.Echap.Models'

__editable__ = True
__description__ = ''
__license__ = 'CeCILL-C'
__url__ = 'http://openalea.gforge.inria.fr'
__alias__ = []
__version__ = '0.8.0'
__authors__ = ''
__institutes__ = ''
__icon__ = ''


__all__ = ['models_nodes_milne_leaf', 'models_nodes_interception_leaf']



models_nodes_milne_leaf = Factory(name='Milne_leaf',
                authors=' (wralea authors)',
                description='Milne model for penetrated pesticide decay',
                category='model',
                nodemodule='models_nodes',
                nodeclass='milne_leaf',
                inputs=[{'interface': IDict, 'name': 'compound_parameters', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'active_dose', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




models_nodes_interception_leaf = Factory(name='Interception_leaf',
                authors=' (wralea authors)',
                description='Caribu adaptor for pesticide interception',
                category='model',
                nodemodule='models_nodes',
                nodeclass='interception_leaf',
                inputs=[{'interface': IStr, 'name': 'product_name', 'value': None, 'desc': 'Commercial product name '}, {'interface': IFloat, 'name': 'dose', 'value': None, 'desc': 'Dose of product in l/ha'}, {'interface': IDict, 'name': 'productsDB', 'value': None, 'desc': ''}, {'interface': IFloat, 'name': 'elevation', 'value': 90, 'desc': 'Angle  from horizontal'}, {'interface': IFloat, 'name': 'azimuth', 'value': 0, 'desc': ''}, {'interface': None, 'name': 'scene', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'doses', 'desc': 'doses of active compound in g/m'}],
                widgetmodule=None,
                widgetclass=None,
               )





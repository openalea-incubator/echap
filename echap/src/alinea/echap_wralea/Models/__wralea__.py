
# This file has been generated at Mon Apr 22 10:38:35 2013

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


__all__ = ['models_nodes_CaribuMicroclimModel', 'models_nodes_PenetratedDecayModel', 'models_nodes_CaribuInterceptModel']



models_nodes_CaribuMicroclimModel = Factory(name='Microclimate_leaf',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='models_nodes',
                nodeclass='CaribuMicroclimModel',
                inputs=[{'interface': IStr, 'name': 'sectors', 'value': '46', 'desc': ''}, {'interface': IFloat, 'name': 'energy', 'value': 1, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'microclimate', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




models_nodes_PenetratedDecayModel = Factory(name='Milne_leaf',
                authors=' (wralea authors)',
                description='Milne model for penetrated pesticide decay',
                category='model',
                nodemodule='models_nodes',
                nodeclass='PenetratedDecayModel',
                inputs=[{'interface': IDict, 'name': 'compound_parameters', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'active_dose', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




models_nodes_CaribuInterceptModel = Factory(name='Interception_leaf',
                authors=' (wralea authors)',
                description='Caribu adaptor for pesticide interception',
                category='model',
                nodemodule='models_nodes',
                nodeclass='CaribuInterceptModel',
                inputs=[{'interface': IDict, 'name': 'productsDB', 'value': None, 'desc': ''}, {'interface': IFloat, 'name': 'elevation', 'value': 90, 'desc': 'Angle  from horizontal'}, {'interface': IFloat, 'name': 'azimuth', 'value': 0, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'doses', 'desc': 'doses of active compound in g/m'}],
                widgetmodule=None,
                widgetclass=None,
               )





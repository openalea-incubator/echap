
# This file has been generated at Mon Apr 29 11:42:34 2013

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


__all__ = ['models_nodes_CaribuMicroclimModel', 'models_nodes_PenetratedDecayModel', 'models_nodes_Meteo', 'models_nodes_CaribuInterceptModel']



models_nodes_CaribuMicroclimModel = Factory(name='Microclimate_leaf',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='models_nodes',
                nodeclass='CaribuMicroclimModel',
                inputs=[{'interface': IStr, 'name': 'sectors', 'value': 46, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'local_meteo', 'desc': ''}],
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




models_nodes_Meteo = Factory(name='Global_meteo',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='models_nodes',
                nodeclass='Meteo',
                inputs=[{'interface': IFileStr, 'name': 'data_file', 'value': None, 'desc': 'Path for read the meteo data file csv'}],
                outputs=[{'interface': None, 'name': 'mean_globalclimate, globalclimate, t_deb', 'desc': ''}],
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





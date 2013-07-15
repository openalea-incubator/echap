
# This file has been generated at Fri May 17 13:44:00 2013

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


__all__ = ['interception_leaf_CaribuInterceptModel', 'microclimate_leaf_MicroclimateLeaf', 'milne_leaf_PenetratedDecayModel']



interception_leaf_CaribuInterceptModel = Factory(name='Interception_Leaf',
                authors=' (wralea authors)',
                description='Caribu adaptor for pesticide interception',
                category='model',
                nodemodule='alinea.echap.interception_leaf',
                nodeclass='CaribuInterceptModel',
                inputs=[{'interface': IDict, 'name': 'productsDB', 'value': None, 'desc': ''}, {'interface': IFloat, 'name': 'elevation', 'value': 90, 'desc': 'Angle  from horizontal'}, {'interface': IFloat, 'name': 'azimuth', 'value': 0, 'desc': ''}, {'interface': IDict, 'name': 'pest_calendar', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'doses', 'desc': 'doses of active compound in g/m'}],
                widgetmodule=None,
                widgetclass=None,
               )




microclimate_leaf_MicroclimateLeaf = Factory(name='Microclimate_Leaf',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='alinea.echap.microclimate_leaf',
                nodeclass='MicroclimateLeaf',
                inputs=[{'interface': IStr, 'name': 'sectors', 'value': '16', 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'local_meteo', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




milne_leaf_PenetratedDecayModel = Factory(name='Milne_Leaf',
                authors=' (wralea authors)',
                description='Milne model for penetrated pesticide decay',
                category='model',
                nodemodule='alinea.echap.milne_leaf',
                nodeclass='PenetratedDecayModel',
                inputs=[{'interface': IDict, 'name': 'compound_parameters', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'active_dose', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )





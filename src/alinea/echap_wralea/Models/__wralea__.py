
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


productDB={'Opus': {'Epoxiconazole': 125}, 'Banko 500': {'Chlorothalonil': 500}}

interception_leaf_CaribuInterceptModel = Factory(name='Interception_Leaf',
                authors=' (wralea authors)',
                description='Caribu adaptor for pesticide interception',
                category='model',
                nodemodule='alinea.echap.interception_leaf',
                nodeclass='CaribuInterceptModel',
                inputs=[{'interface': IDict, 'name': 'productsDB', 'value': productDB, 'desc': ''},
                        {'interface': IFloat, 'name': 'elevation', 'value': 90, 'desc': 'Angle  from horizontal'}, 
                        {'interface': IFloat, 'name': 'azimuth', 'value': 0, 'desc': ''},
                        ],
                outputs=[{'interface': None, 'name': 'interception model'}],
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

microclimate_leafmicroclimate = Factory(name='Leaf Microclimate',
                nodemodule='alinea.echap.microclimate_leaf',
                nodeclass='microclimate_leaf',
               )
__all__.append('microclimate_leafmicroclimate')

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

interception_applications = Factory(name='pesticide applications',
                            nodemodule = 'alinea.echap.interception_leaf',
                            nodeclass='pesticide_applications',
                            )
__all__.append('interception_applications')                            

outputs_leafelement_recorder = Factory(name='LeafElementRecorder',
                            nodemodule = 'alinea.echap.recorder',
                            nodeclass='LeafElementRecorder',
                            )
__all__.append('outputs_leafelement_recorder')   

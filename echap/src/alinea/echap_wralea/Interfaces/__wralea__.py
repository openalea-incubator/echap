
# This file has been generated at Tue Apr 16 11:49:47 2013

from openalea.core import *


__name__ = 'Alinea.Echap.Interfaces'

__editable__ = True
__description__ = ''
__license__ = 'CeCILL-C'
__url__ = 'http://openalea.gforge.inria.fr'
__alias__ = []
__version__ = '0.8.0'
__authors__ = ''
__institutes__ = ''
__icon__ = ''


__all__ = ['interfaces_nodes_pesticide_penetrated_decay', 'interfaces_nodes_local_microclimate', 'interfaces_nodes_pesticide_surfacic_decay', 'interfaces_nodes_pesticide_interception']



interfaces_nodes_pesticide_penetrated_decay = Factory(name='pesticide_penetrated_decay',
                authors=' (wralea authors)',
                description='Update the decay of penetrated doses of pesticide on the MTG',
                category='Unclassified',
                nodemodule='interfaces_nodes',
                nodeclass='pesticide_penetrated_decay',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': 'MTG with penetrated doses of pesticide as property'}, {'interface': None, 'name': 'decay_model', 'value': None, 'desc': 'Model of penetrated pesticide decay'}, {'interface': IStr, 'name': 'label', 'value': 'LeafElement', 'desc': ''}, {'interface': IInt, 'name': 'timestep', 'value': 1, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'g', 'desc': 'Update MTG with decay of the penetrated doses'}],
                widgetmodule=None,
                widgetclass=None,
               )




interfaces_nodes_local_microclimate = Factory(name='local_microclimate',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='interfaces_nodes',
                nodeclass='local_microclimate',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': ''}, {'interface': None, 'name': 'scene', 'value': None, 'desc': ''}, {'interface': None, 'name': 'climate_model', 'value': None, 'desc': ''}, {'interface': IInt, 'name': 'rain', 'value': 0, 'desc': ''}, {'interface': IStr, 'name': 'label', 'value': 'LeafElement', 'desc': ''}],
                outputs=[{'interface': None, 'name': 'g', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




interfaces_nodes_pesticide_surfacic_decay = Factory(name='pesticide_surfacic_decay',
                authors=' (wralea authors)',
                description='Interface between g and the decay model of Pearl',
                category='Unclassified',
                nodemodule='interfaces_nodes',
                nodeclass='pesticide_surfacic_decay',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': 'MTG'}, {'interface': None, 'name': 'decay_model', 'value': None, 'desc': 'Pearl decay model'}, {'interface': IStr, 'name': 'label', 'value': 'LeafElement', 'desc': 'Leaf element of MTG'}, {'interface': IInt, 'name': 'timestep', 'value': 24, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'g', 'desc': 'New MTG with pesticide dose on leaf'}],
                widgetmodule=None,
                widgetclass=None,
               )




interfaces_nodes_pesticide_interception = Factory(name='pesticide_interception',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='interfaces_nodes',
                nodeclass='pesticide_interception',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': ''}, {'interface': None, 'name': 'scene', 'value': None, 'desc': ''}, {'interface': None, 'name': 'interception_model', 'value': None, 'desc': ''}, {'interface': IStr, 'name': 'product_name', 'value': None, 'desc': 'Commercial name of the product'}, {'interface': IFloat, 'name': 'dose', 'value': None, 'desc': 'Dose in l.ha-1'}, {'interface': IStr, 'name': 'label', 'value': 'LeafElement', 'desc': ''}],
                outputs=[{'interface': None, 'name': 'g', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )





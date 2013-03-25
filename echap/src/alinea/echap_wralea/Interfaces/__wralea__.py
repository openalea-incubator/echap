
# This file has been generated at Mon Mar 25 16:35:42 2013

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


__all__ = ['interfaces_nodes_decay_leaf', 'interfaces_nodes_pesticide_surfacic_decay']



interfaces_nodes_decay_leaf = Factory(name='decay_leaf',
                authors=' (wralea authors)',
                description='Update the decay of penetrated doses of pesticide on the MTG',
                category='Unclassified',
                nodemodule='interfaces_nodes',
                nodeclass='decay_leaf',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': 'MTG with penetrated doses of pesticide as property'}, {'interface': IStr, 'name': 'label', 'value': 'LeafElement', 'desc': ''}, {'interface': None, 'name': 'decay_model', 'value': None, 'desc': 'Model of penetrated pesticide decay'}],
                outputs=[{'interface': None, 'name': 'g', 'desc': 'Update MTG with decay of the penetrated doses'}],
                widgetmodule=None,
                widgetclass=None,
               )




interfaces_nodes_pesticide_surfacic_decay = Factory(name='pesticide_surfacic_decay',
                authors=' (wralea authors)',
                description='Interface between g and the decay model of Pearl',
                category='Unclassified',
                nodemodule='interfaces_nodes',
                nodeclass='pesticide_surfacic_decay',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': 'MTG'}, {'interface': None, 'name': 'decay_model', 'value': None, 'desc': 'Pearl decay model'}, {'interface': IStr, 'name': 'label', 'value': 'LeafElement', 'desc': 'Leaf element of MTG'}],
                outputs=[{'interface': None, 'name': 'g', 'desc': 'New MTG with pesticide dose on leaf'}],
                widgetmodule=None,
                widgetclass=None,
               )





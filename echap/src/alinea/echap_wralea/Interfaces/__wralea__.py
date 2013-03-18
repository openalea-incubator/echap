
# This file has been generated at Mon Mar 18 16:27:57 2013

from openalea.core import *


__name__ = 'Alinea.Echap.Interfaces'

__editable__ = True
__description__ = ''
__license__ = 'CeCILL-C'
__url__ = 'http://openalea.gforge.inria.fr'
__alias__ = []
__version__ = '0.8.0'
__authors__ = ''
__institutes__ = None
__icon__ = ''


__all__ = []



intefaces_nodes_pesticide_surfacic_decay = Factory(name='pesticide_surfacic_decay',
                authors=' (wralea authors)',
                description='Interface between g and the Pearl decay model',
                category='Unclassified',
                nodemodule='pesticide_surfacic_decay',
                nodeclass='pesticide_surfacic_decay',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': 'MTG'}, {'interface': None, 'name': 'decay_model', 'value': None, 'desc': 'PearlLeafDecayModel'}, {'interface': IStr, 'name': 'label', 'value': None, 'desc': 'MTG element'}],
                outputs=[{'interface': None, 'name': 'g', 'desc': 'new MTG with pesticide dose'}],
                widgetmodule=None,
                widgetclass=None,
               )
__all__.append('interfaces_nodes_pesticide_surfacic_decay')




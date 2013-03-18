
# This file has been generated at Mon Mar 18 14:42:01 2013

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



interfaces_nodes_pesticide_surfacic_decay = Factory(name='pesticide_surfacic_decay',
                authors=' (wralea authors)',
                description='Interface between g and the Pearl decay model',
                category='Interfaces',
                nodemodule='pesticide_surfacic_decay',
                nodeclass='pesticide_surfacic_decay',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': ''}, {'interface': None, 'name': 'decay_model', 'value': None, 'desc': ''}, {'interface': IStr, 'name': 'label', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'out', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )
__all__.append('interfaces_nodes_pesticide_surfacic_decay')





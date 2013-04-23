
# This file has been generated at Tue Apr 23 17:43:09 2013

from openalea.core import *


__name__ = 'Alinea.Echap.Tests_nodes'

__editable__ = True
__description__ = ''
__license__ = 'CeCILL-C'
__url__ = 'http://openalea.gforge.inria.fr'
__alias__ = []
__version__ = '0.8.0'
__authors__ = ''
__institutes__ = ''
__icon__ = ''


__all__ = ['csv_to_dict_prod_csv_to_dict_prod', 'tests_nodes_plot_pesticide']



csv_to_dict_prod_csv_to_dict_prod = Factory(name='csv_to_dict_prod',
                authors=' (wralea authors)',
                description='Read a products names csv file and convert it to a dict',
                category='Unclassified',
                nodemodule='csv_to_dict_prod',
                nodeclass='csv_to_dict_prod',
                inputs=[{'interface': None, 'name': 'csvname', 'value': None, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'product_dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




tests_nodes_plot_pesticide = Factory(name='plot_pesticide',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='tests_nodes',
                nodeclass='plot_pesticide',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'OUT1', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )





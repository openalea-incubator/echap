
# This file has been generated at Fri Apr 26 16:12:22 2013

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


__all__ = ['tests_nodes_plot_pesticide', 'tests_nodes_update_meteo_date', 'tests_nodes_compounds_from_csv', 'tests_nodes_products_from_csv']



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




tests_nodes_update_meteo_date = Factory(name='update_meteo_date',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='tests_nodes',
                nodeclass='update_meteo_date',
                inputs=[{'interface': None, 'name': 'reader_meteo', 'value': None, 'desc': ''}, {'interface': IInt, 'name': 'timesteps', 'value': 1, 'desc': ''}, {'interface': IStr, 'name': 't_deb', 'value': '2000-10-01 01:00:00', 'desc': ''}],
                outputs=[{'interface': IStr, 'name': 't_deb', 'desc': "String date updated in '%Y-%m-%d %H:%M:%S' format"}],
                widgetmodule=None,
                widgetclass=None,
               )




tests_nodes_compounds_from_csv = Factory(name='compounds_from_csv',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='tests_nodes',
                nodeclass='compounds_from_csv',
                inputs=[{'interface': IFileStr, 'name': 'csvname', 'value': None, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'compound_dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




tests_nodes_products_from_csv = Factory(name='products_from_csv',
                authors=' (wralea authors)',
                description='Read a products names csv file and convert it to a dict',
                category='Unclassified',
                nodemodule='tests_nodes',
                nodeclass='products_from_csv',
                inputs=[{'interface': IFileStr, 'name': 'csvname', 'value': None, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'product_dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )





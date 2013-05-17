
# This file has been generated at Mon May 13 16:37:22 2013

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


__all__ = ['tests_nodes_sum_temp_global', 'tests_nodes_plot_pesticide', 'tests_nodes_compounds_from_csv', 'tests_nodes_interface_meteo', 'tests_nodes_update_meteo_date', 'tests_nodes_products_from_csv']



tests_nodes_sum_temp_global = Factory(name='sum_temp_global',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='alinea.echap.tests_nodes',
                nodeclass='sum_temp_global',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': ''}, {'interface': None, 'name': 'globalclimate', 'value': None, 'desc': 'Pandas dataframe with global climate for the time step'}],
                outputs=[{'interface': None, 'name': 'g', 'desc': 'MTG with the sum of temperature at the plant scale'}],
                widgetmodule=None,
                widgetclass=None,
               )




tests_nodes_plot_pesticide = Factory(name='plot_pesticide',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='alinea.echap.tests_nodes',
                nodeclass='plot_pesticide',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'OUT1', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




tests_nodes_compounds_from_csv = Factory(name='compounds_from_csv',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='alinea.echap.tests_nodes',
                nodeclass='compounds_from_csv',
                inputs=[{'interface': IFileStr, 'name': 'csvname', 'value': None, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'compound_dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




tests_nodes_interface_meteo = Factory(name='interface_meteo',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='alinea.echap.tests_nodes',
                nodeclass='interface_meteo',
                inputs=[{'interface': None, 'name': 'meteo_reader', 'value': None, 'desc': ''}, {'interface': IInt, 'name': 'timestep', 'value': 1, 'desc': ''}, {'interface': IStr, 'name': 't_deb', 'value': '2000-10-01 01:00:00', 'desc': ''}],
                outputs=[{'interface': None, 'name': 'mean_globalclimate', 'desc': 'Pandas dataframe'}, {'interface': None, 'name': 'globalclimate', 'desc': 'Pandas dataframe'}, {'interface': IStr, 'name': 't_deb', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




tests_nodes_update_meteo_date = Factory(name='update_meteo_date',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='alinea.echap.tests_nodes',
                nodeclass='update_meteo_date',
                inputs=[{'interface': IInt, 'name': 'timesteps', 'value': 1, 'desc': ''}, {'interface': IStr, 'name': 't_deb', 'value': '2000-10-01 01:00:00', 'desc': ''}],
                outputs=[{'interface': IStr, 'name': 't_deb', 'desc': "String date updated in '%Y-%m-%d %H:%M:%S' format"}],
                widgetmodule=None,
                widgetclass=None,
               )




tests_nodes_products_from_csv = Factory(name='products_from_csv',
                authors=' (wralea authors)',
                description='Read a products names csv file and convert it to a dict',
                category='Unclassified',
                nodemodule='alinea.echap.tests_nodes',
                nodeclass='products_from_csv',
                inputs=[{'interface': IFileStr, 'name': 'csvname', 'value': None, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'product_dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )





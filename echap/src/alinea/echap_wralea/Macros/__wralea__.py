
# This file has been generated at Wed Nov 13 23:36:59 2013

from openalea.core import *


__name__ = 'Alinea.Echap.Macros'

__editable__ = True
__description__ = ''
__license__ = 'CeCILL-C'
__url__ = 'http://openalea.gforge.inria.fr'
__alias__ = []
__version__ = '0.8.0'
__authors__ = ''
__institutes__ = ''
__icon__ = ''


__all__ = ['_215033328', '_132912976', '_132912944', '_215416464']



_215033328 = CompositeNodeFactory(name='every_rain',
                             description='',
                             category='Unclassified',
                             doc='',
                             inputs=[  {  'desc': '', 'interface': None, 'name': 'time_sequence', 'value': None},
   {  'desc': '', 'interface': None, 'name': 'weather', 'value': None}],
                             outputs=[  {  'desc': '', 'interface': None, 'name': 'values'},
   {  'desc': '', 'interface': None, 'name': 'delays'}],
                             elt_factory={  2: ('alinea.astk', 'rain filter'), 4: ('alinea.astk', 'time_control')},
                             elt_connections={  34429332: (2, 2, 4, 2),
   34429368: ('__in__', 0, 2, 0),
   34429392: (4, 1, '__out__', 1),
   34429404: (4, 0, '__out__', 0),
   34429428: (2, 0, 4, 0),
   34429440: (2, 1, 4, 1),
   34429452: ('__in__', 1, 2, 1)},
                             elt_data={  2: {  'block': False,
         'caption': 'rain filter',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x08F20B50> : "rain filter"',
         'hide': True,
         'id': 2,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 21.210606912550276,
         'posy': 18.665334083044247,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   4: {  'block': False,
         'caption': 'time_control',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x08F20B70> : "time_control"',
         'hide': True,
         'id': 4,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 13.088136466949877,
         'posy': 56.2636414753017,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   '__in__': {  'block': False,
                'caption': 'In',
                'delay': 0,
                'hide': True,
                'id': 0,
                'lazy': True,
                'port_hide_changed': set(),
                'posx': 33.65057628578118,
                'posy': -37.63681792623492,
                'priority': 0,
                'use_user_color': False,
                'user_application': None,
                'user_color': None},
   '__out__': {  'block': False,
                 'caption': 'Out',
                 'delay': 0,
                 'hide': True,
                 'id': 1,
                 'lazy': True,
                 'port_hide_changed': set(),
                 'posx': 33.71948496775704,
                 'posy': 97.54027564092196,
                 'priority': 0,
                 'use_user_color': False,
                 'user_application': None,
                 'user_color': None}},
                             elt_value={  2: [(1, 'None')], 4: [], '__in__': [], '__out__': [(0, 'None'), (1, 'None')]},
                             elt_ad_hoc={  2: {'position': [21.210606912550276, 18.665334083044247], 'userColor': None, 'useUserColor': False},
   4: {'position': [13.088136466949877, 56.2636414753017], 'userColor': None, 'useUserColor': False},
   '__in__': {'position': [33.65057628578118, -37.63681792623492], 'userColor': None, 'useUserColor': False},
   '__out__': {'position': [33.71948496775704, 97.54027564092196], 'userColor': None, 'useUserColor': False}},
                             lazy=True,
                             eval_algo='LambdaEvaluation',
                             )




_132912976 = CompositeNodeFactory(name='every_degreedays',
                             description='',
                             category='time',
                             doc='',
                             inputs=[  {  'desc': '', 'interface': None, 'name': 'time_sequence', 'value': None},
   {  'desc': '', 'interface': None, 'name': 'delay', 'value': None},
   {  'desc': '', 'interface': None, 'name': 'weather', 'value': None}],
                             outputs=[  {  'desc': '', 'interface': None, 'name': 'values'},
   {  'desc': '', 'interface': None, 'name': 'delays'}],
                             elt_factory={  2: ('alinea.astk', 'DegreeDay'),
   3: ('alinea.astk', 'thermal_time filter'),
   4: ('alinea.astk', 'time_control')},
                             elt_connections={  34429356: ('__in__', 2, 3, 1),
   34429368: ('__in__', 1, 3, 3),
   34429380: ('__in__', 0, 3, 0),
   34429392: (4, 1, '__out__', 1),
   34429404: (4, 0, '__out__', 0),
   34429416: (3, 2, 4, 2),
   34429428: (3, 0, 4, 0),
   34429440: (3, 1, 4, 1),
   34429452: (2, 0, 3, 2)},
                             elt_data={  2: {  'block': False,
         'caption': 'DegreeDay',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x08F20C70> : "DegreeDay"',
         'hide': True,
         'id': 2,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 87.25353418395848,
         'posy': -15.190519000709614,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   3: {  'block': False,
         'caption': 'thermal_time filter',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x08F20BF0> : "thermal_time filter"',
         'hide': True,
         'id': 3,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 5.582428903504194,
         'posy': 31.209487169021386,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   4: {  'block': False,
         'caption': 'time_control',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x08F20B70> : "time_control"',
         'hide': True,
         'id': 4,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 13.512348605200884,
         'posy': 66.02052065507483,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   '__in__': {  'block': False,
                'caption': 'In',
                'delay': 0,
                'hide': True,
                'id': 0,
                'lazy': True,
                'port_hide_changed': set(),
                'posx': 33.65057628578117,
                'posy': -38.90945434098794,
                'priority': 0,
                'use_user_color': False,
                'user_application': None,
                'user_color': None},
   '__out__': {  'block': False,
                 'caption': 'Out',
                 'delay': 0,
                 'hide': True,
                 'id': 1,
                 'lazy': True,
                 'port_hide_changed': set(),
                 'posx': 33.71948496775704,
                 'posy': 97.54027564092196,
                 'priority': 0,
                 'use_user_color': False,
                 'user_application': None,
                 'user_color': None}},
                             elt_value={  2: [(0, '0')], 3: [], 4: [], '__in__': [], '__out__': []},
                             elt_ad_hoc={  2: {'position': [87.25353418395848, -15.190519000709614], 'userColor': None, 'useUserColor': False},
   3: {'position': [5.582428903504194, 31.209487169021386], 'userColor': None, 'useUserColor': False},
   4: {'position': [13.512348605200884, 66.02052065507483], 'userColor': None, 'useUserColor': False},
   '__in__': {'position': [33.65057628578117, -38.90945434098794], 'userColor': None, 'useUserColor': False},
   '__out__': {'position': [33.71948496775704, 97.54027564092196], 'userColor': None, 'useUserColor': False}},
                             lazy=True,
                             eval_algo='LambdaEvaluation',
                             )




_132912944 = CompositeNodeFactory(name='adel_mtg3',
                             description='Build a mtg',
                             category='composite',
                             doc='',
                             inputs=[],
                             outputs=[],
                             elt_factory={  2: ('Alep.Test_nodes', 'adel_mtg3'),
   3: ('Alep.Test_nodes', 'leaves_db'),
   4: ('Alep.Test_nodes', 'temp_plot3D'),
   5: ('alinea.adel.simulation', 'RunAdel'),
   6: ('alinea.adel.parameterisation', 'setAdel'),
   7: ('alinea.adel.parameterisation', 'devCsv'),
   8: ('alinea.adel.parameterisation', 'geoLeaf'),
   9: ('alinea.adel.parameterisation', 'geoAxe'),
   10: ('alinea.adel.data', 'axeTCa0N.csv'),
   11: ('alinea.adel.data', 'dimTCa0N.csv'),
   12: ('alinea.adel.data', 'phenTCa0N.csv'),
   13: ('alinea.adel.data', 'earTCa0N.csv'),
   14: ('alinea.adel.data', 'ssi2sen.csv'),
   15: ('alinea.adel.parameterisation', 'freeGeoLeaf'),
   16: ('openalea.flow control', 'annotation'),
   17: ('openalea.data structure', 'int'),
   20: ('openalea.data structure', 'int'),
   21: ('alinea.adel.data', 'leaves_simple.db')},
                             elt_connections={  33967120: (11, 0, 7, 1),
   33967132: (13, 0, 7, 3),
   33967144: (14, 0, 7, 4),
   33967156: (20, 0, 5, 0),
   33967168: (7, 0, 6, 0),
   33967180: (2, 0, 4, 0),
   33967192: (17, 0, 6, 3),
   33967204: (3, 0, 2, 1),
   33967216: (8, 0, 6, 1),
   33967228: (6, 0, 5, 1),
   33967240: (17, 0, 2, 3),
   33967252: (12, 0, 7, 2),
   33967264: (9, 0, 6, 2),
   33967276: (10, 0, 7, 0),
   33967288: (21, 0, 3, 0),
   33967300: (5, 0, 2, 2)},
                             elt_data={  2: {  'block': False,
         'caption': 'adel_mtg3',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x09616710> : "adel_mtg3"',
         'hide': True,
         'id': 2,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 68.86480669392043,
         'posy': 60.841722418900574,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   3: {  'block': False,
         'caption': 'leaves_db',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x09616770> : "leaves_db"',
         'hide': True,
         'id': 3,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': -63.51608384390719,
         'posy': 15.377578193788057,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   4: {  'block': False,
         'caption': 'temp_plot3D',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x09616890> : "temp_plot3D"',
         'hide': True,
         'id': 4,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 68.86480669392043,
         'posy': 145.08410730660907,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   5: {  'block': False,
         'caption': 'RunAdel',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x085CE8F0> : "RunAdel"',
         'hide': True,
         'id': 5,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 83.26019148398382,
         'posy': -26.99086599021677,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   6: {  'block': False,
         'caption': 'setAdel',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x096A7DB0> : "setAdel"',
         'hide': True,
         'id': 6,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 237.1974087457034,
         'posy': -57.873378051322895,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   7: {  'block': False,
         'caption': 'devCsv',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x096A7D90> : "devCsv"',
         'hide': True,
         'id': 7,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 139.48287165494958,
         'posy': -121.12217756718923,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   8: {  'block': False,
         'caption': 'geoLeaf',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x096A7D30> : "geoLeaf"',
         'hide': True,
         'id': 8,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 298.93212587427024,
         'posy': -124.4676999108874,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   9: {  'block': False,
         'caption': 'geoAxe',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x096A7E10> : "geoAxe"',
         'hide': True,
         'id': 9,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 376.26331912335615,
         'posy': -126.41381917635147,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   10: {  'block': False,
          'caption': 'axeTCa0N.csv',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x09579570> : "axeTCa0N.csv"',
          'hide': True,
          'id': 10,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': -15.62384482950334,
          'posy': -265.98501361432903,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   11: {  'block': False,
          'caption': 'dimTCa0N.csv',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x09579650> : "dimTCa0N.csv"',
          'hide': True,
          'id': 11,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': 112.20137920192303,
          'posy': -266.0002373485328,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   12: {  'block': False,
          'caption': 'phenTCa0N.csv',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x09579A50> : "phenTCa0N.csv"',
          'hide': True,
          'id': 12,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': 236.5103345145268,
          'posy': -267.6059088483666,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   13: {  'block': False,
          'caption': 'earTCa0N.csv',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x095795D0> : "earTCa0N.csv"',
          'hide': True,
          'id': 13,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': 373.18705010731117,
          'posy': -272.2350136143291,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   14: {  'block': False,
          'caption': 'ssi2sen.csv',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x095797D0> : "ssi2sen.csv"',
          'hide': True,
          'id': 14,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': 494.9961247671154,
          'posy': -273.29865539439436,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   15: {  'block': False,
          'caption': 'freeGeoLeaf',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x096A7B70> : "freeGeoLeaf"',
          'hide': True,
          'id': 15,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': 439.43705010731117,
          'posy': -179.7350136143291,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   16: {  'factory': '<openalea.core.node.NodeFactory object at 0x0847F1B0> : "annotation"',
          'id': 16,
          'posx': 36.93705010731118,
          'posy': -310.9850136143291,
          'txt': u'Csv Table describing plant development'},
   17: {  'block': False,
          'caption': '3',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x090BC050> : "int"',
          'hide': True,
          'id': 17,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': 496.79182696755316,
          'posy': -62.98501361432909,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   20: {  'block': False,
          'caption': '1000',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x090BC050> : "int"',
          'hide': True,
          'id': 20,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': 17.384370821178123,
          'posy': -95.15212066319657,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   21: {  'block': False,
          'caption': 'leaves_simple.db',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x09579A90> : "leaves_simple.db"',
          'hide': True,
          'id': 21,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': -215.58380639573028,
          'posy': -56.92052441475006,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   '__in__': {  'block': False,
                'caption': 'In',
                'delay': 0,
                'hide': True,
                'id': 0,
                'lazy': True,
                'port_hide_changed': set(),
                'posx': 0,
                'posy': 0,
                'priority': 0,
                'use_user_color': True,
                'user_application': None,
                'user_color': None},
   '__out__': {  'block': False,
                 'caption': 'Out',
                 'delay': 0,
                 'hide': True,
                 'id': 1,
                 'lazy': True,
                 'port_hide_changed': set(),
                 'posx': 0,
                 'posy': 0,
                 'priority': 0,
                 'use_user_color': True,
                 'user_application': None,
                 'user_color': None}},
                             elt_value={  2: [(0, '1')],
   3: [],
   4: [],
   5: [  (  2,
            "{'endLeaf': 1.6, 'stemLeaf': 1.2, 'epsillon': 1e-06, 'startLeaf': -0.4, 'senescence_leaf_shrink': 0.5}")],
   6: [(4, 'None'), (5, 'None'), (6, 'None'), (7, "'random'")],
   7: [],
   8: [(0, '4'), (1, '60'), (2, '10')],
   9: [  (0, '75'),
         (1, '0.0'),
         (2, '0.0'),
         (3, '0.0'),
         (4, '60'),
         (5, '5'),
         (6, '7')],
   10: [  (0, 'PackageData(alinea.adel.data, axeTCa0N.csv)'),
          (1, 'None'),
          (2, 'None')],
   11: [  (0, 'PackageData(alinea.adel.data, dimTCa0N.csv)'),
          (1, 'None'),
          (2, 'None')],
   12: [  (0, 'PackageData(alinea.adel.data, phenTCa0N.csv)'),
          (1, 'None'),
          (2, 'None')],
   13: [  (0, 'PackageData(alinea.adel.data, earTCa0N.csv)'),
          (1, 'None'),
          (2, 'None')],
   14: [  (0, 'PackageData(alinea.adel.data, ssi2sen.csv)'),
          (1, 'None'),
          (2, 'None')],
   15: [  (  0,
             "'\\n\\n# R code for setting Azim and Lindex freely  (reference to index of leaf db)\\n\\ngeoLeaf <- list(\\n\\n        Azim = function(a,n,ntop) {\\n          ifelse(ntop <= 4,\\n                 180 + 60 * runif(1),\\n                 180 + 20 * runif(1))\\n               },\\n\\n        Lindex = function(a,n,ntop,stage) {\\n                 ntop + 1\\n               }\\n        )\\n\\n\\n'")],
   16: [],
   17: [(0, '3')],
   20: [(0, '1000')],
   21: [  (0, 'PackageData(alinea.adel.data, leaves_simple.db)'),
          (1, 'None'),
          (2, 'None')],
   '__in__': [],
   '__out__': []},
                             elt_ad_hoc={  2: {  'position': [68.86480669392043, 60.841722418900574],
         'useUserColor': False,
         'userColor': None},
   3: {  'position': [-63.51608384390719, 15.377578193788057],
         'useUserColor': False,
         'userColor': None},
   4: {  'position': [68.86480669392043, 145.08410730660907],
         'useUserColor': False,
         'userColor': None},
   5: {  'position': [83.26019148398382, -26.99086599021677],
         'useUserColor': False,
         'userColor': None},
   6: {  'position': [237.1974087457034, -57.873378051322895],
         'useUserColor': False,
         'userColor': None},
   7: {  'position': [139.48287165494958, -121.12217756718923],
         'useUserColor': False,
         'userColor': None},
   8: {  'position': [298.93212587427024, -124.4676999108874],
         'useUserColor': False,
         'userColor': None},
   9: {  'position': [376.26331912335615, -126.41381917635147],
         'useUserColor': False,
         'userColor': None},
   10: {  'position': [-15.62384482950334, -265.98501361432903],
          'useUserColor': False,
          'userColor': None},
   11: {  'position': [112.20137920192303, -266.0002373485328],
          'useUserColor': False,
          'userColor': None},
   12: {  'position': [236.5103345145268, -267.6059088483666],
          'useUserColor': False,
          'userColor': None},
   13: {  'position': [373.18705010731117, -272.2350136143291],
          'useUserColor': False,
          'userColor': None},
   14: {  'position': [494.9961247671154, -273.29865539439436],
          'useUserColor': False,
          'userColor': None},
   15: {  'position': [439.43705010731117, -179.7350136143291],
          'useUserColor': False,
          'userColor': None},
   16: {  'color': None,
          'position': [36.93705010731118, -310.9850136143291],
          'rectP2': (-1, -1),
          'text': u'Csv Table describing plant development',
          'textColor': None,
          'visualStyle': 0},
   17: {  'position': [496.79182696755316, -62.98501361432909],
          'useUserColor': False,
          'userColor': None},
   20: {  'position': [17.384370821178123, -95.15212066319657],
          'useUserColor': False,
          'userColor': None},
   21: {  'position': [-215.58380639573028, -56.92052441475006],
          'useUserColor': False,
          'userColor': None},
   '__in__': {  'position': [0, 0], 'useUserColor': True, 'userColor': None},
   '__out__': {  'position': [0, 0], 'useUserColor': True, 'userColor': None}},
                             lazy=True,
                             eval_algo='LambdaEvaluation',
                             )




_215416464 = CompositeNodeFactory(name='every_hours',
                             description='',
                             category='time',
                             doc='',
                             inputs=[  {  'desc': '', 'interface': None, 'name': 'time_sequence', 'value': None},
   {  'desc': '', 'interface': None, 'name': 'delay', 'value': None},
   {  'desc': '', 'interface': None, 'name': 'weather', 'value': None}],
                             outputs=[  {  'desc': '', 'interface': None, 'name': 'values'},
   {  'desc': '', 'interface': None, 'name': 'delays'}],
                             elt_factory={  4: ('alinea.astk', 'time_control'), 5: ('alinea.astk', 'time filter')},
                             elt_connections={  34429332: ('__in__', 1, 5, 1),
   34429344: ('__in__', 0, 5, 0),
   34429392: (4, 1, '__out__', 1),
   34429404: (4, 0, '__out__', 0),
   34429416: ('__in__', 2, 4, 2),
   34429428: (5, 1, 4, 1),
   34429440: (5, 0, 4, 0)},
                             elt_data={  4: {  'block': False,
         'caption': 'time_control',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x08F20B70> : "time_control"',
         'hide': True,
         'id': 4,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 13.512348605200884,
         'posy': 66.02052065507483,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   5: {  'block': False,
         'caption': 'time filter',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x08F20BD0> : "time filter"',
         'hide': True,
         'id': 5,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': -2.969484967757051,
         'posy': 15.69584911528721,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   '__in__': {  'block': False,
                'caption': 'In',
                'delay': 0,
                'hide': True,
                'id': 0,
                'lazy': True,
                'port_hide_changed': set(),
                'posx': 28.135818488518098,
                'posy': -38.06103006448593,
                'priority': 0,
                'use_user_color': False,
                'user_application': None,
                'user_color': None},
   '__out__': {  'block': False,
                 'caption': 'Out',
                 'delay': 0,
                 'hide': True,
                 'id': 1,
                 'lazy': True,
                 'port_hide_changed': set(),
                 'posx': 33.71948496775704,
                 'posy': 97.54027564092196,
                 'priority': 0,
                 'use_user_color': False,
                 'user_application': None,
                 'user_color': None}},
                             elt_value={  4: [], 5: [], '__in__': [], '__out__': [(0, 'None'), (1, 'None')]},
                             elt_ad_hoc={  4: {'position': [13.512348605200884, 66.02052065507483], 'userColor': None, 'useUserColor': False},
   5: {'position': [-2.969484967757051, 15.69584911528721], 'userColor': None, 'useUserColor': False},
   '__in__': {'position': [28.135818488518098, -38.06103006448593], 'userColor': None, 'useUserColor': False},
   '__out__': {'position': [33.71948496775704, 97.54027564092196], 'userColor': None, 'useUserColor': False}},
                             lazy=True,
                             eval_algo='LambdaEvaluation',
                             )





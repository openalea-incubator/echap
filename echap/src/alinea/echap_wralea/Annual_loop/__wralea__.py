
# This file has been generated at Wed Mar 27 14:53:37 2013

from openalea.core import *


__name__ = 'Alinea.Echap.Annual_loop'

__editable__ = True
__description__ = ''
__license__ = 'CeCILL-C'
__url__ = 'http://openalea.gforge.inria.fr'
__alias__ = []
__version__ = '0.8.0'
__authors__ = ''
__institutes__ = ''
__icon__ = ''


__all__ = ['_156861872']



_156861872 = CompositeNodeFactory(name='annual_loop',
                             description='Echap annual_loop',
                             category='composite',
                             doc='',
                             inputs=[],
                             outputs=[],
                             elt_factory={  2: ('Alep.Test_nodes', 'adel_mtg3'),
   3: ('Alep.Test_nodes', 'leaves_db'),
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
   18: ('Alinea.Pearl', 'Pearl'),
   19: ('Alinea.Echap.Interfaces', 'pesticide_surfacic_decay'),
   20: ('openalea.data structure', 'int'),
   21: ('alinea.adel.data', 'leaves_simple.db'),
   22: ('__my package__', 'update_on_leaves'),
   23: ('Alinea.Echap.Interfaces', 'pesticide_penetrated_decay'),
   24: ('openalea.data structure.dict', 'dict'),
   25: ('openalea.data structure.dict', 'dict'),
   26: ('openalea.data structure', 'int'),
   27: ('openalea.data structure', 'int'),
   28: ('Alinea.SimCycle', 'decay_model')},
                             elt_connections={  33445160: (9, 0, 6, 2),
   33445172: (19, 0, 23, 0),
   33445184: (20, 0, 5, 0),
   33445196: (24, 0, 18, 0),
   33445208: (26, 0, 18, 2),
   33445220: (5, 0, 2, 2),
   33445232: (8, 0, 6, 1),
   33445244: (14, 0, 7, 4),
   33445256: (11, 0, 7, 1),
   33445268: (10, 0, 7, 0),
   33445280: (25, 0, 18, 1),
   33445292: (2, 0, 22, 0),
   33445304: (17, 0, 6, 3),
   33445316: (27, 0, 28, 0),
   33445328: (18, 0, 19, 1),
   33445340: (3, 0, 2, 1),
   33445352: (7, 0, 6, 0),
   33445364: (12, 0, 7, 2),
   33445376: (22, 0, 19, 0),
   33445388: (6, 0, 5, 1),
   33445400: (28, 0, 23, 1),
   33445412: (13, 0, 7, 3),
   33445424: (21, 0, 3, 0),
   33445436: (17, 0, 2, 3)},
                             elt_data={  2: {  'block': False,
         'caption': 'adel_mtg3',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x095B48B0> : "adel_mtg3"',
         'hide': True,
         'id': 2,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': -287.67451839383284,
         'posy': 38.12592054350357,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   3: {  'block': False,
         'caption': 'leaves_db',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x095B4910> : "leaves_db"',
         'hide': True,
         'id': 3,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': -420.05540893166045,
         'posy': -7.338223681608948,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   5: {  'block': False,
         'caption': 'RunAdel',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x085C29D0> : "RunAdel"',
         'hide': True,
         'id': 5,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': -273.27913360376937,
         'posy': -49.706667865613774,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   6: {  'block': False,
         'caption': 'setAdel',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x096ABF70> : "setAdel"',
         'hide': True,
         'id': 6,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': -119.34191634204977,
         'posy': -80.58917992671991,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   7: {  'block': False,
         'caption': 'devCsv',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x096ABF50> : "devCsv"',
         'hide': True,
         'id': 7,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': -217.0564534328036,
         'posy': -143.83797944258623,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   8: {  'block': False,
         'caption': 'geoLeaf',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x096ABEF0> : "geoLeaf"',
         'hide': True,
         'id': 8,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': -57.60719921348293,
         'posy': -147.1835017862844,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   9: {  'block': False,
         'caption': 'geoAxe',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x096ABFD0> : "geoAxe"',
         'hide': True,
         'id': 9,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 19.723994035602985,
         'posy': -149.12962105174847,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   10: {  'block': False,
          'caption': 'axeTCa0N.csv',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x09540630> : "axeTCa0N.csv"',
          'hide': True,
          'id': 10,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': -372.16316991725654,
          'posy': -288.70081548972604,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   11: {  'block': False,
          'caption': 'dimTCa0N.csv',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x09540710> : "dimTCa0N.csv"',
          'hide': True,
          'id': 11,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': -244.33794588583015,
          'posy': -288.7160392239298,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   12: {  'block': False,
          'caption': 'phenTCa0N.csv',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x09540B10> : "phenTCa0N.csv"',
          'hide': True,
          'id': 12,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': -120.02899057322638,
          'posy': -290.3217107237636,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   13: {  'block': False,
          'caption': 'earTCa0N.csv',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x09540690> : "earTCa0N.csv"',
          'hide': True,
          'id': 13,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': 16.647725019557996,
          'posy': -294.9508154897261,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   14: {  'block': False,
          'caption': 'ssi2sen.csv',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x09540890> : "ssi2sen.csv"',
          'hide': True,
          'id': 14,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': 138.4567996793622,
          'posy': -296.01445726979136,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   15: {  'block': False,
          'caption': 'freeGeoLeaf',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x096ABD30> : "freeGeoLeaf"',
          'hide': True,
          'id': 15,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': 82.89772501955795,
          'posy': -202.4508154897261,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   16: {  'factory': '<openalea.core.node.NodeFactory object at 0x08485290> : "annotation"',
          'id': 16,
          'posx': -319.6022749804421,
          'posy': -333.7008154897261,
          'txt': u'Csv Table describing plant development'},
   17: {  'block': False,
          'caption': '3',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x090D1050> : "int"',
          'hide': True,
          'id': 17,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': 140.25250187979995,
          'posy': -85.7008154897261,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   18: {  'block': False,
          'caption': 'Pearl',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x09CF03F0> : "Pearl"',
          'hide': True,
          'id': 18,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': -67.04985680711181,
          'posy': 73.94612000473853,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   19: {  'block': False,
          'caption': 'pesticide_surfacic_decay',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x097DB2B0> : "pesticide_surfacic_decay"',
          'hide': True,
          'id': 19,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': -176.09393364868822,
          'posy': 131.63334504031815,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   20: {  'block': False,
          'caption': '1000',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x090D1050> : "int"',
          'hide': True,
          'id': 20,
          'is_in_error_state': False,
          'is_user_application': False,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': -339.15495426657515,
          'posy': -117.8679225385936,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   21: {  'block': False,
          'caption': 'leaves_simple.db',
          'delay': 0,
          'factory': '<openalea.core.data.DataFactory object at 0x09540B50> : "leaves_simple.db"',
          'hide': True,
          'id': 21,
          'lazy': True,
          'port_hide_changed': set([2]),
          'posx': -572.1231314834836,
          'posy': -79.63632629014707,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   22: {  'block': False,
          'caption': 'update_on_leaves',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x094F3FB0> : "update_on_leaves"',
          'hide': True,
          'id': 22,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': -230.72839089230536,
          'posy': 79.89400998061915,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   23: {  'block': False,
          'caption': 'pesticide_penetrated_decay',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x097DB210> : "pesticide_penetrated_decay"',
          'hide': True,
          'id': 23,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': -100.55625394112413,
          'posy': 181.13900538709348,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   24: {  'block': False,
          'caption': 'dict',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x097C2A30> : "dict"',
          'hide': True,
          'id': 24,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': -128.5419294368671,
          'posy': 27.853357339478855,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   25: {  'block': False,
          'caption': 'dict',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x097C2A30> : "dict"',
          'hide': True,
          'id': 25,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': -98.79587439622424,
          'posy': 17.977783447434774,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   26: {  'block': False,
          'caption': '1',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x090D1050> : "int"',
          'hide': True,
          'id': 26,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': 14.484175041664678,
          'posy': 30.285093268935235,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   27: {  'block': False,
          'caption': '1',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x090D1050> : "int"',
          'hide': True,
          'id': 27,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': 67.80135114871125,
          'posy': 83.32214237552466,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   28: {  'block': False,
          'caption': 'decay_model',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x085D9E50> : "decay_model"',
          'hide': True,
          'id': 28,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': -14.484175041664663,
          'posy': 131.674318560588,
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
   5: [  (  2,
            "{'endLeaf': 1.6, 'stemLeaf': 1.2, 'epsillon': 1e-06, 'startLeaf': -0.4, 'senescence_leaf_shrink': 0.5}")],
   6: [(4, 'None'), (5, 'None'), (6, 'None'), (7, "'random'")],
   7: [(1, 'None')],
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
   18: [],
   19: [(2, "'LeafElement'")],
   20: [(0, '1000')],
   21: [  (0, 'PackageData(alinea.adel.data, leaves_simple.db)'),
          (1, 'None'),
          (2, 'None')],
   22: [(1, "'LeafElement'")],
   23: [(2, "'LeafElement'")],
   24: [(0, '{}')],
   25: [(0, '{}')],
   26: [(0, '1')],
   27: [(0, '1')],
   28: [],
   '__in__': [],
   '__out__': []},
                             elt_ad_hoc={  2: {'position': [-287.67451839383284, 38.12592054350357], 'userColor': None, 'useUserColor': False},
   3: {'position': [-420.05540893166045, -7.338223681608948], 'userColor': None, 'useUserColor': False},
   4: {  'position': [-18.59601956445446, 132.23836134723174],
         'useUserColor': False,
         'userColor': None},
   5: {'position': [-273.27913360376937, -49.706667865613774], 'userColor': None, 'useUserColor': False},
   6: {'position': [-119.34191634204977, -80.58917992671991], 'userColor': None, 'useUserColor': False},
   7: {'position': [-217.0564534328036, -143.83797944258623], 'userColor': None, 'useUserColor': False},
   8: {'position': [-57.60719921348293, -147.1835017862844], 'userColor': None, 'useUserColor': False},
   9: {'position': [19.723994035602985, -149.12962105174847], 'userColor': None, 'useUserColor': False},
   10: {'position': [-372.16316991725654, -288.70081548972604], 'userColor': None, 'useUserColor': False},
   11: {'position': [-244.33794588583015, -288.7160392239298], 'userColor': None, 'useUserColor': False},
   12: {'position': [-120.02899057322638, -290.3217107237636], 'userColor': None, 'useUserColor': False},
   13: {'position': [16.647725019557996, -294.9508154897261], 'userColor': None, 'useUserColor': False},
   14: {'position': [138.4567996793622, -296.01445726979136], 'userColor': None, 'useUserColor': False},
   15: {'position': [82.89772501955795, -202.4508154897261], 'userColor': None, 'useUserColor': False},
   16: {'visualStyle': 0, 'position': [-319.6022749804421, -333.7008154897261], 'color': None, 'text': u'Csv Table describing plant development', 'textColor': None, 'rectP2': (-1, -1)},
   17: {'position': [140.25250187979995, -85.7008154897261], 'userColor': None, 'useUserColor': False},
   18: {'position': [-67.04985680711181, 73.94612000473853], 'userColor': None, 'useUserColor': False},
   19: {'position': [-176.09393364868822, 131.63334504031815], 'userColor': None, 'useUserColor': False},
   20: {'position': [-339.15495426657515, -117.8679225385936], 'userColor': None, 'useUserColor': False},
   21: {'position': [-572.1231314834836, -79.63632629014707], 'userColor': None, 'useUserColor': False},
   22: {'position': [-230.72839089230536, 79.89400998061915], 'userColor': None, 'useUserColor': False},
   23: {'position': [-100.55625394112413, 181.13900538709348], 'userColor': None, 'useUserColor': False},
   24: {'position': [-128.5419294368671, 27.853357339478855], 'userColor': None, 'useUserColor': False},
   25: {'position': [-98.79587439622424, 17.977783447434774], 'userColor': None, 'useUserColor': False},
   26: {'position': [14.484175041664678, 30.285093268935235], 'userColor': None, 'useUserColor': False},
   27: {'position': [67.80135114871125, 83.32214237552466], 'userColor': None, 'useUserColor': False},
   28: {'position': [-14.484175041664663, 131.674318560588], 'userColor': None, 'useUserColor': False},
   '__in__': {'position': [0, 0], 'userColor': None, 'useUserColor': True},
   '__out__': {'position': [0, 0], 'userColor': None, 'useUserColor': True}},
                             lazy=True,
                             eval_algo='LambdaEvaluation',
                             )





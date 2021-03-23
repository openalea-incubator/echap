""" Reconstruct a MTG wheat canopy at different treatment dates, save dye interception
    and compare it with field data """
import os
import glob
import pandas
import numpy

from alinea.echap.evaluation_canopy import illuminate_canopies, get_light
from alinea.echap.cache_simulation import get_canopy, cache_simulation_path, \
    get_midribs


def cache_dye_path(tag, rep):
    sim_path = cache_simulation_path(tag, rep)
    path = sim_path / 'dye_interception'
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def dye_aggregation_types(what=None):
    """ aggregating functions for one simulation
    """

    def first_val(x):
        if all(pandas.isnull(x)):
            return x.values[0]
        else:
            return x[x.first_valid_index()]

    types = {'nb_plantes_sim': first_val, 'domain_area': first_val,
             'age': first_val, 'nff': first_val, 'area': numpy.sum,
             'green_area': numpy.sum, 'senesced_area': numpy.sum,
             'id': lambda x: '_'.join(map(str, x)), 'length': numpy.sum,
             'ntop': first_val, 'organ': first_val, 'mature_length': first_val,
             'surfacic_doses_Tartrazine': numpy.mean,
             'deposit_Tartrazine': numpy.sum, 'lifetime': first_val,
             'exposition': first_val, 'hasEar': first_val,
             'light_interception': numpy.mean, 'treatment': first_val}
    if what is None:
        return types
    else:
        return {w: types.get(w, first_val) for w in what}


class LeafElementRecorder:
    def __init__(self):
        self.data = {}
        self.counts = 0

    def record_mtg(self, g, daydate, header={}, label='LeafElement'):
        """
        tentative protocol for recording data during a simulation
        ng tentative debug date 12/12/12
        """
        print(daydate)
        for vid in g:
            n = g.node(vid)
            if n.label is not None:
                if n.label.startswith(label):
                    blade = n.complex()
                    if blade.visible_length > 0:
                        axe = n.complex().complex().complex()
                        header.update({'daydate': daydate,
                                       'plant': n.complex().complex().complex().complex().label,
                                       'axe': axe.label, 'nff': axe.nff,
                                       'hasEar': axe.hasEar, 'metamer': int(
                                ''.join(list(n.complex().complex().label)[7:])),
                                       'organ': n.complex().label,
                                       'ntop': n.complex().ntop,
                                       'mature_length': n.complex().shape_mature_length,
                                       'id': n._vid})
                        self.record(n, header)

    def record(self, node, header):
        if node.area > 0:
            data = {}
            properties = node.properties()
            items = ['length', 'area', 'green_area', 'senesced_area',
                     'light_interception']
            for item in items:
                if item in properties:
                    data[item] = properties[item]

            organ = node.complex()
            properties = organ.properties()
            items = ['exposition', 'lifetime', 'age', 'is_ligulated']
            for item in items:
                if item in properties:
                    data[item] = properties[item]

            data.update(header)
            self.data[self.counts] = data
            self.counts += 1

    def save_records(self, path='./echap_outputs'):
        d = pandas.DataFrame(self.data)
        d.T.to_csv(path)

    def get_records(self):
        d = pandas.DataFrame(self.data)
        df = d.T
        df = df.convert_objects()
        return df


def dye_interception_canopies(variety='Tremie12', nplants=30, tag='reference',
                              rep=1, start=None, stop=None, by=None,
                              at=('T1', 'T2'), reset=False, reset_light=False,
                              reset_build=False, reset_reconstruction=False):


    dd_range = illuminate_canopies(light_tag='spray', z_soil=0, variety=variety,
                                   nplants=nplants, tag=tag, rep=rep,
                                   start=start, stop=stop, by=by, at=at,
                                   reset=reset_light, reset_build=reset_build,
                                   reset_reconstruction=reset_reconstruction)
    head_path = cache_dye_path(tag, rep) / variety.lower() + '_' + str(
        nplants) + 'pl_'

    pattern = head_path + '*.csv'
    done = glob.glob(pattern)
    if len(done) > 0:
        done = [x.split('pl_')[1].split('.')[0] for x in done]
        done = ['-'.join(x.split('_')) for x in done]

    missing = dd_range
    if not reset:
        missing = [d for d in dd_range if d not in done]
        treatments = [tr for tr,d in zip(at, dd_range) if d not in done]

    if len(missing) > 0:
        for daydate, treatment in zip(missing, treatments):
            filename = head_path + '_'.join(daydate.split('-')) + '.csv'
            g = get_canopy(variety=variety, nplants=nplants, daydate=daydate,
                           tag=tag, rep=rep)
            _, light, _ = get_light(light_tag='spray', z_soil=0,
                                    variety=variety, nplants=nplants,
                                    daydate=daydate, tag=tag, rep=rep)
            midribs = get_midribs(variety=variety, nplants=nplants,
                                  daydate=daydate, tag=tag, rep=rep)

            g.add_property('light_interception')
            g.property('light_interception').update(light['Ei'])
            recorder = LeafElementRecorder()
            recorder.record_mtg(g, daydate)
            df = recorder.get_records()
            meta = list(g.property('meta').values())[0]
            # add columns
            df['treatment'] = treatment
            df['surfacic_doses_Tartrazine'] = df['light_interception']
            df['deposit_Tartrazine'] = df['light_interception'] * df['area']
            df['nb_plantes_sim'] = meta['nplants']
            df['domain_area'] = meta['domain_area'] * 10000  # m2 -> cm2
            # aggregate by metamer
            gr = df.groupby(['daydate', 'plant', 'axe', 'metamer'],
                            as_index=False)
            df_leaf = gr.agg(dye_aggregation_types())
            # add midribs
            midribs = midribs.rename(columns={'leaf': 'metamer'})
            midribs['plant'] = ['plant' + str(p) for p in midribs['plant']]
            df_i = df_leaf.merge(midribs)
            df_i.to_csv(filename, index=False)

    return dd_range


def get_dye_interception(variety='Tremie12', nplants=30, daydate='T1',
                         tag='reference', rep=1, reset=False, reset_build=False,
                         reset_reconstruction=False):

    dd_range = dye_interception_canopies(variety=variety, nplants=nplants,
                                         tag=tag, rep=rep, at=[daydate],
                                         reset=reset, reset_build=reset_build,
                                         reset_reconstruction=reset_reconstruction)

    daydate = dd_range[0]
    filename = cache_dye_path(tag, rep) / variety.lower() + '_' + str(
        nplants) + 'pl_' + '_'.join(daydate.split('-')) + '.csv'
    return pandas.read_csv(filename)


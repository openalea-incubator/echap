""" Reconstruct a MTG wheat canopy at different treatment dates, save dye interception 
    and compare it with field data """

import pandas
import numpy
import os
from copy import deepcopy

from alinea.echap.conf_int import conf_int
from alinea.echap.hs_tt import tt_hs_tag, as_daydate
from alinea.echap.evaluation_canopy import canopy_illumination
import alinea.echap.interception_data as idata
from alinea.echap.interception_leaf import InterceptModel, \
    pesticide_applications
from alinea.echap.interfaces import pesticide_interception
from alinea.echap.cache_simulation import cache_analysis_path, get_canopy, \
    get_midribs


def df_interception_path(variety='Tremie12', nplants=30,
                         simulation='reference'):
    filename = 'interception_' + simulation + '_' + variety.lower() + '_' + str(
        nplants) + 'pl.csv'
    return str(cache_analysis_path(simulation) / filename)


def aggregation_types(what=None):
    """ aggregating functions for one simulation
    """

    def first_val(x):
        if all(pandas.isnull(x)):
            return x.values[0]
        else:
            return x[x.first_valid_index()]

    types = {'variety': first_val, 'treatment': first_val,
             'nb_plantes_sim': first_val, 'domain_area': first_val,
             'numero_sim': first_val, 'age': first_val,
             'leaf_emergence': first_val, 'TT': first_val, 'nff': first_val,
             'area': numpy.sum, 'HS': numpy.mean, 'green_area': numpy.sum,
             'senesced_area': numpy.sum,
             'id': lambda (x): '_'.join(map(str, x)), 'length': numpy.sum,
             'ntop': first_val, 'ntop_cur': first_val, 'ntop_lig': first_val,
             'organ': first_val, 'mature_length': first_val,
             'surfacic_doses_Tartrazine': numpy.mean,
             'deposit_Tartrazine': numpy.sum, 'lifetime': first_val,
             'exposition': first_val, 'hasEar':first_val,
             'light_interception':numpy.mean, 'tag_T1': first_val, 'tag_T2': first_val,
             'n_max': first_val, 'nflig': first_val}
    if what is None:
        return types
    else:
        return {w: types[w] for w in what}


def aggregate_by_leaf(df):
    """
    Aggregate interceptioin data by leaf and add colmun 'deposits' (= area * surfacic doses)
    """

    df = df[(df['hasEar'] == 1)]  # do not consider aborting axes

    df = df.convert_objects()
    gr = df.groupby(['daydate', 'plant', 'axe', 'metamer'],
                            as_index=False)
    df_leaf = gr.agg(aggregation_types())
    # strange transforms ?
    # df = df.rename(columns = {'TT':'degree_days', 'plant':'num_plant',
    #                              'nff':'fnl', 'axe':'axis', 'ntop':'num_leaf_top'})
    # df['date'] = df['date'].apply(lambda x: pandas.to_datetime(x, dayfirst=True))
    # df['num_leaf_bottom'] = df['fnl'] - df['num_leaf_top'] + 1
    return df_leaf


def aggregate_by_axe(df_leaf):
    """
    Aggregate  leaf-aggregated interceptioin data by axe
    """

    res = []
    what = ['variety', 'treatment', 'nb_plantes_sim', 'TT', 'HS', 'nff', 'area',
            'daydate', 'green_area', 'senesced_area', 'organ',
            'surfacic_doses_Tartrazine', 'deposit_Tartrazine']

    for name, gr in df_leaf.groupby(['daydate', 'plant', 'axe', 'numero_sim'],
                                    as_index=False):
        df_agg = gr.groupby(['daydate', 'plant', 'axe', 'numero_sim'],
                            as_index=False).agg(aggregation_types(what))
        gr = gr.sort('metamer')
        frac = gr['length'] / gr['mature_length']
        ilig = numpy.max(numpy.where(frac >= 1))
        df_agg['haun_stage'] = gr['metamer'].values[ilig] + frac[
                                                            (ilig + 1):].sum()
        df_agg['leaves'] = '_'.join(map(str, gr['metamer'].drop_duplicates()))
        res.append(df_agg)
    return pandas.concat(res)


def aggregate_by_plant(df_leaf):
    """
    Aggregate  leaf-aggregated interceptioin data by axe
    """

    what = ['variety', 'treatment', 'nb_plantes_sim', 'TT', 'nff', 'area',
            'daydate', 'green_area', 'senesced_area', 'organ', 'domain_area',
            'surfacic_doses_Tartrazine', 'deposit_Tartrazine']
    agg = df_leaf.groupby(['daydate', 'plant', 'numero_sim'], as_index=False).agg(
        aggregation_types(what))
    plant_domain = agg['domain_area'] / agg['nb_plantes_sim']
    agg['fraction_intercepted'] = agg['deposit_Tartrazine'] / plant_domain
    agg['lai'] = agg['area'] / plant_domain
    agg['glai'] = agg['green_area'] / plant_domain
    return agg


def axis_statistics(df_sim, what='deposit_Tartrazine', err=conf_int, axis='MS'):
    data = aggregate_by_axe(df_sim)

    if axis == 'MS':
        data = data[data['axe'] == 'MS']

    sub = data.ix[:, ['treatment', what]]
    agg = sub.groupby('treatment').mean().reset_index()
    errag = sub.groupby('treatment').agg(err).reset_index()
    agg['ybar'] = agg[what]
    agg['yerr'] = errag[what]
    agg = agg.set_index('treatment')
    return agg


def plant_statistics(df_sim, what='deposit_Tartrazine', err=conf_int):
    data = aggregate_by_plant(df_sim)

    sub = data.ix[:, ['treatment', what]]
    agg = sub.groupby('treatment').mean().reset_index()
    errag = sub.groupby('treatment').agg(err).reset_index()
    agg['ybar'] = agg[what]
    agg['yerr'] = errag[what]
    agg = agg.set_index('treatment')
    return agg


def leaf_statistics(df_sim, what='deposit_Tartrazine', err=conf_int,
                    by='ntop_cur', axis='MS'):
    data = df_sim
    if axis == 'MS':
        data = data[data['axe'] == 'MS']
    if not isinstance(what, list):
        what = [what]
    sub = data.ix[:, ['treatment', by] + what]
    agg = sub.groupby(['treatment', by]).mean().reset_index()
    errag = sub.groupby(['treatment', by]).agg(err).reset_index()
    agg['ybar'] = agg[what[0]]
    agg['yerr'] = errag[what[0]]
    if by == 'ntop_cur':
        agg['xbar'] = ['F' + str(int(i)) for i in agg[by]]
    elif by == 'ntop_lig':
        agg['xbar'] = ['Fl' + str(int(i)) for i in agg[by]]
    else:
        agg['xbar'] = ['L' + str(int(leaf)) for leaf in agg['metamer']]
    agg = agg.set_index('treatment')
    return agg


def pdict(value):
    """ create a parameter dict for all echap cultivar with value
    """
    return {k: deepcopy(value) for k in
            ('Mercia', 'Rht3', 'Tremie12', 'Tremie13')}


def simulation_tags():
    tags = {'reference': {'dose': pdict(1e4), 'reconstruction_pars': None}, }
    for shape in (
    'MerciaRht', 'Tremie', 'Soissons', 'Tremie12', 'Tremie13', 'Mercia11',
    'Rht311'):
        tags['shape_' + shape + '_byleafclass'] = {'dose': pdict(1e4),
                                                   'reconstruction_pars': {
                                                       'xy_data': pdict(
                                                           shape + '_byleafclass'),
                                                       'top_leaves': pdict(4)}}
        tags['shape_' + shape] = {'dose': pdict(1e4), 'reconstruction_pars': {
            'xy_data': pdict(shape), 'top_leaves': pdict(0)}}
    return tags


class LeafElementRecorder:
    def __init__(self):
        self.data = {}
        self.counts = 0

    def record_mtg(self, g, daydate, header={}, label='LeafElement'):
        """
        tentative protocol for recording data during a simulation
        ng tentative debug date 12/12/12
        """
        print daydate
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

            # items = ['surfacic_doses', 'penetrated_doses', 'global_efficacy']
            # for item in items:
            #     if item in properties:
            #         for compound in properties[item]:
            #             data['_'.join([item, compound])] = properties[item][
            #                 compound]
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
        return d.T


def get_sim(variety='Tremie12', nplants=30, daydate='T1', tag='reference',
            rep=1, reset=False, reset_reconstruction=False):
    treatment=daydate
    tths = tt_hs_tag(variety, tag)
    daydate = as_daydate(daydate, tths)

    g = get_canopy(variety=variety, nplants=nplants, daydate=daydate, tag=tag,
                   rep=rep, reset=reset,
                   reset_reconstruction=reset_reconstruction)
    _, light, _ = canopy_illumination(variety=variety, nplants=nplants,
                                      daydate=daydate, tag=tag, rep=rep,
                                      light_tag='spray', z_soil=0, reset=reset)
    g.add_property('light_interception')
    g.property('light_interception').update(light['Ei'])
    recorder = LeafElementRecorder()
    recorder.record_mtg(g, daydate)
    midribs = get_midribs(variety=variety, nplants=nplants, daydate=daydate,
                          tag=tag, rep=rep, reset=reset)

    df = recorder.get_records()
    df = df.merge(tths)
    meta = g.property('meta').values()[0]
    # add columns
    df['surfacic_doses_Tartrazine'] = df['light_interception']
    df['deposit_Tartrazine'] = df['light_interception'] * df['area']
    df['variety'] = variety
    df['treatment'] = treatment
    df['nb_plantes_sim'] = meta['nplants']
    df['domain_area'] = meta['domain_area'] * 10000  # m2 -> cm2
    df['numero_sim'] = rep
    # compute n_max, ntop_cur
    gr = df.groupby(['daydate', 'plant', 'axe'], group_keys=False)

    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        if 'is_ligulated' in sub.columns:
            sub['nflig'] = sub['metamer'][sub['is_ligulated'] > 0].max()
        else:
            sub['nflig'] = numpy.nan
        return sub

    df = gr.apply(_fun)
    df['ntop_cur'] = df['n_max'] - df['metamer'] + 1
    df['ntop_lig'] = df['nflig'] - df['metamer'] + 1
    # compute leaf emergence
    df['leaf_emergence'] = df['TT'] - df['age']

    return g, df, midribs


def run_sim(variety='Tremie12', nplants=30, daydate='T1', tag='reference',
            rep=1, reset=False, reset_reconstruction=False):
    g, df, midribs = get_sim(variety=variety, nplants=nplants, daydate=daydate,
                             tag=tag, rep=rep, reset=reset,
                             reset_reconstruction=reset_reconstruction)
    df_i = aggregate_by_leaf(df)
    midribs = midribs.rename(columns={'leaf': 'metamer'})
    midribs['plant'] = ['plant' + str(p) for p in midribs['plant']]
    # add domain area
    return df_i.merge(midribs)


def dye_interception(variety='Tremie12', nplants=30, nrep=1,
                     simulation='reference', treatments=None, reset=False,
                     reset_reconstruction=False):
    path = df_interception_path(variety=variety, nplants=nplants,
                                simulation=simulation)
    repetitions = range(1, nrep + 1)
    if treatments is None:
        treatments = idata.tag_treatments()[variety]['application']
    new_sim = False
    dfint = []

    if os.path.exists(path):
        df_old = pandas.read_csv(path)
        done_reps = list(set(df_old['numero_sim']))
        missing = [rep for rep in repetitions if not rep in done_reps]
        # check all treatments are there
        tocheck = df_old.groupby('numero_sim')
        for i_sim, df in tocheck:
            done = df['treatment'].drop_duplicates().values
            dfint.append(df)
            to_do = [t for t in treatments if not t in done]
            if len(to_do) > 0:
                new_sim = True
                for d in to_do:
                    df_t = run_sim(variety=variety, nplants=nplants, daydate=d,
                                   tag=simulation, rep=i_sim, reset=reset,
                                   reset_reconstruction=reset_reconstruction)
                    dfint.append(df_t)
    else:
        missing = repetitions

    if len(missing) > 0:
        new_sim = True
    for i in missing:
        for d in treatments:
            df_t = run_sim(variety=variety, nplants=nplants, daydate=d,
                           tag=simulation, rep=i, reset=reset,
                           reset_reconstruction=reset_reconstruction)
            dfint.append(df_t)

    df_interception = pandas.concat(dfint).reset_index(drop=True)
    if new_sim:
        df_interception.to_csv(path, index=False)

    df_interception = df_interception[
        df_interception['numero_sim'].isin(repetitions)]
    df_interception = df_interception[
        df_interception['treatment'].isin(treatments)]
    return df_interception

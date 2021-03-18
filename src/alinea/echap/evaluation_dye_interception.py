""" Reconstruct a MTG wheat canopy at different treatment dates, save dye interception 
    and compare it with field data """

import pandas
import numpy
# from copy import deepcopy

from alinea.echap.conf_int import conf_int
from alinea.echap.hs_tt import tt_hs_tag

from alinea.echap.cache_simulation import cache_analysis_path
from alinea.echap.cache_dye_simulation import dye_aggregation_types, \
    dye_interception_canopies, get_dye_interception
from alinea.echap.cache_light_simmulation import tag_to_zenith, lambda0
from alinea.echap.evaluation_canopy import simulate_lai


# def pdict(value):
#     """ create a parameter dict for all echap cultivar with value
#     """
#     return {k: deepcopy(value) for k in
#             ('Mercia', 'Rht3', 'Tremie12', 'Tremie13')}


# def simulation_tags():
#     tags = {'reference': {'dose': pdict(1e4), 'reconstruction_pars': None}, }
#     for shape in (
#     'MerciaRht', 'Tremie', 'Soissons', 'Tremie12', 'Tremie13', 'Mercia11',
#     'Rht311'):
#         tags['shape_' + shape + '_byleafclass'] = {'dose': pdict(1e4),
#                                                    'reconstruction_pars': {
#                                                        'xy_data': pdict(
#                                                            shape + '_byleafclass'),
#                                                        'top_leaves': pdict(4)}}
#         tags['shape_' + shape] = {'dose': pdict(1e4), 'reconstruction_pars': {
#             'xy_data': pdict(shape), 'top_leaves': pdict(0)}}
#     return tags


def decorated(df, variety, tag, rep):
    tths = tt_hs_tag(variety, tag)
    df = df.merge(tths)
    df['leaf_emergence'] = df['TT'] - df['age']
    df['variety'] = variety
    df['rep'] = rep
    df['numero_sim'] = rep

    # compute n_max, ntop_cur
    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        if 'is_ligulated' in sub.columns:
            sub['nflig'] = sub['metamer'][sub['is_ligulated'] > 0].max()
        else:
            sub['nflig'] = numpy.nan
        return sub

    df = df.groupby(['daydate', 'plant', 'axe'], group_keys=False).apply(_fun)
    df['ntop_cur'] = df['n_max'] - df['metamer'] + 1
    df['ntop_lig'] = df['nflig'] - df['metamer'] + 1

    return df


def simulate_dye_interception(variety='Tremie12', nplants=30, rep=1,
                              tag='reference', start=None, stop=None, by=None,
                              at=('T1', 'T2'), reset=False, reset_dye=False,
                              reset_build=False, reset_light=False,
                              reset_reconstruction=False):
    sim_path = cache_analysis_path(tag)
    filename = sim_path / 'interception_' + variety.lower() + '_' + str(
        nplants) + 'pl.csv'

    dd_range = dye_interception_canopies(variety=variety, nplants=nplants,
                                         tag=tag, rep=rep, start=start,
                                         stop=stop, by=by, at=at,
                                         reset=reset_dye,
                                         reset_light=reset_light,
                                         reset_build=reset_build,
                                         reset_reconstruction=reset_reconstruction)
    missing = dd_range
    df = None
    if not reset:
        try:
            df = pandas.read_csv(filename)
            try:
                df = df.set_index('rep').loc[rep, :].reset_index()
                if all(map(lambda x: x in df['daydate'].values, dd_range)):
                    return df.loc[df['daydate'].isin(dd_range), :]
                else:
                    missing = [d for d in dd_range if
                               d not in df['daydate'].values]
            except KeyError:
                pass
        except IOError:
            pass
    new = []
    for d in missing:
        print d
        df_i = get_dye_interception(variety=variety, nplants=nplants, tag=tag,
                                    rep=rep, daydate=d)
        new.append(df_i)
    if len(new) > 0:
        df_new = pandas.concat(new)
        df_new = decorated(df_new, variety=variety, tag=tag, rep=rep)
        if df is not None:
            df = pandas.concat((df, df_new))
        else:
            df = df_new
        df.to_csv(filename, index=False)

    df = df[(df['hasEar'] == 1)]  # do not consider aborting axes
    return df.loc[df['daydate'].isin(dd_range), :]


def dye_interception(variety='Tremie12', nplants=30, nrep=1,
                     simulation='reference', treatments=('T1', 'T2'),
                     reset=False, reset_dye=False, reset_build=False,
                     reset_light=False, reset_reconstruction=False):
    repetitions = range(1, nrep + 1)
    reps = []
    for i_sim in repetitions:
        df = simulate_dye_interception(variety=variety, nplants=nplants,
                                       rep=i_sim, tag=simulation, at=treatments,
                                       reset=reset, reset_dye=reset_dye,
                                       reset_build=reset_build,
                                       reset_light=reset_light,
                                       reset_reconstruction=reset_reconstruction)
        reps.append(df)
    return pandas.concat(reps)


def aggregate_by_axe(df_leaf):
    """
    Aggregate  leaf-aggregated interceptioin data by axe
    """

    res = []
    what = ['treatment', 'nb_plantes_sim', 'nff', 'area',
            'green_area', 'senesced_area', 'organ',
            'surfacic_doses_Tartrazine', 'deposit_Tartrazine']

    for name, gr in df_leaf.groupby(['daydate', 'plant', 'axe', 'numero_sim'],
                                    as_index=False):
        df_agg = gr.groupby(['daydate', 'plant', 'axe', 'numero_sim'],
                            as_index=False).agg(dye_aggregation_types(what))
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

    what = ['variety', 'treatment', 'nb_plantes_sim', 'nff', 'area',
            'green_area', 'senesced_area', 'organ', 'domain_area',
            'surfacic_doses_Tartrazine', 'deposit_Tartrazine']
    agg = df_leaf.groupby(['daydate', 'plant', 'numero_sim'], as_index=False).agg(
        dye_aggregation_types(what))
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


def dye_interception_miller(variety='Tremie12', nplants=30, tag='reference', rep=1,
                          at=('T1','T2')):
    df_sim = simulate_dye_interception(variety=variety, nplants=nplants, tag=tag, rep=rep, at=at)


    df_lai = simulate_lai(variety=variety, nplants=nplants, tag=tag, rep=1,
                          at=at)
    zenith = tag_to_zenith('spray', tag, variety)
    l = lambda0(tag, variety)
    k = l * 0.5 / numpy.cos(numpy.radians(zenith))
    sim = leaf_statistics(df_sim, what='area', by='ntop_cur', axis='MS')
    res = []
    for tr, dat in sim.groupby(sim.index):
        dat['deposit_Tartrazine'] = 0
        area_sum = dat['area'].sum()
        if tr.startswith('T1'):
            coltr = 'tag_T1'
        else:
            coltr = 'tag_T2'
        lai = df_lai.set_index(coltr).loc[tr, 'LAI_tot']
        for i in range(len(dat)):
            lai_cum = [0] + dat['area'].cumsum().values.tolist()
            io = numpy.exp(-k * lai_cum[i] / area_sum * lai)
            lai_leaf = dat['area'].values[i] / area_sum * lai
            dat.ix[i, 'deposit_Tartrazine'] = io * (1 - numpy.exp(-k * lai_leaf)) * dat['area'].values[i]
        dat['ybar'] = dat['deposit_Tartrazine']
        res.append(dat)
    return pandas.concat(res)


def dye_interception_miller_no_layer(variety='Tremie12', nplants=30, tag='reference', rep=1,
                          at=('T1','T2')):
    df_sim = simulate_dye_interception(variety=variety, nplants=nplants, tag=tag, rep=rep, at=at)

    df_lai = simulate_lai(variety=variety, nplants=nplants, tag=tag, rep=1,
                          at=at)
    zenith = tag_to_zenith('spray', tag, variety)
    l = lambda0(tag, variety)
    k = l * 0.5 / numpy.cos(numpy.radians(zenith))
    sim = leaf_statistics(df_sim, what='area', by='ntop_cur', axis='MS')
    res = []
    for tr, dat in sim.groupby(sim.index):
        dat['deposit_Tartrazine'] = 0
        area_sum = dat['area'].sum()
        if tr.startswith('T1'):
            coltr = 'tag_T1'
        else:
            coltr = 'tag_T2'
        lai = df_lai.set_index(coltr).loc[tr, 'LAI_tot']
        intercepted = 1 - numpy.exp(-k * lai)
        for i in range(len(dat)):
            lai_leaf = dat['area'].values[i] / area_sum * lai
            dat.ix[i, 'deposit_Tartrazine'] = lai_leaf / lai * intercepted * dat['area'].values[i]
        dat['ybar'] = dat['deposit_Tartrazine']
        res.append(dat)
    return pandas.concat(res)

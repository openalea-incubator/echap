""" Reconstruct a MTG wheat canopy at different ages and save LAI, cover fraction,
    intercepted radiation and compare it with field data """

import pandas
import numpy
import json
import os

try:
    import cPickle as pickle
except ImportError:
    import pickle

from alinea.adel.astk_interface import AdelWheat
from alinea.caribu.CaribuScene import CaribuScene
# Import for wheat reconstruction
from alinea.echap.hs_tt import tt_hs_tag, as_daydate

# from alinea.adel.postprocessing import ground_cover

from alinea.caribu.light import vecteur_direction, diffuse_source

from alinea.echap.canopy_simulation import (cache_simulation_path, get_canopy,
                                            build_canopies)
from alinea.echap.canopy_data import lai_pai_scan, pai57, ground_cover_data, transmittance_data
from alinea.echap.evaluation_canopy_plot import plot_mean


def aggregate_lai(g, axstat):
    colnames = ['aire du plot', 'Nbr.plant.perplot', 'Nbr.axe.tot.m2', 'ThermalTime', 'LAI_tot', 'LAI_vert',
                 'PAI_tot', 'PAI_vert']
    pstat = AdelWheat.plot_statistics(g,axstat)
    if pstat is None:
        return pandas.DataFrame([numpy.nan] * len(colnames), columns = colnames)
    else:
        return pstat.loc[:, colnames]


def get_lai_properties(g):
    df_axstat = AdelWheat.axis_statistics(g)
    df_lai_tot = aggregate_lai(g, df_axstat)
    df_lai_MS = aggregate_lai(g, df_axstat[df_axstat['axe_id'] == 'MS'])
    df_lai_MS.rename(columns={col: col + '_MS' for col in
                              ['LAI_tot', 'LAI_vert', 'PAI_tot', 'PAI_vert',
                               'Nbr.axe.tot.m2']}, inplace=True)
    df_lai_ferti = aggregate_lai(g, df_axstat[df_axstat['has_ear']])
    df_lai_ferti.rename(columns={col: col + '_ferti' for col in
                                 ['LAI_tot', 'LAI_vert', 'PAI_tot', 'PAI_vert',
                                  'Nbr.axe.tot.m2']}, inplace=True)
    df_lai = pandas.merge(df_lai_tot, df_lai_MS)
    df_lai = pandas.merge(df_lai, df_lai_ferti)
    return df_lai


def simulate_lai(variety='Tremie12', nplants=30, tag='reference', rep=1,
                 start=None, stop=None, by=None, at=None, reset=False,
                 reset_build=False, reset_reconstruction=False):
    sim_path = cache_simulation_path(tag, rep)
    filename = sim_path / 'lai_' + variety.lower() + '_' + str(
        nplants) + 'pl.csv'

    if all(map(lambda x: x is None, [start, stop, by, at])):
        try:
            df = pandas.read_csv(filename)
            return df
        except IOError:
            raise ValueError(
                'No simulation found for ' + variety + ' please specify start, stop, (by) or at')

    print('Check/build canopies..')
    dd_range = build_canopies(variety=variety, nplants=nplants, tag=tag,
                              rep=rep, start=start, stop=stop, by=by, at=at,
                              reset=reset_build,
                              reset_reconstruction=reset_reconstruction)

    print('Compute missing plot statistics...')
    missing = dd_range
    if not reset:
        try:
            df = pandas.read_csv(filename)
            if all(map(lambda x: x in df['daydate'].values, dd_range)):
                return df.loc[df['daydate'].isin(dd_range), :]
            else:
                missing = [d for d in dd_range if d not in df['daydate'].values]
        except IOError:
            pass
    new = []
    for d in missing:
        print d
        g = get_canopy(daydate=d, variety=variety, nplants=nplants, tag=tag,
                       rep=rep, load_geom=False)
        df_lai = get_lai_properties(g)
        df_lai['variety'] = variety
        df_lai['daydate'] = d
        new.append(df_lai)

    df_new = pandas.concat(new)
    tths = tt_hs_tag(variety, tag)
    df_new = df_new.merge(tths)
    df = pandas.concat((df, df_new))
    df.to_csv(filename, index=False)

    return df.loc[df['daydate'].isin(dd_range), :]


def compare_lai(variety='Tremie12', nplants=30, tag='reference', rep=1,
                start=None, stop=None, by=None, at=None, reset=False,
                reset_build=False, reset_reconstruction=False):
    dfsim = simulate_lai(variety, nplants, tag, rep, start, stop, by, at, reset,
                         reset_build, reset_reconstruction)
    ax = plot_mean(dfsim, 'LAI_vert', xaxis='HS', marker=' ')
    plot_mean(dfsim, 'LAI_vert_MS', xaxis='HS', color='r', ax=ax, marker=' ')
    dfobs = lai_pai_scan(variety)
    plot_mean(dfobs, 'GLAI_MS_biomass', xaxis='HS', color='r', ax=ax,
              linestyle=' ')
    plot_mean(dfobs, 'GLAI_biomass', xaxis='HS', ax=ax, linestyle=' ')
    if variety == 'Tremie12':
        plot_mean(dfsim, 'LAI_tot', xaxis='HS', color='g', ax=ax, marker=' ',
                  linestyle='--')
        plot_mean(dfobs, 'LAI_biomass', xaxis='HS', color='g', ax=ax,
                  linestyle=' ')
    return ax


def tag_to_light(tag='zenith'):
    if tag == 'zenith':
        return [(1, (0, 0, -1))]
    elif tag == '57.5':
        return [(1. / 24, vecteur_direction(90 - 57.5, az)) for az in
                range(0, 360, 15)]
    elif tag == 'soc':
        return diffuse_source(46)
    elif tag.startswith('lai2000r'):
        zenith = 7, 23, 38, 53, 68
        r = int(tag.split('r')[1]) - 1
        return [(1. / 24, vecteur_direction(90 - zenith[r], az)) for az in
                range(0, 360, 15)]
    else:
        raise ValueError('Unknown light tag: ' + tag)


def tag_to_zenith(tag='zenith'):
    if tag == 'zenith':
        return 0
    elif tag == '57.5':
        return 57.5
    elif tag == 'soc':
        energy, directions = zip(*diffuse_source(46))
        z = numpy.array(zip(*directions)[2])
        return numpy.mean(numpy.degrees(numpy.arccos(-z)))
    elif tag.startswith('lai2000r'):
        zenith = 7, 23, 38, 53, 68
        r = int(tag.split('r')[1]) - 1
        return zenith[r]
    else:
        raise ValueError('Unknown light tag: ' + tag)


def canopy_illumination(variety='Tremie12', nplants=30, daydate='T1',
                        tag='reference', rep=1, light_tag='zenith', z_soil=0,
                        reset=False, reset_build=False,
                        reset_reconstruction=False):

    tths = tt_hs_tag(variety, tag)
    daydate = as_daydate(daydate, tths)
    sim_path = cache_simulation_path(tag, rep)
    if not os.path.exists(sim_path / 'light'):
        os.makedirs(sim_path / 'light')
    filename = sim_path / 'light' / light_tag + '_z' + str(
        z_soil) + '_' + variety.lower() + '_' + str(nplants) + 'pl_' + '_'.join(
        daydate.split('-')) + '.json'

    if not reset:
        try:
            with open(filename, 'r') as input_file:
                cached = json.load(input_file)
                raw = {w: {int(k): v for k, v in cached['raw'][w].iteritems()}
                       for w in cached['raw']}
                aggregated = {w: {int(k): v for k, v
                                  in cached['aggregated'][w].iteritems()} for
                              w in cached['aggregated']}
            return raw, aggregated, cached['soil']
        except IOError:
            pass
    g = get_canopy(variety=variety, nplants=nplants, daydate=daydate, tag=tag,
                   rep=rep, reset=reset_build, reset_reconstruction=reset_reconstruction)
    meta = AdelWheat.meta_informations(g)
    sources = tag_to_light(light_tag)
    c2u = {v: k for k, v in CaribuScene.units.iteritems()}
    units = c2u.get(meta['convUnit'])
    cscene = CaribuScene(g, sources, pattern=meta['domain'], soil_mesh=1,
                         z_soil=z_soil, scene_unit=units)
    raw, aggregated = cscene.run(direct=True, infinite=True, simplify=True)
    soil, _ = cscene.getSoilEnergy()
    with open(filename, 'w') as output_file:
        saved = {'raw': raw, 'aggregated': aggregated, 'soil': soil}
        json.dump(saved, output_file, sort_keys=True, indent=4,
                  separators=(',', ': '))

    return raw, aggregated, soil


def simulate_light(variety='Tremie12', nplants=30, tag='reference', rep=1,
                   light_tag='zenith', z_soil=0, start=None, stop=None, by=None,
                   at=None, reset=False, reset_build=False,
                   reset_reconstruction=False):
    sim_path = cache_simulation_path(tag, rep)
    filename = sim_path / 'light_interception_' + variety.lower() + '_' + str(
        nplants) + 'pl.csv'
    if all(map(lambda x: x is None, [start, stop, by, at])):
        try:
            df = pandas.read_csv(filename)
            return df
        except IOError:
            raise ValueError(
                'No simulation found for ' + variety +
                ' please specify start, stop, (by) or at')

    print('Check/build canopies..')
    dd_range = build_canopies(variety=variety, nplants=nplants, tag=tag,
                              rep=rep, start=start, stop=stop, by=by, at=at,
                              reset=reset_build,
                              reset_reconstruction=reset_reconstruction)

    print('Compute light interception...')
    sim_tag = 'po_' + light_tag + '_z' + str(z_soil)
    df = None
    done = None
    missing = dd_range
    if not reset:
        try:
            df = pandas.read_csv(filename)
            if sim_tag not in df.columns:
                df[sim_tag] = [numpy.nan] * len(df)
            if all(map(lambda x: x in df['daydate'].values, dd_range)):
                done = df.loc[df['daydate'].isin(dd_range), :]
                missing = []
            else:
                done = df.loc[df['daydate'].isin(dd_range), :]
                missing = [d for d in dd_range if
                           d not in done['daydate'].values]
        except IOError:
            pass

    redo = []
    if done is not None:
        if reset:
            redo = done['daydate'].values
        for d in done['daydate'].values:
            if not pandas.notnull(done.set_index('daydate').loc[d, sim_tag]):
                redo.append(d)

    new = []
    for d in missing + redo:
        print d
        raw, aggregated, soil = canopy_illumination(variety=variety,
                                                    nplants=nplants, daydate=d,
                                                    tag=tag, rep=rep,
                                                    light_tag=light_tag,
                                                    z_soil=0, reset=reset)
        if d in redo:
            df.loc[df['daydate'] == d, sim_tag] = soil
        else:
            df_light = pandas.DataFrame(
                {'rep': rep, 'variety': variety, 'daydate': d, sim_tag: soil},
                index=[0])
            new.append(df_light)

    if len(new) > 0:
        df_new = pandas.concat(new)
        tths = tt_hs_tag(variety, tag)
        df_new = df_new.merge(tths)
        if df is not None:
            df = pandas.concat((df, df_new))
        else:
            df = df_new

    if len(redo + new) > 0:
        df['p1_' + light_tag + '_z' + str(z_soil)] = 1 - df[sim_tag]
        g_miller = 0.5 / numpy.cos(numpy.radians(tag_to_zenith(light_tag)))
        df['lai_' + light_tag + '_z' + str(z_soil)] = - numpy.log(
            df[sim_tag]) / g_miller
        rings = ['lai_lai2000r' + str(i) + '_z0' for i in range(1, 6)]
        if all(map(lambda x: x in df.columns.values, rings)):
            w = 0.034, 0.104, 0.16, 0.218, 0.494
            df['lai_lai2000'] = 0
            for i, ring in enumerate(rings):
                df['lai_lai2000'] += w[i] * df[ring]
        df.to_csv(filename, index=False)

    return df.loc[df['daydate'].isin(dd_range), :]


def compare_po(variety='Tremie12', nplants=30, tag='reference', rep=1,
               light_tag='zenith', z_soil=0, start=None, stop=None, by=None,
               at=None, reset=False, reset_build=False,
               reset_reconstruction=False, ax=None, color='b'):
    dfsim = simulate_light(variety=variety, nplants=nplants, tag=tag, rep=rep,
                           light_tag=light_tag, z_soil=z_soil, start=start,
                           stop=stop, by=by, at=at, reset=reset,
                           reset_build=reset_build,
                           reset_reconstruction=reset_reconstruction)
    sim_tag = 'po_' + light_tag + '_z' + str(z_soil)
    ax = plot_mean(dfsim, sim_tag, xaxis='HS', ax=ax, color=color, marker='')
    dfobs = None
    obs = None
    if light_tag in ['zenith', '57.5']:
        dfobs = ground_cover_data(variety, tag, angle=tag_to_zenith(light_tag))
        obs = 'po'
    elif light_tag == 'soc':
        dfobs = transmittance_data(variety, tag,
                                   start=dfsim['daydate'].sort_values().values[0],
                                   stop=dfsim['daydate'].sort_values().values[-1])
        obs = 'po_' + str(z_soil)
        if obs + '_mean' not in dfobs.columns:
            obs = None
    else:
        pass
    if dfobs is not None and obs is not None:
        plot_mean(dfobs, obs, xaxis='HS', ax=ax, color=color, linestyle='')
    ax.set_yscale('log')
    # get_ticks / set_ticks
    return ax


def compare_all_po(variety='Tremie12', nplants=30, tag='reference', rep=1,
                   start=None, stop=None, by=None, at=None, reset=False,
                   reset_build=False, reset_reconstruction=False):
    ax = compare_po(variety=variety, nplants=nplants, tag=tag, rep=rep,
                    light_tag='57.5', z_soil=0, start=start, stop=stop, by=by,
                    at=at, reset=reset, reset_build=reset_build,
                    reset_reconstruction=reset_reconstruction)
    compare_po(variety=variety, nplants=nplants, tag=tag, rep=rep,
               light_tag='zenith', z_soil=0, start=start, stop=stop, by=by, at=at,
               reset=reset, reset_build=reset_build,
               reset_reconstruction=reset_reconstruction, ax=ax, color='r')
    compare_po(variety=variety, nplants=nplants, tag=tag, rep=rep,
               light_tag='soc', z_soil=0, start=start, stop=stop, by=by, at=at,
               reset=reset, reset_build=reset_build,
               reset_reconstruction=reset_reconstruction, ax=ax, color='lightgreen')
    return ax




# Run and save canopy properties ###################################################################


# def draft_TC(g, adel, zenith, rep, scale = 1):
#     echap_top_camera =  {'type':'perspective', 'distance':200.,
#                          'fov':50., 'azimuth':0, 'zenith':zenith}
#     gc, im, box = ground_cover(g, adel.domain, camera=echap_top_camera, image_width = int(2144 * scale),
#                                 image_height = int(1424*scale), getImages=True, replicate=rep)
#     return gc
#
# def get_cover_fraction_properties(g, adel, nplants, scale = 1):
#     df_TC = pandas.DataFrame(columns = ['TCgreen', 'TCsen', 'TCtot',
#                                          'TCgreen_57', 'TCsen_57', 'TCtot_57'])
#     if nplants <= 30:
#         reps = [1,2]
#     else:
#         reps = [1,1]
#
#     for zenith, rep in zip([0,57],reps):
#         dict_TC = draft_TC(g, adel, zenith, rep, scale)
#         if zenith == 0:
#             suffix = ''
#         elif zenith == 57:
#             suffix = '_57'
#         df_TC.loc[0, 'TCgreen'+suffix] = dict_TC['green']
#         df_TC.loc[0, 'TCsen'+suffix] = dict_TC['senescent']
#         df_TC.loc[0, 'TCtot'+suffix] = sum(dict_TC.values())
#         df_TC.loc[0, 'Gapgreen'+suffix] = 1 - df_TC.loc[0, 'TCgreen'+suffix]
#         df_TC.loc[0, 'Gapsen'+suffix] = 1 - df_TC.loc[0, 'TCsen'+suffix]
#         df_TC.loc[0, 'Gaptot'+suffix] = 1 - df_TC.loc[0, 'TCtot'+suffix]
#     return df_TC
#
#
#
# def run_one_simulation(variety = 'Tremie12', nplants = 30, variability_type = None,
#                         age_range = [400., 2600.], time_steps = [20, 100],
#                         scale_povray = 1., z_levels = [0, 5, 20, 25],
#                         reset = False, reset_data = False, only_lai = False):
#     # Temp
#     if variety in ['Tremie12', 'Tremie13']:
#         pars = reconstruction_parameters()
#         # pars['density_tuning'] = pdict(None)
#         # pars['density_tuning']['Tremie12'] = 0.85
#         # pars['density_tuning']['Tremie13'] = 0.85
#         reconst = EchapReconstructions(reset_data=reset_data, pars=pars)
#     else:
#         reconst = echap_reconstructions(reset=reset, reset_data=reset_data)
#     HSconv = reconst.HS_fit[variety]
#     adel = reconst.get_reconstruction(name=variety, nplants = nplants)
#     ages_1 = numpy.arange(age_range[0], age_range[1], time_steps[0])
#     ages_2 = numpy.arange(age_range[0], age_range[1], time_steps[1])
#     for age in numpy.unique(ages_1.tolist()+ages_2.tolist()):
#         # Get canopy properties
#         g = adel.setup_canopy(age=age)
#         if age in ages_1:
#             df_lai = get_lai_properties(g, adel)
#         if age in ages_2:
#             if only_lai==False:
#                 df_cover = get_cover_fraction_properties(g, adel, nplants, scale_povray)
#                 df_radiation = get_radiation_properties(g, adel, z_levels)
#             else:
#                 cols = ['TCgreen', 'TCsen', 'TCtot', 'TCgreen_57', 'TCsen_57', 'TCtot_57',]
#                 df_cover = pandas.DataFrame([[numpy.nan for col in cols]], columns = cols)
#                 cols = ['LightPenetration_%d' %lev for lev in z_levels]
#                 df_radiation = pandas.DataFrame([[numpy.nan for col in cols]], columns = cols)
#
#         # Group in same dataframe
#         if age == age_range[0]:
#             df_prop = pandas.DataFrame(columns = ['variety', 'HS']+
#                                                  [col for col in df_lai.columns]+
#                                                  [col for col in df_cover.columns]+
#                                                  [col for col in df_radiation.columns])
#         if age in ages_1:
#             df_prop.loc[age, df_lai.columns] = df_lai.loc[0, :]
#         if age in ages_2:
#             df_prop.loc[age, df_cover.columns] = df_cover.loc[0, :]
#             df_prop.loc[age, df_radiation.columns] = df_radiation.loc[0, :]
#         df_prop.loc[age, 'variety'] = variety
#         df_prop.loc[age, 'HS'] = HSconv(age)
#     df_prop.reset_index(drop = True, inplace = True)
#     df_prop = df_prop.applymap(lambda x: numpy.float(x) if numpy.isreal(x) else x)
#     return df_prop



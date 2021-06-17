""" Reconstruct a MTG wheat canopy at different ages and save LAI, cover fraction,
    intercepted radiation and compare it with field data """

import pandas
import numpy
import matplotlib.pyplot as plt


try:
    import pickle as pickle
except ImportError:
    import pickle

from alinea.adel.astk_interface import AdelWheat

# Import for wheat reconstruction
from alinea.echap.hs_tt import tt_hs_tag

# from alinea.adel.postprocessing import ground_cover



from alinea.echap.cache_simulation import (cache_analysis_path, get_canopy,
                                           build_canopies, get_reconstruction)
from alinea.echap.cache_light_simmulation import illuminate_canopies, get_light, tag_to_zenith, lambda0
from alinea.echap.canopy_data import lai_pai_scan, pai57, ground_cover_data, transmittance_data
from alinea.echap.interception_data import petri_data

from alinea.echap.plot_evaluation_canopy import plot_mean


def aggregate_lai(adel, g, axstat):
    colnames = ['aire du plot', 'Nbr.plant.perplot', 'Nbr.axe.tot.m2', 'ThermalTime', 'LAI_tot', 'LAI_vert',
                 'PAI_tot', 'PAI_vert']
    pstat = adel.plot_statistics(g,axstat)
    if pstat is None:
        return pandas.DataFrame([numpy.nan] * len(colnames), columns = colnames)
    else:
        return pstat.loc[:, colnames]


def get_lai_properties(adel, g):
    df_axstat = adel.axis_statistics(g)
    df_lai_tot = aggregate_lai(adel, g, df_axstat)
    df_lai_MS = aggregate_lai(adel, g, df_axstat[df_axstat['axe_id'] == 'MS'])
    df_lai_MS.rename(columns={col: col + '_MS' for col in
                              ['LAI_tot', 'LAI_vert', 'PAI_tot', 'PAI_vert',
                               'Nbr.axe.tot.m2']}, inplace=True)
    df_lai_ferti = aggregate_lai(adel, g, df_axstat[df_axstat['has_ear']])
    df_lai_ferti.rename(columns={col: col + '_ferti' for col in
                                 ['LAI_tot', 'LAI_vert', 'PAI_tot', 'PAI_vert',
                                  'Nbr.axe.tot.m2']}, inplace=True)
    df_lai = pandas.merge(df_lai_tot, df_lai_MS)
    df_lai = pandas.merge(df_lai, df_lai_ferti)
    return df_lai


def simulate_lai(variety='Tremie12', nplants=30, tag='reference', rep=1,
                 start=None, stop=None, by=None, at=None, reset=False,
                 reset_build=False, reset_reconstruction=False):

    sim_path = cache_analysis_path(tag)
    filename = sim_path / 'lai_' + variety.lower() + '_' + str(
        nplants) + 'pl.csv'

    print('Check/build canopies..')
    dd_range = build_canopies(variety=variety, nplants=nplants, tag=tag,
                              rep=rep, start=start, stop=stop, by=by, at=at,
                              reset=reset_build,
                              reset_reconstruction=reset_reconstruction)

    print('Compute missing plot statistics...')
    missing = dd_range
    df = None
    if not reset:
        try:
            df = pandas.read_csv(filename)
            try:
                df = df.set_index('rep').loc[rep,:].reset_index()
                if all([x in df['daydate'].values for x in dd_range]):
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
        print(d)
        adel = get_reconstruction(variety=variety, nplants=nplants, tag=tag,
                                  rep=rep)
        g = get_canopy(daydate=d, variety=variety, nplants=nplants, tag=tag,
                       rep=rep, load_geom=False)
        df_lai = get_lai_properties(adel, g)
        df_lai['variety'] = variety
        df_lai['daydate'] = d
        df_lai['rep'] = rep
        new.append(df_lai)
    if len(new) > 0:
        df_new = pandas.concat(new)
        tths = tt_hs_tag(variety, tag)
        df_new = df_new.merge(tths)
        if df is not None:
            df = pandas.concat((df, df_new))
        else:
            df = df_new
        df.to_csv(filename, index=False)

    return df.loc[df['daydate'].isin(dd_range), :]


def compare_lai(variety='Tremie12', nplants=30, tag='reference', rep=1,
                start=None, stop=None, by=None, at=None, reset=False,
                reset_build=False, reset_reconstruction=False):
    dfsim = simulate_lai(variety, nplants, tag, rep, start, stop, by, at, reset,
                         reset_build, reset_reconstruction)
    ax = plot_mean(dfsim, 'LAI_vert', xaxis='HS', color='g', marker='', xlabel='Haun Stage', ylabel='LAI')
    plot_mean(dfsim, 'LAI_vert_MS', xaxis='HS', color='r', ax=ax, marker=' ')
    dfobs = lai_pai_scan(variety)
    plot_mean(dfobs, 'GLAI_MS_biomass', xaxis='HS', color='r', ax=ax,
              linestyle=' ')
    plot_mean(dfobs, 'GLAI_biomass', xaxis='HS', ax=ax, linestyle=' ', color='g')
    if variety == 'Tremie12':
        plot_mean(dfsim, 'LAI_tot', xaxis='HS', color='b', ax=ax, marker=' ',
                  linestyle='--')
        plot_mean(dfobs, 'LAI_biomass', xaxis='HS', color='b', ax=ax,
                  linestyle=' ')
    return ax


def simulate_po_light(light_tag='zenith', z_soil=0, variety='Tremie12',
                      nplants=30, tag='reference', rep=1, start=None, stop=None,
                      by=None, at=None, reset=False, reset_light=False,
                      reset_build=False, reset_reconstruction=False):
    sim_path = cache_analysis_path(tag)
    filename = sim_path / 'light_po_' + variety.lower() + '_' + str(
        nplants) + 'pl.csv'
    sim_tag = 'po_' + light_tag + '_z' + str(z_soil)

    print('Check/illuminate canopies..')
    dd_range = illuminate_canopies(light_tag=light_tag, z_soil=z_soil,
                                   variety=variety, nplants=nplants, tag=tag,
                                   rep=rep, start=start, stop=stop, by=by,
                                   at=at, reset=reset_light,
                                   reset_build=reset_build,
                                   reset_reconstruction=reset_reconstruction)

    print('Compute po_light statistics...')
    df = None
    done = None
    missing = dd_range
    if not reset:
        try:
            df = pandas.read_csv(filename)
            if sim_tag not in df.columns:
                df[sim_tag] = [numpy.nan] * len(df)
            try:
                df = df.set_index('rep').loc[rep, :].reset_index()
                if all([x in df['daydate'].values for x in dd_range]):
                    done = df.loc[df['daydate'].isin(dd_range), :]
                    missing = []
                else:
                    done = df.loc[df['daydate'].isin(dd_range), :]
                    missing = [d for d in dd_range if
                               d not in done['daydate'].values]
            except KeyError:
                pass
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
        print(d)
        raw, aggregated, soil = get_light(variety=variety, nplants=nplants,
                                          daydate=d, tag=tag, rep=rep,
                                          light_tag=light_tag, z_soil=0)
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
        zenith = tag_to_zenith(light_tag, tag, variety)
        g_miller = 0.5 / numpy.cos(numpy.radians(zenith))
        df['lai_' + light_tag + '_z' + str(z_soil)] = - numpy.log(
            df[sim_tag]) / g_miller
        rings = ['lai_lai2000r' + str(i) + '_z0' for i in range(1, 6)]
        if all([x in df.columns.values for x in rings]):
            w = 0.034, 0.104, 0.16, 0.218, 0.494
            df['lai_lai2000'] = 0
            for i, ring in enumerate(rings):
                df['lai_lai2000'] += w[i] * df[ring]
        df.to_csv(filename, index=False)

    return df.loc[df['daydate'].isin(dd_range), :]


def po_miller(light_tag='zenith', variety='Tremie12', nplants=30,
              tag='reference', rep=1, start=None, stop=None, by=None, at=None,
              reset_build=False, reset_reconstruction=False, l=None):
    df = simulate_lai(variety=variety, nplants=nplants, tag=tag, rep=rep,
                       start=start, stop=stop, by=by, at=at,
                       reset_build=reset_build,
                       reset_reconstruction=reset_reconstruction)
    zenith = tag_to_zenith(light_tag, tag, variety)
    if l is None:
        l = lambda0(tag, variety)
    k = l * 0.5 / numpy.cos(numpy.radians(zenith))
    sim_tag = 'po_' + light_tag + '_z0'
    df[sim_tag] = numpy.exp(-k * df['LAI_tot'])
    df['p1_' + light_tag + '_z0'] = 1 - df[sim_tag]
    return df


def compare_po(model='adel', variety='Tremie12', nplants=30, tag='reference', rep=1,
               light_tag='zenith', z_soil=0, start=None, stop=None, by=None,
               at=None, reset=False, reset_light=False, reset_build=False,
               reset_reconstruction=False, ax=None, color='b'):
    if model == 'adel':
        dfsim = simulate_po_light(variety=variety, nplants=nplants, tag=tag, rep=rep,
                               light_tag=light_tag, z_soil=z_soil, start=start,
                               stop=stop, by=by, at=at, reset=reset,
                               reset_light=reset_light,reset_build=reset_build,
                               reset_reconstruction=reset_reconstruction)
    elif model == 'miller':
        dfsim = po_miller(variety=variety, nplants=nplants, tag=tag, rep=rep,
                               light_tag=light_tag, start=start,
                               stop=stop, by=by, at=at, reset_build=reset_build,
                               reset_reconstruction=reset_reconstruction)
    elif model == 'miller_raw':
        dfsim = po_miller(variety=variety, nplants=nplants, tag=tag, rep=rep,
                          light_tag=light_tag, start=start, stop=stop, by=by,
                          at=at, reset_build=reset_build,
                          reset_reconstruction=reset_reconstruction, l=1)
    else:
        raise ValueError('unknown model: ' + model)

    sim_tag = 'po_' + light_tag + '_z' + str(z_soil)
    ax = plot_mean(dfsim, sim_tag, xaxis='HS', ax=ax, color=color, marker='')
    dfobs = None
    obs = None
    if light_tag in ['zenith', '57.5']:
        dfobs = ground_cover_data(variety, tag, angle=tag_to_zenith(light_tag))
        obs = 'po'
    elif light_tag == 'spray':
        dfobs = petri_data(variety, tag, 'soil')
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


def compare_all_po(model='adel',variety='Tremie12', nplants=30, tag='reference', rep=1, top_petri=False,
                   start=None, stop=None, by=None, at=None, reset=False,
                   reset_build=False, reset_reconstruction=False):
    ax = compare_po(model=model, variety=variety, nplants=nplants, tag=tag, rep=rep,
                    light_tag='57.5', z_soil=0, start=start, stop=stop, by=by,
                    at=at, reset=reset, reset_build=reset_build,
                    reset_reconstruction=reset_reconstruction)
    compare_po(model=model, variety=variety, nplants=nplants, tag=tag, rep=rep,
               light_tag='zenith', z_soil=0, start=start, stop=stop, by=by, at=at,
               reset=reset, reset_build=reset_build,
               reset_reconstruction=reset_reconstruction, ax=ax, color='r')
    compare_po(model=model, variety=variety, nplants=nplants, tag=tag, rep=rep,
               light_tag='soc', z_soil=0, start=start, stop=stop, by=by, at=at,
               reset=reset, reset_build=reset_build,
               reset_reconstruction=reset_reconstruction, ax=ax, color='green')
    compare_po(model=model, variety=variety, nplants=nplants, tag=tag, rep=rep,
               light_tag='spray', z_soil=0, start=start, stop=stop, by=by, at=at,
               reset=reset, reset_build=reset_build,
               reset_reconstruction=reset_reconstruction, ax=ax, color='m')

    if top_petri:
        dfobs = petri_data(variety, tag, 'top')
        if dfobs is not None:
            plot_mean(dfobs, 'po', xaxis='HS', ax=ax, color='darkviolet', linestyle='')


    xlabel = 'Haun Stage'
    ylabel = 'Gap Fraction'
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)

    return ax


def check_lai57(variety='Tremie12', nplants=30, tag='reference', rep=1,
                start=None, stop=None, by=None, at=None, reset=False,
                reset_build=False, reset_reconstruction=False):
    df = simulate_lai(variety, nplants, tag, rep, start, stop, by, at, reset,
                      reset_build, reset_reconstruction)
    po = simulate_po_light('57.5')
    df = df.merge(po)
    fig, ax = plt.subplots()

    ax.plot(list(range(7)), list(range(7)), color='b', linestyle='--')
    ax.plot(df['LAI_tot'], df['LAI_vert'], color='g', linestyle='--')
    ax.plot(df['LAI_tot'], df['PAI_tot'], color='m', linestyle='--')
    ax.plot(df['LAI_tot'], df['lai_57.5_z0'], marker='d',markersize=7,color='m')

    ax.set_xlim((0, 6))
    ax.set_ylim((0, 6))
    xlabel = 'LAI'
    ylabel = 'effective LAI'
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)

    return ax


def check_lai2000(variety='Tremie12', nplants=30, tag='reference', rep=1,
                start=None, stop=None, by=None, at=None, reset=False,
                reset_build=False, reset_reconstruction=False):
    df = simulate_lai(variety, nplants, tag, rep, start, stop, by, at, reset,
                      reset_build, reset_reconstruction)
    po = simulate_po_light()
    df = df.merge(po)
    fig, ax = plt.subplots()

    ax.plot(list(range(7)), list(range(7)), color='b', linestyle='--')
    ax.plot(df['LAI_tot'], df['LAI_vert'], color='g', linestyle='--')
    ax.plot(df['LAI_tot'], df['PAI_tot'], color='m', linestyle='--')
    ax.plot(df['LAI_tot'], df['lai_lai2000'], marker='d',markersize=7,color='b')

    ax.set_xlim((0, 6))
    ax.set_ylim((0, 6))
    xlabel = 'LAI'
    ylabel = 'effective LAI'
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)

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



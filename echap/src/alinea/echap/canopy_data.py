"""Canopy level data"""

import pandas
import numpy
import scipy.stats

from alinea.echap.hs_tt import tt_hs_tag, derived_data_path
from alinea.echap.weather_data import par_transmittance

import alinea.echap
from openalea.deploy.shared_data import shared_data



def conf_int(lst, perc_conf=95):
    """
    Confidence interval - given a list of values compute the square root of
    the variance of the list (v) divided by the number of entries (n)
    multiplied by a constant factor of (c). This means that I can
    be confident of a result +/- this amount from the mean.
    The constant factor can be looked up from a table, for 95pcent confidence
    on a reasonable size sample (>=500) 1.96 is used.
    """

    n, v = len(lst), numpy.var(lst, ddof=1)
    c = scipy.stats.t.interval(perc_conf * 1.0 / 100, n - 1)[1]

    return numpy.sqrt(v / n) * c

#
# To do : include in compil nb plant/rep and generate from notation stem area
#
def scan_data():
    data_file = shared_data(alinea.echap, 'architectural_measurements/Compil_scan.csv')
    df = pandas.read_csv(data_file, decimal='.', sep=',')
    df['daydate'] = pandas.to_datetime(df['prelevement'].values, dayfirst=True).strftime('%Y-%m-%d')
    return df
#
# Biomass scanned plants / biomass microplot
#
def plot_biomass_data():
    t12 = {'daydate': ['2012-03-09'] * 3 + ['2012-04-11'] * 3 + ['2012-05-09'] * 3,
           'rep': [1, 2, 3] * 3,
           'nb_plant_scanned': [10] * 9,
           'dry_mass_scanned_blade': [2.3912, 1.4955, 1.8367] + [6.1498, 7.252,7.895]
                                       + [11.16599, 11.44891, 10.1178],
           'dry_mass_scanned_stem': [1.5412, 0.97, 1.1083] + [10.1299, 10.3606, 11.2501]
                                       + [30.60269, 34.3479, 29.9004],
           'nb_plant_plot': [123, 142, 141] + [129, 117, 135] + [110, 129, 113],
           'dry_mass_plot': [46.5, 51.8, 44.4] + [180.5, 185.0, 179.2] + [438.1, 484.4, 418.5],
           'plot_area': [5 * 0.143 * 0.6] * 9
           }
    df12 = pandas.DataFrame(t12)
    df12['variety'] = 'Tremie12'
    df12['dry_mass_scanned'] = df12['dry_mass_scanned_blade'] + df12['dry_mass_scanned_stem']
    t13 = {'daydate': ['2013-04-22'] * 3 + ['2013-05-03'] * 3,
           'rep': [1, 2, 3] * 2,
           'nb_plant_scanned': [15, 13, 15] + [15, 15, 15],
           'dry_mass_scanned': [14.95, 17.12, 15.89] + [None] * 3,
           'nb_plant_plot': [118, 116, 128] + [98, 130, 126],
           'nb_axe_plot': [625, 615, 678] + [None] * 3,
           'dry_mass_plot': [133.1, 136.8, 147.2] + [None] * 3, # dry mass not available for scanned plants [180.5, 185.0, 179.2],
           'plot_area': [5 * 0.15 * 0.6] * 6}
    df13 = pandas.DataFrame(t13)
    df13['variety'] = 'Tremie13'
    return pandas.concat([df12, df13], axis=0)


# Green pixels on  vertical and obliques images
def image_analysis_data():
    genos = ['MerciaRht3', 'Tremie12', 'Tremie13']
    paths = [shared_data(alinea.echap) / 'canopy_data' / var + '_vertical_images.csv' for var in genos]
    data = [pandas.read_csv(path, decimal=',', sep=';') for path in paths]
    df = pandas.concat(data)
    df['daydate'] = pandas.to_datetime(df['date'].values,
                                       dayfirst=True).strftime('%Y-%m-%d')
    df['zenith_angle'] = 0
    df['po'] = (100 - df['pcent_veg']) / 100.
    paths = [shared_data(alinea.echap) / 'canopy_data' / var + '_inclined_images.csv'
     for var in genos]
    data = [pandas.read_csv(path, decimal=',', sep=';') for path in paths]
    df57 = pandas.concat(data)
    df57['daydate'] = pandas.to_datetime(df57['date'].values,
                                         dayfirst=True).strftime('%Y-%m-%d')
    df57['zenith_angle'] = 57.5
    df57['po'] = (100 - df57['pcent_veg']) / 100.
    what = ['daydate', 'variety', 'zenith_angle', 'po']

    return pandas.concat((df[what], df57[what]), axis=0)

#derived / cached data

def ground_cover_data(variety='Mercia', tag='reference', angle=0):
    df = None
    path = derived_data_path(tag) / 'ground_cover_aggregated.csv'
    try:
        df = pandas.read_csv(path)
    except IOError:
        df = image_analysis_data()
        df['ground_cover'] = 1 - df['po']
        df = df.groupby(['variety', 'daydate', 'zenith_angle']).agg([numpy.mean, conf_int])
        df.columns = ['_'.join(c) for c in df.columns]
        df = df.reset_index()
        df.to_csv(path, index=False)

    df = df.loc[(df['variety'] == variety) & (df['zenith_angle'] == angle),:]
    tths = tt_hs_tag(variety, tag)

    return df.merge(tths)


def transmittance_data(variety='Mercia', tag='reference', reset=False):
    tths = tt_hs_tag(variety, tag)
    path = derived_data_path(tag) / 'transmittance_aggregated.csv'
    if not reset:
        try:
            df = pandas.read_csv(path)
            return df.merge(tths)
        except IOError:
            pass
    df = par_transmittance(variety)
    df = df.loc[df['daydate'].isin(tths['daydate']),:]
    df.to_csv(path, index=False)
    return df.merge(tths)


def pai57(variety='Mercia', tag='reference'):
    df = None
    path = derived_data_path(tag) /  'pai_aggregated.csv'
    try:
        df = pandas.read_csv(path)
    except IOError:
        df = image_analysis_data()
        g_miller = 0.5 / numpy.cos(numpy.radians(df['zenith_angle']))
        df['pai'] = - numpy.log(df['po']) / g_miller
        df = df.groupby(['variety', 'daydate', 'zenith_angle']).agg([numpy.mean, conf_int])
        df.columns = ['_'.join(c) for c in df.columns]
        df = df.reset_index()
        df.to_csv(path, index=False)

    df = df.loc[(df['variety'] == variety) & (df['zenith_angle'] == 57.5),:]
    tths = tt_hs_tag(variety, tag)

    return df.merge(tths)

def lai_pai_scan(variety='Tremie12', tag='reference', reset=False):
    """ Lai/Pai estimates from scanned leaves / biomass data"""
    df = None
    path = derived_data_path(tag) / 'lai_pai_scan.csv'
    try:
        if reset:
            raise IOError
        df = pandas.read_csv(path)
    except IOError:
        bm = plot_biomass_data()
        scan = scan_data().loc[:, (
        'variety', 'daydate', 'rep', 'N', 'id_Axe', 'A_bl', 'A_bl_green',
        'stem_half_area')]
        scan = scan.loc[scan['daydate'].isin(bm['daydate']),:]
        areas = scan.groupby(('daydate','rep')).agg(sum).reset_index()
        areas = areas.drop(['stem_half_area', 'N'],axis=1)
        stem_areas_ms = scan.groupby(('daydate', 'rep', 'N')).agg(
            numpy.mean).reset_index().groupby(('daydate', 'rep')).agg(
            sum).reset_index().loc[:,
                        ('daydate', 'rep', 'stem_half_area')].rename(
            columns={'stem_half_area': 'A_stem_MS'})
        areas_axe = scan.groupby(('daydate','rep', 'id_Axe')).agg(sum).reset_index()
        areas_ms = areas_axe.loc[areas_axe['id_Axe'] == 'MB', (
        'daydate', 'rep', 'A_bl', 'A_bl_green')].rename(
            columns={'A_bl_green': 'A_bl_green_MS', 'A_bl': 'A_bl_MS'})
        df = bm.merge(areas)
        df = df.merge(areas_ms)
        df = df.merge(stem_areas_ms)
        frac_scanned = df['nb_plant_scanned'] / df['nb_plant_plot']
        frac_scanned_drymass = df['dry_mass_scanned'] / df['dry_mass_plot']
        df['correction_factor'] = frac_scanned / frac_scanned_drymass
        df['stem_leaf_ratio_biomass'] = df['dry_mass_scanned_stem'] / df[
            'dry_mass_scanned_blade']
        df['stem_leaf_ratio_area'] = df['A_stem_MS'] / df[
            'A_bl_MS']
        df['GLAI_density'] = df['A_bl_green'] * 1e-4 / frac_scanned / df['plot_area']
        df['GLAI_biomass'] = df['A_bl_green'] * 1e-4 / frac_scanned_drymass / df['plot_area']
        df['GLAI_MS_density'] = df['A_bl_green_MS'] * 1e-4 / frac_scanned / df['plot_area']
        df['GLAI_MS_biomass'] = df['A_bl_green_MS'] * 1e-4 / frac_scanned_drymass / df['plot_area']
        # for LAI /A_bl, keep 09-03-2012 only, as only this sampling includes dried leaves
        df['LAI_density'] = numpy.where(df['daydate'] == '2012-03-09',
                                        df['A_bl'] * 1e-4 / frac_scanned / df['plot_area'], numpy.nan)
        df['LAI_biomass'] = numpy.where(df['daydate'] == '2012-03-09',
                                        df['A_bl'] * 1e-4 / frac_scanned_drymass / df['plot_area'],
                                        numpy.nan)
        df['plant_density'] = df['nb_plant_plot'] / df['plot_area']
        df['greeness'] = df['A_bl_green'] / df['A_bl']
        df.to_csv(path, index=False)
        #

    df = pandas.concat((df.loc[:,('variety', 'daydate')], df.ix[:, 'stem_leaf_ratio_area':]), axis=1)
    df = df.groupby(['variety', 'daydate']).agg([numpy.mean, conf_int])
    df.columns = ['_'.join(c) for c in df.columns]
    df = df.reset_index()
    df = df.loc[df['variety'] == variety, :]
    tths = tt_hs_tag(variety, tag)

    return df.merge(tths)


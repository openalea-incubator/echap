""" isolated Thermal/time / HS fit mmodule to be used by other modules"""
import pandas
import numpy
import os
import json
import scipy.stats as stats

from alinea.astk.Weather import linear_degree_days
from alinea.adel.plantgen_extensions import HaunStage
from alinea.echap.weather_data import get_weather
from openalea.deploy.shared_data import shared_data
import alinea.echap


def cache_reconstruction_path(tag):
    path = shared_data(alinea.echap,share_path='../../share') / 'cache' / 'reconstructions' / tag
    if not os.path.exists(str(path)):
        os.makedirs(str(path))
    return path


def derived_data_path(tag=None):
    if tag is None:
        path = shared_data(alinea.echap, share_path='../../share') / 'cache' / 'derived_data'
    else:
        path = shared_data(alinea.echap, share_path='../../share') / 'cache' / 'derived_data' / tag
    if not os.path.exists(str(path)):
        os.makedirs(str(path))
    return path


def hs_fit_parameters(tag='reference', reset=False):
    path = cache_reconstruction_path(tag) / 'hs_fit_parameters.json'

    if not reset:
        try:
            with open(path, 'r') as input_file:
                cached = json.load(input_file)
            return cached
        except IOError:
            pass
    pars = {}
    #
    # Haun Stage = f(TT), convergence between axis
    # ---------------------------------------------
    # delay between emergence of flag leaf on mainstem and flag leaf emergence
    # on cohorts (60 is for Maxwell that has quite large desynchronisation,
    # 30 may be more realistic)
    pars['dTT_cohort'] = {k: {'first': 30, 'increment': 10} for k in
                          ('Mercia', 'Rht3', 'Tremie12', 'Tremie13')}
    with open(path, 'w') as output_file:
        json.dump(pars, output_file, sort_keys=True, indent=4,
                  separators=(',', ': '))

    return pars


def tt_lin(variety='Mercia', reset=False):
    """ compute the daydate <-> TT mapping from sowing to harvest"""

    sowing = {'Mercia': '2010-10-05', 'Rht3': '2010-10-05',
              'Tremie12': '2011-10-21', 'Tremie13': '2012-10-29'}
    harvest = {'Mercia': '2011-06-20', 'Rht3': '2011-06-20',
               'Tremie12': '2012-06-19', 'Tremie13': '2013-06-09'}
    filename = {'Mercia': 'MerciaRht3_TTlin_sowing_harvest.csv',
                'Rht3': 'MerciaRht3_TTlin_sowing_harvest.csv',
                'Tremie12': 'Tremie12_TTlin_sowing_harvest.csv',
                'Tremie13': 'Tremie13_TTlin_sowing_harvest.csv'}

    path = derived_data_path() / filename[variety]
    if not reset:
        try:
            return pandas.read_csv(path)
        except IOError:
            pass
    w = get_weather(variety)
    tt = linear_degree_days(w.data, start_date=sowing[variety])
    seq = pandas.date_range(sowing[variety], harvest[variety], tz='UTC')
    df = pandas.DataFrame(
        {'daydate': seq.strftime('%Y-%m-%d'), 'TT': tt.reindex(seq).values})
    df.to_csv(path, index=False)

    return df


def TT_lin():
    def _read(label):
        df = tt_lin(label)
        df['label'] = label
        return df
    data = [_read(k) for k in ('Mercia', 'Rht3', 'Tremie12', 'Tremie13')]
    return pandas.concat(data)


def Pheno_data(pheno_dict={},
               sources=['archi_tagged', 'archi_sampled', 'symptom_tagged'],
               count=0):
    """ Phenological data (HS, GL, SSI) found for
        - Treated tagged plants followed for architectural notations
        - Treated tagged plant followed for symptom notation (SSI/GL)
        - some destructive samples (scan /silhouettes data

        Details in architectural_measurements R pre-processing scripts
    """

    def dateparse(x):
        return pandas.datetime.strptime(x, '%d/%m/%Y')
    src = sources[count]
    count += 1
    filename = 'Compil_Pheno_treated_' + src + '.csv'
    filepath = str(
        shared_data(alinea.echap,share_path='../../share') / 'architectural_measurements' / filename)
    df = pandas.read_csv(filepath, sep=',', decimal='.')
    df['Date'] = df['Date'].apply(dateparse)
    df['daydate'] = df.set_index('Date').index.strftime('%Y-%m-%d')
    df_TT = TT_lin()
    df = df.merge(df_TT)
    pheno_dict[src] = df
    if count < len(sources):
        return Pheno_data(pheno_dict, sources, count)
    else:
        return pheno_dict


def application_tag(variety, daydate, which='T1'):
    """ generate tag relative to pesticide application date

    :param variety: variety name
    :param daydate: a list of daydate strings
    :return: tags
    """

    sowing = {'Mercia': '2010-10-05', 'Rht3': '2010-10-05',
              'Tremie12': '2011-10-21', 'Tremie13': '2012-10-29'}
    t1 = {'Mercia': '2011-04-19', 'Rht3': '2011-04-19',
          'Tremie12': '2012-04-11', 'Tremie13': '2013-04-25'}
    t2 = {'Mercia': '2011-05-11', 'Rht3': '2011-05-11',
          'Tremie12': '2012-05-09', 'Tremie13': '2013-05-17'}

    if which == 'T1':
        origin = t1[variety]
    elif which == 'T2':
        origin = t2[variety]
    else:
        origin = sowing[variety]

    delta = (pandas.to_datetime(daydate) - pandas.to_datetime(origin)).days
    tags = numpy.where(delta == 0, which, numpy.where(delta > 0, [which + '+' + str(x) for x in delta], [which + '-' + str(abs(x)) for x in delta]))
    return tags


def linreg_df(x, y):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return pandas.Series(
        {'slope': slope, 'intercept': intercept, 'r2': r_value ** 2,
         'std_err': std_err})


def fit_hs(tag='reference'):
    """ Linear regression for phyllochron data
    """
    parameters = hs_fit_parameters(tag)
    data = Pheno_data()
    tagged = data['archi_tagged']
    sampled = data['archi_sampled']
    dTT_cohort = parameters.get('dTT_cohort')

    # HS of mean plant
    g = tagged.groupby('label')

    def _fit(df):
        dat = df.loc[df['HS'] < df['nff'], ('TT', 'HS')].dropna()
        nff = df.groupby('N').agg('mean')['nff'].mean()
        n = df.groupby('N').agg('mean')['nff'].count()
        res = linreg_df(dat['TT'], dat['HS'])
        res['nff'] = nff
        res['n'] = n
        return res

    hs_ms = g.apply(_fit)
    # complete fit if destructive sample available
    # for Tremie 12, keep only first sample
    # for Tremie 13, sampled plants are late / tagged plant : not used
    gd = sampled.groupby('label')
    for lab in gd.groups:
        #        if lab == 'Tremie12':
        dat = g.get_group(lab)
        dat = dat.loc[dat['HS'] < dat['nff'], ('TT', 'HS')].dropna()
        datd = gd.get_group(lab)
        datd = datd.loc[
            (datd['HS'] < datd['nff']) | (numpy.isnan(datd['nff'])), (
                'TT', 'HS')].dropna()
        # datd = datd.loc[datd['HS'] < 10,:]
        datreg = pandas.concat((dat, datd))
        res = linreg_df(datreg['TT'], datreg['HS'])
        hs_ms.loc[lab, 'intercept':'std_err'] = res
    # TTem per label
    TTem = - hs_ms['intercept'] / hs_ms['slope']
    # Flag leaf delay per nff
    TT_mean_flag = (hs_ms['nff'] - hs_ms['intercept']) / hs_ms['slope']
    g = tagged.groupby(['label', 'nff'])
    hs_nff = g.apply(_fit)
    TT_nff_flag = (hs_nff['nff'] - hs_nff['intercept']) / hs_nff['slope']
    TT_flag = TT_nff_flag.reset_index().rename(columns={0: 'flag'}).merge(
        TT_mean_flag.reset_index().rename(columns={0: 'mean_flag'}))
    TT_flag['dTT_flag'] = TT_flag['flag'] - TT_flag['mean_flag']
    TT_flag = TT_flag.merge(
        hs_nff['n'].reset_index().rename(columns=({0: 'n'})))
    TT_flag = TT_flag.merge(
        hs_ms['nff'].reset_index().rename(columns=({'nff': 'mean_nff'})))
    TT_flag = TT_flag.merge(
        TTem.reset_index().rename(columns=({0: 'mean_TTem'})))
    dTTnff = TT_flag.groupby('label').apply(
        lambda x: numpy.average(x['dTT_flag'] * 1. / (x['nff'] - x['mean_nff']),
                                weights=x['n']))
    # residual variabilityy of emergenece between plants
    # mean slopes per nff:
    dat = TT_flag.set_index(['label', 'nff'], drop=False)
    a_nff = dat['nff'] / (dat['flag'] - dat['mean_TTem'])
    # estimate dTTEm per plant wyth modeled slope
    g = tagged.groupby(['label', 'N'])

    def _TTem(df):
        res = None
        dat = df.loc[df['HS'] < df['nff'], ('TT', 'HS')].dropna()
        nff = df['nff'].values[0]
        label = df['label'].values[0]
        if not numpy.isnan(nff):
            a = a_nff[(label, nff)]
            b = numpy.mean(dat['HS'] - a * dat['TT'])
            res = -1. * b / a
        return res

    TTem_p = g.apply(_TTem)
    # standard deviation
    std_TTem = TTem_p.reset_index().groupby('label').agg('std').loc[:, 0]

    # HS fits
    hs_fits = {
        k: HaunStage(1. / hs_ms['slope'][k], TTem[k], std_TTem[k],
                              hs_ms['nff'][k], dTTnff[k] / hs_ms['slope'][k],
                              dTT_cohort[k]) for k in dTT_cohort}
    return hs_fits


def HS_fit(tag='reference', reset=False):
    """ Handle fitted phyllochron persitent object
    """
    dir_path = cache_reconstruction_path(tag)
    if not reset:
        try:
            hs_fits = {}
            for k in ('Mercia', 'Rht3', 'Tremie12', 'Tremie13'):
                file_path = dir_path / 'HSfit_' + k + '.json'
                hs_fits[k] = HaunStage.load(file_path)
            return hs_fits
        except IOError:
            pass
    hs_fits = fit_hs(tag)
    for k, hs in hs_fits.items():
        file_path = dir_path / 'HSfit_' + k + '.json'
        hs.dump(file_path)

    return hs_fits


def tt_hs_tag(variety='Mercia', tag='reference'):
    """A multiple time index"""

    df = None
    path = derived_data_path(tag) / variety + '_TT_HS_tag.csv'
    try:
        df = pandas.read_csv(path)
    except IOError:
        hsfit = HS_fit(tag)
        df = tt_lin(variety)
        df['HS'] = hsfit[variety](df['TT'])
        for t in ('T1', 'T2'):
            df['tag_' + t] = application_tag(variety, df['daydate'].values, t)
        df.to_csv(path, index=False)

    return df


def as_daydate(daydate, tths):
    if str(daydate).startswith('T1'):
        return tths.set_index('tag_T1')['daydate'][daydate]
    elif str(daydate).startswith('T2'):
        return tths.set_index('tag_T2')['daydate'][daydate]
    else:
        return daydate


def daydate_range(variety, tag, start=None, stop=None, by=None, at=None):
    tths = tt_hs_tag(variety, tag)
    if at is None:
        if start is None:
            start = tths['daydate'][0]
        else:
            start = as_daydate(start, tths)
        if stop is None:
            stop = tths['daydate'][-1]
        else:
            stop = as_daydate(stop, tths)
        at = tths.set_index('daydate').ix[start:stop:by,].index.values.tolist()
    else:
        at = [as_daydate(x, tths) for x in at]

    return at
import numpy
import pandas
import scipy.stats
import matplotlib.pyplot as plt
from openalea.deploy.shared_data import shared_data
import alinea.echap

#from alinea.echap.evaluation_canopy import get_lai_sim, get_light_sim


def get_file_name(variety = 'Tremie12', nplants = 30, group = None):
    if group is None:
        return 'leaf_'+variety.lower() + '_' + str(nplants) + 'pl.csv'
    else:
        assert group in ['sum', 'mean'], ValueError("group unknown, try 'sum' or 'mean'")
        return 'leaf_'+variety.lower() + '_' + str(nplants) + 'pl_' + group + '.csv'

def get_file_path(variety = 'Tremie12', nplants = 30, group = None):
    filename = get_file_name(variety = variety, nplants = nplants, group = group)
    return str(shared_data(alinea.echap)/'architectural_simulations'/filename)


def get_simu_results(variety = 'Tremie12', nplants = 30, group = None):
    file_path = get_file_path(variety = variety, nplants = nplants, group = group)
    df = pandas.read_csv(file_path)
    df['necro'] = 100*df['senesced_area']/df['area'].replace({ 0 : numpy.inf })
    df['necro_tot'] = 100*df['senesced_area']/df['area'].replace({ 0 : numpy.inf })
    return df

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
    ci = numpy.sqrt(v / n) * c
    if numpy.isnan(ci):
        ci = 0
    return ci


def plot_mean(data, variable='LAI_vert', xaxis='ThermalTime', marker='d', markersize=9,
              empty_marker=False, linestyle='-', color='b', title=None,
              xlabel=None, ylabel=None, xlims=None, ylims=None, ax=None, logy=False):
    var_mean = '_'.join((variable, 'mean'))
    var_err = '_'.join((variable, 'conf_int'))
    df = data.copy(deep=True)
    if var_mean not in df.columns:
        if variable not in df.columns:
            raise KeyError('variable ' + variable + ' not found in data!')
        else:
            if logy:
                df[variable] = numpy.log(df[variable])
            df = df.groupby(xaxis).agg([numpy.mean, conf_int])
            df.columns = ['_'.join(c) for c in df.columns]
            df = df.reset_index()

    if ax is None:
        fig, ax = plt.subplots()
    if empty_marker:
        markerfacecolor = 'none'
    else:
        markerfacecolor = color

    df = df.loc[:, [xaxis, var_mean, var_err]].dropna()
    if df[var_err].sum() > 0:
        ax.errorbar(df[xaxis], df[var_mean].values, yerr=df[var_err].values,
                    marker=marker, linestyle=linestyle, color=color,
                    markerfacecolor=markerfacecolor, markeredgecolor=color, markersize=markersize)
    else:
        ax.plot(df[xaxis], df[var_mean].values,
                    marker=marker, linestyle=linestyle, color=color,
                    markerfacecolor=markerfacecolor, markeredgecolor=color, markersize=markersize)
    if title is not None:
        ax.set_title(title, fontsize=18)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=18)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=18)
    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)

    return ax




# deprecated
def plot_sum(data, variable='LAI_vert', xaxis='ThermalTime', marker='d',
             linestyle='-', color='b', title=None, xlabel=None, ylabel=None,
             xlims=None, ylims=None, ax=None, return_ax=False):
    if variable in data.columns:
        if ax == None:
            fig, ax = plt.subplots()
        df = data[pandas.notnull(data.loc[:, variable])]

        df_sum = df.groupby([xaxis, 'num_plant']).sum()
        df_sum = df_sum.reset_index()
        df_sum = df_sum.groupby(xaxis).mean()

        ax.plot(df_sum.index, df_sum[variable], color=color, marker=marker,
                linestyle=linestyle)
        if title is not None:
            ax.set_title(title, fontsize=18)
        if xlabel is not None:
            ax.set_xlabel(xlabel, fontsize=18)
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=18)
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)
        if return_ax == True:
            return ax


# Get observation data for canopy properties #######################################################





def get_all_obs(variety='Tremie12', origin_lai_data='biomass'):
    df_lai = get_lai_obs(variety=variety, origin=origin_lai_data)
    df_TC = get_TC_obs(variety=variety)
    df_rad = get_radiation_obs(variety=variety)
    df_all = pandas.concat([df_lai, df_TC, df_rad])
    return df_all.reset_index(drop=True)


def plot_sim_obs(df_sim, df_obs=None, variable='LAI_vert', xaxis='HS',
                 title=None, xlabel='Haun Stage', ylabel='LAI',
                 colors=['b', 'r'], markers=['d', 'o'], linestyles=['-', '--'],
                 error_bars=[True, True], xlims=None, ylims=None, legend=True,
                 ax=None):
    if ax == None:
        fig, ax = plt.subplots()
    plot_mean(df_sim, variable=variable, xaxis=xaxis, color=colors[0],
              marker=markers[0], linestyle=linestyles[0], title=title,
              xlabel=xlabel, ylabel=ylabel, xlims=xlims, ylims=ylims, ax=ax)
    if df_obs is not None:
        plot_mean(df_obs, variable=variable, xaxis=xaxis, color=colors[1],
                  marker=markers[1], linestyle=linestyles[1], title=title,
                  xlabel=xlabel, ylabel=ylabel, xlims=xlims, ylims=ylims, ax=ax)

    if legend == True:
        ax.legend(['Simulated', 'Observed'], loc='center left',
                  bbox_to_anchor=(1, 0.5))





def test_architecture_canopy_single(variety='Tremie12', nplants=30, nrep=1,
                                    fig_size=(10, 10), color='b', title=None,
                                    axs=None, tight_layout=False):
    if axs.all() == None:
        fig, axs = plt.subplots(2, 2, figsize=fig_size)
    df_sim = get_simu_results(variety=variety, nplants=nplants)
    df_obs = get_all_obs(variety=variety, origin_lai_data='biomass')
    plot_sim_obs(df_sim, df_obs, variable='TCgreen_57', xaxis='HS',
                 xlabel='Haun stage', ylabel='Cover fraction (oblique view)',
                 colors=[color, color], markers=['', 'o'],
                 linestyles=['-', '--'], error_bars=[False, True], legend=False,
                 ax=axs[0][0])
    plot_sim_obs(df_sim, df_obs, variable='TCgreen', xaxis='HS',
                 xlabel='Haun stage', ylabel='Cover fraction (vertical view)',
                 colors=[color, color], markers=['', 'o'],
                 linestyles=['-', '--'], error_bars=[False, True], legend=True,
                 ax=axs[0][1])
    plot_sim_obs(df_sim, df_obs, variable='LAI_vert', xaxis='HS',
                 xlabel='Haun stage', ylabel='Green Leaf Area Index',
                 colors=[color, color], markers=['', 'o'],
                 linestyles=['-', '--'], error_bars=[False, True], legend=False,
                 ax=axs[1][0])
    # Temp
    if not 'LightInterception_0' in df_sim.columns:
        df_sim['LightInterception_0'] = 1 - df_sim['LightPenetration_0']
    plot_sim_obs(df_sim, df_obs, variable='LightInterception_0', xaxis='HS',
                 xlabel='Haun stage',
                 ylabel='Intercepted fraction of radiation',
                 colors=[color, color], markers=['', ''],
                 linestyles=['-', '--'], error_bars=[False, False],
                 legend=False, ax=axs[1][1])
    letters = iter(['a', 'b', 'c', 'd'])
    for ax in axs.flat:
        ax.annotate(next(letters), xy=(0.05, 0.85), xycoords='axes fraction',
                    fontsize=18)
    if title is not None:
        ax.set_title(title, fontsize=18)
    if tight_layout == True:
        plt.tight_layout()


def test_architecture_canopy(
        varieties=['Mercia', 'Rht3', 'Tremie12', 'Tremie13'], nplants=30,
        nrep=1, fig_size=(15, 15), color='b', title=None, tight_layout=False):
    fig, axs = plt.subplots(2, 2, figsize=fig_size)
    colors = {'Mercia': 'r', 'Rht3': 'g', 'Tremie12': 'b', 'Tremie13': 'm'}
    proxy = []
    for variety in varieties:
        color = colors[variety]
        test_architecture_canopy_single(variety=variety, nplants=nplants,
                                        nrep=nrep, fig_size=fig_size,
                                        color=color, title=title, axs=axs,
                                        tight_layout=tight_layout)
        proxy += [
            plt.Line2D((0, 1), (0, 0), color=color, marker='', linestyle='-')]
    axs[1][1].legend(proxy, varieties, loc='center left',
                     bbox_to_anchor=(1, 0.5))


def compare_lai_by_axis(nplants=30, nrep=1, fig_size=(10, 10), color='b',
                        title=None, tight_layout=False, only_lai=False):
    fig, axs = plt.subplots(2, 2, figsize=fig_size)
    linestyles = {'total': '-', 'MS': '--', 'ferti': '-.'}
    varieties = ['Mercia', 'Rht3', 'Tremie12', 'Tremie13']
    colors = {'Mercia': 'r', 'Rht3': 'g', 'Tremie12': 'b', 'Tremie13': 'm'}
    axs_iter = iter(axs.flat)
    proxys = []
    labels = []
    for variety in varieties:
        df_sim = get_simu_results(variety=variety, nplants=nplants, nrep=nrep,
                                  only_lai=only_lai)
        df_obs = get_lai_obs(variety=variety, origin='biomass')
        ax = next(axs_iter)
        color = colors[variety]
        proxys += [
            plt.Line2D((0, 1), (0, 0), color=color, marker='', linestyle='-')]
        labels += [variety]
        plot_sim_obs(df_sim, df_obs, variable='LAI_vert', xaxis='HS',
                     xlabel='Haun stage', ylabel='Green Leaf Area Index',
                     colors=[color, color], markers=['', 'o'],
                     linestyles=[linestyles['total'], ''],
                     error_bars=[False, True], legend=False, ylims=[0, 10],
                     ax=ax)
        plot_sim_obs(df_sim, df_obs, variable='LAI_vert_MS', xaxis='HS',
                     xlabel='Haun stage', ylabel='Green Leaf Area Index',
                     colors=[color, color], markers=['', 'o'],
                     linestyles=[linestyles['MS'], ''],
                     error_bars=[False, True], legend=False, ylims=[0, 10],
                     ax=ax)
        plot_sim_obs(df_sim, df_obs, variable='LAI_vert_ferti', xaxis='HS',
                     xlabel='Haun stage', ylabel='Green Leaf Area Index',
                     colors=[color, color], markers=['', 'o'],
                     linestyles=[linestyles['ferti'], ''],
                     error_bars=[False, True], legend=False, ylims=[0, 10],
                     ax=ax)
    proxys += [plt.Line2D((0, 1), (0, 0), color='k', marker='', linestyle='-'),
               plt.Line2D((0, 1), (0, 0), color='k', marker='', linestyle='--'),
               plt.Line2D((0, 1), (0, 0), color='k', marker='', linestyle='-.')]
    labels += ['Total', 'MS', 'Ferti']
    axs[0][1].legend(proxys, labels, loc='center left', bbox_to_anchor=(1, 0.5))

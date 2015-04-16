# -*- coding: latin1 -*- 

""" Analysis and plotting of pandas DataFrames complying with the format generated 
    by readers in 'septo_data_reader'. """
    
from alinea.echap.disease.septo_data_reader import *

# Imports for statistical data description
import scipy
import scipy.stats as st
import scikits.bootstrap as bootstrap
from math import ceil
from scipy.stats.mstats import normaltest

# Imports for plotting
from datetime import datetime, timedelta, date
from matplotlib.dates import MonthLocator, DateFormatter
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from matplotlib._png import read_png
from matplotlib.offsetbox import TextArea, AnnotationBbox, OffsetImage
from scipy.interpolate import InterpolatedUnivariateSpline 

# Description of data ##############################################################################
def get_basic_stat(df, column='severity', by_leaf=True):
    """ Calculate basic statistical indicators on given column in data 
        (mean, std deviation, distribution) """
    if by_leaf==True:
        df_ = df.copy()
        change_index(df_, new_index = ['datetime', 'num_leaf_top'])
        df_ = df_.groupby(level=[0,1]).describe()
        return df_[column].unstack('num_leaf_top')
    else:
        df = df.groupby(level=0).describe()
        return df[column].unstack()

def get_std(df, column = 'severity'):
    """ Calculate standard deviation on given column in data """
    stats = get_basic_stat(data, column='severity', by_leaf=True)
    return stats.iloc[stats.index.get_level_values(level = 1) == 'std']
    
def get_mean_one_leaf(df, variable = 'severity', xaxis = 'degree_days',
                          num_leaf = 1, from_top = True):
    """ Get average of argument variable on argument leaf over all plants in canopy """
    if from_top == True:
        df = df[df['num_leaf_top'] == num_leaf]
    else:
        df = df[df['num_leaf_bottom'] == num_leaf]
    if xaxis in ['datetime', 'degree_days']:
        return df.groupby(xaxis).mean()[variable]
    elif xaxis in ['age_leaf', 'age_leaf_lig', 'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg']:
        df_mean = df.groupby('degree_days').mean()[variable]
        df_dates = get_df_dates_xaxis(df, xaxis)
        df_mean.index -= df_dates.loc[num_leaf, xaxis].astype(int)
        return df_mean
        
def get_mean(df, column='severity', xaxis = 'datetime', by_leaf=True, from_top = True):
    """ Get average of argument variable on all leaves over all plants measured """
    if by_leaf==True:
        df_ = df.copy()
        if from_top == True:
            num_leaf = 'num_leaf_top'
        else:
            num_leaf = 'num_leaf_bottom'
        if xaxis in ['datetime', 'degree_days']:
            change_index(df_, new_index = [xaxis, num_leaf])
            df_ = df_.groupby(level=[0,1]).mean()
            return df_[column].unstack()
        elif xaxis in ['age_leaf', 'age_leaf_lig', 'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg']:
            dfs = []
            for lf in set(df_[num_leaf]):
                df_mean_lf = get_mean_one_leaf(df_, variable = column, xaxis = xaxis, num_leaf = lf, from_top = from_top)
                df_mean_lf = df_mean_lf.to_frame()
                df_mean_lf.columns = [lf]
                dfs.append(df_mean_lf)
            return pd.concat(dfs, axis = 1)
    else:
        df = df.groupby(level=0).mean()
        return pd.DataFrame(df[column], columns=[column])
    
def get_standard_deviation(df, column='severity', by_leaf=True):
    """ Calculate lower and upper bounds of standard deviation around mean of column in data """
    stats = get_basic_stat(df, column=column, by_leaf=by_leaf)
    df_low = stats['mean']-stats['std']
    df_low[df_low<0]=0.
    df_high = stats['mean']+stats['std']
    df_high[df_high>100]=100
    return df_low, df_high
    
def get_error_margin_one_leaf(df, column='severity', xaxis = 'datetime', 
                                num_leaf = 1, by_leaf=True, from_top = True):
    """ Get 95% confidence interval of argument variable 
            on argument leaf over all plants in canopy """
    df_ = df.copy()
    if from_top == True:
        df_ = df_[df_['num_leaf_top'] == num_leaf]
    else:
        df_ = df_[df_['num_leaf_bottom'] == num_leaf]
    if xaxis in ['datetime', 'degree_days']:
        change_index(df_, new_index = ['degree_days'])
        df_ = df_[column]
        return df_.groupby(level=0).apply(lambda x: 2*st.sem(x[~np.isnan(x)]))
    elif xaxis in ['age_leaf', 'age_leaf_lig', 'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg']:
        change_index(df_, new_index = ['degree_days'])
        df_ = df_[column]
        df_ = df_.groupby(level=0).apply(lambda x: 2*st.sem(x[~np.isnan(x)]))
        df_dates = get_df_dates_xaxis(df, xaxis)
        df_.index -= df_dates.loc[num_leaf, xaxis].astype(int)
        return df_
    
def get_error_margin(df, column='severity', xaxis = 'datetime', by_leaf=True, from_top = True):
    """ Get 95% confidence interval of argument variable on all leaves over all plants """
    df_ = df.copy()
    if from_top == True:
        num_leaf = 'num_leaf_top'
    else:
        num_leaf = 'num_leaf_bottom'
    if by_leaf==True:
        if xaxis in ['datetime', 'degree_days']:
            change_index(df_, new_index = [xaxis, num_leaf])
            df_ = df_[column]
            df_ = df_.groupby(level=[0,1]).apply(lambda x: 2*st.sem(x[~np.isnan(x)]))
            return df_.unstack()
        elif xaxis in ['age_leaf', 'age_leaf_lig', 'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg']:
            dfs = []
            for lf in set(df_[num_leaf]):
                df_lf = get_error_margin_one_leaf(df_, column = column, xaxis = xaxis, num_leaf = lf, from_top = from_top)
                df_lf = df_lf.to_frame()
                df_lf.columns = [lf]
                dfs.append(df_lf)
            return pd.concat(dfs, axis = 1)
    else:
        df_ = df_.groupby(level=[0]).apply(lambda x: bootstr(x[~np.isnan(x)]))
        return df_
        
def bootstr(x):
    """ Calculate bootstrap sample """
    if len(x) == 0:
        return np.array([np.nan, np.nan])
    elif all(x.values[0] == item for item in x):
        return np.array([x.values[0], x.values[0]])
    else:
        try:
            return bootstrap.ci(data=x, statfunction=scipy.mean)
        except:
            if np.mean(x) == np.max(x):
                return(np.array([np.mean(x), np.mean(x)]))
            else:
                return np.array([np.nan, np.nan])

def get_bootstrap_error_margin_one_leaf(df, column='severity', xaxis = 'datetime', 
                                            num_leaf = 1, by_leaf=True, from_top = True):
    """ Get bootstraped 95% confidence interval of argument variable 
            on argument leaf over all plants in canopy """
    df_ = df.copy()
    if from_top == True:
        df_ = df_[df_['num_leaf_top'] == num_leaf]
    else:
        df_ = df_[df_['num_leaf_bottom'] == num_leaf]
    if xaxis in ['datetime', 'degree_days']:
        change_index(df_, new_index = ['degree_days'])
        df_ = df_[column]
        return df_.groupby(level=0).apply(lambda x: bootstr(x[~np.isnan(x)]))
    elif xaxis in ['age_leaf', 'age_leaf_lig', 'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg']:
        change_index(df_, new_index = ['degree_days'])
        df_ = df_[column]
        df_ = df_.groupby(level=0).apply(lambda x: bootstr(x[~np.isnan(x)]))
        df_dates = get_df_dates_xaxis(df, xaxis)
        df_.index -= df_dates.loc[num_leaf, xaxis].astype(int)
        return df_

def get_bootstrap_error_margin(df, column='severity', xaxis = 'datetime', by_leaf=True, from_top = True):
    """ Get 95% confidence interval of argument variable on all leaves over all plants """
    df_ = df.copy()
    if from_top == True:
        num_leaf = 'num_leaf_top'
    else:
        num_leaf = 'num_leaf_bottom'
    if by_leaf==True:
        if xaxis in ['datetime', 'degree_days']:
            change_index(df_, new_index = [xaxis, num_leaf])
            df_ = df_[column]
            df_ = df_.groupby(level=[0,1]).apply(lambda x: bootstr(x[~np.isnan(x)]))
            return df_.unstack()
        elif xaxis in ['age_leaf', 'age_leaf_lig', 'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg']:
            dfs = []
            for lf in set(df_[num_leaf]):
                df_lf = get_bootstrap_error_margin_one_leaf(df_, column = column, xaxis = xaxis, num_leaf = lf, from_top = from_top)
                df_lf = df_lf.to_frame()
                df_lf.columns = [lf]
                dfs.append(df_lf)
            return pd.concat(dfs, axis = 1)
    else:
        df_ = df_.groupby(level=[0]).apply(lambda x: bootstr(x[~np.isnan(x)]))
        return df_
        
def get_confidence_interval(df, weather, column='severity', by_leaf=True, 
                                xaxis = 'datetime', add_ddays = True):
    """ Calculate lower and upper bounds of 95% confidence interval
        around mean of column in data versus given xaxis """
    
    def separate_low_high(x, position=0):
        return x[position] if ~np.isnan(np.sum(x)) else np.nan
        
    df_mean = get_mean(df, column=column, by_leaf=by_leaf, xaxis = xaxis)
    df_count = table_count_notations(df, weather, variable = column, xaxis = xaxis, add_ddays = add_ddays)
    if any(x<20 for x in df_count.values.flat if x!='-' and not np.isnan(x)):
        df_low_high = get_bootstrap_error_margin(df, column=column, xaxis = xaxis, by_leaf = by_leaf)
        if not by_leaf:
            df_low_high = pd.DataFrame(df_low_high, xaxis = 'datetime')
        return df_low_high.applymap(lambda x: separate_low_high(x, 0)), df_low_high.applymap(lambda x: separate_low_high(x, 1))
    else:
        df_err = get_error_margin(df, column=column, by_leaf=by_leaf, xaxis = xaxis)
        df_low = df_mean-df_err
        df_low[df_low<0] = 0.
        df_high = df_mean+df_err
        df_high[df_high>100.] = 100.
    return df_low, df_high

def test_normality(df):
    df_normal_test = pd.DataFrame(index = ['p-values', 'normality'])
    for lf in df.columns[:4]:
        z, pval = normaltest(df[lf].values[np.where(df[lf].values==df[lf].values)])
        if pval<0.05:
            df_normal_test.loc[:,lf] = (pval, 'not normal')
        else:
            df_normal_test.loc[:,lf] = (pval, 'normal')
    return df_normal_test
    
# Save tables to latex #############################################################################
class TaggedFloat(float):
    """ Add a tag to a float """
    def __new__(self, value, tag):
        return float.__new__(self, value)
    
    def __init__(self, value, tag):
        self.tag = tag

def tagged_age_leaf_notations(data, weather, adel, nff=None):
    """ Tag particular notations to customize them in LaTeX """
    
    def tag_empty(x):
        if x != '-':
            return [TaggedFloat(x, '')]
    
    df_age = age_leaf_notations(data, weather, adel, nff=nff)
    df_sen = sen_leaf_notations(data, weather, adel, nff=nff)
    df_tag = df_age.applymap(tag_empty)
    for i_row, row in df_tag.iterrows():       
        for j_col in df_tag.ix[i_row].index:
            if df_age.ix[i_row, j_col]!='-' and df_age.ix[i_row, j_col] < 350:
                df_tag.ix[i_row, j_col] = [TaggedFloat(df_tag.ix[i_row, j_col][0], 'too_young')] # include in list to avoid reformating by panda
            if df_sen.ix[i_row, j_col]!='-' and df_sen.ix[i_row, j_col] > 0:
                df_tag.ix[i_row, j_col] = [TaggedFloat(df_tag.ix[i_row, j_col][0], 'senescent')]
    df_tag[df_tag==[None]] = '-'
    return df_tag

def format_latex_bold(x):
    return '\\textbf{'+str(int(x))+'}'

def format_latex_it(x):
    return '\\textit{'+str(int(x))+'}'

def bold_if_tagged(x):
    """ Save in bold or italic LaTeX text if value is a TaggedFloat """
    x = x[0]
    if x.tag == 'too_young':
        return format_latex_it(x)
    elif x.tag == 'senescent':
        return format_latex_bold(x)
    else:
        return str(int(x))
        
def save_table_to_latex(data, filename, format_model=None):
    with open(filename, "w") as f:
        txt = ('\documentclass[a4paper,reqno,11pt]{amsart}' + '\n' + '\n' +
                '\usepackage[french]{babel}' + '\n' +
                '\usepackage{booktabs}' + '\n' +
                '\\title{number of observations}' + '\n'+ '\n' + 
                '\\begin{document}' + '\n'+ '\n' + 
                '\\begin{tabular}{'+
                ''.join('l'*len(data.index.names))+
                ''.join('r'*len(data.columns))+ '}\n' +
                '\\toprule' + '\n' + '{} & &' + '& '.join(str(col) for col in data.columns) + '\\\\' + '\n' +
                data.index.names[0] + ' & ' + data.index.names[1] + ''.join(' & '*len(data.columns))+ '\\\\' + '\n' +
                '\midrule' + '\n')
        for i, row in data.iterrows():
            if format_model == None:
                txt+=str(i[0])[:10] + ' & ' + str(i[1]) + ' & ' + ' & '.join([str(int(x) if x!='-' else x) for x in row.values])+ '\\\\' + '\n'
            else:
                txt+=str(i[0])[:10] + ' & ' + str(i[1]) + ' & ' + ' & '.join([format_model(x) if x!='-' else x for x in row.values])+ '\\\\' + '\n'
            
        txt+='\\bottomrule' + '\n' + '\end{tabular}' + '\n' + '\n' +  '\end{document}'
        f.write(txt)

# Plotting functions for individual samples of 1 leaf ##############################################
def form_int(x, pos):
    return int(x)
    
def change_zero_sample(data, df, variable = 'severity', xaxis = 'degree_days'):
    """ Interpolate data of given variable between notations 
        to find the last date when variable = 0 """
    indx = np.nonzero(df[variable][~np.isnan(df[variable])])[0]
    if sum(df[variable]==0)>0 and len(indx)>=2:
        first_date_zero = df[xaxis][df[variable]==0].iloc[0]
        df_insert = pd.DataFrame(index = [df.index[0]], columns = df.columns)
        df_insert[variable] = 0.
        x_data = df.iloc[indx[:2]][variable].values
        y_data_xaxis = df.iloc[indx[:2]][xaxis].values
        s_xaxis = InterpolatedUnivariateSpline(x_data, y_data_xaxis, k=1)
        if s_xaxis(0.) > first_date_zero:
            df_insert[xaxis] = s_xaxis(0.)
            if xaxis != 'degree_days':
                y_data_dd = df.iloc[indx[:2]]['degree_days'].values
                s_xaxis = InterpolatedUnivariateSpline(x_data, y_data_dd, k=1)
                df_insert['degree_days'] = s_xaxis(0.)
            #df_insert['datetime'] = np.argmin(np.abs(data['degree_days'] - s_xaxis(0.)))
            return pd.concat([df, df_insert]).sort('datetime')
        else:
            return df
        
def plot_by_leaf_sample(data, weather, variable = 'severity', 
                        leaf = 1, plant = None,  xaxis = 'degree_days', variety = 'Tremie 12',
                        return_df = False, display_fit_function = None, ax = None, 
                        fixed_color = None, change_zero = False, marker = None, 
                        linestyle = '-', xlims = None, ylims = None, 
                        xlabel = None, ylabel = None):
    """ Plot given leaf (numbered from top) as individual samples for all plants or for given plant """
    try:
        import mpld3
        mpld3.disable_notebook()
    except:
        pass
    
    # Rename xaxis according to column names
    if not xaxis in ['degree_days', 'age_leaf', 'age_leaf_lig',
                        'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg']:
        raise ValueError("xaxis unknown: try 'degree_days', 'age_leaf', 'age_leaf_lig', 'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg'")

    # Isolate data for particular leaf
    df = data.copy()
    change_index(df, ['plant'])
    df = df[df['num_leaf_top']==leaf]
    
    # Set color cycle
    if ax == None:
        fig, ax = plt.subplots(figsize=(8,6))
    colors = ax._get_lines.color_cycle
    
    # Isolate data for particular leaf of particular plant
    if plant == None:
        plant = set(df.index)
        title = variety + ' all leaves #' + str(leaf)
    elif not is_iterable(plant):
        plant = [plant]
        title = variety + ' leaf ' + str(leaf) + ' on plant ' + str(plant)
    else:
        plant = [p for p in plant if p in set(df.index)]
        title = variety + ' leaf ' + str(leaf) + ' on plant ' + str([p for p in plant])

    for pl in plant:     
        df_ = df.loc[pl,:]
        if change_zero == True:
            df_ = change_zero_sample(data, df_, variable = variable, xaxis = xaxis)

        if fixed_color == None:
            color = next(colors)
        else:
            color = fixed_color
            
        # Draw plot
        ax.plot(df_[xaxis], df_[variable], color = color, marker = marker, linestyle = linestyle)
        ax.set_title(title, fontsize = 18)

        # Display particular fit
        if display_fit_function is not None:
            try:
                popt, pcov = find_fit(df_[xaxis], df_[variable], function_name = display_fit_function)
                if display_fit_function == 'logi':
                    ax.plot(min(df_[xaxis]) + np.arange(0,700), logi(np.arange(0,700), *popt), color = color)
                elif display_fit_function == 'gomp':
                    ax.plot(min(df_[xaxis]) + np.arange(0,700), gomp(np.arange(0,700), *popt), color = color)
            except:
                pass
    
    # Customize plot
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize = 18)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize = 18)
    
    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)
    else:
        ax.set_ylim([0, 100])
            
    # Display particular data table
    if return_df == True and plant is not None and len(plant)==1:
        return df.loc[plant]
        
def plot_by_leaf_samples_by_fnl(data, weather, variable = 'septo_green', 
                                leaf = 2, xaxis = 'age_emergence', 
                                marker = None, xlims = None, ylims = None):
    """ Plot given leaf (numbered from top) as individual samples
            for all plants or for given plant grouped by final leaf number """
    data_grouped = group_by_fnl(data)
    fig, ax = plt.subplots(figsize = (8,6))
    colors = iter(['b', 'r', 'k'])
    legend = []
    handles = []
    for fnl, df in data_grouped.iteritems():
        color = next(colors)
        plot_by_leaf_sample(df, weather, variable = variable, leaf = leaf,
                            xaxis = xaxis, ax = ax, fixed_color = color,
                            marker = marker, xlims = xlims, ylims = ylims)
        handles += [plt.Line2D([],[], color = color, linewidth = 2)]
        legend += ['%d' % fnl]
    ax.legend(handles, legend, title = 'FNL', loc='best', fontsize=12)

# Plotting functions for group of leaves ###########################################################
def plot_grouped_leaves(data, weather, variable='severity', degree_days=True,
                        alpha=0.05, fig_size=(8,6), title_suffix='', ylims = [0, 100],
                        title = None, return_descriptor = False, comparison_table=None, 
                        minimum_sample_size = 5, ax = None, fixed_color = None, marker = None, 
                        linestyle = '-', filling_contour = True, display_box = True):
    """ Plot mean of given variable over leaves in argument on all plants """
    # Format DataFrame to get mean and standard error
    df_mean = get_mean(data, column=variable, by_leaf=False, xaxis = 'datetime')
    df_low, df_high = get_confidence_interval(data, weather, column=variable, by_leaf=False)
    
    if return_descriptor == True:
        df_mean.columns = ['mean_severity']
        df_interval = (df_high - df_low)/2
        df_interval.columns = ['error_margin']
        descriptor = df_mean.join(df_interval)
        descriptor = add_index_ddays(descriptor, weather)
    
    if degree_days == True:
        idx = 'Degree days'
        (df_mean, df_high, df_low) = map(lambda x: add_index_ddays(x, weather).reset_index(level=0, drop = True),
                                         (df_mean, df_high, df_low))
        if comparison_table is not None and not 'Degree days' in comparison_table.index.names:
            comparison_table = add_index_ddays(comparison_table, weather)
            comparison_table.reset_index(level=0, drop = True, inplace=True)
    else:
        idx = 'Date'
    
    # Take out samples with size < minimum_sample_size
    count = table_count_notations(data, weather, variable = variable, add_ddays = degree_days)
    count = count[(minimum_sample_size<=count) & (count<=1000)]
    count = count.dropna(how='all')
    df_mean = df_mean.ix[count.index.get_level_values(idx),:]
    df_high = df_high.ix[count.index.get_level_values(idx),:]
    df_low = df_low.ix[count.index.get_level_values(idx),:]
    if return_descriptor == True:
        descriptor = descriptor[descriptor.index.get_level_values(idx).isin(count.index.get_level_values(idx))]

    # Plot mean and standard error
    if ax == None:
        fig, ax = plt.subplots(figsize=(8,6))
        
    if fixed_color == None:
        color = 'k'
        filling_color = 'grey'
    else:
        color = fixed_color
        filling_color = fixed_color

    if marker == None:
        marker = 'o'
    
    filling_linewidth = 1. if filling_contour == True else 0.
    
    if degree_days==True:
        h1 = ax.fill_between([ind for ind in df_low.index], df_low.values.flatten(), df_high.values.flatten(), 
                             linestyle='dashed', linewidth = filling_linewidth, color = color, facecolor = filling_color, alpha = 0.3)
    else:
        h1 = ax.fill_between(df_low.index, df_low.values.flatten(), df_high.values.flatten(), 
                             linewidth = filling_linewidth, color = color, facecolor = filling_color, alpha = 0.3)
    h2 = ax.plot(df_mean.index, df_mean, color = color, marker = marker, linestyle = linestyle)
    
    # Format DataFrame to draw boxplots
    sev = pd.DataFrame(data[variable])
    sev['subindex'] = sev.groupby(level=0).cumcount()
    # sev.set_index('subindex', append=True, inplace=True)
    # sev.index = sev.index.droplevel(1)
    change_index(sev, ['datetime', 'subindex'])
    sev = sev.unstack('datetime')[variable]
    sev = sev.ix[:,count.index.get_level_values('Date')]
    if degree_days==True:
        sev.rename(columns={sev.columns[col]:df_mean.index[col] for col in range(len(sev.columns))}, inplace=True)    
        
    # Plot boxplots for each date of observation
    if display_box == True:
        if degree_days==True:
            bp = sev.boxplot(ax=ax, positions=sev.columns, widths=10)
        else:
            locs = [loc.toordinal() for loc in df_mean.index]
            bp = sev.boxplot(ax=ax, positions=locs)
        if fixed_color is not None:
            plt.setp(bp['boxes'], color=fixed_color)
            plt.setp(bp['medians'], color=fixed_color)
            plt.setp(bp['whiskers'], color=fixed_color)
            plt.setp(bp['fliers'], color=fixed_color)
            plt.setp(bp['caps'], color=fixed_color)

    # Format plot: axis, legend and title
    if degree_days == True:
        xlims = [df_mean.index[0]-20, df_mean.index[-1]+20]
        ax.set_xlabel('Degree-days', fontsize=20)
        formatter = FuncFormatter(form_int)
        ax.xaxis.set_major_formatter(FuncFormatter(formatter))
    else:
        xlims = [df_mean.index[0]-timedelta(days=3), df_mean.index[-1]+timedelta(days=3)]
        month_fmt = DateFormatter('%b-%d')
        ax.xaxis.set_major_formatter(month_fmt)
        
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    proxy = [plt.Line2D((0,1),(0,0), color='k', marker='o', linestyle='-'), 
             plt.Rectangle((0,0), 0,0, facecolor=h1.get_facecolor()[0])]
    ax.legend(proxy, ['Mean', 'Confidence\n interval'], bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
    txtbox = TextArea('         Distribution', minimumdescent=False)
    ab1 = AnnotationBbox(txtbox, xy=(1.14, 0.75), xycoords="axes fraction", frameon=True, fontsize=32)
    ax.add_artist(ab1)
    img = read_png('box.png')
    imagebox = OffsetImage(img, zoom=.25)
    ab2 = AnnotationBbox(imagebox, xy=(1.06, 0.75), xycoords="axes fraction", frameon=False)
    ax.add_artist(ab2)
    if variable=='septo_green':
        ylab = 'Septoria on green (in % of leaf area)'
    elif variable=='severity':
        ylab = 'Severity (in % of leaf area)'
    else:
        ylab = variable+' (in % of leaf area)'
    yl = ax.set_ylabel(ylab, fontsize=20)
    num_leaves = np.unique(np.array(data['num_leaf_top']))
    leaves = 'F%d' % num_leaves[0] + ' to ' + 'F%d' % num_leaves[-1] if num_leaves[0]!=num_leaves[-1] else 'F%d' % num_leaves[0]
    if title == None:
        title = data['variety'][0]+'_'+str(data.index[0].year)+': '+leaves+title_suffix
    ax.set_title(title, fontsize=18)
    
    if comparison_table is not None:
        trues = comparison_table[comparison_table==True] * df_high[comparison_table==True] + 5
        h3 = ax.plot(trues.index+15, trues, 'k* ', markersize = 10)
    
    if return_descriptor==True:
        return descriptor
        
def plot_top_leaves(data, weather, leaves=[1, 2, 3, 4], variable='severity',
                    degree_days=True, alpha=0.05, fig_size=(8,6), title_suffix='', 
                    comparison_table=None, return_descriptor = False):
    """ Plot given variable averaged on 4 top leaves """
    # Keep only data values for leaves 1 to 4
    data_top = data[data['num_leaf_top'].isin(leaves)]
    # Plot
    df = plot_grouped_leaves(data_top, weather = weather, 
                                variable = variable, degree_days = degree_days,
                                alpha = alpha, fig_size = fig_size, 
                                title_suffix = title_suffix, comparison_table = comparison_table,
                                return_descriptor = return_descriptor)
    if return_descriptor==True:
        return df
        
def plot_lower_leaves(data, weather, leaves=[7, 8], variable='severity', 
                      degree_days=True, alpha=0.05, fig_size=(8,6), title_suffix='', 
                      comparison_table=None, return_descriptor = False):
    """ Plot given variable averaged on lower leaves """
    # Keep only data values for leaves 1 to 4
    data_top = data[data['num_leaf_top'].isin(leaves)]
    # Plot
    df = plot_grouped_leaves(data_top, weather = weather, 
                                variable = variable, degree_days = degree_days,
                                alpha = alpha, fig_size = fig_size, 
                                title_suffix = title_suffix, comparison_table = comparison_table,
                                return_descriptor = return_descriptor)
    if return_descriptor==True:
        return df

# Plotting functions for separate leaves layers averaged ###########################################
def get_data_of_interest(data, weather, variable = 'severity', xaxis = 'degree_days',
                          leaves = range(1,7), minimum_sample_size = 5, error_bars = False):
    """ Get average of given variable for separate leaves in argument versus xaxis in argument.
        If error_bars = True then also get 95% confidence interval """
    df_mean = get_mean(data, column=variable, by_leaf=True, xaxis = xaxis)
    if error_bars == True:
        df_low, df_high = get_confidence_interval(data, weather, column=variable,
                                                    xaxis = xaxis, by_leaf=True)
    
    # Count notations
    df_count = table_count_notations(data, weather, variable = variable, 
                                        xaxis = xaxis, add_ddays = True)
    df_count = df_count[(minimum_sample_size<=df_count) & (df_count<=1000)]
    
    # Drop columns of low leaves
    for col in df_mean.columns:
        if col not in leaves:
            (df_mean, df_count) = map(lambda x: x.drop(col, 1), (df_mean, df_count))
            if error_bars == True:
                (df_low, df_high) = map(lambda x: x.drop(col, 1), (df_low, df_high))
            
    # Take out data points with not enough individuals in sample
    mask = df_count.applymap(lambda x: np.isnan(x))
    df_mean = df_mean[~mask]
    if error_bars == True:
        (df_low, df_high) = map(lambda x: x[~mask], (df_low, df_high))
        return df_mean, df_low, df_high
    else:
        return df_mean
        
def set_color_cycle(ax, df_mean, with_brewer = True):
    try:
        if with_brewer == True:
            import brewer2mpl
            #ax.set_color_cycle(brewer2mpl.get_map('Paired', 'qualitative', max(3, len(df_mean.columns)+1)).mpl_colors)
            ax.set_color_cycle(brewer2mpl.get_map('Dark2', 'qualitative', max(3, len(df_mean.columns)+1)).mpl_colors)
        else:
            ax.set_color_cycle = ax._get_lines.color_cycle
    except:
        ax.set_color_cycle = ax._get_lines.color_cycle
    return ax

def get_fig_lims(df_mean, degree_days = True):
    if degree_days == True:
        xlims = [df_mean.index[0]-20, df_mean.index[-1]+20]
    else:
        xlims = [df_mean.index[0]-timedelta(days=3), df_mean.index[-1]+timedelta(days=3)]
    ylims = [0, 100]
    return xlims, ylims

def custom_ylabel(ax, variable):
    if variable=='septo_green':
        ylab = 'Septoria on green (in % of leaf area)'
    elif variable=='severity':
        ylab = 'Severity (in % of leaf area)'
    else:
        ylab = variable+' (in % of leaf area)'
    ax.set_ylabel(ylab, fontsize=20)
    return ax

def custom_xlabel(ax, degree_days):
    if degree_days == True:
        ax.set_xlabel('Degree-days', fontsize=20)
    else:
        month_fmt = DateFormatter('%b-%d')
        ax.xaxis.set_major_formatter(month_fmt)
    ax.tick_params(labelsize=12)
    return ax

def set_pointer(fig, points, df_mean, degree_days):
    """ Add pointer tooltip on graphs in ipython notebook"""
    try:
        import mpld3
        mpld3.enable_notebook()
        for leaf in df_mean.columns:
            if degree_days == True:
                labels = [(str(round(t,0))+': '+str(round(df_mean[leaf][t],1))+' %') for t in df_mean.index[~np.isnan(df_mean.ix[:,leaf])]]
            else:
                labels = [(date.strftime('%d-%b')+': '+str(round(df_mean[leaf][date],1))+' %') for date in df_mean.index[~np.isnan(df_mean.ix[:,leaf])]]
            tooltip = mpld3.plugins.PointLabelTooltip(points[leaf-1], labels)
            mpld3.plugins.connect(fig, tooltip)
    except:
        pass
        
def plot_df_mean(ax, df_mean, color = None, marker='d', markersize=8, 
                 linestyle='-', linewidth=2, df_low = None, df_high = None):
    """ Plot average of given variable for separate leaves in 
            df_mean columns versus xaxis in df_mean index """
    points = []
    for leaf in df_mean.columns:
        x_data = df_mean.index[~np.isnan(df_mean.ix[:,leaf])]
        y_data = df_mean.ix[:,leaf][~np.isnan(df_mean.ix[:,leaf])]
        options = dict(marker = marker, markersize = markersize, 
                        linestyle = linestyle, linewidth = linewidth)
        if color is None:
            if df_low is None:
                points += ax.plot(x_data, y_data, **options)
            else:
                points += ax.errorbar(x_data, y_data, 
                                      yerr = [y_data - df_low.ix[:,leaf][~np.isnan(df_mean.ix[:,leaf])], 
                                              df_high.ix[:,leaf][~np.isnan(df_mean.ix[:,leaf])]-y_data], 
                                      **options)
        else:
            if df_low is None:
                points += ax.plot(x_data, y_data, color = color, **options)
            else:
                points += ax.errorbar(x_data, y_data, 
                                      yerr = [y_data - df_low.ix[:,leaf][~np.isnan(df_mean.ix[:,leaf])], 
                                              df_high.ix[:,leaf][~np.isnan(df_mean.ix[:,leaf])]-y_data], 
                                      color = color, **options)
    return points, ax
    
def plot_by_leaf(data, weather, variable='severity', xaxis = 'degree_days', leaves = range(1,7), 
                 fig_size=(8,6), pointer=True, title_suffix='', minimum_sample_size = 5, 
                 linestyle = '-', marker = 'd', ax = None, fig = None, jump_colors = 0,
                 return_df_mean = False, xlims = None, ylims = None, title = None, 
                 error_bars = False, with_brewer = True, xlabel = True):
    """ Plot average of given variable for separate leaves in argument versus xaxis in argument """
    # Extract data
    if error_bars == True:
        df_mean, df_low, df_high = get_data_of_interest(data = data, weather = weather, 
                                                        xaxis = xaxis, variable = variable, 
                                                        leaves = leaves,
                                                        minimum_sample_size = minimum_sample_size,
                                                        error_bars = error_bars)
    else:
        df_mean = get_data_of_interest(data = data, weather = weather, xaxis = xaxis,
                                       variable = variable,
                                       leaves = leaves, 
                                       minimum_sample_size = minimum_sample_size, 
                                       error_bars = error_bars)
        df_low = None
        df_high = None
    
    # Set plot
    if ax == None:
        fig, ax = plt.subplots(figsize=(8,6))
    ax = set_color_cycle(ax = ax, df_mean = df_mean, with_brewer = with_brewer)
    colors = ax._get_lines.color_cycle
    for jump in range(jump_colors):
        next(colors)        
        
    # Draw plot
    points, ax = plot_df_mean(ax, df_mean, marker = marker, 
                              markersize = 8, linestyle = linestyle, 
                              linewidth = 2, df_low = df_low, df_high = df_high)
    
    # Customize plot
    leg = ax.legend(df_mean.columns, title='Leaf number', loc='best', fontsize=12)
    plt.setp(leg.get_title(), fontsize=14)
    ax = custom_ylabel(ax = ax, variable = variable)
    if xlabel == True:
        ax = custom_xlabel(ax, degree_days = True)
    if ylims is not None:
        ax.set_ylim(ylims)
    if xlims is None:
        xlims, ylims_ = get_fig_lims(df_mean = df_mean, degree_days = True)
    else:
        _ , ylims_ = get_fig_lims(df_mean = df_mean, degree_days = True)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims_)
    if ylims is not None:
        ax.set_ylim(ylims)
    ax.grid(alpha=0.5)
    if title is None:
        title = data['variety'][0]+'_'+str(data.index[-1].year)+title_suffix
    else:
        title += title_suffix
    ax.set_title(title, fontsize=18)
    if pointer==True:
        set_pointer(fig = fig, points = points, df_mean = df_mean, degree_days = True)
        
    if return_df_mean == True:
        return df_mean
        
def plot_by_leaf_by_fnl(data, weather, variable='severity', xaxis = 'degree_days',
                        leaves = range(1,7), fig_size=(8,6), xlims = None, ylims = None,
                        title_suffix='_control', minimum_sample_size = 5, 
                        with_brewer = True, error_bars = False, fig = None, ax = None,
                        xlabel = False, jumps = {11:0, 12:0, 13:0}):
    """ Plot mean of given variable for separate leaves in argument versus xaxis in argument 
        grouped by final leaf number """
    # Group data by fnl
    data_grouped = group_by_fnl(data)
    
    # Plot
    if ax is None:
        fig, ax = plt.subplots(figsize = fig_size)
    linestyles = iter(['-', '--', ':'])
    markers = iter(['o', 'd', '*'])
    labels = []
    for fnl, df in data_grouped.iteritems():
        (marker, linestyle) = map(lambda x: next(x), (markers, linestyles))
        df_mean = plot_by_leaf(df, weather, leaves = leaves, variable = variable, pointer = False,
                                title_suffix=title_suffix, xaxis = xaxis, ax = ax, fig = fig, 
                                linestyle = linestyle, marker = marker, jump_colors = jumps[fnl],
                                return_df_mean = True, xlims = xlims, ylims = ylims, 
                                with_brewer = with_brewer, xlabel = xlabel, 
                                minimum_sample_size = minimum_sample_size, 
                                error_bars = error_bars)
        labels += ['L '+str(lf)+': FNL '+str(fnl) for lf in df_mean.columns]

    # Customize legend
    color_list = [next(ax._get_lines.color_cycle) for col in range(7)]
    leg = ax.legend(ax._get_legend_handles(), labels, loc='center left', bbox_to_anchor=(1, 0.5))
    labels, handles = zip(*sorted(zip(labels, ax._get_legend_handles()), 
                            key=lambda t: color_list.index(t[1].get_c())))
    leg = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))

def compare_plot_by_leaf_two_datasets(data1, data2, weather, variable='severity',
                                      xaxis = 'degree_days', leaves = range(1,7),
                                      labels = ['data1', 'data2'], xlims = None, ylims = None,
                                      minimum_sample_size = 5, error_bars = False,
                                      title_suffix='_control', markers = ['d', '*'],
                                      linestyles = ['-', '--'], ax = None):
    if ax == None:
        fig, ax = plt.subplots(1, figsize = (8,6))
    plot_by_leaf(data1, weather, leaves = leaves, variable = variable, xaxis = xaxis, pointer=False, 
                 title_suffix = title_suffix, xlims = xlims, ylims = ylims, ax = ax,
                 minimum_sample_size = minimum_sample_size, error_bars = error_bars,
                 marker = markers[0], linestyle = linestyles[0])
    plot_by_leaf(data2, weather, leaves = leaves, variable = variable, xaxis = xaxis, pointer=False, 
                 title_suffix = title_suffix, xlims = xlims, ylims = ylims, ax = ax,
                 minimum_sample_size = minimum_sample_size, error_bars = error_bars,
                 marker = markers[1], linestyle = linestyles[1])
                 
    labels = ['F%d ' %lf +labels[0] for lf in leaves]+['F%d ' %lf +labels[1] for lf in leaves]
    color_list = [next(ax._get_lines.color_cycle) for lf in range(7)]
    leg = ax.legend(ax._get_legend_handles(), labels, loc='center left', bbox_to_anchor=(1, 0.5))
    labels, handles = zip(*sorted(zip(labels, ax._get_legend_handles()), 
                           key=lambda t: color_list.index(t[1].get_c())))
    leg = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    
def compare_plot_by_leaf_two_variables(data, weather, 
                                        variable1 = 'septo_green', variable2 = 'severity',
                                        xaxis = 'degree_days', leaves = range(1,7),
                                        labels = ['septo_green', 'severity'], 
                                        xlims = None, ylims = None,
                                        markers = ['d', '*'], linestyles = ['-', '--'],
                                        minimum_sample_size = 5, error_bars = False,
                                        title_suffix='_control', ax = None):
    if ax == None:
        fig, ax = plt.subplots(1, figsize = (8,6))
    plot_by_leaf(data, weather, leaves = leaves, variable = variable1, xaxis = xaxis, pointer=False, 
                 title_suffix = title_suffix, xlims = xlims, ylims = ylims, ax = ax,
                 minimum_sample_size = minimum_sample_size, error_bars = error_bars,
                 marker = markers[0], linestyle = linestyles[0])
    plot_by_leaf(data, weather, leaves = leaves, variable = variable2, xaxis = xaxis, pointer=False, 
                 title_suffix = title_suffix, xlims = xlims, ylims = ylims, ax = ax,
                 minimum_sample_size = minimum_sample_size, error_bars = error_bars,
                 marker = markers[1], linestyle = linestyles[1])
                 
    labels = ['F%d ' %lf +labels[0] for lf in leaves]+['F%d ' %lf +labels[1] for lf in leaves]
    color_list = [next(ax._get_lines.color_cycle) for lf in range(7)]
    leg = ax.legend(ax._get_legend_handles(), labels, loc='center left', bbox_to_anchor=(1, 0.5))
    labels, handles = zip(*sorted(zip(labels, ax._get_legend_handles()), 
                           key=lambda t: color_list.index(t[1].get_c())))
    leg = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    
def plot_confidence_and_boxplot(data, weather, variable='severity', xaxis = 'degree_days', 
                                leaves=range(1,7), fig_size=(15, 10), title_suffix='', 
                                comparison_table = None, minimum_sample_size = 5, 
                                return_fig = False, xlims = None, ylims = None, fig = None, 
                                axs = None, marker = 'o', fixed_color = None, alpha=0.3, 
                                linestyle = '-', filling_contour = True, tight_layout = True, 
                                title = True, xlabel = True, display_error_margin = True,
                                display_box = True, forced_date_leaves = None):
    """ Plot mean of given variable, with distribution and confidence interval
        for separate leaves in argument versus xaxis in argument """
    from matplotlib.dates import MonthLocator, DateFormatter
    from math import ceil
    try:
        import mpld3
        mpld3.disable_notebook()
    except:
        pass
    
    # All severity data and mean, by date and by leaf
    data_ = data.copy()
    change_index(data_, ['datetime', 'num_leaf_top'])
    df_sev = pd.DataFrame(data_[variable])
    df_mean = get_mean(data_, column=variable, by_leaf=True, xaxis = 'datetime')
    df_low, df_high = get_confidence_interval(data_, weather, column=variable, by_leaf=True)
    if xaxis != 'datetime':
        degree_days = True
        (df_mean, df_high, df_low) = map(lambda x: add_index_ddays(x, weather).reset_index(level=0, drop = True),
                                         (df_mean, df_high, df_low))
        
        if comparison_table is not None and not 'Degree days' in comparison_table.index.names:
            comparison_table = add_index_ddays(comparison_table, weather)
            comparison_table.reset_index(level=0, drop = True, inplace=True)
    else:
        degree_days = False
    
    if comparison_table is not None:
        signif_diffs = comparison_table[comparison_table==True] * df_high[comparison_table==True] + 5
        
    # Count notations
    count = table_count_notations(data_, weather, variable = variable, add_ddays = degree_days)
    count = count[(minimum_sample_size<=count) & (count<=1000)]
    
    # Prepare figure
    nb_leaves = min(df_mean.shape[1], len(leaves))
    if axs == None:
        fig, axs = plt.subplots(int(ceil(nb_leaves/2.)), 2, figsize=fig_size)
        if len(axs.flat)>nb_leaves:
            fig.delaxes(axs.flat[-1])
    
    for ax, leaf in zip([ax for i,ax in enumerate(axs.flat) if i<len(leaves)], leaves):
        # Get severity data and mean for chosen leaf
        sev = df_sev.xs(leaf, level=1)
        mn = df_mean[leaf]
        ind_enough_data = np.where(count[leaf]>=minimum_sample_size)
        sev = sev[sev.index.isin(count.index.get_level_values('Date')[ind_enough_data[0]])]
        notnans = np.where(mn.values==mn.values) and ind_enough_data[0]  
        
        # Reshape dataframe to draw boxplots from severity data
        sev['subindex'] = sev.groupby(level=0).cumcount()
        sev.set_index('subindex', append=True, inplace=True)
        sev = sev.unstack('datetime')[variable]
        if degree_days==True:
            sev.rename(columns={sev.columns[col]:mn.iloc[notnans].index[col] for col in range(len(sev.columns))}, inplace=True)
            
        # Set plot parameters
        if fixed_color == None:
            color = 'k'
            filling_color = 'grey'
        else:
            color = fixed_color
            filling_color = fixed_color

        if marker == None:
            marker = 'o'
        
        filling_linewidth = 1. if filling_contour == True else 0.
        
        # Draw figure
        ax.annotate('Leaf %d' % leaf, xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
        if len(notnans>0):
            if xaxis != 'datetime':
                # Get data for x axis
                x_data = np.array([ind for ind in df_low.iloc[notnans,df_low.columns==leaf].index])
                if xaxis in ['age_leaf', 'age_leaf_lig', 'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg']:
                    df_dates = pd.DataFrame(index = set(data.num_leaf_top))
                    if forced_date_leaves is None:
                        date_name = get_date_name_event(xaxis)
                        df_dates[xaxis] = map(lambda x: np.mean(data[date_name][data.num_leaf_top==x]), 
                                              df_dates.index)
                
                    else:
                        date_name = forced_date_leaves.columns[-1]
                        df_dates[xaxis] = map(lambda x: np.mean(forced_date_leaves[date_name][forced_date_leaves.num_leaf_top==x]), 
                                              df_dates.index)
                    x_data -= int(df_dates.loc[leaf, xaxis])
                    
                # Draw plot
                if display_error_margin == True:
                    h1 = ax.fill_between(x_data, 
                                         df_low.iloc[notnans, df_low.columns==leaf].iloc[:,0], 
                                         df_high.iloc[notnans, df_high.columns==leaf].iloc[:,0], 
                                         linestyle='dashed', linewidth = filling_linewidth,
                                         color = color, facecolor = filling_color, alpha = alpha)
                h2 = ax.plot(x_data, mn.iloc[notnans], color = color, 
                             marker = marker, linestyle = linestyle)

                # Customize and draw boxplot
                if display_box == True:
                    if fixed_color is not None:
                        bp = sev.boxplot(ax = ax, positions = x_data, widths = 10)
                        plt.setp(bp['boxes'], color=fixed_color)
                        plt.setp(bp['medians'], color=fixed_color)
                        plt.setp(bp['whiskers'], color=fixed_color)
                        plt.setp(bp['fliers'], color=fixed_color)
                        plt.setp(bp['caps'], color=fixed_color)
                    else:
                        sev.boxplot(ax=ax, positions=x_data, widths=10)
                if xlabel == True:
                    ax.set_xlabel('Degree-days', fontsize=16)
                formatter = FuncFormatter(form_int)
                ax.xaxis.set_major_formatter(FuncFormatter(formatter))
                if xlims == None:
                    xlims = [x_data[0]-20, x_data[-1]+20]
            else:
                # Draw plot
                if display_error_margin == True:
                    h1 = ax.fill_between(df_low.iloc[notnans,df_low.columns==leaf].index, 
                                         df_low.iloc[notnans,df_low.columns==leaf].iloc[:,0], 
                                         df_high.iloc[notnans,df_high.columns==leaf].iloc[:,0], 
                                         linestyle='dashed', color = color, 
                                         facecolor = filling_color, alpha = alpha)
                h2 = ax.plot(mn.index[notnans], mn.iloc[notnans], color = color, 
                                marker = marker, linestyle = linestyle)
                
                ll = [loc.toordinal() for loc in mn.index]
                locs = [ll[ind] for ind in notnans]
                
                # Customize and draw boxplot
                if display_box == True:
                    if fixed_color is not None:
                        bp = sev.boxplot(ax = ax, positions = locs, widths = 10)
                        plt.setp(bp['boxes'], color=fixed_color)
                        plt.setp(bp['medians'], color=fixed_color)
                        plt.setp(bp['whiskers'], color=fixed_color)
                        plt.setp(bp['fliers'], color=fixed_color)
                        plt.setp(bp['caps'], color=fixed_color)
                    else:
                        sev.boxplot(ax=ax, positions=locs, widths=10)
                month_fmt = DateFormatter('%b-%d')
                ax.xaxis.set_major_formatter(month_fmt)
                xlims = [df_mean.index[0]-timedelta(days=3), df_mean.index[-1]+timedelta(days=3)]

            # Set legend
            if len(axs)>1:
                if (isinstance(axs[0], np.ndarray) and ax==axs[0][1]) or (not 
                    isinstance(axs[0], np.ndarray) and ax==axs[1]):
                    proxy = [plt.Line2D((0,1),(0,0), color=color, 
                                marker=marker, linestyle=linestyle)]
                    labels = ['Mean']
                    if display_error_margin == True:
                        proxy += [plt.Rectangle((0,0), 0,0, facecolor=h1.get_facecolor()[0])]
                        labels += ['Confidence\n interval']
                    if fig_size==(8,6):
                        ax.legend(proxy, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                        locbox = (1.1, 0.4)
                        locim = (1.12, 0.4)
                    else:
                        ax.legend(proxy, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
                        locbox = (1.04, 0.5)
                        locim = (1.05, 0.5)

                    txtbox = TextArea('        Distribution', minimumdescent=False)
                    ab1 = AnnotationBbox(txtbox, xy=locbox, xycoords="axes fraction", 
                                         frameon=True, fontsize=32, box_alignment=(0., 0.5))
                    ax.add_artist(ab1)
                    try:
                        img = read_png('box.png')
                        imagebox = OffsetImage(img, zoom=.25)
                        ab2 = AnnotationBbox(imagebox, xy=locim, 
                                             xycoords="axes fraction", frameon=False)         
                        ax.add_artist(ab2)
                    except:
                        pass
        
        # Customize
        ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)
        else:
            ax.set_ylim([0, 100])
        if ax in axs.flat[::2]:
            if variable=='septo_green':
                ylab = 'Septoria on green (in %)'
            elif variable=='severity':
                ylab = 'Severity (in %)'
            else:
                ylab = variable+' (in %)'
            yl = ax.set_ylabel(ylab, fontsize=14)
        if len(notnans)<1:
            fig.delaxes(ax)
            
        if comparison_table is not None:
            h3 = ax.plot(signif_diffs.index+15, signif_diffs[leaf], 'k* ', markersize = 8)

    if title == True:
        plt.text(0.5, 0.98, data['variety'][0]+'_'+str(data.index[-1].year)+title_suffix, 
                fontsize=18, transform=fig.transFigure, horizontalalignment='center')
    
    if tight_layout == True:
        fig.tight_layout(rect=(0,0,0.98,0.98))
    
    if return_fig == True:
        return df_mean, df_high, df_low, fig, axs
    else:
        return df_mean, df_high, df_low
        
def plot_confidence_and_boxplot_by_fnl(data, weather, variable='severity', xaxis = 'degree_days',
                                       leaves=range(1,7), fig_size=(15, 10), title_suffix='', 
                                       minimum_sample_size = 5, xlims = None, ylims = None):
    """ Plot mean of given variable, with distribution and confidence interval
        for separate leaves in argument versus xaxis in argument, plants grouped by FNL """
    # Group data by fnl
    data_grouped = group_by_fnl(data)
    
    # Plot
    colors = iter(['b', 'r', 'g'])
    markers = iter(['o', 'd', '*'])
    linestyles = iter(['-', '--', ':'])
    labels = []
    proxy = []
    count = 0
    alpha = 0.2
    for fnl, df in data_grouped.iteritems():
        (color, marker, linestyle) = map(lambda x: next(x), (colors, markers, linestyles))
        if count == 0:
            (fig, axs) = map(lambda x: None, range(2))
            count += 1
        df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(df, weather, leaves = leaves, variable = variable,
                                                                        xaxis = xaxis, title_suffix='_control', fig_size = fig_size,
                                                                        xlims = xlims, ylims = ylims, return_fig = True, fig = fig,
                                                                        axs = axs, marker = marker, fixed_color = color, alpha = alpha,
                                                                        linestyle = linestyle, filling_contour = False,
                                                                        minimum_sample_size = minimum_sample_size)
        labels += ['FNL '+str(fnl)]
        proxy += [plt.Rectangle((0,0), 0,0, facecolor=color, alpha = alpha)]
    
    # Set legend
    axs[1][1].legend(proxy, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    
def compare_confidence_and_boxplot_green_sev(data, weather, xaxis = 'degree_days', 
                                             leaves=range(1,7), fig_size=(15, 10), 
                                             title_suffix='', minimum_sample_size = 5, 
                                             xlims = None, ylims = None, display_necro = False, 
                                             display_green = True):
    """ Compare plots of mean, distribution and confidence interval 
            for variables 'septo_green' and 'severity' on separate leaves in data """
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data, weather, leaves = leaves, variable = 'septo_green',
                                                                        xaxis = xaxis, title_suffix='_control', fig_size = fig_size,
                                                                        xlims = xlims, ylims = ylims, return_fig = True, 
                                                                        marker = 'o', fixed_color = 'g', alpha = 0.2,
                                                                        linestyle = '-', filling_contour = False)
    labels = ['Septoria on green']
    proxy = [plt.Rectangle((0,0), 0,0, facecolor='g', alpha = 0.2)]
    
    if display_necro == True:
        df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data, weather, leaves = leaves, variable = 'necro',
                                                                    xaxis = xaxis, title_suffix='_control', fig_size = fig_size,
                                                                    xlims = xlims, ylims = ylims, return_fig = True, fig = fig, 
                                                                    axs = axs, marker = 'o', fixed_color = 'k', alpha = 0.2,
                                                                    linestyle = '--', filling_contour = False,
                                                                    display_error_margin = False, display_box = False)
        labels.insert(0, 'Apical senescence')
        proxy.insert(0, plt.Line2D((0,1),(0,0), color = 'k', marker ='o', linestyle ='--'))
        
    if display_green == True:
        df = data.copy()
        df['green'] = 100 - df['necro']
        df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data, weather, leaves = leaves, variable = 'green',
                                                                    xaxis = xaxis, title_suffix='_control', fig_size = fig_size,
                                                                    xlims = xlims, ylims = ylims, return_fig = True, fig = fig,
                                                                    axs = axs, marker = 'o', fixed_color = 'k', alpha = 0.2,
                                                                    linestyle = '--', filling_contour = False,
                                                                    display_error_margin = False, display_box = False)
        labels.insert(0, 'Non-senescent')
        proxy.insert(0, plt.Line2D((0,1),(0,0), color = 'k', marker ='o', linestyle ='--'))
        
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data, weather, leaves = leaves, variable = 'severity',
                                                                        xaxis = xaxis, title_suffix='_control', fig_size = fig_size,
                                                                        xlims = xlims, ylims = ylims, return_fig = True, fig = fig, 
                                                                        axs = axs, marker = 'o', fixed_color = 'r', alpha = 0.2,
                                                                        linestyle = '-', filling_contour = False)
    labels += ['Severity']
    proxy += [plt.Rectangle((0,0), 0,0, facecolor='r', alpha = 0.2)]   
        
    axs[1][1].legend(proxy, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    
    for ax in axs.flat:
        ax.set_ylabel('Percentage of leaf surface', fontsize = 14)
        
def plot_confidence_and_boxplot_with_global(data, data_global, weather, xaxis = 'degree_days',
                                            variable = 'severity', leaves=range(1,7), 
                                            fig_size=(15, 10), title_suffix='', title = '',
                                            minimum_sample_size = 5, xlims = None, ylims = None,
                                            display_necro = False, display_green = True):
    """ Plot mean of given variable, with distribution and confidence interval
        for separate leaves versus xaxis. 
        
        Data from global notations are plotted against measurements on tagged plants.
    """
    fig, axs = plt.subplots(int(ceil(len(leaves)/2.)), 2, figsize=fig_size)
    
    if xaxis in ['age_leaf',  'age_leaf_lig', 'age_leaf_vs_flag_lig', 'age_leaf_vs_flag_emg']:
        forced_date_leaves = data[['num_leaf_top', get_date_name_event(xaxis)]]
    else:
        forced_date_leaves = None
    
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data_global, weather, leaves = range(1,min(max(leaves),5)),
                                                                     variable = variable, xaxis = xaxis, 
                                                                     title_suffix = title_suffix, title = title,
                                                                     fig_size = fig_size, xlims = xlims, ylims = None,
                                                                     return_fig = True, fig = fig, axs = axs, marker = 'o', 
                                                                     fixed_color = 'r', alpha = 0.2, linestyle = '', 
                                                                     filling_contour = False, display_error_margin = False,
                                                                     display_box = False, forced_date_leaves = forced_date_leaves)
    
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data, weather, leaves = leaves, variable = variable,
                                                                        xaxis = xaxis, title_suffix = title_suffix, title = title,
                                                                        fig_size = fig_size, xlims = xlims, ylims = ylims, 
                                                                        return_fig = True, fig = fig, axs = axs, 
                                                                        minimum_sample_size = minimum_sample_size)
    labels = ['Weekly\n measurements', 'Consistency\n measurements']
    proxy = [plt.Line2D((0,1),(0,0), color = 'k', marker ='o', linestyle ='-'), 
             plt.Line2D((0,1),(0,0), color = 'r', marker ='o', linestyle ='')]
                
    try:
        axs[1][1].legend(proxy, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    except:
        pass
        
def compare_confidence_and_boxplot_two_datasets(data1, data2, weather, variable = 'severity',
                                                xaxis = 'degree_days', leaves=range(1,5),
                                                fig_size=(15, 10), minimum_sample_size = 5., 
                                                labels = ['data1', 'data2'], 
                                                xlims = None, ylims = None, 
                                                display_error_margin = [True, True],
                                                display_box = [True, True],
                                                forced_date_leaves = None):
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data1, weather, leaves = leaves, variable = variable,
                                                                        xaxis = xaxis, title = False, fig_size = fig_size,
                                                                        xlims = xlims, ylims = ylims, return_fig = True, 
                                                                        marker = 'o', fixed_color = 'g', alpha = 0.2,
                                                                        linestyle = '-', filling_contour = False,
                                                                        display_error_margin = display_error_margin[0],
                                                                        display_box = display_box[0],
                                                                        minimum_sample_size = minimum_sample_size, 
                                                                        forced_date_leaves = forced_date_leaves)
    proxy = [plt.Rectangle((0,0), 0,0, facecolor='g', alpha = 0.2)]
        
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data2, weather, leaves = leaves, variable = variable,
                                                                        xaxis = xaxis, title = False, fig_size = fig_size,
                                                                        xlims = xlims, ylims = ylims, return_fig = True,
                                                                        fig = fig, axs = axs,
                                                                        marker = 'o', fixed_color = 'r', alpha = 0.2,
                                                                        display_error_margin = display_error_margin[1],
                                                                        display_box = display_box[1],
                                                                        linestyle = '-', filling_contour = False, 
                                                                        minimum_sample_size = minimum_sample_size,
                                                                        forced_date_leaves = forced_date_leaves)
    proxy += [plt.Rectangle((0,0), 0,0, facecolor='r', alpha = 0.2)]   
        
    axs[1][1].legend(proxy, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    
# Fit of curves ####################################################################################
from scipy.optimize import curve_fit

def gomp(x, b, c):
    """ Gompertz law """
    #a = 100 upper asymptote
    #b = negative = x axis displacement
    #c = negative = growth rate
    return 100*(np.exp(b*(np.exp(c*x))))

def logi(x, b, c):
    """ Logistic law """
    #a = 100 upper asymptote
    #b = negative = x axis displacement
    #c = negative = growth rate
    return 100/(1+np.exp(-c*(x-b)))

def find_fit(degree_days, severity, function_name = 'logi'):
    ddays = degree_days - min(degree_days)
    if function_name == 'logi':
        return curve_fit(logi, ddays, severity, p0=[400, 0.1])
    elif function_name == 'gomp':
        return curve_fit(gomp, ddays, severity, p0=[-500, -0.01])
    else:
        print 'Fitting function unknown'
        
# Extraction of indicators #########################################################################
def get_date_threshold(data, weather, variable = 'severity', xaxis = 'degree_days', from_top = True,
                       threshold = 5, fit_function_name = None):
    """ Get date at which variable overpass given threshold """
    df = data.copy()
    change_index(df, ['plant'])
    
    if from_top == True:
        num_leaf = 'num_leaf_top'
    else:
        num_leaf = 'num_leaf_bottom'

    df_dates = pd.DataFrame()
    for lf in set(df[num_leaf]):
        for pl in set(df[df[num_leaf]==lf].index):
            df_ = df[df[num_leaf]==lf].loc[pl,:]
            if is_iterable(df_[variable]) and any(df_[variable] < threshold) and any(df_[variable] > threshold):
                if fit_function_name is None :
                    # Extract date when threshold is overpassed from interpolation between notations
                    df_2 = pd.Series([np.nan for i in np.arange(int(min(df_[xaxis]))-100,
                                                                int(max(df_[xaxis]))+100)],
                                                                index = np.arange(int(min(df_[xaxis]))-100,
                                                                int(max(df_[xaxis]))+100))
                    df_2[df_[xaxis][pd.notnull(df_[xaxis])].astype(int)] = df_[variable][pd.notnull(df_[xaxis])].values
                    df_2.interpolate(inplace=True)
                    df_dates.loc[pl,lf] = np.argmax(df_2>threshold)
                else:
                    try:
                        # Try and find best fit
                        popt, pcov = find_fit(df_[xaxis], df_[variable], 
                                                function_name = fit_function_name)
                        # Extract date when threshold is overpassed
                        if fit_function_name == 'logi':
                            df_dates.loc[pl,lf] = min(df_[xaxis]) + np.argmax(logi(np.arange(1000),*popt)>threshold)
                        elif fit_function_name == 'gomp':
                            df_dates.loc[pl,lf] = min(df_[xaxis]) + np.argmax(gomp(np.arange(1000),*popt)>threshold)
                    except:
                        pass
    return df_dates
    
def error_margin_threshold(df_dates):
    """ Get error margin of date at which threshold is overpassed """
    mini = df_dates.apply(lambda x: bootstr(x[~np.isnan(x).values])[0] 
                          if len(x)<20 
                          else x[~np.isnan(x).values].mean() - 2*st.sem(x[~np.isnan(x).values]))
    maxi = df_dates.apply(lambda x: bootstr(x[~np.isnan(x).values])[1]
                          if len(x)<20 
                          else x[~np.isnan(x).values].mean() + 2*st.sem(x[~np.isnan(x).values]))
    return mini, maxi

def std_dev_threshold(df_dates):
    """ Get standard deviation of date at which threshold is overpassed """
    mini = df_dates.mean() - df_dates.std()
    maxi = df_dates.mean() + df_dates.std()
    return mini, maxi

def plot_mean_thr_date_by_leaf(df_dates, ax = None, legend = None,
                               marker = 'd', linestyle = '', color = 'b',
                               fig_size=(8,6), invert_xaxis = True, error_bars = True, 
                               error_method = error_margin_threshold, 
                               leaf_dates = None, ylims = None):
    """ Plot date at which variable overpass given threshold """
    df = df_dates.mean()
    if ax == None:
        fig, ax = plt.subplots(1,1, figsize=fig_size)
    if error_bars == True:
        mini, maxi = error_method(df_dates)
        ax.errorbar(df.index, df.values, yerr = [df-mini, maxi-df], color = color, 
                    marker = marker, linestyle = linestyle, markersize = 8)
    else:
        ax.plot(df.index, df.values, color = color, marker = marker, 
                linestyle = linestyle, markersize = 8)
    
    # Display leaf dates if given
    if leaf_dates is not None:
        leaf_dates = leaf_dates.loc[df.index]
        ax.plot(leaf_dates.index, leaf_dates.values, marker = "^", color = 'g', linestyle = '')
    
    # Customize
    if ylims != None:
        ax.set_ylim(ylims)
    if invert_xaxis == True:
        ax.invert_xaxis()
    ax.set_xlabel('Leaf number', fontsize = 18)
    ax.set_ylabel('Threshold date (Cd)', fontsize = 18)
    if legend is not None:
        ax.legend(legend, loc = 'best')
        
def plot_mean_thr_date_by_leaf_by_fnl(data, weather, variable = 'severity',
                                      xaxis = 'degree_days', from_top = True, threshold = 3, 
                                      fit_function_name = None, fig_size=(8,6), 
                                      add_mean = True, error_method = error_margin_threshold):
    """ Plot date at which variable overpass given threshold from plants grouped by FNL """
    data_grouped = group_by_fnl(data)
    fig, ax = plt.subplots(figsize = fig_size)
    colors = iter(['b', 'g', 'r'])
    legend = []
    for fnl, df in data_grouped.iteritems():
        df_dates = get_date_threshold(df, weather, variable = variable, xaxis = xaxis,
                                      from_top = from_top, threshold = threshold,
                                      fit_function_name = fit_function_name)
        plot_mean_thr_date_by_leaf(df_dates, ax = ax, color = next(colors), 
                                   invert_xaxis = False, error_method = error_method)
        legend += ['%d' % fnl]
    
    if add_mean == True:
        df_dates = get_date_threshold(data, weather, variable = variable, xaxis = xaxis,
                                      from_top = from_top, threshold = threshold,
                                      fit_function_name = fit_function_name)
        plot_mean_thr_date_by_leaf(df_dates, ax = ax, marker = None, linestyle = '--', color = 'k', 
                                   invert_xaxis = False, error_method = error_method)
        legend += ['Mean']
    
    if from_top == True:
        ax.invert_xaxis()
    ax.legend(legend, title = 'FNL', loc='best', fontsize=12)
    
def compare_thr_dates_with_fit(data, weather, threshold = 5, fig_size=(8,6)):
    """ Compare plots of date at which variable overpass given threshold obtained 
        with linear, gompertz or logistic fit of data """
    df_dates = get_date_threshold(data, weather, threshold = threshold, fit_function_name = None)
    df_dates_logi = get_date_threshold(data, weather, threshold = threshold, 
                                        fit_function_name = 'logi')
    df_dates_gomp = get_date_threshold(data, weather, threshold = threshold, 
                                        fit_function_name = 'gomp')
    fig, ax = plt.subplots(1,1, figsize=fig_size)
    ax.plot(df_dates.mean().index, df_dates.mean().values)
    ax.plot(df_dates_logi.mean().index, df_dates_logi.mean().values)
    ax.plot(df_dates_gomp.mean().index, df_dates_gomp.mean().values)
    ax.invert_xaxis()
    ax.set_xlabel('Leaf number', fontsize = 18)
    ax.set_ylabel('Threshold date (Cd)', fontsize = 18)
    ax.legend(['Linear', 'Logistic', 'Gompertz'])
    
def get_value_at_age(data, weather, variable = 'severity', age = 300, fit_function_name = None):
    """ Get value of variable by leaf at age in argument """
    df = data.copy()
    change_index(df, ['plant'])
    
    df_dates = pd.DataFrame()
    for lf in set(df['num_leaf_top']):
        for pl in set(df[df['num_leaf_top']==lf].index):
            df_ = df[df['num_leaf_top']==lf].loc[pl,:]
            if is_iterable(df_[variable]) and any(df_.age_leaf < age) and any(df_.age_leaf > age):
                if fit_function_name is None :
                    # Extract severity at age from interpolation between notations
                    df_2 = pd.Series([np.nan for i in np.arange(0,2000)], index = np.arange(0,2000))
                    df_2.loc[np.round(df_['age_leaf'])] = df_[variable].values
                    df_2.interpolate(inplace=True)
                    df_dates.loc[pl,lf] = df_2[age]
                else:
                    try:
                        # Try and find best fit
                        popt, pcov = find_fit(df_.age_leaf, df_[variable], 
                                                function_name = fit_function_name)
                        # Extract severity at age
                        if fit_function_name == 'logi':
                            df_dates.loc[pl,lf] = logi(age,*popt)
                        elif fit_function_name == 'gomp':
                            df_dates.loc[pl,lf] = gomp(age, *popt)
                    except:
                        pass
    return df_dates

def plot_mean_value_at_age_by_leaf(df_value_at_age, ax = None, age = 300, 
                                   legend = None, marker = 'd', linestyle = '', 
                                   color = 'b', fig_size=(8,6), invert_xaxis = True, 
                                   error_bars = True, error_method = error_margin_threshold, 
                                   ylims = None):
    """ Plot value of variable by leaf at age in argument """
    df = df_value_at_age.mean()
    if ax == None:
        fig, ax = plt.subplots(1,1, figsize=fig_size)
    if error_bars == True:
        mini, maxi = error_method(df_value_at_age)
        ax.errorbar(df.index, df.values, yerr = [df-mini, maxi-df], color = color, 
                    marker = marker, linestyle = linestyle, markersize = 8)
    else:
        ax.plot(df.index, df.values, color = color, marker = marker, 
                linestyle = linestyle, markersize = 8)

    # Customize
    if ylims != None:
        ax.set_ylim(ylims)
    if invert_xaxis == True:
        ax.invert_xaxis()
    ax.set_xlabel('Leaf number', fontsize = 18)
    ax.set_ylabel('Severity at '+str(age)+' Cd', fontsize = 18)
    if legend is not None:
        ax.legend(legend, loc = 'best')
        
def plot_mean_value_at_age_by_leaf_by_fnl(data, weather, variable = 'severity', age = 300, 
                                          fit_function_name = None, fig_size=(8,6),
                                          add_mean = True, error_bars = True, 
                                          error_method = error_margin_threshold, ylims = None):
    """ Plot value of variable by leaf at age in argument, from plants grouped by FNL """
    data_grouped = group_by_fnl(data)
    fig, ax = plt.subplots(figsize = fig_size)
    colors = iter(['b', 'g', 'r'])
    legend = []
    for fnl, df in data_grouped.iteritems():
        df_value_at_age = get_value_at_age(df, weather, variable = variable,
                                            age = age, fit_function_name = fit_function_name)
        plot_mean_thr_date_by_leaf(df_value_at_age, ax = ax, 
                                   color = next(colors), invert_xaxis = False, 
                                   error_bars = error_bars, error_method = error_method, 
                                   ylims = ylims)
        legend += ['%d' % fnl]
    
    if add_mean == True:
        df_value_at_age = get_value_at_age(data, weather, variable = variable, age = age,
                                            fit_function_name = fit_function_name)
        plot_mean_thr_date_by_leaf(df_value_at_age, ax = ax, marker = None, 
                                    linestyle = '--', color = 'k', invert_xaxis = False)
        legend += ['Mean']
    
    ax.invert_xaxis()
    ax.legend(legend, title = 'FNL', loc='best', fontsize=12)
    
def get_nb_risk_by_leaf(adel, weather, age_inf = 0, age_sup = 500,
                        phyllochron_added = 0, unit = 'hours'):
    """ Get number of events at epidemiological risk during leaf life span """
    df_dates = get_leaf_dates(adel, event='tip', nff = None)
    df_risk = pd.DataFrame(index = df_dates.index, columns = ['nb_events'])
    for lf, date in df_dates.iterrows():
        df_ = weather.data.copy()
        df_['age_leaf'] = df_['degree_days'] - date.values
        df_ = df_[(df_['age_leaf']>age_inf) & (df_['age_leaf']<age_sup)]
        if unit == 'hours':
            df_risk.loc[lf] = df_['septo_infection_risk_with_event'].sum()
        elif unit == 'events':
            df_risk.loc[lf] = sum(map(lambda x: 0 if x<0 or x==np.nan else x, 
                                    df_['septo_infection_risk_with_event'].diff())[1:])
        else:
            raise ValueError("Unit unknown")
    return df_risk

def get_speed(data, weather, variable = 'severity', from_top = True):
    """ Get speed of curve to reach value max  """
    df = data.copy()
    change_index(df, ['plant'])

    if from_top == True:
        num_leaf = 'num_leaf_top'
    else:
        num_leaf = 'num_leaf_bottom'

    df_speed = pd.DataFrame()
    for lf in set(df[num_leaf]):
        for pl in set(df[df[num_leaf]==lf].index):
            df_ = df[df[num_leaf]==lf].loc[pl,:]
            if is_iterable(df_[variable]):
                # Extract date of first positive value from interpolation between notations
                change_zero_sample(data, df_, variable = variable, xaxis = 'degree_days')
                first_date_zero = df_['degree_days'][df_[variable]==0].iloc[0]
                # Extract date of first reach of max value 
                # (TODO from interpolation between notations)
                max_value = max(df_[variable])
                first_date_max = df_['degree_days'][df_[variable]==max_value].iloc[0]
                delta_date = first_date_max - first_date_zero
                df_speed.loc[pl,lf] = max_value/delta_date if delta_date>0. else 0.
    return df_speed
    
def get_max_value(data, weather, variable = 'severity', from_top = True):
    """ Get value max of curve on each sample """
    df = data.copy()
    change_index(df, ['plant'])

    if from_top == True:
        num_leaf = 'num_leaf_top'
    else:
        num_leaf = 'num_leaf_bottom'

    df_max = pd.DataFrame()
    for lf in set(df[num_leaf]):
        for pl in set(df[df[num_leaf]==lf].index):
            df_ = df[df[num_leaf]==lf].loc[pl,:]
            if is_iterable(df_[variable]):
                df_max.loc[pl,lf] = max(df_[variable])
    return df_max
    
# Plots related to weather-disease interactions ####################################################
import matplotlib.gridspec as gridspec
from copy import deepcopy

def form_tick(x, pos):
    t = date.fromordinal(int(x))
    return t.strftime('%b')+'\n'+str(int(weather.get_variable('degree_days', t)))
    
def plot_septo_infection_risk(weather, start_date="2010-10-15 12:00:00", 
                              leaf_dates=None, xaxis = 'degree_days', 
                              ax = None, xlims=None, only_with_event=True, 
                              title = None, xlabel = True, ylabel = False, 
                              display_first_event = True, wetness_duration_min = 10.,
                              temp_min = 10., temp_max = 30.):
    """ Plot events with risk of infection by septoria as in Alep model """
    if ax == None:
        fig, ax = plt.subplots(figsize=(8,1.8))
    
    if not 'degree_days' in weather.data.columns:
        weather.check(varnames=['degree_days'], 
                      models={'degree_days':linear_degree_days},
                      start_date=start_date, base_temp=0., max_temp=30.)
    if not 'wetness' in weather.data.columns:
        weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    if not 'septo_infection_risk' in weather.data.columns:
        weather.check(varnames=['septo_infection_risk'],
                      models={'septo_infection_risk':add_septoria_infection_risk}, 
                      wetness_duration_min = wetness_duration_min, 
                      temp_min = temp_min, temp_max = temp_max)
        
    if xaxis == 'degree_days':
        index = weather.data.degree_days
    elif xaxis == 'date':
        index = weather.data.index  
    if only_with_event==True:
        if not 'septo_risk_with_event' in weather.data.columns:
            weather.check(varnames=['septo_risk_with_event'],
                          models={'septo_risk_with_event':add_septoria_risk_with_event}, 
                          wetness_duration_min = wetness_duration_min, 
                          temp_min = temp_min, temp_max = temp_max)
        ax.vlines(index, [0], weather.data.septo_risk_with_event, 'k', alpha=0.5)
    else:
        ax.vlines(index, [0], weather.data.septo_infection_risk, 'k', alpha=0.5)
        
    if leaf_dates is not None:
        col_name = leaf_dates.columns[0]
        for lf, row in leaf_dates.iterrows():
            ax.annotate('', xy=(row[col_name], 0.7), xycoords='data',
                        xytext=(row[col_name], 1.), 
                        arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='g'))
            
        if display_first_event == True:
            df = pd.DataFrame(np.zeros(len(weather.data)), index = index)
            for lf, row in leaf_dates.iterrows():
                ind_first = weather.data.degree_days[(weather.data['septo_risk_with_event']==1) & (weather.data['degree_days']>=row[col_name])].iloc[0]
                df.loc[ind_first,:] = 1
            ax.vlines(index, [0], df, 'r', linewidth = 1)
    
    ax.set_yticks([])
    ax.set_ylim([0,1])
    if title == None:
        ax.set_title(str(weather.data.index[0].year)+'-'+str(weather.data.index[-1].year))
    elif title != '':
        ax.set_title(title)
    if xlabel == True:
        ax.set_xlabel('Degree days', fontsize=16)
    if ylabel == True:
        ax.set_ylabel('Infection risk events')
    if xaxis == 'date':
        formatter = FuncFormatter(form_tick)
        ax.xaxis.set_major_formatter(FuncFormatter(formatter))
    if xlims!=None:
        ax.set_xlim(xlims)       
    return ax
    
def plot_rain_and_temp(weather, xaxis = 'degree_days', leaf_dates=None, 
                       ax = None, xlims=None, ylims_rain = None, ylims_temp = None,
                       title = None, xlabel = True, arrowstyle = '->', arrow_color = 'g'):
    """ Plot rain intensity (mm/h) and temperature (°Cd) """
    if ax == None:
        fig, ax = plt.subplots(figsize=(8,1.8))
        
    if xaxis == 'degree_days':
        index = weather.data.degree_days
    elif xaxis == 'date':
        index = weather.data.index  
    
    ax.bar(index, weather.data.rain, width = 1, color = 'b')
    ax2 = ax.twinx()
    ax2.plot(index, weather.data.temperature_air, 'r', alpha = 0.5)
    
    if leaf_dates is not None:
        col_name = leaf_dates.columns[0]
        for lf, row in leaf_dates.iterrows():
            ax.annotate('', xy=(row[col_name], 4.), xycoords='data',
                        xytext=(row[col_name], 5.), 
                        arrowprops=dict(arrowstyle=arrowstyle, 
                                        connectionstyle="arc3", 
                                        color=arrow_color))
    
    if title == None:
        ax.set_title(str(weather.data.index[0].year)+'-'+str(weather.data.index[-1].year))
    elif title != '':
        ax.set_title(title)
    if xlabel == True:
        ax.set_xlabel('Degree days', fontsize=16)
    if xaxis == 'date':
        formatter = FuncFormatter(form_tick)
        ax.xaxis.set_major_formatter(FuncFormatter(formatter))
    if xlims!=None:
        ax.set_xlim(xlims)
    if ylims_rain != None:
        ax.set_ylim(ylims_rain)
    if ylims_temp != None:
        ax2.set_ylim(ylims_temp)
        
    ax.set_ylabel('Rain (mm/h)')
    ax2.set_ylabel('Temperature (C)')
    return ax
    
def plot_relative_humidity(weather, xaxis = 'degree_days', 
                              ax = None, xlims=None, ylims = None, linestyle = '-', color = 'b',
                              title = None, xlabel = True):
    """ Plot relative humidity (%) """
    if ax == None:
        fig, ax = plt.subplots(figsize=(8,1.8))
        
    if xaxis == 'degree_days':
        index = weather.data.degree_days
    elif xaxis == 'date':
        index = weather.data.index  
    
    ax.plot(index, weather.data.relative_humidity, color = color, linestyle = linestyle)
    
    if title == None:
        ax.set_title(str(weather.data.index[0].year)+'-'+str(weather.data.index[-1].year))
    elif title != '':
        ax.set_title(title)
    if xlabel == True and xaxis == 'degree_days':
        ax.set_xlabel('Degree days', fontsize=16)
    if xaxis == 'date':
        formatter = FuncFormatter(form_tick)
        ax.xaxis.set_major_formatter(FuncFormatter(formatter))
    if xlims!=None:
        ax.set_xlim(xlims)
    if ylims != None:
        ax.set_ylim(ylims)
        
    ax.set_ylabel('Relative humidity (%)')
    return ax
    
def plot_daily_relative_humidity(weather, xaxis = 'degree_days', 
                                  ax = None, xlims=None, ylims = None, linestyle = '-', color = 'b',
                                  title = None, xlabel = True):
    
    df = weather.data.resample('D', how='mean')

    if ax == None:
        fig, ax = plt.subplots(figsize=(8,1.8))
        
    if xaxis == 'degree_days':
        index = df.degree_days
    elif xaxis == 'date':
        index = df.index  
    
    ax.plot(index, df.relative_humidity, color = color, linestyle = linestyle)
    
    if title == None:
        ax.set_title(str(df.index[0].year)+'-'+str(df.index[-1].year))
    elif title != '':
        ax.set_title(title)
    if xlabel == True and xaxis == 'degree_days':
        ax.set_xlabel('Degree days', fontsize=16)
    if xaxis == 'date':
        formatter = FuncFormatter(form_tick)
        ax.xaxis.set_major_formatter(FuncFormatter(formatter))
    if xlims!=None:
        ax.set_xlim(xlims)
    if ylims != None:
        ax.set_ylim(ylims)
        
    ax.set_ylabel('Relative humidity (%)')
    return ax
   
def plot_with_infection_risk(data, weather, leaves = range(1,2), 
                              variable='severity', xaxis = 'degree_days',
                              minimum_sample_size = 5, xlims = [1150, 2050], ylims = [0, 100], 
                              error_bars = False, return_ax = False, marker = 'd', with_brewer = True):
    """ Plot disease evolution with weather conditions underneath """
    weather_ = deepcopy(weather)
    try:
        import mpld3
        mpld3.disable_notebook()
    except:
        pass

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(3, 1, height_ratios=[5, 1, 1])

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])

    plot_by_leaf(data, weather_, leaves = leaves,  variable = variable, 
                 pointer = False, title_suffix = '_control', xaxis = xaxis, 
                 ax = ax1, xlims = xlims, ylims = ylims, with_brewer = with_brewer, xlabel = False, 
                 minimum_sample_size = minimum_sample_size, error_bars = error_bars, 
                 marker = marker)
    ax1.legend(ax1._get_legend_handles(), leaves, title = 'Leaf number',
                loc='center left', bbox_to_anchor=(1, 0.5))

    if xaxis in ['age_leaf', 'age_leaf_lig']:
        df_dates = get_df_dates_xaxis(data, xaxis)
        assert len(leaves) == 1, "Can't plot with 'age_leaf' in xaxis if more than 1 leaf at a time"
        date = df_dates.loc[leaves[0], xaxis]
        weather_.data['degree_days'] -= date
        df_dates -= date
    elif xaxis == 'age_leaf_vs_flag_lig':
        df_dates = get_df_dates_xaxis(data, xaxis)
        date = df_dates.loc[leaves[0], xaxis]
        weather_.data['degree_days'] -= date
        df_dates = get_df_dates_xaxis(data, 'age_leaf_lig')
        df_dates -= date
    elif xaxis == 'age_leaf_vs_flag_emg':
        df_dates = get_df_dates_xaxis(data, xaxis)
        date = df_dates.loc[leaves[0], xaxis]
        weather_.data['degree_days'] -= date
        df_dates = get_df_dates_xaxis(data, 'age_leaf')
        df_dates -= date
    elif xaxis == 'degree_days':
        df_dates = get_df_dates_xaxis(data, 'age_leaf')
        
    plot_septo_infection_risk(weather = weather_, leaf_dates = df_dates, ax = ax2,
                              xlims = xlims, title = '', xlabel = False, ylabel = True)
    plot_rain_and_temp(weather = weather_, ax = ax3, 
                       ylims_rain = [0., 5.], ylims_temp = [0., 30.],
                       xlims = xlims, title = '')
    if return_ax == True:
        return [ax1, ax2, ax3]

def plot_with_weather(data, weather, leaves = range(1,2), 
                      variable='severity', xaxis = 'degree_days',
                      minimum_sample_size = 5, xlims = [1150, 2050], ylims = [0, 100], 
                      error_bars = False, return_ax = False, marker = 'd', with_brewer = True):
    """ Plot disease evolution with weather conditions underneath """
    weather_ = deepcopy(weather)
    try:
        import mpld3
        mpld3.disable_notebook()
    except:
        pass

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[5, 2])

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    plot_by_leaf(data, weather_, leaves = leaves,  variable = variable, 
                 pointer = False, title_suffix = '_control', xaxis = xaxis, 
                 ax = ax1, xlims = xlims, ylims = ylims, with_brewer = with_brewer, xlabel = False, 
                 minimum_sample_size = minimum_sample_size, error_bars = error_bars, 
                 marker = marker)
    ax1.legend(ax1._get_legend_handles(), leaves, title = 'Leaf number',
                loc='center left', bbox_to_anchor=(1, 0.5))

    if xaxis in ['age_leaf', 'age_leaf_lig']:
        df_dates = get_df_dates_xaxis(data, xaxis)
        assert len(leaves) == 1, "Can't plot with 'age_leaf' in xaxis if more than 1 leaf at a time"
        date = df_dates.loc[leaves[0], xaxis]
        weather_.data['degree_days'] -= date
        df_dates -= date
    elif xaxis == 'age_leaf_vs_flag_lig':
        df_dates = get_df_dates_xaxis(data, xaxis)
        date = df_dates.loc[leaves[0], xaxis]
        weather_.data['degree_days'] -= date
        df_dates = get_df_dates_xaxis(data, 'age_leaf_lig')
        df_dates -= date
    elif xaxis == 'age_leaf_vs_flag_emg':
        df_dates = get_df_dates_xaxis(data, xaxis)
        date = df_dates.loc[leaves[0], xaxis]
        weather_.data['degree_days'] -= date
        df_dates = get_df_dates_xaxis(data, 'age_leaf')
        df_dates -= date
    elif xaxis == 'degree_days':
        df_dates = get_df_dates_xaxis(data, 'age_leaf')
            
    plot_rain_and_temp(weather = weather_, leaf_dates = df_dates, ax = ax2,
                       ylims_rain = [0., 5.], ylims_temp = [0., 30.],
                       xlims = xlims, title = '')
    if return_ax == True:
        return [ax1, ax2]

def plot_with_weather_by_fnl(data, weather, variable = 'severity',
                                      leaves = range(1, 5), xaxis = 'degree_days', 
                                      minimum_sample_size = 5, fig_size = (10, 10),
                                      xlims = None, error_bars = False, return_ax = False):
    """ Plot disease evolution with weather conditions underneath """
    data_grouped = group_by_fnl(data)
    try:
        import mpld3
        mpld3.disable_notebook()
    except:
        pass

    fig = plt.figure(figsize=fig_size)
    height_ratios = [5] + [1]
    gs = gridspec.GridSpec(2, 1, height_ratios=[5, 2])

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    
    plot_by_leaf_by_fnl(data, weather, leaves = leaves,  variable = variable, 
                        title_suffix = '', xaxis = xaxis, fig = fig,
                        ax = ax1, xlims = xlims, with_brewer = True, xlabel = False, 
                        minimum_sample_size = minimum_sample_size, error_bars = error_bars)
    
    arrow_styles = iter(['->', '-|>'])
    arrow_colors = iter(['g', 'k'])
    for fnl, df in data_grouped.iteritems():
        weather_ = deepcopy(weather)
        if xaxis in ['age_leaf', 'age_leaf_lig']:
            df_dates = get_df_dates_xaxis(df, xaxis)
            assert len(leaves) == 1, "Can't plot with 'age_leaf' in xaxis if more than 1 leaf at a time"
            date = df_dates.loc[leaves[0], xaxis]
            weather_.data['degree_days'] -= date
            df_dates -= date
        elif xaxis == 'age_leaf_vs_flag_lig':
            df_dates = get_df_dates_xaxis(df, xaxis)
            date = df_dates.loc[leaves[0], xaxis]
            weather_.data['degree_days'] -= date
            df_dates = get_df_dates_xaxis(df, 'age_leaf_lig')
            df_dates -= date
        elif xaxis == 'age_leaf_vs_flag_emg':
            df_dates = get_df_dates_xaxis(df, xaxis)
            date = df_dates.loc[leaves[0], xaxis]
            weather_.data['degree_days'] -= date
            df_dates = get_df_dates_xaxis(df, 'age_leaf')
            df_dates -= date
        elif xaxis == 'degree_days':
            df_dates = get_df_dates_xaxis(df, 'age_leaf')

        plot_rain_and_temp(weather = weather_, leaf_dates = df_dates, ax = ax2, 
                           ylims_rain = [0., 5.], ylims_temp = [0., 30.],
                           xlims = xlims, title = '', xlabel = True, 
                           arrowstyle = next(arrow_styles), arrow_color = next(arrow_colors))
    if return_ax == True:
        return axs
        
def plot_with_infection_risk_by_fnl(data, weather, variable = 'severity',
                                      leaves = range(1, 5), xaxis = 'degree_days', 
                                      minimum_sample_size = 5, fig_size = (10, 15),
                                      xlims = None, error_bars = False, return_ax = False):
    """ Plot disease evolution with weather conditions underneath """
    data_grouped = group_by_fnl(data)
    try:
        import mpld3
        mpld3.disable_notebook()
    except:
        pass

    fig = plt.figure(figsize=fig_size)
    height_ratios = [5] + [1 for i,k in enumerate(data_grouped.keys())] + [1]
    gs = gridspec.GridSpec(2+len(data_grouped.keys()), 1, height_ratios=height_ratios)

    axs = [fig.add_subplot(gs[i]) for i in range(len(gs.get_height_ratios()))]
    iter_axs = iter(axs)
    plot_by_leaf_by_fnl(data, weather, leaves = leaves,  variable = variable, 
                    title_suffix = '', xaxis = xaxis, fig = fig,
                    ax = iter_axs.next(), xlims = xlims, with_brewer = True, xlabel = False, 
                    minimum_sample_size = minimum_sample_size, error_bars = error_bars)
    
    for fnl, df in data_grouped.iteritems():
        weather_ = deepcopy(weather)
        if xaxis in ['age_leaf', 'age_leaf_lig']:
            df_dates = get_df_dates_xaxis(df, xaxis)
            assert len(leaves) == 1, "Can't plot with 'age_leaf' in xaxis if more than 1 leaf at a time"
            date = df_dates.loc[leaves[0], xaxis]
            weather_.data['degree_days'] -= date
            df_dates -= date
        elif xaxis == 'age_leaf_vs_flag_lig':
            df_dates = get_df_dates_xaxis(df, xaxis)
            date = df_dates.loc[leaves[0], xaxis]
            weather_.data['degree_days'] -= date
            df_dates = get_df_dates_xaxis(df, 'age_leaf_lig')
            df_dates -= date
        elif xaxis == 'age_leaf_vs_flag_emg':
            df_dates = get_df_dates_xaxis(df, xaxis)
            date = df_dates.loc[leaves[0], xaxis]
            weather_.data['degree_days'] -= date
            df_dates = get_df_dates_xaxis(df, 'age_leaf')
            df_dates -= date
        elif xaxis == 'degree_days':
            df_dates = get_df_dates_xaxis(df, 'age_leaf')
            
        plot_septo_infection_risk(weather = weather_, leaf_dates = df_dates, ax = iter_axs.next(),
                                  xlims = xlims, title = 'FNL %d' % fnl,
                                  xlabel = False, ylabel = True)
    
    plot_rain_and_temp(weather = weather_, ax = iter_axs.next(), 
                       ylims_rain = [0., 5.], ylims_temp = [0., 30.],
                       xlims = xlims, title = '', xlabel = True)
    if return_ax == True:
        return axs

def dates_to_julian(data, sowing_date="2010-10-15"):
    indexes = np.array([ind.strftime('%Y-%m-%d') for ind in data.index])
    ind_sowing = np.where(indexes==sowing_date)[0][0]
    data['Julian days'] = [(data.index[t]-data.index[ind_sowing]).days for t in range(len(data.index))]
    return data.set_index('Julian days')

def compare_degree_days(sowing_dates = ["2010-10-15", "2011-10-21", "2012-10-29"]):
    """ Compare evolution of degree days between years """
    sowing_dates = {int(sd[:4])+1:sd for sd in sowing_dates}
    data = {yr:[] for yr in sowing_dates.iterkeys()}
    for yr, date in sowing_dates.iteritems():
        weather = read_weather(*format_date(date, str(yr)+"-08-01"))
        df = pd.DataFrame(weather.data.degree_days).ix[::24]
        df.columns = [str(yr-1)+'-'+str(yr)]
        data[yr] = dates_to_julian(df, sowing_dates[yr])
    
    iterval = data.itervalues()
    df = next(iterval)
    while True:
        try:
            df = df.join(next(iterval))
        except StopIteration:
            break
    
    fig, ax = plt.subplots(figsize=(8,6))
    df.plot(ax=ax, linewidth=1.5)
    ax.set_ylabel('Degree days from sowing date', fontsize=18)
    ax.set_xlabel('Julian days from sowing date', fontsize=18)
    
# Comparison functions for two varieties in 2011 ###################################################
def compare_contingency_2011(data_mercia, data_rht3, weather):
    dict_df_count = {}
    df_count_mercia = table_count_notations(data_mercia, weather, variable = 'septo_green')
    df_count_mercia.columns = pd.MultiIndex.from_product([['Mercia'], df_count_mercia.columns])
    dict_df_count['mercia'] = df_count_mercia
    df_count_rht3 = table_count_notations(data_rht3, weather, variable = 'septo_green')
    df_count_rht3.columns = pd.MultiIndex.from_product([['Rht3'], df_count_rht3.columns])
    dict_df_count['rht3'] = df_count_rht3    
    df_empty = pd.DataFrame(['/' for i in range(len(df_count_mercia))], 
                            index = df_count_mercia.index, 
                            columns=pd.MultiIndex.from_product([['/'], ['/']]))
    return pd.concat([pd.concat([df, df_empty], axis=1) for df in dict_df_count.itervalues()],
                      axis = 1)

def get_comparison_table_2011(data_mercia, data_rht3, num_leaves = range(1,5), by_leaf = False):
    """ Create a comparison table by testing significative differences between vaireties """
    mercia_2011 = data_mercia.copy()
    rht3_2011 = data_rht3.copy()
    change_index(mercia_2011, new_index = ['datetime', 'num_leaf_top'])
    change_index(rht3_2011, new_index = ['datetime', 'num_leaf_top'])
    
    def one_step_comparison(leaves = range(1,5)):
        m_2011 = mercia_2011[[leaf in leaves for date, leaf in mercia_2011.index]]
        r_2011 = rht3_2011[[leaf in leaves for date, leaf in rht3_2011.index]]
        df = pd.DataFrame(index = m_2011.index.levels[0], columns = ['septo_green'])
        for date_ind in m_2011.index.levels[0]:
            treatment1 = m_2011.ix[date_ind][['septo_green']]
            treatment1.reset_index(drop=True, inplace=True)
            treatment1 = treatment1.dropna()

            treatment2 = r_2011.ix[date_ind][['septo_green']]
            treatment2.reset_index(drop=True, inplace=True)
            treatment2 = treatment2.dropna()

            # Do treatments differ ? (non parametric test)
            z_stat, p_val = st.ranksums(treatment1, treatment2)
            df.ix[date_ind,:] = True if p_val <= 0.05 else False   
        return df

    if by_leaf == True:
        comp_by_leaf = pd.concat([one_step_comparison([lf]) for lf in num_leaves], axis = 1)
        comp_by_leaf.columns = num_leaves
        return comp_by_leaf
    else:
        return one_step_comparison(num_leaves)
        
def save_comparison_table_2011(data_mercia, data_rht3, weather, num_leaves = range(1,5)):
    def add_column_level(df, column_name='Mercia'):
        df.columns = pd.MultiIndex.from_product([[column_name], df.columns])
        return df
    
    def get_dfs(data, num_leaves):
        data_ = data.copy()
        change_index(data_, ['datetime', 'num_leaf_top'])
        df_sev = pd.DataFrame(data_['septo_green'])
        df_mean = get_mean(data_, column='septo_green', by_leaf=True)
        df_low, df_high = get_confidence_interval(data_, weather, column='septo_green', by_leaf=True)
        df_mean, df_high, df_low = map(lambda x: x.loc[:,num_leaves], (df_mean, df_high, df_low))
        return map(lambda x: x[~np.isnan(x.sum(axis = 1))], (df_mean, df_high, df_low))
    
    mercia_2011 = data_mercia.copy()
    rht3_2011 = data_rht3.copy()
    df_mean_m, df_high_m, df_low_m = get_dfs(mercia_2011, num_leaves)
    df_mean_r, df_high_r, df_low_r = get_dfs(rht3_2011, num_leaves)
    stat_comp = get_comparison_table_2011(mercia_2011, rht3_2011, 
                                            num_leaves = num_leaves, by_leaf = True)
    (stat_comp, df_mean_m, df_high_m, df_low_m, df_mean_r, df_high_r, df_low_r) = map(lambda x: add_index_ddays(x, weather), 
                                                                                      (stat_comp, df_mean_m, df_high_m, df_low_m,
                                                                                       df_mean_r, df_high_r, df_low_r))

    df_conf_m = df_mean_m - df_low_m
    df_conf_r = df_mean_r - df_low_r
    (df_mean_m, df_conf_m) = map(lambda x: add_column_level(x, 'Mercia'), (df_mean_m, df_conf_m))
    (df_mean_r, df_conf_r) = map(lambda x: add_column_level(x, 'Rht3'),(df_mean_r, df_conf_r))

    data = pd.concat([df_mean_m, df_mean_r], axis = 1)
    df_conf = pd.concat([df_conf_m, df_conf_r], axis = 1)

    nb_leaves = len(num_leaves)
    nb_index = len(stat_comp.index.names)

    with open('figures/comparison_2011/comparison_table.tex', "w") as f:
        txt = ('\documentclass[a4paper,reqno,11pt]{amsart}' + '\n' + '\n' +
                '\usepackage[landscape]{geometry}' + '\n' +
                '\usepackage[french]{babel}' + '\n' +
                '\usepackage{booktabs}' + '\n' +
                '\\title{septoria on green}' + '\n'+ '\n' + 
                '\\begin{document}' + '\n'+ '\n' + 
                '\\begin{tabular}{'+ ''.join('l'*nb_index) + ''.join('r'*nb_leaves) + 'c' + ''.join('r'*nb_leaves)+ '}\n' +
                '\\toprule' + '\n' +
                '&  \multicolumn{1}{c}{} & \multicolumn{' + str(nb_leaves) + '}{c}{'+ data.columns.levels[0][0] + '} & \phantom{abc}'+
                ' & \multicolumn{' + str(nb_leaves) + '}{c}{'+ data.columns.levels[0][1] + '}' + '\\\\' + '\n' +
                '\cmidrule{'+ str(nb_index+1) + '-' + str(nb_index+nb_leaves)+ '}' + '\n' +
                '\cmidrule{'+ str(nb_index+nb_leaves+2) + '-' + str(nb_index+1+2*nb_leaves)+ '}' + '\n' +
                data.index.names[0] + ' & ' + data.index.names[1]+ ' & ' + ' & '.join(str(col) for col in data.columns.levels[1]) + 
                ' & '+ ' & ' + ' & '.join(str(col) for col in data.columns.levels[1]) + '\\\\' + '\n' +
                '\midrule' + '\n')
        for i, row in data.iterrows():
            txt += str(i[0])[:10] + ' & ' + str(i[1])
            for j, x in enumerate(row.values[:nb_leaves]):
                txt += ' &' + ''.join(' $\\ast$'*stat_comp.ix[i,j+1]) + '%.2f' % x + '$\pm$' + '%.2f' % df_conf.ix[i,j] if x!='-' else x
            txt += ' & '
            for j, x in enumerate(row.values[nb_leaves:]):
                txt += ' &' + ''.join(' $\\ast$'*stat_comp.ix[i,j+1]) + '%.2f' % x + '$\pm$' + '%.2f' % df_conf.ix[i,nb_leaves+j] if x!='-' else x
            txt += '\\\\' + '\n'

        txt+='\\bottomrule' + '\n' + '\end{tabular}' + '\n' + '\n' +  '\end{document}'
        f.write(txt)
    
    df_empty = pd.DataFrame(['/' for i in range(len(df_mean_m))], index = df_mean_m.index, 
                            columns=pd.MultiIndex.from_product([['/'], ['/']]))
    return pd.concat([df_mean_m, df_empty, df_mean_r], axis = 1)
    
def compare_confidence_and_boxplot_mercia_rht3(data_mercia, data_rht3, weather, xaxis = 'degree_days', leaves=range(1,5),
                                               fig_size=(15, 10), minimum_sample_size = 5., xlims = None):
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data_mercia, weather, leaves = leaves, variable = 'septo_green',
                                                                        xaxis = xaxis, title = False, fig_size = fig_size,
                                                                        xlims = xlims, return_fig = True, 
                                                                        marker = 'o', fixed_color = 'g', alpha = 0.2,
                                                                        linestyle = '-', filling_contour = False,
                                                                        minimum_sample_size = minimum_sample_size)
    labels = ['Mercia']
    proxy = [plt.Rectangle((0,0), 0,0, facecolor='g', alpha = 0.2)]
        
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data_rht3, weather, leaves = leaves, variable = 'septo_green',
                                                                        xaxis = xaxis, title = False, fig_size = fig_size,
                                                                        xlims = xlims, return_fig = True, fig = fig, axs = axs,
                                                                        marker = 'o', fixed_color = 'r', alpha = 0.2,
                                                                        linestyle = '-', filling_contour = False, 
                                                                        minimum_sample_size = minimum_sample_size)
    labels += ['Rht3']
    proxy += [plt.Rectangle((0,0), 0,0, facecolor='r', alpha = 0.2)]   
        
    axs[1][1].legend(proxy, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    
def compare_plot_with_infection_risk_mercia_rht3(data_mercia, data_rht3, weather, 
                                          leaves = range(1, 5), xaxis = 'degree_days', 
                                          minimum_sample_size = 5, xlims = None, error_bars = False):

    axs = plot_with_infection_risk(data_mercia, weather, leaves = leaves, variable = 'septo_green', 
                                    xaxis = xaxis, xlims = xlims, error_bars = error_bars,
                                    return_ax = True)
    
    weather_ = deepcopy(weather)
    
    plot_by_leaf(data_rht3, weather, leaves = leaves,  variable = 'septo_green', 
             pointer = False, title = '', title_suffix = '', xaxis = xaxis, 
             ax = axs[0], xlims = xlims, with_brewer = True, xlabel = False, 
             minimum_sample_size = minimum_sample_size, error_bars = error_bars, marker = 'o', linestyle = '--')
    
    # Customize legend
    ax = axs[0]
    labels = ['L '+str(lf)+': Mercia' for lf in leaves] + ['L '+str(lf)+': Rht3' for lf in leaves]
    color_list = [next(ax._get_lines.color_cycle) for col in range(7)]
    leg = ax.legend(ax._get_legend_handles(), labels, loc='center left', bbox_to_anchor=(1, 0.5))
    labels, handles = zip(*sorted(zip(labels, ax._get_legend_handles()), key=lambda t: color_list.index(t[1].get_c())))
    leg = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    
# Comparison functions for 2012 vs. 2013 analysis ##################################################
def compare_confidence_and_boxplot_2012_2013(data_12, data_13, weather_12, weather_13, 
                                             xaxis = 'degree_days', variable = 'severity', 
                                             leaves=range(1,5), fig_size=(15, 10), 
                                             minimum_sample_size = 5., xlims = None):
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data_12, weather_12, leaves = leaves, variable = variable,
                                                                    xaxis = xaxis, title_suffix='_control', fig_size = fig_size,
                                                                    xlims = xlims, return_fig = True, 
                                                                    marker = 'o', fixed_color = 'g', alpha = 0.2,
                                                                    linestyle = '-', filling_contour = False)
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data_13, weather_13, leaves = leaves, variable = variable,
                                                                        xaxis = xaxis, title_suffix='_control', fig_size = fig_size,
                                                                        xlims = xlims, return_fig = True, fig = fig, axs = axs,
                                                                        marker = 'o', fixed_color = 'r', alpha = 0.2,
                                                                        linestyle = '-', filling_contour = False)
    labels = ['2012', '2013']
    proxy = [plt.Rectangle((0,0), 0,0, facecolor='g', alpha = 0.2), plt.Rectangle((0,0), 0,0, facecolor='r', alpha = 0.2)]

    axs[1][1].legend(proxy, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    
def compare_plot_with_infection_risk_2012_2013(data_12, data_13, weather_12, weather_13,
                                                variable = 'severity', leaves = range(1, 5),
                                                xaxis = 'degree_days', minimum_sample_size = 5,
                                                xlims = None, fig_size = (10,15),
                                                error_bars = False, return_ax = False):
    weather_12_ = deepcopy(weather_12)
    weather_13_ = deepcopy(weather_13)
    try:
        import mpld3
        mpld3.disable_notebook()
    except:
        pass
    
    fig = plt.figure(figsize=fig_size)
    gs = gridspec.GridSpec(5, 1, height_ratios=[5, 1, 1, 1, 1])

    axs = [fig.add_subplot(gs[i]) for i,_ in enumerate(gs)]
    iter_axs = iter(axs)
    
    ax = iter_axs.next()
    plot_by_leaf(data_12, weather_12_, leaves = leaves,  variable = variable, 
                 pointer = False, title_suffix = '_control', xaxis = xaxis, 
                 ax = ax, xlims = xlims, with_brewer = True, xlabel = False, 
                 error_bars = error_bars, minimum_sample_size = minimum_sample_size)
    plot_by_leaf(data_13, weather_13_, leaves = leaves,  variable = variable, 
                 pointer = False, title_suffix = '_control', xaxis = xaxis, 
                 ax = ax, xlims = xlims, with_brewer = True, xlabel = False, 
                 error_bars = error_bars, minimum_sample_size = minimum_sample_size, 
                 marker = 'o', linestyle = '--')
    # Customize legend
    ax = axs[0]
    labels = ['L '+str(lf)+': Mercia' for lf in leaves] + ['L '+str(lf)+': Rht3' for lf in leaves]
    color_list = [next(ax._get_lines.color_cycle) for col in range(7)]
    leg = ax.legend(ax._get_legend_handles(), labels, loc='center left', bbox_to_anchor=(1, 0.5))
    labels, handles = zip(*sorted(zip(labels, ax._get_legend_handles()), key=lambda t: color_list.index(t[1].get_c())))
    leg = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))

    # Other plots    
    for weather_, df in zip([weather_12_, weather_13_], [data_12, data_13]):
        if xaxis in ['age_leaf', 'age_leaf_lig']:
            df_dates = get_df_dates_xaxis(df, xaxis)
            assert len(leaves) == 1, "Can't plot with 'age_leaf' in xaxis if more than 1 leaf at a time"
            date = df_dates.loc[leaves[0], xaxis]
            weather_.data['degree_days'] -= date
            df_dates -= date
        elif xaxis == 'age_leaf_vs_flag_lig':
            df_dates = get_df_dates_xaxis(df, xaxis)
            date = df_dates.loc[leaves[0], xaxis]
            weather_.data['degree_days'] -= date
            df_dates = get_df_dates_xaxis(df, 'age_leaf_lig')
            df_dates -= date
        elif xaxis == 'age_leaf_vs_flag_emg':
            df_dates = get_df_dates_xaxis(df, xaxis)
            date = df_dates.loc[leaves[0], xaxis]
            weather_.data['degree_days'] -= date
            df_dates = get_df_dates_xaxis(df, 'age_leaf')
            df_dates -= date
        elif xaxis == 'degree_days':
            df_dates = get_df_dates_xaxis(df, 'age_leaf')
            
        plot_septo_infection_risk(weather = weather_, leaf_dates = df_dates, ax = iter_axs.next(),
                                  xlims = xlims, title = weather_.data.index[-1].year,
                                  xlabel = False, ylabel = True)
    
        plot_rain_and_temp(weather = weather_, ax = iter_axs.next(), 
                           ylims_rain = [0., 5.], ylims_temp = [0., 30.],
                           xlims = xlims, title = '', xlabel = True)
                       
    if return_ax == True:
        return axs
        
# Comparison functions for reference treatment vs. control #########################################
def compare_confidence_and_boxplot_ref_ctrl(data_ctrl, data_ref, weather, xaxis = 'degree_days', variable = 'severity', 
                                            leaves=range(1,5), fig_size=(15, 10), minimum_sample_size = 5., xlims = None):
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data_ctrl, weather, leaves = leaves, variable = variable,
                                                                    xaxis = xaxis, title_suffix='_control', fig_size = fig_size,
                                                                    xlims = xlims, return_fig = True, 
                                                                    marker = 'o', fixed_color = 'r', alpha = 0.2,
                                                                    linestyle = '-', filling_contour = False)
    df_mean, df_high, df_low, fig, axs = plot_confidence_and_boxplot(data_ref, weather, leaves = leaves, variable = variable,
                                                                        xaxis = xaxis, title_suffix='_control', fig_size = fig_size,
                                                                        xlims = xlims, return_fig = True, fig = fig, axs = axs,
                                                                        marker = 'o', fixed_color = 'g', alpha = 0.2,
                                                                        linestyle = '-', filling_contour = False)
    labels = ['Control', 'Reference']
    proxy = [plt.Rectangle((0,0), 0,0, facecolor='r', alpha = 0.2), plt.Rectangle((0,0), 0,0, facecolor='g', alpha = 0.2)]

    axs[1][1].legend(proxy, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    
def compare_plot_with_infection_risk_ref_ctrl(data_ctrl, data_ref, weather, 
                                               leaves = range(1, 5), xaxis = 'degree_days', 
                                               minimum_sample_size = 5, xlims = None,
                                               error_bars = False):

    axs = plot_with_infection_risk(data_ctrl, weather, leaves = leaves, variable = 'severity',
                                    xaxis = xaxis, xlims = xlims, 
                                    minimum_sample_size = minimum_sample_size,
                                    error_bars = error_bars, return_ax = True)
    
    plot_by_leaf(data_ref, weather, leaves = leaves,  variable = 'severity', 
             pointer = False, title = '', title_suffix = '', xaxis = xaxis, 
             ax = axs[0], xlims = xlims, with_brewer = True, xlabel = False, 
             minimum_sample_size = minimum_sample_size, error_bars = error_bars, marker = 'o', linestyle = '--')
    
    # Customize legend
    ax = axs[0]
    labels = ['L '+str(lf)+': Control' for lf in leaves] + ['L '+str(lf)+': Reference' for lf in leaves]
    color_list = [next(ax._get_lines.color_cycle) for col in range(7)]
    leg = ax.legend(ax._get_legend_handles(), labels, loc='center left', bbox_to_anchor=(1, 0.5))
    labels, handles = zip(*sorted(zip(labels, ax._get_legend_handles()), key=lambda t: color_list.index(t[1].get_c())))
    leg = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
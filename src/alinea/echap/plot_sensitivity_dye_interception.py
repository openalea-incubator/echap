import alinea.echap.architectural_data as archidb
import pandas
import numpy


import alinea.echap.interception_data as idata


import matplotlib.pyplot as plt
plt.ion()

run_plot=False
# __________________ deprecated code from steph (to be checked !!

# sensibilite aux nombres de plantes et de simulation
# ici : choix de faire 5 simulations a 30, 60, 100 et 200 plantes
def sensibilite_nplants(var_lst=['Mercia', 'Rht3', 'Tremie12', 'Tremie13'],
                        axis='MS', csv=True):  # axis = MS or all
    df_all = pandas.DataFrame()
    for var in var_lst:
        df_sim = pandas.DataFrame()
        for stade in ['T1', 'T2']:
            if var == 'Mercia':
                nplants_lst = [43, 94, 149, 279]  # soit n=31, 61, 93, 205
            elif var == 'Rht3':
                nplants_lst = [49, 100, 150, 285]  # soit n=31, 58, 107, 201
            elif var == 'Tremie12':
                nplants_lst = [35]  # seulement n=25 pour Tremie12 car bug
            elif var == 'Tremie13':
                nplants_lst = [42, 82, 130, 235]  # soit n=28, 60, 100, 200
            for n in nplants_lst:
                x = 1
                while x <= 5:
                    npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n,
                                                 axis=axis, to_csv=False)
                    print('var = ' + var + ' stade = ' + stade + ' - nplants = ' + str(
                        npl) + ' - simulation num = ' + str(
                        x) + '/5 ok')  # verif etat avancement
                    dfmoy['var'] = var
                    dfmoy['dim'] = stade
                    dfmoy['nb_plantes_sim'] = npl
                    dfmoy['numero_sim'] = str(x)
                    df_sim = df_sim.append(dfmoy)
                    x += 1
        df_sim_gr = df_sim.groupby(
            ['var', 'HS', 'dim', 'nb_plantes_sim', 'numero_sim',
             'ntop_cur']).mean()
        df_sim_gr = df_sim_gr.reset_index()
        if csv == True:
            df_sim.to_csv('sensibilite_' + axis + '_all_' + var + '.csv')
            df_sim_gr.to_csv('sensibilite_' + axis + '_mean_' + var + '.csv')
        df_all = df_all.append(df_sim_gr)

    df_all = df_all.sort(['var', 'HS', 'ntop_cur'])
    df_all.to_csv('analyse_' + axis + '.csv', index=False)
    return df_all


# simulation interception des talles et des MB (4 valeurs des simulations cote a cote pour les maquettes)
def plot_diff(varieties=['Mercia']):
    for var in varieties:
        file = simulation_diff(var_lst=[var], axis='all')
        file.to_csv('plot_diff_' + var + '.csv', index=False)
        file = file[file['area'] > 0]
        file.to_csv('plot_diff_withoutArea0_' + var + '.csv', index=False)
        df_MS = file[file['axe'] == 'MS']
        df_MS = df_MS.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).mean()
        # df_MS_count = df_MS.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).count()
        # df_MS = df_MS.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).sum()
        df_MS = df_MS.reset_index();
        df_MS = df_MS[df_MS['ntop_cur'] <= 5]
        df_all = file
        df_all = df_all.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).mean()
        # df_all_count = df_all.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).count()
        # df_all = df_all.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).sum()
        df_all = df_all.reset_index();
        df_all = df_all[df_all['ntop_cur'] <= 5]
        df_talles = file[file['axe'] != 'MS']
        df_talles = df_talles.groupby(
            ['HS', 'nb_plantes_sim', 'ntop_cur']).mean()
        # df_talles_count = df_talles.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).count()
        # df_talles = df_talles.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).sum()
        df_talles = df_talles.reset_index();
        df_talles = df_talles[df_talles['ntop_cur'] <= 5]
        df_T1 = file[file['axe'] == 'T1']
        df_T1 = df_T1.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).mean()
        # df_T1_count = df_T1.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).count()
        # df_T1 = df_T1.groupby(['HS', 'nb_plantes_sim', 'ntop_cur']).sum()
        df_T1 = df_T1.reset_index();
        df_T1 = df_T1[df_T1['ntop_cur'] <= 5]

        bar_width = 0.2;
        opacity = 0.4
        if var == 'Mercia':
            val = [[0, 9.74], [1, 12.8]]
        elif var == 'Rht3':
            val = [[0, 9.15], [1, 12.48]]
        elif var == 'Tremie12':
            val = [[0, 10.98], [1, 12.63]]
        else:
            val = [[0, 8.7], [1, 11.04]]
        fig, axes = plt.subplots(nrows=1, ncols=2)
        for x, date in val:
            df_MS_date = df_MS[df_MS['HS'] == date]
            df_all_date = df_all[df_all['HS'] == date]
            df_talles_date = df_talles[df_talles['HS'] == date]
            df_T1_date = df_T1[df_T1['HS'] == date]

            n_groups = len(df_MS_date)
            index = numpy.arange(n_groups)

            rects1 = axes[x].bar(index, df_MS_date['deposit_Tartrazine'],
                                 bar_width, alpha=opacity, color=['c'])
            rects2 = axes[x].bar(index + bar_width,
                                 df_all_date['deposit_Tartrazine'], bar_width,
                                 alpha=opacity, color=['m'])
            rects3 = axes[x].bar(index + 2 * bar_width,
                                 df_talles_date['deposit_Tartrazine'],
                                 bar_width, alpha=opacity, color=['g'])
            rects4 = axes[x].bar(index + 3 * bar_width,
                                 df_T1_date['deposit_Tartrazine'], bar_width,
                                 alpha=opacity, color=['y'])
            '''rects1 = axes[x].bar(index, df_MS_date['area'], bar_width, alpha=opacity, color=['c'])
            rects2 = axes[x].bar(index + bar_width, df_all_date['area'], bar_width, alpha=opacity, color=['m'])
            rects3 = axes[x].bar(index + 2*bar_width, df_talles_date['area'], bar_width, alpha=opacity, color=['g'])
            rects4 = axes[x].bar(index + 3*bar_width, df_T1_date['area'], bar_width, alpha=opacity, color=['y'])'''

            # Mise en forme
            axes[x].set_ylim(0, 6)
            axes[x].set_xlim(0, len(df_MS_date))
            axes[x].set_xticks(index + bar_width)
            df_all_date['label'] = df_all_date['ntop_cur'].astype(str)
            axes[x].set_xticklabels(df_all_date['label'].tolist(), rotation=90,
                                    fontsize='small')
            if x == 0:
                axes[x].set_xlabel('ntop')
            axes[x].text(0.4, 5.5, var + ' - HS = ' + str(date),
                         bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10},
                         fontsize=12)

            fig.suptitle('Simulation somme surfaces des talles et des MB',
                         fontsize=10)
            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        axes[x].legend((rects1[0], rects2[0], rects3[0], rects4[0]), (
        'Interception MB', 'Interception MB + talles', 'Interception talles',
        'Interception T1'), bbox_to_anchor=[1.10, 1.12], prop={'size': 14})


# effet dimension (effet dimension _T2/T1)
def simulation_dimension(name='Tremie12', n_sim=5, n_plt=200, ajust=True):
    if name == 'Tremie12':
        n = 30
    else:
        n = n_plt

    # nff moyen
    adel = Reconstructions.get_reconstruction(name=name, nplants=n_plt)
    df_phenT = adel.phenT()
    nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
    df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
    mean = df_phenT.mean()
    hs_moyen = mean['nff']
    # HS T2 - nff moyen
    if ajust == True:
        lst = [[name, str(hs_moyen)]]
    else:
        lst = [[name, 'T2']]

    df_sim = pandas.DataFrame()
    x = 1
    while x <= n_sim:
        for var, stade in lst:
            dim_lst = [[0.5, '/2'], [0.75, '/1.25'], [0.9, '/1.1'], [1, '1'],
                       [1.1, 'x1.1'], [1.25, 'x1.25'], [1.5, 'x1.5']]
            for dim, dim_label in dim_lst:
                npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n,
                                             axis='MS', dimension=dim,
                                             to_csv=False)
                print('var = ' + var + ' stade = ' + dim_label + ' - nbre de plantes sim = ' + str(
                    npl))

                dfmoy['var'] = var;
                dfmoy['dim'] = dim_label
                df_sim = df_sim.append(dfmoy)
        x += 1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'dim', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    return df_sim_gr


def plot_dimension(plot1=True, plot2=False, plot3=True,
                   varieties=['Tremie12', 'Tremie13'], n_sim=1, n_plt=30,
                   ajust=True):
    '''
    Plot 1 =  Barplot observation / simulation avec ntop_cur en abscisse et (dfmoy['deposit_Tartrazin']  versus obs 'moyenne' par genotype et par date interception
    Plot 2 = Barplot de controle 'area' versus scan data (avec la colone area de dfmoy)
    Plot 3 = Barplot observation/simulation de deposit/area et obs/area

    !!! On plotte seulement les 3 feuilles du haut !!!
    '''

    # base
    df_scan_tremie12 = archidb.treatment_scan('Tremie12');
    df_scan_tremie12['var'] = 'Tremie12'
    df_scan_tremie13 = archidb.treatment_scan('Tremie13');
    df_scan_tremie13['var'] = 'Tremie13'
    df_scan = df_scan_tremie12.append(df_scan_tremie13)
    df_obs = idata.dye_interception()
    df_sim = pandas.DataFrame()
    for n in varieties:
        df_sim_var = simulation_dimension(name=n, n_sim=n_sim, n_plt=n_plt,
                                          ajust=ajust)
        df_sim = df_sim.append(df_sim_var)
    df_sim.to_csv('dimension.csv')
    # df_sim = 'dim_tremie13.csv'
    # df_sim = pandas.read_csv(df_sim, decimal='.', sep=',')

    # PLOT 1
    if plot1 is True:
        for var in varieties:
            df_obs_var = df_obs[df_obs['name'] == var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var'] == var]
            # df_obs_var.HS=map(str,df_obs_var.HS)
            # df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(
                df_obs_var.ix[:, ['HS', 'ntop_cur', 'mean']], how='outer')
            df_all = df_all[df_all['ntop_cur'] <= 3]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            df_all['HS-nffmean'] = df_all['HS'] - hs_moyen

            # plot
            bar_width = 0.4;
            opacity = 0.4
            if ajust == False:
                if var == 'Mercia':
                    val = [[0, 9.74], [1, 12.8]]
                elif var == 'Rht3':
                    val = [[0, 9.15], [1, 12.48]]
                elif var == 'Tremie12':
                    # val = [[0,10.98],[1,12.63]]
                    val = [[0, 12.63]]
                else:
                    # val = [[0,8.7],[1,11.04]]
                    val = [[0, 11.04]]
            else:
                val = [[0, round(hs_moyen, 2)]]

            fig, ax0 = plt.subplots(nrows=1, ncols=1)
            for x, date in val:
                df_all.HS = list(map(float, df_all.HS))
                df_all['HS'] = numpy.round(df_all['HS'], decimals=2)
                df_fin = df_all[df_all['HS'] == date]
                # index_dim pour classer les dim dans l ordre souhaite
                df_fin['index_dim'] = df_fin['dim']
                df_fin.ix[df_fin.dim.isin(['/2']), 'index_dim'] = 1
                df_fin.ix[df_fin.dim.isin(['/1.25']), 'index_dim'] = 2
                df_fin.ix[df_fin.dim.isin(['/1.1']), 'index_dim'] = 3
                df_fin.ix[df_fin.dim.isin(['1']), 'index_dim'] = 4
                df_fin.ix[df_fin.dim.isin(['x1.1']), 'index_dim'] = 5
                df_fin.ix[df_fin.dim.isin(['x1.25']), 'index_dim'] = 6
                df_fin.ix[df_fin.dim.isin(['x2']), 'index_dim'] = 7
                # rects2 pour ne tracer les obs que quand dim = normale
                df_fin['mean_dim'] = df_fin['mean']
                df_fin.ix[~df_fin.dim.isin(['1']), 'mean_dim'] = 0
                # ---
                df_fin = df_fin.sort(['ntop_cur', 'HS', 'index_dim'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)

                rects1 = ax0.bar(index, df_fin['deposit_Tartrazine'], bar_width,
                                 alpha=opacity,
                                 color=['c', 'c', 'c', 'c', 'c', 'c', 'c', 'm',
                                        'm', 'm', 'm', 'm', 'm', 'm', 'g', 'g',
                                        'g', 'g', 'g', 'g', 'g'])
                rects2 = ax0.bar(index + bar_width, df_fin['mean_dim'],
                                 bar_width, alpha=opacity, color='y')
                # Mise en forme
                ax0.set_ylim(0, 8)
                ax0.set_xlim(0, len(df_fin))
                ax0.set_xticks(index + bar_width)
                df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                df_fin['label'] = df_fin['dim'].astype(str) + ' - ' + df_fin[
                    'metamer'].astype(str)
                ax0.set_xticklabels(df_fin['label'].tolist(), rotation=90,
                                    fontsize='small')
                ax0.set_xlabel('dim/ntop')
                ax0.text(0.4, 7.5, var + ' - HS = ' + str(date),
                         bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10},
                         fontsize=12)

                fig.suptitle(
                    'Sim [deposit_Tartrazine] / Obs [mean] par ntop_cur',
                    fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
                # ax0.legend((rects1[0], rects1[7], rects1[14], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})

    # PLOT 2
    if plot2 is True:
        for var in varieties:
            if var == 'Tremie12' or var == 'Tremie13':
                df_scan_var = df_scan[df_scan['var'] == var]
                df_sim_var = df_sim[df_sim['var'] == var]
                df_scan_var.HS = list(map(str, df_scan_var.HS))
                df_sim_var.HS = list(map(str, df_sim_var.HS))
                df_all = df_sim_var.merge(
                    df_scan_var.ix[:, ['HS', 'ntop_cur', 'Area A_bl']],
                    how='outer')
                df_all = df_all[df_all['ntop_cur'] <= 3]
                # plot
                bar_width = 0.4;
                opacity = 0.4
                fig, axes = plt.subplots(nrows=1, ncols=2)
                if var == 'Mercia':
                    val = [[0, 9.74], [1, 12.8]]
                elif var == 'Rht3':
                    val = [[0, 9.15], [1, 12.48]]
                elif var == 'Tremie12':
                    val = [[0, 10.98], [1, 12.63]]
                else:
                    val = [[0, 8.7], [1, 11.04]]
                for x, date in val:
                    df_all.HS = list(map(float, df_all.HS))
                    df_fin = df_all[df_all['HS'] == date]
                    # index_dim
                    df_fin['index_dim'] = df_fin['dim']
                    df_fin.ix[df_fin.dim.isin(['/2']), 'index_dim'] = 1
                    df_fin.ix[df_fin.dim.isin(['/1.5']), 'index_dim'] = 2
                    df_fin.ix[df_fin.dim.isin(['1']), 'index_dim'] = 3
                    df_fin.ix[df_fin.dim.isin(['x1.5']), 'index_dim'] = 4
                    df_fin.ix[df_fin.dim.isin(['x2']), 'index_dim'] = 5
                    # rects2
                    df_fin['Area A_bl_dim'] = df_fin['Area A_bl']
                    df_fin.ix[~df_fin.dim.isin(['1']), 'Area A_bl_dim'] = 0
                    # ---
                    df_fin = df_fin.sort(['ntop_cur', 'HS', 'index_dim'])
                    n_groups = len(df_fin)
                    index = numpy.arange(n_groups)
                    rects1 = axes[x].bar(index, df_fin['Area A_bl_dim'],
                                         bar_width, alpha=opacity, color='y')
                    rects2 = axes[x].bar(index + bar_width, df_fin['area'],
                                         bar_width, alpha=opacity,
                                         color=['c', 'c', 'c', 'c', 'c', 'm',
                                                'm', 'm', 'm', 'm', 'g', 'g',
                                                'g', 'g', 'g'])

                    # Mise en forme
                    axes[x].set_ylim(0, 50)
                    axes[x].set_xlim(0, len(df_fin))
                    axes[x].set_xticks(index + bar_width)
                    df_fin['metamer'] = numpy.round(df_fin['metamer'],
                                                    decimals=1)
                    df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                    df_fin['label'] = df_fin['dim'].astype(str) + ' - ' + \
                                      df_fin['metamer'].astype(str)
                    axes[x].set_xticklabels(df_fin['label'].tolist(),
                                            rotation=90, fontsize='small')
                    if x == 0:
                        axes[x].set_xlabel('dim/ntop');
                        axes[x].set_ylabel('area (cm2)')
                    axes[x].text(0.4, 48, var + ' - HS = ' + str(date),
                                 bbox={'facecolor': '#FCF8F8', 'alpha': 0.6,
                                       'pad': 10}, fontsize=12)

                    fig.suptitle('Surface scan/sim par ntop_cur', fontsize=10)
                    fig.subplots_adjust(left=0.1, right=0.9, top=0.9,
                                        bottom=0.2)
                axes[x].legend((rects1[0], rects2[0], rects2[5], rects2[10]),
                               ('Scan', 'Sim F1', 'Sim F2', 'Sim F3'),
                               bbox_to_anchor=[1.10, 1.12], prop={'size': 14})

    # PLOT 3
    if plot3 is True:
        for var in varieties:
            df_obs_var = df_obs[df_obs['name'] == var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var'] == var]
            # df_obs_var.HS=map(str,df_obs_var.HS)
            # df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(
                df_obs_var.ix[:, ['HS', 'ntop_cur', 'mean']], how='outer')
            df_all = df_all[df_all['ntop_cur'] <= 3]
            df_all['tartrazine/area'] = df_all['deposit_Tartrazine'] / df_all[
                'area']
            df_all['mean/area'] = df_all['mean'] / df_all['area']
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            # plot
            bar_width = 0.4;
            opacity = 0.4
            if ajust == False:
                if var == 'Mercia':
                    # val = [[0,9.74],[1,12.8]]
                    val = [[0, 12.8]]
                elif var == 'Rht3':
                    # val = [[0,9.15],[1,12.48]]
                    val = [[0, 12.48]]
                elif var == 'Tremie12':
                    # val = [[0,10.98],[1,12.63]]
                    val = [[0, 12.63]]
                else:
                    # val = [[0,8.7],[1,11.04]]
                    val = [[0, 11.04]]
            else:
                val = [[0, round(hs_moyen, 2)]]

            fig, ax0 = plt.subplots(nrows=1, ncols=1)
            for x, date in val:
                df_all.HS = list(map(float, df_all.HS))
                df_all['HS'] = numpy.round(df_all['HS'], decimals=2)
                df_fin = df_all[df_all['HS'] == date]
                # index_dim
                df_fin['index_dim'] = df_fin['dim']
                df_fin.ix[df_fin.dim.isin(['/2']), 'index_dim'] = 1
                df_fin.ix[df_fin.dim.isin(['/1.25']), 'index_dim'] = 2
                df_fin.ix[df_fin.dim.isin(['/1.1']), 'index_dim'] = 3
                df_fin.ix[df_fin.dim.isin(['1']), 'index_dim'] = 4
                df_fin.ix[df_fin.dim.isin(['x1.1']), 'index_dim'] = 5
                df_fin.ix[df_fin.dim.isin(['x1.25']), 'index_dim'] = 6
                df_fin.ix[df_fin.dim.isin(['x2']), 'index_dim'] = 7
                # rects2
                df_fin['mean_dim'] = df_fin['mean/area']
                df_fin.ix[~df_fin.dim.isin(['1']), 'mean_dim'] = 0
                # ---
                df_fin = df_fin.sort(['ntop_cur', 'HS', 'index_dim'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)

                rects1 = ax0.bar(index, df_fin['tartrazine/area'], bar_width,
                                 alpha=opacity,
                                 color=['c', 'c', 'c', 'c', 'c', 'c', 'c', 'm',
                                        'm', 'm', 'm', 'm', 'm', 'm', 'g', 'g',
                                        'g', 'g', 'g', 'g', 'g'])
                rects2 = ax0.bar(index + bar_width, df_fin['mean_dim'],
                                 bar_width, alpha=opacity, color='y')
                # Mise en forme
                ax0.set_ylim(0, 0.5)
                ax0.set_xlim(0, len(df_fin))
                ax0.set_xticks(index + bar_width)
                df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                df_fin['label'] = df_fin['dim'].astype(str) + ' - ' + df_fin[
                    'metamer'].astype(str)
                ax0.set_xticklabels(df_fin['label'].tolist(), rotation=90,
                                    fontsize='small')
                ax0.set_xlabel('dim/ntop')
                ax0.text(0.4, 0.46, var + ' - HS = ' + str(date),
                         bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10},
                         fontsize=12)

                fig.suptitle(
                    'Sim [deposit_Tartrazine/area] / Obs [mean/area] par ntop_cur',
                    fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
                # ax0.legend((rects1[0], rects1[7], rects1[14], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )

    return df_scan, df_obs, df_sim


def courbe_dimension(varieties=['Tremie12', 'Tremie13']):
    # courbe interception par surface a partir du fichier de sortie de plot_dimension
    df_sim = 'dimension.csv'
    df_sim = pandas.read_csv(df_sim, decimal='.', sep=',')

    for var in varieties:
        df_sim_var = df_sim[df_sim['var'] == var]
        for ntop in [1, 2, 3]:
            df_sim_var_ntop = df_sim_var[df_sim_var['ntop'] == ntop]
            if var == 'Tremie12':
                picto = 'o'
            elif var == 'Tremie13':
                picto = '^'
            if ntop == 1:
                color = 'c'
            elif ntop == 2:
                color = 'm'
            elif ntop == 3:
                color = 'g'
            if ntop == 1:
                plt.plot(df_sim_var_ntop['area'],
                         df_sim_var_ntop['deposit_Tartrazine'], picto + color,
                         label=var)
            else:
                plt.plot(df_sim_var_ntop['area'],
                         df_sim_var_ntop['deposit_Tartrazine'], picto + color,
                         label='_nolegend_')
        plt.ylim(0, 8);
        plt.xlim(0, 50)
        plt.xlabel('surface');
        plt.ylabel('interception')
        plt.legend(numpoints=1, bbox_to_anchor=(1., 1.), prop={'size': 9})


# effet densite  (efffet density T2- T1)
def simulation_density(name='Tremie12', n_sim=5, n_plt=200, ajust=True):
    if name == 'Tremie12':
        n = 30
    else:
        n = n_plt

    # nff moyen
    adel = Reconstructions.get_reconstruction(name=name, nplants=n_plt)
    df_phenT = adel.phenT()
    nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
    df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
    mean = df_phenT.mean()
    hs_moyen = mean['nff']

    if ajust == True:
        lst = [[name, str(hs_moyen)]]
    else:
        lst = [[name, 'T2']]

    df_sim = pandas.DataFrame()
    x = 1
    while x <= n_sim:
        for var, stade in lst:
            dens_lst = [0.1, 0.5, 0.75, 0.9, 1, 1.1, 1.25, 1.5, 1.9]
            for dens in dens_lst:
                npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n,
                                             axis='MS', to_csv=False,
                                             density=dens)
                print('var = ' + var + ' stade = ' + stade + ' density = ' + str(
                    dens) + ' - nbre de plantes sim = ' + str(npl))
                dfmoy['var'] = var
                dfmoy['density'] = dens
                df_sim = df_sim.append(dfmoy)
        x += 1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'density', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    return df_sim_gr


def plot_density(plot1=True, plot2=False, plot3=True, plot4=False,
                 varieties=['Tremie12', 'Tremie13'], n_sim=1, n_plt=30,
                 ajust=True):
    '''
    Plot 1 =  Barplot observation / simulation avec ntop_cur en abscisse et (dfmoy['deposit_Tartrazin']  versus obs 'moyenne' par genotype et par date interception
    Plot 2 = Barplot de controle 'area' versus scan data (avec la colone area de dfmoy)
    !!! plot 2 = PAS DE SCAN POUR MERCIA ET RHT3 donc pas de graph genere (test si plot2=True, ne posera pas de pb)
    Plot 3 = Barplot observation/simulation de deposit/area et obs/area
    Plot 4 = 2 graphes supplementaires recapitulatifs avec en x = age des feuilles (en dd depuis la ligulation, quelque soit le numero des feuilles) ET Ntop (2 graphes donc) et en y = interception feuille par cm2
    '''

    # base
    df_scan_tremie12 = archidb.treatment_scan('Tremie12');
    df_scan_tremie12['var'] = 'Tremie12'
    df_scan_tremie13 = archidb.treatment_scan('Tremie13');
    df_scan_tremie13['var'] = 'Tremie13'
    df_scan = df_scan_tremie12.append(df_scan_tremie13)
    df_obs = idata.dye_interception()
    df_sim = pandas.DataFrame()
    for n in varieties:
        df_sim_var = simulation_density(name=n, n_sim=n_sim, n_plt=n_plt,
                                        ajust=ajust)
        df_sim = df_sim.append(df_sim_var)
    df_sim.to_csv('density.csv')
    # test
    # df_sim = 'density_synth_tremie12.csv'
    # df_sim = pandas.read_csv(df_sim, decimal='.', sep=',')

    # PLOT 1
    if plot1 is True:
        for var in varieties:
            df_obs_var = df_obs[df_obs['name'] == var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var'] == var]
            # df_obs_var.HS=map(str,df_obs_var.HS)
            # df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(
                df_obs_var.ix[:, ['HS', 'ntop_cur', 'mean']], how='outer')
            df_all = df_all[df_all['ntop_cur'] <= 3]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            # plot
            bar_width = 0.4;
            opacity = 0.4
            if ajust == False:
                if var == 'Mercia':
                    # val = [[0,9.74],[1,12.8]]
                    val = [[0, 12.8]]
                elif var == 'Rht3':
                    # val = [[0,9.15],[1,12.48]]
                    val = [[0, 12.48]]
                elif var == 'Tremie12':
                    # val = [[0,10.98],[1,12.63]]
                    val = [[0, 12.63]]
                else:
                    # val = [[0,8.7],[1,11.04]]
                    val = [[0, 11.04]]
            else:
                val = [[0, round(hs_moyen, 2)]]

            fig, ax0 = plt.subplots(nrows=1, ncols=1)
            for x, date in val:
                df_all.HS = list(map(float, df_all.HS))
                df_all['HS'] = numpy.round(df_all['HS'], decimals=2)
                df_fin = df_all[df_all['HS'] == date]
                # mean_dens pour dessiner les obs seulement qd densite=1
                df_fin['mean_dens'] = df_fin['mean']
                df_fin.ix[~df_fin.density.isin([1]), 'mean_dens'] = 0
                # ---
                df_fin = df_fin.sort(['ntop_cur', 'HS', 'density'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)

                rects1 = ax0.bar(index, df_fin['deposit_Tartrazine'], bar_width,
                                 alpha=opacity,
                                 color=['c', 'c', 'c', 'c', 'c', 'c', 'c', 'c',
                                        'c', 'm', 'm', 'm', 'm', 'm', 'm', 'm',
                                        'm', 'm', 'g', 'g', 'g', 'g', 'g', 'g',
                                        'g', 'g', 'g', 'r', 'r', 'r', 'r', 'r',
                                        'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b',
                                        'b', 'b', 'b', 'b', 'b'])
                rects2 = ax0.bar(index + bar_width, df_fin['mean_dens'],
                                 bar_width, alpha=opacity, color='y')
                # Mise en forme
                ax0.set_ylim(0, 8)
                ax0.set_xlim(0, len(df_fin))
                ax0.set_xticks(index + bar_width)
                df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                df_fin['label'] = df_fin['density'].astype(str) + ' - ' + \
                                  df_fin['metamer'].astype(str)
                ax0.set_xticklabels(df_fin['label'].tolist(), rotation=90,
                                    fontsize='small')
                if x == 0:
                    ax0.set_xlabel('densite/ntop')
                ax0.text(0.4, 7.5, var + ' - HS = ' + str(date),
                         bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10},
                         fontsize=12)

                fig.suptitle(
                    'Sim [deposit_Tartrazine] / Obs [mean] par ntop_cur pour ' + var,
                    fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
                # ax0.legend((rects1[0], rects1[9], rects1[18], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})

    # PLOT 2 - !!! PAS DE SCAN POUR MERCIA ET RHT3 donc pas de graph genere
    if plot2 is True:
        for var in varieties:
            if var == 'Tremie12' or var == 'Tremie13':
                df_scan_var = df_scan[df_scan['var'] == var]
                df_sim_var = df_sim[df_sim['var'] == var]
                df_scan_var.HS = list(map(str, df_scan_var.HS))
                df_sim_var.HS = list(map(str, df_sim_var.HS))
                df_all = df_sim_var.merge(
                    df_scan_var.ix[:, ['HS', 'ntop_cur', 'Area A_bl']],
                    how='outer')
                df_all = df_all[df_all['ntop_cur'] <= 5]
                # plot
                bar_width = 0.4;
                opacity = 0.4
                fig, axes = plt.subplots(nrows=1, ncols=2)
                if var == 'Tremie12':
                    val = [[0, 10.98], [1, 12.63]]
                else:
                    val = [[0, 8.7], [1, 11.04]]
                for x, date in val:
                    df_all.HS = list(map(float, df_all.HS))
                    df_fin = df_all[df_all['HS'] == date]
                    # Area A_bl dens pour dessiner les obs seulement qd densite=1
                    df_fin['Area A_bl dens'] = df_fin['Area A_bl']
                    df_fin.ix[~df_fin.density.isin([1]), 'Area A_bl dens'] = 0
                    df_fin = df_fin.sort(['ntop_cur', 'HS', 'density'])
                    n_groups = len(df_fin)
                    index = numpy.arange(n_groups)

                    rects1 = axes[x].bar(index, df_fin['Area A_bl dens'],
                                         bar_width, alpha=opacity, color='y')
                    rects2 = axes[x].bar(index + bar_width, df_fin['area'],
                                         bar_width, alpha=opacity,
                                         color=['c', 'c', 'c', 'c', 'c', 'c',
                                                'c', 'm', 'm', 'm', 'm', 'm',
                                                'm', 'm', 'g', 'g', 'g', 'g',
                                                'g', 'g', 'g', 'r', 'r', 'r',
                                                'r', 'r', 'r', 'r', 'b', 'b',
                                                'b', 'b', 'b', 'b', 'b'])

                    # Mise en forme
                    axes[x].set_ylim(0, 30)
                    axes[x].set_xlim(0, len(df_fin))
                    axes[x].set_xticks(index + bar_width)
                    df_fin['metamer'] = numpy.round(df_fin['metamer'],
                                                    decimals=1)
                    df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                    df_fin['label'] = df_fin['density'].astype(str) + ' - ' + \
                                      df_fin['TT'].astype(str) + ' - ' + df_fin[
                                          'metamer'].astype(str)
                    axes[x].set_xticklabels(df_fin['label'].tolist(),
                                            rotation=90, fontsize='x-small')
                    if x == 0:
                        axes[x].set_xlabel('densite/TT/ntop');
                        axes[x].set_ylabel('area (cm2)')
                    axes[x].text(0.4, 28.5, var + ' - HS = ' + str(date),
                                 bbox={'facecolor': '#FCF8F8', 'alpha': 0.6,
                                       'pad': 10}, fontsize=12)
                    fig.suptitle('Surface scan/sim par ntop_cur pour ' + var,
                                 fontsize=10)
                    fig.subplots_adjust(left=0.1, right=0.9, top=0.9,
                                        bottom=0.2)
                axes[x].legend((rects1[0], rects2[0], rects2[7], rects2[14]),
                               ('Scan', 'Sim F1', 'Sim F2', 'Sim F3'),
                               bbox_to_anchor=[1.10, 1.12], prop={'size': 14})
            else:
                print('Pas de scan pour variete : ' + var)

    # PLOT 3
    if plot3 is True:
        for var in varieties:
            df_obs_var = df_obs[df_obs['name'] == var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var'] == var]
            # df_obs_var.HS=map(str,df_obs_var.HS)
            # df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(
                df_obs_var.ix[:, ['HS', 'ntop_cur', 'mean']], how='outer')
            df_all = df_all[df_all['ntop_cur'] <= 3]
            df_all['tartrazine/area'] = df_all['deposit_Tartrazine'] / df_all[
                'area']
            df_all['mean/area'] = df_all['mean'] / df_all['area']
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            # plot
            bar_width = 0.4;
            opacity = 0.4
            if ajust == False:
                if var == 'Mercia':
                    # val = [[0,9.74],[1,12.8]]
                    val = [[0, 12.8]]
                elif var == 'Rht3':
                    # val = [[0,9.15],[1,12.48]]
                    val = [[0, 12.48]]
                elif var == 'Tremie12':
                    # val = [[0,10.98],[1,12.63]]
                    val = [[0, 12.63]]
                else:
                    # val = [[0,8.7],[1,11.04]]
                    val = [[0, 11.04]]
            else:
                val = [[0, round(hs_moyen, 2)]]

            fig, ax0 = plt.subplots(nrows=1, ncols=1)
            for x, date in val:
                df_all.HS = list(map(float, df_all.HS))
                df_all['HS'] = numpy.round(df_all['HS'], decimals=2)
                df_fin = df_all[df_all['HS'] == date]
                # mean_dens pour dessiner les obs seulement qd densite=1
                df_fin['mean/area_dens'] = df_fin['mean/area']
                df_fin.ix[~df_fin.density.isin([1]), 'mean/area_dens'] = 0
                # ---
                df_fin = df_fin.sort(['ntop_cur', 'HS', 'density'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)

                rects1 = ax0.bar(index, df_fin['tartrazine/area'], bar_width,
                                 alpha=opacity,
                                 color=['c', 'c', 'c', 'c', 'c', 'c', 'c', 'c',
                                        'c', 'm', 'm', 'm', 'm', 'm', 'm', 'm',
                                        'm', 'm', 'g', 'g', 'g', 'g', 'g', 'g',
                                        'g', 'g', 'g', 'r', 'r', 'r', 'r', 'r',
                                        'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b',
                                        'b', 'b', 'b', 'b', 'b'])
                rects2 = ax0.bar(index + bar_width, df_fin['mean/area_dens'],
                                 bar_width, alpha=opacity, color='y')
                # Mise en forme
                ax0.set_ylim(0, 0.5)
                ax0.set_xlim(0, len(df_fin))
                ax0.set_xticks(index + bar_width)
                df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                df_fin['label'] = df_fin['density'].astype(str) + ' - ' + \
                                  df_fin['metamer'].astype(str)
                ax0.set_xticklabels(df_fin['label'].tolist(), rotation=90,
                                    fontsize='x-small')
                ax0.set_xlabel('densite/ntop')
                ax0.text(0.4, 0.46, var + ' - HS = ' + str(date),
                         bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10},
                         fontsize=12)
                fig.suptitle(
                    'Sim [deposit_Tartrazine/area] / Obs [mean/area] par ntop_cur pour ' + var,
                    fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
                # ax0.legend((rects1[0], rects1[9], rects1[18], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )

    # PLOT 4
    if plot4 is True:
        for var in varieties:
            df_obs_var = df_obs[df_obs['name'] == var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var'] == var]
            df_obs_var.HS = list(map(str, df_obs_var.HS))
            df_sim_var.HS = list(map(str, df_sim_var.HS))
            df_all = df_sim_var.merge(df_obs_var.ix[:, ['HS', 'ntop_cur', 'mean']], how='outer')
            df_all = df_all[df_all['ntop_cur'] <= 5]
            # plot
            bar_width = 0.4;
            opacity = 0.4
            if var == 'Mercia':
                val = [[0, 9.74], [1, 12.8]]
            elif var == 'Rht3':
                val = [[0, 9.15], [1, 12.48]]
            elif var == 'Tremie12':
                val = [[0, 10.98], [1, 12.63]]
            else:
                val = [[0, 8.7], [1, 11.04]]
            fig, axes = plt.subplots(nrows=1, ncols=2)
            for x, date in val:
                df_all.HS = list(map(float, df_all.HS))
                df_fin = df_all[df_all['HS'] == date]
                # mean_dens pour dessiner les obs seulement qd densite=1
                df_fin['mean_dens'] = df_fin['mean']
                df_fin.ix[~df_fin.density.isin([1]), 'mean_dens'] = 0
                # ---
                df_fin = df_fin.sort(['ntop_cur', 'leaf_emergence'],
                                     ascending=[1, 0])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)

                rects1 = axes[x].bar(index, df_fin['deposit_Tartrazine'],
                                     bar_width, alpha=opacity,
                                     color=['c', 'c', 'c', 'c', 'c', 'm', 'm',
                                            'm', 'm', 'm', 'g', 'g', 'g', 'g',
                                            'g', 'r', 'r', 'r', 'r', 'r', 'b',
                                            'b', 'b', 'b', 'b'])
                rects2 = axes[x].bar(index + bar_width, df_fin['mean_dens'],
                                     bar_width, alpha=opacity, color='y')
                # Mise en forme
                axes[x].set_ylim(0, 8)
                axes[x].set_xlim(0, len(df_fin))
                axes[x].set_xticks(index + bar_width)
                df_fin['leaf_emergence'] = numpy.round(df_fin['leaf_emergence'],
                                                       decimals=1)
                df_fin['label'] = df_fin['leaf_emergence'].astype(str) + ' - ' + \
                                  df_fin['density'].astype(str)
                axes[x].set_xticklabels(df_fin['label'].tolist(), rotation=90,
                                        fontsize='small')
                if x == 0:
                    axes[x].set_xlabel('leaf_since_emergence/densite')
                axes[x].text(0.4, 7.4, var + ' - HS = ' + str(date),
                             bbox={'facecolor': '#FCF8F8', 'alpha': 0.6,
                                   'pad': 10}, fontsize=12)

                fig.suptitle(
                    'Sim [deposit_Tartrazine] / Obs [mean] par ntop_cur pour ' + var,
                    fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            axes[x].legend((rects1[0], rects1[5], rects1[10], rects2[0]),
                           ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'),
                           bbox_to_anchor=[1.10, 1.12], prop={'size': 14})

    return df_scan, df_obs, df_sim


def courbe_density(varieties=['Tremie12', 'Tremie13']):
    # courbe interception par nbr epi/plante*densite imposee a partir du fichier de sortie de plot_density
    df_sim = 'density.csv'
    df_sim = pandas.read_csv(df_sim, decimal='.', sep=',')

    for var in varieties:
        df_sim_var = df_sim[df_sim['var'] == var]
        for ntop in [1, 2, 3]:
            df_sim_var_ntop = df_sim_var[df_sim_var['ntop'] == ntop]
            if var == 'Tremie12':
                picto = 'o'
                df_sim_var_ntop['earsm2xdens'] = 728 * df_sim_var_ntop[
                    'density']
            elif var == 'Tremie13':
                picto = '^'
                df_sim_var_ntop['earsm2xdens'] = 530 * df_sim_var_ntop['density']
            if ntop == 1:
                color = 'c'
            elif ntop == 2:
                color = 'm'
            elif ntop == 3:
                color = 'g'
            if ntop == 1:
                plt.plot(df_sim_var_ntop['earsm2xdens'],
                         df_sim_var_ntop['deposit_Tartrazine'], picto + color,
                         label=var)
            else:
                plt.plot(df_sim_var_ntop['earsm2xdens'],
                         df_sim_var_ntop['deposit_Tartrazine'], picto + color,
                         label='_nolegend_')
        plt.ylim(0, 8)
        plt.xlabel('nb epis/plante * densite');
        plt.ylabel('interception')
        plt.legend(numpoints=1, bbox_to_anchor=(1., 1.), prop={'size': 9})


# effet HS _T2 / T1
def simulation_HS(name='Tremie12', n_sim=5, n_plt=200, ajust=True):
    if name == 'Tremie12':
        n = 30
    else:
        n = n_plt

    # nff moyen
    adel = Reconstructions.get_reconstruction(name=name, nplants=n_plt)
    df_phenT = adel.phenT()
    nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
    df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
    mean = df_phenT.mean()
    hs_moyen = mean['nff']
    # HS T2 - nff moyen
    date, hsT2 = HS_applications[name]['T2']
    # date a simuler
    if ajust == True:
        lst = [[name, str(hs_moyen - 2)], [name, str(hs_moyen - 1.5)],
               [name, str(hs_moyen - 1)], [name, str(hs_moyen - 0.5)],
               [name, str(hs_moyen)], [name, str(hs_moyen + 0.5)],
               [name, str(hs_moyen + 1)], [name, str(hs_moyen + 1.5)],
               [name, str(hs_moyen + 2)], [name, str(hs_moyen + 2.5)],
               [name, str(hs_moyen + 3)], [name, str(hs_moyen + 3.5)],
               [name, str(hs_moyen + 4)], [name, str(hs_moyen + 4.5)],
               [name, str(hs_moyen + 5)]]
    else:
        # lst = [[name, str(hsT2-0.8)], [name, str(hsT2-0.6)], [name, str(hsT2-0.4)], [name, str(hsT2-0.2)], [name, str(hsT2)], [name, str(hsT2+0.2)], [name, str(hsT2+0.4)], [name, str(hsT2+0.6)], [name, str(hsT2+0.8)]]
        lst = [[name, str(hsT2 - 2)], [name, str(hsT2 - 1.5)],
               [name, str(hsT2 - 1)], [name, str(hsT2 - 0.5)],
               [name, str(hsT2)], [name, str(hsT2 + 0.5)],
               [name, str(hsT2 + 1)], [name, str(hsT2 + 1.5)],
               [name, str(hsT2 + 2)], [name, str(hsT2 + 2.5)],
               [name, str(hsT2 + 3)], [name, str(hsT2 + 3.5)],
               [name, str(hsT2 + 4)], [name, str(hsT2 + 4.5)],
               [name, str(hsT2 + 5)]]

    df_sim = pandas.DataFrame()
    x = 1
    while x <= n_sim:
        for var, stade in lst:
            npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n,
                                         axis='MS', to_csv=False)
            print('var = ' + var + ' stade = ' + stade + ' - nbre de plantes sim = ' + str(
                npl))
            dfmoy['var'] = var
            df_sim = df_sim.append(dfmoy)
        x += 1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    return df_sim_gr


def plot_HS(plot1=True, plot2=False, plot3=True,
            varieties=['Tremie12', 'Tremie13'], n_sim=1, n_plt=30, ajust=False):
    '''
    Plot 1 =  Barplot observation / simulation avec ntop_cur en abscisse et (dfmoy['deposit_Tartrazin']  versus obs 'moyenne' par genotype et par date interception
    Plot 2 = Barplot de controle 'area' versus scan data (avec la colone area de dfmoy)
    Plot 3 = Barplot observation/simulation de deposit/area et obs/area
    '''

    # base
    df_scan_tremie12 = archidb.treatment_scan('Tremie12');
    df_scan_tremie12['var'] = 'Tremie12'
    df_scan_tremie13 = archidb.treatment_scan('Tremie13');
    df_scan_tremie13['var'] = 'Tremie13'
    df_scan = df_scan_tremie12.append(df_scan_tremie13)
    df_obs = idata.dye_interception()
    df_sim = pandas.DataFrame()
    for n in varieties:
        df_sim_var = simulation_HS(name=n, n_sim=n_sim, n_plt=n_plt,
                                   ajust=ajust)
        df_sim = df_sim.append(df_sim_var)
    df_sim.to_csv('HS.csv')
    # test
    # df_sim = 'HS_synth_tremie13.csv'
    # df_sim = pandas.read_csv(df_sim, decimal='.', sep=',')

    # PLOT 1
    if plot1 is True:
        for var in varieties:
            df_obs_var = df_obs[df_obs['name'] == var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var'] == var]
            # df_obs_var.HS=map(str,df_obs_var.HS)
            # df_sim_var.HS=map(str,df_sim_var.HS)
            df_all = df_sim_var.merge(
                df_obs_var.ix[:, ['HS', 'ntop_cur', 'mean']], how='outer')
            df_all = df_all[df_all['ntop_cur'] <= 3]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            # plot
            bar_width = 0.4;
            opacity = 0.4
            if ajust == False:
                if var == 'Mercia':
                    # val = [[0,9.74],[1,12.8]]
                    val = [[0, 12.8]]
                elif var == 'Rht3':
                    # val = [[0,9.15],[1,12.48]]
                    val = [[0, 12.48]]
                elif var == 'Tremie12':
                    # val = [[0,10.98],[1,12.63]]
                    val = [[0, 12.63]]
                elif var == 'Tremie13':
                    # val = [[0,8.7],[1,11.04]]
                    val = [[0, 11.04]]
            else:
                val = [[0, round(hs_moyen, 2)]]
            fig, ax0 = plt.subplots(nrows=1, ncols=1)
            for x, date in val:
                df_all.HS = list(map(float, df_all.HS))
                df_all['HS'] = numpy.round(df_all['HS'], decimals=2)
                df_fin = df_all[df_all['HS'] == date]
                # for f in [-0.8,-0.6,-0.4,-0.2,+0.2,+0.4,+0.6,+0.8]:
                for f in [-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4,
                          4.5, 5]:
                    df_fin = df_fin.append(
                        df_all[df_all['HS'] == round(date + f, 2)])
                df_fin = df_fin.sort(['ntop_cur', 'HS'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)

                if date == 10.98 or date == 8.7 or date == 9.15 or date == 9.74:
                    rects1 = ax0.bar(index, df_fin['deposit_Tartrazine'],
                                     bar_width, alpha=opacity,
                                     color=['c', 'c', 'c', 'c', 'c', 'm', 'm',
                                            'm', 'm', 'm', 'g', 'g', 'g', 'g',
                                            'g', 'r', 'r', 'r', 'r', 'r', 'b',
                                            'b', 'b', 'b', 'b'])
                else:
                    if ajust == True:
                        rects1 = ax0.bar(index, df_fin['deposit_Tartrazine'],
                                         bar_width, alpha=opacity,
                                         color=['c', 'c', 'c', 'c', 'c', 'c',
                                                'c', 'c', 'c', 'c', 'c', 'c',
                                                'c', 'c', 'c', 'm', 'm', 'm',
                                                'm', 'm', 'm', 'm', 'm', 'm',
                                                'm', 'm', 'm', 'm', 'm', 'm',
                                                'g', 'g', 'g', 'g', 'g', 'g',
                                                'g', 'g', 'g', 'g', 'g', 'g',
                                                'g', 'g', 'g', 'r', 'r', 'r',
                                                'r', 'r', 'r', 'r', 'r', 'r',
                                                'r', 'r', 'r', 'r', 'r', 'r',
                                                'b', 'b', 'b', 'b', 'b', 'b',
                                                'b', 'b', 'b', 'b', 'b', 'b',
                                                'b', 'b', 'b'])
                    else:
                        rects1 = ax0.bar(index, df_fin['deposit_Tartrazine'],
                                         bar_width, alpha=opacity,
                                         color=['c', 'c', 'c', 'c', 'c', 'c',
                                                'c', 'c', 'c', 'c', 'c', 'c',
                                                'c', 'c', 'c', 'm', 'm', 'm',
                                                'm', 'm', 'm', 'm', 'm', 'm',
                                                'm', 'm', 'm', 'm', 'm', 'm',
                                                'g', 'g', 'g', 'g', 'g', 'g',
                                                'g', 'g', 'g', 'g', 'g', 'g',
                                                'g', 'g', 'g', 'r', 'r', 'r',
                                                'r', 'r', 'r', 'r', 'r', 'r',
                                                'r', 'r', 'r', 'r', 'r', 'r',
                                                'b', 'b', 'b', 'b', 'b', 'b',
                                                'b', 'b', 'b', 'b', 'b', 'b',
                                                'b', 'b', 'b'])

                rects2 = ax0.bar(index + bar_width, df_fin['mean'], bar_width,
                                 alpha=opacity, color='y')
                # Mise en forme
                ax0.set_ylim(0, 8)
                ax0.set_xlim(0, len(df_fin))
                ax0.set_xticks(index + bar_width)
                df_fin['metamer'] = numpy.round(df_fin['metamer'], decimals=1)
                df_fin['lig'] = (df_fin['HS'] - 0.3 - df_fin['metamer']) / \
                                HSconv[var].a_cohort
                df_fin['lig'] = numpy.round(df_fin['lig'], decimals=2)
                ax0.set_xticklabels(df_fin['lig'].tolist(), rotation=90,
                                    fontsize='small')
                ax0.set_xlabel('age depuis ligulation')
                ax0.text(0.4, 7.5, var + ' - HS = ' + str(date),
                         bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10},
                         fontsize=14)
                fig.suptitle(
                    'Sim [deposit_Tartrazine] / Obs [mean] par ntop_cur pour ' + var,
                    fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
                # ax0.legend((rects1[0], rects1[15], rects1[30], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})

    # PLOT 2
    if plot2 is True:
        for var in varieties:
            if var == 'Tremie12' or var == 'Tremie13':
                df_scan_var = df_scan[df_scan['var'] == var]
                df_sim_var = df_sim[df_sim['var'] == var]
                df_scan_var.HS = list(map(str, df_scan_var.HS))
                df_sim_var.HS = list(map(str, df_sim_var.HS))
                df_all = df_sim_var.merge(
                    df_scan_var.ix[:, ['HS', 'ntop_cur', 'Area A_bl']],
                    how='outer')
                df_all = df_all[df_all['ntop_cur'] <= 5]
                # plot
                bar_width = 0.4;
                opacity = 0.4
                fig, axes = plt.subplots(nrows=1,
                                         ncols=2)  # fig, axes = plt.subplots(nrows=1, ncols=len(df_all['HS'].unique()))
                if var == 'Tremie12':
                    val = [[0, 10.98], [1, 12.63]]
                elif var == 'Tremie13':
                    val = [[0, 8.7], [1, 11.04]]
                for x, date in val:
                    df_all.HS = list(map(float, df_all.HS))
                    df_fin = df_all[df_all['HS'] == date]
                    if x == 0:
                        for f in [-0.4, -0.2, 0.2, 0.4]:
                            df_fin = df_fin.append(
                                df_all[df_all['HS'] == round(date + f, 2)])
                    else:
                        for f in [-0.4, -0.2, 0.2, 0.4, 0.5, 1, 1.5, 2, 2.5]:
                            df_fin = df_fin.append(
                                df_all[df_all['HS'] == round(date + f, 2)])
                    df_fin = df_fin.sort(['ntop_cur', 'HS'])
                    n_groups = len(df_fin)
                    index = numpy.arange(n_groups)
                    rects1 = axes[x].bar(index, df_fin['Area A_bl'], bar_width,
                                         alpha=opacity, color='y')

                    if date == 10.98 or date == 8.7:
                        rects2 = axes[x].bar(index + bar_width, df_fin['area'],
                                             bar_width, alpha=opacity,
                                             color=['c', 'c', 'c', 'c', 'c',
                                                    'm', 'm', 'm', 'm', 'm',
                                                    'g', 'g', 'g', 'g', 'g',
                                                    'r', 'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b', 'b'])
                    else:
                        rects2 = axes[x].bar(index + bar_width, df_fin['area'],
                                             bar_width, alpha=opacity,
                                             color=['c', 'c', 'c', 'c', 'c',
                                                    'c', 'c', 'c', 'c', 'c',
                                                    'c', 'c', 'c', 'c', 'c',
                                                    'm', 'm', 'm', 'm', 'm',
                                                    'm', 'm', 'm', 'm', 'm',
                                                    'm', 'm', 'm', 'm', 'm',
                                                    'g', 'g', 'g', 'g', 'g',
                                                    'g', 'g', 'g', 'g', 'g',
                                                    'g', 'g', 'g', 'g', 'g', 'r', 'r', 'r', 'r', 'r',
                                                    'r', 'r', 'r', 'r', 'r',
                                                    'r', 'r', 'r', 'r', 'r',
                                                    'b', 'b', 'b', 'b', 'b',
                                                    'b', 'b', 'b', 'b', 'b',
                                                    'b', 'b', 'b', 'b', 'b'])

                    # Mise en forme
                    axes[x].set_ylim(0, 30)
                    axes[x].set_xlim(0, len(df_fin))
                    axes[x].set_xticks(index + bar_width)
                    df_fin['metamer'] = numpy.round(df_fin['metamer'],
                                                    decimals=1)
                    df_fin['TT'] = numpy.round(df_fin['TT'], decimals=1)
                    df_fin['label'] = df_fin['HS'].astype(str) + ' - ' + df_fin[
                        'metamer'].astype(str)
                    axes[x].set_xticklabels(df_fin['label'].tolist(),
                                            rotation=90, fontsize='small')
                    if x == 0:
                        axes[x].set_xlabel('HS/ntop');
                        axes[x].set_ylabel('area (cm2)')
                    axes[x].text(0.4, 28.5, var + ' - HS = ' + str(date),
                                 bbox={'facecolor': '#FCF8F8', 'alpha': 0.6,
                                       'pad': 10}, fontsize=14)
                    fig.suptitle('Surface scan/sim par ntop_cur pour ' + var,
                                 fontsize=10)
                    fig.subplots_adjust(left=0.1, right=0.9, top=0.9,
                                        bottom=0.2)
                axes[x].legend((rects2[0], rects2[15], rects2[30], rects1[0]),
                               ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'),
                               bbox_to_anchor=[1.10, 1.12], prop={'size': 14})
            else:
                print('Pas de scan pour variete : ' + var)

    # PLOT 3
    if plot3 is True:
        for var in varieties:
            df_obs_var = df_obs[df_obs['name'] == var]
            df_obs_var['ntop_cur'] = df_obs_var['N feuille']
            df_sim_var = df_sim[df_sim['var'] == var]
            df_obs_var.HS = list(map(str, df_obs_var.HS))
            df_sim_var.HS = list(map(str, df_sim_var.HS))
            df_all = df_sim_var.merge(
                df_obs_var.ix[:, ['HS', 'ntop_cur', 'mean', 'sd']], how='outer')
            df_all = df_all[df_all['ntop_cur'] <= 3]
            df_all['tartrazine/area'] = df_all['deposit_Tartrazine'] / df_all[
                'area']
            df_all['mean/area'] = df_all['mean'] / df_all['area']
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            # plot
            bar_width = 0.4;
            opacity = 0.4
            if ajust == False:
                if var == 'Mercia':
                    # val = [[0,9.74],[1,12.8]]
                    val = [[0, 12.8]]
                elif var == 'Rht3':
                    # val = [[0,9.15],[1,12.48]]
                    val = [[0, 12.48]]
                elif var == 'Tremie12':
                    # val = [[0,10.98],[1,12.63]]
                    val = [[0, 12.63]]
                else:
                    # val = [[0,8.7],[1,11.04]]
                    val = [[0, 11.04]]
            else:
                val = [[0, round(hs_moyen, 2)]]
            fig, ax0 = plt.subplots(nrows=1, ncols=1)
            for x, date in val:
                df_all.HS = list(map(float, df_all.HS))
                df_all['HS'] = numpy.round(df_all['HS'], decimals=2)
                df_fin = df_all[df_all['HS'] == date]
                # for f in [-0.8,-0.6,-0.4,-0.2,+0.2,+0.4,+0.6,+0.8]:
                for f in [-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4,
                          4.5, 5]:
                    df_fin = df_fin.append(
                        df_all[df_all['HS'] == round(date + f, 2)])
                df_fin = df_fin.sort(['ntop_cur', 'HS'])
                n_groups = len(df_fin)
                index = numpy.arange(n_groups)

                # rects1 = ax0.bar(index, df_fin['tartrazine/area'], bar_width, alpha=opacity, color=['c','m','g','r','b'])
                if date == 10.98 or date == 8.7 or date == 9.15 or date == 9.74:
                    rects1 = ax0.bar(index, df_fin['tartrazine/area'],
                                     bar_width, alpha=opacity,
                                     color=['c', 'c', 'c', 'c', 'c', 'm', 'm',
                                            'm', 'm', 'm', 'g', 'g', 'g', 'g',
                                            'g', 'r', 'r', 'r', 'r', 'r', 'b',
                                            'b', 'b', 'b', 'b'])
                else:
                    if ajust == True:
                        rects1 = ax0.bar(index, df_fin['tartrazine/area'],
                                         bar_width, alpha=opacity,
                                         color=['c', 'c', 'c', 'c', 'c', 'c',
                                                'c', 'c', 'c', 'c', 'c', 'c',
                                                'c', 'c', 'c', 'm', 'm', 'm',
                                                'm', 'm', 'm', 'm', 'm', 'm',
                                                'm', 'm', 'm', 'm', 'm', 'm',
                                                'g', 'g', 'g', 'g', 'g', 'g',
                                                'g', 'g', 'g', 'g', 'g', 'g',
                                                'g', 'g', 'g', 'r', 'r', 'r',
                                                'r', 'r', 'r', 'r', 'r', 'r',
                                                'r', 'r', 'r', 'r', 'r', 'r',
                                                'b', 'b', 'b', 'b', 'b', 'b',
                                                'b', 'b', 'b', 'b', 'b', 'b',
                                                'b', 'b', 'b'])
                    else:
                        rects1 = ax0.bar(index, df_fin['tartrazine/area'],
                                         bar_width, alpha=opacity,
                                         color=['c', 'c', 'c', 'c', 'c', 'c',
                                                'c', 'c', 'c', 'c', 'c', 'c',
                                                'c', 'c', 'c', 'm', 'm', 'm',
                                                'm', 'm', 'm', 'm', 'm', 'm',
                                                'm', 'm', 'm', 'm', 'm', 'm',
                                                'g', 'g', 'g', 'g', 'g', 'g',
                                                'g', 'g', 'g', 'g', 'g', 'g',
                                                'g', 'g', 'g', 'r', 'r', 'r',
                                                'r', 'r', 'r', 'r', 'r', 'r',
                                                'r', 'r', 'r', 'r', 'r', 'r',
                                                'b', 'b', 'b', 'b', 'b', 'b',
                                                'b', 'b', 'b', 'b', 'b', 'b',
                                                'b', 'b', 'b'])

                rects2 = ax0.bar(index + bar_width, df_fin['mean/area'],
                                 bar_width, alpha=opacity, color='y')
                # Mise en forme
                ax0.set_ylim(0, 0.5)
                ax0.set_xlim(0, len(df_fin))
                ax0.set_xticks(index + bar_width)
                df_fin['lig'] = (df_fin['HS'] - 0.3 - df_fin['metamer']) / \
                                HSconv[var].a_cohort
                df_fin['lig'] = numpy.round(df_fin['lig'], decimals=2)
                ax0.set_xticklabels(df_fin['lig'].tolist(), rotation=90,
                                    fontsize='small')
                if x == 0:
                    ax0.set_xlabel('age depuis ligulation')
                ax0.text(0.4, 0.46, var + ' - HS = ' + str(date),
                         bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10},
                         fontsize=14)
                fig.suptitle(
                    'Sim [deposit_Tartrazine/area] / Obs [mean/area] par ntop_cur pour ' + var,
                    fontsize=10)
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
                # ax0.legend((rects1[0], rects1[15], rects1[30], rects2[0]), ('Sim F1', 'Sim F2', 'Sim F3', 'Obs'), bbox_to_anchor=[1.10, 1.12], prop={'size':14} )

    return df_scan, df_obs, df_sim


# comparer ntop_cur=1 pour des ages de feuilles identiques (Fameuse figure corinne)
def sim_HScom(var='Mercia', stade='T1', axis='MS'):
    df_sim = pandas.DataFrame()
    x = 1
    while x <= 1:
        if var == 'Tremie12':
            n = 30
        else:
            n = 30
        npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n, axis=axis,
                                     to_csv=False)
        print('var = ' + var + ' stade = ' + stade + ' - nbre de plantes sim = ' + str(
            npl) + ' - nbr simul = ' + str(x) + '/5')
        dfmoy['var'] = var
        if stade == 'T2':
            dfmoy['stade'] = 'T2_' + var
        else:
            dfmoy['stade'] = stade
        df_sim = df_sim.append(dfmoy)
        x += 1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'stade', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    return df_sim_gr


def plot_HScom(varieties=['Mercia', 'Rht3', 'Tremie12', 'Tremie13'],
               appdate_lst=['T1', 'T2'], delta_opt=True, axis='MS', nplants=30,
               plot1=True, plot2=True):  # delta = FALSE (6) ou TRUE
    '''
    Comparer les F1 T1 et F1 T2 pour les 4 varietes au HS de Mercia, au HS de Rht3, au HS de Tremie12 et au HS de Tremie13
    Methode : (ex. pour Mercia, T1)
        delta_HS_Mercia = HS T1 Mercia - nff moyen de Mercia
    Pour les autres varietes, faire :
        HS = nff moyen + delta_HS_Mercia

    Plot 1 = Tartrazine
    Plot 2 = Tartrazine / area
    '''
    df_all = pandas.DataFrame()
    # Ref Mercia
    for name in varieties:
        for appdate in appdate_lst:
            # HS
            date, hs = HS_applications[name][appdate]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=name,
                                                      nplants=nplants)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            # calcul du delta
            # delta = hs - mean['nff']
            if delta_opt == True:
                delta = hs - mean['nff']
            else:
                delta = 6
            # verif
            print('nff moyen ref ' + name + ' = ' + str(mean['nff']))
            print('delta ref ' + name + ' = ' + str(delta))
            # simulation
            # hs_ref = hs+delta
            df_ref = sim_HScom(var=name, stade=str(hs), axis=axis)
            df_ref['dim'] = appdate + '_' + name
            df_ref['delta'] = delta
            df_all = df_all.append(df_ref)

            if name == 'Mercia':
                name_other = ['Rht3', 'Tremie12', 'Tremie13']
            elif name == 'Rht3':
                name_other = ['Mercia', 'Tremie12', 'Tremie13']
            elif name == 'Tremie12':
                name_other = ['Mercia', 'Rht3', 'Tremie13']
            elif name == 'Tremie13':
                name_other = ['Mercia', 'Rht3', 'Tremie12']
            for nameo in name_other:
                # nff moyen
                adel = Reconstructions.get_reconstruction(name=nameo,
                                                          nplants=nplants)
                df_phenT2 = adel.phenT()
                nffs = df_phenT2.set_index('plant').groupby(level=0).count()[
                    'n']
                df_phenT2['nff'] = [nffs[v] for v in df_phenT2['plant']]
                mean_other = df_phenT2.mean()
                # HS
                HS_other = mean_other['nff'] + delta
                # verif
                print('nff moyen ' + nameo + ' = ' + str(mean_other['nff']))
                print('HS other ' + nameo + ' = ' + str(HS_other))
                # simulation
                df_other = sim_HScom(var=nameo, stade=str(HS_other), axis=axis)
                df_other['dim'] = appdate + '_' + name
                df_all = df_all.append(df_other)

    if plot1 == True:
        bar_width = 0.2;
        opacity = 0.4
        fig, axes = plt.subplots(nrows=1, ncols=2)

        for x, date in [[0, 'T1'], [1, 'T2']]:

            df_sim = df_all[df_all['ntop_cur'] == 1]

            df_sim_HSMercia = df_sim[df_sim['dim'] == date + '_Mercia'];
            df_sim_HSMercia = df_sim_HSMercia.reset_index()

            df_sim_HSRht3 = df_sim[df_sim['dim'] == date + '_Rht3'];
            df_sim_HSRht3 = df_sim_HSRht3.sort(['var']);
            df_sim_HSRht3 = df_sim_HSRht3.reset_index()

            df_sim_HSTremie12 = df_sim[df_sim['dim'] == date + '_Tremie12'];
            df_sim_HSTremie12 = df_sim_HSTremie12.sort(['var']);
            df_sim_HSTremie12 = df_sim_HSTremie12.reset_index()

            df_sim_HSTremie13 = df_sim[df_sim['dim'] == date + '_Tremie13'];
            df_sim_HSTremie13 = df_sim_HSTremie13.sort(['var']);
            df_sim_HSTremie13 = df_sim_HSTremie13.reset_index()

            n_groups = len(df_sim_HSMercia)
            index = numpy.arange(n_groups)

            rects1 = axes[x].bar(index, df_sim_HSMercia['deposit_Tartrazine'],
                                 bar_width, alpha=opacity, color=['r'])
            rects2 = axes[x].bar(index + bar_width,
                                 df_sim_HSRht3['deposit_Tartrazine'], bar_width,
                                 alpha=opacity, color=['g'])
            rects3 = axes[x].bar(index + 2 * bar_width,
                                 df_sim_HSTremie12['deposit_Tartrazine'],
                                 bar_width, alpha=opacity, color=['b'])
            rects4 = axes[x].bar(index + 3 * bar_width,
                                 df_sim_HSTremie13['deposit_Tartrazine'],
                                 bar_width, alpha=opacity, color=['m'])

            # Mise en forme
            axes[x].set_ylim(0, 6.5)
            axes[x].set_xlim(0, len(df_sim_HSMercia))
            axes[x].set_xticks(index + bar_width)
            # Label au dessus des colonnes
            df_sim_HSMercia['leaf_emergence'] = numpy.round(
                df_sim_HSMercia['leaf_emergence'], decimals=1)
            df_sim_HSRht3['leaf_emergence'] = numpy.round(
                df_sim_HSRht3['leaf_emergence'], decimals=1)
            df_sim_HSTremie12['leaf_emergence'] = numpy.round(
                df_sim_HSTremie12['leaf_emergence'], decimals=1)
            df_sim_HSTremie13['leaf_emergence'] = numpy.round(
                df_sim_HSTremie13['leaf_emergence'], decimals=1)
            df_sim_HSMercia['HS'] = numpy.round(df_sim_HSMercia['HS'],
                                                decimals=1)
            df_sim_HSRht3['HS'] = numpy.round(df_sim_HSRht3['HS'], decimals=1)
            df_sim_HSTremie12['HS'] = numpy.round(df_sim_HSTremie12['HS'],
                                                  decimals=1)
            df_sim_HSTremie13['HS'] = numpy.round(df_sim_HSTremie13['HS'],
                                                  decimals=1)
            df_sim_HSMercia['label'] = df_sim_HSMercia['HS'].astype(
                str) + ' / ' + df_sim_HSMercia['leaf_emergence'].astype(str)
            df_sim_HSRht3['label'] = df_sim_HSRht3['HS'].astype(str) + ' / ' + \
                                     df_sim_HSRht3['leaf_emergence'].astype(str)
            df_sim_HSTremie12['label'] = df_sim_HSTremie12['HS'].astype(
                str) + ' / ' + df_sim_HSTremie12['leaf_emergence'].astype(str)
            df_sim_HSTremie13['label'] = df_sim_HSTremie13['HS'].astype(
                str) + ' / ' + df_sim_HSTremie13['leaf_emergence'].astype(str)

            list1 = df_sim_HSMercia['label']
            list2 = df_sim_HSRht3['label']
            list3 = df_sim_HSTremie12['label']
            list4 = df_sim_HSTremie13['label']

            def autolabel(rects):
                # attach some text labels
                if rects == rects1:
                    list = list1
                elif rects == rects2:
                    list = list2
                elif rects == rects3:
                    list = list3
                elif rects == rects4:
                    list = list4
                for rect in rects:
                    height = rect.get_height()
                    axes[x].text(rect.get_x() + rect.get_width() / 2.,
                                 1.05 * height, list[rect.get_x()], ha='center',
                                 va='bottom', fontsize=8, rotation=90)

            autolabel(rects1)
            autolabel(rects2)
            autolabel(rects3)
            autolabel(rects4)
            # Label x
            df_label = df_sim_HSMercia.append(df_sim_HSRht3)
            df_label = df_label.append(df_sim_HSTremie12)
            df_label = df_label.append(df_sim_HSTremie13)
            delta = df_label['delta'].dropna()
            axes[x].set_xticklabels(delta.tolist(), rotation=90,
                                    fontsize='small')
            if x == 0:
                axes[x].set_xlabel('delta_HS')
            axes[x].text(0.4, 6.2, 'treatment = ' + str(date),
                         bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10},
                         fontsize=12)

            fig.suptitle('', fontsize=10)
            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        axes[x].legend((rects1[0], rects2[0], rects3[0], rects4[0]),
                       ('Mercia', 'Rht3', 'Tremie12', 'Tremie13'),
                       bbox_to_anchor=[1.10, 1.12], prop={'size': 14})

    if plot2 == True:
        bar_width = 0.2;
        opacity = 0.4
        fig, axes = plt.subplots(nrows=1, ncols=2)

        for x, date in [[0, 'T1'], [1, 'T2']]:

            df_sim = df_all[df_all['ntop_cur'] == 1]
            df_sim['deposit_Tartrazine/area'] = df_sim['deposit_Tartrazine'] / \
                                                df_sim['area']

            df_sim_HSMercia = df_sim[df_sim['dim'] == date + '_Mercia'];
            df_sim_HSMercia = df_sim_HSMercia.reset_index()

            df_sim_HSRht3 = df_sim[df_sim['dim'] == date + '_Rht3'];
            df_sim_HSRht3 = df_sim_HSRht3.sort(['var']);
            df_sim_HSRht3 = df_sim_HSRht3.reset_index()

            df_sim_HSTremie12 = df_sim[df_sim['dim'] == date + '_Tremie12'];
            df_sim_HSTremie12 = df_sim_HSTremie12.sort(['var']);
            df_sim_HSTremie12 = df_sim_HSTremie12.reset_index()

            df_sim_HSTremie13 = df_sim[df_sim['dim'] == date + '_Tremie13'];
            df_sim_HSTremie13 = df_sim_HSTremie13.sort(['var']);
            df_sim_HSTremie13 = df_sim_HSTremie13.reset_index()

            n_groups = len(df_sim_HSMercia)
            index = numpy.arange(n_groups)

            rects1 = axes[x].bar(index,
                                 df_sim_HSMercia['deposit_Tartrazine/area'],
                                 bar_width, alpha=opacity, color=['r'])
            rects2 = axes[x].bar(index + bar_width,
                                 df_sim_HSRht3['deposit_Tartrazine/area'],
                                 bar_width, alpha=opacity, color=['g'])
            rects3 = axes[x].bar(index + 2 * bar_width,
                                 df_sim_HSTremie12['deposit_Tartrazine/area'],
                                 bar_width, alpha=opacity, color=['b'])
            rects4 = axes[x].bar(index + 3 * bar_width,
                                 df_sim_HSTremie13['deposit_Tartrazine/area'],
                                 bar_width, alpha=opacity, color=['m'])

            # Mise en forme
            axes[x].set_ylim(0, 0.6)
            axes[x].set_xlim(0, len(df_sim_HSMercia))
            axes[x].set_xticks(index + bar_width)
            # Label au dessus des colonnes
            df_sim_HSMercia['leaf_emergence'] = numpy.round(
                df_sim_HSMercia['leaf_emergence'], decimals=1)
            df_sim_HSRht3['leaf_emergence'] = numpy.round(
                df_sim_HSRht3['leaf_emergence'], decimals=1)
            df_sim_HSTremie12['leaf_emergence'] = numpy.round(
                df_sim_HSTremie12['leaf_emergence'], decimals=1)
            df_sim_HSTremie13['leaf_emergence'] = numpy.round(
                df_sim_HSTremie13['leaf_emergence'], decimals=1)
            df_sim_HSMercia['HS'] = numpy.round(df_sim_HSMercia['HS'],
                                                decimals=1)
            df_sim_HSRht3['HS'] = numpy.round(df_sim_HSRht3['HS'], decimals=1)
            df_sim_HSTremie12['HS'] = numpy.round(df_sim_HSTremie12['HS'],
                                                  decimals=1)
            df_sim_HSTremie13['HS'] = numpy.round(df_sim_HSTremie13['HS'],
                                                  decimals=1)

            df_sim_HSMercia['label'] = df_sim_HSMercia['HS'].astype(
                str) + ' / ' + df_sim_HSMercia['leaf_emergence'].astype(str)

            df_sim_HSRht3['label'] = df_sim_HSRht3['HS'].astype(str) + ' / ' + \
                                     df_sim_HSRht3['leaf_emergence'].astype(str)
            df_sim_HSTremie12['label'] = df_sim_HSTremie12['HS'].astype(
                str) + ' / ' + df_sim_HSTremie12['leaf_emergence'].astype(str)

            df_sim_HSTremie13['label'] = df_sim_HSTremie13['HS'].astype(
                str) + ' / ' + df_sim_HSTremie13['leaf_emergence'].astype(str)

            list1 = df_sim_HSMercia['label']
            list2 = df_sim_HSRht3['label']
            list3 = df_sim_HSTremie12['label']
            list4 = df_sim_HSTremie13['label']

            def autolabel(rects):
                # attach some text labels
                if rects == rects1:
                    list = list1
                elif rects == rects2:
                    list = list2
                elif rects == rects3:
                    list = list3
                elif rects == rects4:
                    list = list4
                for rect in rects:
                    height = rect.get_height()
                    axes[x].text(rect.get_x() + rect.get_width() / 2.,
                                 1.05 * height, list[rect.get_x()], ha='center',
                                 va='bottom', fontsize=8, rotation=90)

            autolabel(rects1)
            autolabel(rects2)
            autolabel(rects3)
            autolabel(rects4)
            # Label x
            df_label = df_sim_HSMercia.append(df_sim_HSRht3)
            df_label = df_label.append(df_sim_HSTremie12)
            df_label = df_label.append(df_sim_HSTremie13)
            delta = df_label['delta'].dropna()
            axes[x].set_xticklabels(delta.tolist(), rotation=90,
                                    fontsize='small')

            if x == 0:
                axes[x].set_xlabel('delta_HS')
            axes[x].text(0.2, 0.56, 'treatment = ' + str(date),
                         bbox={'facecolor': '#FCF8F8', 'alpha': 0.6, 'pad': 10},
                         fontsize=12)

            fig.suptitle('', fontsize=10)
            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        axes[x].legend((rects1[0], rects2[0], rects3[0], rects4[0]),
                       ('Mercia', 'Rht3', 'Tremie12', 'Tremie13'),
                       bbox_to_anchor=[1.10, 1.12], prop={'size': 14})

    df_sim.to_csv('HS_com_synth_all.csv')
    return df_sim


# efficacy_T2/T1
# plot tartrazine/area par HS (HSnff-6 a HSnff+6) pour ntop=1 a3
def simulation_efficacy(name='Mercia', hs=12, n_sim=5, n_plt=200, axis='MS'):
    def frange(x, y, jump):
        while x < y:
            yield x
            x += jump

    hs_lst = frange(hs - 6, hs + 6, ((float((hs + 6) - (hs - 6))) / 20))
    lst = []
    for hs_ in hs_lst:
        lst.append([name, str(hs_)])

    if name == 'Tremie12':  # cas particulier de tremie12 n=30
        n = 30
    else:
        n = n_plt

    df_sim = pandas.DataFrame()
    x = 1
    while x <= n_sim:
        for var, stade in lst:
            npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n,
                                         axis=axis, to_csv=False)
            print('var = ' + var + ' stade = ' + stade + ' - nbre de plantes sim = ' + str(
                npl))
            dfmoy['var'] = var
            df_sim = df_sim.append(dfmoy)
        x += 1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'ntop']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    return df_sim_gr


def plot_efficacy(varieties=['Mercia', 'Rht3', 'Tremie12', 'Tremie13'], n_sim=1,
                  n_plt=30, axis='MS', plot_tartrazine=True,
                  plot_intercept=True, plot_cov=True, plot_protect=True,
                  plot_global=True, csv=True, ntop_lst=[1, 2, 3], data=True):
    '''
    !!! Tous les plots sont declines en 2 versions :
    HS et HS-nff moyen dans une fourchette de -6 a +6 par pas de temps calcule d'environ 0.5
    plot tartrazine : simulation 'deposit_Tratrazine'
    plot intercept : simulation 'interception efficacy'
    RAPPEL : interception efficacy = deposit_Tartrazine/area
    plot cov : simulation 'coverage efficacy'
    coverage efficacy = interception efficacy * exposition
    plot protect : simulation 'protection efficacy'
    protection efficacy = (1-lifetime) * coverage efficacy
    plot global : simulation 'somme protection efficacy'
    '''
    df_sim = pandas.DataFrame()
    for var in varieties:
        # nff moyen
        adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
        df_phenT = adel.phenT()
        nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
        df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
        mean = df_phenT.mean()
        hs_moyen = mean['nff']
        # sim
        df_sim_var = simulation_efficacy(name=var, hs=hs_moyen, n_sim=n_sim,
                                         n_plt=n_plt, axis=axis)
        df_sim_var['deposit_Tartrazine/area'] = df_sim_var[
                                                    'deposit_Tartrazine'] / \
                                                df_sim_var['area']
        df_sim_var['coverage_efficacy'] = df_sim_var[
                                              'deposit_Tartrazine/area'] * \
                                          df_sim_var['exposition']
        df_sim_var['protection_efficacy'] = (1 - df_sim_var['lifetime']) * \
                                            df_sim_var['coverage_efficacy']
        df_sim_var['HS-nffmean'] = df_sim_var['HS'] - hs_moyen
        df_sim = df_sim.append(df_sim_var)
    if csv == True:
        df_sim.to_csv('efficacy_data.csv')

    # obs
    df_obs = idata.dye_interception()
    '''df_obs.rename(columns={'name':'var'}, inplace=True)
    df_obs.rename(columns={'N feuille':'ntop_cur'}, inplace=True)
    df_obs.rename(columns={'HS':'HS_obs'}, inplace=True)'''

    if plot_tartrazine == True:
        # x = hs
        plt.figure()
        for var in varieties:
            if var == 'Mercia':
                color = 'r'
            elif var == 'Rht3':
                color = 'g'
            elif var == 'Tremie12':
                color = 'b'
            elif var == 'Tremie13':
                color = 'm'
            df_sim_tartra = df_sim[df_sim['var'] == var]
            # plot
            for ntop in ntop_lst:
                df_sim_ntop = df_sim_tartra[df_sim_tartra['ntop'] == ntop]
                # gestion de la legende
                if var == 'Mercia':
                    label = 'true F' + str(ntop)
                else:
                    label = '_nolegend_'
                if ntop == 1:
                    picto = '-'
                elif ntop == 2:
                    picto = '--'
                elif ntop == 3:
                    picto = '-.'
                elif ntop == 4:
                    picto = ':'
                # plot par ntop
                plt.plot(df_sim_ntop['HS'], df_sim_ntop['deposit_Tartrazine'],
                         picto + color, label=label)
            # fleche hsT2
            date, hsT2 = HS_applications[var]['T2']
            plt.annotate('', xy=(hsT2, 0), xytext=(hsT2, -0.7),
                         arrowprops=dict(color=color, arrowstyle="->",
                                         connectionstyle="arc3"))
            # observation
            if data == True:
                df_obs_var = df_obs[df_obs['name'] == var]
                df_obs_var = df_obs_var[df_obs_var['treatment'] == 'T2']
                df_obs_var['x'] = hsT2
                for n_feuille in [1, 2, 3, 4]:
                    df_obs_var_f = df_obs_var[
                        df_obs_var['N feuille'] == n_feuille]
                    if n_feuille == 1:
                        colorf = 'c'
                    elif n_feuille == 2:
                        colorf = 'm'
                    elif n_feuille == 3:
                        colorf = 'g'
                    elif n_feuille == 4:
                        colorf = 'r'
                    plt.plot(df_obs_var_f['x'], df_obs_var_f['mean'],
                             'o' + colorf)

        # mise en forme
        plt.xlabel('haun stage');
        plt.ylabel('deposit_tartrazine')
        plt.grid(False)  # enlever la grille
        # plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})

        # x = hs - nff moyen
        plt.figure()
        for var in varieties:
            if var == 'Mercia':
                color = 'r'
            elif var == 'Rht3':
                color = 'g'
            elif var == 'Tremie12':
                color = 'b'
            elif var == 'Tremie13':
                color = 'm'
            df_sim_tartra2 = df_sim[df_sim['var'] == var]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            # plot
            for ntop in ntop_lst:
                df_sim_ntop = df_sim_tartra2[df_sim_tartra2['ntop'] == ntop]
                # gestion de la legende
                if var == 'Mercia':
                    label = 'true F' + str(ntop)
                else:
                    label = '_nolegend_'
                if ntop == 1:
                    picto = '-'
                elif ntop == 2:
                    picto = '--'
                elif ntop == 3:
                    picto = '-.'
                elif ntop == 4:
                    picto = ':'
                # plot par ntop
                plt.plot(df_sim_ntop['HS-nffmean'],
                         df_sim_ntop['deposit_Tartrazine'], picto + color,
                         label=label)
            # fleche hsT2 - nff moyen
            hs_new = hsT2 - hs_moyen
            plt.annotate('', xy=(hs_new, 0), xytext=(hs_new, -0.7),
                         arrowprops=dict(color=color, arrowstyle="->",
                                         connectionstyle="arc3"))
            # observation
            if data == True:
                df_obs_var = df_obs[df_obs['name'] == var]
                df_obs_var = df_obs_var[df_obs_var['treatment'] == 'T2']
                df_obs_var['x'] = hs_new
                for n_feuille in [1, 2, 3, 4]:
                    df_obs_var_f = df_obs_var[
                        df_obs_var['N feuille'] == n_feuille]
                    if n_feuille == 1:
                        colorf = 'c'
                    elif n_feuille == 2:
                        colorf = 'm'
                    elif n_feuille == 3:
                        colorf = 'g'
                    elif n_feuille == 4:
                        colorf = 'r'
                    plt.plot(df_obs_var_f['x'], df_obs_var_f['mean'],
                             'o' + colorf)
        # mise en forme
        plt.xlabel('haun stage - nff moyen');
        plt.ylabel('deposit_tartrazine')
        plt.grid(False)  # enlever la grille
        # plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})

    if plot_intercept == True:
        # x = hs
        plt.figure()
        for var in varieties:
            if var == 'Mercia':
                color = 'r'
            elif var == 'Rht3':
                color = 'g'
            elif var == 'Tremie12':
                color = 'b'
            elif var == 'Tremie13':
                color = 'm'
            df_sim_intercept = df_sim[df_sim['var'] == var]
            # plot
            for ntop in ntop_lst:
                df_sim_ntop = df_sim_intercept[df_sim_intercept['ntop'] == ntop]
                # gestion de la legende
                if var == 'Mercia':
                    label = 'true F' + str(ntop)
                else:
                    label = '_nolegend_'
                if ntop == 1:
                    picto = '-'
                elif ntop == 2:
                    picto = '--'
                elif ntop == 3:
                    picto = '-.'
                elif ntop == 4:
                    picto = ':'
                # plot par ntop
                plt.plot(df_sim_ntop['HS'],
                         df_sim_ntop['deposit_Tartrazine/area'], picto + color,
                         label=label)
            # fleche hsT2
            date, hsT2 = HS_applications[var]['T2']
            plt.annotate('', xy=(hsT2, 0), xytext=(hsT2, -0.03),
                         arrowprops=dict(color=color, arrowstyle="->",
                                         connectionstyle="arc3"))
        # mise en forme
        plt.xlabel('haun stage');
        plt.ylabel('interception efficacy')
        plt.grid(False)  # enlever la grille
        # plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})

        # x = hs - nffmoyen
        plt.figure()
        for var in varieties:
            if var == 'Mercia':
                color = 'r'
            elif var == 'Rht3':
                color = 'g'
            elif var == 'Tremie12':
                color = 'b'
            elif var == 'Tremie13':
                color = 'm'
            df_sim_intercept2 = df_sim[df_sim['var'] == var]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            # plot
            for ntop in ntop_lst:
                df_sim_ntop = df_sim_intercept2[
                    df_sim_intercept2['ntop'] == ntop]
                # gestion de la legende
                if var == 'Mercia':
                    label = 'true F' + str(ntop)
                else:
                    label = '_nolegend_'
                if ntop == 1:
                    picto = '-'
                elif ntop == 2:
                    picto = '--'
                elif ntop == 3:
                    picto = '-.'
                elif ntop == 4:
                    picto = ':'
                # plot par ntop
                plt.plot(df_sim_ntop['HS-nffmean'],
                         df_sim_ntop['deposit_Tartrazine/area'], picto + color,
                         label=label)
            # fleche hsT2 - nff moyen
            hs_new = hsT2 - hs_moyen
            plt.annotate('', xy=(hs_new, 0), xytext=(hs_new, -0.03),
                         arrowprops=dict(color=color, arrowstyle="->",
                                         connectionstyle="arc3"))
        # mise en forme
        plt.xlabel('haun stage - nff moyen');
        plt.ylabel('interception efficacy')
        plt.grid(False)  # enlever la grille
        # plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})

    if plot_cov == True:
        # x = hs
        plt.figure()
        for var in varieties:
            if var == 'Mercia':
                color = 'r'
            elif var == 'Rht3':
                color = 'g'
            elif var == 'Tremie12':
                color = 'b'
            elif var == 'Tremie13':
                color = 'm'
            df_sim_cov = df_sim[df_sim['var'] == var]
            # plot
            for ntop in ntop_lst:
                df_sim_ntop = df_sim_cov[df_sim_cov['ntop'] == ntop]
                # gestion de la legende
                if var == 'Mercia':
                    label = 'true F' + str(ntop)
                else:
                    label = '_nolegend_'
                if ntop == 1:
                    picto = '-'
                elif ntop == 2:
                    picto = '--'
                elif ntop == 3:
                    picto = '-.'
                elif ntop == 4:
                    picto = ':'
                # plot par ntop
                plt.plot(df_sim_ntop['HS'], df_sim_ntop['coverage_efficacy'],
                         picto + color, label=label)
            # fleche hsT2
            date, hsT2 = HS_applications[var]['T2']
            plt.annotate('', xy=(hsT2, 0), xytext=(hsT2, -0.7),
                         arrowprops=dict(color=color, arrowstyle="->",
                                         connectionstyle="arc3"))
        # mise en forme
        plt.xlabel('haun stage');
        plt.ylabel('coverage_efficacy')
        plt.grid(False)  # enlever la grille
        # plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})

        # x = hs - nff moyen
        plt.figure()
        for var in varieties:
            if var == 'Mercia':
                color = 'r'
            elif var == 'Rht3':
                color = 'g'
            elif var == 'Tremie12':
                color = 'b'
            elif var == 'Tremie13':
                color = 'm'
            df_sim_cov2 = df_sim[df_sim['var'] == var]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            # plot
            for ntop in ntop_lst:
                df_sim_ntop = df_sim_cov2[df_sim_cov2['ntop'] == ntop]
                # gestion de la legende
                if var == 'Mercia':
                    label = 'true F' + str(ntop)
                else:
                    label = '_nolegend_'
                if ntop == 1:
                    picto = '-'
                elif ntop == 2:
                    picto = '--'
                elif ntop == 3:
                    picto = '-.'
                elif ntop == 4:
                    picto = ':'
                # plot par ntop
                plt.plot(df_sim_ntop['HS-nffmean'],
                         df_sim_ntop['coverage_efficacy'], picto + color,
                         label=label)
            # fleche hsT2 - nff moyen
            hs_new = hsT2 - hs_moyen
            plt.annotate('', xy=(hs_new, 0), xytext=(hs_new, -0.03),
                         arrowprops=dict(color=color, arrowstyle="->",
                                         connectionstyle="arc3"))
        # mise en forme
        plt.xlabel('haun stage - nff moyen');
        plt.ylabel('coverage_efficacy')
        plt.grid(False)  # enlever la grille
        # plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})

    if plot_protect == True:
        # x = hs
        plt.figure()
        for var in varieties:
            if var == 'Mercia':
                color = 'r'
            elif var == 'Rht3':
                color = 'g'
            elif var == 'Tremie12':
                color = 'b'
            elif var == 'Tremie13':
                color = 'm'
            df_sim_protect = df_sim[df_sim['var'] == var]
            # plot
            for ntop in ntop_lst:
                df_sim_ntop = df_sim_protect[df_sim_protect['ntop'] == ntop]
                # gestion de la legende
                if var == 'Mercia':
                    label = 'true F' + str(ntop)
                else:
                    label = '_nolegend_'
                if ntop == 1:
                    picto = '-'
                elif ntop == 2:
                    picto = '--'
                elif ntop == 3:
                    picto = '-.'
                elif ntop == 4:
                    picto = ':'
                # plot par ntop
                plt.plot(df_sim_ntop['HS'], df_sim_ntop['protection_efficacy'],
                         picto + color, label=label)
            # fleche hsT2
            date, hsT2 = HS_applications[var]['T2']
            plt.annotate('', xy=(hsT2, 0), xytext=(hsT2, -0.7),
                         arrowprops=dict(color=color, arrowstyle="->",
                                         connectionstyle="arc3"))
        # mise en forme
        plt.xlabel('haun stage');
        plt.ylabel('protection_efficacy')
        plt.grid(False)  # enlever la grille
        # plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})

        # x = hs - nff moyen
        plt.figure()
        for var in varieties:
            if var == 'Mercia':
                color = 'r'
            elif var == 'Rht3':
                color = 'g'
            elif var == 'Tremie12':
                color = 'b'
            elif var == 'Tremie13':
                color = 'm'
            df_sim_protect2 = df_sim[df_sim['var'] == var]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            # plot
            for ntop in ntop_lst:
                df_sim_ntop = df_sim_protect2[df_sim_protect2['ntop'] == ntop]
                # gestion de la legende
                if var == 'Mercia':
                    label = 'true F' + str(ntop)
                else:
                    label = '_nolegend_'
                if ntop == 1:
                    picto = '-'
                elif ntop == 2:
                    picto = '--'
                elif ntop == 3:
                    picto = '-.'
                elif ntop == 4:
                    picto = ':'
                # plot par ntop
                plt.plot(df_sim_ntop['HS-nffmean'],
                         df_sim_ntop['protection_efficacy'], picto + color,
                         label=label)
            # fleche hsT2 - nff moyen
            hs_new = hsT2 - hs_moyen
            plt.annotate('', xy=(hs_new, 0), xytext=(hs_new, -0.03),
                         arrowprops=dict(color=color, arrowstyle="->",
                                         connectionstyle="arc3"))
        # mise en forme
        plt.xlabel('haun stage - nff moyen');
        plt.ylabel('protection_efficacy')
        plt.grid(False)  # enlever la grille
        # plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})

    if plot_global == True:
        # x = hs
        plt.figure()
        for var in varieties:
            if var == 'Mercia':
                color = 'r'
            elif var == 'Rht3':
                color = 'g'
            elif var == 'Tremie12':
                color = 'b'
            elif var == 'Tremie13':
                color = 'm'
            df_sim_protect = df_sim[df_sim['var'] == var]
            # traitement dataframe
            df_sim_protect_all = pandas.DataFrame()
            for ntop in ntop_lst:
                df_sim_ntop = df_sim_protect[df_sim_protect['ntop'] == ntop]
                df_sim_protect_all = df_sim_protect_all.append(df_sim_ntop)
            df_sim_protect_all = df_sim_protect_all.groupby(['HS']).sum()
            df_sim_protect_all = df_sim_protect_all.reset_index()  # gestion de la legende
            if var == 'Mercia':
                label = 'Somme des protections efficacy'
            else:
                label = '_nolegend_'
            # plot
            plt.plot(df_sim_protect_all['HS'],
                     df_sim_protect_all['protection_efficacy'], '-' + color,
                     label=label)
            # fleche hsT2
            date, hsT2 = HS_applications[var]['T2']
            plt.annotate('', xy=(hsT2, 0), xytext=(hsT2, -0.7),
                         arrowprops=dict(color=color, arrowstyle="->",
                                         connectionstyle="arc3"))
        # mise en forme
        plt.xlabel('haun stage');
        plt.ylabel('global_efficacy')
        plt.grid(False)  # enlever la grille
        # plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})

        # x = hs - nff moyen
        plt.figure()
        for var in varieties:
            if var == 'Mercia':
                color = 'r'
            elif var == 'Rht3':
                color = 'g'
            elif var == 'Tremie12':
                color = 'b'
            elif var == 'Tremie13':
                color = 'm'
            df_sim_protect = df_sim[df_sim['var'] == var]
            # nff moyen
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            # HS T2 - nff moyen
            date, hsT2 = HS_applications[var]['T2']
            # traitement dataframe
            df_sim_protect_all = pandas.DataFrame()
            for ntop in ntop_lst:
                df_sim_ntop = df_sim_protect[df_sim_protect['ntop'] == ntop]
                df_sim_protect_all = df_sim_protect_all.append(df_sim_ntop)
            df_sim_protect_all = df_sim_protect_all.groupby(
                ['HS-nffmean']).sum()
            df_sim_protect_all = df_sim_protect_all.reset_index()  # gestion de la legende
            if var == 'Mercia':
                label = 'Somme des protections efficacy'
            else:
                label = '_nolegend_'
            # plot
            plt.plot(df_sim_protect_all['HS-nffmean'],
                     df_sim_protect_all['protection_efficacy'], '-' + color,
                     label=label)
            # fleche hsT2 - nff moyen
            hs_new = hsT2 - hs_moyen
            plt.annotate('', xy=(hs_new, 0), xytext=(hs_new, -0.03),
                         arrowprops=dict(color=color, arrowstyle="->",
                                         connectionstyle="arc3"))
        # mise en forme
        plt.xlabel('haun stage - nff moyen');
        plt.ylabel('global_efficacy')
        plt.grid(False)  # enlever la grille
        # plt.legend(numpoints=1, bbox_to_anchor=(1.1, 1.1), prop={'size': 9})

    return df_sim


# plot data publi contre simulation au nff moyen ou HS   (comp base sim_obs T1_T2)
def simulation_nffmoyen(var='Mercia', n_sim=5, n_plt=200, hs=12,
                        axis='MS'):  # axis = MS or all
    if var == 'Tremie12':
        nplants_lst = 30
    else:
        nplants_lst = n_plt

    x = 1
    df_sim = pandas.DataFrame()
    while x <= n_sim:
        npl, dfmoy, dfsd = treatment(name=var, sim=str(hs), nplants=n_plt,
                                     axis=axis, to_csv=False)
        print('var = ' + var + ' stade = ' + str(hs) + ' - nplants = ' + str(
            npl))
        dfmoy['var'] = var
        df_sim = df_sim.append(dfmoy)
        x += 1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    return df_sim_gr


def plot_data_simnffmoyen(varieties=['Mercia', 'Rht3', 'Tremie12', 'Tremie13'],
                          n_sim=1, n_plt=30, axis='MS', sim_date='T2',
                          decobs=['Mercia', 'Rht3', 'Tremie12']):
    '''Plot1 = obs
    Plot2 = obs/area
    sim_date is one of T1 T2 delta'''
    df_sim = pandas.DataFrame()
    # obs + preparation de obs
    df_obs = idata.dye_interception(decobs)
    df_obs.rename(columns={'name': 'var'}, inplace=True)
    df_obs.rename(columns={'N feuille': 'ntop_cur'}, inplace=True)
    df_obs.rename(columns={'HS': 'HS_obs'}, inplace=True)
    # fichier de simulation
    for var in varieties:

        # if ajust=False, attention au pb de decalage de feuille
        if sim_date == 'T1':
            date, hsr = HS_applications[var]['T1']
            df_obs = df_obs[df_obs['treatment'] == 'T1']
        elif sim_date == 'T2':
            date, hsr = HS_applications[var]['T2']
            df_obs = df_obs[df_obs['treatment'] == 'T2']

            # sim
        if sim_date == 'delta':
            # nff moyen
            if var == 'Mercia':
                delta = -0.7
            elif var == 'Rht3':
                delta = -1
            elif var == 'Tremie12':
                delta = -0.5
            elif var == 'Tremie13':
                delta = -0.4
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            hs_moyen = hs_moyen + delta
            df_sim_var = simulation_nffmoyen(var=var, hs=hs_moyen, n_sim=n_sim,
                                             n_plt=n_plt, axis=axis)
            df_obs = df_obs[df_obs['treatment'] == 'T2']
        else:
            df_sim_var = simulation_nffmoyen(var=var, hs=hsr, n_sim=n_sim,
                                             n_plt=n_plt, axis=axis)
        df_sim_var['deposit_Tartrazine/area'] = df_sim_var[
                                                    'deposit_Tartrazine'] / \
                                                df_sim_var['area']
        # df_sim_var['HS-nffmean'] = df_sim_var['HS'] - hs_moyen
        df_sim = df_sim.append(df_sim_var)
    # merge des 2 tableaux obs et sim
    df_all = df_sim.merge(df_obs)
    # column tartrazine/area
    df_all['meanObs/area'] = df_all['mean'] / df_all['area']

    # plot
    bar_width = 0.3;
    opacity = 0.4
    val = [[0, 'mean', 'deposit_Tartrazine'],
           [1, 'meanObs/area', 'deposit_Tartrazine/area']]
    fig, axes = plt.subplots(nrows=1, ncols=2)

    for x, col_data, col_sim in val:
        # on ne veut que les 4 feuilles du haut
        df_all = df_all[df_all['ntop_cur'] <= 4]

        n_groups = len(df_all)
        index = numpy.arange(n_groups)

        # publi
        if x == 0:
            rects1 = axes[x].bar(index, df_all[col_data], bar_width,
                                 alpha=opacity, color='y', yerr=df_all['IC'],
                                 error_kw=dict(ecolor='k'))
        else:
            rects1 = axes[x].bar(index, df_all[col_data], bar_width,
                                 alpha=opacity, color='y')
        # sim a hs nff moyen + delta
        rects2 = axes[x].bar(index + bar_width, df_all[col_sim], bar_width,
                             alpha=opacity,
                             color=['r', 'r', 'r', 'r', 'g', 'g', 'g', 'g', 'b',
                                    'b', 'b', 'b', 'm', 'm', 'm', 'm'])

        # Mise en forme
        df_all['HS'] = numpy.round(df_all['HS'], decimals=1)
        label = ['', 'Mercia', '', 'Rht3', '', 'Tremie12', '', 'Tremie13']
        axes[x].set_xticklabels(label, rotation=90, fontsize='small')
        if x == 0:
            axes[x].set_xlabel('var')
            axes[x].set_ylabel('deposit_tartrazine')
            axes[x].set_ylim(0, 8)
        else:
            axes[x].set_ylabel('deposit_tartrazine/area')
            axes[x].set_ylim(0, 0.5)

        fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)

        def autolabel(rects):
            # attach some text labels
            if rects == rects2:
                list = ['F1', 'F2', 'F3', 'F4']
            else:
                print('pb')
            for rect in rects:
                height = rect.get_height()
                if int(rect.get_x()) < 4:
                    if x == 0:
                        print(int(rect.get_x()))
                        axes[x].text(rect.get_x() + rect.get_width() / 2. - 0.1,
                                     1.05 * height + 0.2,
                                     list[int(rect.get_x())], ha='center',
                                     va='bottom', fontsize=8)

        autolabel(rects2)

        '''axes[x].text(0.2, 6.1, ''+str(date), bbox={'facecolor':'#FCF8F8', 'alpha':0.6, 'pad':10}, fontsize=12)

    axes[x].legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('F1', 'F2', 'F3', 'F4'), bbox_to_anchor=[1.10, 1.12], prop={'size':14})'''


# plot data publi contre simulation T2 et simulation au nff moyen  (data-sim T2 - sim delta)
def simulation_sim_simnffmoyen(var='Tremie12', n_sim=5, n_plt=200,
                               axis='MS'):  # axis = MS or all
    if var == 'Tremie12':
        nplants_lst = 30
        delta = -0.5
    else:
        nplants_lst = n_plt
        delta = -0.4

    # HS T2
    date, hsT2 = HS_applications[var]['T2']
    # nff moyen + delta
    adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
    df_phenT = adel.phenT()
    nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
    df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
    mean = df_phenT.mean()
    hs_moyen = mean['nff']
    hs_moyen = hs_moyen + delta
    # date a simuler
    lst = [[var, 'T2'], [var, str(hs_moyen)]]

    df_sim = pandas.DataFrame()
    x = 1
    while x <= n_sim:
        for var, stade in lst:
            npl, dfmoy, dfsd = treatment(name=var, sim=stade, nplants=n_plt,
                                         axis='MS', to_csv=False)
            print('var = ' + var + ' stade = ' + stade + ' - nbre de plantes sim = ' + str(
                npl))
            dfmoy['var'] = var
            df_sim = df_sim.append(dfmoy)
        x += 1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    return df_sim_gr


def plot_data_sim_simnffmoyen(varieties=['Tremie12', 'Tremie13'], n_sim=1,
                              n_plt=30, axis='MS', plot1=True, plot2=True):
    '''Plot1 = obs publi
    Plot2 = obs/area'''
    # obs + preparation de obs
    df_obs = idata.dye_interception()
    df_obs.rename(columns={'name': 'var'}, inplace=True)
    df_obs.rename(columns={'N feuille': 'ntop_cur'}, inplace=True)
    df_obs.rename(columns={'HS': 'HS_obs'}, inplace=True)
    # simulation
    df_sim = pandas.DataFrame()
    for var in varieties:
        df_sim_var = simulation_sim_simnffmoyen(var=var, n_sim=n_sim,
                                                n_plt=n_plt, axis=axis)
        df_sim_var['deposit_Tartrazine/area'] = df_sim_var[
                                                    'deposit_Tartrazine'] / \
                                                df_sim_var['area']
        df_sim = df_sim.append(df_sim_var)
    # merge des 2 tableaux obs et sim
    df_obs = df_obs[df_obs['treatment'] == 'T2']
    df_all = df_sim.merge(df_obs)
    # column tartrazine/area
    df_all['meanObs/area'] = df_all['mean'] / df_all['area']

    df_all['HS'] = numpy.round(df_all['HS'], decimals=2)

    if plot1 == True:
        for var in varieties:
            bar_width = 0.2;
            opacity = 0.4
            fig, ax0 = plt.subplots(nrows=1, ncols=1)
            # on ne veut que les 3 feuilles du haut
            df_all_var = df_all[df_all['var'] == var]
            df_all_var = df_all_var[df_all_var['ntop_cur'] <= 4]
            date, hsT2 = HS_applications[var]['T2']
            df_all_var_T2 = df_all_var[df_all_var['HS'] == hsT2]
            if var == 'Mercia':
                delta = -0.7
            elif var == 'Rht3':
                delta = -1
            elif var == 'Tremie12':
                delta = -0.5
            elif var == 'Tremie13':
                delta = -0.4
            # nff moyen + delta
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            hs_moyen = hs_moyen + delta
            df_all_var_delta = df_all_var[
                df_all_var['HS'] == round(hs_moyen, 2)]

            n_groups = len(df_all_var_T2)
            index = numpy.arange(n_groups)

            # publi
            rects1 = ax0.bar(index, df_all_var_T2['mean'], bar_width,
                             alpha=opacity, color='y', yerr=df_all_var_T2['IC'],
                             error_kw=dict(ecolor='k'))
            # sim T2
            rects2 = ax0.bar(index + bar_width,
                             df_all_var_T2['deposit_Tartrazine'], bar_width,
                             alpha=opacity, color='r')
            # sim T2+delta
            rects3 = ax0.bar(index + 2 * bar_width,
                             df_all_var_delta['deposit_Tartrazine'], bar_width,
                             alpha=opacity, color='b')

            # Mise en forme
            label = ['']
            ax0.set_xticklabels(label, rotation=90, fontsize='small')
            ax0.text(0.3, -0.3, '1', fontsize='small')
            ax0.text(1.3, -0.3, '2', fontsize='small')
            ax0.text(2.3, -0.3, '3', fontsize='small')
            ax0.text(3.3, -0.3, '4', fontsize='small')

            def autolabel(rects):
                # attach some text labels
                if rects == rects1:
                    list = ['publication', '', '', '']
                elif rects == rects2:
                    list = ['simulation T2', '', '', '']
                elif rects == rects3:
                    list = ['simulation nff moyen + delta', '', '', '']
                for rect in rects:
                    height = rect.get_height()
                    if rects == rects1:
                        ax0.text(rect.get_x() + rect.get_width() / 2.,
                                 1.05 * height + 0.7, list[int(rect.get_x())],
                                 ha='center', va='bottom', fontsize=8,
                                 rotation=90)
                    else:
                        ax0.text(rect.get_x() + rect.get_width() / 2.,
                                 1.05 * height + 0.2, list[int(rect.get_x())],
                                 ha='center', va='bottom', fontsize=8,
                                 rotation=90)

            autolabel(rects1)
            autolabel(rects2)
            autolabel(rects3)

            plt.ylim(ymax=7)
            # ax0.set_xlabel('ntop_cur')
            ax0.set_ylabel('deposit_tartrazine ' + var)

            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)

    if plot2 == True:
        for var in varieties:
            bar_width = 0.2;
            opacity = 0.4
            fig, ax1 = plt.subplots(nrows=1, ncols=1)
            # on ne veut que les 3 feuilles du haut
            df_all_var = df_all[df_all['var'] == var]
            df_all_var = df_all_var[df_all_var['ntop_cur'] <= 4]
            date, hsT2 = HS_applications[var]['T2']
            df_all_var_T2 = df_all_var[df_all_var['HS'] == hsT2]
            if var == 'Mercia':
                delta = -0.7
            elif var == 'Rht3':
                delta = -1
            elif var == 'Tremie12':
                delta = -0.5
            elif var == 'Tremie13':
                delta = -0.4
            # nff moyen + delta
            adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
            df_phenT = adel.phenT()
            nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
            df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
            mean = df_phenT.mean()
            hs_moyen = mean['nff']
            hs_moyen = hs_moyen + delta
            df_all_var_delta = df_all_var[
                df_all_var['HS'] == round(hs_moyen, 2)]

            n_groups = len(df_all_var_T2)
            index = numpy.arange(n_groups)

            # publi tartrazine/area
            rects1 = ax1.bar(index, df_all_var_T2['meanObs/area'], bar_width,
                             alpha=opacity, color='y')
            # sim T2
            rects2 = ax1.bar(index + bar_width,
                             df_all_var_T2['deposit_Tartrazine/area'],
                             bar_width, alpha=opacity, color='r')
            # sim T2+delta
            rects3 = ax1.bar(index + 2 * bar_width,
                             df_all_var_delta['deposit_Tartrazine/area'],
                             bar_width, alpha=opacity, color='b')

            # Mise en forme
            label = ['']
            ax1.set_xticklabels(label, rotation=90, fontsize='small')
            ax1.text(0.3, -0.01, '1', fontsize='small')
            ax1.text(1.3, -0.01, '2', fontsize='small')
            ax1.text(2.3, -0.01, '3', fontsize='small')
            ax1.text(3.3, -0.01, '4', fontsize='small')

            def autolabel(rects):
                # attach some text labels
                if rects == rects1:
                    list = ['publication', '', '', '']
                elif rects == rects2:
                    list = ['simulation T2', '', '', '']
                elif rects == rects3:
                    list = ['simulation nff moyen + delta', '', '', '']
                for rect in rects:
                    height = rect.get_height()
                    print(height)
                    ax1.text(rect.get_x() + rect.get_width() / 2.,
                             1.05 * height, list[int(rect.get_x())],
                             ha='center', va='bottom', fontsize=8, rotation=90)

            autolabel(rects1)
            autolabel(rects2)
            autolabel(rects3)

            plt.ylim(ymax=0.28)
            # ax1.set_xlabel('ntop_cur')
            ax1.set_ylabel('deposit_tartrazine/area ' + var)

            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)


# plot courbure
def simulation_courbure(var='Tremie12', n_sim=1, n_plt=30,
                        axis='MS'):  # axis = MS or all
    if var == 'Tremie12':
        nplants_lst = 30
        delta = -0.5
    else:
        nplants_lst = n_plt
        delta = -0.4

    # HS T2
    date, hsT2 = HS_applications[var]['T2']
    # nff moyen + delta
    adel = Reconstructions.get_reconstruction(name=var, nplants=n_plt)
    df_phenT = adel.phenT()
    nffs = df_phenT.set_index('plant').groupby(level=0).count()['n']
    df_phenT['nff'] = [nffs[v] for v in df_phenT['plant']]
    mean = df_phenT.mean()
    hs_moyen = mean['nff']
    hs_moyen = hs_moyen + delta
    # date a simuler
    lst = [[var, str(hs_moyen)]]

    df_sim = pandas.DataFrame()
    x = 1
    while x <= n_sim:
        for var, stade in lst:
            for nlim in [1, 2, 3, 4, 5, 6]:
                if var == 'Tremie12':
                    npl, dfmoy, dfsd = treatment(name=var, sim=stade,
                                                 nplants=n_plt, axis='MS',
                                                 to_csv=False,
                                                 nlim_Tremie12=nlim)
                elif var == 'Tremie13':
                    npl, dfmoy, dfsd = treatment(name=var, sim=stade,
                                                 nplants=n_plt, axis='MS',
                                                 to_csv=False,
                                                 nlim_Tremie13=nlim)
                print('var = ' + var + ' stade = ' + stade + ' - nbre de plantes sim = ' + str(
                    npl))
                dfmoy['var'] = var
                dfmoy['nlim'] = nlim
                df_sim = df_sim.append(dfmoy)
        x += 1
    df_sim_gr = df_sim.groupby(['var', 'HS', 'nlim', 'ntop_cur']).mean()
    df_sim_gr = df_sim_gr.reset_index()
    return df_sim_gr


def plot_courbure(varieties=['Tremie12', 'Tremie13'], n_sim=1, n_plt=30,
                  axis='MS', plot1=True, plot2=True):
    '''Plot1 = obs publi
    Plot2 = obs/area'''
    # simulation
    df_sim = pandas.DataFrame()
    for var in varieties:
        df_sim_var = simulation_courbure(var=var, n_sim=n_sim, n_plt=n_plt,
                                         axis=axis)
        df_sim_var['deposit_Tartrazine/area'] = df_sim_var[
                                                    'deposit_Tartrazine'] / \
                                                df_sim_var['area']
        df_sim = df_sim.append(df_sim_var)
    df_all = df_sim
    df_all['HS'] = numpy.round(df_all['HS'], decimals=2)

    if plot1 == True:
        for var in varieties:
            bar_width = 0.1;
            opacity = 0.4
            fig, ax0 = plt.subplots(nrows=1, ncols=1)
            # on ne veut que les 3 feuilles du haut
            df_all_var = df_all[df_all['var'] == var]
            df_all_var = df_all_var[df_all_var['ntop_cur'] <= 3]

            df_all_var_nlim1 = df_all_var[df_all_var['nlim'] == 1]
            df_all_var_nlim2 = df_all_var[df_all_var['nlim'] == 2]
            df_all_var_nlim3 = df_all_var[df_all_var['nlim'] == 3]
            df_all_var_nlim4 = df_all_var[df_all_var['nlim'] == 4]
            df_all_var_nlim5 = df_all_var[df_all_var['nlim'] == 5]
            df_all_var_nlim6 = df_all_var[df_all_var['nlim'] == 6]

            n_groups = len(df_all_var_nlim1)
            index = numpy.arange(n_groups)

            rects1 = ax0.bar(index, df_all_var_nlim1['deposit_Tartrazine'],
                             bar_width, alpha=opacity, color='y')
            rects2 = ax0.bar(index + bar_width,
                             df_all_var_nlim2['deposit_Tartrazine'], bar_width,
                             alpha=opacity, color='r')
            rects3 = ax0.bar(index + 2 * bar_width,
                             df_all_var_nlim3['deposit_Tartrazine'], bar_width,
                             alpha=opacity, color='b')
            rects4 = ax0.bar(index + 3 * bar_width,
                             df_all_var_nlim4['deposit_Tartrazine'], bar_width,
                             alpha=opacity, color='g')
            rects5 = ax0.bar(index + 4 * bar_width,
                             df_all_var_nlim5['deposit_Tartrazine'], bar_width,
                             alpha=opacity, color='m')
            rects6 = ax0.bar(index + 5 * bar_width,
                             df_all_var_nlim6['deposit_Tartrazine'], bar_width,
                             alpha=opacity, color='c')

            # Mise en forme
            label = ['']
            ax0.set_xticklabels(label, rotation=90, fontsize='small')
            ax0.text(0.3, -0.4, '1', fontsize='small')
            ax0.text(1.3, -0.4, '2', fontsize='small')
            ax0.text(2.3, -0.4, '3', fontsize='small')

            def autolabel(rects):
                # attach some text labels
                if rects == rects1:
                    list = ['nlim=1', '', '']
                elif rects == rects2:
                    list = ['nlim=2', '', '']
                elif rects == rects3:
                    list = ['nlim=3', '', '']
                elif rects == rects4:
                    list = ['nlim=4', '', '']
                elif rects == rects5:
                    list = ['nlim=5', '', '']
                elif rects == rects6:
                    list = ['nlim=6', '', '']
                for rect in rects:
                    height = rect.get_height()
                    if rects == rects1:
                        ax0.text(rect.get_x() + rect.get_width() / 2.,
                                 1.05 * height + 0.2, list[int(rect.get_x())],
                                 ha='center', va='bottom', fontsize=8,
                                 rotation=90)
                    else:
                        ax0.text(rect.get_x() + rect.get_width() / 2.,
                                 1.05 * height + 0.2, list[int(rect.get_x())],
                                 ha='center', va='bottom', fontsize=8,
                                 rotation=90)

            autolabel(rects1)
            autolabel(rects2)
            autolabel(rects3)
            autolabel(rects4)
            autolabel(rects5)
            autolabel(rects6)

            plt.ylim(ymax=12)
            # ax0.set_xlabel('ntop_cur')
            ax0.set_ylabel('deposit_tartrazine ' + var)

            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)

    if plot2 == True:
        for var in varieties:
            bar_width = 0.1;
            opacity = 0.4
            fig, ax0 = plt.subplots(nrows=1, ncols=1)
            # on ne veut que les 3 feuilles du haut
            df_all_var = df_all[df_all['var'] == var]
            df_all_var = df_all_var[df_all_var['ntop_cur'] <= 3]

            df_all_var_nlim1 = df_all_var[df_all_var['nlim'] == 1]
            df_all_var_nlim2 = df_all_var[df_all_var['nlim'] == 2]
            df_all_var_nlim3 = df_all_var[df_all_var['nlim'] == 3]
            df_all_var_nlim4 = df_all_var[df_all_var['nlim'] == 4]
            df_all_var_nlim5 = df_all_var[df_all_var['nlim'] == 5]
            df_all_var_nlim6 = df_all_var[df_all_var['nlim'] == 6]

            n_groups = len(df_all_var_nlim1)
            index = numpy.arange(n_groups)

            rects1 = ax0.bar(index, df_all_var_nlim1['deposit_Tartrazine/area'],
                             bar_width, alpha=opacity, color='y')
            rects2 = ax0.bar(index + bar_width,
                             df_all_var_nlim2['deposit_Tartrazine/area'],
                             bar_width, alpha=opacity, color='r')
            rects3 = ax0.bar(index + 2 * bar_width,
                             df_all_var_nlim3['deposit_Tartrazine/area'],
                             bar_width, alpha=opacity, color='b')
            rects4 = ax0.bar(index + 3 * bar_width,
                             df_all_var_nlim4['deposit_Tartrazine/area'],
                             bar_width, alpha=opacity, color='g')
            rects5 = ax0.bar(index + 4 * bar_width,
                             df_all_var_nlim5['deposit_Tartrazine/area'],
                             bar_width, alpha=opacity, color='m')
            rects6 = ax0.bar(index + 5 * bar_width,
                             df_all_var_nlim6['deposit_Tartrazine/area'],
                             bar_width, alpha=opacity, color='c')

            # Mise en forme
            label = ['']
            ax0.set_xticklabels(label, rotation=90, fontsize='small')
            ax0.text(0.3, -0.02, '1', fontsize='small')
            ax0.text(1.3, -0.02, '2', fontsize='small')
            ax0.text(2.3, -0.02, '3', fontsize='small')

            def autolabel(rects):
                # attach some text labels
                if rects == rects1:
                    list = ['nlim=1', '', '']
                elif rects == rects2:
                    list = ['nlim=2', '', '']
                elif rects == rects3:
                    list = ['nlim=3', '', '']
                elif rects == rects4:
                    list = ['nlim=4', '', '']
                elif rects == rects5:
                    list = ['nlim=5', '', '']
                elif rects == rects6:
                    list = ['nlim=6', '', '']
                for rect in rects:
                    height = rect.get_height()
                    if rects == rects1:
                        ax0.text(rect.get_x() + rect.get_width() / 2.,
                                 1.05 * height + 0.02, list[int(rect.get_x())],
                                 ha='center', va='bottom', fontsize=8,
                                 rotation=90)
                    else:
                        ax0.text(rect.get_x() + rect.get_width() / 2.,
                                 1.05 * height + 0.02, list[int(rect.get_x())],
                                 ha='center', va='bottom', fontsize=8,
                                 rotation=90)

            autolabel(rects1)
            autolabel(rects2)
            autolabel(rects3)
            autolabel(rects4)
            autolabel(rects5)
            autolabel(rects6)

            plt.ylim(ymax=0.5)
            # ax0.set_xlabel('ntop_cur')
            ax0.set_ylabel('deposit_tartrazine/area ' + var)

            fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
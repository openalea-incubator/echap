import pandas 
from new_annual_loop_pest import pesticide_loop, repartition_at_application,repartition_at_applicationArch,get_reconstruction
import pesticide_protocol as proto
import numpy as np
import pylab as p
import matplotlib.pyplot as plt
from alinea.echap.weather_data import *

from itertools import cycle, islice
from pylab import *

# import pandas
# from openalea.deploy.shared_data import shared_data

# from alinea.astk.Weather import Weather
# from alinea.astk.TimeControl import *
# from alinea.caribu.caribu_star import rain_and_light_star
# from alinea.alep.protocol import infect
# from alinea.echap.interception_leaf import pesticide_applications
# from alinea.echap.microclimate_leaf import microclimate_leaf
# from alinea.echap.recorder import LeafElementRecorder
# from alinea.echap.interfaces import record as do_record #to avoid confusion with numpy record

# from macros_annual_loop import setup_canopy, update_lesions, update_pesticides, dispersion, contamination, pesticide_intercept
# import alinea.septo3d

# from alinea.echap.tests_nodes import plot_pesticide

#g, recorder = pesticide_loop(meteo_file='meteo00-01.txt', applications=proto.Mercia2011_treat1, start="2011-04-11", periods=24, freq='H', TB=0, delayH=1, delayDD=15)
# recorder.save_records('test.csv')   
# recorder.save_plot(what='area',t = 'dd',  by='ntop', axe = 'MS', plant='all', fig_name='test.png')   


# test
"""
g, recorder = pesticide_loop(meteo_file='meteo00-01.txt', start="2001-04-25", periods=50, freq='H', TB=0, delayH=1, delayDD=15, applications=proto.Test2001_treat1)
recorder.save_records('testMAI2014.csv')   
recorder.save_plot(what='surfacic_doses_Epoxiconazole', t = 'iter',  by='ntop', axe = 'MS', plant='all', fig_name='test2.png')   
recorder.plt.show()
"""

# Appel a la fonction repartition_at_application
# df = repartition_at_application(appdate = '2011-04-19', dose = 0.5, age = 1166)
# print 'df.columns after = ', df.columns

# -----------------------------------------------------------------------------------------------------------------------------------
# VERSION SIMPLE (prend une date TT en entree)
def g_constr(name='Mercia',nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0):

    # Appel a la fonction get_reconstruction
    pgen,adel,domain, domain_area, convUnit, nplants = get_reconstruction(name=name, nplants=nplants, nsect=nsect, seed=seed, sample=sample, as_pgen=as_pgen, dTT_stop=dTT_stop)
    # Appel a la fonction repartition_at_application lié à l'architecture
    #age=1166
    age = 1500
    g = adel.setup_canopy(age)
    #df = repartition_at_applicationArch(appdate = '2011-04-19', dose = 40, g=g)
    df = repartition_at_applicationArch(appdate = '2011-05-11', dose = 40, g=g)
    return df

def treatment(name='Mercia',nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0): 

    df = g_constr(name, nplants, nsect, seed, sample, as_pgen, dTT_stop)

    gr=df.groupby(['plant', 'date', 'axe'], group_keys=False)
    
    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        return sub
        
    data = gr.apply(_fun)    
    data['ntop_cur'] = data['n_max'] - data['metamer'] + 1    
    data = data.convert_objects() # conversion des objets en types numpy
    
    #test groupby by multiple apres avoir fait calcul sur data
    """data = data[data.ntop_cur <= 4]
    f = {'surfacic_doses_Epoxiconazole':['sum','mean']}
    dftest = fdf.groupby(['plant', 'ntop_cur', 'date', 'axe']).agg(f)
    """

    data = data[data.ntop_cur <= 4]
    gr=data.groupby(['plant', 'ntop_cur', 'date', 'axe'], as_index=False)
    dmp=gr.sum() # dmp contient, pour chaque quadruplet ('plant', 'ntop_cur', 'date', 'axe'), 
                 # la somme des area, green_area, id, length, metamer, ntop, penetrated_doses_Epoxiconazole, 
                 # senesced_area, surfacic_doses_Epoxiconazole, n_max
     
    gr = dmp.groupby(['ntop_cur','plant','date'],as_index=False)
    # calcul, pour chaque triplet ('ntop_cur','plant','date'), de la moyenne et de l ecart type de surfacic_doses_Epoxiconazole
    grtest = gr['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std' : np.std})

    gr2 = dmp.groupby(['ntop_cur','date'],as_index=False)
    # calcul, pour chaque doublon ('ntop_cur','date'), de la moyenne et de l ecart type de surfacic_doses_Epoxiconazole
    grtest4 = gr2['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std' : np.std})
    
    return grtest4

#plot des donnees obs
def plot(name='Mercia'):

    if name is 'Mercia' or 'Rht3':
        data_file = 'ObsData2011.csv'

    df_obs = pandas.read_csv(data_file, sep=';')
    df_obs = df_obs[df_obs['feuille']<5]
    df_obs['var_stade'] = df_obs['variete'] + "_" + df_obs['stade']
    
    cols = df_obs.var_stade.value_counts().shape[0]
    fig, axes = plt.subplots(1, cols, figsize=(8, 8))

    for x, var_stade in enumerate(df_obs.var_stade.value_counts().index.values):
        data = df_obs[(df_obs['var_stade'] == var_stade)]
        data = data.groupby(['dose','feuille']).colorant.sum()
        print (data)
        print type(data.index)
        left = [k[0] for k in enumerate(data)]
        right = [k[1] for k in enumerate(data)]

        #gestion des 4 couleurs differentes des barplot
        my_colors = list(islice(cycle(['b', 'r', 'g', 'y']), None, 4))
        axes[x].bar(left,right,label="%s" % (var_stade), color=my_colors)
        axes[x].set_xticks(left, minor=False)
        # modification de la legende axe des abcisses
        data = data.reset_index()
        axes[x].set_xticklabels(data['dose'].values)
        #axe des ordonnees
        axes[x].set_ylim(0, 25)

        axes[x].legend(loc='best')
        axes[x].grid(True)
        #changement couleur autour graph
        fig.patch.set_facecolor('#FCF8F8')
        #changement couleur fond graph
        axes[x].patch.set_facecolor('#D1D1D1')
        #titre
        fig.suptitle('Dose par Feuille par Stade par Variete', fontsize=15)

    axes[x].text(12, 20, 'feuille 1', style='italic', bbox={'facecolor':'b', 'alpha':0.5, 'pad':10})
    axes[x].text(12, 19, 'feuille 2', style='italic', bbox={'facecolor':'r', 'alpha':0.5, 'pad':10})
    axes[x].text(12, 18, 'feuille 3', style='italic', bbox={'facecolor':'g', 'alpha':0.5, 'pad':10})
    axes[x].text(12, 17, 'feuille 4', style='italic', bbox={'facecolor':'y', 'alpha':0.5, 'pad':10})
    plt.show()

# -----------------------------------------------------------------------------------------------------------------------------------
# IDEM mais gestion d'un range de date  
def g_constr_age(name='Mercia',nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, dose=40, age=1166):

    # Conversion age (TT) en date ('aaaa-mm-jj')
    from alinea.astk.TimeControl import thermal_time
    if name is 'Mercia' or 'Rht3':
        data_file = 'METEO_stationINRA_20102011.csv'
        # correspondance entre date et TT
        met = Boigneville_2010_2011()
        seq = pandas.date_range(start="2010-10-15", end="2011-06-20", freq='H')
        bid = thermal_time(seq, met.data)
        # groupe par date
        seq = pandas.date_range(start="2010-10-15", end="2011-06-20")
        bid = bid[seq]
        #renommer colonnes
        bid = bid.reset_index()
        bid.columns = ['datetime','TT']
        bid['datetime'] = bid['datetime'].map(lambda x: x.strftime('%Y-%m-%d'))
        bid['TT'] = bid['TT'].astype(int)
        # trouver la correspondance en notre age TT et la date dont a besoin repartition_at_applicationArch()
        da=0; appdate=0
        while appdate==0:
            if da<len(bid):
                if bid['TT'][da]==age:
                    appdate=bid['datetime'][da]
                else:
                    da = da+1
            if da==len(bid):
                da=0; age=age+1
    # Appel a la fonction get_reconstruction
    pgen,adel,domain, domain_area, convUnit, nplants = get_reconstruction(name=name, nplants=nplants, nsect=nsect, seed=seed, sample=sample, as_pgen=as_pgen, dTT_stop=dTT_stop)
    # Appel a la fonction repartition_at_application lié à l'architecture
    g = adel.setup_canopy(age)
    df = repartition_at_applicationArch(appdate = appdate, dose = 40, g = g)
    return df

def treatment_age(name='Mercia',nplants=30, nsect=3, seed=1, sample='sequence', as_pgen=False, dTT_stop=0, dose=40): 

    dd = range(1000,1400,300)
    #df_sim = [g_constr_age(name, nplants, nsect, seed, sample, as_pgen, dTT_stop, dose, age) for age in dd]
    
    # lancer la reconstruction pour chaque age du range et concatener les dataframe de sortie
    dl = len(dd); dlc = 1
    df = g_constr_age(name, nplants, nsect, seed, sample, as_pgen, dTT_stop, dose, dd[0])
    while dlc < dl :
        df_other = g_constr_age(name, nplants, nsect, seed, sample, as_pgen, dTT_stop, dose, dd[dlc])
        dlc = dlc+1
        df = df.append(df_other, ignore_index=True)

    gr = df.groupby(['plant', 'date', 'axe'], group_keys=False)
    
    def _fun(sub):
        sub['n_max'] = sub['metamer'].max()
        return sub
    
    data = gr.apply(_fun)    
    data['ntop_cur'] = data['n_max'] - data['metamer'] + 1    
    data = data.convert_objects() # conversion des objets en types numpy
    
    #test groupby by multiple apres avoir fait calcul sur data
    """data = data[data.ntop_cur <= 4]
    f = {'surfacic_doses_Epoxiconazole':['sum','mean']}
    dftest = fdf.groupby(['plant', 'ntop_cur', 'date', 'axe']).agg(f)
    """

    data = data[data.ntop_cur <= 4]
    gr=data.groupby(['plant', 'ntop_cur', 'date', 'axe'], as_index=False)
    dmp=gr.sum() # dmp contient, pour chaque quadruplet ('plant', 'ntop_cur', 'date', 'axe'), 
                 # la somme des area, green_area, id, length, metamer, ntop, penetrated_doses_Epoxiconazole, 
                 # senesced_area, surfacic_doses_Epoxiconazole, n_max
     
    gr = dmp.groupby(['ntop_cur','plant','date'],as_index=False)
    # calcul, pour chaque triplet ('ntop_cur','plant','date'), de la moyenne et de l ecart type de surfacic_doses_Epoxiconazole
    grtest = gr['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std' : np.std})

    gr2 = dmp.groupby(['ntop_cur','date'],as_index=False)
    # calcul, pour chaque doublon ('ntop_cur','date'), de la moyenne et de l ecart type de surfacic_doses_Epoxiconazole
    grtest4 = gr2['surfacic_doses_Epoxiconazole'].agg({'surfacic_doses_Epoxiconazole_mean' : np.mean, 'surfacic_doses_Epoxiconazole_std' : np.std})

    return grtest4   											 
	
#grtest = dmp.groupby('ntop_cur','date').apply(lambda x: 
#             Series({'r': (x.y + x.z).sum() / x.z.sum(), 
#                     's': (x.y + x.z ** 2).sum() / x.z.sum()}))

#gr=dmp.groupby(['ntop_cur','date'])
#dms=gr.agg(numpy.mean)


#NG ajout ecart type
#dmserr=gr.agg(numpy.std)
#dms=dms.reset_index()
#dmserr=dmserr.reset_index()

#plt.figure()
#plt = grtest.plot('ntop_cur','surfacic_doses_Epoxiconazole_mean',kind="bar")
#plt.show()

'''mise en forme sous format barplot
fig = p.figure()
ax = fig.add_subplot(1,1,1)
#N = len(ntop_cur)

#err = [1.2, 1.5, 2.5, 1.2]
#ax.bar('ntop_cur', 'n_max', facecolor='#777777',
#       align='center', yerr=err, ecolor='black')
ax.bar(dms.ntop_cur, dms.n_max, facecolor='#777777',
       align='center', ecolor='black')
#Create a y label
ax.set_ylabel('Intercepted dose (?/?)')
 
# Create a title, in italics
ax.set_title('Intercepted dose, by ntop_cur',fontstyle='italic')
 
# This sets the ticks on the x axis to be exactly where we put
# the center of the bars.
ax.set_xticks(dms.ntop_cur)
 
# Labels for the ticks on the x axis.  It needs to be the same length
# as y (one label for each bar)
group_labels = ['1', '2',
                 '3', '4']
 
# Set the x tick labels to the group_labels defined above.
ax.set_xticklabels(group_labels)
 
# Extremely nice function to auto-rotate the x axis labels.
# It was made for dates (hence the name) but it works
# for any long x tick labels
fig.autofmt_xdate()
 
p.show()

#
# OK autre solution pour les barplots
#

#definition du nombre de groupe : correspond ici à la dose de pesticide epandu
n_groups = 3

fig, ax = plt.subplots()

#index = np.arange(n_groups)
bar_width = 1

opacity = 0.4
error_config = {'ecolor': '0.3'}

#for i in range (0,n_groups):
rects1 = plt.bar(grtest2.ntop_cur[0], grtest2.surfacic_doses_Epoxiconazole_mean[0], bar_width,
                 alpha=opacity,
                 color='b',
                 yerr=grtest2.surfacic_doses_Epoxiconazole_std[0],
                 error_kw=error_config,
                 label='leaf 1')




                 
plt.xlabel('Dose 0.5')
plt.ylabel('Quantite retenue')
plt.title('Quantite de pesticide par feuille')
plt.xticks(grtest2.ntop_cur, (' '))
plt.legend()

plt.tight_layout()
plt.show()

#solution 3 : probleme avec non reconnaissance DataFrame

#df2 = DataFrame(dms.n_max[1], columns=['1', '2', '3', '4'])
plt.figure();
#grplot = DataFrame(grtest2.surfacic_doses_Epoxiconazole_mean,columns=['1','2','3','4'])
grplot = grtest2.surfacic_doses_Epoxiconazole_mean
grplot.plot(kind='bar',yerr=grtest2.surfacic_doses_Epoxiconazole_std,grid=False)
plt.show()

"""autre solution"""
fig = plt.figure()
ax1 = plt.subplot(111,ylabel='A')
#ax2 = gcf().add_axes(ax1.get_position(), sharex=ax1, frameon=False, ylabel='axes2')
#ax2 =ax1.twinx()
#ax2.set_ylabel('B')
ax1.bar(grtest2.ntop_cur.values,grtest2.surfacic_doses_Epoxiconazole_mean.values, width =0.4, color ='g', align = 'center')
#ax2.bar(drtest2.index,df.B.values, width = 0.4, color='r', align = 'edge')
ax1.legend(['A'], loc = 'upper left')
#ax2.legend(['B'], loc = 'upper right')
fig.show()

N = 1
ind = np.arange(N)  # the x locations for the groups dans echap = nombre de doses = N
width = 0.2    # the width of the bars

fig = plt.figure()
ax = fig.add_subplot(111)

yvals = [grtest2.surfacic_doses_Epoxiconazole_mean[0]]
rects1 = ax.bar(ind, yvals, width, color='r')
zvals = [grtest2.surfacic_doses_Epoxiconazole_mean[1]]
rects2 = ax.bar(ind+width, zvals, width, color='g')
kvals = [grtest2.surfacic_doses_Epoxiconazole_mean[2]]
rects3 = ax.bar(ind+width*2, kvals, width, color='b')
mvals =  [grtest2.surfacic_doses_Epoxiconazole_mean[3]]
rects4 = ax.bar(ind+width*3, mvals, width, color='y',yerr=grtest2.surfacic_doses_Epoxiconazole_std[3])

ax.set_ylabel('Depot pesticides par etage foliaire')
ax.set_xticks(ind+width)
ax.set_xticklabels( ('etage 1', 'etage 2', 'etage 3','etage 4') )
#ax.legend( (rects1[0], rects2[0], rects3[0],rects4[0]), ('1', '2', '3','4') )
ax.legend( (yvals, zvals, kvals,mvals), ('1', '2', '3','4') )

def autolabel(rects):
    for rect in rects:
        h = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%int(h),
                ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
autolabel(rects4)

plt.show()

'''
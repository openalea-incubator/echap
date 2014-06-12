from new_annual_loop_pest import pesticide_loop, repartition_at_application
import pesticide_protocol as proto
import numpy as np
import pylab as p
import matplotlib.pyplot as plt

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
g, recorder = pesticide_loop(meteo_file='meteo00-01.txt', start="2001-04-25", periods=50, freq='H', TB=0, delayH=1, delayDD=15, applications=proto.Test2001_treat1)
recorder.save_records('testMAI2014.csv')   
recorder.save_plot(what='surfacic_doses_Epoxiconazole', t = 'iter',  by='ntop', axe = 'MS', plant='all', fig_name='test2.png')   
recorder.plt.show()




# Function repartition_at_application : appeler ensuite ici pour test appel
"""appdate = '2011-04-19' ; dose = 0.5 ; age = 1166
print '\n\nrepartition_at_application 2 !!\n\n'
from macros_annual_loop import setup_canopy
from alinea.echap.recorder import LeafElementRecorder
recorder = LeafElementRecorder()
g, adel, domain, domain_area, convUnit, nplants = setup_canopy(age=1166)
applications= """date,dose, product_name
%s 10:00:00, %f, Opus
"""%(appdate, dose)
application_data = pesticide_applications(applications)
g,_=pesticide_intercept(g, application_data)
do_record(g, application_data, recorder)
df_1 =  recorder.get_records()
print 'df_1.columns = ', df_1.columns
"""

# Appel a la fonction repartition_at_application
df = repartition_at_application(appdate = '2011-04-19', dose = 1.0, age = 1166)
print 'df.columns after = ', df.columns


from itertools import cycle, islice
from pylab import *


gr=df.groupby(['plant', 'date', 'axe'], group_keys=False)
def _fun(sub):
    sub['n_max'] = sub['metamer'].max()
    return sub
data = gr.apply(_fun)    
data['ntop_cur'] = data['n_max'] - data['metamer'] + 1    



data.to_csv('data.csv', index=False) # cc

data = pandas.read_csv('data.csv')

#test groupby by multiple apres avoir fait calcul sur data
"""data = data[data.ntop_cur <= 4]
f = {'surfacic_doses_Epoxiconazole':['sum','mean']}
dftest = fdf.groupby(['plant', 'ntop_cur', 'date', 'axe']).agg(f)
"""

data = data[data.ntop_cur <= 4]
gr=data.groupby(['plant', 'ntop_cur', 'date', 'axe'], as_index=False)
# dmp=gr.agg([numpy.sum])
dmp=gr.sum()
dmp.to_csv('dmpSum.csv') # cc


#test apply pour effectuer mean et std
grtest =dmp.groupby(['ntop_cur','plant','date'],as_index=False).apply(lambda subf: 
	[subf['surfacic_doses_Epoxiconazole'].mean()])
	
grtest = dmp.groupby(['ntop_cur','plant','date'],as_index=False)	
surfacic_doses_Epoxiconazole_mean = grtest['surfacic_doses_Epoxiconazole'].mean()
surfacic_doses_Epoxiconazole_std = grtest['surfacic_doses_Epoxiconazole'].std()
	
#grtest = dmp.groupby('ntop_cur','date').apply(lambda x: 
#             Series({'r': (x.y + x.z).sum() / x.z.sum(), 
#                     's': (x.y + x.z ** 2).sum() / x.z.sum()}))

grtestdf=pandas.DataFrame(grtest)
grtesterr =dmp.groupby(['ntop_cur','plant','date']).apply(lambda subf: 
	[subf['surfacic_doses_Epoxiconazole'].std()])

#gr=dmp.groupby(['ntop_cur','date'])
#dms=gr.agg(numpy.mean)


#NG ajout ecart type
#dmserr=gr.agg(numpy.std)
#dms=dms.reset_index()
#dmserr=dmserr.reset_index()

grtest.plot('ntop_cur','surfacic_doses_Epoxiconazole',kind="bar")

"""mise en forme sous format barplot
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
"""
#
# OK autre solution pour les barplots
#

#definition du nombre de groupe : correspond ici à la dose de pesticide epandu
#n_groups = 2

fig, ax = plt.subplots()

#index = np.arange(n_groups)
bar_width = 1

opacity = 0.4
#error_config = {'ecolor': '0.3'}

#for i in range (0,n_groups):
rects1 = plt.bar(d.ntop_cur[0], dms.n_max[0], bar_width,
                 alpha=opacity,
                 color='b',
                 yerr=dmserr.n_max[0],
#                error_kw=error_config,
                 label='leaf 1')

rects2 = plt.bar(dms.ntop_cur[5], dms.n_max[5], bar_width,
                 alpha=opacity,
                 color='r',
               yerr=dmserr.n_max[5],
 #                error_kw=error_config,
                 label='leaf 2')

rects3 = plt.bar(dms.ntop_cur[7], dms.n_max[7], bar_width,
                 alpha=opacity,
                 color='green',
                yerr=dmserr.n_max[7],
#                 error_kw=error_config,
                 label='leaf 3')


                 
plt.xlabel('Dose 0.5')
plt.ylabel('Quantite retenue')
plt.title('Quantite de pesticide par feuille')
plt.xticks(dms.ntop_cur, (' '))
plt.legend()

plt.tight_layout()
plt.show()

"""
"""
#solution 3 : probleme avec non reconnaissance DataFrame

#df2 = DataFrame(dms.n_max[1], columns=['1', '2', '3', '4'])
plt.figure();
dms.n_max.plot(kind='bar');
plt.show()
"""
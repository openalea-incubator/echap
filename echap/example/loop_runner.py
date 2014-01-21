from new_annual_loop_pest import pesticide_loop, repartition_at_application
import pesticide_protocol as proto
import numpy

#g, recorder = pesticide_loop(meteo_file='meteo00-01.txt', applications=proto.Mercia2011_treat1, start="2011-04-11", periods=24, freq='H', TB=0, delayH=1, delayDD=15)
# recorder.save_records('test.csv')   
# recorder.save_plot(what='area',t = 'dd',  by='ntop', axe = 'MS', plant='all', fig_name='test.png')   


# test
g, recorder = pesticide_loop(meteo_file='meteo00-01.txt', start="2001-04-25", periods=50, freq='H', TB=0, delayH=1, delayDD=15, applications=proto.Test2001_treat1)
#recorder.save_records('test.csv')   
#recorder.save_plot(what='area',t = 'iter',  by='ntop', axe = 'MS', plant='all', fig_name='test.png')   


df = repartition_at_application(appdate = '2011-04-19', dose = 0.5, age = 1166)
from itertools import cycle, islice
from pylab import *

gr=df.groupby(['plant', 'date', 'axe'], group_keys=False)
def _fun(sub):
    sub['n_max'] = sub['metamer'].max()
    return sub
data = gr.apply(_fun)    
data['ntop_cur'] = data['n_max'] - data['metamer'] + 1    

data = data[data.ntop_cur <= 4]
gr=data.groupby(['plant', 'ntop_cur', 'date', 'axe'])
dmp=gr.agg(numpy.sum)
dmp = dmp.reset_index()
gr=dmp.groupby(['ntop_cur','date'])
dms=gr.agg(numpy.mean)
dms=dms.reset_index()
dms.plot(x,y,kind="bar")


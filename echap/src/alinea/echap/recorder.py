import pandas
import numpy
from pylab import *

class LeafElementRecorder:
    def __init__(self):
        self.data = {}
        self.counts = 0
        
    def record(self, node, header):
        if node.area > 0:
            data = {}
            properties = node.properties()
            items = ['length', 'area', 'green_area', 'senesced_area']
            for item in items:
                if item in properties:
                    data[item] = properties[item]
                    
            items = ['surfacic_doses', 'penetrated_doses', 'global_efficacy']
            for item in items:
                if item in properties:
                    for compound in properties[item]:
                        data['_'.join([item,compound])] = properties[item][compound]
            
            if 'lesions' in properties:
                les = properties['lesions'][0]
                data.update(les.Lesions.EtatS())
           
            data.update(header)
            self.data[self.counts] = data
            self.counts += 1
        
    def save_records(self, path='./echap_outputs'):
        d = pandas.DataFrame(self.data)
        d.T.to_csv(path)
        
    def get_records(self):
        d = pandas.DataFrame(self.data)
        return d.T
        
    def plot(self, what='area',t = 'dd',  by='ntop', axe = 'MS', plant='all'):
        from itertools import cycle, islice
        d = self.get_records()
        d = d[d[what].notnull()]
        if plant == 'all':
            dm=d[d.axe==axe]
        else:
            dm=d[(d.axe==axe) & (d.plant==plant)]
        gr=dm.groupby(['plant', by, t])
        dmp=gr.agg(numpy.sum)
        dmp = dmp.reset_index()
        gr=dmp.groupby([by,t])
        dms=gr.agg(numpy.mean)
        dms=dms.reset_index()
        gr=dms.groupby(by)
        colors = list(islice(cycle(['k', 'r', 'g', 'b', 'y', 'c','m']),None,len(gr)))
        for i,group in enumerate(gr):
            group[1].plot(t,what,color=colors[i])
       

    def save_plot(self, what='area',t = 'dd',  by='ntop', axe = 'MS', plant='all', fig_name='test.png'):
        from itertools import cycle, islice
        d = self.get_records()
        d = d[d[what].notnull()]
        if plant == 'all':
            dm=d[d.axe==axe]
        else:
            dm=d[(d.axe==axe) & (d.plant==plant)]
        gr=dm.groupby(['plant', by, t])
        dmp=gr.agg(numpy.sum)
        dmp = dmp.reset_index()
        gr=dmp.groupby([by,t])
        dms=gr.agg(numpy.mean)
        dms=dms.reset_index()
        gr=dms.groupby(by)
        colors = list(islice(cycle(['k', 'r', 'g', 'b', 'y', 'c','m']),None,len(gr)))
        for i,group in enumerate(gr):
            df=group[1]
            plot(df[t].values, df[what].values, label=what)
        legend( loc='upper right', numpoints = 1 )
        savefig(fig_name)           
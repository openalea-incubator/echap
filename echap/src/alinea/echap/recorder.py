import pandas
import numpy

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
                    
            items = ['surfacic_doses', 'penetrated_doses']
            for item in items:
                if item in properties:
                    for compound in properties[item]:
                        data['_'.join([item,compound])] = properties[item][compound]
            
            if 'lesions' in properties:
                les = properties['lesions'][0]
                data.update(les.EtatS())
           
            data.update(header)
            self.data[self.counts] = data
            self.counts += 1
        
    def save_records(self, path='./echap_outputs'):
        d = pandas.DataFrame(self.data)
        d.T.to_csv(path)
        
    def get_records(self):
        d = pandas.DataFrame(self.data)
        return d.T
        
    def plotkin(self, what='area', n='ntop'):
        from itertools import cycle, islice
        d = self.get_records()
        dm=d[(d.axe=='MS') & (d.plant=='plant1')]
        gr=dm.groupby(['plant', n,'dd'])
        dmp=gr.agg(numpy.sum)
        dmp = dmp.reset_index()
        gr=dmp.groupby([n,'dd'])
        dms=gr.agg(numpy.mean)
        dms=dms.reset_index()
        gr=dms.groupby(n)
        colors = list(islice(cycle(['k', 'r', 'g', 'b', 'y', 'c','m']),None,len(gr)))
        for i,group in enumerate(gr):
            group[1].plot('dd',what,color=colors[i])
       
        
import pandas

class LeafElementRecorder:
    def __init__(self):
        self.data = {}
        self.counts = 0
        
    def record(self, node, header):
        if node.area > 0:
            data = {}
            properties = node.properties()
            items = ['area', 'green_area']
            for item in items:
                if item in properties:
                    data[item] = properties[item]
            items = ['surfacic_doses', 'penetrated_doses']
            for item in items:
                if item in properties:
                    for compound in properties[item]:
                        data['_'.join([item,compound])] = properties[item][compound]
           
            data.update(header)
            self.data[self.counts] = data
            self.counts += 1
        
    def save_records(self, path='./echap_outputs'):
        d = pandas.DataFrame(self.data)
        d.T.to_csv(path)
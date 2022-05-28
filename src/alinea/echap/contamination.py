# rapporter de septo3D dispersion.alep.interface

import random
from functools import reduce

def as_deposits(g,DU_stock, ud):
    """ returns deposits property build from deposited DU stock and ud property returned by septo3D
    """
    total_ud = int(round(sum(ud.values())))
    deposited_DU = min(len(DU_stock), total_ud)
    #print 'emited = %d, deposits = %d'%(total_emited, canopy_deposits)
 
    if deposited_DU>0:
        # first repartition of DUS along int(du)
        deposits = {}
        for k,v in ud.items():
            deposits[k] = []
            for i in range(v):
                deposits[k].append(DU_stock.pop(0))
            deposited_DU -= int(v)

        # repartition of remaining DUs, priority given to elements with highest deposits
        sorted_id = sorted(ud, key=ud.get)#last id is the bigest, that will be poped first
        while deposited_DU > 0:
            id = sorted_id.pop()
            deposits[id].append(DU_stock.pop(0))
            deposited_DU -= 1
    else:
        deposits = {}
    return deposits

class SimpleSoilInoculum:
    """ class useful for debug """
    def __init__(self, fsol = 0.5, DU_generator = None, sporulating_fraction = 0.01, pFA = 6.19e7, pNoSpo = 0.36, dh = 0.1, nAng=10, domain_area = 1, convUnit=0.01):
        """
        DU_generator : a class with 'create_stock(N)' method that can produce a list of DU of a given contaminant
        sporulating_fraction : fraction of soil that is sporulating 
        """
        self.opt = {'reference_surface' : domain_area, 'dh': dh, 'nAng': nAng, 'convUnit':convUnit}
        self.areaspo = sporulating_fraction
        self.pFA = pFA
        self.pNoSpo = pNoSpo
        self.generator = DU_generator
        self.fsol = 0.5
        
    def emission(self, g, weather_data):
        stock = []
        rain = weather_data[['rain']].sum().item()
        if rain > 0 :
            Ip = rain / len(weather_data)
            Isol = self.fsol * Ip
            eclins = self.pFA * self.pNoSpo * self.areaspo * Isol * self.opt['reference_surface']
            stock = self.generator.create_stock(eclins)
            print('eclins sol: %d'%(len(stock)))
        return stock

class SimpleContamination:
    """
    Class used for debug 
    """
    
    def __init__(self, fracud = 0.1, dh = 0.1, nAng=10, 
                 domain=None, domain_area=1, convUnit=0.01):
        """
        Setup model
        - reference surface is the soil surface (m2) occupied by plants (inverse of density)
        - dh is the heigh of layer
        - convUnit is the proportionality between scene units and meter
        """
        self.opt = {'reference_surface' : domain_area, 'dh': dh, 'nAng': nAng, 'convUnit':convUnit}
        self.domain = domain
        self.fracud = fracud
        
    def contaminate(self, g, DU, weather_data,label='LeafElement'):
        vids = [n for n in g if g.label(n).startswith(label)]
        areas = g.property('area')
        vids = [vid for vid in vids if vid in g.property('geometry') if areas[vid]>0.]

        n = len(vids)
        d = int(len(DU) / n)
        if d < 1:
            ud = {k:1 for k in random.sample(vids, len(DU))}
        else:
            ud = {k:d for k in vids}
        deposits = as_deposits(g,DU, ud)
        return deposits

class SimpleTransport:
    """
    Class wraping splash model of septo3D for alep generic model
    """
    
    def __init__(self, dh = 0.1, nAng=10, domain = None,
                 domain_area=1, convUnit=0.01, wash = True):
        """
        Setup splash model
        - reference surface is the soil surface (m2) occupied by plants (inverse of density)
        - dh is the heigh of layer
        - convUnit is the proportionality between scene units and meter
        """
        self.opt = {'reference_surface' : domain_area, 'dh': dh, 'nAng': nAng, 'convUnit':convUnit}
        self.wash = wash
        self.domain = domain
        
    def disperse(self, g, DU, weather_data):
        vids = [n for n in g if g.label(n).startswith('LeafElement')]
        areas = g.property('area')
        vids = [vid for vid in vids if vid in g.property('geometry') if areas[vid]>0.]

        emited_DU = reduce(lambda x,y: x + y, list(DU.values()))
        n = len(vids)
        d = int(len(emited_DU) / n)
        if d < 1:
            ud = {k:1 for k in random.sample(vids, len(emited_DU) )}
        else:
            ud = {k:d for k in vids}
        
        deposits = as_deposits(g,emited_DU, ud)
                    
        return deposits

""" Reconstruction of Wheat for Boigneville data
"""

from alinea.adel.AdelR import devCsv
from openalea.deploy.shared_data import shared_data
import alinea.echap
from alinea.adel.stand.stand import agronomicplot
from alinea.adel.astk_interface import AdelWheat


def devT_Mercia():
    axeT = shared_data(alinea.echap, 'Mercia_axeT.csv')
    dimT = shared_data(alinea.echap, 'Mercia_dimT.csv')
    earT = shared_data(alinea.echap, 'Mercia_earT.csv')
    phenT = shared_data(alinea.echap, 'Mercia_phenT.csv')
    ssisenT = shared_data(alinea.echap, 'ssi2sen.csv')
    
    return devCsv(axeT,dimT,phenT,earT,ssisenT)
    

def Mercia_2010(age = 300, nsect=3, length = 0.2, width=0.15, sowing_density=150, plant_density=150, inter_row=0.125, seed=1, sample='random'):
    devT = devT_Mercia()
    nplants, positions, domain, domain_area, convUnit = agronomicplot(length=length, 
                                                            width=width, 
                                                            sowing_density=sowing_density, 
                                                            plant_density=plant_density,
                                                            inter_row=inter_row)
    adel = AdelWheat(nplants=nplants, positions = positions, nsect=nsect, devT=devT, seed= seed, sample=sample)
    g = adel.setup_canopy(age)
    return g, adel, domain, domain_area, convUnit, nplants
    
    
age_at_application = {'Mercia_2010': {'2011-04-19': 1166,
                                      '2011-05-11' : 1500}}
                                
                                   
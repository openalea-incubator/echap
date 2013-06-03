from alinea.echap.wheat_mtg import *

from alinea.echap.interception_leaf import *
from alinea.echap.interfaces import pesticide_interception

    
##################################### test interception    
    
def test_intercept_elevation():
    g90 = adel_mtg()
    g45 = adel_mtg()   
    g1 = adel_mtg()   

    interception_model_90 = CaribuInterceptModel(elevation=90.)
    g90 = pesticide_interception(g90, interception_model_90, product_name='Opus', dose=1.5)
    print g90.property('surfacic_doses')

    interception_model_45 = CaribuInterceptModel(elevation=45.)
    g45 = pesticide_interception(g45, interception_model_45, product_name='Opus', dose=1.5)
    print g45.property('surfacic_doses')

    interception_model_1 = CaribuInterceptModel(elevation=1.)
    g1 = pesticide_interception(g1, interception_model_1, product_name='Opus', dose=1.5)
    print g1.property('surfacic_doses')



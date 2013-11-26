from macros_annual_loop import setup_canopy
from alinea.echap.recorder import LeafElementRecorder

g, adel, domain, domain_area, convUnit, nplants = setup_canopy(length=0.2, nsect=3)


def test_setup():
    ages = range(100,2500, 100)
    recorder = LeafElementRecorder()
    for a in ages:
        g=adel.setup_canopy(a)
        header={'dd':a}
        for vid in g:
            if g.label(vid).startswith('LeafElement'):
                n = g.node(vid)
                header.update({'plant' : n.complex().complex().complex().complex().label,
                                'axe' : n.complex().complex().complex().label,
                                'metamer' : int(''.join(list(n.complex().complex().label)[7:])),
                                'organ' : n.complex().label,
                                'ntop' : n.complex().ntop,
                                'id' : n._vid
                         })
 
                recorder.record(n,header)
        adel.plot(g)
    return g, recorder

def test_grow():
    g=adel.setup_canopy(100)
    recorder = LeafElementRecorder()
    for i in range(24):
        g=adel.grow_dd(g,100)
        header={'dd':adel.canopy_age}
        for vid in g:
            if g.label(vid).startswith('LeafElement'):
                n = g.node(vid)
                header.update({'plant' : n.complex().complex().complex().complex().label,
                                'axe' : n.complex().complex().complex().label,
                                'metamer' : int(''.join(list(n.complex().complex().label)[7:])),
                                'organ' : n.complex().label,
                                'ntop' : n.complex().ntop,
                                'id' : n._vid
                         })
 
                recorder.record(n,header)

        adel.plot(g)
    return g, recorder
        

    
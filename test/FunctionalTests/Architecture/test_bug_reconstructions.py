from alinea.echap.architectural_reconstructions import EchapReconstructions
variety = 'Tremie12'
nplants = 10
nsect = 7
disc_level = 5
reconst = EchapReconstructions()
adel = None
tries = 0.
while adel is None:
    tries += 1
    print('='*80)
    print(tries)
    print('='*80)
    try:
        adel = reconst.get_reconstruction(name=variety, nplants = nplants, nsect = nsect, 
                                            disc_level = disc_level, aspect = 'line')
        print('='*80)
        print("C EST BON ")
        print('='*80)
    except:
        pass
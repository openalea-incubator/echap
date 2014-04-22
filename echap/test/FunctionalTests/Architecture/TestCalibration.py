import pandas
import numpy
import os
import matplotlib.pyplot as plt

from alinea.adel.WheatTillering import WheatTillering
import alinea.echap.architectural_data as archidb

from multiprocessing import Pool,freeze_support
import itertools

plt.ion()

def test_axis_dynamics():
    pgen = archidb.Mercia_2010_plantgen()
    tdata = archidb.Tillering_data_2010()
    
    # fit by Mariem
    primary_proba = pgen['decide_child_axis_probabilities']
    plant_density = pgen['plants_number']
    ears_density = pgen['ears_density']
    ears_per_plant = float(ears_density) / plant_density
    nff = 12
    m = WheatTillering(primary_tiller_probabilities=primary_proba, ears_per_plant = ears_per_plant, nff=nff)
    
    #newfit
    obs = tdata['tillering']
    obs = obs[obs.Var=='Mercia']
    primary_proba={'T2': obs['MB'].values + obs['TC'].values}#handle damages to MB and simplify TC
    for w in ('T1','T3','T4','T5'):
        primary_proba[w] = obs[w].values
    ears_per_plant = obs['MB'].values + obs['TT'].values
    nff = float(obs['Nff'].values)
    new_m = WheatTillering(primary_tiller_probabilities=primary_proba, ears_per_plant = ears_per_plant, nff=nff)
    
    # todo : add comparison to primary emission
    
    #test plantgen fit&
    
    fit = m.axis_dynamics(plant_density = plant_density)
    
    new_fit = new_m.axis_dynamics(plant_density = float(tdata['ear_density_at_harvest']['Mercia']) / new_m.ears_per_plant)
    fit.plot('HS', ['total','primary','others'], style=['--r','--g','--b'])
    new_fit.plot('HS', ['total','primary','others'],style=['-r','-g','-b'])
    plt.plot([1,13],[tdata['plant_density_at_emergence']['Mercia'],tdata['ear_density_at_harvest']['Mercia']] ,'ob')

def setAdel(TT_stop, var):
    from alinea.adel.plantgen.plantgen_interface import gen_adel_input_data,plantgen2adel
    from alinea.adel.stand.stand import agronomicplot
    from alinea.adel.astk_interface import AdelWheat
    from alinea.adel.AdelR import devCsv
    pgen = archidb.Mercia_2010_plantgen()
    #pgen = archidb.Tremie_2011_plantgen()
    tdata = archidb.Tillering_data_Mercia_Rht3_2010_2011()
    #tdata = archidb.Tillering_data_Tremie1_2011_2012()
    obs = tdata['tillering']
    obs = obs[obs.Var==var]
    #primary_proba={'T2': obs['MB'].values + obs['TC'].values}#handle damages to MB and simplify TC
    primary_proba={}
    for w in ('T1','T2','T3','T4'):
        primary_proba[w] = obs[w].values
    ears_per_plant = obs['MB'].values + obs['TT'].values
    nff = float(obs['Nff'].values)
    new_m = WheatTillering(primary_tiller_probabilities=primary_proba, ears_per_plant = ears_per_plant, nff=nff)
    #update avec les probas primares de Tillering data
    pgen['decide_child_axis_probabilities'] = primary_proba
    # replace plant_number par eardensity / ear_per_plant, avec earsperplant =tillering data['MB'].values + Tillering data['TT'].values
    pgen['plants_number'] = float(tdata['ear_density_at_harvest'][var]) / new_m.ears_per_plant
    #
    pgen['delais_TT_stop_del_axis'] -= TT_stop
    tx = pgen['dynT_user'].a_cohort[0]
    axeT_, dimT_, phenT_, _, _, _, _, _, _, _, _ = gen_adel_input_data(**pgen)
    axeT, dimT, phenT = plantgen2adel(axeT_, dimT_, phenT_)
    devT = devCsv(axeT, dimT, phenT)

    nplants, positions, domain, domain_area, convUnit = agronomicplot(length=0.5, 
                                                            width=0.5, 
                                                            sowing_density=pgen['plants_density'], 
                                                            plant_density=pgen['plants_number'],
                                                            inter_row=0.125)
    adel = AdelWheat(nplants=nplants, positions = positions, nsect=10, devT=devT, seed= 1, sample='random')
    return adel, nplants, positions, domain, domain_area, convUnit
    

    
def sample_adel():
    adel, nplants, positions, domain, domain_area, convUnit = setAdel(0, 'Mercia')
    return adel, domain

def simLAI(TT_stop, var):
    
    from alinea.adel.postprocessing import axis_statistics, plot_statistics 
    
    adel, nplants, positions, domain, domain_area, convUnit = setAdel(TT_stop,var)
    dd = range(0,2000,100)
    outs = [adel.get_exposed_areas(g, convert=True) for g in (adel.setup_canopy(age) for age in dd)]
    new_outs = [df for df in outs if not df.empty]
    out = reduce(lambda x,y: pandas.concat([x,y],ignore_index=True), new_outs)
    axstat = axis_statistics(out, domain_area, convUnit)    
    res =  plot_statistics(axstat, nplants, domain_area)
    return res, tx
    #res['HS'] = res.ThermalTime*tx
    #print res.plot('HS','PAI_vert',color='b')

def compare_LAI():
    sim = test_adel()
    #Mercia
    #obs = pandas.DataFrame({'TT':[1137,1281,1796],'PAI':[3.27,4.03,4.05]})
    #Rht3
    obs = pandas.DataFrame({'TT':[1137,1281,1796],'PAI':[3.40,3.96,3.45]})
    #Tremie2012
    #obs = pandas.DataFrame({'TT':[923,1260,1550],'PAI':[1.1,4.31,4.43]})
    #Tremie2013
    #obs = pandas.DataFrame({'TT':[917.5,1018.3,1377,1612.25],'PAI':[1.73,2.7,2.78,1.95]})
    #sim.plot('ThermalTime','PAI_vert')
    #obs.plot('TT','PAI',style='o')
    sim['HS'] = sim.ThermalTime*tx
    sim.plot('HS','PAI_vert',color='g')
    obs['HS'] = obs.TT*tx
    obs.plot('HS','PAI',style='o',color='r')

def draft_TC(adel, domain,thermal_time, pov_dirpath='.'):
    from alinea.adel.povray import povray
    from alinea.adel.mtg_interpreter import plot3d   
    
    g =adel.setup_canopy(thermal_time)
    # colors
    colors_def = {'green':[0, 255, 0], 'red' : [255, 0, 0]}
    greeness = g.property('is_green')
    colors = {k:colors_def['green'] if greeness[k] else colors_def['red'] for k in greeness}
    # scene
    scene = plot3d(g, colors=colors)
    # plot3D
    # povray
    pov_filename = os.path.join(pov_dirpath, '%s%s%s' % ('scene_', thermal_time, '.pov'))
    povray_image_filepath, stand_box_image_filepath = povray.povray(scene, pov_filename, 200.0, 50, 4288, 2848, domain, 0, 0, 'perspective', False, 'povray')
    # genFilepath
    from alinea.adel.povray import post_processing
    outpath = './counts.txt'
    post_processing.count_pixels(povray_image_filepath, stand_box_image_filepath, colors_def.values(), outpath)
    return povray_image_filepath, stand_box_image_filepath, colors_def
    
def color_count(image,RGBcolor = (0,255,0)):
    import cv2
    BGRcolor = numpy.array(RGBcolor[::-1])
    res = cv2.inRange(image, BGRcolor, BGRcolor)
    return res.sum() / 255.
    
def count_pixels(povray_image_filepath, stand_box_image_filepath, colors_def):
    import cv2
    im = cv2.imread(povray_image_filepath)
    box = cv2.imread(stand_box_image_filepath)
    mask = numpy.uint8(box[:,:,0] / box.max() * 255)
    total = mask.sum() / 255.
    masked=cv2.bitwise_and(im,im,mask=mask)
    return {k:color_count(masked,v) / total for k,v in colors_def.iteritems()}
    
    
    

    
    
def func_star(a_b):
    return test_adel(*a_b)
def main():
    p = Pool()
    TT_stop=[0,600,400,200]
    #var = "Mercia"
    var = "Rht3"
    #var = "Tremie"
    p.map(func_star, itertools.izip(TT_stop, itertools.repeat(var)))
# if __name__=="__main__":
    # freeze_support()
    # main()


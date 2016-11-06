import json
import os
import glob
import numpy
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.light import vecteur_direction, diffuse_source
from alinea.echap.cache_simulation import (cache_simulation_path, get_canopy,
                                           build_canopies)


def light_simulation_parameters(tag='reference', reset=False):
    path = cache_simulation_path(tag) / 'light_simulation_parameters.json'
    if not reset:
        try:
            with open(path, 'r') as input_file:
                cached = json.load(input_file)
            return cached
        except IOError:
            pass
    pars={}
    genos = ('Mercia', 'Rht3', 'Tremie12', 'Tremie13')
    pars['spraying_angle'] = {k: 17 for k in genos}
    pars['lambda0'] = {k: 0.6 for k in genos}
    with open(path, 'w') as output_file:
        json.dump(pars, output_file, sort_keys=True, indent=4,
                  separators=(',', ': '))
    #
    return pars


def lambda0(tag='reference', variety='Tremie12'):
    pars = light_simulation_parameters(tag)
    return pars['lambda0'][variety]


def tag_to_light(tag='zenith', sim_tag='reference', variety='Tremie12'):
    pars = light_simulation_parameters(sim_tag)
    spraying_angle = pars['spraying_angle'][variety]
    if tag == 'zenith':
        return [(1, (0, 0, -1))]
    elif tag == '57.5':
        return [(1. / 24, vecteur_direction(90 - 57.5, az)) for az in
                range(0, 360, 15)]
    elif tag == 'spray':
        return [(1. / 24, vecteur_direction(90 - spraying_angle, az)) for az in
                range(0, 360, 15)]
    elif tag == 'soc':
        return diffuse_source(46)
    elif tag.startswith('lai2000r'):
        zenith = 7, 23, 38, 53, 68
        r = int(tag.split('r')[1]) - 1
        return [(1. / 24, vecteur_direction(90 - zenith[r], az)) for az in
                range(0, 360, 15)]
    else:
        raise ValueError('Unknown light tag: ' + tag)


def tag_to_zenith(tag='zenith', sim_tag='reference', variety='Tremie12'):
    pars = light_simulation_parameters(sim_tag)
    spraying_angle = pars['spraying_angle'][variety]
    if tag == 'zenith':
        return 0
    elif tag == '57.5':
        return 57.5
    elif tag == 'spray':
        return spraying_angle
    elif tag == 'soc':
        energy, directions = zip(*diffuse_source(46))
        z = numpy.array(zip(*directions)[2])
        return numpy.mean(numpy.degrees(numpy.arccos(-z)))
    elif tag.startswith('lai2000r'):
        zenith = 7, 23, 38, 53, 68
        r = int(tag.split('r')[1]) - 1
        return zenith[r]
    else:
        raise ValueError('Unknown light tag: ' + tag)


def cache_light_path(tag, rep):
    sim_path = cache_simulation_path(tag, rep)
    path = sim_path / 'light'
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def illuminate_canopies(light_tag='zenith', z_soil=0, variety='Tremie12',
                        nplants=30, tag='reference', rep=1, start=None,
                        stop=None, by=None, at=None, reset=False,
                        reset_build=False, reset_reconstruction=False):
    print('Check/build canopies..')
    dd_range = build_canopies(variety=variety, nplants=nplants, tag=tag,
                              rep=rep, start=start, stop=stop, by=by, at=at,
                              reset=reset_build,
                              reset_reconstruction=reset_reconstruction)
    print('check/compute canopy illumination')
    head_path = cache_light_path(tag, rep) / light_tag + '_z' + str(
        z_soil) + '_' + variety.lower() + '_' + str(nplants) + 'pl_'

    pattern = head_path + '*.json'
    done = glob.glob(pattern)
    if len(done) > 0:
        done = map(lambda x: x.split('pl_')[1].split('.')[0], done)
        done = map(lambda x: '-'.join(x.split('_')), done)

    missing = dd_range
    if not reset:
        missing = [d for d in dd_range if d not in done]

    if len(missing) > 0:
        for daydate in missing:
            filename = head_path + '_'.join(daydate.split('-')) + '.json'
            g = get_canopy(variety=variety, nplants=nplants, daydate=daydate,
                           tag=tag, rep=rep)
            meta = g.property('meta').values()[0]
            sources = tag_to_light(light_tag, sim_tag=tag, variety=variety)
            c2u = {v: k for k, v in CaribuScene.units.iteritems()}
            units = c2u.get(meta['convUnit'])
            cscene = CaribuScene(g, sources, pattern=meta['domain'],
                                 soil_mesh=1, z_soil=z_soil, scene_unit=units)
            raw, aggregated = cscene.run(direct=True, infinite=True,
                                         simplify=True)
            soil, _ = cscene.getSoilEnergy()
            with open(filename, 'w') as output_file:
                saved = {'raw': raw, 'aggregated': aggregated, 'soil': soil}
                json.dump(saved, output_file, sort_keys=True, indent=4,
                          separators=(',', ': '))

    return dd_range


def get_light(light_tag='zenith', z_soil=0, variety='Tremie12', nplants=30,
              daydate='T1', tag='reference', rep=1, reset=False, reset_build=False,
              reset_reconstruction=False):
    dd_range = illuminate_canopies(light_tag=light_tag, z_soil=z_soil,
                                   variety=variety, nplants=nplants, tag=tag,
                                   rep=rep, at=[daydate], reset=reset, reset_build=reset_build,
                                   reset_reconstruction=reset_reconstruction)
    daydate = dd_range[0]
    filename = cache_light_path(tag, rep) / light_tag + '_z' + str(
        z_soil) + '_' + variety.lower() + '_' + str(nplants) + 'pl_' + '_'.join(
        daydate.split('-')) + '.json'
    with open(filename, 'r') as input_file:
        cached = json.load(input_file)
        raw = {w: {int(k): v for k, v in cached['raw'][w].iteritems()} for w in
               cached['raw']}
        aggregated = {
            w: {int(k): v for k, v in cached['aggregated'][w].iteritems()} for w
            in cached['aggregated']}
    return raw, aggregated, cached['soil']

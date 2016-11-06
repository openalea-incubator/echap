import os
import glob
import pandas
# try:
#     import cPickle as pickle
# except ImportError:
#     import pickle

from alinea.adel.astk_interface import AdelWheat

import alinea.echap
from openalea.deploy.shared_data import shared_data
from alinea.echap.hs_tt import tt_hs_tag, daydate_range
from alinea.echap.architectural_reconstructions import echap_reconstructions


def cache_reconstruction_path(tag):
    path = shared_data(alinea.echap) / 'cache' / 'reconstructions' / tag
    if not os.path.exists(str(path)):
        os.makedirs(str(path))
    return path


def cache_simulation_path(tag, rep=None):
    path = shared_data(alinea.echap) / 'cache' / 'simulations' / tag
    if not os.path.exists(str(path)):
        os.makedirs(str(path))
    if rep is not None:
        path = shared_data(alinea.echap) / 'cache' / 'simulations' / tag / 'rep_' + str(rep)
        if not os.path.exists(str(path)):
            os.makedirs(str(path))
    return path


def cache_analysis_path(tag):
    path = shared_data(alinea.echap) / 'cache' / 'analysis' / tag
    if not os.path.exists(str(path)):
        os.makedirs(str(path))
    return path


def cache_canopy_path(tag, rep):
    sim_path = cache_simulation_path(tag, rep)
    path = sim_path / 'canopy'
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def get_reconstruction(variety='Tremie12', nplants=30, tag='reference', rep=1,
                       reset=False, reset_reconstruction=False):
    """ Return a new adel wheat instance"""
    # sim_path = cache_simulation_path(tag, rep)
    #
    # filename = sim_path / 'adel_' + variety.lower() + '_' + str(nplants) + 'pl.pckl'
    # if not reset:
    #     try:
    #         with open(filename) as saved:
    #             return pickle.load(saved)
    #     except IOError:
    #         pass
    echap = echap_reconstructions(tag, reset = reset_reconstruction)
    adel = echap.get_reconstruction(name=variety, nplants=nplants)
    # TO DO : make (part of ?) adel picklble !
    # with open(str(filename), 'w') as output:
    #     pickle.dump(adel, output)
    return adel


def build_canopies(variety='Tremie12', nplants=30, tag='reference', rep=1,
                   start=None, stop=None, by=None, at=None, reset=False,
                   reset_reconstruction=False):

    head_path = cache_canopy_path(tag, rep) / variety.lower() + '_' + str(
            nplants) + 'pl_'

    pattern = head_path + '*.pckl'
    done = glob.glob(pattern)
    if len(done) > 0:
        done = map(lambda x: x.split('pl_')[1].split('.')[0], done)
        done = map(lambda x: '-'.join(x.split('_')), done)

    if all(map(lambda x: x is None, [start, stop, by, at])):
        if len(done) > 0:
            dd_range = done
        else:
            raise ValueError(
                'No simulation found for ' + variety + ' please specify start, stop, (by) or at')
    else:
        dd_range = daydate_range(variety, tag, start, stop, by, at)

    missing = dd_range
    if not reset:
        missing = [d for d in dd_range if d not in done]

    if len(missing) > 0:
        adel = get_reconstruction(variety=variety, nplants=nplants, tag=tag,
                                  rep=rep,
                                  reset_reconstruction=reset_reconstruction)
        tths = tt_hs_tag(variety, tag)
        for d in missing:
            print d
            basename = head_path + '_'.join(d.split('-'))
            age = tths.set_index('daydate')['TT'][d]
            g = adel.setup_canopy(age=age)
            adel.save(g, basename=str(basename))
            midribs = adel.get_midribs(g)
            midribs.to_csv(basename + '_midribs.csv', index=False)

    return dd_range


def get_canopy(variety='Tremie12', nplants=30, daydate='T1', tag='reference',
               rep=1, load_geom=True, reset=False, reset_reconstruction=False):
    # check/build
    dd_range = build_canopies(variety=variety, nplants=nplants, tag=tag, rep=rep,
                   at=[daydate], reset=reset,
                   reset_reconstruction=reset_reconstruction)
    # load cached canopy
    daydate = dd_range[0]
    basename = cache_canopy_path(tag, rep) / variety.lower() + '_' + str(
        nplants) + 'pl_' + '_'.join(daydate.split('-'))

    return AdelWheat.load(basename=str(basename), load_geom=load_geom)


def get_midribs(variety='Tremie12', nplants=30, daydate='T1', tag='reference',
               rep=1, reset=False, reset_reconstruction=False):
    # check/build
    dd_range = build_canopies(variety=variety, nplants=nplants, tag=tag, rep=rep,
                   at=[daydate], reset=reset,
                   reset_reconstruction=reset_reconstruction)
    # load cached canopy
    daydate = dd_range[0]
    basename = cache_canopy_path(tag, rep) / variety.lower() + '_' + str(
        nplants) + 'pl_' + '_'.join(daydate.split('-'))
    return pandas.read_csv(basename + '_midribs.csv')

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 13:58:41 2013

@author: lepse
"""
from operator import itemgetter
from itertools import groupby
from numpy import exp,log

from alinea.echap.pesticide_data import milne_parameters
from functools import reduce

class SimcyclePesticide(Exception): pass       


def _dose_norm(dose, dose_max_ha):
    """ normalise doses(g.m-2) """
    dn = 0
    if dose_max_ha > 0 :
        dn = float(dose * 1e4) / dose_max_ha
    return dn


def _dose_response(dose, dose_max_ha, A, K): 
    return A * (1 - exp(-K * _dose_norm(dose,dose_max_ha)))


def _dose_equivalent(dose_other, dmax_ref, dmax_other, A_ref, A_other, K_ref, K_other):
    dneq = 0 
    dn_other = _dose_norm(dose_other,dmax_other)
    if A_ref > 0 and A_other > 0 :
        #print((float(A_other)/float(A_ref)) * (1 - exp(- K_other * dn_other)))
        dneq = - (1./K_ref) * log(1 - (float(A_other)/A_ref) * (1 - exp(- K_other * dn_other)))  
    return dneq * dmax_ref / 1e4


def add_doses(pref,pother, effect='protectant'):
    """ add equivalent dose of other in pref."""
    deq = 0
    if effect == 'protectant':
        deq = _dose_equivalent(pother['dose'], pref['dose_max_ha'], pother['dose_max_ha'], pref['Ap'],pother['Ap'], pref['Kp'],pother['Kp'])
    elif effect == 'eradicant':
        deq = _dose_equivalent(pother['dose'], pref['dose_max_ha'], pother['dose_max_ha'], pref['Ae'],pother['Ae'], pref['Ke'],pother['Ke'])
    else :
        raise SimcyclePesticide('Unknown effect: %s'%effect)
    pref['dose'] += deq
    return pref


def global_efficacy(doses, compound_parameters=[]):
    """
    compute global efficacy of a mixture of coumpounds
    
    :Parameters:
      - `doses` : A dict of ('compound' : doses) items giving the amount of product 'compound' present on the surface of the leaves (g.m-2)
      - `coumpound_parameters` : A dict of parameters dict. Parameters are :
          - `dose_max_ha`(float) - Maximum recommended dose of coumpound for a single dose application (g.ha-1)
          - `type_code` (int) - Code for the mode of action of coumpound
          - `Ap` and `Kp` (floats) - Parameters for dose/response curve for protectant effect of compound
          - `Ae` and `Ke` (floats) - Parameters for dose/response curve for eradicant effect of compound
          - `decay_rate` (float, [0,1]) - Decay rate of the active substance over time

  :Returns:
      - `efficacy` : A dict with the following items :
          - `protectant` (float, [0,1]) - Protectant efficacy of the active subsance. Affects the number of successful infections. No infection will occur when `protectant` = 1. 
          - `eradicant` (float, [0,1]) - Eradicant efficacy of the active subsance. Affects the fungal development. Fungal development will be stopped when `eradicant` = 1.
          
    """
    efficacy = {'protectant': None, 'eradicant': None}
    ptable = compound_parameters
    d = {}
    for k in doses:
        try:
            d[k]=ptable[k]
        except KeyError:
            raise SimcyclePesticide('product %s not found in parameter dict'%k)
        d[k].update(dose=doses[k])
    # Sort products in order to group them by mode of action with 'groupby'
    d = sorted(list(d.values()),key= itemgetter('type_code')) 
    d = [list(g) for k,g in groupby(d,key= itemgetter('type_code'))]
    
    A = {'eradicant':'Ae','protectant':'Ap'}
    K = {'eradicant':'Ke','protectant':'Kp'}
    for mode in ('eradicant','protectant'):
        # Sort products inside the lists to get reference product (highest A) in first position
        ds = [sorted(g, key=itemgetter(A[mode]), reverse=True) for g in d]
        # Reduce lists with equivalent doses computation (additive effect)
        dc = [reduce(lambda x,y: add_doses(x,y,effect=mode),g) for g in ds]
        # 1 - efficacy computation for reference products with updated doses
        effs = [(1-_dose_response(p['dose'], p['dose_max_ha'],p[A[mode]], p[K[mode]])) for p in dc]
        #print effs
        # Application of the multiplicative effect
        efficacy[mode] = 1 - reduce(lambda x,y: x * y, effs)
    return efficacy


def protectant_eradicant_efficacy(surfacic_doses, penetrated_doses, compound_parameters=[]):
    """ Implement pesticide efficacy_model using simcycle pesticide model """
    protectant = 0
    eradicant = 0
    if len(surfacic_doses) > 0:
        protectant = global_efficacy(surfacic_doses, compound_parameters)['protectant']
    if len(penetrated_doses) > 0:
        eradicant = global_efficacy(penetrated_doses, compound_parameters)['eradicant']
    efficacy = {"protectant":protectant, 'eradicant':eradicant}
    return efficacy


class PesticideEfficacyModel(object):
    """ Adaptor for PesticideEfficacy model compliying echap pesticide_efficacy protocol"""
    def __init__(self, compound_parameters = milne_parameters):
        self.compound_parameters = compound_parameters
    def efficacy(self, surfacic_doses, penetrated_doses):
        efficacy = protectant_eradicant_efficacy(surfacic_doses, penetrated_doses, self.compound_parameters)
        return efficacy

def as_doses(doses, name='Chlorothalonil'):
    return [{name:d} for d in doses]
    


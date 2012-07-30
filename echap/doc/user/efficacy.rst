
WP 2: Effect of pesticide on Infectious Cycle
#############################################


Description
===========

This module allows: 

- to compute globaly decay of dose on leaves as a function of time
- to compute the global efficacy on infection and lesion development rate  of a set of compound doses fungii
- to compute an erosion factor of compounds as a function of doses sprayed per year

It is based on the paper : A model of the effect of fungicides on disease-induced yield loss, for use in wheat disease management decision support system. Annals of Applied Biology. 151, p.113-125 by Milne et al., 2007

Interfaces
==========

The module comes with 3 functions

.. function:: dose_decay(decay_rate, initial_dose, time_step)

  :Parameters:
      - `decay_rate` (float, [0,1]) - Decay rate of the active substance over time
      - `initial_dose` (float) - Mass of product intercepted by the surface of the leaves (g.m-2)
      - `time_step` (float) - Number of days since the product application
      
  :Returns:
      - `current_dose` (float) - Updated mass of product on the surface of the leaves (g.m-2)
  
.. function:: global_efficacy(doses,coumpound_parameters)
  
  :Parameters:
      - `doses` : A dict of ('compound' : doses) items giving the amount of product 'compound' present on the surface of the leaves (g.m-2)
      - `coumpound_parameters` : A dict of ('compound' : (parameters)) items. Parameters are :
          - `dose_max`(float) - Maximum recommended dose of coumpounds for a single dose application (g.ha-1)
          - `action_mode_class` (int) - Code for the mode of action of coumpounds
          - `Ap` and `Kp` (floats) - Parameters for dose/protectant response curve of compounds
          - `Ae` and `Ke` (floats) - Parameters for dose/eradicant response curve of compounds

  :Returns:
      - `efficacy` : A dict with the following items :
          - `protectant` (float, [0,1]) - Protectant efficacity of the active subsance. Affects the number of successful infections. No infection will occur when `protectant` = 1. 
          - `eradicant` (float, [0,1]) - Eradicant efficacity of the active subsance. Affects the fungal development. Fungal development will be stopped when `eradicant` = 1.
          
.. TODO:: 
  - use Milne et al to implement additive and multiplicative composition of products for efficacy
  - get/parameterise parameters for dose response curves ofcompounds, using global simulation of efficacy
  - allows infectious cycle (WP1) to respond to efficacy
  - `erosion_factor` (list of float) - Scaling factor for efficacy of compounds linked to erosion of efficacy due to evolution of resistance of strains

.. function:: erode_products(erosion_factor, compound_parameters)

    :Parameters:
      - `erosion_factor` (float, [0,1]) - Scaling factor for efficacy of compounds linked to erosion of efficacy due to evolution of resistance of strains
      - `coumpound_parameters` : A dict of ('compound' : (parameters)) items. Parameters are :
          - `dose_max`(float) - Maximum recommended dose of coumpounds for a single dose application (g.ha-1)
          - `action_mode_class` (int) - Code for the mode of action of coumpounds
          - `Ap` and `Kp` (floats) - Parameters for dose/protectant response curve of compounds
          - `Ae` and `Ke` (floats) - Parameters for dose/eradicant response curve of compounds
          
    :Returns:
      - `new_pars` : Same dict as `coumpound_parameters` with updated `Ap` and `Ae` linked to erosion of efficacy due to evolution of resistance of strains
      
.. function:: erosion_factor(cumulative_doses,current_erosion_factor,erosion_type)
  
  :Parameters:
      - `cumulative doses` (float) - Masses of product that acted on fungi during the year(g.lesionsurface ???)
      - `current_erosion_factor` (float) - Scaling factor for efficacy of compounds linked to erosion of efficacy due to evolution of resistance of strains at the begiging of the year.
      - `erosion_type` (string) : type of erosion (slow, rapid...) or erosion=f(t) function

  :Returns:
      - `new_erosion_factor` (float, [0,1]) - Updated scaling factor for efficacy of compounds linked to erosion of efficacy due to evolution of resistance of strains
            


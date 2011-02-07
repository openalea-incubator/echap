
WP 2: Effect of pesticide on Infectious Cycle
#############################################


Description
===========

This module allows: 

- to compute globaly decay of dose on leaves as a function of time
- to compute the global efficacy on infection and lesion development rate  of a set of compound doses fungii
- to compute an erosion factor of compounds as a function of doses sprayed per year

It is based on Milne etal.

Interfaces
==========

The module comes with 3 functions

.. function:: global_efficacy(doses,compounds,erosion_factor,dose_response_curves)
  
  :Parameters:
      - `doses` (list of float) - Masses of product present at the leaf surface the leaf (g.m-2)
      - `coumpounds` (list of string) - Compound name or compound ntype of acttion
      - `erosion factor` (list of float) - Scaling factor for efficacy of compounds linked to erosion of efficacy dueto evolution of resistance of strains.
      - 'dose_response_curves' (dict) : parameters for dose/response curve of compounds

  :Returns:
      - 'eradicant_efficacy' (float, [0,1]) - efficacy of the set of products for reducing probability of infectionUpdated mass of product at the end of the day (g.m-2) 
      - 'curative_efficacy' (float [0,1]) - efficacy of the set of products for reducing the rate of fungal development
            

.. TODO:: 
  - use Milne et al to implement additive and multiplicative composition of products for efficacy
  - get/parameterise parameters for dose response curves ofcompounds, using global simulation of efficacy
  - allows infectious cycle (WP1) to respond to efficacy

.. function:: erosion_factor(cumulative_doses,current_erosion,erosion_type)
  
  :Parameters:
      - `cumulative doses` (float) - Masses of product taht acted on fungi during the year(g.lesionsurface ???)
      - `curent_erosion factor` (float) - Scaling factor for efficacy of compounds linked to erosion of efficacy dueto evolution of resistance of strains at the begiging of the year.
      - 'erosion_type' (string) : type of erosion (slow, rapid...) or erosion=f(t) function

  :Returns:
      - 'new_erosion' (float, [0,1]) - 
            


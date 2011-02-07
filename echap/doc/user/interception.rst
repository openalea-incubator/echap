.. _echap_interception:


WP 3: Interception
############################

:Version: |version|
:Release: |release|
:Date: |today|

.. .. seealso:: :ref:`echap_dispersion_reference`.


Description
=============

Intercept raindrops or droplets of pesticide with Monte Carlo integration.


- Inputs:

    - Rain R [mm/h] -> to build distribution of velocity and diameter of raindropd
    - Spray mixture : distributions of velocities and diameters for different nozzles + concentration of pesticide [dose l/ha] (+ coef of drift loss) :TODO 
    
    - 3D canopy structure ie list of primitives:
      * one cylinder per axis,
      * list of quad for leaves
      * quad of soil 
      * id primitives with plants and MTG  
  
- Ouput:
    - list of same order of primitive with eg quantities of liquid or physics parameters of impact  
    

Tutorials and Examples
=======================

:Todo: Tutorials


:Example: colorant de saucisses


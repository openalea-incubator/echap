.. _echap_interception:


WP 3: Interception
############################

:Version: |version|
:Release: |release|
:Date: |today|

.. .. seealso:: :ref:`echap_dispersion_reference`.


Description
=============

Intercept raindrops or droplets of pesticide.


We consider that spray is homogeneous within a ramp, ie that arrangement of nozzles is such that it ensure this homogeneity.

we consider that rain and pesticide are vertical.

This will combine the computation of the fraction intercepted using a ligth interception model (fractalysis or caribu)
and an empirical model insprired from Ulbricht et al. that gives probability distribution function (pdf) of droplets per unit time and unit surface(pdf, drops/m2/s).
and associated distribution of velocity.
As a result, the model will compute physical statistics (Weber, kinetic energy, Oh) and associated fraction of product deposited, lost by runoff and lost by splash.

For pesticide the model will probably be simplified to fraction deposited = 1

In this case pdf could be use to estimate the surface touched by the pesticide (that will increase concentration locally)

Two strategies are possible: 

- Either we consider that pdfs are the same throughout the canopy, ie that the period of time of the interception event is long enougth to ensure an unbiased pdf 
on all elementary surfaces, even those that receive a small amount of rain/product.

- Either we consider that pdf can be locally variable, due to sampling, and that this variability should be simulated explicitely. 
In that case we sould include the sampling in the model


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


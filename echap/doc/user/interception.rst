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
    - Spray mixture : distributions of velocities and diameters for different nozzles + concentration of pesticide [dose l/ha] (+ coef of drift loss)
    - 3D canopy structure
  
- Ouput:

    - spatial distribution of liquids (rather rain or spray) within the canopy ie, liquids on plant organs and soil    


Tutorials and Examples
=======================

:Todo: Tutorials


:Example: colorant de saucisses

Dataflow
==========


.. .. dataflow:: Alinea.Echap.Concept - Annual loop
..    :width: 50%

..	Conceptual dataflow simulating one year experiment.

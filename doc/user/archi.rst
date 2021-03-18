.. _echap_archi:


WP1: Architecture Simulation
############################

:Version: |version|
:Release: |release|
:Date: |today|
:Authors: Christian Fournier, Camille Chambon, Mariem Abichou, Bruno Andrieu, Corinne Robert
:Target: users, developers and administrators


Objects manipulated
===================

MTG representing the plant with dimensions, tissue state, parameters (rate senescence)


ADEL
====

Adel will provide: 

  * Parameterisation procedure that allows to generate ADEL parameter tables from field experiment databases
  * a module/function that is able to update canopy sate / 3D structure (the MTG) during a day, using the structutre at previous step, adel tables and temperature.
  * a visualisation routine to see plants in 3D
  * Adaptors to Pearl, to generate cropgrowth parameters from adel parameterisation tables




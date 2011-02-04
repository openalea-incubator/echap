.. _echap_archi:


WP1: Architecture Simulation
############################

:Version: |version|
:Release: |release|
:Date: |today|
:Authors: Christian Fournier, Camille Chambon, Mariem Abichou, Bruno Andrieu
:Target: users, developers and administrators


Objects manipulated
===================

MTG representing the plant with dimensions, tissue state, parameters (rate senescence…)


ADEL
====

Plant/canopy model initialisation as an MTG = f(architectural field data)
Plant/canopy growth (updated MTG) = f(current MTG,  Temperature, dt)
Plant/canopy senescence (updated MTG + current progress = f(current MTG, Temperature, dt)
3D canopy (PlantGL scene graph)= f(current MTG)


Tutorials and Examples
=======================
:TODO: 


Dataflow
==========
:TODO: 



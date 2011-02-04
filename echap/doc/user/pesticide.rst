.. _echap_pesticide:


WP 3: Pesticide decay
############################

:Version: |version|
:Release: |release|
:Date: |today|
:Authors: Pierre BENOIT, add Erik Van den BERG
:Target: users, developers and administratorss

.. .. seealso:: :ref:`echap_pesticide_reference`.


General introduction
====================

Pesticide decay is designed to simulate amounts of fungicide residues on the aerial parts (leaves) of wheat plants as a function of time. It describes elementary processes of phodegradtion, penetration, volatilization and wash off. Output is a concentration/amount as a function of time for each plant element.

The model is coupled : 
with a light model adapted to field conditions (Caribu model), in order to be able to calculate accurately the
light interception and the photodegradation on phytoelements during the crop cycle
with an interception module 
with a infection cycle module




Description of Pesticide decay's inputs
=======================================

Pesticide decay has several kinds of inputs available for the user:
 * inputs characterising the pesticide type
 * inputs characterising the initial pesticide concentration after application (interception)
 * inputs characterising the microclimate at the leaf scale (temperature, radiation)
 * inputs characterising the degree of infection of leaves
 

.. list-table::
    :widths: 10 50 50
    :header-rows: 1

    * - Column
      - Description
      - From
    * - **Pesticide** 
      - the name of the active ingredient
      - 
    * - **c0** 
      - the intial concentration after application
      - Interception
    * - **temp** 
      - the temperature at the leaf scale
      - Caribu  
    * - **uv** 
      - the UV radiations at the leaf scale
      - Caribu
    * - **is_infected** 
      - indicator for infected status
      - Septo3D



Description of Pesticide decay's outputs
========================================
Pesticide amount on MTG =f(time; T°,light, green/necrotic) from description of individual processes : photodegradation, penetration, volatilization, wash off

Pesticide has one output available for the user:

 * ouputs characterising the amounts of fungicide residues on each leaf sector as a function of time


Tutorials and Examples
=======================


Function
----------

Pesticide amount on MTG

.. math::
        	
    f(T^0, light, green/necrotic)

.. math::

    A(t) =A0-(A_PEN + A_DEG + A_VOL + A_WAS)

.. math::
    
    A0 initial amount on MTG


Pesticide penetration

.. math::

    A_PEN (t) =A0 exp (-k_PEN.t)

Pesticide photodegradation

.. math::

    A_DEG (t) =A0 exp (-k_DEG.t)

Pesticide volatilization

.. math::

    Volat Flux ref PEARL

Pesticide wash off

.. math::

    Wash off ref PEARL
    


Dataflow
========

.. dataflow:: Alinea.Echap pesticide_global
    :width: 100%

    The global dataflow associated with WP 3



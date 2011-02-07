.. _echap_pesticide:


WP 3: Pesticide leaf
############################

:Version: |version|
:Release: |release|
:Date: |today|
:Authors: Pierre BENOIT, add Alex Ter HALLE, Erik Van den BERG, Carole BEDOS
:Target: users, developers and administratorss

.. .. seealso:: :ref:`echap_pesticide_reference`.


General introduction
====================

Pesticide leaf is designed to simulate amounts of fungicide residues on the aerial parts (leaves) of wheat plants as a function of time. It describes elementary processes of phodegradtion, penetration, volatilization and wash off. Output is a concentration/amount as a function of time for each plant element.

The model is coupled : 
with a light model adapted to field conditions (Caribu model), in order to be able to calculate accurately the
light interception and the photodegradation on phytoelements during the crop cycle
with an interception module 
with a infection cycle module



Description of Pesticide leaf's inputs
=======================================

Pesticide decay has several kinds of inputs available for the user:
 * inputs characterising the pesticide type
 * inputs characterising the initial pesticide concentration after application (interception)
 * inputs characterising the microclimate at the leaf scale (temperature, radiation)
 * inputs characterising the degree of infection of leaves
 

.. list-table::
    :widths: 10 50 50
    :header-rows: 1

    * - Parameters
      - Description
      - From
    * - **Pesticide** 
      - the name of the active ingredient
      - 
    * - **c0** 
      - the intial concentration after application
      - Interception
    * - **leaf_temp** 
      - the temperature at the leaf scale
      - Caribu  
    * - **leaf_uv** 
      - the UV radiations at the leaf scale
      - Caribu
    * - **leaf_water** 
      - the UV radiations at the leaf scale
      - Caribu
    * - **is_infected** 
      - indicator for infected status
      - Septo3D

            

.. TODO:: 
    * Parametrisation of kinetic parameters according to local meteo 
    - DT50 photodegradation as a function of UV
        - from Alexandra experimental data for several fungicides
 
    - DT50 penetration as a function of temperature
        - from litterature data data see for instance Baur et Schönherr, 1995 Chemosphere, 30, 1331-1340

    * Check processes description in SURFATM Model

    * Check possible description of the distribution of concentrations (heterogenous vs homogeneous)

::

Function
----------

Pesticide amount on leaf sector

.. math::
        	
    f(T^0, light, green/necrotic)

.. math::

    A(t) =A_0-(A_{PEN} + A_{DEG} + A_{VOL} + A_{WAS})

.. .. math::
    
    :math:`A0` initial amount on MTG


Pesticide penetration

.. math::

    A_{PEN} (t) =A_0 exp (-k_{PEN} . t)

Pesticide photodegradation

.. math::

    A_{DEG} (t) =A_0 exp (-k_{DEG}.t)

Pesticide volatilization

.. .. math::

    Volat Flux ref PEARL

Pesticide wash off

.. .. math::

    Wash off ref PEARL


Description of Pesticide leaf 's outputs
========================================

Pesticide leaf has one output available for the user:

 * ouputs characterising the amounts of fungicide residues on each leaf sector as a function of time

.. function:: pest_leaf(pesticide, dose,leaf_temp,leaf_uv, leaf_water, is_infected)
  

  :Returns:
      - 'dose' (float) - Updated mass of product at the end of the day (g.m-2) 
      - 'volat' (float) - Mass of product volatilised (g)
      - 'washoff' (float) - Mass of product washed off (g)
      - 'penetrated' (float) - Mass of product penetrated (g)
      - 'photodegraded' (float) - Mass of product photodegraded (g)

.. TODO:: 
    * check possible link with WP2
	- integrate Surface and Penetrated Amount to link with dose-response curves cf Neil Paveley
 
    * link with WP3 Environmental Impact : check possible input of pesticide to the soil pesticide Module from Wash-off Amount - Cumulated amounts over a season ?

    * link with WP3 Environmental Impact :  check possible input of pesticide to the soil pesticide Module from Leaf Penetrated Amount - Cumulated amounts over a season ?

::


Links with other Modules
========================
The input of Pesticide leaf is to be checked to convert interception output (possibly volumes) as masses on leaf elements

The output of Pesticide leaf are to be checked for consistency  : 
    - with the input of environmental impact module

    - with the input of pesticide efficiency



RoadMap
=======

- Feb-June 2011 : complete information on 1st order equation parameters and dependecies from local meteo parameters and/or spatial distribution of intercepted product 
- Interact with the Pearl leaf inclusion see Pearl doc case 1
- Store data flow for parameter estimation either from experimental data or culculation
- Check consistency with WP2 decay function
- Plan calibration of decay function from data



Tutorials and Examples
=======================


Dataflow
========

.. dataflow:: Alinea.Echap pesticide_global
    :width: 100%

    The global dataflow associated with WP 3



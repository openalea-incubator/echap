
WP 3: Use of PRZM
############################

:Version: |version|
:Release: |release|
:Date: |today|


Description
=============

PRZM allows to compute the fate of pesticides on plant leaves and loss to the environment.

In the context of the project, it could be used in two cases : 

- case 1 : In this case the outputs of the annual loop are used to compute a DT50 of dissipation in the crop that will force the decay in the vegetation.

- case 2 (alternative to case 1) : to run PRZM independantly of other module, in one run, and get directly the fate of pesticide at the canopy scale level + soil.


Main input parameters
=====================

Whatever the case, przm needs to be parameterised with tables decribing : 

:Pesticides and metabolites:

 * Total number of pesticide applications occurring at different dates(1 to 50)
 * number of pesticide(s) in the simulation (1 to 3)
 * names of pesticide(s)
 * target application date (day, month, year)
 * chemical application method (soil applied, soil incorporation depth (cm), interception based on crop canopy)
 * target application rate of the pesticide(s) (kg ha-1)
 * application efficiency (fraction)
 * spray drift (fraction)
 * day, month when first half-life begins, number of days after begin of first half-life to swap from first to second bi-phase half-life.

:Soil:

 * Total depth of soil core (cm)
 * bulk density (kg L-1)
 * monthly values of soil surface albedo (12 values)
 * reflectivity of soil surface to longwave radiation (fraction)
 * height of wind speed measurement above the soil surface (m)
 * average monthly values of bottom boundary soil temperatures in degrees Celsius (12 values)
 * total number of horizons
 * thickness of the horizon (cm)
 * thickness of compartments in the horizon (cm)
 * organic carbon in the horizon (percent)
 * sand content in the horizon (%)
 * clay content in the horizon (%). 

:Climate:

 * Precipitation (cm day-1)
 * Pan evaporation data (cm day-1)
 * Average daily temperature (°C)
 * Wind speed (cm s-1)
 * Solar radiation (Langley)

:CropGrowth:

 * date of emergence, maturity, harvest
 * maximum rooting depth
 * maximum canopy height at maturation date
 * maximum interception storage of the crop
 * maximum aeral coverage of the canopy (see page 4-19 of Users manual)

Interface of use in case 1
==========================

Here PRZM is used in plant+soil mode, for all a season, with forced decay of products on plants.

One run per product.

.. function:: przm_forced_dissipation_on_crop(substance_parameters, spray_scenario, daily_meteo, soil parameters,lumped Coefficient,CropGrowth_parameters)

  
  :Parameters:
      - see above parameter tables
      - Lumped coefficient : to be fitted on output of daily loop
      - Substance parameters : see above
      - Spray scenario : date/dose of application, application efficiency, drift...
      - Meteo : see variable list above (daily time step)
      - Soil parameters (see PRZM users manual pages 4-24 to 4-29: capacity approach with water contents at field capacity and wilting point)
  
  :Returns:
      - mass balance for the soil (other terms already computed in case 1)
      - mass balance can be detailed at the level of every siingle soil layer or sumarised    

	
.. note::  in the case we only use the output for soil, a precice control of fraction of product with reduce exposure may be skipped, as we will not use the "plant" degradation/volatilisation outputs



Interface of use in case 3
==========================

Same as case 1, but with "Calculated" vegetation and including an estimation of the fraction of product in the reduce exposure class as input.

PRZM is used in plant+soil mode, for all a season.
One run per product.

RoadMap
=======

- Feb-June 2011 : complete inclusion of PRZM (call/return from python) / test of model at leaf scale. extract/depose input from/to mtg
- check consistency with WP2 decay function
- Plan calibration of decay function from data


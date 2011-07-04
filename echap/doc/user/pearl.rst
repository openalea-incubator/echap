
WP 3: Use of Pearl
############################

:Version: |version|
:Release: |release|
:Date: |today|


Description
=============

Pearl allows to compute the fate of pesticides on plant leaves and loss to the environment.

In the context of the project, it could be used in three cases : 

- case 1 : to compute localy (on small leaf parts) the fate of products at a daily (or hourly) time step. In this case the model is run with micrometeorological data computed at leaf level and becomes an alternative to other pesticide decay functions in the within day loop.


- case 2 : to compute the fate of pesticide in the whole system, taking into account the outputs of the detailed simulation to compute a DT50 of dissipation in the crop (option Lumped in pearl)

- case 3 : to run pearl independantly of other module, in one run, and get directly the fate of pesticide at the canopy scale level + soil. Vegetation is in mode 'Calculated'

Use in case 1
=============

The idea is to configure pearl inputs for a canopy consisting of a single leaf element. Pearl is also run without soil module. Model wil run 1 day

The model can be run for all leaf elements, or for different classes of leaf elements experiencing similar micloclimatic conditions.

One model run is needed per product.

.. note:: One important question will be to determine if we need the 'reduction factor' for exposure or if we are able to drive te reduction of exposure with microclimatic variables directly (light, rain interception, wind speed) .


The model can then be called with the following function: 


.. function:: pearl_leaf(dose,substance_parameters,local_meteo, TransportResistance)

  
  :Parameters:
      - `dose` (float) - Mass of product present at the leaf surface the leaf (g.m-2)
      - `coumpound_parameters` (dict) - list of parameters the physico-chemical properties gouverning decay equation (Phothodegradation, wash-off, see pearl doc)
      - `local_meteo` (dict) - list of hourly meteorological variables, either in the air if we use reduce exposure coefficient or value at the leaf scale.
      - transportResistance (dict): type of computation for resistance for boundary layer : either laminar, or ..see doc pearl, 2.2)
      

  :Returns:
      - 'dose' (float) - Updated mass of product at the end of the day (g.m-2) 
      - 'volat' (float) - Mass of product volatilised (g)
      - 'washoff' (float) - Mass of product washed off (g)
      - 'penetrated' (float) - Mass of product penetrated (g)
      - 'photodegraded' (float) - Mass of product photodegraded (g)
            

.. TODO:: 
    - PlantOnly.prl, configured 
	* to run one day  with date application == date debut
        * Calculated OptDspCrp_CompoundId
     -modify for every run PlantOnly.prl with : 
	- dose as:
	::
	    table Applications
	    xxxx AppCrpUsr dose_in_kg.ha-1 1.0 
	    end_table
	- compound parameters are : MolMass (g.mol-1), DT50PenCrp (d), DT50TraCrp (d), RadGloRef (W.m-2), FacWasCrp (m-1), TemRefDif (°C), CofDifAirRef (m2.d-1). May be also coefficient of reduction of parameters for less exposed leaves if needed (to decide).
        - option OptTraRes depending on bondary layer conditionn (laminar,...)	 

     - generate meteo file for one day, with local or global variables : radiation (kJ.m-2.h-1),Air temperature (°C), air Humidity (kPa), WindSpeed (m.s-1), Rainfall/Intercepted rain (mm)   

      - compute mass balance from pearl output (*.sum)
::



Interface of use in case 2
==========================

Here Pearl is used in plant+soil mode, for all a season, with forced decay of products on plants.

One run per product.

.. function:: pearl_forced_dissipation_on_crop(substance_parameters,spray_scenario,hourly_meteo, TransportResistance, soilparameters,lumpedCoefficient,CropGrowth_parameters,irrigation,interception_product)

  
  :Parameters:
      - CropGrowth (CrpPar,...): date of emergence, date of harvest, LAI(t), Rooting__depth(t), cultural coefficient for transpiration(t), hauteur(t), parameter for water extraction
      - irrigation ??
      - Lumped coefficient : to be fitted on output of daily loop at leaf scale
      - substance parameters : see above
      -spray scenario : date/dose of application
      - meteo : see variable list above (hourly time step)
      - Soil parameters (see pearl doc)
      - transportResistance (dict): type of computation for resistance for boundary layer : either laminar, or ..see doc pearl, 2.2)
      -  the fraction of product intercepted per spray event       

  :Returns:
      - mass balance for the soil (other terms already computed in case 1)
            

.. note::  in the case we only use the output for soil, a precice control of fraction of product with reduce exposure may be skipped, as we will not use the "plant" degradation/volatilisation outputs



Interface of use in case 3
==========================

Same as case 2, but with "Calculated" vegetation and including an estimation of the fraction of product in the reduce exposure class as input.

Here Pearl is used in plant+soil mode, for all a season.
One run per product.

The other echap module may be used to parameterise : 
-  the fraction of product intercepted per spray event 
- the fraction of product belonging to the reduce exposure class (kept constant during one run : to be kept as this ??)
- CropGrowth parameters...



RoadMap
=======

- Feb-June 2011 : complete inclusion of pearl (call/return from python) / test of model at leaf scale. extract/depose input from/to mtg
- check run time in a scenario with 5000 ? leaf element
- check consistency with WP2 decay function
- Plan calibration of decay fuinction from data



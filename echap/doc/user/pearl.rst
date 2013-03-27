
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



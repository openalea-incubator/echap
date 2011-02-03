.. _echap_scenarii:

.. ..test
WP6 : Scenarii
############################

:Version: |version|
:Release: |release|
:Date: |today|



Description
=============
What we do in this WP6 ?

Run *Big model* for validation against experimental data ? This is a question. 

Run *Big model* on a range of combination of genotypes and fungicide applications strategies, each time for a sequence  of climatic conditions representing 5 (?) years  to generate a database that will be used to analyse the environmental impact and productivity impacts of genotypes* fungicides stratrégies


Analyse the data base to 
-	identify general trends exploratory statistics such as ACP and regression analysis
-	identify appropriate genotypes * application strategies that lead to the best compromises between productivity and environmental impact

First proposal for the output table to be analysed
Variables describing the conditions
scenario definition : identifier + other variables (type of fungicide, total dose...)
genotype definition : identifier + other variables (height,stature,..)
(Choice of descriptive variables will impact the type of analysis
one year climate (variables describing the global climate - cumulated rain, mean temperature...)
Variables describing the results
Yield (or yield loss ?)
quantity of spores
quantity of fungicide emitted in the atmosphere
quantity of fungicide in the soil
		qt stored in the soil
		wash off 
		if secondary products : same (stored and wash off)
quantity of fungicide absorbed by the plant
quantity of fungicide photodegradated



:TODO: 
Define explicitely the list of input data that have to be given for running one year simulation
Are there some input depending on output of previous year ?

Define explicitely the list of output data that have to be stored as output of a one year simulation. 
 
The consortium must decide whether the whole model will be validated against experimental data (eg from Boigneville experiment) 

Define the time step  at which climatic data must be provided for “big model” simulations (daily ? hourly ?). Define the method for estimating data at this time step from the data existing in meteorological databases.

Define the meteorological database that will be used in the scenaris (what should be the range of conditions that are represented ? eg only grignon, multiple sites ?)
Get some examples of these meteorological database 

:MISCEALLENOUS:
Will it be necessary for the project to store daily output data ? 
 

Scenario
=======================

.. .. image:: annual_loop.png


.. .. dataflow:: Alinea.Echap.Concept - Annual loop
..    :width: 50%
..
..	Conceptual dataflow simulating one year experiment.

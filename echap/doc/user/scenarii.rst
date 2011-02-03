. _echap_scenarii:

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


Analyse the output data base to 
-	identify general trends by exploratory statistics such as ACP and regression analysis
-	identify appropriate genotypes * application strategies that lead to the best compromises between productivity and environmental impact. This will rise the difficulty of ranking according to several criteria. 

For this we may think to associate a positive or negative value to each variable (impact, yield, erosion of pesticide efficiency) and rank according to the cumulated value but we do not know if this is valuable.

An alternative way is to identify questions than can be answered unequivocally for instance “create a subset of the database according to one column, rank lines according to a column” an example is “amongst all scenariis with yield more than 90% of the maximum, which is the one with lower pesticide leach ?”
- 


First proposal for the output table to be analysed

One line per year of simulation. Each line consists in

Variables describing the conditions
scenario definition : identifier + variables describing the scenario (type of fungicide, total dose...)
genotype definition : identifier + variables describing the genotype (height, stature,..)
(Choice of descriptive variables will impact the type of analysis
one year climate (identifier + variables describing globally the climate - cumulated rain, mean temperature...)

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
pesticide efficiency erosion resulting from this year



first lines for data analysis

the table with yearly variables will be used for exploratory analysis and for building a table with variable integrated over 5? years that will serve for ranking genotytpe + fungicide application strategies




:TODO: 
Define explicitely the list of input data that have to be given for running one year simulation with Big Model
Are there some input depending on output of previous year ? ( at least : erosion of pesticide efficiency)

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


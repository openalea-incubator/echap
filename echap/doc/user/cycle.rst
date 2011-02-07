.. _echap_cycle:


WP1: Infectious cycle
#####################

:Version: |version|
:Release: |release|
:Date: |today|
:Authors: C Robert, C Fournier, C Pradal


This module will be based on Septo3D (Robert et al) model, adapted to respond to pesticide efficacy (see WP2)

Objects manipulated
===================

  * Lesions (on plant or soil) representing one cohorte (same day, same zone) with : Sector, effectif, cycle parameters (latency duration, growth rate,state in terms of tissue : chlorotic area, sporoluating area,)
  * dissemination units (UD) : droplets with spores that will produce one lesion


Interfaces
=========

Septo3D will provide: 

  * a soil inoculum model that allows to set-up the epidemics
  * germination model that predict lesions creation from UD, humidity at leaf surface, temperature, and global efficacy of a pesticide mixture (surfacic ?)
  * UD production model/lesions emptying that compute UD incorporated in droplets from rain drops impact,and lesion state
  * UD washing model, as a function of rain drops impacts
  * senescence response model of lesions
  * lesion development model as a function of temperature and global efficacy of (penetrated) pesticide mixture




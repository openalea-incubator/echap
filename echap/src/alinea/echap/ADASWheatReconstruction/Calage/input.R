#
#              Calage ADEL a partir de RDatabase
#              Fev 2011
#
#
#fonctions generiques importee
#
source("FunCalage.R")
source("FunConvMeteo.R")
#
#import db
#
rm(gsdb,dimdb,lhd,diseasedb)
attach("../RDatabase/.RData")
gsdb <- get("gsdb")
dimdb <- get("dimdb")
lhd <- get("lhd")
diseasedb <- get("diseasedb")
detach()
#
# Ajout donnees manuelles
#
latitude <- list(HM=54.11,RM=52.13)

#
#              Calcul sommes temperature
#              Fev 2011
#
#
#
# 1. Convertir les données en tableaux homogènes:
#************************************************
#
# Attendu : metdb, une liste de 4 elements (RM99,RM00,...), chaque element etant un dataframe : 
# DatePosix,Annee,Jour,HeureTU,Tair, allant du 1er septembre an n au 31 aout an n+1. 
# RM99 = du 1er sept 98 au 31 aout 99, etc...
# Verifier la convention 0h/24h pour minuit : le jour j demarre a 1h et se termine a 24h
#
# Pour Tair:
# - a partir de données horaires, on fait la moyenne (Tmin+Tmax)/2
# - a partir de donnees journaliere, on applique la fonction PartonLogan avec Tmin, Tmax et la latidude (cf input.R).
# - on bouche les trous avec approx
#
# Fonction pour t'inspirer dans cette etape :
#
# Mise en forme d'un fichier meteo Arvalis(datm) avec des dates posix, et un minuit avec une autre convention:
#
#mefMeteo <- function(datm) {
#  #print("coucou")
#  #lecture composant par composant 
#  a <- as.numeric(format(datm$Date,"%Y"))
#  j <- as.numeric(format(datm$Date,"%j"))
#  h <- as.numeric(format(datm$hhmm,"%H"))
#  m <- as.numeric(format(datm$hhmm,"%M"))
#  #arrondi a l'heure pres pour gestion lecture R 2.7 avec bug (1400 est lue 1359 !)
#  newh <- round(h+m/60)
#  j[h==23 & newh==24] <- j[h==23 & newh==24]-1
#  h <- newh
#  #homogeneisation des notation 0h en h=24
#  j[h==0] <- j[h==0]-1
#  h[h==0] <- 24
#  #cas jour 0
#  j[j==0] <- 365
#  p <- datm[,"Pluies(mm)"]
#  #homogeneisation filtrage
#  p <- ifelse(p==0.5,0.2,p)
#  dat <- data.frame(An=a,
#                    Jour=j,
#                    hhmm=h*100,
#                    PAR=round(datm[,"Rj (j/cm²)"]/3600/2.174*1e5*.48,2),
#                    Tair=datm[,"T (°C)"],
#                    HR=datm[,"HR(%)"],
#                    Vent=round(datm[,"Vent(Km/h)"]*1000/3600,3),
#                    Pluie=p)
#  #on vire minuit du jour d'avant
#  dat <- dat[-1,]
#  cbind(dat,Jsim=rep(seq(unique(dat$Jour)),rep(24,length(unique(dat$Jour))))[1:nrow(dat)])
#}
##
#meteodb <- lapply(datm,mefMeteo)
#  
##exemple bouchage trou
#approx(met$numJ+met$hhmm/100/24,met$Tair,met$numJ+met$hhmm/100/24)$y

###############################################

constrHMMeteo <- function(HMj, year) {
    
    HMjSelect <- subset(HMj, 
                        Date >= as.POSIXct(strptime(paste(year-1, "-09-01", sep=""),"%Y-%m-%d")) & Date <= as.POSIXct(strptime(paste(year, "-08-31", sep=""),"%Y-%m-%d")),
                        c(Date,Maximum,Minimum))
    
    annee <- rep(as.numeric(format(HMjSelect$Date,"%Y")),each=24)
    jour <- as.numeric(format(HMjSelect$Date,"%j"))
    
    pl <- PartonLogan(jour, HMjSelect$Minimum, HMjSelect$Maximum, latitude$HM)

    HM99 <- data.frame(datePosix=HMjSelect$Date,
                       annee=annee,
                       jour=pl$numJ,
                       heure=pl$hTU,
                       Tair=pl$T)
    
}


constrRM99Meteo <- function(RMj) {
    
    RMjSelect <- subset(RMj,
                        (year == 1998 & month >= 9 & day >= 1) | (year == 1999 & month <= 8 & day <= 31),
                        c(year,month,day,max,min))    
    
    datePosix <- as.POSIXct(strptime(paste(RMjSelect$year, "-", RMjSelect$month, "-", RMjSelect$day, "-", sep=""),"%Y-%m-%d"))
    datePosix <- rep(datePosix,each=24)
                
    pl <- PartonLogan(RMjSelect$day, RMjSelect$min, RMjSelect$max, latitude$RM)
    
    annee <- rep(RMjSelect$year,each=24)
    
    RM99 <- data.frame(datePosix=datePosix,
                       annee=annee,
                       jour=pl$numJ,
                       heure=pl$hTU,
                       Tair=pl$T)
               
}


constrRM00Meteo <- function(RMh) {
    
    RM00 <- NULL
    
    RMhSelect <- subset(RMh,
                        date >= as.POSIXct(strptime("1999-09-01","%Y-%m-%d")) & date <= as.POSIXct(strptime("2000-08-31","%Y-%m-%d")),
                        c(date,time,max.hourly.temp,min.hourly.temp))
    
    erroneousDates <- paste(RMhSelect$date, RMhSelect$time/100) # contains duplicates
    incompleteDates <- unique(erroneousDates) # after removing duplicates, 2 dates are missing
    
    elementsToRemove <- grep(2,table(erroneousDates))
    
    RMhSelect <- RMhSelect[-c(elementsToRemove),]
       
    # construct a complete dates sequence
    completeDates <- seq(as.POSIXct("1999-09-01 00:00:00", "GMT"), as.POSIXct("2000-08-31 23:00:00", "GMT"), by="1 hour")
    completeHours <- as.character(format(completeDates,"%H"))
    completeHours <- as.character(as.numeric(completeHours)+1)
    completeDates <- as.character(format(completeDates,"%Y-%m-%d"))
    completeDates <- paste(completeDates, completeHours)    
    
    elementToInsert <- grep(FALSE, is.element(completeDates, incompleteDates)) # it lacks "2000-08-25 1" data
   
    RMhSelect <- RMhSelect[c(1:8616, 8616, 8617:nrow(RMhSelect)),]
    RMhSelect[8617,] <- list(as.POSIXct("2000-08-25"),as.numeric("100"),as.numeric("14.1"),as.numeric("12.865"))
                  
    annee <- as.numeric(format(RMhSelect$date,"%Y"))
    jour <- as.numeric(format(RMhSelect$date,"%j"))
    RMhSelect$time <- as.integer(RMhSelect$time/100)
    Tair <- (RMhSelect$min.hourly.temp + RMhSelect$max.hourly.temp)/2
    
    RM00 <- data.frame(datePosix=RMhSelect$date,
                       annee=annee,
                       jour=jour,
                       heure=RMhSelect$time,
                       Tair=Tair)            
               
}


meteoHomogene <- function(meteodb) {
    HMj <- meteodb$HMj
    RMh <- meteodb$RMh
    RMj <- meteodb$RMj
    
    HM99 <- constrHMMeteo(HMj,1999)
    HM00 <- constrHMMeteo(HMj,2000)
    RM99 <- constrRM99Meteo(RMj)
    RM00 <- constrRM00Meteo(RMh)
    
    list(RM99=RM99, RM00=RM00, HM99=HM99, HM00=HM00)
}

metdb <- meteoHomogene(meteodb)


#
#
# 2. Calcul des sommes de temperature
#************************************
#
# 2.1 : Ajout, dans les tables metdb, de 2 colones : dsT et dsTc (= incréments de temps thermiques ou de temps thermique compensé) par heure
#

repTlin <- function(Th,Tb=0) max(0,Th-Tb) * 1 / 24

loiT <- function(TK, k, EaR, DSR, DHR) {
    # TK : temperature en Kelvin
    k * TK * exp(-EaR/TK) / (1 + exp(DSR - DHR / TK))
}

repTc <- function(Tair, TCref=12, k=3.8088e10, EaR=8899.6, DSR=68, DHR=20736) {
    TKref = 273 + TCref
    loiTref = loiT(TKref, k, EaR, DSR, DHR)
    sT=Tair
    for (i in 1:length(Tair)) {
        TK = 273 + Tair[i]
        sT[i] = loiT(TK, k, EaR, DSR, DHR) * TCref / loiTref * 1 / 24
    }
    sT
}

# repTc <- function(T,Ea....cf fonction C ci dessous)

#metdb <- lapply(metdb,function(dat) {
#            dat$dsT = sapply(dat$Tair,repTlin)
#            dat$dsTc = sapply(dat$Tair, repTc)}
#)

metdb$HM99$dsT <- sapply(metdb$HM99$Tair,repTlin)
metdb$HM00$dsT <- sapply(metdb$HM00$Tair,repTlin)
metdb$RM99$dsT <- sapply(metdb$RM99$Tair,repTlin)
metdb$RM00$dsT <- sapply(metdb$RM00$Tair,repTlin)
      
metdb$HM99$dsTc <- sapply(metdb$HM99$Tair, repTc)
metdb$HM00$dsTc <- sapply(metdb$HM00$Tair, repTc)
metdb$RM99$dsTc <- sapply(metdb$HM99$Tair, repTc)
metdb$RM00$dsTc <- sapply(metdb$RM00$Tair, repTc)

# pseudo code C pour repTc :
#double loiT(double TK,double k, double EaR, double DSR, double DHR) {
#    TK : temperature en Kelvin
#    #define TCref 12
#    #define k 3.8088e10
#    #define EaR 8899.6
#    #define DSR 68
#    #define DHR 20736
#
#    return(k * TK * exp(-EaR/TK) / (1 + exp(DSR - DHR / TK)));
#}

#float repTc(float Tair[],float TCref, double k, double EaR, double DSR, double DHR) {
#  double TKref = (double) (273. + TCref);
#  double TK; 
#  double loiTref = repTc(TKref,k,EaR,DSR,DHR);
#  double sT=Tair;
#  int i;
#  for (i = 0; i < nTair; i++) {
#    TK = (double) (273. + Th[i]);
#    sT[i] = loiT(TK,k,EaR,DSR,DHR) * TCref / loiTref * 1 / 24;
#  }
#  return sT;
#}

#
#            Bibliotheque modele meteo h a partir de meteoJ
#
#
#
#
# Fonction PartonLogan, pour estimer T horraire a partir de T journaliere
#
PartonLogan <- function(numJ,Tmin,Tmax,latitude = 55,param="air150cm") {
  # formule Parton & Logan, AFM 1981
  #
  #numJ, Tmin, Tmax : consecutive day time series  for julian day number, minimal and max temperature
  #latitude = latitude in degrees
  # param = named list of 3 coefficient (a,b,c) OR string indicating the type of temperature record ("air150cm","air10cm","soilsurface","soil-10cm")
  # a,b,c are for :
  # a = time lag coeff for  max Temperature after noon (h)
  # b = coefficient for temperature decrease at night
  # c = timeLag for the minimum temperature after sunrise (h)
  #
  # returns data.frame(hTU Tair)

# Select a,,b,c
paramref <- list("air150cm" = list(a = 1.86,b = 2.2, c = - 0.17),
                 "air10cm" = list(a = 1.52, b =  2, c = -0.18),
                 "soilsurface" = list(a = 0.5, b = 1.81, c = 0.49),
                 "soil-10cm" = list(a = 0.45, b = 2.28, c = 1.83))
                 
if (length(param) == 1)
  p <- paramref[[param]]
else
  p <- param
a <- p$a; b <- p$b; c <- p$c
#
  phi <- latitude / 180 * pi
  delta <- 0.4014 * sin( 2 * pi * (numJ - 77) / 365)
  t2 <- (-tan(phi)) * tan(delta)
  t1 <- sqrt(1 - t2^2)
  #amplitude angulaire du/des jour(s)
  ahr <- atan2(t1,t2)
  daylength <- ahr / pi * 24
  sunrise <- 12 - daylength / 2
  sunset <- 12 + daylength / 2
  # heure Tmin
  hmin <- sunrise + c
  # Temperature sunset
  Tsunset <- Tmin + (Tmax - Tmin) * sin (pi * (sunset - hmin) / (daylength + 2 * a))
  #Temperature sunset Jour d'avant
  Tsunsetb <- c(Tsunset[1],Tsunset[-length(Tsunset)])
  #Tmin jour d'apres
  Tmina <- c(Tmin[-1],Tmin[length(Tmin)])
  #
  h <- 1:24
  T <- sapply(seq(along=numJ),function(i)
              ifelse(h < sunrise[i],
                     Tmin[i] + (Tsunsetb[i] - Tmin[i]) * exp (-b * (24 - sunset[i] + h) / (24 - daylength[i])),
                     ifelse(h > sunset[i],
                            Tmina[i] + (Tsunset[i] - Tmina[i]) * exp (-b * (h - sunset[i]) / (24 - daylength[i])),
                            Tmin[i] + (Tmax[i] - Tmin[i]) * sin (pi * (h - hmin[i]) / (daylength[i] + 2 * a))
                            )
                     )
              ,simplify=FALSE)
  data.frame(numJ=rep(numJ,rep(24,length(numJ))),hTU=rep(h,length(numJ)),T=unlist(T))
}

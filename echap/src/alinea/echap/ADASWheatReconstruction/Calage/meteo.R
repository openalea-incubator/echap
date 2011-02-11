#
#              Calcul sommes temperature
#              Fev 2011
#
#
#
# (0). Rameuter (dans input.R) la meteodb:
#
# 1. Convertir les données en tableaux homogènes:
#************************************************
#
# Attendu : metdb, une liste de 4 elements (RM99,RM00,...), chaque element etant un dataframe : DatePosix,An,Jour,hhmm (heureTU),Tair, allant du 1er septembre an n au 31 aout an n+1. RM99 = du 1er sept 98 au 31 aout 99, etc...
#verifier la convention 0h/24h pour minuit : le jour j demarre a 1h et se termine a 24h
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
mefMeteo <- function(datm) {
  #print("coucou")
  #lecture composant par composant 
  a <- as.numeric(format(datm$Date,"%Y"))
  j <- as.numeric(format(datm$Date,"%j"))
  h <- as.numeric(format(datm$hhmm,"%H"))
  m <- as.numeric(format(datm$hhmm,"%M"))
  #arrondi a l'heure pres pour gestion lecture R 2.7 avec bug (1400 est lue 1359 !)
  newh <- round(h+m/60)
  j[h==23 & newh==24] <- j[h==23 & newh==24]-1
  h <- newh
  #homogeneisation des notation 0h en h=24
  j[h==0] <- j[h==0]-1
  h[h==0] <- 24
  #cas jour 0
  j[j==0] <- 365
  p <- datm[,"Pluies(mm)"]
  #homogeneisation filtrage
  p <- ifelse(p==0.5,0.2,p)
  dat <- data.frame(An=a,
                    Jour=j,
                    hhmm=h*100,
                    PAR=round(datm[,"Rj (j/cm²)"]/3600/2.174*1e5*.48,2),
                    Tair=datm[,"T (°C)"],
                    HR=datm[,"HR(%)"],
                    Vent=round(datm[,"Vent(Km/h)"]*1000/3600,3),
                    Pluie=p)
  #on vire minuit du jour d'avant
  dat <- dat[-1,]
  cbind(dat,Jsim=rep(seq(unique(dat$Jour)),rep(24,length(unique(dat$Jour))))[1:nrow(dat)])
}
#
meteodb <- lapply(datm,mefMeteo)
  
              #exemple bouchage trou
approx(met$numJ+met$hhmm/100/24,met$Tair,met$numJ+met$hhmm/100/24)$y
#
#
# 2. Calcul des sommes de temperature
#************************************
#
# 2.1 : Ajout, dans les tables metdb, de 2 colones : dsT et dsTc (= incréments de temps thermiques ou de temps thermique compensé) par heure
#
#
repTlin <- function(Th,Tb=0) max(0,Th-Tb) * 1 / 24
#
repTc <- function(T,Ea....cf fonction C ci dessous)
#
# on doit pouvoir alors faire directement :
#
metdb <- lapply(metdb,function(dat) {
  dat$dsT = sapply(dat$Tair,repTlin)
  dat$dsTc = sapply(dat$Tair, repTc)}
                )
# pseudo code C pour repTc :
double loiT(double TK,double k, double EaR, double DSR, double DHR) {
  TK : temperature en Kelvin
  #define TCref 12
#define k 3.8088e10
#define EaR 8899.6
#define DSR 68
#define DHR 20736

  return(k * TK * exp(-EaR/TK) / (1 + exp(DSR - DHR / TK)));
    }

float repTc(float Tair[],float TCref, double k, double EaR, double DSR, double DHR) {
  double TKref = (double) (273. + TCref);
  double TK; 
  double loiTref = repTc(TKref,k,EaR,DSR,DHR);
  double sT=Tair;
  int i;
  for (i = 0; i < nTair; i++) {
    TK = (double) (273. + Th[i]);
    sT[i] = loiT(TK,k,EaR,DSR,DHR) * TCref / loiTref * 1 / 24;
  }
  return sT;
}

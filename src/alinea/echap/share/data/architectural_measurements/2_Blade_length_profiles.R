#
#
#                          Preprocessing 2: estimate blade length profile
#                          (per nff and for mean plant)
#
#
# Estimate ligulation times per plant
#
# pentes : groupes par nff
phenM <- lapply(phendb, function(dat) {res <- aggregate(dat$HS, list(nff=dat$nff,TT=dat$TT), mean, na.rm=TRUE);res$HS=res$x;res})
fitM <- lapply(phenM, function(dat) lapply(split(dat,dat$nff), hsfit))
slopes <- lapply(fitM,function(fit) lapply(fit, function(f) f$coeff[2]))
# TTo per nff
TTo <- lapply(fitM,function(fit) lapply(fit, function(f) -f$coeff[1]/f$coeff[2]))
# TT_ligulation per leaf per plant
TTem <- sapply(genos,function(g) do.call('rbind',lapply(split(phendb[[g]],phendb[[g]]$N), function(dat) TTleaf(dat,slopes[[g]][[as.character(dat$nff[1])]]))), simplify=FALSE)
#
# Phyllochones moyens
#
phenMM <- lapply(phendb, function(dat) aggregate(dat[,c('nff','HS')], list(TT=dat$TT), mean, na.rm=TRUE))
fitMM <- lapply(phenM, hsfit)
slopeM <- lapply(fitMM,function(fit) fit$coeff[2])
TToM <- lapply(fitMM,function(fit) -fit$coeff[1]/fit$coeff[2])
#
#Construct Lblade = f(TTem) profiles
#
Lblade <- sapply(genos, function(g) merge(dimP[[g]], TTem[[g]]), simplify=FALSE)
#
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(Lblade, function(dim) {
  plot(c(0,1600),c(0,30),type='n')
  lapply(split(dim,dim$N), function(x) points(x$TTlig,x$L,col=x$N,pch=16))
  lines(smooth.spline(dim$TTlig,dim$L,df=8))
})
#
# Raw fit : shape bof
#
# final Strategy: use Tremie 13 shape, align data along x using position of the max + divide by phyllochronic difference, then linear fit of yref/yobs with position to get smooth deformation of original profile that fit the data
#
# demo effet normalisation x
#
par(mfrow=c(1,1))
plot(c(-10,5),c(0,1.2),type='n')
sapply(seq(Lblade), function(ig) {
  dim <- Lblade[[ig]]
  a <- slopeM[[ig]]
  spl <- smooth.spline(dim$TTlig,dim$L,df=8)
  ymax <- max(spl$y)
  xmax <- spl$x[which.max(spl$y)]
  lapply(split(dim,dim$N), function(x) points((x$TTlig - xmax) * a ,x$L/ymax,col=ig,pch=16))
  lines(smooth.spline((dim$TTlig - xmax) * a,dim$L/ymax,df=6),col=ig)
})
#
# new fit
#
# add destructive data + Leaf 1 of Tremie13 for completing first ranks of Tremie12
#
scans12 <- scandb$Tremie12
scans12 <- scans12[scans12$id_Axe=='MB' & scans12$stat < 3,1:9]
Lb12 <- data.frame(N=NA,nff=NA,rank=scans12$rank, L=scans12$lmax, TTlig=TToM$Tremie12 + scans12$rank / slopeM$Tremie12)
#
par(mfrow=c(1,1))
plot(c(0,1600),c(0,30),type='n')
points(Lblade$Tremie12$TTlig, Lblade$Tremie12$L, col=2, pch=16)
points(Lb12$TTlig, Lb12$L, col=1,pch=16)
# keep rank 5,6,7,8
Lb12 <- Lb12[Lb12$rank < 8,]
#
Lbladec <- Lblade
Lbladec$Tremie12 <- rbind(Lblade$Tremie12, Lb12)
spldb <- lapply(Lbladec, function(dim) smooth.spline(dim$TTlig,dim$L,df=8))
splref <- spldb$Tremie13
aref <- slopeM$Tremie13
def <- c(TRUE, TRUE, TRUE, FALSE)
names(def) <- genos
fitLdb <- sapply(genos, function(g) fitLspl(spldb[[g]], slopeM[[g]], splref, aref, deform=def[g]), simplify=FALSE)
#
par(mfrow=c(2,2),mar=c(4,4,1,1))
sapply(genos, function(g) {
  dim <- Lblade[[g]]
  plot(c(0,1600),c(0,30),type='n')
  lapply(split(dim,dim$N), function(x) points(x$TTlig,x$L,col=x$N,pch=16))
  if (g=='Tremie12')
    points(Lb12$TTlig, Lb12$L, pch=16, col=8)
  lines(fitLdb[[g]]$x, fitLdb[[g]]$y,col=2)
})
#
# Profils moyen par nff et globaux
#
LbM <- sapply(genos, function(g) {
  fitp <- fitM[[g]]
  fitL <- fitLdb[[g]]
  res <- fitp
  for (nff in names(fitp)) {
    ranks = seq(as.numeric(nff))
    TT = (ranks - fitp[[nff]]$coeff[1]) / fitp[[nff]]$coeff[2]
    res[[nff]] <- data.frame(ranks = ranks, L=approx(fitL$x, fitL$y,TT)$y)
  }
  res}, simplify=FALSE)
#
LbMM <- sapply(genos, function(g) {
  fitp <- fitMM[[g]]
  dim <- dimPMM[[g]]
  fitL <- fitLdb[[g]]
  nff <- max(dim$rank)
  ranks = seq(as.numeric(nff))
  TT = (ranks - fitp$coeff[1]) / fitp$coeff[2]
  data.frame(ranks = ranks, L=approx(fitL$x, fitL$y,TT)$y)
  }, simplify=FALSE)
#
#
par(mfrow=c(2,2),mar=c(4,4,1,1))
sapply(genos, function(g) {
  dim <- dimP[[g]]
  Lb <- LbM[[g]]
  plot(c(0,15),c(0,30),type='n')
  nffs <- names(Lb)
  for (i in seq(nffs)) {
    m <- Lb[[nffs[i]]]
    d <- dim[dim$nff==as.numeric(nffs[i]),]
    if (nrow(d) > 0)
      lapply(split(d,d$N), function(x) points(x$rank,x$L,col=i,pch=16,type='p'))
    lines(m$rank,m$L,col=i)
    lines(LbMM[[g]]$rank,LbMM[[g]]$L,col=4,lty=2)
  }
})

#
#
#                          Preprocessing 2: estimate blade length profile =f(TTem)
#
#
# Estimate ligulation times per plant
#
# pentes : groupes par nff
phenM <- lapply(phendb, function(dat) {res <- aggregate(dat$HS, list(nff=dat$nff,TT=dat$TT), mean, na.rm=TRUE);res$HS=res$x;res})phyl <- sapply(split(dat,dat$N), function(x) ifelse(nrow(x) >=2, lsfit(x$TT,x$HS / x$nff)$coeff[2],NA))
fitM <- lapply(phenM, function(dat) lapply(split(dat,dat$nff), hsfit))
slopes <- lapply(fitM,function(fit) lapply(fit, function(f) f$coeff[2]))
TTo <- lapply(fitM,function(fit) mean(sapply(fit, function(f) -f$coeff[1]/f$coeff[2])))
# TT_ligulation per leaf per plant
TTem <- sapply(genos,function(g) do.call('rbind',lapply(split(phendb[[g]],phendb[[g]]$N), function(dat) TTleaf(dat,slopes[[g]][[as.character(dat$nff[1])]]))), simplify=FALSE)
#
Lblade <- sapply(genos, function(g) merge(dimP[[g]], TTem[[g]]), simplify=FALSE)
#
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(Lblade, function(dim) {
  plot(c(0,1600),c(0,30),type='n')
  lapply(split(dim,dim$N), function(x) points(x$TTlig,x$L,col=x$N,pch=16))
  lines(smooth.spline(dim$TTlig,dim$L,df=6))
})
#
par(mfrow=c(1,1))
plot(c(-10,5),c(0,1.2),type='n')
sapply(seq(Lblade), function(ig) {
  dim <- Lblade[[ig]]
  tto <- TTo[[ig]]
  a <- mean(unlist(slopes[[ig]]))
  spl <- smooth.spline((dim$TTlig-tto)*a,dim$L,df=6)
  ymax <- max(spl$y)
  xmax <- spl$x[which.max(spl$y)]
  lapply(split(dim,dim$N), function(x) points((x$TTlig-tto)*a - xmax,x$L/ymax,col=ig,pch=16))
})

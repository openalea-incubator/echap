#
#      shared function for echap pipeline
#
# xy -> sphi transform
#
sphi_shapes <- function(sh,df=3,ds=0.02) {
  smid <- (sh$s[-nrow(sh)] + sh$s[-1]) / 2
  phi <- atan2(diff(sh$yr),diff(sh$xr))
  # HERE USE LINEAR INTERPOLATION TO BE CONSISTENT WITH TINO
  #fit <- smooth.spline(smid,phi,df=df)
  fit <- splinefun(smid,phi)
  news <- seq(0,1,ds)
  newsmid <- (news[-length(news)] + news[-1]) / 2 * max(sh$s)
  #newphi <- predict(fit,newsmid)$y
  newphi <- fit(newsmid)
  na.omit(data.frame(inerv=sh$inerv[1],phiS=sh$phiS[1],smid=newsmid, ds=ds * max(sh$s), phi = newphi, smid_n = newsmid / max(sh$s), ds_n = ds))
}
#
# transform from s,phi to x,y
#
sphi2xy <- function(ds,phi) {
  dx <- ds*cos(phi)
  dy <- ds*sin(phi)
  data.frame(x=c(0,cumsum(dx)), y=c(0,cumsum(dy)), s= c(0,cumsum(ds)))
}
#
# get descriptors
#
get_desc <- function(fit) {
  #print(fit$inerv[1])
  maxs <- (max(fit$smid) + fit$ds[1] / 2)
  xy <- sphi2xy(fit$ds,fit$phi)
  imax <- min(which.max(xy$y))
  Pa <- xy$s[imax] / max(xy$s)
  rs <- fit$smid / maxs
  phiS = fit$phiS[1]
  phiL <- abs(90 - mean(fit$phi[rs < 0.1]) / pi * 180)
  phi0 <- phiS + phiL
  C <- sum(diff(fit$phi)) /pi * 180
  Cabs <- sum(abs(diff(fit$phi))) / pi * 180
  proj <- max(xy$x) / maxs
  hproj <- (max(xy$x) - min(xy$x)) / maxs
  vproj <- (max(xy$y) - min(xy$y)) /maxs
  ymax <- max(xy$y) / maxs
  ymin <- min(xy$y) / maxs
  last <- nrow(xy)
  xtip <- xy$x[last] / maxs
  ytip <- xy$y[last] / maxs
  data.frame(inerv=fit$inerv[1],phiS=phiS,phi0=phi0,phiL=phiL,C=C,Cabs=Cabs,smax=maxs, Pa=Pa,proj=proj, ymax=ymax,ymin=ymin,xtip=xtip,ytip=ytip,hproj=hproj,vproj=vproj)
}
#
# creating shape descriptor tables
#
get_shapes <- function(sphi_fits,blades,topleaves=4) {
  desc <- lapply(sphi_fits,get_desc)
  desc <- do.call("rbind",desc)
  shapes <- blades
  shapes <- merge(shapes,desc)
  shapes$pos <- ifelse(shapes$ranktop <=topleaves, 2,1)
  shapes$psen <- 1 - shapes$Agr/shapes$A
  # grouping is rank,variety, harvest for lower leaves and ranktop, variety, harvest for upper leaves
  shapes$rankg <- ifelse(shapes$ranktop <=topleaves, shapes$ranktop, shapes$rank)
  shapes
}
#
# fit interpolated x(s,age),y(s,age)
#
sxy_interpolation <- function(sxy,shapes,at=seq(0,12,0.5),nzero=1000) {
  xymed <- NULL
  s <- do.call('cbind',lapply(sxy,function(fit) fit$s))
  x <- do.call('cbind',lapply(sxy,function(fit) fit$x))
  y <- do.call('cbind',lapply(sxy,function(fit) fit$y))
  #
  smed <- apply(s,1,median)
  #
  xmed <- try(t(apply(x,1,function(xage) predict(smooth.spline(shapes$age,xage,df=3),at)$y)),TRUE)
  ymed <- try(t(apply(y,1,function(yage) predict(smooth.spline(shapes$age,yage,df=3),at)$y)),TRUE)
  hproj_med <- try(predict(smooth.spline(shapes$age,shapes$hproj,df=4),at)$y, TRUE)
  vproj_med <- try(predict(smooth.spline(shapes$age,shapes$vproj,df=4),at)$y, TRUE)
  #
  if (!inherits(xmed, "try-error") & !inherits(ymed, "try-error"))
    xymed <- sapply(seq(ncol(xmed)),function(icol) {
      x <- xmed[,icol]
      y <- ymed[,icol]
      #x <- x / (max(x) - min(x)) * hproj_med[icol]
      #y <- y / (max(y) - min(y)) * vproj_med[icol]
      ds <- sqrt(diff(x)^2+diff(y)^2)
      phi <- atan2(diff(y),diff((x)))
      dsn <- ds / sum(ds)
      dx <- dsn*cos(phi)
      dy <- dsn*sin(phi)
      s <- c(0,cumsum(dsn))
      data.frame(x=c(0,cumsum(dx)), y=c(0,cumsum(dy)), s= s)
      },simplify=FALSE)
  xymed
}
#
# sphi interpolation
#
dphi_interpolation <- function(sphi,shapes,at=seq(0,12,0.5)) {
  xymed <- NULL
  phi <- do.call('cbind',lapply(sphi,function(fit) fit$phi))
  ds <-  do.call('cbind',lapply(sphi,function(fit) fit$ds_n))
  dsmed <- apply(ds,1,median)
  #
  # initial median shapes : force dphi to continuous decrease
  phi_i <- phi[1,]
  phi_i_med <- try(predict(smooth.spline(shapes$age,phi_i,df=3),at)$y, TRUE)
  dphi_abs<- apply(phi,2,function(x) x - x[1])
  dphi_abs_med <- try(t(apply(dphi_abs,1,function(x) predict(smooth.spline(shapes$age,x,df=3),at)$y)),TRUE)
  dphi_med <- apply(dphi_abs_med,2,function(x) pmin(0,diff(x)))
  if (!inherits(dphi_abs_med, "try-error") & !inherits(phi_i_med, "try-error")) {
    phi_med <- apply(rbind(phi_i_med,dphi_med),2,cumsum)
    xymed <- sapply(seq(ncol(phi_med)),function(icol) sphi2xy(dsmed,phi_med[,icol]),simplify=FALSE)
  }
  xymed
}
#
# deprecated multicriteria fit
#

#
# fit phi_at_s=f(age) to get smooth interpolations at regular intervals.
#
# Filter data age >5 for lower leaves
# get data:
#shapeg <- split(shapedb, g)[['Tremie']]
#sphig <- split(sphidb, g)[['Tremie']]
#sphi <- split(sphig, shapeg$pos)[[2]]
#shapes <- split(shapeg, shapeg$pos)[[2]]
dphi_interpolation <- function(sphi,shapes,at=seq(0,12,0.5)) {
  xymed <- NULL
  phi <- do.call('cbind',lapply(sphi,function(fit) fit$phi))
  ds <-  do.call('cbind',lapply(sphi,function(fit) fit$ds_n))
  dsmed <- apply(ds,1,median)
  #
  # initial median shapes : force dphi to continuous decrease
  phi_i <- phi[1,]
  phi_i_med <- try(predict(smooth.spline(shapes$age,phi_i,df=3),at)$y, TRUE)
  dphi_abs<- apply(phi,2,function(x) x - x[1])
  dphi_abs_med <- try(t(apply(dphi_abs,1,function(x) predict(smooth.spline(shapes$age,x,df=3),at)$y)),TRUE)
  dphi_med <- apply(dphi_abs_med,2,function(x) pmin(0,diff(x)))
  #
  # targets
  xtip_med <-  try(predict(smooth.spline(shapes$age,shapes$xtip,df=4),at)$y, TRUE)
  ytip_med <-  try(predict(smooth.spline(shapes$age,shapes$ytip,df=4),at)$y, TRUE)
  proj_med <- try(predict(smooth.spline(shapes$age,shapes$proj,df=4),at)$y, TRUE)
  ymax_med <- try(predict(smooth.spline(shapes$age,shapes$ymax,df=4),at)$y, TRUE)
  ymin_med <- try(predict(smooth.spline(shapes$age,shapes$ymin,df=4),at)$y, TRUE)
  if (!inherits(dphi_abs_med, "try-error") & !inherits(phi_i_med, "try-error") & !inherits(proj_med, "try-error") & !inherits(ymin_med, "try-error") & !inherits(ymax_med, "try-error")) {
    #original : use phi__i median as base angle
    #phi_med <- apply(rbind(phi_i_med,dphi_med),2,cumsum)
    #xymed <- sapply(seq(ncol(phi_med)),function(icol) sphi2xy(dsmed,phi_med[,icol]),simplify=FALSE)
    #alternative : choose phi so that proj is okay
    #
    #adjust dphi
    dphi_med_opt <- dphi_med
    for (icol in seq(ncol(dphi_med))) {
      targets <- c(proj_med[icol],ymax_med[icol],ymin_med[icol])
      targets <- c(xtip_med[icol],ytip_med[icol])
      scfactor <- seq(0.7,1.3,length.out=50)
      crit <- sapply(scfactor, function(sc) {
        phim <- cumsum(c(phi_i_med[icol],dphi_med[,icol]*sc))
        xy <- sphi2xy(dsmed,phim)
        proj <- max(xy$x) / max(xy$s)
        ymax <- max(xy$y) / max(xy$s)
        ymin <- min(xy$y) / max(xy$s)
        xtip <- xy$x[nrow(xy)] / max(xy$s)
        ytip <- xy$y[nrow(xy)] / max(xy$s)
        c(proj,ymax,ymin)
        c(xtip,ytip)
      },simplify=FALSE)
      scores <- sapply(crit,function(x) sum(abs(x - targets)))
      best <- scfactor[which.min(scores)]
      dphi_med[,icol] <- best * dphi_med[,icol]
    }
    #adjust phi__i
    phi_i_proj <- sapply(seq(ncol(dphi_med)),function(icol) {
      targets <- c(proj_med[icol],ymax_med[icol],ymin_med[icol])
      targets <- c(xtip_med[icol],ytip_med[icol])
      cphi <- seq(0.7* abs(phi_i_med[icol]),1.3 * abs(phi_i_med[icol]),length.out=50)
      crit <- sapply(cphi, function(phi) {
        phim <- cumsum(c(phi,dphi_med[,icol]))
        xy <- sphi2xy(dsmed,phim)
        proj <- max(xy$x) / max(xy$s)
        ymax <- max(xy$y) / max(xy$s)
        ymin <- min(xy$y) / max(xy$s)
        xtip <- xy$x[nrow(xy)] / max(xy$s)
        ytip <- xy$y[nrow(xy)] / max(xy$s)
        c(proj,ymax,ymin)
         c(xtip,ytip)
      },simplify=FALSE)
      scores <- sapply(crit,function(x) sum(abs(x - targets)))
      best <- cphi[which.min(scores)]
      best
    })
    phi_med <- apply(rbind(phi_i_proj,dphi_med),2,cumsum)
    xymed <- sapply(seq(ncol(phi_med)),function(icol) sphi2xy(dsmed,phi_med[,icol]),simplify=FALSE)
  }
  xymed
}

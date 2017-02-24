#
#            Compute median leaf interpolated trajectories for leaf shape=f(age)
#
#
source('funShapes.R')
#
bl_e <-  read.csv('leafshape_metainfo_Boigneville_2012_2013.csv', sep=',', dec='.')
nervs_e <-  read.csv('leafshape_xydata_Boigneville_2012_2013.csv', sep=',', dec='.')
bl_ij <-  read.csv('leafshape_metainfo_Boigneville_2011.csv', sep=',', dec='.')
nervs_ij <-  read.csv('leafshape_xydata_Boigneville_2011.csv', sep=',', dec='.')
bl <-  read.csv('leafshape_metainfo_Grignon2010.csv', sep=',', dec='.')
nervs <-  read.csv('leafshape_xydata_Grignon2010.csv', sep=',', dec='.')
# merge indices
imax <- max(bl$inerv)
bl_e$inerv <- bl_e$inerv + imax
nervs_e$inerv <- nervs_e$inerv + imax
#
imax <- max(bl_e$inerv)
bl_ij$inerv <- bl_ij$inerv + imax
bl_ij$label <- paste(bl_ij$label,'11',sep='')
nervs_ij$inerv <- nervs_ij$inerv + imax
#merge
nervs_e <- split(nervs_e,nervs_e$inerv)
nervs_ij <- split(nervs_ij,nervs_ij$inerv)
nervs <- split(nervs,nervs$inerv)
nervsdb <-c(nervs,nervs_e,nervs_ij)
#
cols <- c('inerv','label','Source','plant','rank','ranktop','age','lmax','wmax','A','Agr','mass','stat','phiP')
bldb <- rbind(bl[,cols],bl_e[,cols], bl_ij[,cols])
bldb$variety = as.numeric(as.factor(bldb$label))
bldb$harvest = as.numeric(as.factor(bldb$Source))
#
# sphi spline fit
#
sphidb <- lapply(nervsdb, sphi_shapes,df=10)
#normalised leaves
sxydb <- lapply(sphidb,function(fit) sphi2xy(fit$ds_n,fit$phi))
#
# check fits
#
par(mfrow=c(4,5),mar=c(1,1,1,1))
ids <- sample(seq(nervsdb),20)
mapply(function(dat,fit) {plot(dat$xr,dat$yr);lines(sphi2xy(fit$ds,fit$phi),col=2)}, nervsdb[ids],sphidb[ids])
#
# shapes descriptors and groups
#
shapedb <- get_shapes(sphidb, bldb,topleaves=4)
groups <- c('MerciaRht','MerciaRht','MerciaRht', 'MerciaRht','Soissons','Tremie', 'Tremie','MerciaRht','MerciaRht')
names(groups) <- c('Mercia','Rht2','Rht10','Rht3','Soissons','Tremie12','Tremie13','Mercia11','Rht311')
shapedb$group <- groups[as.character(shapedb$label)]
litpos <- ifelse(shapedb$pos > 1,'upper','lower')
shapedb$leafgroup <- paste(shapedb$group,litpos,sep='_')

  median_leaf <- function(sxy,inervs) {
    x <- do.call('cbind',lapply(sxy,function(x) x$x))
     y <- do.call('cbind',lapply(sxy,function(x) x$y))
    xmed <- apply(x,1,median)
    ymed <- apply(y,1,median)
    dx <- x - matrix(xmed,nrow=nrow(x),ncol=ncol(x),byrow=FALSE)
    dy <- y - matrix(ymed,nrow=nrow(y),ncol=ncol(y),byrow=FALSE)
    err <- abs(dx) + abs(dy)
    inervs[which.min(apply(err,2,sum))]
  }
#
# compute median leaf per ntop and per source
#
groups <- split(shapedb,list(shapedb$Source, shapedb$label), drop=TRUE)
#
fits <- lapply(groups,function(sh) {
  lapply(split(sh,sh$ranktop), function(x) {
  sxy <- sxydb[x$inerv]
  ileaf <- median_leaf(sxy, x$inerv)
  sxydb[[ileaf]]})
})
# median age of leaves in fits
medage <- lapply(groups,function(sh) {
  lapply(split(sh,sh$ranktop), function(x) median(x$age))
})
# horizontal projection
medhp <- lapply(fits,function(fit) {
  lapply(fit, function(x) max(x$x) - min(x$x))
})



#
# median leaf reduction
# (keep only median leaves for homogeneous ageclass)
#
require(zoo)
groups <- split(shapedb, shapedb$leafgroup)
sh <- groups[[1]]



agedisp <- function(sh, nmed=7) {
  sh <- sh[order(sh$age),]
  sxy <- sxydb[sh$inerv]
 #elementary ageclass
  ages <- seq(0,13,0.1)
 # number of leaves to form a median
  agemin <- rollapply(sh$age,nmed,min)
  agemax <- rollapply(sh$age,nmed,max)
  agemax-agemin
}
#
medagedisp <- sapply(groups, function(x) agedisp(x,nmed=5))
par(mfrow=c(3,3),mar=c(4,4,1,1))
mapply(function(x,name) hist(x,main=name),medagedisp,sapply(groups, function(x) x$label[1]))
#
# conc : take O.15 max dispersion class for median leaves
#
leafzero <- sxy[[1]]
leafzero$x <- 0
leafzero$y <- leafzero$s
#
shzero <- sh[1,]
shzero$age <- 0
shzero$hproj <- 0
shzero$vproj <- 1
shzero$inerv <- -1


median_shapes <- function(sh, nmed=9, maxdisp=0.2) {
}

par(mfrow=c(6,3),mar=c(4,4,1,1))
for (g in 1:6) {
  sh <- groups[[g]]
  nmed=7
  
  sh <- sh[order(sh$age),]
  sxy <- sxydb[sh$inerv]
  
  agemin <- rollapply(sh$age,nmed,min)
  agemax <- rollapply(sh$age,nmed,max)
  agemean <- rollapply(sh$age,nmed,mean)
  disp <- agemax - agemin
  maxdisp <- round(median(disp),1)
  sel <- disp <= maxdisp
                                        #
  medleaf <- rollapply(seq(sxy), nmed, function(isel) median_leaf(sxy[isel], sh$inerv[isel]))
                                        #xmed <- rollapply(seq(sxy), nmed, function(isel) apply(do.call('cbind',lapply(sxy[isel],function(x) x$x)),1,median))
                                        #ymed <- rollapply(seq(sxy), nmed, function(isel) apply(do.call('cbind',lapply(sxy[isel],function(x) x$y)),1,median))
    # smed <- rollapply(seq(sxy), nmed, function(isel) apply(do.call('cbind',lapply(sxy[isel],function(x) x$s)),1,median))
                                        #
  age <- agemean[sel]
  #xmed <- xmed[sel,]
                                        #ymed <- ymed[sel,]
                                        #smed <- smed[sel,]
                                        #
  medleaves <- sxydb[medleaf[sel]]
  medsh <- shapedb[medleaf[sel],]
                                        #
  ages <- seq(0,max(age),0.5)
#add leaf zero
  medz <- c(lapply(seq(nmed), function(x) leafzero),medleaves)
  shz <- rbind(do.call('rbind',lapply(seq(nmed), function(x) shzero)),medsh)
  dxymed <- sxy_interpolation(medz,shz,at=ages)
 # dxymed <- sxy_interpolation(medleaves,medshsh,at=ages)
  dxymedraw <- sxy_interpolation(sxy,sh,at=ages)
  jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  maxage <- 6
  agerange <- seq(0,maxage,0.5)
  palette(jet(length(agerange)))
  plot(c(-.3,1), c(-.3,1),type='n',xlab=names(groups)[g])
  for (i in seq(age)) {
    sxymed <- medleaves[[i]]
    lines(sxymed$x, sxymed$y, col=findInterval(age[i],agerange))
    text(sxymed$x[51], sxymed$y[51], round(age[i],1),cex=0.6)
  }

plot(c(-.3,1), c(-.3,1),type='n')
sapply(seq(dxymed),function(j) {lines(dxymed[[j]],col=findInterval(ages[j],agerange))})


hp <- sapply(medleaves,function(leaf) max(leaf$x) -min(leaf$x))
hpr <- sapply(sxy,function(leaf) max(leaf$x) -min(leaf$x))
  palette('default')
  plot(age,hp,col=2,pch=16,ylim=c(0,1),xlim=c(0,maxage))
  lines(smooth.spline(age,hp,df=6),col=2)
  points(sh$age,hpr,pch=16,col=c(8,1)[sh$pos])
  lines(smooth.spline(sh$age,hpr,df=6),col=1)
  hpm <- sapply(dxymed, function(leaf) max(leaf$x) -min(leaf$x))
  points(ages,hpm,col=3,pch=16)
  lines(smooth.spline(ages,hpm,df=6),col=3)
  points(ages,sapply(dxymedraw, function(leaf) max(leaf$x) -min(leaf$x)),col=4,pch=16)
}
#
  #ageclass <- cut(sh$age, ages)

#
# trajectories
#
ages <- seq(-2,12,0.5)
#shapedb$age <- pmax(min(ages)+0.01,shapedb$age)
#

#
# All ranks
#
dxymedM <- sxy_interpolation(sxydb,shapedb,at=ages)
dxymedgM <- mapply(sxy_interpolation,split(sxydb, shapedb$group), split(shapedb,shapedb$group), MoreArgs=list(at=ages),SIMPLIFY=FALSE)
dxymedlabM <- mapply(sxy_interpolation,split(sxydb, shapedb$label), split(shapedb,shapedb$label), MoreArgs=list(at=ages),SIMPLIFY=FALSE)
# subset T1 fit


# per leaf pos
dxymed <- mapply(sxy_interpolation,split(sxydb, shapedb$pos), split(shapedb,shapedb$pos), MoreArgs=list(at=ages),SIMPLIFY=FALSE)
dxymedg <- mapply(function(sxyg, shapeg) mapply(sxy_interpolation,split(sxyg, shapeg$pos), split(shapeg, shapeg$pos), MoreArgs=list(at=ages), SIMPLIFY=FALSE),
                  split(sxydb,shapedb$group), split(shapedb, shapedb$group),SIMPLIFY=FALSE)
#
#
# plot trajectories
#
jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
palette(jet(length(ages)))
par(mfcol=c(3,5), mar=c(4,4,1,1))
#
# Pooled per pos
leg <- c('lower leaves','upper leaves')
for (i in seq(dxymed)) {
  plot(c(0,1.3),c(-0.3,1),type='n',xlab=leg[i],ylab='interpolated median leaf',cex.lab=1.5)
  sapply(seq(dxymed[[i]]),function(j) {lines(dxymed[[i]][[j]],col=j)})
  if (i==2)
    text(1,0.8,'Pooled')
}
# poooled all ranks
plot(c(0,1.3),c(-0.3,1),type='n',xlab='All leaves',ylab='interpolated median leaf',cex.lab=1.5)
sapply(seq(dxymedM),function(j) {lines(dxymedM[[j]],col=j)})
# per genotype group per pos
for (g in names(dxymedg)) {
  dxy <- dxymedg[[g]]
  for (i in seq(dxy)) {
    plot(c(0,1.3),c(-0.3,1),type='n',xlab=leg[i],ylab='interpolated median leaf',cex.lab=1.5)
    sapply(seq(dxy[[i]]),function(j) {lines(dxy[[i]][[j]],col=j)})
    if (i==2)
    text(1,0.8,g)
  }
  # all ranks per grenotype group
  plot(c(0,1.3),c(-0.3,1),type='n',xlab="All leaves",ylab='interpolated median leaf',cex.lab=1.5)
  sapply(seq(dxymedgM[[g]]),function(j) {lines(dxymedgM[[g]][[j]],col=j)})
}
#
# fits per label
#
jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
palette(jet(length(ages)))
par(mfrow=c(3,3), mar=c(4,4,1,1))
for (g in names(dxymedlabM)) {
  dxy <- dxymedlabM[[g]]
  plot(c(0,1.3),c(-0.3,1),type='n',xlab=leg[i],ylab='interpolated median leaf',cex.lab=1.5)
  sapply(seq(dxy),function(j) {lines(dxy[[j]],col=j)})
  text(1,0.8,g)
}
#
#
# Check relation with age / median trajectorie projection
#
shapeT <- shapedb[shapedb$group=='MerciaRht11',c('inerv','label','Source','rank','ranktop','age','pos','hproj','vproj')]
nervT <- do.call('rbind',split(nervsdb,shapedb$group)[['MerciaRht11']])[,c('inerv','s','xr','yr')]
nervp <- merge(nervT,shapeT)
bins <- c(-1,ages)
nervp$ageclass <- cut(nervp$age, bins)
symb <- c(1,16)
palette(jet(length(ages)))
plot(nervp$age,nervp$hproj,col=as.numeric(nervp$ageclass),pch=symb[nervp$pos])
palette('default')
plot(nervp$age,nervp$hproj,col=as.numeric(as.factor(nervp$Source)),pch=16)
t13 <- nervp[nervp$label=='Tremie13',]
palette('default')
points(t13$age,t13$hproj,col=1,pch='^')
lapply(split(nervp,nervp$pos),function(x) lines(smooth.spline(x$age,x$hproj,df=4),col=8))
#
medleaf <- trajfit[['Tremie']]
medp <- do.call('rbind',lapply(split(medleaf,list(medleaf$age,medleaf$lindex),drop=TRUE),function(leaf) data.frame(age=leaf$age[1],lindex=leaf$lindex[1],proj=max(leaf$x))))
lapply(split(medp,medp$lindex),function(x) lines(x$age,x$proj,col=2))
#
# export
#
loc <- c('Grignon2010','Grignon2010','Boigneville2012_2013')
names(loc) <- c('MerciaRht','Soissons','Tremie')
#by leafclass
trajfit <- lapply(dxymedg, function(dxy) {
  xya <- lapply(dxy,function(traj) do.call('rbind',mapply(function(xys,age) data.frame(age=age, x=xys$x,y=xys$y), traj,ages,SIMPLIFY=FALSE)))
  xyai <- do.call('rbind',mapply(function(xy,index) {xy$lindex=index;xy},xya,as.numeric(names(xya)),SIMPLIFY=FALSE))
})
mapply(function(fit,label) {
  fn <- paste('median_leaf_trajectories_',label, '_',loc[label],'_byleafclass.csv',sep='')
  write.csv(fit, fn, row.names=FALSE)
  },trajfit, names(trajfit))
# by genotypic group only
trajfit <- lapply(dxymedgM, function(dxy) {
  xya <- do.call('rbind',mapply(function(xys,age) data.frame(age=age, x=xys$x,y=xys$y), dxy,ages,SIMPLIFY=FALSE))
  xya$lindex=1
  xya})
mapply(function(fit,label) {
  fn <- paste('median_leaf_trajectories_',label, '_',loc[label],'.csv',sep='')
  write.csv(fit, fn, row.names=FALSE)
  },trajfit, names(trajfit))
# by label
trajfit <- lapply(dxymedlabM, function(dxy) {
  xya <- do.call('rbind',mapply(function(xys,age) data.frame(age=age, x=xys$x,y=xys$y), dxy,ages,SIMPLIFY=FALSE))
  xya$lindex=1
  xya})
#
mapply(function(fit,label) {
  fn <- paste('median_leaf_trajectories_',label,'.csv',sep='')
  write.csv(fit, fn, row.names=FALSE)
  },trajfit, names(trajfit))
# 
#
#
#
#
#
#
# deprecated sphi interpolation fits
#
ages <- seq(0,12,0.5)
dxymed <- mapply(dphi_interpolation,split(sphidb, shapedb$pos), split(shapedb,shapedb$pos), MoreArgs=list(at=ages),SIMPLIFY=FALSE)
dxyfit <- list(trajs=dxymed,ages=ages)
#
# fit per 'genotype'
#
groups <- c('MerciaRht','MerciaRht', 'MerciaRht', 'MerciaRht','Soissons','Tremie', 'Tremie')
names(groups) <- c('Mercia','Rht10','Rht3','Soissons','Tremie12','Tremie13')
g <- groups[shapedb$label]
dxymedg <- mapply(function(sphig, shapeg) mapply(dphi_interpolation,split(sphig, shapeg$pos), split(shapeg, shapeg$pos), MoreArgs=list(at=ages), SIMPLIFY=FALSE),
                  split(sphidb,g), split(shapedb, g),SIMPLIFY=FALSE)
dxyfitg <- lapply(dxymedg, function(dxy) list(trajs=dxy,ages=ages))

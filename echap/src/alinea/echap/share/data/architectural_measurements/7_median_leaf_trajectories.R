#
#            Compute median leaf interpolated trajectories for leaf shape=f(age)
#
#
source('funShapes.R')
#
bl_e <-  read.csv('leafshape_metainfo_Boigneville_2012_2013.csv', sep=',', dec='.')
nervs_e <-  read.csv('leafshape_xydata_Boigneville_2012_2013.csv', sep=',', dec='.')
bl <-  read.csv('leafshape_metainfo_Grignon2010.csv', sep=',', dec='.')
nervs <-  read.csv('leafshape_xydata_Grignon2010.csv', sep=',', dec='.')
# merge indices
imax <- max(bl$inerv)
bl_e$inerv <- bl_e$inerv + imax
nervs_e$inerv <- nervs_e$inerv + imax
#merge
nervs_e <- split(nervs_e,nervs_e$inerv)
nervs <- split(nervs,nervs$inerv)
nervsdb <-c(nervs,nervs_e)
#
cols <- c('inerv','label','Source','plant','rank','ranktop','age','lmax','wmax','A','Agr','mass','stat','phiP')
bldb <- rbind(bl[,cols],bl_e[,cols])
bldb$variety = as.numeric(as.factor(bldb$label))
bldb$harvest = as.numeric(as.factor(bldb$Source))
#
# sphi spline fit
#
sphidb <- lapply(nervsdb, sphi_shapes,df=10)
sxydb <- lapply(sphidb,function(fit) sphi2xy(fit$ds_n,fit$phi))
#
# check fits
#
par(mfrow=c(4,5),mar=c(1,1,1,1))
ids <- sample(seq(nervsdb),20)
mapply(function(dat,fit) {plot(dat$xr,dat$yr);lines(sphi2xy(fit$ds,fit$phi),col=2)}, nervsdb[ids],sphidb[ids])
#
# shapes descriptors
#
shapedb <- get_shapes(sphidb, bldb,topleaves=4)
groups <- c('MerciaRht','MerciaRht', 'MerciaRht', 'MerciaRht','Soissons','Tremie', 'Tremie')
names(groups) <- c('Mercia','Rht10','Rht3','Soissons','Tremie12','Tremie13')
shapedb$group <- groups[shapedb$label]

#
# trajectories
#
ages <- seq(0,12,0.5)
#
# pooled data
#
dxymed <- mapply(sxy_interpolation,split(sxydb, shapedb$pos), split(shapedb,shapedb$pos), MoreArgs=list(at=ages),SIMPLIFY=FALSE)
# fit per 'genotype'
dxymedg <- mapply(function(sxyg, shapeg) mapply(sxy_interpolation,split(sxyg, shapeg$pos), split(shapeg, shapeg$pos), MoreArgs=list(at=ages), SIMPLIFY=FALSE),
                  split(sxydb,shapedb$group), split(shapedb, shapedb$group),SIMPLIFY=FALSE)
#
trajfit <- lapply(dxymedg, function(dxy) {
  xya <- lapply(dxy,function(traj) do.call('rbind',mapply(function(xys,age) data.frame(age=age, x=xys$x,y=xys$y), traj,ages,SIMPLIFY=FALSE)))
  xyai <- do.call('rbind',mapply(function(xy,index) {xy$lindex=index;xy},xya,as.numeric(names(xya)),SIMPLIFY=FALSE))
})
#
# export
#
loc <- c('Grignon2010','Grignon2010','Boigneville2012_2013')
names(loc) <- c('MerciaRht','Soissons','Tremie')
mapply(function(fit,label) {
  fn <- paste('median_leaf_trajectories_',label, '_',loc[label],'.csv',sep='')
  write.csv(fit, fn, row.names=FALSE)
  },trajfit, names(trajfit))  
#
# plot trajectories
#
jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
palette(jet(length(ages)))

#
palette(jet(length(ages)))
par(mfcol=c(2,4), mar=c(4,4,1,1))
leg <- c('lower leaves','upper leaves')
for (i in seq(dxymed)) {
  plot(c(0,1.3),c(-0.3,1),type='n',xlab=leg[i],ylab='interpolated median leaf',cex.lab=1.5)
  sapply(seq(dxymed[[i]]),function(j) {lines(dxymed[[i]][[j]],col=j)})
  if (i==2)
    text(1,0.8,'Pooled')
}
#
for (g in names(dxymedg)) {
  dxy <- dxymedg[[g]]
  for (i in seq(dxy)) {
    plot(c(0,1.3),c(-0.3,1),type='n',xlab=leg[i],ylab='interpolated median leaf',cex.lab=1.5)
    sapply(seq(dxy[[i]]),function(j) {lines(dxy[[i]][[j]],col=j)})
    if (i==2)
    text(1,0.8,g)
  }
}
#
# Check relation with age / median trajectorie projection
#
shapeT <- shapedb[shapedb$group=='Tremie',c('inerv','rank','ranktop','age','pos','hproj','vproj')]
nervT <- do.call('rbind',split(nervsdb,shapedb$group)[['Tremie']])[,c('inerv','s','xr','yr')]
nervp <- merge(nervT,shapeT)
bins <- c(-1,ages)
nervp$ageclass <- cut(nervp$age, bins)
# raw shapes
palette(jet(length(ages)))
plot(nervp$x,nervp$y,col=as.numeric(nervp$ageclass),pch=16,cex=0.5)
#
# hproj
#
symb <- c(1,16)
palette(jet(length(ages)))
plot(nervp$age,nervp$hproj,col=as.numeric(nervp$ageclass),pch=symb[nervp$pos])
palette('default')
lapply(split(nervp,nervp$pos),function(x) lines(smooth.spline(x$age,x$hproj,df=4),col=8))
#
medleaf <- trajfit[['Tremie']]
medp <- do.call('rbind',lapply(split(medleaf,list(medleaf$age,medleaf$lindex),drop=TRUE),function(leaf) data.frame(age=leaf$age[1],lindex=leaf$lindex[1],proj=max(leaf$x))))
lapply(split(medp,medp$lindex),function(x) lines(x$age,x$proj,col=2))


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

#                 Compute shapes of organ profiles using tagged data + MaxwellREF
#
#
# add phyllochronic rank to dimensions of tagged plants
# phyllochronic rank = rank_phy = phylochronic time of mean axis at leaf ligulation = TT_leaf_ligulation / mean_phyllochron
# xn = phyllochronic rank - phyllochronic_rank of longuest leaf
#
profdb <- dimtdb
for (g in names(profdb)) {
  dat <- merge(dimtdb[[g]],TTem[[g]])
  rphy <- dat$TTlig * slopeM[[g]]
  dat$rank_phy <- rphy
  dat <- dat[!is.na(rphy),]
  mat <- na.omit(dat[,c('rank_phy','Lb')])
  spl <- smooth.spline(mat$rank_phy, mat$Lb,df=10)
  ymax <- max(spl$y)
  xmax <- spl$x[which.max(spl$y)]
  dat$xn <- dat$rank_phy - xmax
  profdb[[g]] <- dat
}
#
# add Maxwell data
#
maxdb <- read.csv('../Maxwell_dimdb.csv', sep=',',dec='.')
newnames <- list(L_blade='Lb',W_blade='Wb',L_sheath='Ls', L_internode='Li',W_sheath='Ws',Hcol='Hc')
colnames(maxdb)[match(names(newnames),colnames(maxdb))] <- unlist(newnames)
maxdb$TTlig <- maxdb$TT
maxdb$N <-  maxdb$id_dim %/% 1000
maxdb$Source='Maxwell_REF'
maxdb$Ab <- as.numeric(NA)
maxdb$Li[maxdb$rank < 6] <- 0
maxdb$Ws[maxdb$rank < 6] <- median(na.omit(maxdb$Ws))
dat <- maxdb[,c('N','nff','rank','Lb','Source','Ab','Wb','Ls','Li','Hc','TTlig','Ws','rank_phy')]
dat <- dat[!is.na(dat$rank_phy),]
mat <- na.omit(dat[,c('rank_phy','Lb')])
spl <- smooth.spline(mat$rank_phy, mat$Lb,df=8)
ymax <- max(spl$y)
xmax <- spl$x[which.max(spl$y)]
dat$xn <- dat$rank_phy - xmax
profdb[['Maxwell']] <- dat
#
ym <- list(Lb=30,Wb=3,Ls=30,Ws=1,Li=30,Hc=80)
par(mfrow=c(2,3), mar=c(4,4,1,1))
for (w in names(ym)) {
  plot(c(0,15),c(0,ym[[w]]),xlab='phyllochronic rank',ylab=w,type='n')
  sapply(seq(profdb),function(i) points(profdb[[i]]$rank_phy,profdb[[i]][[w]],col=i,pch=c(16,1,15)[as.factor(profdb[[i]]$nff)]))
}
#
#ym <- list(Lb=1.5,Wb=0.05,Ls=1.5,Ws=0.05,Li=1.5,Hc=4)
par(mfrow=c(2,3), mar=c(4,4,1,1))
for (w in names(ym)) {
  plot(c(-10,5),c(0,ym[[w]]),xlab=' normalised phyllochronic rank',ylab=w,type='n')
  sapply(seq(profdb),function(i) points(profdb[[i]]$xn,profdb[[i]][[w]],col=i,pch=c(16,1,15)[as.factor(profdb[[i]]$nff)]))
}
#
# Compute scales for indivudual plants using maxwell profiles except for blade where Tremie 13 is more sharp
#
fits <- sapply(names(ym), function(w) {
  if (w != 'Lb')
    mat <- na.omit(profdb$Maxwell[,c('xn',w)])
  else
     mat <- na.omit(profdb$Tremie13[,c('xn',w)])
  smooth.spline(mat$xn, mat[[w]],df=8)},simplify=FALSE)
#
scaledb <- lapply(profdb,function(dat) do.call('rbind',lapply(split(dat,dat$N),function(mat) {
  for (w in names(ym)) {
    sc <- NA
    matw <- na.omit(mat[mat[[w]] > 0,c('xn',w)])
    #matw <- matw[matw$xn > 0,]
    if (nrow(matw) >= 1)
      sc <- median(matw[[w]] / predict(fits[[w]],matw$xn)$y)
    mat[[w]] <- sc
   }
  mat})))
#check
par(mfrow=c(2,3), mar=c(4,4,1,1))
for (w in names(ym)) {
  plot(c(-10,5),c(0,ym[[w]]),xlab=' normalised phyllochronic rank',ylab=w,type='n')
  sapply(seq(profdb),function(i) points(profdb[[i]]$xn,profdb[[i]][[w]]/scaledb[[i]][[w]],col=i,pch=16))
  lines(fits[[w]],col=2)
}
#
# Final fits
#
pdb <- do.call('rbind',mapply(function(prof,sc,label) {
  dat <- prof
  for (w in names(ym))
    dat[[w]] <- dat[[w]] / sc[[w]]
  dat$label <- label
  dat},profdb,scaledb,names(profdb),SIMPLIFY=FALSE))
#
fits <- sapply(names(ym), function(w) {
  mat <- na.omit(pdb[,c('xn',w)])
  smooth.spline(mat$xn, mat[[w]],df=8)},simplify=FALSE)
#check
par(mfrow=c(2,3), mar=c(4,4,1,1))
for (w in names(ym)) {
  plot(c(-10,5),c(0,ym[[w]]),xlab=' normalised phyllochronic rank',ylab=w,type='n')
  points(pdb$xn,pdb[[w]],col=as.factor(pdb$label),pch=16)
  lines(fits[[w]],col=2,lwd=2)
  lines(predict(fits[[w]],seq(-15,5,0.1)),col=2)
}
#
# export profiles, normalised by value at xn=0
#
xn <- seq(-13,4.5,0.1)
cols <- list(L_blade='Lb',W_blade='Wb',L_sheath='Ls', L_internode='Li',W_sheath='Ws',Hcol='Hc')
profiles <- data.frame(xn=xn)
for (w in names(cols))
  profiles[[w]] <- predict(fits[[cols[[w]]]], xn)$y
# check negative values
profiles$L_internode[1:max(which(profiles$L_internode < 0))] <- 0
#flaten Ws and add W_internode
profiles$W_sheath <- mean(na.omit(profiles$W_sheath))
profiles$W_internode <- profiles$W_sheath
#
write.csv(profiles, '../Dimension_profiles.csv',row.names=FALSE)

#

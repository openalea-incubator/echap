#
#
#           visual check/correc leaf height data
#
#
#
#
#corrections manuelles
#
#RM00, data mercia rht2 plot 117,3 = data mercia rht3 plot 115,3
#
id <- lhrb$RM00$PLOT*10+lhrb$RM00$OBSNO
p1173 <- id==1173
p1153 <- id==1153
v1173 <- lhrb$RM00[p1173,]
v1153 <- lhrb$RM00[p1153,]
lhrb$RM00[p1153,14:29] <- v1173[,14:29]
lhrb$RM00[p1173,14:29] <- v1153[,14:29]
#
sel1 <- (id == 1151 | id == 1152 | id == 1153) & lhrb$RM00$JULIAN == 123
sel2 <- (id == 1171 | id == 1172 | id == 1173) & lhrb$RM00$JULIAN == 123
lhrb$RM00[sel1,15:29] <- NA
lhrb$RM00[sel2,15:29] <- NA
#
#
fr <- list("RM99" = c(3,4),"HM99"=c(3,5),"HM00" = c(3,6),"RM00"=c(3,6))
#
#view lb
#
kinlb <- function(g,what="BL") {
  dat <- lhd[[g]]
  par(mfrow=fr[[g]],mar=c(1,1,0,0),oma=c(1,1,4,0))
  for (i in levels(dat$SEEDNAME) ) {
    bid <- dat[dat$SEEDNAME == i,]
    plot(bid$JULIAN,bid[,what],col=bid$LN,pch=16,cex=.6,xlim=c(100,180),ylim=c(0,100))
    text(120,80,i,cex=.8)
  }
  mtext(paste(fich[[g]],ifelse(what=="BL","Collar heights","Leaf length")),3,line=2,outer=TRUE)
}
#
kinlb("RM99")
kinlb("RM99","LL")
kinlb("RM00")
kinlb("HM99")
kinlb("HM00")
#
#plante a plante
kinlbp <- function(g,w,what="BL") {
  dat <- lhd[[g]]
  par(mfrow=c(3,3),mar=c(2,2,1,1))
  i <- levels(dat$SEEDNAME)[w]
  bid <- dat[dat$SEEDNAME == i,]
  id <- bid$PLOT*10+bid$OBSNO
  for (r in unique(id)) {
    tmp <- bid[id==r,]
    plot(tmp$JULIAN,tmp[,what],col=tmp$LN,xlim=c(80,180),ylim=c(0,100))
    text(120,80,paste(i,"id",r,"row",paste(unique(tmp$ROW),collapse=",")),cex=.6)
  }
}
#
kinlbp("RM99",12,"LL")

#
#
#essais geometrie
#
#
geom <- function(g,mature=FALSE) {
  dat <- lhd[[g]]
  par(mfrow=fr[[g]],mar=c(1,1,0,0),oma=c(1,1,4,0))
  for (i in levels(dat$SEEDNAME) ) {
    bid <- dat[dat$SEEDNAME == i,]
    plot(c(-50,50),c(0,100), type="n") 
    text(0,0,i,cex=.8)
    id <- bid$PLOT*10+bid$OBSNO
    for (r in unique(id)) {
      tmp <- bid[id==r,]
      for (f in na.omit(unique(tmp$LN)))
        if (length(na.omit(tmp$LN[tmp$LN == f])) > 0) {
        sel <- tmp$LN == f
        y <- tmp$BL[sel]
        if (mature)
          y <- ifelse(y < max(y,na.rm=T)- 5,NA,y)
        dy <- ifelse(is.na(tmp$ML[sel]),tmp$TL[sel]-tmp$BL[sel],tmp$ML[sel]-tmp$BL[sel])
        dy2 <- ifelse(is.na(tmp$ML[sel]),0,tmp$ML[sel]-tmp$BL[sel])
        lseg <- tmp$LL[sel] * ifelse(is.na(tmp$ML[sel]),1,0.5)
        dx <- lseg * cos(asin(dy / lseg))
        dx2 <- ifelse(dy2==0,0,lseg * cos(asin(dy2 / lseg)))
        segments(0,y,ifelse(f%%2 == 0,1,-1)*dx,y+dy,col=f)
        segments(ifelse(f%%2 == 0,1,-1)*dx,y+dy,ifelse(f%%2 == 0,1,-1)*(dx+dx2),y+dy-dy2,col=f)
      }
    }
  }
  mtext(paste(fich[[g]],"Plant Stature"),3,line=2,outer=TRUE)
}
#
geom("RM99",TRUE)
geom("HM99",TRUE)
geom("RM00",TRUE)
geom("HM00",TRUE)
#
geomp <- function(g,w,mature=FALSE) {
  dat <- lhd[[g]]
  par(mfrow=c(3,3),mar=c(2,2,1,1))
  i <- levels(dat$SEEDNAME)[w]
  bid <- dat[dat$SEEDNAME == i,]
  id <- bid$PLOT*10+bid$OBSNO
  for (r in unique(id)) {
    plot(c(-50,50),c(0,100), type="n")
    tmp <- bid[id==r,]
    for (f in na.omit(unique(tmp$LN)))
      if (length(na.omit(tmp$LN[tmp$LN == f])) > 0){
      sel <- tmp$LN == f
      y <- tmp$BL[sel]
      if (mature)
        y <- ifelse(y < max(y,na.rm=T)- 5,NA,y)
      dy <- ifelse(is.na(tmp$ML[sel]),tmp$TL[sel]-tmp$BL[sel],tmp$ML[sel]-tmp$BL[sel])
      dy2 <- ifelse(is.na(tmp$ML[sel]),0,tmp$ML[sel]-tmp$BL[sel])
      lseg <- tmp$LL[sel] * ifelse(is.na(tmp$ML[sel]),1,0.5)
      dx <- lseg * cos(asin(dy / lseg))
      dx2 <- ifelse(dy2==0,0,lseg * cos(asin(dy2 / lseg)))
      segments(0,y,ifelse(f%%2 == 0,1,-1)*dx,y+dy,col=f)
      segments(ifelse(f%%2 == 0,1,-1)*dx,y+dy,ifelse(f%%2 == 0,1,-1)*(dx+dx2),y+dy-dy2,col=f)
    }
  }
}
#
geomp("RM99",11,TRUE)
#
# Leaf stage
#
lstage <- function(g) {
  dat <- lhd[[g]]
  par(mfrow=fr[[g]],mar=c(1,1,0,0),oma=c(1,1,4,0))
  for (i in levels(dat$SEEDNAME) ) {
    bid <- dat[dat$SEEDNAME == i,]
    plot(c(80,180),c(-6,0), type="n") 
    text(150,-5,i,cex=.8)
    id <- bid$PLOT*10+bid$OBSNO
    for (r in seq(unique(id))) {
      tmp <- bid[id==unique(id)[r],]
      llmax <- by(tmp,tmp$LN,function(mat) max(mat$LL,na.rm=T))
      maxl <- tmp$LN
      for (i in seq(maxl))
        if (!is.na(maxl[i])) 
          maxl[i] <- llmax[[as.character(maxl[i])]]
      tmp <- cbind(tmp,dec=tmp$LL/maxl)
      gs <- by(tmp,tmp$JULIAN,function(mat) {
        mat <- mat[order(mat$LN,decreasing=TRUE),]
        -mat$LN[max(which(!is.na(mat$BL)))] + sum(na.omit(mat$dec[seq(mat$LN) > max(which(!is.na(mat$BL)))]))
      })
      points(as.numeric(names(gs)),unlist(gs),col=r,pch=16)
    }
  }
  mtext(paste(fich[[g]],"LeafStage"),3,line=2,outer=TRUE)
}
#
#
lstage("HM99")
lstage("RM00")
#
#press book
#
pdf("ArchiADAS.pdf",paper="a4r", width=0, height=0)
for (g in names(fich)) {
  lstage(g)
  kinlb(g)
  geom(g,TRUE)
}
dev.off()

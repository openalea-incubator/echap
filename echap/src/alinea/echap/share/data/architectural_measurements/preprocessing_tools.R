#
#
#  Preprocessing functions
#
# read scan files
#
readScanned <- function(prefix, scanfiles, name=NULL) {
  scans <- sapply(scanfiles, function(x) {
    df <- read.table(paste(prefix,'_scan_',x,'.txt',sep=''),sep='\t', dec='.',header=TRUE)
    df$Source <- x
    df$N <- df$plant
    if ('pcent_green'%in%colnames(df))
      df$A_bl_green <- df$A_bl * df$pcent_green / 100
    df},simplify=FALSE)
  cols <- c('Source', 'prelevement','N','id_Axe','rank','lmax','wmax','A_bl','A_bl_green','stat', grep('^w[0-9]',colnames(scans[[1]]),value=TRUE))
  do.call('rbind',lapply(scans,function(x) x[,cols[cols%in%colnames(x)]]))
}
#
# read silhouette files
#
readCurv <- function(prefix, scanfiles, name=NULL) {
  curvs <- sapply(scanfiles, function(x) {
    df <- read.table(paste(prefix,'_curvature_',x,'.csv',sep=''),sep=';', dec=',',header=TRUE)
    df$Source <- x
    df$N <- df$ID_Plant
    if ('ID_Metamer' %in% colnames(df))
      df$rank <- df$ID_Metamer
    if ('ID_Metamer_top' %in% colnames(df))
      df$relative_ranktop <- df$ID_Metamer_top
    df},simplify=FALSE)
  do.call('rbind',curvs)
}
#
# Extration leaf length, hcol from silhouettes
#
HcLbCurv <- function(curv) {
  do.call('rbind',lapply(split(curv,list(curv$Source,curv$N),drop=TRUE), function(pl) {
    #print(paste(pl$Source[1],pl$N[1]))
    stem <- na.omit(data.frame(x=unlist(pl[pl$Organ==0 & pl$XY==0,grep('^Pt',colnames(pl))]),y=unlist(pl[pl$Organ==0 & pl$XY==1,grep('^Pt',colnames(pl))])))
    phiP <- mean(atan2(diff(stem$y),diff(stem$x)))
    ranks <- pl[pl$Organ > 0 & pl$XY==0,'rank']
    xbase <- pl[pl$Organ > 0 & pl$XY==0,'Pt1'] - stem$x[1]
    ybase <- pl[pl$Organ > 0 & pl$XY==1,'Pt1'] - stem$y[1]
    #corrective rotation to align with vertical
    theta <- pi / 2 - phiP
    xybase <- list(x = xbase * cos(theta) - sin(theta) * ybase, y = sin(theta) * xbase + cos(theta) * ybase)
    # corrective factor for insertion heigth taking into account curved stems
    lstem <- sum(sqrt(diff(stem$x)^2 + diff(stem$y)^2))
    sc <- lstem / max(stem$y)
    # leaf lengths
    lpl <- pl[pl$Organ > 0,]
    lb <- sapply(split(lpl,lpl$rank), function(mat) {
      xy <- na.omit(t(mat[,grep("Pt",colnames(mat))]))
      sum(sqrt(diff(xy[,1])^2+diff(xy[,2])^2))
    })
    data.frame(Source=pl$Source[1], N=pl$N[1],rank=ranks, Hcol=xybase$y * sc, Lb=lb)
  }))
}
#
# read notations files
#
readNotations <- function(pref, notationfiles, name=NULL) {
  notations <- sapply(notationfiles,function(x) {
    res <- read.table(paste(pref,'_notations_',x,'.txt',sep=''),sep='\t', dec=',',header=TRUE)
    # complete Nflig and  nff
    nflig <- rep(NA,nrow(res))
    nff <- rep(NA,nrow(res))
    if ('nff' %in% colnames(res))
      nff <- res$nff
    if ('Nflig' %in% colnames(res)) {
      nflig <- res$Nflig
      if ('Nfvis' %in% colnames(res) & !('nff' %in% colnames(res)))
        nff <- ifelse(res$Nfvis == 0, res$Nflig, NA)
    }
    res$Nflig <- as.numeric(nflig)
    res$nff <- as.numeric(nff)
    # if sheath_length present + (internode_length measured or rank <=6): compute hins_Fx = Lstem_x + sheath length
    if (length(grep('sheath_length', colnames(res))) > 0) {
      #print(paste(res$Date[1],res$var[1],collapse=' '))
      shcols <- colnames(res)[grep('sheath_length', colnames(res))]
      numphy <- sapply(shcols, function(x) as.numeric(strsplit(x,'_F')[[1]][2]))
      if (length(grep('internode_length', colnames(res))) > 0) {
        encols <- colnames(res)[grep('internode_length', colnames(res))]
        numen <- sapply(encols, function(x) as.numeric(strsplit(x,'_F')[[1]][2]))
        en <- res[,grep('internode_length', colnames(res))]
        sel <- seq(nrow(en))[apply(en,1,function(e) any(!is.na(e)))]
        first <- sapply(sel, function(irow) min(which(!is.na(unlist(en[irow,])))))
        for (i in seq(length(first)))
          en[sel[i],1:max(1,(first[i] - 1))] <- 0
        for (i in numphy) {
          stem <- sapply(seq(nrow(en)), function(irow) sum(en[irow,1:match(i,numen)]))
          res[[paste('Hins_F',i,sep='')]] <- stem + res[[paste('sheath_length_F',i,sep='')]]
        }
      } else {#no internodes
        for (i in numphy[numphy <= 7])
          res[[paste('Hins_F',i,sep='')]] <- res[[paste('sheath_length_F',i,sep='')]]
      }
    }      
    # if Nflig and Hcol present split into Hcol_Fx, and keep Hcol only for Nflig unknonwn
    if (length(na.omit(res$Nflig)) > 0 & 'Hcol' %in% colnames(res)) {
      for (i in sort(unique(na.omit(res$Nflig))))
        res[[paste('Hcol_F',i,sep='')]] <- ifelse(res$Nflig==i,res$Hcol,NA)
      res$Hcol <- ifelse(is.na(res$Nflig), res$Hcol, NA)
    }
    # if number of elongated internode known, set non elongated internode length to zero
    if ('nb_elongated_internode' %in% colnames(res) & length(na.omit(res$nff)) > 0) {
      last_null = res$nff - res$nb_elongated_internode
      for (i in seq(max(last_null)))
        res[[paste('internode_length_F',i,sep='')]] <- ifelse(i <= last_null,0,NA)
    }
    if ('first_elongated_internode' %in% colnames(res)) {
      last_null = res$first_elongated_internode - 1
      for (i in seq(max(last_null)))
        res[[paste('internode_length_F',i,sep='')]] <- ifelse(i <= last_null,0,NA)
    }
    #filter data (Mercai/Rht3)
    if (!is.null(name))
      res<- res[res$var==name,]
    res},simplify=FALSE)
  notations
}
#
# read reference maladie
#
readSymTagged <- function(symfile, TT) {
  symtagged <- read.csv2(symfile)
  symtagged$Date <- sapply(as.character(symtagged$date),function(x) paste(strsplit(x,split='/',fixed=TRUE)[[1]][c(2,1,3)],collapse='/'))
  symtagged <- merge(symtagged,TT)
}
#
# read suivis + TT
#
readTagged <- function(taggedfile,TT, name=NULL) {
# do not consider cut leaves(to be counted later on ?)
  tagged <- read.table(taggedfile,sep='\t', dec=',',header=TRUE, na.strings=c('NA','c'))
  if (!any(grep('fsen',colnames(tagged)))) {
    mat <- tagged[grep('fvert',colnames(tagged))]
    mat[,] <- NA
    colnames(mat) <- sub('fvert','fsen',colnames(mat))
    tagged <- cbind(tagged,mat)
  }
  # comptue fvert/fsen if missing
  tagged <- do.call('rbind',lapply(split(tagged,tagged$Date),function(pl) {
    fvert <- pl[grep('fvert',colnames(pl))]
    fsen <- pl[grep('fsen',colnames(pl))]
    if (all(is.na(fsen)))
      fsen <- 100 - fvert
    if (all(is.na(fvert)))
      fvert <- 100 - fsen
    pl[,grep('fvert',colnames(pl))] <- fvert
    pl[,grep('fsen',colnames(pl))] <- fsen
    pl}))
  tagged <- merge(tagged,TT)
  tagged$numJ <- as.numeric(format(as.Date(as.character(tagged$Date),"%d/%m/%Y"),"%j"))
  #filter data (Mercai/Rht3)
  if (!is.null(name))
    tagged <- tagged[tagged$var==name,]
  tagged
}
#
dimTagged <- function(dat) {
  first <- min(grep('Lg',colnames(dat)))
  last <- max(grep('Lg',colnames(dat)))
  nobs <- apply(dat[,first:last],1,function(x) length(na.omit(x)))
  Flig <- dat[!is.na(dat$Nflig) & nobs > 0,1:last]
  for (i in seq(nrow(Flig))) {
    if (is.na(Flig$nff[i])) {#plant is dead during measurements
      Flig[i,(first + Flig$Nflig[i]):last] <- NA
    } else {
      if (Flig$Nflig[i] < Flig$nff[i])
        Flig[i,(first + Flig$Nflig[i]):last] <- NA
    }
  }
  dimT <- aggregate(Flig[,c(grep('nff',colnames(dat)),first:last)],list(N=Flig$N),mean,na.rm=TRUE)
  dimT <- do.call('rbind',lapply(split(dimT,dimT$N),function(x) {
    x <- data.frame(x)
    n <- length(x) - 2
    mat <- na.omit(cbind(1:n, unlist(x)[3:(2 + n)]))
    data.frame(N=x$N[1],nff=x$nff[1],rank=mat[,1],L=mat[,2])}))
  dimT
}
#
final_length <- function(obs, dimM, dimMM) {  
  nff <- obs$nff[1]
  ref <- dimMM
  if (!is.na(nff))
    ref <- na.omit(dimM[dimM$nff == nff,])
  yobs <- obs$L
  yref <- ref$x[match(obs$rank, ref$rank)]
  scale <- mean(yobs/yref, na.rm=TRUE)
  res <- data.frame(rank=seq(min(ref$rank), max(ref$rank)))
  res$L <- approx(ref$rank, ref$x * scale, res$rank)$y
  res$L[match(obs$rank,res$rank)] <- obs$L
  res$nff <- obs$nff[1]
  res$N <- obs$N[1]
  res[,c('N','nff','rank','L')]
}
#
HaunStage <- function(d,dimf) {
  lb <- d[grep('Lg',colnames(d))]
  hs <- NA
  if (ncol(lb) > 0) {
    if (!all(is.na(lb)) & !is.na(d$nff) & !is.na(d$Nflig)) {
      hs <- d$Nflig
      if (hs < d$nff) {
        rg <- d$Nflig + (1:d$Nfvis)
        frac <- NA
        if  (all(rg%in%dimf$rank))
          frac <- sum(lb[rg]/dimf$L[match(rg,dimf$rank)])  
        hs = d$Nflig + frac
      }
      
    }
  }
  if (!is.na(d$nff) & !is.na(d$Nflig))
    if (d$Nflig >= d$nff)
      hs <- d$Nflig
  hs
}
#
SSI <- function(d) {
  fsen <- d[grep('fsen',colnames(d))]
  ssi <- NA
  if (!all(is.na(fsen)) & !is.na(d$nff) & !is.na(d$Nflig)) {
    last <- max(which(!is.na(fsen))) 
    if (fsen[last] == 0 | last==d$nff | last==d$Nflig)# a growing leaf is considered to be sen=zero if ot noted
      ssi <- sum(fsen[1:last]) / 100
  }
  ssi
}
# range of rank to sum for fractional ssi part
SSImax <- function(d) {
  fsen <- d[grep('fsen',colnames(d))]
  last <- max(which(!is.na(fsen)))
  ssiM <- NA
  if (fsen[last] == 0 | last==d$nff | last==d$Nflig)
    ssiM <- last
  ssiM
}
#
SSImin <- function(d) {
  fsen <- d[grep('fsen',colnames(d))]
  min(which(fsen<100))
}
#
#
pheno <- function(tagged, dimP) {
  phen <- do.call('rbind',lapply(split(tagged,tagged$N),function(pl) {
    dimf <- dimP[dimP$N==pl$N[1],]
    pl$HS <- sapply(seq(nrow(pl)), function(i) HaunStage(pl[i,],dimf))
    pl$SSI <- sapply(seq(nrow(pl)), function(i) SSI(pl[i,]))
    pl$GL <- pl$HS - pl$SSI
    #pl$SSIfirst <- sapply(seq(nrow(pl)), function(i) SSImin(pl[i,]))
    #pl$SSIlast <- sapply(seq(nrow(pl)), function(i) SSImax(pl[i,]))
    pl
  }))
  phen[,-c(grep('Lg',colnames(tagged)), grep('fvert',colnames(tagged)), grep('fsen',colnames(tagged)))]
}
#
hsfit <- function(phen) {
  dat <- phen[phen$HS < phen$nff & !is.na(phen$HS),]
  fit <- NULL
  if (length(unique(dat$TT)) > 1)
    fit <- lsfit(dat$TT,dat$HS)
}
#
fitLspl <- function(spl, a, splref, aref, deform=FALSE) {
  xmax <- spl$x[which.max(spl$y)]
  ymax <- max(spl$y)
  xref <- splref$x[which.max(splref$y)]
  yref <- max(splref$y)
  #
  xfit <- seq(0,1500,10)
  xnorm <- (xfit - xmax) * a / aref + xref
  yfit <- predict(splref,xnorm)$y / yref * ymax
  if (deform) {
    # deformation for first leaves (before max)
    xobs <- spl$data$x[spl$data$x < xmax]
    yobs <- spl$data$y[spl$data$x < xmax]
    xp <- (xobs - xmax) * a / aref + xref
    yp <- predict(splref,xp)$y / yref * ymax
    def <- smooth.spline(xobs,yobs / yp,df=3)
    coef <- predict(def,xfit)$y
    coef[xfit <= min(spl$x)] <- predict(def,min(spl$x))$y
    yfit <- yfit * ifelse(xfit < xmax, coef, 1)
  }
  list(x=xfit, y=yfit)
}
#
# pheno on scanned plant
#
pheno_scan <- function(pl, dim) {
  #print(paste(pl$prelevement[1], pl$N[1]))
  hs <- NA
  if (!is.na(pl$Nfvis[1])) {
    hs <- pl$Nflig[1]
    frac = 0
    if (pl$Nfvis[1] > 0) {
      vis <- pl[pl$rank > pl$Nflig,]
      if (nrow(vis) >= pl$Nfvis[1]) {
        lig <- pl[pl$rank <= pl$Nflig & pl$stat < 3,]
        sc <- 1
        if (nrow(lig) > 0)
          sc <- mean(lig$lmax / dim$L[match(lig$rank,dim$ranks)])
        frac = sum(vis$lmax / sc / dim$L[match(vis$rank,dim$ranks)])
      } else {
        hs <- NA
      }
    }
  }
  ssi <- NA
  p <-pl[!is.na(pl$rank) & !is.na(pl$A_bl_green) & !is.na(pl$A_bl),]
  if (nrow(p) > (max(p$rank) - min(p$rank))) {
    fsen <- 1 - p$A_bl_green / p$A_bl
    if (fsen[which.min(p$rank)] > 0.3)
      ssi <- min(p$rank) - 1 + sum(fsen)
  }
  data.frame(Source=pl$Source[1], Date = pl$prelevement[1], N=pl$N[1], nff=pl$nff[1], Nflig=pl$Nflig[1], Nfvis=pl$Nfvis[1], HS=hs + frac, SSI=ssi, GL=hs+frac-ssi)
}
#
pheno_ssi <- function(ssitable) {
  phen <- do.call('rbind',lapply(split(ssitable,ssitable$N),function(pl) {
    pl$HS <- NA
    ranks <- seq(grep('fsen',colnames(pl)))
    pl$SSI <- sapply(seq(nrow(pl)), function(i) {
      fsen <- pl[i,grep('fsen',colnames(pl))]
      lim <- max(which(fsen >= 100))
      ranks[lim] -1 + sum(fsen[lim:length(fsen)] / 100)
      })
    pl$GL <- pl$HS - pl$SSI
    pl
  }))
  phen[,-c(grep('Lg',colnames(phen)), grep('fvert',colnames(phen)), grep('fsen',colnames(phen)))]
}
#
pheno_symptom <- function(symptom) { 
  phen <- do.call('rbind',lapply(split(symptom, list(symptom$Date, symptom$plant),drop=TRUE),function(pl) {
  nleaf <- length(seq(min(pl$num_leaf_bottom),max(pl$num_leaf_bottom)))
  ssinec <- NA
  ssiap <- NA
  hs <- max(pl$num_leaf_bottom)
  if (hs < pl$fnl[1])
    hs <- NA
  if (nrow(pl) >= nleaf) {
    ssinec <- min(pl$num_leaf_bottom) - 1 + sum(pl$necro_tot) / 100
    ssiap <- min(pl$num_leaf_bottom) - 1 + sum(pl$necro) / 100
  }
  data.frame(Date=pl$Date[1], var='Tremie', Rep=pl$bloc[1], N=pl$plant[1], Axis=pl$axis[1], Nflig=pl$fnl[1], Nfvis=0, nff=pl$fnl[1], TT=pl$TT[1], HS=hs, SSI=ssinec, SSIap=ssiap, GL=hs - ssinec, GLap=hs - ssiap)
}))
}
#
dim_notations <- function(not, what='sheath_length')  {
  sources <- names(not)
  sources <- sources[sapply(sources, function(x) length(grep(what,colnames(not[[x]]))) > 0)]
  res <- NULL  
  if (length(sources) > 0) {
    res <- do.call('rbind', lapply(sources, function(s) {
      #print(s)
      nt <- not[[s]]
      cols <- grep(what,colnames(nt))
      sel <- seq(nrow(nt))
      if ('axe' %in% colnames(nt))
        sel <- nt$axe=='MB'
      lg <- nt[sel,cols]
      lg$N=nt$N[sel]
      numphy <- sapply(colnames(nt)[cols], function(x) as.numeric(strsplit(x,'_F')[[1]][2]))
      dat <- do.call('rbind', lapply(split(lg, lg$N), function(x) data.frame(Source=s, N=x$N,rank=numphy, L=unlist(x[1,seq(length(cols))]))))
      dat[!is.na(dat$L),]
    }))
  }
  res
}
#
plant_notations <- function(not)  {
  what <- c('nff','Nflig','Daxe_mm','Hcol','dh_ped','nb_elongated_internode','lped','Wped_mm','H_node','first_elongated_internode')
  columns <- c('Source','N',what)
  sources <- names(not)
  sources <- sources[sapply(sources, function(x) any(what%in%colnames(not[[x]])))]
  res <- NULL
  if (length(sources) > 0) {
    res <- do.call('rbind', lapply(sources, function(s) {
      nt <- not[[s]]
      cols <- what[what%in%colnames(nt)]
      sel <- seq(nrow(nt))
      if ('axe' %in% colnames(nt))
        sel <- nt$axe=='MB'
      
      dat <- nt[sel,cols]
      if (is.null(dim(dat))) {
        data <- dat
        dat <- data.frame(Source=rep(s,length(dat)))
        dat[[cols]] <- data
      }
      
      dat$Source <- s
      dat$N <- nt$N[sel]
      dat$nff <- nt$nff[sel]
      dat$Nflig <- nt$Nflig[sel]
      for (w in what[!what %in% cols])
        dat[[w]] <- NA
      dat <- dat[,match(columns,colnames(dat))]
      dat
    }))
  }
  res
}
#
# view dimension from notation
#
view_dim <- function(data, what='L',ylim=c(0,30)) {
  par(mfrow=c(2,2),mar=c(4,4,1,1))
  lapply(data, function(dim) {
    plot(c(0,15),ylim,type='n')
    if (!is.null(dim))
      lapply(split(dim,list(dim$Source,dim$N)), function(x) points(x$rank,x[[what]],col=x$N,pch=16,type='b'))
  })
}
#
# add dimension data from a notation source for tagged plants
#
add_dimt <- function(db, notsource, as='Ls'){
  if (is.null(notsource)) {
    db[[as]] <- as.numeric(NA)
  } else {
    dat <- notsource[grep('tagged', notsource$Source),]
    if (nrow(dat) <=0) {
      db[[as]] <- as.numeric(NA)
    } else {
      dat[[as]] <- dat$L
      #add nff
      nff <- db[,c('N','nff')]
      dat <- merge(dat, nff[!duplicated(nff),])
      dat <- dat[,c('N','nff','rank',as)]
      db <- merge(db,dat,all=TRUE)
    }
  }
  db
}
#
# add dimension data from a notation source for sampled non-tagged plants
#
add_dims <- function(db, newsource, as='Ls') {
  for (g in names(newsource)) {
    if (is.null(newsource[[g]])) {
      if (!is.null(db[[g]]))
        if (!as%in%names(db[[g]]))
          db[[g]][[as]] <- NA
    } else {
      new <- newsource[[g]]
      new[[as]] <- new$L
      dat <- merge(db[[g]], new[,c('Source','N','rank',as)],all=TRUE)
      db[[g]] <- dat
    }
  }
  db
}
#
# deprecated ?
rssi_patternT <- function(n,nf,ssisenT=data.frame(ndel=1:4,rate1=0.07, dssit1=c(0,1.2,2.5,3),dssit2=c(1.2,2.5,3.7,4)),hasEar=TRUE) {
  ndelsen <- max(ssisenT$ndel)
  pattern <- list(t=c(-1, 0),p=c(0,1))
  if (hasEar & n > (nf - ndelsen)) {
    idel <- n - (nf - ndelsen)
    t0 <- -idel
    t1 <- t0 + ssisenT$dssit1[idel]
    t2 <-  min(t0 + ssisenT$dssit2[idel],nf - n)
    if (nf < ndelsen) {
      t0 <- -nf
      t1 <- min(nf - n, max(t1,t0))#nf - n is complete senescence of last leaf
      t2 <- min(nf - n, max(t2,t1))
      }
    p1 <- ssisenT$rate[idel] * (t1 - t0)
    pattern <- list(t=c(t0,t1,t2),p=c(0,p1,1))
  }
  pattern
}
#
ssi_table <- function(r1=.1, ndel=3) {
  table <- matrix(0,ncol=ndel,nrow=ndel)
  table[1,] <- c(rep(r1,ndel-1),1-(ndel-1)*r1)
  for (i in 2:ndel) {
    if ((ndel-i) >= 1)
      table[i,1:(ndel-i)] <- r1
    table[i,ndel-i+2] <- 1 - sum(table[1:(i-1),(ndel-i+2)])
    table[i,ndel-i+1] <- 1 - sum(table[i,])
  }
  table
}
#
rssi_pattern <- function(n,nf,hasEar=TRUE,pars=list(r1=0.07,ndelsen=3)) {
  pattern <- list(t=c(-1, 0),p=c(0,1))
  ndel <- min(pars$ndelsen,nf)
  if (ndel > 1 & hasEar & (nf - n) < pars$ndel) {
    table <- ssi_table(r1=pars$r1,ndel=ndel)
    t <- ((nf - ndel):nf) - n
    p <- cumsum(c(0,table[nf - n + 1,]))
    pattern <- list(t=t,p=p)
  }
  pattern
}
#
#
#
#
# fit HS
#
#
# HS model (Jessica/Mariem):
# - hypothese convergence moyenne phase 2 + synchro moyenne changement de phase => HS/nff definit le phyllo en phase 2, Variabilite inter-plante = yorr
#
fiths <- function(phen) {
  dat <- phen[phen$HS < phen$nff & !is.na(phen$HS),]
  phyl <- sapply(split(dat,dat$N), function(x) ifelse(nrow(x) >=2, lsfit(x$TT,x$HS / x$nff)$coeff[2],NA))
  phylM <- mean(na.omit(phyl))
  do.call('rbind',lapply(split(dat,dat$N),function(x) data.frame(N=x$N[1], var=x$var[1],Trt=x$Trt[1],Rep=x$Rep[1], Axis=x$Axis[1],nff=x$nff[1], phyl = 1 / phylM / x$nff[1],a=phylM*x$nff[1], b=mean(x$HS - phylM*x$nff*x$TT))))
}
#
fitTT0 <- function(phen, slope) {
   dat <- phen[phen$HS < phen$nff & !is.na(phen$HS),]
   b <- mean(dat$HS - slope * dat$TT)
   -b/slope
}
#
TTleaf <- function(dat, slope) {
  res <- NULL
  TT0 <- fitTT0(dat, slope)
  nff <- dat$nff[1]
  if (!is.na(nff)) {
    rank <- seq(nff)
    res <- data.frame(N=dat$N[1],nff=nff,rank=rank,TTlig=TT0 + rank / slope)
  }
  res
}
#
TTphy <- function(fitp) {
  n <- seq(fitp$nff)
  hs <- sapply(n, function(x) (x-fitp$b)/fitp$a)
  lig <- hs + 0.3 * fitp$phyl
  tip <- hs - 1.3 * fitp$phyl
  data.frame(nff=fitp$nff[1],rank=n,ranktop = fitp$nff[1] - n + 1,HS=hs,TTem=tip,TTlig=lig)
}

                                       #
# flatten tagged files (for guillaume)
#
flattaged <- function(tagged) {
  do.call('rbind',lapply(split(tagged,seq(nrow(tagged))),function(pl) {
    res <- NULL
    #print(pl[,1:5])
    lg <- pl[grep('Lg_',colnames(pl))]
    if (length(lg) > 0) {
      lg <- na.omit(cbind(seq(ncol(lg)),unlist(lg)))
      if (nrow(lg) > 0) {
        fv <-  pl[grep('fvert_',colnames(pl))]
        lastv <- max(which(!is.na(fv)))
        if (fv[lastv] >= 100)
          fv[lastv:length(fv)] <- 100
        firstv <-  min(which(!is.na(fv)))
        if (fv[firstv] <= 0)
          fv[1:firstv] <- 0
        fv <- cbind(seq(ncol(fv)),unlist(fv))
        res <- data.frame(Date=pl$Date[1],var=pl$var[1],Trt=pl$Trt[1],Rep=pl$Rep[1], N=pl$N[1], Axis=pl$Axis[1],Nflig=pl$Nflig[1],Nfvis=pl$Nfvis[1],nff=pl$nff[1],num_leaf_bottom = lg[,1],length=lg[,2],green_length=lg[,2]*fv[fv[,1]%in%lg[,1],2]/100, senesced_length=lg[,2]*(1 -fv[fv[,1]%in%lg[,1],2]/100))
      }
    }
    res
  }))
}

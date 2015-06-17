#
#
#  Preprocessing functions
#
# read scan files
#
readScanned <- function(prefix, scanfiles, name=NULL) {
  scans <- sapply(scanfiles, function(x) {
    df <- read.table(paste(prefix,'_',x,'.txt',sep=''),sep='\t', dec='.',header=TRUE)
    df$Source <- x
    df$N <- df$plant
    df},simplify=FALSE)
  cols <- c('Source', 'prelevement','N','id_Axe','rank','lmax','wmax','A_bl','A_bl_green','stat','pcent_green', 'Nflig', grep('^w[0-9]',colnames(scans[[1]]),value=TRUE))
  do.call('rbind',lapply(scans,function(x) x[,cols[cols%in%colnames(x)]]))
}
#
# read silhouette files
#
readCurv <- function(prefix, scanfiles, name=NULL) {
  curvs <- sapply(scanfiles, function(x) {
    df <- read.table(paste(prefix,'_',x,'.csv',sep=''),sep=';', dec=',',header=TRUE)
    df$Source <- x
    df$N <- df$ID_Plant
    df$rank <- df$ID_Metamer
    df},simplify=FALSE)
  do.call('rbind',curvs)
}
#
# read notations files
#
readNotations <- function(prefix, notationfiles, name=NULL) {
  notations <- sapply(notationfiles,function(x) {
    res <- read.table(paste(prefix,'_notations_',x,'.txt',sep=''),sep='\t', dec=',',header=TRUE)
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
  #print(paste(pl$prelevement[1], pl$plant[1]))
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
  data.frame(Source=pl$Source[1], Date = pl$prelevement[1], N=pl$plant[1], nff=pl$nff[1], Nflig=pl$Nflig[1], Nfvis=pl$Nfvis[1], HS=hs + frac, SSI=ssi, GL=hs+frac-ssi)
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
dim_notations <- function(not, what='sheath_length', as='sheath')  {
  sources <- names(not)
  sources <- sources[sapply(sources, function(x) length(grep(what,colnames(not[[x]]))) > 0)]
  res <- NULL  
  if (length(sources) > 0) {
    res <- do.call('rbind', lapply(sources, function(s) {
      nt <- not[[s]]
      cols <- grep(what,colnames(nt))
      sel <- seq(nrow(nt))
      if ('axe' %in% colnames(nt))
        sel <- nt$axe=='MB'
      lg <- nt[sel,cols]
      lg$N=nt$N[sel]
      lg$nff = NA
      lg$Nflig = NA
      if ('Nflig' %in% colnames(nt)) {
        lg$Nflig = nt$Nflig[sel]
        if ('Nfvis' %in% colnames(nt)) {
          lg$nff = ifelse(nt$Nfvis[sel] == 0, nt$Nflig[sel], NA)
        } else {
          lg$nff = lg$Nflig#if Nfvis absent, nff=Nflig
          }
      }
      numphy <- sapply(colnames(nt)[cols], function(x) as.numeric(strsplit(x,'_F')[[1]][2]))
      dat <- do.call('rbind', lapply(split(lg, lg$N), function(x) data.frame(Source=s, N=x$N, nff=x$nff, organ=as,rank=numphy, L=unlist(x[1,seq(length(cols))]))))
      dat[!is.na(dat$L),]
    }))
  }
  res
}
#
plant_notations <- function(not)  {
  what <- c('Daxe_mm','Hcol','dh_ped','nb_elongated_internode','lped','Wped_mm','H_node','H_last_col')
  columns <- c('Source','N','nff','Nflig',what)
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
      dat$nff <- NA
      dat$Nflig <- NA
      if ('Nflig' %in% colnames(nt)) {
        dat$Nflig = nt$Nflig[sel]
        if ('Nfvis' %in% colnames(nt)) {
          dat$nff = ifelse(nt$Nfvis[sel] == 0, nt$Nflig[sel], NA)
        } else {
          dat$nff = lg$Nflig#if Nfvis absent, nff=Nflig
        }
      }
      for (w in what[!what %in% cols])
        dat[[w]] <- NA
      dat <- dat[,match(columns,colnames(dat))]
      dat
    }))
  }
  res
}
#
#
#
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

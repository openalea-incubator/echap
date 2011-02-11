#
#              Lecture fichier ADAS pour reconstruction ADEL
#              Nov 2009
#
#
#fonctions generiques importee
#
source("FunCalage.R")
#
fich <- list(RM99="RM 1999",RM00="RM 2000",HM99="HMo 1999",HM00="HMo 2000")
#
narestore <- function(x, nacode=-1) {
    res=x
    if (class(x)[1]=="numeric") {
        res=as.numeric(ifelse(x == -1,NA,x))
    }
    res
}
#
#
#Leaf size and GS
#
readGSls <- function(what,dir="../XLData/",end=" leaf size and crop GS.xls",ong=NULL) {
  print(what)
  file <- paste(dir,what,end,sep="")
  dat <- readExcell(file, "GrowthStages")
  assd <- grep(glob2rx("ASS*JD"),colnames(dat))
  for (i in min(assd):ncol(dat))
    dat[,i] <- as.numeric(narestore(dat[,i]))
  assv <- sapply(assd,function(c) length(na.omit(dat[,c])) > 0)
  assd <- assd[assv]
  #
  if (is.null(ong)) {
    ch <- odbcConnectExcel(file)
    ong <- sqlTables(ch)$TABLE_NAME
    odbcClose(ch)
  }
  ongi <- grep(glob2rx("ASS*"),ong)
  dim <- readExcell(file, ong[1],asis=TRUE)
  for (i in ongi[-1]) {
    print(paste("read onglet",ong[i]))
    dim <- rbind(dim,readExcell(file, ong[i],asis=TRUE)[,1:ncol(dim)])
  }
  ld <- grep(glob2rx("L*"),colnames(dim))
  for (i in ld)
    dim[,i] <- as.numeric(narestore(dim[,i]))
  dim <- dim[!is.na(dim[,1]),]
  gshead <- cbind(dat[,1:(min(assd)-1)])
  gs <- NULL
  for (i in seq(assd))
    gs <- rbind(gs,cbind(gshead,ASS = i, JD = dat[,assd[i]], GSMIN = dat[,assd[i]+1], GSMAX = dat[,assd[i]+2]))
  list(gs=gs,dim=dim)
}
#
gsdim <- lapply(fich[-4],readGSls)
gsdim$HM00 <- readGSls(fich$HM00,ong=paste("ASS",1:13,"$",sep=""))
#
gsdb <- lapply(gsdim,function(l) l$gs)
dimdb <- lapply(gsdim,function(l) l$dim)
#
#leaf heigth
#
# hauteurs absolues : onglet 8 : data missing, onglet 15 : autre site !
readlha <- function(what,dir="../XLData/",end=" leaf height.xls",ong=paste("ASS",(1:17)[c(-8,-15)],"$",sep="")) {
  print(what)
  file <- paste(dir,what,end,sep="")
  lh <- readExcell(file, ong[1],asis=TRUE)
  ld <- grep(glob2rx("L*"),colnames(lh))
  ihead <- 1:(min(ld)-1)
  lhd <- NULL
  #
  for (i in seq(ong)) {
    print(ong[i])
    lh <- readExcell(file, ong[i],asis=TRUE)
    ld <- grep(glob2rx("L*"),colnames(lh))
    lh <- lh[!is.na(lh[,1]),1:max(ld)]
    for (n in 1:6) {
      w <- grep(as.character(n),colnames(lh))
      if (length(w) > 0)
        lhd <- rbind(lhd,cbind(lh[,ihead],LN = n, BL = narestore(lh[,w[1]]), ML = narestore(lh[,w[2]]), TL = narestore(lh[,w[3]]), LL = narestore(lh[,w[4]])))
    }
  }
  lhd
}
#
lhda <- readlha("HMo 1999")
#
# hauteurs relatives
#
readlhr <- function(what,ong,dir="../XLData/",end=" leaf height.xls") {
  print(what)
  lhd <- NULL
  file <-  paste(dir,what,end,sep="")
  for (i in seq(ong)) {
    print(ong[i])
    lh <- readExcell(file, ong[i],asis=TRUE)
    if (length(grep("HMo",what)) > 0)
      ld <- grep("[MTBL]L[-+]",colnames(lh))
    else
      ld <- grep("[MTBL]L$",colnames(lh))
    lh <- lh[!is.na(lh[,1]),1:max(ld)]
    for (c in ld)
      lh[,c] <- narestore(lh[,c])
    lhd <- rbind(lhd,lh)
  }
  #homogeneisation seedname RM00
  if (what == fich$RM00) {
    l <- levels(lhd$SEEDNAME)
    l[l=="MERCIA "] <- "MERCIA"
    l[l=="MERCIA PPD1"] <- "MERCIA Ppd1"
    l[l=="MERCIA PPD2"] <- "MERCIA Ppd2"
    l[l=="MERCIA RHT 2"] <- "MERCIA Rht2"
    l[l=="MERCIA RHT2"] <- "MERCIA Rht2"
    l[l=="MERCIA RHT3"] <- "MERCIA Rht3"
    l[l=="MERCIA Rht 3"] <- "MERCIA Rht3"
    levels(lhd$SEEDNAME) <- l
  }
  lhd
}
#
onglhr <- list(HM00=paste("ASS",1:14,"$",sep=""),RM99=paste("ASS",1:23,"$",sep=""),RM00=paste("ASS",1:20,"$",sep=""))
lhrb <- onglhr
for (g in names(lhr))
  lhrb[[g]] <- readlhr(fich[[g]],onglhr[[g]])
#
#Transformation tag -> numero depuis le haut (ramene RM99 et RM00 au cas HM00)
#
tagseq <- list(o=1,g=2,r=3,b=4,y=5,o2=6)
#
col2num <- function(dat) {
  numb <- sapply(dat$TAGCOL,function(c) if (is.na(c))NA else tagseq[[c]])
  numt <- rep(NA,length(numb))
  id <- dat$PLOT*10+dat$OBSNO
  for (i in unique(id)) {
    sel <- (id == i)
    #correction tagage regressif
    tag <- numb[sel]
    tag[which(diff(tag) < 0)] <- tag[which(diff(tag) < 0)] - 1
    #conversion
    if (length(na.omit(tag)) > 0) 
      numt[sel] <- max(tag,na.rm=T) + 1 - tag
  }
  dat$TAGCOL <- numt
  dat
}
#
# Correction erreur rajout tag flag leaf par Julie : lorsque tag = 1, c'est en fait (sauf exception) tag = 2 (conformement au protocole). Avant cette date, T+1BL est une copie de TBL, a effacer pour avoir une trace de l'apparition flag leaf
#attention aux feuilles numerotee 0 !
#
corflag <- function(dat) {
  #homogeneisation noms de colones
  colnames(dat)[15:29] <- c(paste(rep(c("T-2","T-1"),c(3,3)),rep(c("BL","ML","TL"),2),sep=""),
                            paste(rep(c("T","T+1"),c(4,4)),rep(c("BL","ML","TL","LL"),2),sep=""),
                            "T+2TL")
  #detection/correction (RM99) de lignes ou tag flag = 1 en vrai (en general le jour le l'apparition)
  trueflag <- is.na(dat$"T+1BL") & dat$TAGCOL == 1
  #effacement T+1BL pour autre que flag, si TBL=T+1BL (julie a recopie tbl dans t+1bl ?)
  dat$"T+1BL"[dat$TBL == dat$"T+1BL" & dat$TAGCOL != 1] <- NA
  # correction tag =1 est en fait tag = 2 (sauf trueflag)
  dat$TAGCOL <- ifelse(dat$TAGCOL == 1 & !trueflag, 2, dat$TAGCOL)
  dat
}
#
#mise en forme pour lhd: to do garder une trace du numero de la derniere ligulee (avec tagcol de lhrb + ask jilian pour ldha)
#
dec <- c(-2,-1,0,1,2)
blc <- c(15,18,21,25,29)
#
meflhr <- function(dat) {
  lhd <- NULL
  numt <- dat$TAGCOL
  for (l in seq(dec)) {
    if (dec[l] < 0)
      lhd <- rbind(lhd,cbind(dat[,1:14],LN=numt-dec[l],BL=dat[,blc[l]],ML=dat[,blc[l]+1],TL=dat[,blc[l]+2],LL=NA))
    else if (dec[l] < 2)
      lhd <- rbind(lhd,cbind(dat[,1:14],LN=numt-dec[l],BL=dat[,blc[l]],ML=dat[,blc[l]+1],TL=dat[,blc[l]+2],LL=dat[,blc[l]+3]))
    else
      lhd <- rbind(lhd,cbind(dat[,1:14],LN=numt-dec[l],BL=NA,ML=NA,TL=dat[,blc[l]],LL=NA))
  }
  lhd
}
#
lhr <- lhrb
for (g in c("RM00","RM99"))
  lhr[[g]] <- col2num(lhrb[[g]])
lhr <- lapply(lhr,corflag)
lhd <- c(lapply(lhr,meflhr),list(HM99=lhda))

#
#                 disease data
#

# funtion to read HMo 1999 disease data.xls and rearrange the data
readDisease <- function(what,dir="../XLData/",end=" disease data.xls",ong=NULL) {
    print(what)    
    file <- paste(dir,what,end,sep="") # concate input parameters
    require(RODBC)
    if (is.null(ong)) {
        # detect table names automatically
        ch <- odbcConnectExcel(file)
        ong <- sqlTables(ch)$TABLE_NAME
        odbcClose(ch)
    }
    # select the interesting tables
    ongi <- grep("^(ASS|ass)[0-9]{1,2}.$",ong)
    # read the first table to fixe which columns have to be filled
    disease <- readExcell(file, ong[1],asis=TRUE)
    # read all the other tables
    for (i in ongi) {
        print(paste("read onglet",ong[i]))
        # put each table at the end of the precedent one ; the coherence between column names and colum numbers is also checked 
        disease <- rbind(disease,readExcell(file, ong[i],asis=TRUE)[,1:ncol(disease)])
    }
    # replace "-1" by "NA"
    disease = data.frame(lapply(disease,narestore))
    # keep only the set data (i.e. those which are not nan)
    disease[!is.na(disease[,1]),]
   
}

# launch readDisease for all files and store the results list in diseasedb
diseasedb <- lapply(fich, readDisease)


#
#                 meteo data
#


# funtion to read HMo met data.xls and rearrange the data
readHMoMeteo <- function(what="HMo",dir="../XLData/",end=" met data.xls",ong=NULL) {
    print(what)    
    file <- paste(dir,what,end,sep="")
    require(RODBC)
    # detect table names automatically
    if (is.null(ong)) {
        ch <- odbcConnectExcel(file)
        ong <- sqlTables(ch)$TABLE_NAME
        odbcClose(ch)
    }
    print(paste("read onglet",ong))
    meteo <- readExcell(file, ong, asis=TRUE)
    # find the indexes of the "Date" columns
    dateColsIndex <- grep("Date", colnames(meteo))
    # find the indexes of the "Wind" columns
    windColsIndex <- grep("Wind", colnames(meteo))
    res <- NULL
    # for each element in dateColsIndex...
    for (i in seq(dateColsIndex)) {
      # get the value for each column between dateColsIndex[i] and windColsIndex[i] and for each line except the 2 first one)
      newdat <- meteo[-1:-2,dateColsIndex[i]:windColsIndex[i]]
      if (i > 1)
        # keep the column names read at the first loop pass
        colnames(newdat) <- colnames(res)
      # put each table at the end of the precedent one ; the coherence between column names and colum numbers is also checked 
      res <- rbind(res,newdat)
    }
    # replace "-9999" by "NA"
    res = data.frame(lapply(res,narestore,nacode=-9999))
    # replace "6999" by "NA"
    res = data.frame(lapply(res,narestore,nacode=6999))
    # keep only the set data (i.e. those which are not nan)
    res[!is.na(res[,1]),]
   
}

metHMj <- readHMoMeteo()

#
#               RM met data.xls
#

# funtion to read "New logger hourly" table from "RM met data.xls" and rearrange the data
readRMhMeteo <- function(what="RM",dir="../XLData/",end=" met data.xls",ong=NULL) {
    print(what)
    file <- paste(dir,what,end,sep="")
    require(RODBC)
    # detect table names automatically
    if (is.null(ong)) {
        ch <- odbcConnectExcel(file)
        ong <- sqlTables(ch)$TABLE_NAME
        odbcClose(ch)
    }
    # select the interesting table
    ongi <- grep("^New logger hourly$",ong)
    print(paste("read onglet",ongi))
    meteo <- readExcell(file, ongi, asis=TRUE)
    # find the indexes of the "Date" columns
    dateColsIndex <- grep("Date", colnames(meteo))
    # find the indexes of the "Wind" columns
    windColsIndex <- grep("Wind", colnames(meteo))
    res <- NULL
    # for each element in dateColsIndex...
    for (i in seq(dateColsIndex)) {
      # get the value for each column between dateColsIndex[i] and windColsIndex[i] and for each line except the 2 first one)
      newdat <- meteo[-1:-2,dateColsIndex[i]:windColsIndex[i]]
      if (i > 1)
        # keep the column names read at the first loop pass
        colnames(newdat) <- colnames(res)
      # put each table at the end of the precedent one ; the coherence between column names and colum numbers is also checked 
      res <- rbind(res,newdat)
    }
    # replace "-9999" by "NA"
    res = data.frame(lapply(res,narestore,nacode=-9999))
    # replace "6999" by "NA"
    res = data.frame(lapply(res,narestore,nacode=6999))
    # keep only the set data (i.e. those which are not nan)
    res[!is.na(res[,1]),]
   
}
# new logger hourly table
metRMh <- readRMhMeteo()

# funtion to read "Old logger" table from "RM met data.xls" and rearrange the data
metRMj


# merge "RM met data" and "HMo met data"
meteodb <- list(HMj=metHMj, RMh=metRMh, RMj=metRMj)

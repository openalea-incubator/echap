#
#
#           Compilations donees dimensions
#
#
# Lblade plantes tagged
#
Lb_tagged <- lapply(tagged, dimTagged)
#
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(Lb_tagged, function(dim) {
  plot(c(0,15),c(0,30),type='n')
  lapply(split(dim,dim$N), function(x) points(x$rank,x$L,col=x$N,pch=16,type='b'))
})
#
# dimensions per rank
#
Lsheath <- lapply(notdb, function(x) dim_notations(x,'sheath_length', 'sheath'))
Linternode <- lapply(notdb, function(x) dim_notations(x,'internode_length', 'internode'))
Hcol <- lapply(notdb, function(x) dim_notations(x,'Hcol_', 'col'))
#
# TO DO : filter growing organ
#
#
#view sheath
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(Lsheath, function(dim) {
  plot(c(0,15),c(0,30),type='n')
  if (!is.null(dim))
    lapply(split(dim,dim$N), function(x) points(x$rank,x$L,col=x$N,pch=16,type='b'))
})
#
#view internodes
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(Linternode, function(dim) {
  plot(c(0,15),c(0,30),type='n')
  if (!is.null(dim))
    lapply(split(dim,dim$N), function(x) points(x$rank,x$L,col=x$N,pch=16,type='b'))
})
#
#view hcol
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(Hcol, function(dim) {
  plot(c(0,15),c(0,30),type='n')
  if (!is.null(dim))
    lapply(split(dim,dim$N), function(x) points(x$rank,x$L,col=x$N,pch=16,type='b'))
})

#
# plant level variables
plant_dim <- lapply(notdb, function(x) plant_notations(x)) 

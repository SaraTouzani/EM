nk <-
function(tkx){
  k<-ncol(tkx)
  nk<-list(0, k)
  for (i in 1:k) {
    nk[i]<-sum(tkx[,i])
  }
  return(unlist(nk))
}

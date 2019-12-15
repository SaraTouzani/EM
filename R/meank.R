meank <-
function(Nk,tkx,Xquanti,Y){
  K<-length(Nk)
  p<-ncol(Xquanti)
  moyenne<-matrix(0, ncol = p, nrow =K )
  moyenne[,1]<-unique(Y)
  for (k in 1:K) {
    for (j in 1:p) {
      moyenne[k,j]<-sum((tkx[,k]*Xquanti)[,j])/Nk[k]
    }
  }
  moyenne<-cbind(sort(unique(Y)),moyenne)
  return(moyenne)
}

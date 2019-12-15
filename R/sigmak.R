sigmak <-
function(Nk,tkx,Xquanti,Meank){
  K<-length(Nk)
  n<-nrow(Xquanti)
  p<-ncol(Xquanti)
  sigma<-as.list(seq_len(K))
  for (k in 1:K) {
    sig<-matrix(0,ncol = p , nrow = p)
    for (i in 1:n) {
      sig<-sig + ((tkx[i,k]) * t(Xquanti[i,]- Meank[k,-1])%*%as.matrix(Xquanti[i,]-Meank[k,-1])) 
    }
    sigma[[k]]<-sig/Nk[k]
  }
  return(sigma)
  
}

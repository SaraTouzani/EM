Bic <-
function(f,Y=NULL,nbr_par=NULL,Xquali=NULL,Xquanti=NULL,K=length(unique(Y))) {
  nbr_par<-0
  if(!is.null(Xquali)){
    p<-ncol(Xquali)
    
    for (j in 1:p) {
      nbr_par<-nbr_par+length(unique(Xquali[,j])) -1
    }
    nbr_par<-nbr_par*K
  }
  if(!is.null(Xquanti)){
    p<-ncol(Xquanti)
    nbr_par<- K*(p*(p+1)/2+p+1)-1 + nbr_par    
  }

  return(-2*sum(log(f)) + nbr_par*log(length(Y)))
}

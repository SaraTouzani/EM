xijh <-
function(Xquali){
  if(is.null(dim(Xquali))){
    p<-1
    Xijh<-vector("list",1)
    Xquali<-as.matrix(Xquali)
  }
  else{
    p<-length(colnames(Xquali))
    Xijh<-as.list(colnames(Xquali))
  }
  for (j in 1:p) {
    
    Xih<-matrix(0, ncol = length (levels(Xquali[,j])), nrow = nrow(Xquali))
    modalite_ij<-matrix("",  nrow = nrow(Xquali))
    h_index<-1
    for (h in sort(levels(Xquali[,j]))) {
      for (i in 1:nrow(Xquali)) {
        if((Xquali[i,j]==h)){
          Xih[i,h_index]<-1
          modalite_ij[i]<-h
        }
        else{
          Xih[i,h_index]<-0
        }
        
      }
      h_index<-h_index+1
      
    }
    Xijh[[j]]<-data.frame(proba=Xih,
                          modalite=modalite_ij)
    
  }
  
  return(Xijh)
}

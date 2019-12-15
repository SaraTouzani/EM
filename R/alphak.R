alphak <-
function(tkx,Xijh,Y,Xquali){
  
  alpha_k<-as.list(unique(Y))
  
  for (k in 1:length(alpha_k)) {
    variables<-as.list(seq_len(length(Xijh)))
    for (j in 1:length(variables)) {
      
      X<-vector("list",length(levels(Xquali[,j])))
      h_index<-1
      for (h in sort(levels(Xquali[,j]))) {
        s<-0
        for (i in 1:nrow(Xquali)) {
          s<-s+ tkx[i,k]*Xijh[[j]][i  ,-ncol(Xijh[[j]])][h_index]
        }
        X[[h_index]]<-s
        h_index<-h_index+1
      }
      variables[[j]]<-data.frame(proba=unlist(X)/sum(unlist(X))
                                 #,Var=variables[[j]],
                                # modalite=unlist(attributes(table( Xquali_k[[k]][,j] )/ nrow(Xquali_k[[k]]))$dimnames)
                                )
      
    }
    alpha_k[[k]]<-variables
    
  }
  return(alpha_k)
  
  
}

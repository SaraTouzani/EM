alphai <-
function(Y, Xquali){
  #repartition des observations de variables qualitatives suivant
  Xquali_k<-as.list(unique(Y))
  if(is.null(dim(Xquali))) { # une seule variable 
    j<-1
    for(i in unique(Y)) {
      Xquali_k[[j]]<-data.frame(MAT = Xquali[Y==i])
      j<-j+1
    }
    
    alpha_k<-as.list(unique(Y))
    for (k in 1:length(alpha_k)) {
      variables<-vector("list",1)
      X<-as.numeric(table( Xquali_k[[k]] )/ nrow(Xquali_k[[k]]))
      variables[[1]]<-data.frame(proba=X,
                                 modalite=unlist(attributes(table( Xquali_k[[k]] )/ nrow(Xquali_k[[k]]))$dimnames))
      for (h in unique(Xquali)) {
        if(!(h %in% variables[[1]]$modalite)){
          mod_miss<-data.frame(proba=0,
                               modalite=h)
          variables[[1]]<-rbind( variables[[1]],mod_miss)
        }
      }
      variables[[1]]$modalite<-as.character( variables[[1]]$modalite)
      variables[[1]]<- variables[[1]][order( variables[[1]]$modalite , decreasing = F),]
      
      alpha_k[[k]]<- variables
      
    } 
  }
  else{
    j<-1
    for(i in unique(Y)) {
      Xquali_k[[j]]<-data.frame(MAT = Xquali[Y==i,])
      j<-j+1
    }
    alpha_k<-as.list(unique(Y))
    for (k in 1:length(alpha_k)) {
      variables<-as.list(colnames(Xquali))
      for (j in 1:length(variables)) {
        X<-as.numeric(table( Xquali_k[[k]][,j] )/ nrow(Xquali_k[[k]]))
        variables[[j]]<-data.frame(proba=X,
                                   Var=variables[[j]],
                                   modalite=unlist(attributes(table( Xquali_k[[k]][,j] )/ nrow(Xquali_k[[k]]))$dimnames))
        for (h in unique(Xquali[,j])) {
          if(!(h %in% variables[[j]]$modalite)){
            mod_miss<-data.frame(proba=0,
                                 Var=unique(variables[[j]]$Var),
                                 modalite=h)
            variables[[j]]<-rbind(variables[[j]],mod_miss)
          }
          
        }
        variables[[j]]$modalite<-as.character(variables[[j]]$modalite)
        variables[[j]]<-variables[[j]][order(variables[[j]]$modalite , decreasing = F),]
        
        
      }
      alpha_k[[k]]<-variables
      
    }  
    
  }
  return(alpha_k)
}

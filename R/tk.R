tk <-
function(fkx,pk){
   K<-length(pk)
   n<-length(fkx)
   tk<-matrix(0,ncol = K,nrow = n)
   for (i in 1:n) {
     tk[i,]<-(as.numeric(unlist(fkx[[i]]))*pk)/(as.numeric(unlist(fkx[[i]]))%*%pk)[1,1] 
     
   }
   return(tk)
 }

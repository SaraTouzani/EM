fx <-
function(Pk,fkx,Y){
   f<-matrix(0, nrow = length(Y) )
   for (i in 1:length(Y)) {
     f[i]<- as.numeric(unlist(fkx[[i]]))%*%Pk
   }
   
   return(f)
 }

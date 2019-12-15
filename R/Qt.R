Qt <-
function(tkx,fkx,Pk){
  n<-nrow(tkx)
  K<-ncol(tkx)
  somme<-0
    for (i in 1:n) {
#      somme<-somme + as.numeric(unlist(fkx[[i]])%*%tkx[i,]) 
      somme<-somme+((log1p(Pk)+log1p(unlist(fkx[[i]])))%*%tkx[i,])
    }
   return(somme)
}

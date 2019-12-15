sigma_init <-
function(data,Y){
  sigma<-as.list(unique(Y))
  j<-1
  data$class<-Y
  for(i in unique(Y)) {
    Mat_cov<-(cov(data[data$class==i, -match("class", names(data))][,]))
    sigma[[j]]<-data.frame(MAT = Mat_cov)
    j<-j+1
  }
  return(sigma)
}

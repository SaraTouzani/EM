init_class <-
function(data,k){
  p<-vector("list",1)
  p[[1]][1:k]<- (1/k)
  return(sample(k,nrow(data), replace = TRUE, prob =unlist(p) ))
}

Plot_quanti <-
function(Dataquanti,Em_results,all=TRUE){
  
  tkx <-Em_results[[length(Em_results)]][["part_MAP"]]
  Y<-apply(tkx, MARGIN = 1 ,which.max)
  if(all){
    pairs.panels(Dataquanti,gap=0, bg=Y,pch = 21)
  }
  else{
    i<-1
    while(i< ((ncol(Dataquanti))-6) ){
      pairs.panels(Dataquanti[,i:(i+5)],gap=0, bg=Y,pch = 21)
      i<-i+6
    }
    pairs.panels(Dataquanti[,i:ncol(Dataquanti)],gap=0, bg=Y,pch = 21)
    
  }

  
  }

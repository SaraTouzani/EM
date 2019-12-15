Em_quali <-
function( Xquali,  alpha_ki=NULL,Xijh_i=NULL, pki=NULL,Y=NULL, methode="init",e=0.1,nbr_iteration=NULL,k) {
  #  Yobs<-Y[!is.na(Y)]
  #  Ymis<-Y[is.na(Y)]
  #  n<-length(c(Yobs,Ymis))
  #  o<-length(Yobs)
  if(methode=="kmodes"){
    Y<-kmodes(Xquali,k)$cluster
  }
  else{
    Y<-init_class(iris,k)
  }
  if(is.null(alpha_ki)){
    alpha_ki<-alphai(Y,Xquali)
  }
  if(is.null(Xijh_i)){
    Xijh_i<-xijh(Xquali)
  }
  if(is.null(pki)){
    pki<-as.numeric(table(Y)/length(Y))
  }
  Pki <- pki
  Alphaki <- alpha_ki
  p=length(colnames(Xquali))
  n<-nrow(Xquali)
  fkx<-fk_quali(Y,Xijh_i,alpha_ki, p=p, n = n)  
  #fkx<-fk_quali(Y,Xijh_i,alpha_ki)
  f<-fx(Pki,fkx,Y)
  tkx<-tk(fkx,Pki)  
  li <- l(f)
  Qi<-Qt(tkx,fkx,Pki)
  BIC<-Bic(f,Y,Xquali = Xquali)
  ICL<-BIC-sum(tkx*log(tkx))
  i <- 1
#  teta_init<-list(iteration=i,
#                  Proba=Pki,
#                  Alphak = Alphaki,
#                  log_like=li,
#                  Q=Qi,
#                  Bic=BIC,
#                  ICL=ICL,
#                  proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
#                  part_MAP=tkx  
#  )
  #plot(Xquanti$Sepal.Length,Xquanti$Sepal.Width, col = Y, pch = 19)
  #Y[1]<-which.max(tkx[1,])
  teta<-list()
  K<-k
 # teta[[1]]<-teta_init
#  x11()
  repeat {
    #layout(matrix(1:6, ncol =K ,nrow=2))
 #   vars<-colnames(Xquali)
#    #par(mfrow=c(length(Xquali),1))
 #   for (j in 1:length(Xquali)) {
  #    XQ_K<-NULL
   #   for(k in unique(Y)){
    #      XQ_K<-rbind(XQ_K,table(Xquali[Y==k,j]))
          #pie(table(Xquali[Y==k,j]),main =paste(c("Variable:", vars[j]), collapse=" ") )
         # text(-1, 0, paste(c("  Freq:", table(Xquali[Y==k,j])), collapse=" "), col = "black")
     # }
      #barplot(XQ_K,beside = FALSE,legend.text = TRUE,xlab = paste(c("Variable:", vars[j]), collapse=" ") )
    #}

    Y<-apply(tkx, MARGIN = 1 ,which.max)
    Nk<-nk(tkx)
    Pk<-pk(Nk,Xquali)
    Alphak <- alphak(tkx,Xijh_i,Y,Xquali)
    #   tk <- tk(pk, alphak)
    fkx<-fk_quali(Y,Xijh_i,Alphak, p=p, n = n)
    f<-fx(Pk,fkx,Y)
    lf<-l(f)
    tkx<-tk(fkx,Pk)  
    Qf<-Qt(tkx,fkx,Pk)
    BIC<-Bic(f,Y,Xquali = Xquali)
    ICL<-BIC-sum(tkx*log(tkx))
    teta[[i]]<-list(iteration=i,
                      Proba=Pk,
                      Alphak = Alphak,
                      log_like=lf,
                      Q=Qf,
                      Bic=BIC,
                      ICL=ICL,
                      proba_individu=fkx,
                      part_MAP=tkx          
    )
    if(!is.null(nbr_iteration)){
      if(i>=nbr_iteration) break
    }
    if (abs(li - lf) < e)  {
      break
    }
    else{

      Pki <- Pk
      Alphaki <- Alphak
      
      li <- lf
      i<-i+1
    }
    
    
  }
  
  # Alphak = Alphak   ,               #   tk <- tk(pk, alphak)
  #fkx<-fk(Xquanti, Meank, Sigk, Y, Xijh_i , Alphak, p=ncol(Xquali) )
  #f<-fx(Pk,fkx,Y)
  return(teta)
}

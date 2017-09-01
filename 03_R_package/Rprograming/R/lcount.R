lcount <-
function(lst,obj,index=FALSE){
  
  if (index!=FALSE){
    suma<-sum(sapply(lst[[index]],function(x){length(as.character(which(x==obj)))}))
  
  }
  else if (index==FALSE) {
    
    suma=sum(sapply(lst,function(x){length(as.character(which(x==obj)))}))
    
  }
  return(suma)
}

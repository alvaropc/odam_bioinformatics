lremove <-
function(lst,obj,index=FALSE){
  
  if (index>0){
    lst[[index]]<-unlist(sapply(lst[[index]],function(x){x[x!=obj]}))
  }
  else if (index==FALSE){
    lst<-sapply(lst,function(x){x[!x==obj]})
  }
  return(lst)
}

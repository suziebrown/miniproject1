## generates sequence of colours e.g. for "image" plots
## colour sequence goes from "basecol" to white
## rev=TRUE reverses order of colour sequence

suzie.colors <- function(n,basecol="purple",rev=FALSE){
  #stopifnot(n>0,n<257)
  cols<-character(n)
  basecol <- paste(as.hexmode(col2rgb(basecol)),sep="",collapse="")
  for (i in 1:(n-1)){
    alpha <- as.hexmode((i*256)%/%n)
    if (nchar(paste(alpha))==1){alpha<-paste("0",alpha,sep="")}
    cols[i]<-toupper(paste("#",basecol,alpha,sep=""))
  }
  cols[n]<-toupper(paste("#",basecol,"ff",sep=""))
  if(rev==TRUE){rev(cols)}
  else{cols}
}

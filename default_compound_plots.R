source("K_calcs_Johnson_OS.R")

TSSolPlot<-function(compound){
S_list<-c(0,17.5,32,33,34,35,36)
plot(T_LOTS,KH(compound,T_LOTS,max(S_list)),type="n",xlab="Temperature / Celcius", ylab="KH / dimensionless (gas/liquid)")
linetype<-1
legendlist<-NULL
for(S in S_list){
 lines(T_LOTS,KH(compound,T_LOTS,S),lty=linetype)
 linetype<-linetype+1
 legendlist<-c(legendlist,paste("S=",S,sep=""))
}

legend("topleft",legendlist,lty=seq(from=1,to=length(S_list),by=1))
}

TSSolPlot<-function(compound){
S_list<-c(0,17.5,32,33,34,35,36)
plot(T_LOTS,KH_Molar_per_atmosphere(compound,T_LOTS,min(S_list)),type="n",xlab="Temperature / Celcius", ylab="KH / dimensionless (gas/liquid)")
linetype<-1
legendlist<-NULL
for(S in S_list){
 lines(T_LOTS,KH_Molar_per_atmosphere(compound,T_LOTS,S),lty=linetype)
 linetype<-linetype+1
 legendlist<-c(legendlist,paste("S=",S,sep=""))
}

legend("topright",legendlist,lty=seq(from=1,to=length(S_list),by=1))
}

TSSchmidtPlot<-function(compound){
S_list<-c(0,17.5,32,33,34,35,36)
plot(T_LOTS,schmidt(compound,T_LOTS,max(S_list)),type="n",xlab="Temperature / Celcius", ylab="Schmidt number (dimensionless)")
linetype<-1
legendlist<-NULL
for(S in S_list){
 lines(T_LOTS,schmidt(compound,T_LOTS,S),lty=linetype)
 linetype<-linetype+1
 legendlist<-c(legendlist,paste("S=",S,sep=""))
}

legend("topleft",legendlist,lty=seq(from=1,to=length(S_list),by=1))
}

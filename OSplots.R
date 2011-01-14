# plots for OS numerical scheme paper

source("K_calcs_Johnson_OS.R")

#Make smoother curves by having more winspeed values

WS_LOTS <- c(0.1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)
T_LOTS <- c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)



f01 <- function(){
	#compare various ka parameterizations	
	compound <- "CHI3"
	NDuce<-360000*Duce_ka(compound,WS_LOTS)
	NLiss<-360000*Liss_ka(WS_LOTS)
	NMackay<-360000*MackayYeun_ka(compound,WS_LOTS, 15)
	NShahin<-360000*Shahin_ka(compound,WS_LOTS,15)	
	compound2 <- "O2"
	CDuce<-360000*Duce_ka(compound2,WS_LOTS)
	CLiss<-360000*Liss_ka(WS_LOTS)
	CMackay<-360000*MackayYeun_ka(compound2,WS_LOTS, 15)
	CShahin<-360000*Shahin_ka(compound2,WS_LOTS,15)
	xrange<-c(0,30)
	yrange<-c(0,36000)
	
	pdf("f01.pdf", paper="special", width=16.6, height=8.3)

	par(mfrow=c(1,2), mar=c(6,6,2,2))
	plot(xrange,yrange,type="n",
		xaxs="i", yaxs="i", lwd=2,tcl=0.5, cex.axis=1.3, cex.lab=1.8,
		xlab=expression(paste("wind-speed / m ", s^{-1})),
		ylab=expression(paste(k[a]," / cm ", hr^{-1}))
	)
	
	lines(WS_LOTS,CLiss,lty=1,lwd=3, col="black")
	lines(WS_LOTS,CMackay,lty=4,lwd=3, col="firebrick3")
	lines(WS_LOTS,CDuce,lty=5,lwd=3, col="firebrick3")
	lines(WS_LOTS,CShahin,lty=3,lwd=4, col="firebrick3")
	#lines(COAREWINDS,CCOARE,lty=6,lwd=2, col="firebrick3")
	lines(WS_LOTS,NMackay,lty=4,lwd=3, col="forestgreen")
	lines(WS_LOTS,NDuce,lty=5,lwd=3, col="forestgreen")
	lines(WS_LOTS,NShahin,lty=3,lwd=4, col="forestgreen")
	#lines(COAREWINDS,NCOARE, lty=6, lwd=2,col="forestgreen")

	legend(1,34500,c("Liss 1973","Mackay and Yeun 1983","Duce et al. 1991","Shahin et al. 2002","O2","CHI3"),lty=c(1,4,5,3,1,1),col=c("black","black","black","black","firebrick3","forestgreen"))
	xrange<-c(0,7)
	yrange<-c(0,5400)

	plot(xrange,yrange,type="n",
		xaxs="i", yaxs="i", lwd=2,tcl=0.5, cex.axis=1.3, cex.lab=1.8,
		xlab=expression(paste("wind-speed / m ", s^{-1})),
		ylab=expression(paste(k[a]," / cm ", hr^{-1}))
	)
	
	lines(WS_LOTS,CLiss,lty=1,lwd=3, col="black")
	lines(WS_LOTS,CMackay,lty=4,lwd=3, col="firebrick3")
	lines(WS_LOTS,CDuce,lty=5,lwd=3, col="firebrick3")
	lines(WS_LOTS,CShahin,lty=3,lwd=4, col="firebrick3")
	#lines(COAREWINDS,CCOARE,lty=6,lwd=2, col="firebrick3")
	lines(WS_LOTS,NMackay,lty=4,lwd=3, col="forestgreen")
	lines(WS_LOTS,NDuce,lty=5,lwd=3, col="forestgreen")
	lines(WS_LOTS,NShahin,lty=3,lwd=4, col="forestgreen")
	#lines(COAREWINDS,NCOARE, lty=6, lwd=2,col="forestgreen")

	#title(expression(paste("Gas phase transfer velocity, ", k[a]," with wind speed")))
	dev.off()
}


f02<-function(T=15){
    #compare Duce etal 1991 ka with schmidt number approach
    pdf("f02.pdf", paper="special", width=16.6,height=8.3)
    #COAREWINDS<-c(0.52589,1.13269,2.1845,3.1553,4.2071,6.6343,8.9401,12.6618,17.4757,21.2783,25.0)
    compound <-"O2"
    compound2 <-"CHI3"
    Ca<-Duce_standard<-Duce_ka(compound,WS_LOTS)*360000
    Cb<-Duce_schmidt<-Duce_ka_with_Sc(compound,WS_LOTS,1,T)*360000
    Cc<-Duce_schmidt_Smith_Cd<-Duce_ka_with_Sc(compound,WS_LOTS,2,T)*360000
    #Cd<-Duce_schmidt_Largepond_Cd<-Duce_ka_with_Sc(compound,WS_LOTS,3,T)*100
    Ce<-Mackay_Yeun<-MackayYeun_ka(compound,WS_LOTS,T)*360000
    #Cf<-c(0.0597,0.1215,0.2274,0.3259,0.4412,0.6905,0.9376,1.3680,1.9851,2.4845,2.9983)

    Na<-Duce_standard<-Duce_ka(compound2,WS_LOTS)*360000
    Nb<-Duce_schmidt<-Duce_ka_with_Sc(compound2,WS_LOTS,1,T)*360000
    Nc<-Duce_schmidt_Smith_Cd<-Duce_ka_with_Sc(compound2,WS_LOTS,2,T)*360000
    #Nd<-Duce_schmidt_Largepond_Cd<-Duce_ka_with_Sc(compound2,WS_LOTS,3,T)*100
    Ne<-Mackay_Yeun<-MackayYeun_ka(compound2,WS_LOTS,T)*360000
    #Nf<-c(0.0761,0.1548,0.2897,0.4152,0.5621,0.8797,1.1945,1.7429,2.5291,3.1654,3.8200)
    xrange<-c(0,30)
    yrange<-c(0,20000)
    par(mfrow=c(1,2), mar=c(6,6,2,2))

    plot(xrange,yrange,type="n",xaxs="i",yaxs="i", cex.axis=1.3, cex.lab=1.8,
        lwd=2,tcl=0.5,xlab=expression(paste("wind speed / m ", s^{-1})),
        ylab=expression(paste(k[a]," / cm ", hr^{-1}))
    )
    lines(WS_LOTS,Ca,lty=1,lwd=2, col="firebrick3")
    lines(WS_LOTS,Cb,lty=3,lwd=2, col="firebrick3")
    lines(WS_LOTS,Cc,lty=2,lwd=2, col="firebrick3")
    #lines(WS_LOTS,Cd,lty=3,lwd=3, col="firebrick3")
    lines(WS_LOTS,Ce,lty=4,lwd=2, col="firebrick3")
    #lines(COAREWINDS,Cf,lty="92",lwd=2, col="firebrick3")
    lines(WS_LOTS,Na,lty=1,lwd=2, col="forestgreen")
    lines(WS_LOTS,Nb,lty=3,lwd=2, col="forestgreen")
    lines(WS_LOTS,Nc,lty=2,lwd=2, col="forestgreen")
    #lines(WS_LOTS,Nd,lty=3,lwd=3, col="forestgreen")
    lines(WS_LOTS,Ne,lty=4,lwd=2, col="forestgreen")
    #lines(COAREWINDS,Nf,lty="92",lwd=2, col="forestgreen")
        
    legend(1,19400,c("Duce/MW","Duce/Sc","Duce/Smith","Mackay and Yeun 1983","O2","CHI3"),lty=c("solid","11","44","1343","solid","solid"),col=c("black","black","black","black","firebrick3","forestgreen"))
    xrange<-c(0,7)
    yrange<-c(0,3600)

 plot(xrange,yrange,type="n",xaxs="i",yaxs="i", cex.axis=1.3, cex.lab=1.8,
        lwd=2,tcl=0.5,xlab=expression(paste("wind speed / m ", s^{-1})),
        ylab=expression(paste(k[a]," / cm ", hr^{-1}))
    )
    lines(WS_LOTS,Ca,lty=1,lwd=2, col="firebrick3")
    lines(WS_LOTS,Cb,lty=3,lwd=2, col="firebrick3")
    lines(WS_LOTS,Cc,lty=2,lwd=2, col="firebrick3")
    #lines(WS_LOTS,Cd,lty=3,lwd=3, col="firebrick3")
    lines(WS_LOTS,Ce,lty=4,lwd=2, col="firebrick3")
    #lines(COAREWINDS,Cf,lty="92",lwd=2, col="firebrick3")
    lines(WS_LOTS,Na,lty=1,lwd=2, col="forestgreen")
    lines(WS_LOTS,Nb,lty=3,lwd=2, col="forestgreen")
    lines(WS_LOTS,Nc,lty=2,lwd=2, col="forestgreen")
    #lines(WS_LOTS,Nd,lty=3,lwd=3, col="forestgreen")
    lines(WS_LOTS,Ne,lty=4,lwd=2, col="forestgreen")
    #lines(COAREWINDS,Nf,lty="92",lwd=2, col="forestgreen")
    dev.off()
}



# u_* (friction velocity) is related to windspeed through the drag coefficient, C_D. C_D = (u*/u)^2. 
# rearranging, u* = sqrt(C_d)*u

# in Duce et al, C_D is constant at 1.3e-3

Duce_u_star<-function(u){
	u_star<-u*sqrt(1.3e-03)
	u_star
} 

# Mackay and Yeun use the parameterisation of Smith (1980), J. Phys. Oceanog,10,709

Smith_u_star<-function(u){
	u_star<-u*sqrt(1e-4*(6.1+0.63*u))
	u_star
}

# Large and Pond (1981), J. Phys Oceanog 11, 324 also present one:

Large_Pond_u_star<-function(u){
	C_d<-ifelse(u<11,1.2,0.49 + (0.065*u))
	u_star=u*sqrt(C_d*1e-3)
	u_star
}


f03<-function(T=15){
    #compare new duce term with Jeffery and Mackay and Yeun
    pdf("f03.pdf", paper="special", width=16.6,height=8.3)
    compound <-"O2"
    compound2 <-"CHI3"
    Cc<-(new_ka(compound,WS_LOTS,T)-1e-3)*360000
    Cd<-Jeffrey_ka(compound,WS_LOTS,T)*360000
    Ce<-Mackay_Yeun<-MackayYeun_ka(compound,WS_LOTS,T)*360000

    Nc<-(new_ka(compound2,WS_LOTS,T)-1e-3)*360000
    Nd<-Jeffrey_ka(compound2,WS_LOTS,T)*360000
    Ne<-Mackay_Yeun<-MackayYeun_ka(compound2,WS_LOTS,T)*360000
    xrange<-c(0,30)
    yrange<-c(0,20000) 
    par(mfrow=c(1,2), mar=c(6,6,2,2))

    plot(xrange,yrange,type="n",xaxs="i",yaxs="i", cex.axis=1.3, cex.lab=1.8,
        lwd=2,tcl=0.5,xlab=expression(paste("wind speed / m ", s^{-1})),
        ylab=expression(paste(k[a]," / cm ", hr^{-1}))
    )

    lines(WS_LOTS,Cc,lty=1,lwd=2, col="firebrick3")
    lines(WS_LOTS,Cd,lty="3732", lwd=3, col="firebrick3")
    lines(WS_LOTS,Ce,lty=4,lwd=2,col="firebrick3")

    lines(WS_LOTS,Nc,lty=1,lwd=2, col="forestgreen")
    lines(WS_LOTS,Nd,lty="3732", lwd=3, col="forestgreen")
    lines(WS_LOTS,Ne,lty=4, lwd=2, col="forestgreen")
 
    legend(1,19400,c("Duce et al. (1991) with Smith(1980) CD","Jeffrey et al. (2010)","Mackay and Yeun (1983)","O2","CHI3"),lty=c("solid","3732","1343","solid","solid"),col=c("black", "black","black","firebrick3","forestgreen"))

    xrange<-c(0,7)
    yrange<-c(0,3600) 

    plot(xrange,yrange,type="n",xaxs="i",yaxs="i", cex.axis=1.3, cex.lab=1.8,
        lwd=2,tcl=0.5,xlab=expression(paste("wind speed / m ", s^{-1})),
        ylab=expression(paste(k[a]," / cm ", hr^{-1}))
    )
    lines(WS_LOTS,Cc,lty=1,lwd=2, col="firebrick3")
    lines(WS_LOTS,Cd,lty="3732", lwd=3, col="firebrick3")
    lines(WS_LOTS,Ce,lty=4,lwd=2,col="firebrick3")

    lines(WS_LOTS,Nc,lty=1,lwd=2, col="forestgreen")
    lines(WS_LOTS,Nd,lty="3732", lwd=3, col="forestgreen")
    lines(WS_LOTS,Ne,lty=4, lwd=2, col="forestgreen")
 
    dev.off()
}
 

#fairall 2003 Cd10N from graph:
Fairall_winds<-c(0.485829959514, 0.890688259109, 1.25506072874,	1.65991902834, 2.30769230769, 3.15789473684, 4.08906882591, 5.78947368421, 8.17813765182, 9.91902834008, 12.2672064777, 17.5708502024, 18.1376518219, 24.9392712551, 30)
Fairall_Cd<-c(0.000553888035739,0.000405855487525,0.000344054765561,0.000294544384835,0.000257199696855,0.000244360902256,0.000280733332004,0.000402838994037,0.000586097205879,0.00072049320915,	0.000989983247243,0.00162711154544,0.00168833888435,0.00229991424184,0.00275) 

#compare all this with the Jeffrey et al 2010 implementation of the NOAA COARE algorithm

#first with CD as and input
Jeffrey_with_CD<-function(compound,u,T=15,Cd){
	von_kar<-0.4
	#Cd<-1e-4*(6.1+0.63*u)
	Sc<-Sc_air(compound,T)
	ra<-13.3*sqrt(Sc) + (Cd^(-0.5)) - 5 + log(Sc)/(2*von_kar)
	u_star<-sqrt(Cd*(u^2))
	u_star/ra
}

#now applying Smith(1980) Cd
Jeffrey_ka<-function(compound,u,T=15){
	von_kar<-0.4
	Cd<-1e-4*(6.1+0.63*u)
	Sc<-Sc_air(compound,T)
	ra<-13.3*sqrt(Sc) + (Cd^(-0.5)) - 5 + log(Sc)/(2*von_kar)
	u_star<-sqrt(Cd*(u^2))
	u_star/ra
}


#compare Jeffrey implementation with Fairall 2003 Cd with Jeffery/Smith

f04<-function(){
	pdf("f04.pdf", paper="special",width=16.6, height=8.3)
	counter<-1
	Jeffrey_fairall<-NULL
	for(i in Fairall_winds){
		Jeffrey_fairall<-c(Jeffrey_fairall,360000*Jeffrey_with_CD("CO2",i,15,Fairall_Cd[counter]))
		counter<-counter+1
	}
	Jeffrey_smith<-360000*Jeffrey_ka("CO2",Fairall_winds,15)
	#plot(Fairall_winds,Jeffrey_fairall,type="l",xlim=c(0,7),ylim=c(0,1))
	par(mfrow=c(1,2), mar=c(6,6,2,2))
	plot(c(0,30),c(0,20000),type="n",xaxs="i",yaxs="i", cex.axis=1.3, cex.lab=1.8,
        lwd=2,tcl=0.5,xlab=expression(paste("wind speed / m ", s^{-1})),
        ylab=expression(paste(k[a]," / cm ", hr^{-1}))
    	)
	lines(Fairall_winds,Jeffrey_fairall,lty=3, lwd=6)
	lines(Fairall_winds,Jeffrey_smith, lty="3732", lwd=2)
	lines(WS_LOTS,(new_ka("CO2",WS_LOTS,15)-1e-3)*360000,lty=1, lwd=2)
	legend(2,19400,c(expression(paste("Jeffrey et al. (2010) [COARE] with Smith (1980) ", C[D])),expression(paste("Jeffrey et al. (2010) [COARE] with Fairall et al. (2003) ", C[D])),expression(paste("Duce et al. (1991) with Smith (1980) ", C[D]))),lty=c("3237","11","solid"),lwd=c(2,4,2))

	plot(c(0,7),c(0,3600),type="n",xaxs="i",yaxs="i", cex.axis=1.3, cex.lab=1.8,
        lwd=2,tcl=0.5,xlab=expression(paste("wind speed / m ", s^{-1})),
        ylab=expression(paste(k[a]," / cm ", hr^{-1}))
    	)
	lines(Fairall_winds,Jeffrey_fairall,lty=3,lwd=6)
	lines(Fairall_winds,Jeffrey_smith, lty="3732", lwd=2)
	lines(WS_LOTS,(new_ka("CO2",WS_LOTS,15)-1e-3)*360000,lty=1,lwd=2)
	dev.off()
}


f05<-function(){
	data<-read.csv("coarecomp.dat")
	Nightingale<-360000*kw("CO2",15,WS_LOTS,35,normalize=660)
	Woolf<-Woolf97kw("CO2",15,WS_LOTS,35,normalize=660)*360000
	WannLT<-0.39*(WS_LOTS^2)
	WannST<-0.31*(WS_LOTS^2)
	McGillis<-3.3+0.026*(WS_LOTS^3)
 	xrange<-c(0,20)
	yrange<-c(0,100)
	pdf("f05.pdf", paper="special", width=16.6, height=8.3)
	
	par(mfrow=c(1,2), mar=c(6,6,2,2))	
	plot(xrange,yrange,type="n",
		xaxs="i", yaxs="i", lwd=2,tcl=0.5, cex.axis=1.3, cex.lab=1.8,
		xlab=expression(paste("wind speed / m.", s^{-1})),
		ylab=expression(paste(k[w_660]," / cm.", hr^{-1}))
	)
	lines(WS_LOTS,Nightingale,lty=2,lwd=4,col="dodgerblue")
	lines(WS_LOTS,Woolf,lty=1,lwd=3,col="firebrick")
	lines(WS_LOTS,WannLT,lty=2,lwd=4,col="forestgreen")
	lines(WS_LOTS,WannST,lty=4,lwd=4,col="forestgreen")
	#liss + merlivat
	lines(data[,35],data[,36],lty=2,lwd=4,col="darkslateblue")
	#Jeffrey coare
	lines(data[,33],data[,34],lty=1,lwd=3,col="darkblue")
	#gasex98 tuned COARE
	lines(data[,3],data[,4],lty=1,lwd=3,col="olivedrab4")
	#Mcgillis cubic relationship
	lines(WS_LOTS,McGillis,lty=2,lwd=4,col="darkgoldenrod4")
	#nightingale data
	points(data[,37],data[,38],pch=3,col="black")
	#gasex 98 observations
	points(data[,1],data[,2],pch=2,cex=1.5,col="black")
	#blomquist DMS observaitons
	points(data[,7],data[,8],pch=15,cex=0.8)	
	#Huebert DMS observations
	points(data[,21],data[,22],pch=1)
	#Miller knorr a + b observations
	points(data[,31],data[,32],pch=7)
	#Watston 1991 dual tracer
	points(data[,41],data[,42],pch=16)
	#Kholouiski 1995 radon data averaged
	points(c(10,10),c(22.1,15.8),pch=23,cex=1.6,bg="black")
	#Sweeney 2007 global estimate
	points(6.9,14.6,pch=21,cex=2,bg="goldenrod")
	#Wanninkhof 1992 global estimate
	points(7.4,21.9,pch=21,cex=2,bg="forestgreen")	
	#Mueller 2007 global estimate
	points(6.9,17.1,pch=21,cex=2,bg="burlywood4")
	legend(0.5,98,c("Nightingale et al. (2000)",
			"Woolf (1997)",
			"Wanninkhof (1992) Long-term",
			"Wanninkhof (1992) Short-term",
			"Liss and Merlivat (1983)", 
			"Jeffrey et al., (2010) Tuned NOAA COARE", 
			"Hare et al., (2004) GasEx '98 tuned COARE algorithm",
			"McGillis et al., (2001) GasEx '98 empirical fit (cubic)",  
			"Nightingale et al., (2000) dual tracer",
			"Hare et al., (2004) GasEx '98 CO2 fluxes",
			"Blomquist et al., (2006) DMS fluxes", 
			"Huebert et al., (2010) DMS fluxes", 
			"Miller et al., (2009) CO2 fluxes",
			"Watson et al., (1991) dual tracer",
			"Kholouiski (1995) averaged radon data",
			"Sweeney et al., (2007) 14C global estimate",
			"Wanninkhof (1992) 14C global estimate",
			"Müller et al., (2008) 14C global estimate"
			),
		lty=c(2,1,2,4,2,1,1,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1),
		pch=c(-1,-1,-1,-1,-1,-1,-1,-1,3,2,15,1,7,16,23,21,21,21),
		col=c(	"dodgerblue",
			"firebrick",
			"forestgreen",
			"forestgreen", 
			"darkslateblue",
			"darkblue",
			"olivedrab4",
			"darkgoldenrod4",
			"black",
			"black",
			"black",
			"black",
			"black",
			"black",
			"black",
			"black",
			"black",
			"black"),
		lwd=c(3,2,3,3,3,2,2,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1),
		pt.bg=c("white","white","white","white","white","white","white","white","white","white","white","white","white","black","black","goldenrod","forestgreen","burlywood4"),
		bty="n")
	xrange<-c(0,7)
	yrange<-c(0,30)
	plot(xrange,yrange,type="n",
		xaxs="i", yaxs="i", lwd=2,tcl=0.5, cex.axis=1.3, cex.lab=1.8,
		xlab=expression(paste("wind speed / m.", s^{-1})),
		ylab=expression(paste(k[w_660]," / cm.", hr^{-1}))
	)
	lines(WS_LOTS,Nightingale,lty=2,lwd=4,col="dodgerblue")
	lines(WS_LOTS,Woolf,lty=1,lwd=3,col="firebrick")
	lines(WS_LOTS,WannLT,lty=2,lwd=4,col="forestgreen")
	lines(WS_LOTS,WannST,lty=4,lwd=4,col="forestgreen")
	#liss + merlivat
	lines(data[,35],data[,36],lty=2,lwd=4,col="darkslateblue")
	#Jeffrey coare
	lines(data[,33],data[,34],lty=1,lwd=3,col="darkblue")
	#gasex98 tuned COARE
	lines(data[,3],data[,4],lty=1,lwd=3,col="olivedrab4")
	#Mcgillis cubic relationship
	lines(WS_LOTS,McGillis,lty=2,lwd=4,col="darkgoldenrod4")
	#nightingale data
	points(data[,37],data[,38],pch=3,col="black")
	#gasex 98 observations
	points(data[,1],data[,2],pch=2,cex=1.5,col="black")
	#blomquist DMS observaitons
	points(data[,7],data[,8],pch=15,cex=0.8)	
	#Huebert DMS observations
	points(data[,21],data[,22],pch=1)
	#Miller knorr a + b observations
	points(data[,31],data[,32],pch=7)
	#Watston 1991 dual tracer
	points(data[,41],data[,42],pch=16)
	#Kholouiski 1995 radon data averaged
	points(c(10,10),c(22.1,15.8),pch=23,cex=1.6,bg="black")
	#Sweeney 2007 global estimate
	points(6.9,14.6,pch=21,cex=2,bg="goldenrod")
	#Wanninkhof 1992 global estimate
	points(7.4,21.9,pch=21,cex=2,bg="forestgreen")	
	#Mueller 2007 global estimate
	points(6.9,17.1,pch=21,cex=2,bg="burlywood4")
	
 dev.off()
}

f06<-function(){
	C<-c(7.66E-004,-7.16E-004,6.04E-004,-2.94E-004,3.96E-004,5.68E-004,3.50E-004,8.30E-004,8.50E-004,8.69E-004,8.80E-004,9.19E-004,7.61E-004,3.26E-004,4.04E-004,4.14E-004,6.28E-004,6.37E-004,7.81E-004,7.98E-004,7.46E-004,8.83E-004,6.80E-004,6.36E-004,5.96E-004,3.64E-004,3.96E-004,5.58E-004,6.60E-004,8.03E-004)

	lnKH<-c(2.82,-14.6,-1.36,-11.2,-5.84,-5.76,-6.6,0.31,4.51,3.38,4.17,3.45,2.14,-6.6,-5.84,-7.11,-1.99,0.31,-1.59,-1.3,3.45,5.14,0.53,0.16,-2.5,-4.61,-4.26,-3.89,-3.67,-3.2
)

	testlnK<-(-15:5)

	fitC<-(7.3353282561828962e-04+(3.3961477466551352e-05*testlnK) +(-2.4088830102075734E-06*testlnK^2)+(1.5711393120941302E-07*testlnK^3))


	pdf("f06.pdf",paper="special",width=8.3,height=8.3)
		plot(c(-16,6),c(-1e-3,1e-3),type="n",xaxs="i",yaxs="i",ylim=c(-1e-3,1e-3),xlab=expression(paste("ln",K[H])),ylab=expression(theta),tcl=0.5,cex.axis=1.3, cex.lab=1.8)
		points(lnKH,C)
		points(-7.291,-5.26E-005,pch=20,ylim=c(-1e-3,1e-3))
		lines(testlnK,fitC)
	dev.off()
}


f07_old<-function(T=15){
 pdf("f07.pdf", width=8.3, height=8.3)
 DMS_data<-read.table("huebertDMS.dat")
 par(mar=c(6,6,2,2))
 plot(c(0,15),c(0,60),type="n",tcl=-0.5, cex.axis=1.3, cex.lab=1.8, xaxs="i",yaxs="i",xlab="u10 / ms-1",ylab="kw660 or Kw(660) / cm hr-1")


 #liss+merlivat
 data<-read.csv("coarecomp.dat")
 lines(data[,35],data[,36],lty=3,lwd=3)

 #NOAA COARE DMS with bubbles from Huebert 2010 
 #cubic fit to DOGEE data from Mingxi Yang (Barry Huebert's group)
 coef_a<-2.2177674432959229E+00
 coef_b<-6.5708740314358582E-01
 coef_c<-5.5949520530126154E-02
 coef_d<-1.3651758186355976E-03
 x_data<-data[,"HuebertCOAREDMSbubws"]
 y_data<-coef_a + coef_b*x_data + (coef_c*(x_data^2)) + (coef_d*(x_data^3))

 lines(x_data,y_data,lwd=3)

 #NOAA COARE CO2 with bubbles from Huebert 2010
 lines(read.table("COARECO2Bubbles.dat"),lty=2,lwd=3)

 lines(WS_LOTS,360000*Nightingkw("DMS",T,WS_LOTS,35,normalize=660),lty=3, lwd=3,col="firebrick")
 lines(WS_LOTS,360000*Wannkw("DMS",T,WS_LOTS,35,normalize=660),lty=3,lwd=3,col="darkslateblue")
 lines(WS_LOTS,360000*Woolf97kw("DMS",T,WS_LOTS,35,normalize=660),lty=3,lwd=3,col="forestgreen")
 total_K_N2000<-360000/((1/(kw("DMS",T,WS_LOTS,35,normalize=660)))+(1/(KH("DMS",T,35)*ka("DMS",WS_LOTS,T)))) 
 total_K_W1997<-360000/((1/(Woolf97kw("DMS",T,WS_LOTS,35,normalize=660)))+(1/(KH("DMS",T,35)*ka("DMS",WS_LOTS,T))))
total_K_W1997_CO2<-360000/((1/(Woolf97kw("CO2",T,WS_LOTS,35,normalize=660)))+(1/(KH("CO2",T,35)*ka("CO2",WS_LOTS,T))))
 lines(WS_LOTS,total_K_N2000,lty=4, lwd=4, col="firebrick")
 lines(WS_LOTS,total_K_W1997,lty=4, lwd=4, col="forestgreen")
 lines(WS_LOTS,total_K_W1997_CO2,lty=5, lwd=4, col="forestgreen")
 #Huebert 2010 DMS
 points(DMS_data,pch=21)

 #Hare et al gasex CO2 data
 points(data[,1],data[,2],pch=25,col="black",bg="black")

 legend(0.5,58.5,bty="n",c("Heubert et al 2010 eddy covariance DMS derived kw660","Hare et al 2004 eddy covariance CO2 derived kw660","Liss and Merlivat 1986 kw660","Wanninkhof 1992 kw660","Nightingale et al 2000 kw660","Woolf 1997 bubble-mediated kw660 for DMS","NOAA COARE DMS kw660 with bubbles (from Huebert et al 2010)","NOAA COARE CO2 kw660 with bubbles (from Huebert et al 2010)","This scheme Kw with Nightingale et al 2000 kw660","This scheme Kw with Woolf 1997 DMS kw660","This scheme Kw with Woolf 1997 CO2 kw660"),pch=c(21,25,-1,-1,-1,-1,-1,-1,-1,-1,-1),lty=c(0,0,3,3,3,3,1,2,4,4,5),lwd=c(0,0,3,3,3,3,2,3,3,3,3),col=c("black","black","black","darkslateblue","firebrick","forestgreen","black","black","firebrick","forestgreen","forestgreen"),bg=c("white","black","white","white","white","white","white","white","white","white","white"))
 dev.off()
}

f07<-function(T=15){
 pdf("f07.pdf", width=8.3, height=8.3)
 DMS_data<-read.table("huebertDMS.dat")
  par(mar=c(6,6,2,2))
 plot(c(0,15),c(0,60),type="n",tcl=-0.5, cex.axis=1.3, cex.lab=1.8, xaxs="i",yaxs="i",xlab=expression(paste("wind-speed / m ", s^{-1})),ylab=expression(paste(k[w][660]," or ",K[w][(660)]," / cm h",r^{-1})))



 #NOAA COARE DMS with bubbles from Huebert 2010
 #cubic fit to DOGEE data from Mingxi Yang (Barry Huebert's group)
 coef_a<-2.2177674432959229E+00
 coef_b<-6.5708740314358582E-01
 coef_c<-5.5949520530126154E-02
 coef_d<-1.3651758186355976E-03
 data<-read.csv("coarecomp.dat")
 x_data<-data[,"HuebertCOAREDMSbubws"]
 y_data<-coef_a + coef_b*x_data + (coef_c*(x_data^2)) + (coef_d*(x_data^3))

 lines(x_data,y_data,lwd=3)

 #NOAA COARE CO2 with bubbles from Huebert 2010
 lines(read.table("COARECO2Bubbles.dat"),lty=2,lwd=3)

 lines(WS_LOTS,360000*Nightingkw("DMS",T,WS_LOTS,35,normalize=660),lty=3, lwd=3,col="firebrick")

 lines(WS_LOTS,360000*Woolf97kw("DMS",T,WS_LOTS,35,normalize=660),lty=1,lwd=3,col="lightblue")

 total_K_W1997<-360000/((1/(Woolf97kw("DMS",T,WS_LOTS,35,normalize=660)))+(1/(KH("DMS",T,35)*ka("DMS",WS_LOTS,T))))
 total_K_W1997_CO2<-360000/((1/(Woolf97kw("CO2",T,WS_LOTS,35,normalize=660)))+(1/(KH("CO2",T,35)*ka("CO2",WS_LOTS,T))))

 lines(WS_LOTS,total_K_W1997,lty=1, lwd=4, col="forestgreen")
 lines(WS_LOTS,total_K_W1997_CO2,lty=2, lwd=4, col="forestgreen")
 #Huebert 2010 DMS
 points(DMS_data,pch=21)

 #Hare et al gasex CO2 data
 data<-read.csv("coarecomp.dat") 
 points(data[,1],data[,2],pch=25,col="black",bg="black")

 legend(0.5,58.5,bty="n",c("Heubert et al 2010 eddy covariance DMS derived kw660","Hare et al 2004 eddy covariance CO2 derived kw660","Nightingale et al 2000 kw660","Woolf 1997 bubble-mediated kw660 for DMS","NOAA COARE DMS kw660 with bubbles (from Huebert et al 2010)","NOAA COARE CO2 kw660 with bubbles (from Huebert et al 2010)","Johnson 2010 Kw with Woolf 1997 DMS kw660","Johnson 2010 Kw with Woolf 1997 CO2 kw660"),pch=c(21,25,-1,-1,-1,-1,-1,-1),lty=c(0,0,3,1,1,2,1,2),lwd=c(1,0,3,3,2,3,3,3),col=c("black","black","firebrick","lightblue","black","black","forestgreen","forestgreen"),pt.bg=c("white","black","white","white","white","white","white","white"))
 dev.off()
}


f08<-function(){
 source("NH3comp.r")
 pdf("f08.pdf",width=8.3,height=8.3)
 par(mar=c(6,6,2,2))
 plot(newflux,oldflux,pch=20,xaxs="i",yaxs="i",ylab=expression(paste("Johnson et al. 2008 flux / pmol ",m^{-2}," ",s^{-1})),
        xlab=expression(paste("This scheme flux / pmol ",m^{-2}," ",s^{-1})),xlim=c(-60,20),ylim=c(-60,20),tcl=0.5,cex.axis=1.3, cex.lab=1.8,)
 #arrows(newflux, oldminflux, newflux, oldmaxflux, code=3, angle=90, length=0.1)
 oneone<-(-60:20) 
 lines(oneone,oneone)
 dev.off()
}

f01()
f02()
f03()
f04()
f05()
f06()
f07()
f08()


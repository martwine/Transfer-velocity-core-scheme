source("K_calcs_Johnson_OS.R")

# Chemical enhancement calculations based on Hoover and Berkshire (1969) equation



#############################################
#Calculating chemical enhancement to kw
###############################################


# get mass boundarly layer dpeth from windspeed for given compound
z_from_wind<- function(compound,ws,T,S){
	#calculate z in cm given compound and conditions
	diff<-0.5*(diff_HM(compound,T,S)+diff_WC(compound,T,S))
	diff/(kw(compound,T,ws,S)*100)
}

alpha_kw <- function (compound,ws,T,S,khyd=1e-100,tau=1.0000000001){
        #khyd is pseudo first-order rate constant of reaction in water
	# z is apparent mass boundary layer depth (see Liss and Slater 1974)
        z <- z_from_wind(compound,T,ws,S)
	x <- z*sqrt(khyd*tau/Diff)
	alpha <- (tau)/((tau-1) + (tanh(x)/x))
    	alpha
}

#define a range of tau (minus 1) for generic calculation 
tau_minus_one_list <- c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100)

#a_test <- function(khyd,tauminusone,Diff,ws,T=15, S=35){
#	tau <- tauminusone + 1
#	z <- z_from_wind("O2",ws,T,S)
#	x <- z*sqrt(khyd*tau/Diff)
#	(tau*x)/((tauminusone*x) + tanh(x))
#}


# this function retuns the residual between a given target value of alpha and the actual one calculated 
# given the input params - for use in the uniroot solver below, where all but khyd are contsrained.
a_general <- function(khyd,tauminusone,Diff,z,target_value){
	tau <- tauminusone + 1
	x <- z*sqrt(khyd*tau/Diff)
	target_value-(tau/(tauminusone + (tanh(x)/x)))
}

# use an iterative solver to find the value of khyd for a given enhancement factor
# Diffusivity is compound specific and needs to be fed into this function appropriately
solve_a_general_for_rate <- function(tauminusone,Diff,z,target){
	uniroot(a_general,c(1e-300,1e300),tol=1e-100,tauminusone=tauminusone,Diff=Diff,z=z,target_value=target,maxiter=10000)$root
}




#generic table of k values for a given threshold enhancement over a range of tau and z. This is not compound-specific by default, as it isn't particularly compound-sensitive (at the order of magnitude level, which appropriate for this analysis). Use O2 as default typical average gas.

zeds <- function(ws=WS, compound="O2", T=15, S=35){
# get z values in cm from standard windspeeds, by default using O2 as model gas (after Liss and Slater 1974)
	z_from_wind(compound,ws,T,S)
}

khyd_for_given_enhancement_table<-function(target, Diff=1e-5, ws=WS, T=15, S=35, tauminusonelist=tau_minus_one_list){
	#use t and WS (or user-specified vectors of temp and windspeed) to calculate the total transfer velocity (Kw) for a given compound over these ranges of conditions
	zedds<-zeds(ws,T=T,S=S)	
	tabler<-NULL
	#need to modify tau-minus-one list so max enhancement = tau/tau-1 can always be satisfied
	taucutter<-NULL
	for (t in tauminusonelist) {
		if (((t+1)/t)>target){
			taucutter<-c(taucutter,t)
		}
	}	
	tauminusonelist<-tauminusonelist[1:length(taucutter)]
	for(i in tauminusonelist){
		rowtemp<-NULL
		for(j in zedds){
			rowtemp<-c(rowtemp,solve_a_general_for_rate(i,Diff,j,target))
		}
	tabler<-rbind(tabler,rowtemp)	
	}
	output_table<-data.frame(tabler,row.names=(tauminusonelist))
	names(output_table)<-ws
	output_table
}


#####################################################################
#   Calculate chemical enhancement threshold values in gas phase    #
#####################################################################

alpha_ka <- function (compound,ws,T,S,katm=1e-100,tau=1.00000000001){
        #katm is rate constant of reaction in gas phase
	Diff <- D_air(compound,T)
        z <- Diff/(100*ka(compound,ws,T))
	x <- z*sqrt(katm*tau/Diff)
	alpha <- (tau)/((tau-1) + (tanh(x)/x))
    	alpha
}

zg_from_wind<- function(compound,ws,T){
	#calculate z in cm given compound and conditions
	D_air(compound,T)/(ka(compound,ws,T)*100)
}


zgeds <- function(ws=WS, compound="O2", T=15){
# get z values in cm from standard windspeeds using O2 as model gas (after Liss and Slater 1974)
	zg_from_wind(compound,ws,T)
}


#generic table of k (rate constant) values for a given threshold enhancement over a range of tau and z. This is not compound-specific by default, as it isn't particularly compound-sensitive (at the order of magnitude level, which appropriate for this analysis). Use O2 as default typical average gas.

khyd_for_given_gas_phase_enhancement_table<-function(target, Diff=0.1, ws=WS, T=15, tauminusonelist=tau_minus_one_list){
	#use t and WS (or user-specified vectors of temp and windspeed) to calculate the total transfer velocity (Kw) for a given compound over these ranges of conditions
	zedds<-zgeds(ws,T=T)	
	tabler<-NULL
	#need to modify tau-minus-one list so max enhancement = tau/tau-1 can always be satisfied
	taucutter<-NULL
	for (t in tauminusonelist) {
		if (((t+1)/t)>target){
			taucutter<-c(taucutter,t)
		}
	}	
	tauminusonelist<-tauminusonelist[1:length(taucutter)]
	for(i in tauminusonelist){
		rowtemp<-NULL
		for(j in zedds){
			rowtemp<-c(rowtemp,solve_a_general_for_rate(i,Diff,j,target))
		}
	tabler<-rbind(tabler,rowtemp)	
	}
	output_table<-data.frame(tabler,row.names=(tauminusonelist))
	names(output_table)<-ws
	output_table
}




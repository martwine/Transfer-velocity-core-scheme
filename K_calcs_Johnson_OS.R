# R program for calculating the liquid phase and gas phase transfer velocities for any volatile compound given the relevant physical and chemical data. For scientific details see Johnson, M.T., A numerical scheme to calculate the temperature and salinity dependent air-water transfer velocity for any gas. Ocean Science, 2010. Questions / comments gladly received: martin.johnson@uea.ac.uk.

## first of all read in the data file with each of the gases of interest
# note: if Vb is specified, then this value is used in diffusion coefficient calculations rather than the Schroeder method
# KH and tVar shoudl be taken from Rolf Sander's compilation of Henry's lay constants (http://www.mpch-mainz.mpg.de/~sander/res/henry.html) if available

compounds<-read.table("compounds.dat", header = TRUE, row.names = 1)


#let's define a range of temperature and windspeed over which to calculate transfer velocities etc by default

t = c(-5,0,5,10,15,20,25,30,35)
WS = c(0.1,0.5,1,2,5,10,20,30)


##########################################################
#              MOLAR VOLUME AT BOILING POINT             #
##########################################################
 


Vb <- function(compound){
	## Calculate the molar volume at boiling point, Vb, using the Schroeder method, or just take the overrulling value from compounds.dat if specified
	# note: for compounds containing elements other than C,H,O,N,S, Cl, Br, F, and I, the schroeder method won't work so Vb must be specified
	# db, tb and rings should be the number of double and triple bonds, and ring features in the moelcule respectively. 
	ringval<-ifelse(compounds[compound,"rings"]>0,-7,0)
	ifelse(compounds[compound,"Vb"]>0,compounds[compound,"Vb"],7*(compounds[compound,"C"]+compounds[compound,"H"]+compounds[compound,"O"]+compounds[compound,"N"]+compounds[compound,"db"])+14*compounds[compound,"tb"]+31.5*compounds[compound,"Br"]+24.5*compounds[compound,"Cl"]+10.5*compounds[compound,"F"]+38.5*compounds[compound,"I"]+21*compounds[compound,"S"]+32*compounds[compound,"Se"]+ringval)	

}


##########################################################
##		DRAG coefficient C_D			##
##########################################################

smithCD<-function(u){
	1e-4*(6.1+0.63*u)
}

C_D<-function(u,formulation=1){
	if (formulation==1) {smithCD(u)} else {-999}
}

############################################################
############################################################
##      LIQUID PHASE TRANSFER VELOCITY CALCULATIONS       ##
############################################################
############################################################


####################################################
###       Viscosity and density of water         ###
####################################################

n_0 <- function(T){
	#Calculate the viscosity of pure water in cP (mPa s) according to LaLiberte 2007 (doi:10.1021/je0604075). Requires T in celcius
    	(T+246)/(137.37+(5.2842*T)+(0.05594*(T^2)))
}



n_sw <- function(T,S,return_relative_viscotiy=0){
	#calculate the dynamic viscosity of seawater in cP using the viscosity model / mixing rule approach of LaLiberte 2007
	#to return relative viscosity (i.e. viscosity scaling factor) for comparison with other studies, pass n_sw a final argument of "1"
	#after T and S e.g. n_sw(10,35,1)
    	#read in a data file containing the mass fraction of each constituent solute in seawater per salinity unit in per mil and the coefficients determined by LaLiberte 2007 for each salt
    	sw_cmf<-read.table("sw_constituent_mass_fractions.dat", header = TRUE, row.names = 1)
    	#w_i_ln_n_i_tot is the sum of the product of the mass fraction and the natural log of the viscosities contributed by each solute individually (see LaLiberte 2007, equation 8)
    	w_i_ln_n_i_tot <-0
    	#sum up the mass fractions to get water mass fraction
    	w_i_tot<-0
    	for (salt in row.names(sw_cmf)){
       		w_i <- sw_cmf[salt,"massfraction"]*S/1000
       		w_i_tot<-w_i_tot + w_i
    	}
    	for (salt in row.names(sw_cmf)){
       		w_i <- sw_cmf[salt,"massfraction"]*S/1000
       		v1 <- sw_cmf[salt,"v1"]
       		v2 <- sw_cmf[salt,"v2"]
       		v3 <- sw_cmf[salt,"v3"]
       		v4 <- sw_cmf[salt,"v4"]
       		v5 <- sw_cmf[salt,"v5"]
       		v6 <- sw_cmf[salt,"v6"]
       		#wi_tot is used here as eq 12 in LaLiberte requires (1-w_w), which is equivalent. Using this term seems initially counterintuitive - one might expect to use w_i here for each solute individually. However,  "at any given solute concentration, the solute viscosity increases with the total solute content" so 1-w_w (or w_i_tot) is the correct term - pers. comm. from Marc LaLiberte
       		n_i<-(exp(((v1*w_i_tot^v2)+v3)/((v4*T) + 1)))/((v5*(w_i_tot^v6))+1)
       		w_i_ln_n_i_tot <- w_i_ln_n_i_tot + (w_i*log(n_i))
    	}
    	ln_n_m<- (1-w_i_tot)*log(n_0(T))+w_i_ln_n_i_tot
    	if (return_relative_viscotiy==0) exp(ln_n_m) else
        	(exp(ln_n_m)/n_0(T))
}

n_sw_hardy <- function(T,n_zero=1){
	#Non-salinity dependent Hardy  1953 method for calculating viscosity - not used in calculations, purely for comparison with Laliberte method.
	# η = η0 * K/[1+0.03338.T + 0.00018325.T2]
	# use n_zero = 1 to use Laliberte n_0 or n_zero=0 for Hardy n_0
    	if (n_zero==0) n0<- 1.787 else
        	n0 <- n_0(T)
    	(n0*1.052/(1+(0.03338*T)+(0.00018325*(T^2))))
}

p_sw <- function(T,S){
	## Calculate the density of seawater in kg/m3 according to millero and Poisson (1981)
    	#coefficients for millero calc
    	A <- 0.824493-(0.0040899*T)+(0.000076438*(T^2))-(0.00000082467*(T^3))+(0.0000000053875*(T^4))
    	B <- -0.00572466+(0.00010277*T)-(0.0000016546*(T^2))
    	C <- 0.00048314
    	# density of pure water
    	p_o <- 999.842594+(0.06793952*T)-(0.00909529*(T^2))+(0.0001001685*(T^3))-(0.000001120083*(T^4))+(0.000000006536332*(T^5))
    	# return salinity-and-temperature-dependent density of seawater (at atmospheric pressure)
    	(p_o+(A*S)+(B*(S^(3/2)))+(C*S))
} 

v_sw <- function(T,S) {
	#calculate kinmatic viscosity of seawater in cm/s for Schmidt number calculation
    	# dynamic viscosity - convert to S.I. (Kg/m/s)
    	n = n_sw(T,S)/1000
    	# density already in S.I.
    	p = p_sw(T,S)
    	#multiply by 10000 to go from m2/s to cm2/s
    	10000*n/p
}

#################################################
#all diffusion coefficients calculated in cm^2/s#
#################################################

diff_HM <- function(compound,T,S){
	#Hayduk and Minhas (1982) diffusion coefficient calculation
    	EpsilonStar <- (9.58/Vb(compound))-1.12
    	1.25e-8*(Vb(compound)^(-0.19)-0.292)*((T+273.15)^1.52)*((n_sw(T,S))^EpsilonStar)
}

schmidt_HM <- function(compound,T,S){
	#calculate schmidt number from HM diffusivity
	(v_sw(T,S))/diff_HM(compound,T,S)
}

diff_HL <- function(compound,T,S){
	#calculate diffudivity by Hayduk and Laudie (1974) method
	#NOTE - only T dependence from n_sw calculation
	13.26e-5/(((n_sw(T,S))^1.4)*(Vb(compound)^0.589))
}

schmidt_HL <- function(compound,T,S){
	#calculate schmidt number from HM diffusivity
    	(v_sw(T,S))/diff_HL(compound,T,S)
}

diff_WC <- function(compound,T,S){
	#Wilkie and Chang (1955) diffusion coefficient
    	# associaction factor of solvent (2.6 in the case of water according to Poling 2001; although Wanninkhof suggests 2.26)
    	phi <- 2.6
    	((T+273.15)*7.4e-8*(phi*18.01)^0.5)/((n_sw(T,S))*(Vb(compound)^0.6))
}

schmidt_WC <- function(compound,T,S){
	#calculate schmidt number from WC diffusivity
    	(v_sw(T,S))/diff_WC(compound,T,S)
}

schmidt <- function(compound,T,S){
	#calculate mean schmidt number in water
    	mean_diff<-0.5*(diff_WC(compound,T,S)+diff_HM(compound,T,S))
        v_sw(T,S)/mean_diff
}

schmidt_out<-function(compound,T,S,schmidt=0){
	#output schmidt number T dependence for a list of T
	#calculated with mean schmidt number by default, or specify 1 for Hayduck and Minhas or 2 for Wilkie and Chang
    	if (schmidt==0) schmidt(compound,T,S) else
    		if (schmidt==1) schmidt_HM(compound,T,S) else
    			schmidt_WC(compound,T,S)
}



#########################################################
###   liquid phase transfer velcocity calculations    ###
###           all values returned in m/s              ###
#########################################################


Wannkw <- function (compound,T,u,S,normalize=0,schmidt_formulation=0){
	#calculate kw transfer velocity in m/s according to Wanninkhof 1992
	#specifying a 'normalize' value allows the calculation of the transfer velocity for a user-specified value of schmidt number
    	if (normalize!=0) schmidt_number<-normalize else
    	schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
    	(0.31*u^2*((schmidt_number/660)^(-0.5)))/(100*3600)
}

Lissmerkw<- function(compound,T,u,S,normalize=0,schmidt_formulation=0){
		k600<-ifelse(u<3.6,0.17*u,
			ifelse(u<13,(2.85*u)-9.65,(5.9*u)-49.3))
	    	if (normalize!=0) schmidt_number<-normalize else
        		schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
		ifelse(u<3.6,(k600*((schmidt_number/600)^(-0.66)))/(100*3600),(k600*((schmidt_number/600)^(-0.5)))/(100*3600))
}

Lissmer_approx_kw <- function(compound,T,u,S,normalize=0,schmidt_formulation=0){
	#calculate kw in m/s according to liss and merlivat 
	# note k600 not k660
	# in order to make calculations easier, an exponential relationship has been fitted to the Liss and Merlivat Curve)
    	if (normalize!=0) schmidt_number<-normalize else
        	schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
    	(0.075*u^2.25*((schmidt_number/600)^(-0.5))/(100*3600))
}

Sweeneykw <- function (compound,T,u,S,normalize=0,schmidt_formulation=0){
	#calculate kw according to Sweeney 2007
    	if (normalize!=0) schmidt_number<-normalize else
    		schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
    	(0.27*u^2*((schmidt_number/660)^(-0.5)))/(100*3600)
}

Nightingkw <- function(compound,T,u,S,normalize=0,schmidt_formulation=0){
	# empirical fit to dual tracer data by Nightingale et al 2000 (GBC)
	# note k600 not k660 for their study
    	if (normalize!=0) schmidt_number<-normalize else
        	schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
    	(((0.222*u^2)+0.333*u)*(schmidt_number/600)^(-0.5))/(100*3600)
}

Woolf97kw <- function(compound,T,u,S,normalize=0,wb_formulation=0,schmidt_formulation=0){
	if (normalize!=0) schmidt_number<-normalize else
        	schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
    	(((0.222*u^2)+0.333*u)*(schmidt_number/600)^(-0.5))/(100*3600)
	# whitecapping component included in this physically based model of gas exchange
	# W_b is whitecap coverage
	#different whitecapping paramaterisations could be applied here:
	#0 is woolf default from Monahan..
	if (wb_formulation==0) W_b<-3.84e-6 * (u^(3.41))
	
	# K_b is whitecap-dependent (ostensibly bubble-related) component of total k_w
	# beta is equal to liquid over gas unitless KH using water T for both phases (equal to 1/KH in this scheme).  
		
	beta<-1/KH(compound,T,S)
	denominator<-(beta*(1+((14*beta*(schmidt_number^(-0.5)))^(-1/1.2)))^1.2)	
	K_b = 2540*W_b/denominator
	# K_o is non-whitecapping friction velocity relationship component of k_w
	u_star<-u*sqrt(C_D(u))
	K_o<-1.57e-4*u_star*((600/schmidt_number)^0.5)
	((K_o*360000) + K_b)/360000
}

McGilliskw<-function(compound,T,u,S,normalize=0,schmidt_formulation=0){
	#empirical fit to GasEx data by McGillis et al (JGR 2001)
	schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
	((3.3+(0.026*(u^3)))*((schmidt_number/600)^(-0.5)))/(100*3600)
}

kw <- function(compound,T,u,S,normalize=0,schmidt_formulation=0) {
	#calculate mean kw. set normalize to schmidt number value e.g. 660 to compare with other kw curves. Doesn't matter what compound or value of T is passed to the function (unless using the Woolf1997 bubble scheme), only u and schmidt value to normalise to, but you have to pass it some values otherwise R will break.
    	sf<-schmidt_formulation
    	nmlz<-normalize
    	Nightingkw(compound,T,u,S,nmlz,sf)
	#Woolf97kw(compound,T,u,S,normalize=nmlz,schmidt_formulation=sf)
	#Wannkw(compoundT,u,S,normalize=nmlz,schmidt_formulation=sf)
}



########################################################
#              HENRY'S LAW COEFFICIENTS                #
########################################################

KH0 <- function(compound,T=25){
	#Calculate Henry's law constant at a given T in pure water according to Sander (1999)
    	12.2/((273.15+T)*(compounds[compound,"KH"])*exp((compounds[compound,"tVar"])*((1/(T+273.15))-(1/298.15))))
}


Ks<-function(compound){
	theta = (7.3353282561828962e-04 + (3.3961477466551352e-05*log(KH0(compound))) + (-2.4088830102075734E-06*(log(KH0(compound)))^2) + (1.5711393120941302E-07*(log(KH0(compound)))^3))
	theta*log(Vb(compound))
}

K_H_factor <-function(compound,S){
	#calculate salinity-dependent salting-out scaling factor for KH0 (see manuscript)
    	10^(Ks(compound)*S)
}

KH <- function(compound,T,S){
	#Calculate gas-over-liquid unitless Henry's law constant at a given T and Salinity
    	KH0(compound,T)*K_H_factor(compound,S)
}

KH_Molar_per_atmosphere <- function(compound,T,S){
	# calculate the Henry's law constant in M/atm from Sander data
	# applying the salting out factor above
    	(compounds[compound,"KH"]*exp((compounds[compound,"tVar"])*((1/(T+273.15))-(1/298.15))))/K_H_factor(compound,S)
}





########################################################
#       GAS PHASE TRANSFER VELOCITY CALCULATIONS       #
########################################################


Tucker_D_air <- function(compound,T){
	#calculate diffusivity in air in cm2/sec
	#M_a is molar weight of air
	M_a <- 28.97
	M_b <- compounds[compound,"mw"]
	M_r <- (M_a + M_b)/(M_a*M_b)
	#assume 1ATM
	P <- 1
	#assume molar volume air is 20.1 cm3/mol
	V_a <- 20.1	
	(0.001*((T+273.15)^1.75)*sqrt(M_r))/(P*((V_a^(1/3))+(Vb(compound)^(1/3))))^2
}

D_air <- function(compound,T){
	Tucker_D_air(compound,T)
} 

n_air <- function(T){
	# dynamic viscosity of saturated air according to Tsiligiris 2008
	SV_0 = 1.715747771e-5
	SV_1 = 4.722402075e-8
	SV_2 = -3.663027156e-10
	SV_3 = 1.873236686e-12
	SV_4 = -8.050218737e-14
	
	# in N.s/m^2 (Pa.s)
	u_m = SV_0+(SV_1*T)+(SV_2*T^2)+(SV_3*T^3)+(SV_4*T^4)
	u_m
}

p_air <- function(T){
	# density of saturated air according to Tsiligiris 2008 in kg/m^3
	SD_0 = 1.293393662
	SD_1 = -5.538444326e-3
	SD_2 = 3.860201577e-5
	SD_3 = -5.2536065e-7
	p = SD_0+(SD_1*T)+(SD_2*T^2)+(SD_3*T^3)
	p
}

v_air <- function(T) {
	#calculate kinmatic viscosity of air in cm2/s for Schmidt number calculation
    	# dynamic viscosity 
    	n = n_air(T)
    	# density 
    	p = p_air(T)
    	#multiply by 10000 to go from m2/s to cm2/s
    	10000*n/p
}

Sc_air <- function (compound,T){
	#calculate the schmidt number of a given gas in air
	v_air(T)/D_air(compound,T)
}


####################################################################
####      all ka paramterisations return values in m/s          ####
####################################################################

Duce_ka <- function(compound,u){
	#calculate ka gas phase transfer velocity according to Duce et al 1991
    	u/(770+(45*((compounds[compound,"mw"])^(1/3))))
}

Duce_ka_with_Sc<-function(compound,u,C_d_form,T=15){
	# use C_d = 1 > Duce; C_d = 2 > Smith; C_d = 3 > Large Pond
	if (C_d_form==1) u_star<-Duce_u_star(u)
	if (C_d_form==2) u_star<-Smith_u_star(u)
	if (C_d_form==3) u_star<-Large_Pond_u_star(u)
	ra_turb <-u/(u_star)^2
	ra_diff <- (5/u_star)*((Sc_air(compound,T))^(2/3))
	1/(ra_turb + ra_diff)
	
}


Liss_ka <- function(u){
	#calculate ka according to Liss 1973 (non compound-specific)
	(0.005+(0.21*u))/100
}

MackayYeun_ka <- function(compound,u,T)
	#ka according to Mackay and Yeun 1983
	1e-3 + (46.2e-5*sqrt(6.1+(0.63*u))*u*(Sc_air(compound,T))^-0.67)

Shahin_ka <- function(compound,u,T){
	#calculate transfer velocity from gas phase diffusivity in air according to Shahin et al 2002 in m/s
	#checked against shahin data - they find a ka of ~3 cm/s at 6m/s wind - approx 5 times that of Duce!)- see Johnson 2010 (final article) for an explanation
	(sqrt(D_air(compound,T))*((0.98*u)+1.26))/100
}

Jeffrey_ka<-function(compound,u,T){
	#using smith(1980) Cd
	von_kar<-0.4
	Cd<-(1e-4*(6.1+0.63*u))
	Sc<-Sc_air(compound,T)
	ra<-13.3*sqrt(Sc) + (Cd^(-0.5)) - 5 + log(Sc)/(2*von_kar)
	u_star<-u*sqrt(C_D(u))
	u_star/ra
}

new_ka<-function(compound,u,T){
#this is the ka propsed in Johnson 2010 (OS discussions paper), which has been superceded by a modified version of Jeffrey 2010 (scheme_ka below) for the resubmitted  paper which will hopefully be published in OS
 u_star<-u*sqrt(C_D(u))
 ra_turb <-u/(u_star)^2
 ra_diff <-(5/u_star)*((Sc_air(compound,T))^(2/3))
 katest<-1e-3 +1/(ra_turb + ra_diff)
 ifelse(is.na(katest),1e-3,katest)
}

Jeffrey_modified_ka<-function(compound,u,T){
	(1e-3+Jeffrey_ka(compound,u,T))
}

ka <- function(compound,u,T){
   #use new formulation, based on Jeffrey et al. (2010), by default
   Jeffrey_modified_ka(compound,u,T)
   
   #new_ka(compound,u,T)
   #Jeffrey_ka(compound,u,T)
   #Shahin_ka(compound,u,T)
   #MackayYeun_ka(compound,u,T)
   #Duce_ka_with_sc(compound,u,T)
   #Duce_ka(compound,u)
   #Liss_ka(u)
}

#####################################################################
#  RESISTANCE TO TRANSFER, RG/RL and TOTAL TRANSFER VECLOCITIES     #
#####################################################################

rl <- function(compound,T,u,S,schmidt_formulation=0){
	#calculate rl (=1/kw) after Liss and Slater 1974
    	1/kw(compound,T,u,S,0,schmidt_formulation)
}

rl_prime <- function(compound,T,u,S,schmidt_formulation=0){
	#rl_prime = KH/kw
    	KH(compound,T,S)/kw(compound,T,u,S,0,schmidt_formulation)
}

rg <- function(compound,T,u,S){
	#calculate rg (=1/KH*ka) after Liss and Slater
    	1/(KH(compound,T,S)*ka(compound,u,T))
}

rg_prime <- function(compound,u){
	#rg_prime = 1/ka 
    	1/ka(compound,u,T)
}

rg_by_rl <- function(compound,T,u,S){
	#calculate rg/rl to give relative contributions of the two phases to the total transfer velocity
    	rg(compound,T,u,S)/rl(compound,T,u,S)
}

rgp_by_rlp <- function(compound,T,u,S){
	#calculate rg_prime/rl_prime to give relative contributions of the two phases to the total transfer velocity.
	#if this and rg_by_rl aren't equal for a given gas and set of conditions, something's gone very wrong
    	rg_prime(compound,u)/rl_prime(compound,T,u,S)
}

Kw <- function(compound,T,u,S){
	#calculate total transfer velocity (Kw) after Liss and Slater
    	1/(rg(compound,T,u,S)+rl(compound,T,u,S))
}

Ka <- function(compound,T,u,S){
	#calculate total transfer velocity (Ka) after Liss and Slater
    	1/(rg_prime(compound,u)+rl_prime(compound,T,u,S))
}

#Sanity checking
Sanity <- function(){
cat("Sanity checking - make sure the numbers coming out are approximately what would be expected")
cat("They won't be exactly the same as the de-facto values as calculation methods, Henry's law and other input values etc will not be identical")
cat("\n\n")
cat("first of all, what about the Schmidt numbers for CO2 at 20 Celcius in freshwater and seawater? \n - the commonly quoted values are 600 and 660 respectively, but these are very sensitive so within +/- 30 of these values is fine") 
cat(paste("\nfreshwater: ", schmidt("CO2",20,0)))
cat(paste("\nseawater: ", schmidt("CO2",20,35)))
cat("\n\n")
cat("Hopefully that was OK. How about the values of some transfer velocities?")
cat("\n")
cat("for each formulation, the expected values at 4,10, and 16 m/s and schmidt number of 600 are given, and then the values calculated by this model \n")
cat("1: Liss and Merlivat (1986) \n")
cat("             4m/s    10m/s    16m/s \n")
cat("expected:    2       19       45.5  cm/hr \n")
cat(paste("calculated: ", round(360000*Lissmerkw("CO2",20,4,0,normalize=600),1),"   ", round(360000*Lissmerkw("CO2",20,10,0,normalize=600),1), "   ", round(360000*Lissmerkw("CO2",20,16,0,normalize=600),1)))

cat("\n \n 2: Wanninkhof(1992) \n")
cat("             4m/s    10m/s    16m/s \n")
cat("expected:    5.2     32.5     83   cm/hr \n")
cat(paste("calculated: ", round(360000*Wannkw("CO2",20,4,0,normalize=600),1),"   ", round(360000*Wannkw("CO2",20,10,0,normalize=600),1), "   ", round(360000*Wannkw("CO2",20,16,0,normalize=600),1)))

cat("\n \n 3: Nightingale et al (1999) \n")
cat("             4m/s    10m/s    16m/s \n")
cat("expected:    4.9     25.5     62   cm/hr \n")
cat(paste("calculated: ", round(360000*Nightingkw("CO2",20,4,0,normalize=600),1),"   ", round(360000*Nightingkw("CO2",20,10,0,normalize=600),1), "   ", round(360000*Nightingkw("CO2",20,16,0,normalize=600),1),"\n"))



}


#THE END


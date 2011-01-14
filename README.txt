Transfer velocity scheme - core module

author Martin Johnson (martin.johnson@uea.ac.uk)

file list:

K_calcs_Johnson_OS.R 
is the core transfer velocity file program as detailed in 
Johnson, M. T. A numerical Scheme to calculate the temperature and salinity dependent air-water transfer velocity for any gas, Ocean Science, 6, 913-932).

compounds.dat
provides input data for a list of compounds. Please extend as you see fit.

compounds.notes 
notes data sources for compounds in compounds.dat where appropriate

sw_constituent_mass_fractions.dat
the ionic make-up of average seawater, for viscosity calculations in core scheme. Actually the scheme isn't very sensitive to this - could be all NaCl and it wouldn't really make any difference. Still, it's in now so it might as well stay...

OSplots.R
code to reproduce all the plots in the Ocean Science manuscript where the core scheme is presented (http://www.ocean-sci.net/6/913/2010/os-6-913-2010.pdf)

coarecomp.dat
A whole load of kw-u10 relationships and tunes COARE models, all taken from papers or publicly available datasets. If there's anything here that shoudn't be public I apologise sincerely and will take it down immediately. These are mostly cribbed from graphs in papers, so they're good enough for making plots, but I wouldn't rely on them for numerical applications - better to go back to the original works. 

COARECO2Bubbles.dat
CO2 eddy covariance data from huebert et al 2010 (see OS manuscript)

huebertDMS.dat
DMS eddy covariance data from Huebert et al 2010

NH3comparison.csv
data from Johnson et al 2008. GBC paper on ammonia fluxes, calculated the 'old' way and using the Johnson 2010 scheme - for comparison. Used in figure plotting in OSplots.R


# CS Calculation

import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import scipy.interpolate
# import pdb; bp=pdb.set_trace

# setting the constants ##############

Mx=98.08;    #g/mol   # molar mass of H2SO4
Mair=28.965; #g/mol   #  - " -        air
Pr=1.0;   #atm        # pressure
Dair=19.7; # ??       # diffusion volume of air
Dx=51.96;  # ??       # diffusion volume of H2SO4
k = 8314.7;           # actually R
pi = np.pi

def CS_general(Dp,N,T,alpha,yesmatrix=[]):
	# function [CS,CS_prime,varargout] = CS_general(Dp,N,T,alpha,yesmatrix)
	#
	#
	# Version 0.1, last update: 02.7.2001
	#
	# This function calculates the total condensation sink
	# value and also the CS size distribution if yesmatrix is
	# nonzero.
	# Input:
	# * Dp		: diameter vector
	# * N 		: (absolute) concentration vector corrseponding to Dp
	#             NOTE: NOT dN/dlogDp !!!!!!!!
	# * T			: temperature
	# * alpha	: sticking coefficient
	# * yesmatrix: [optional] output the cs size distribution
	#
	# Output:
	# * CS		: the condensation sink [s^-1]
	# * CS_prime: CS / (4*pi*Dif)	[cm^-2]
	#
	# CS calculated as per Pirjola: Effects of aerosol dynamics...
	# (1998), J. Aerosol. Sci

	####################################################

	R = Dp/2;   # diameter to radius


	# The diffusion coefficient [cm^2/s]# from Reid et al.
	Temp = T;

	Dif = (0.001 * (Temp**1.75)*np.sqrt( (1/Mair)+(1/Mx))) / (Pr*(Dair**(1/3)+Dx**(1/3))**2);

	# lambda [m]#
	lam=3*(np.sqrt( (pi*Mx)/(8*k*Temp) )) * Dif *1e-4;

	# the knudsen number [dimensionless]#
	knud=lam/R;

	# beta, the correction coefficient [dimensionless]#
	beta=(knud+1)/((0.377*knud)+1+(4/(3*alpha))*(knud**2)+(4/(3*alpha))*knud);

	# the total CS calculation loop, cs in [1/s] #
	# CS_prime is given in cm^-2;
	# done this way to avoid total NaN:ing if only some values are NaN

	CS=0;
	CS_prime=0;
	CS_m = np.zeros(len(N))
	for j in range(len(N)):
		if (~np.isnan(N[j])):
			CS=CS + (4*pi*Dif)*N[j]*beta[j]*R[j]*1e2; # 1e2 comes for units
			CS_prime = CS_prime + N[j]*beta[j]*R[j]*1e2; # 1e2 comes for units

			CS_m[j] = 4*pi*Dif*N[j]*beta[j]*R[j]*1e2; # for the size distribution
			CS_m_prime = N[j]*beta[j]*R[j]*1e2; # for the size distribution

	# this outputs the size distribution if yesmatrix is given and nonzero

	if yesmatrix!=[]:
		if yesmatrix!=0:
			varargout = CS_m;
			varargout = varargout.append(CS_m_prime);
			return CS, CS_prime,varargout#, dCS_prime
	else:
		return CS, CS_prime#, dCS_prime

def CS_general_dlog(Dp,dN,T,alpha):
	# function [dCS, dCS_prime] = CS_general(Dp,dN,T,alpha)
	#
	#
	# Version 0.5, last update: 04.02.2004
	#(c) Miikka Dal Maso
	#  calculates dCS/dlogDp if dN/dlogDp is known;
	#
	# Input:
	# * Dp		: diameter vector
	# * N 		: particle concentration vector as dN/dlogDp
	# * T		: temperature
	# * alpha	: sticking coefficient (usually using 1.0)
	#
	# Output:
	# * dCS		: the condensation sink [s^-1]
	# * sCS_prime: CS / (4*pi*Dif)	[cm^-2]
	#
	# CS calculated as per Pirjola: Effects of aerosol dynamics...
	# (1998), J. Aerosol. Sci
	#
	# notes: n_CS(log Dp) = 4.*pi.*D.*beta(Dp).*(Dp./2).*n_N(log Dp)
	# so in this formula is the RADIUS!! not diameter even if the
	# size distribution is given as a function of diameter


	##########################

	R = Dp/2;   # diameter to radius

	# The diffusion coefficient [cm^2/s]# from Reid et al.
	Temp = T;
	Dif = (0.001 * (Temp**1.75)*np.sqrt( (1/Mair)+(1/Mx))) / (Pr*(Dair**(1/3)+Dx**(1/3))**2);

	# lambda [m]#
	lam=3*(np.sqrt( (pi*Mx)/(8.*k*Temp) )) * Dif *1e-4;

	# the knudsen number [dimensionless]#
	knud=lam/R;

	# beta, the correction coefficient [dimensionless]#
	beta=(knud+1)/((0.377*knud)+1+(4/(3.*alpha))*(knud**2)+(4/(3.*alpha))*knud);

	# the total CS calculation loop, cs in [1/s] #
	# CS_prime is given in cm^-2;
	# done this way to avoid total NaN:ing if only some values are NaN
	dCS = np.zeros(len(dN))

	for j in range(len(dN)):
		if (~np.isnan(dN[j])):
			dCS[j] = 4.*pi*Dif*dN[j]*beta[j]*R[j]*1e2; # for the size distribution
			dCS_prime = dN[j]*beta[j]*R[j]*1e2; # for the size distribution
		else:
			dCS[j] = 0; # this is not very good; however, the error is not too bad

	return dCS#, dCS_prime


def GF_gamma_lauri(dp_new):
# gamma as a function of Dp
	RH_norm=0.9; # [1:99]/100;

	dp = np.array([15, 20, 35, 50, 73, 109, 166, 264])*1e-9;
	GF_l= np.array([1.15,1.12,1.16,1.15,1.17,1.17,1.15,1.14]); # 1.16
	GF_i=np.array([1.31,1.25,1.22,1.24,1.32,1.36,1.32,1.46]);
	GF_m=np.array([1.38,1.32,1.36,1.32,1.37,1.46,1.53,1.59]); # 1.35

	osuus_l=np.array([ 61.5,57.2,39.8,38.7,46.1,53.5,38.3,28.0])/100; #61.5
	osuus_i=np.array([0.9,0.7,1.9,2.1,3.0,1.2,0.5,0.4])/100;
	osuus_m=np.array([ 37.6,42.1,58.3,59.2,50.9,45.3,61.2,71.6])/100; #37.5
	GF_mean=GF_l*osuus_l+GF_i*osuus_i+GF_m*osuus_m;

	gamma=np.log(GF_mean)/np.log(1-RH_norm);
	p=np.polyfit(dp,gamma,1);


	gamma_out=p[0]*np.minimum(dp_new,np.ones(len(dp_new))*2.9e-7)+p[1];
	#print(gamma_out)
	return gamma_out


def Dlog_to_N_vec(Dp, dN):

	lkm=len(Dp)
	vali=np.ones(lkm)
	vpist=np.zeros(lkm+1)

	for bi in range(1,lkm):
		vpist[bi]=Dp[bi-1]+0.5*(Dp[bi]-Dp[bi-1]);

	vpist[lkm]=Dp[lkm-1]+(Dp[lkm-1]-vpist[lkm-1])
	vpist[0]=Dp[0]-(vpist[0]-Dp[0])
	logpist=np.log10(vpist)
	vali[:]=np.diff(logpist) # dlogDp  jonkin verran approksimoiden

	N = dN*vali

	return Dp, N


def GFparam(Dp,dN,RH):
	# function[outDp, outdN] = GFparam(Dp,dN,RH)
	# Creates a humid size distribution for Hyytiälä using the parametrization
	# of Lauri Laakso et al, ACP, Page(s) 1933-1943. SRef-ID: 1680-7324/acp/2004-4-1933, 2004.
	# Rh over 99 are set to 99.
	# (c) Miikka Dal Maso


	dN[dN<0] = 0

	RH = min(RH,99)

	gamma = GF_gamma_lauri(Dp);

	GF = (1-(RH/100))**gamma;
	tr,N = Dlog_to_N_vec(Dp,dN);

	GF_Dp = Dp*GF;
	N_cum = np.cumsum(N);

	# with dense interpolation

	dens = np.logspace(np.log10(min(Dp)/10),np.log10(max(Dp)*10),100); # dense
	dgamma = GF_gamma_lauri(dens);
	densGF =   (1-(RH/100))**(dgamma);
	densN = sc.interpolate.interp1d(Dp,dN, fill_value=0, bounds_error=False,)(dens);

	# padding at the start and end to avoid NaNs

	densN[densN!=densN] = 0

	# cumulative distribution

	tr, densN_real = Dlog_to_N_vec(dens,densN);
	densN_real_cum = np.cumsum(densN_real);
	GFdens = dens*densGF;

	# getting the dN back

	newdN = np.diff(densN_real_cum)/np.diff(np.log10(dens));
	newdNGF = np.diff(densN_real_cum)/np.diff(np.log10(GFdens));

	newdens = np.zeros(len(newdN))
	newGFdens = np.zeros(len(newdN))
	for i in range(len(newdN)):
		newdens[i] = np.sqrt(dens[i]*dens[i+1])
		newGFdens[i] = np.sqrt(GFdens[i]*GFdens[i+1])

	outDp = np.logspace(np.log10(min(Dp)),np.log10(max(GF_Dp)*2),len(Dp));

	outdN = sc.interpolate.interp1d(newGFdens,newdNGF)(outDp);
	# checking the number matching:
	oldNtot = integrate_distribution(Dp,dN,1e-9,100e-6);
	newNtot = integrate_distribution(outDp,outdN,1e-9,100e-6);
	if oldNtot != 0:
		badness = (abs(oldNtot-newNtot)/oldNtot)
	else:
		badness=0
	if badness>0.05:
		print('Big np.difference in concentrations! Badness = #3.3g##\n',badness)
	# if this number gets big, something is wrong with the integration!
	return (outDp, outdN)

def integrate_distribution(Dp,dN,dmin,dmax):
	if dmin<min(Dp):
		dmin=min(Dp)
	if dmax>max(Dp):
		dmax=max(Dp);
	DpX = np.append(np.arange(dmin,dmax,1e-10),dmax);
	dNX = sc.interpolate.interp1d(Dp,dN, fill_value=0, bounds_error=False)(DpX)
	lDpX = len(DpX)
	return np.trapz(dNX, np.log10(DpX))

def log_to_lin(diam, sum):
	log = np.log10(diam)
	FF = (log[2:]-log[1:-1])/2 + (log[1:-1]-log[0:-2])/2
	dm_nconc_f = sum.copy()
	dm_nconc_f[:,1:-1] = dm_nconc_f[:,1:-1]*FF
	dm_nconc_f[:,0] = dm_nconc_f[:,0]*(log[1]-log[0])
	dm_nconc_f[:,-1] = dm_nconc_f[:,-1]*(log[-1]-log[-2])
	return dm_nconc_f

def filter_nan(mod_data, obs_data):
	"""
	this functions removed the data  from modeled and observed data
	whereever the observed data contains nan

	this is used by all other functions, otherwise they will produce nan as
	output
	"""
	data = np.array([mod_data.flatten(), obs_data.flatten()])
	data = np.transpose(data)
	data = data[~np.isnan(data).any(1)]
	return (data[:, 0], data[:, 1])


def CS_calc_day(v,RH,temp,meth):
	# calculating the condensation sink for a
	# dmps data file v
	# three methods are possible:
	# meth = 1: calculating by converting v to 'real'
	#           concentrations N(i), assuming the particles
	#           are all of size Dp(i), where i is the
	#           channel number
	# meth = 2: calculating the CS in dCS/dlogDp and integrating
	#           over the distribution (better)
	# meth = 3: using L. Laakso's parametrisation for Hyytiälä
	#           to correct v to ambient hygroscopicity
	# input
	# RH        : a 2-column vector with time and rh (1..100)
	#           in the columns [tim(:) rh(:)]
	# temp      : a 2-column vector with time and temperature
	#           in the columns [tim(:) temp(:)] temp is in celsius!
	#
	#           NOTE: RH and temp can also be scalar (1x1); this means
	#           a constant value is used for the whole day
	#           YOU MUST GIVE A VALUE FOR TEMP! RH is obligatory if meth = 3
	# v         : a matrix obtained by loading the DMPS datafile
	#           to the matlab workspace
	V = v.copy()
	V[1:,2:] = log_to_lin(v[0,2:], v[1:,2:])
	ro,co = np.shape(V);

	tim = v[1:,0];

	if meth>2:
		if all(~np.isnan(np.ravel(RH))):
			if len(RH)>1:
				rhi = sc.interpolate.interp1d(RH[:,0],RH[:,1],kind='linear',fill_value ='extrapolate')(tim);
			else:
				rhi = np.ones(np.shape(tim))*RH;
		else:
			rhi = ones(size(tim))*float('NaN')

	if all(~np.isnan(np.ravel(temp))):
		if len(temp)>1:
			tempi = sc.interpolate.interp1d(temp[:,0],temp[:,1],kind='linear',fill_value ='extrapolate')(tim)+273.15;
		else:
			tempi = np.ones(size(tim))*temp+273.15;

	else:
		tempi = np.ones(len(tim))*float('NaN');

	# checkin if temp might be in Kelvin:
		if any(tempi>373):
			print('Serious strangeness found in CS_calc_day:\n')
			print('Your ambient temperature is above the boiling point of water (T = #5.1f C)! Aborting...\n',max(temp));

	# RH you have to take care yourself... :-)
	CS = np.zeros(len(tim))
	if meth == 1:
	#  method 1: calculating the old way
		for i in range(len(tim)):
			Dp = V[0,2:];
			N  = V[i+1,2:];
			CS[i],_ = CS_general(Dp,N,tempi[i],1.0);

	elif meth ==2:
	#  method 2: using dCS/dlogDp
		for i in range(len(tim)):
			Dp = V[0,2:];
			dN  = v[i+1,2:];
			dCS = CS_general_dlog(Dp,dN,tempi[i],1.0);
			CS[i] = integrate_distribution(Dp,dCS,1e-9,50e-6);

	else:
	# method 3: correcting for growth factor
		Dp = v[0,2:];
		for i in range(len(tim)):
			dN  = v[i+1,2:];
			GF_Dp, GF_dN = GFparam(Dp,dN,rhi[i]);
			# if i==14:
			dCS = CS_general_dlog(GF_Dp,GF_dN,tempi[i],1.0);
			CS[i] = integrate_distribution(GF_Dp,dCS,1e-9,50e-6);

	return tim, CS

day = '180401'
CSfile = 'SMEAR_CS_T.dat'
RHfile = 'SMEAR_RH.dat'
TEMPfile = 'SMEAR_TEMP_168.dat'
selected_days = ([
'180401',
# '180402',
# '180403',
# '180404',
# '180405',
# '180406',
# '180407',
# '180408',
# '180409',
# '180410',
# '180411',
# '180412',
# '180413',
# '180414',
# '180415',
# '180416',
# '180417',
# '180418',
# '180419',
# '180420',
# '180421',
# '180422',
# '180423',
# '180424',
# '180425',
# '180426',
# '180427',
# '180428',
# '180429',
# '180430'
])
# plt.ion()
total_cs = []
total_time = []
for day in selected_days:
	# times need to be in same format
	v = np.genfromtxt('/home/pecl/04-MALTE/dMalte/Malte_in/Box/April2018/filled_dmps/dm'+day+'.sum')
	RH = np.genfromtxt('/home/pecl/04-MALTE/dMalte/Malte_in/Box/April2018/PC'+day+'/'+RHfile)[:144,:]
	temp = np.genfromtxt('/home/pecl/04-MALTE/dMalte/Malte_in/Box/April2018/PC'+day+'/'+TEMPfile)[:144,:]
	RH[:,1] = 80
	time,dd = CS_calc_day(v,RH,temp,3)
	total_cs = np.append(total_cs, dd)
	total_time = np.append(total_time, time)
	out = np.zeros((len(time)+2, 2))
	out[:144,0] = time
	out[:144,1] = dd
	out[144,0] = time[0]+1.
	out[144,1] = dd[-1]
	# np.savetxt('/home/pecl/04-MALTE/dMalte/Malte_in/Box/April2018/PC'+day+'/'+CSfile, out, fmt='%20.12e', delimiter = '  ')
	plt.plot(time,dd)
	plt.grid()

plt.show()

#

#!/usr/bin/python

#	'Copyright' 2012 Kamil Mielczarek (kamil.m@utdallas.edu), University of Texas at Dallas
#	Modifications:
#		10/2012 - Matlab code was ported to python, which is free and readily accessible to all :)
#                   3/2015 Improved accuracy of fundamental constants, fix plotting bugs on some platforms
#
#	Installation Notes:
#		Requires use of the 2.x Branch of python and
#		Matplotlib 	: http://matplotlib.org/
#		Scipy		: http://www.scipy.org/
#		Numpy		: http://numpy.scipy.org/
#                  Pylab             : wiki.scipy.org/PyLab
#
#		A free turn-key solution for windows/mac users is available via Enthought in the 'EPD Free'
#		package which is available via:
#		 http://www.enthought.com/products/epd_free.php
#
#		Linux users need to consult their distributions repositories for available packages


# Original MATLAB CODE WAS :
# Copyright 2010 George F. Burkhard, Eric T. Hoke, Stanford University

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

# This program calculates the field profile, exciton generation profile
# and generated current from the wavelength dependent complex indicies of
# refraction in devices using the transfer matrix method. It assumes the light source
# light source is in an n=1 environment (air) and that the first layer is
# a thick superstrate, so that incoherent reflection from the air/1st layer
# interface is taken into account before the coherent interference is
# calculated in the remaining layers. If there is no thick superstrate,
# input 'Air' as the first layer so that the incoherent reflection calculates
# to 0.
# The program
# also returns the calculated short circuit current for the device in
# mA/cm^2 if the device had an IQE of 100%.

# The procedure was adapted from J. Appl. Phys Vol 86, No. 1 (1999) p.487
# and JAP 93 No. 7 p. 3693.

# George Burkhard and Eric Hoke February 5, 2010
# When citing this work, please refer to:
#
# G. F. Burkhard, E. T. Hoke, M. D. McGehee, Adv. Mater., 22, 3293.
# Accounting for Interference, Scattering, and Electrode Absorption to Make
# Accurate Internal Quantum Efficiency Measurements in Organic and Other
# Thin Solar Cells

# Modifications:
# 3/3/11 Parastic absorption (parasitic_abs) is calculated and made
# accessable outside of script.

## USER CONFIGURATION

# Build up the device, light enters from left to right, the layer names
# correspond to the csv filenames in the 'matDir' ignoring the material prefix.
# example : you have a material file named 'nk_P3HT.csv' the layer name would be 'P3HT'
# Layer thicknesses are in nm.
layers 			= ['SiO2' , 'ITO' , 'PEDOT' , 'P3HTPCBM_BHJ' , 'Ca' , 'Al']
thicknesses		= [0 , 110 , 35  , 220 , 7 , 200]

plotGeneration	= True		# Make generation plot , True/False
activeLayer		= 3			# indexing starts from 0

lambda_start	= 350	# build a wavelength range to calculate over, starting wavelength (nm)
lambda_stop		= 800	# final wavelength (nm)
lambda_step		= 1		# wavelength step size

plotWavelengths	= [400 , 500 , 600]		# Wavelengths to plot |E|^2 for, in nm.

x_step=1.0 # grid spacing of device cross section simulation in nm 
# this is where the n/k optical constants for each material are stored.
# file format is : CSV (comma separated values)
# each materials file is prefixed with 'nk_' and the first 'matHeader' number of lines are skipped.
matDir			= 'matdata'	# materials data
matPrefix		= 'nk_'		# materials data prefix
matHeader		= 1				# number of lines in header

##
## START OF CODE. DO NOT EDIT FURTHER.
##

from os.path import join,isfile
from scipy.interpolate import interp1d
import numpy as np
import pylab as pl

# helper function
lambdas			= np.arange(lambda_start,lambda_stop+lambda_step,lambda_step,float)
t				= thicknesses

def openFile(fname):
	"""
	opens files and returns a list split at each new line
	"""
	fd = []
	if isfile(fname):
		fn = open(fname, 'r')
		fdtmp = fn.read()
		fdtmp = fdtmp.split('\n')
		# clean up line endings
		for f in fdtmp:
			f = f.strip('\n')
			f = f.strip('\r')
			fd.append(f)
		# make so doesn't return empty line at the end
		if len(fd[-1]) == 0:
			fd.pop(-1)
	else:
		print("%s Target is not a readable file" % fname)
	return fd

def get_ntotal(matName,lambdas):
	fname = join(matDir,'%s%s.csv' % (matPrefix,matName))
	fdata = openFile(fname)[matHeader:]
	# get data from the file
	lambList	= []
	nList		= []
	kList		= []
	for l in fdata:
		wl , n , k = l.split(',')
		wl , n , k = float(wl) , float(n) , float(k)
		lambList.append(wl)
		nList.append(n)
		kList.append(k)
	# make interpolation functions
	int_n	= interp1d(lambList,nList)
	int_k	= interp1d(lambList,kList)
	# interpolate data
	kintList	= int_k(lambdas)
	nintList	= int_n(lambdas)
	# make ntotal
	ntotal = []
	for i,n in enumerate(nintList):
		nt = complex(n,kintList[i])
		ntotal.append(nt)
	return ntotal

def I_mat(n1,n2):
	# transfer matrix at an interface
	r = (n1-n2)/(n1+n2)
	t = (2*n1)/(n1+n2)
	ret = np.array([[1,r],[r,1]],dtype=complex)
	ret = ret / t
	return ret

def L_mat(n,d,l):
	# propagation matrix
	# n = complex dielectric constant
	# d = thickness
	# l = wavelength
	xi = (2*np.pi*d*n)/l
	L = np.array( [ [ np.exp(complex(0,-1.0*xi)),0] , [0,np.exp(complex(0,xi))] ] )
	return L


# Constants
h 	= 	6.62606957e-34 	# Js Planck's constant
c	= 	2.99792458e8	# m/s speed of light
q	=	1.60217657e-19	#C electric charge

# PROGRAM

# initialize an array
n = np.zeros((len(layers),len(lambdas)),dtype=complex)

# load index of refraction for each material in the stack
for i,l in enumerate(layers):
	ni = np.array(get_ntotal(l,lambdas))
	n[i,:] = ni

# calculate incoherent power transmission through substrate

T_glass = abs((4.0*1.0*n[0,:])/((1+n[0,:])**2))
R_glass = abs((1-n[0,:])/(1+n[0,:]))**2

# calculate transfer marices, and field at each wavelength and position
t[0] 		= 0
t_cumsum	= np.cumsum(t)
x_pos		= np.arange((x_step/2.0),sum(t),x_step)
# get x_mat
comp1	= np.kron(np.ones( (len(t),1) ),x_pos)
comp2	= np.transpose(np.kron(np.ones( (len(x_pos),1) ),t_cumsum))
x_mat 	= sum(comp1>comp2,0) 	# might need to get changed to better match python indices

R		= lambdas*0.0
T		= lambdas*0.0
E		= np.zeros( (len(x_pos),len(lambdas)),dtype=complex )

# start looping
for ind,l in enumerate(lambdas):
	# calculate the transfer matrices for incoherent reflection/transmission at the first interface
	S = I_mat(n[0,ind],n[1,ind])
	for matind in np.arange(1,len(t)-1):
		mL = L_mat( n[matind,ind] , t[matind] , lambdas[ind] )
		mI = I_mat( n[matind,ind] , n[matind+1,ind])
		S  = np.asarray(np.mat(S)*np.mat(mL)*np.mat(mI))
	R[ind] = abs(S[1,0]/S[0,0])**2
	T[ind] = abs((2/(1+n[0,ind])))/np.sqrt(1-R_glass[ind]*R[ind])
	# good up to here
	# calculate all other transfer matrices
	for material in np.arange(1,len(t)):
		xi = 2*np.pi*n[material,ind]/lambdas[ind]
  		dj = t[material]
		x_indices	= np.nonzero(x_mat == material)
		x			= x_pos[x_indices]-t_cumsum[material-1]               
		# Calculate S_Prime
		S_prime		= I_mat(n[0,ind],n[1,ind])
		for matind in np.arange(2,material+1):
			mL = L_mat( n[matind-1,ind],t[matind-1],lambdas[ind] )
			mI = I_mat( n[matind-1,ind],n[matind,ind] )
			S_prime  = np.asarray( np.mat(S_prime)*np.mat(mL)*np.mat(mI) )
		# Calculate S_dprime (double prime)
		S_dprime	= np.eye(2)
		for matind in np.arange(material,len(t)-1):
			mI	= I_mat(n[matind,ind],n[matind+1,ind])
			mL	= L_mat(n[matind+1,ind],t[matind+1],lambdas[ind])
			S_dprime = np.asarray( np.mat(S_dprime) * np.mat(mI) * np.mat(mL) )
		# Normalized Electric Field Profile
		num = T[ind] * (S_dprime[0,0] * np.exp( complex(0,-1.0)*xi*(dj-x) ) + S_dprime[1,0]*np.exp(complex(0,1)*xi*(dj-x)))
		den = S_prime[0,0]*S_dprime[0,0]*np.exp(complex(0,-1.0)*xi*dj) + S_prime[0,1]*S_dprime[1,0]*np.exp(complex(0,1)*xi*dj)
		E[x_indices,ind] = num / den

# overall Reflection from device with incoherent reflections at first interface
Reflection = R_glass+T_glass**2*R/(1-R_glass*R)

#
# plot |E|^2 vs Position in device for wavelengths specified
#
fig1 = pl.figure(1)
fig1.clf()

ax1 = fig1.add_subplot(111)
ax1.set_ylim(ymin=0)
ax1.set_title('E-Field Intensity in Device')
ax1.set_ylabel('Normalized |E|$^2$Intensity')
ax1.set_xlabel('Position in Device (nm)')


# |E|^2 Plot
for ind,wvln in enumerate(plotWavelengths):
	xData = x_pos
	# get index
	nind  = np.where(lambdas==wvln)[0]
	yData = abs(E[:,nind])**2
	label = '%s nm' % wvln
	ax1.plot(xData,yData,label=label)
# Layer Bars
for matind in np.arange(2,len(t)+1):
	xpos = (t_cumsum[matind-2]+t_cumsum[matind-1])/2.0
	ax1.axvline(np.sum(t[0:matind]),linestyle=':')
	ax1.text(xpos,0.01,layers[matind-1],va='bottom',ha='center')
ax1.legend(loc='upper right')
fig1.show()

# Absorption coefficient in 1/cm
a = np.zeros( (len(t),len(lambdas)) )
for matind in np.arange(1,len(t)):
	a[matind,:] = ( 4 * np.pi * np.imag(n[matind,:]) ) / ( lambdas * 1.0e-7 )

#
# Plot normalized intensity absorbed / cm3-nm at each position and wavelength
# as well as the total reflection expected from the device
#
fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(111)
ax2.set_title('Fraction of Light Absorbed or Reflected')
ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('Light Intensity Fraction')
Absorption = np.zeros( (len(t),len(lambdas)) )
for matind in np.arange(1,len(t)):
	Pos 		= np.nonzero(x_mat == matind)
	AbsRate 	= np.tile( (a[matind,:] * np.real(n[matind,:])),(len(Pos),1)) * (abs(E[Pos,:])**2)
	Absorption[matind,:] = np.sum(AbsRate,1)*lambda_step*1.0e-7
	ax2.plot(lambdas,Absorption[matind,:],label=layers[matind])
ax2.plot(lambdas,Reflection,label='Reflection')
ax2.legend(loc='upper right',ncol=3)
fig2.show()

if plotGeneration:
	# load and interpolate AM1.5G Data
	am15_file = join(matDir,'AM15G.csv')
	am15_data = openFile(am15_file)[1:]
	am15_xData = []
	am15_yData = []
	for l in am15_data:
		x,y = l.split(',')
		x,y = float(x),float(y)
		am15_xData.append(x)
		am15_yData.append(y)
	am15_interp = interp1d(am15_xData,am15_yData,'linear')
	am15_int_y  = am15_interp(lambdas)
	#
	ActivePos = np.nonzero(x_mat == activeLayer)
	tmp1	= (a[activeLayer,:]*np.real(n[activeLayer,:])*am15_int_y)
	Q	 	= np.tile(tmp1,(np.size(ActivePos),1))*(abs(E[ActivePos,:])**2)
	# Exciton generatoion are
	Gxl		= (Q*1.0e-3)*np.tile( (lambdas*1.0e-9) , (np.size(ActivePos),1))/(h*c)
	if len(lambdas) == 1:
		lambda_step = 1
	else:
		lambda_step = (sorted(lambdas)[-1] - sorted(lambdas)[0])/(len(lambdas) - 1)
	Gx		= np.sum(Gxl,2)*lambda_step

	# plot
	fig3 = pl.figure(3)
	fig3.clf()
	ax3 = fig3.add_subplot(111)
	ax3.set_title('Generation Rate in Device')
	ax3.set_xlabel('Position in Device (nm)')
	ax3.set_ylabel('Generation Rate/sec-cm$^3$')
	ax3.set_xlim(0,t_cumsum[-1])
	ax3.plot(x_pos[ActivePos[0]] , Gx[0])
	# Layer Bars
	for matind in np.arange(2,len(t)+1):
		xpos = (t_cumsum[matind-2]+t_cumsum[matind-1])/2.0
		ax3.axvline(np.sum(t[0:matind]),linestyle=':')
		ax3.text(xpos,sorted(Gx[0])[0],layers[matind-1],va='bottom',ha='center')
	fig3.show()

	# calculate Jsc
	Jsc = np.sum(Gx)*x_step*1.0e-7*q*1.0e3

	# calculate parasitic absorption
	parasitic_abs = (1.0 - Reflection - Absorption[activeLayer,:])

	print Jsc

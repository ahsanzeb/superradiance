
import numpy as np
import matplotlib
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.axes as axx
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.ticker import MaxNLocator,FixedLocator,IndexLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import StrMethodFormatter, NullFormatter
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib.ticker import ScalarFormatter
from matplotlib import colors, ticker, cm
import matplotlib.colors as mcol
import matplotlib.patches as patches
from qutip import *


#plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rc('text', usetex=True)


wv=0.2;

color = ["b","g","r","k","m","y","o"]
lss= ['--','--','-','--','-'];
dashes = [[2,2],[10,4],[1000,1],[10,2]]

lwl = 1.5 ; # line width , loop
ms=1.5;


# column title label size:
ttlfs= 18;
# tics label size:
lblsz = 16;
# x label font size
xlfs = 20;
# y label font size
ylfs = 20;


#****************************************
#	fonts of text and legends
lgf=11
lgf2=14

#****************************************
#	make colors list
#****************************************
#start_time = 0
#end_time = 7
# Generate some dummy data.
#tim = range(start_time,end_time+1)
ls1='--k'
ls2='--k'

tim=range(5);

# Make a user-defined colormap.
# cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["b","r"])
#cm1 = cm.get_cmap('Blues', 10) 
#cm1 = cm.get_cmap('coolwarm', 7)
b2r = ['#4444FF','#AA0000']
b2r = ['#FFBABA','#FF0000']

b2r = ['#fdcc8a','#e34a33']
b2r = ['#FFDDDD','#FF0000']
b2r = ['red','yellow']

cm1 = mcol.LinearSegmentedColormap.from_list('mycol',b2r)
#cm1 = cm.get_cmap('cool', 7)

# Make a normalizer that will map the time values from
# [start_time,end_time+1] -> [0,1].
cnorm = mcol.Normalize(vmin=min(tim),vmax=max(tim))
# Turn these into an object that can be used to map time values to colors and
# can be passed to plt.colorbar().
cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
cpick.set_array([])
# use as:
#for y, t in zip(ydat,tim):
#    A.plot(xdat,y,color=cpick.to_rgba(t))
#
#plt.colorbar(cpick,label="Time (seconds)")





#****************************************
# ----------------------------------------------------------------------
nxmax= 200
nymax= 150
xi=-2; xf=4
yi=-2.5; yf=2.5;

xvec = np.linspace(xi,xf, nxmax)
yvec = np.linspace(yi,yf, nymax)
dp=yvec[1]-yvec[0]; # dp
ndel=3; # temporarily use delta for nev=3 eigenstates, modify if a set of delta vals
nds=ndel;
M=31;
wrs = np.array([0.5,1.0,1.5,2.0]);
nwrs=len(wrs);

ds = np.array([0]); #np.linspace(-2,2, ndel)
ms =np.array([1]); #[1,2,5,10,100]
rhoex=ms/1.0;
nms = len(ms);

print('wrs = ', wrs)
print('ms = ', ms)

# ----------------------------------------------------------------------
def checkfile(filename):
	import os
	if os.path.exists(filename):
		print("The files "+filename+" exist")
		print("Overwrite? y/n")
		overwrite = input('===> ')
		if(overwrite=='y'):
			os.remove(filename)
		else:
			print("Stopping.... ")
			exit()
	else:
		print("Creating output files ... ")
	return
# ----------------------------------------------------------------------

# check if old files are present and you dont want to remove them...
checkfile("parameters-f.dat")



# ----------------------------------------------------------------------
files=["n-30/dmfield.dat"]
# ----------------------------------------------------------------------
print('reading data files')
filename=files[0];
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split() ] )
a=np.array(a);
dmdn = a.reshape((nms,nwrs, ndel, M, M))
# ----------------------------------------------------------------------
Pxdn = np.zeros((nms,nwrs, nds,nxmax))
for im in range(nms):
	print('rho to wigner function: working on m = ',ms[im])
	for iwr in range(nwrs):
		for idel in range(ndel):
			rho=Qobj(dmdn[im,iwr,idel,:,:])
			W = wigner(rho, xvec, yvec)
			Pxdn[im,iwr,idel,:] = np.sum(W,axis=0) * dp # integrate over p
# ----------------------------------------------------------------------


def FitSHO(i,j,x,y):
	from scipy.optimize import curve_fit
	from scipy import asarray as ar,exp
	def gaus(x,a,x0,sigma):
		return a*exp(-(x-x0)**2/(2*sigma**2))

	# For normalised gaussian
	#n = len(x)                          #the number of data
	#mean = sum(x*y)/n                   #note this correction
	#sigma = sum(y*(x-mean)**2)/n        #note this correction

	# for unnormalised gaussian
	mean = sum(x * y) / sum(y)
	sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))

	popt,pcov = curve_fit(gaus,x,y,p0=[max(y),mean,sigma])

	#print('popt = ',popt)

	grid[i,j].plot(x,y,'b-',linewidth=0.5,label='data')
	grid[i,j].plot(x,gaus(x,*popt),'r:',linewidth=0.5,label='fit')

	return popt

# ----------------------------------------------------------------------


def FitTwoSHOs(i,j,x,y):
	from scipy.optimize import curve_fit
	from scipy import asarray as ar,exp
	def gaus(x,a,x0,sigma):
		return a*exp(-(x-x0)**2/(2*sigma**2))

	def biexp(x,a1,a2,x01,x02,sigma):
		z = a1*exp(-(x-x01)**2/(2*sigma**2));
		z += a2*exp(-(x-x02)**2/(2*sigma**2))
		return z


	# for unnormalised gaussian
	mean = sum(x * y) / sum(y)
	sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))

	p0 = [max(y),max(y),0,mean,sigma];
	popt,pcov = curve_fit(biexp,x,y,p0=p0)

	#print('popt = ',popt)

	grid[i,j].plot(x,y,'b-',linewidth=0.5,label='data')
	grid[i,j].plot(x,biexp(x,*popt),'r:',linewidth=0.5,label='fit')

	return

# ----------------------------------------------------------------------
def OneSHO():
	sqrt2 = 1/np.sqrt(2);
	par = np.zeros((nms,nwrs, nds, 2, 4))
	fdn=open('parameters-f.dat','ab')
	for im in range(nms):
		print('working on m = ',ms[im])
		for iwr in range(nwrs):
			for idel in range(ndel):
				par[im,iwr,idel,0,0:3] = FitSHO(iwr,0,xvec,Pxdn[im,iwr,idel,:])	
				# displacements normalise lam+-
				par[im,iwr,idel,0,1] = par[im,iwr,idel,0,1]*sqrt2
				# frequencies nu+- from sigma+-
				par[im,iwr,idel,0,3] = np.reciprocal(par[im,iwr,idel,0,2]**2)/2
				x = par[im,iwr,idel,0,:].reshape((1,4)) # dn
				np.savetxt(fdn,x, fmt='%15.10f', delimiter=',')
	fdn.close()
	return par
# ----------------------------------------------------------------------



def plot(im):
	for iwr in range(nwrs):
		#  par[iwr,id,0:3] = [a, mean, sigma]
		grid2[0].plot(ds,par[im,iwr,:,0,1]	,'-',linewidth=3,color=cpick.to_rgba(iwr))
		grid2[0].plot(ds,par[im,iwr,:,1,1],':',linewidth=3,color=cpick.to_rgba(iwr))
		grid2[1].plot(ds,par[im,iwr,:,0,3]	,'-',linewidth=3,color=cpick.to_rgba(iwr))
		grid2[1].plot(ds,par[im,iwr,:,1,3],':',linewidth=3,color=cpick.to_rgba(iwr))

	grid2[0].yaxis.set_minor_formatter(NullFormatter())
	grid2[0].yaxis.set_major_formatter(ScalarFormatter())
	grid2[1].yaxis.set_minor_formatter(NullFormatter())
	grid2[1].yaxis.set_major_formatter(ScalarFormatter())

	#grid2[1].set_ylim([0.82,1.01])
	#grid2[1].set_yticks([0.85,0.9,0.95,1.0])
	#grid2[1].set_yticklabels(['$0.80$','$0.85$','$0.90$','$0.95$','$1.00$'])

	grid2[0].set_ylabel(r'$\lambda_{\pm}/\lambda$')
	grid2[1].set_ylabel(r'$\nu_{\pm}/\nu$')
	return
# ----------------------------------------------------------------------


cm2 = mcol.LinearSegmentedColormap.from_list('mycol',['red','blue'])
cnorm2 = mcol.Normalize(vmin=0,vmax=nwrs-1)
cpick2 = cm.ScalarMappable(norm=cnorm2,cmap=cm2)
cpick2.set_array([])
lwl = 0.5 ; # line width

fig, grid = plt.subplots(nrows=nwrs, ncols=2)

# call the function OneSHO() to calc all parameters:
par = OneSHO();

plt.savefig('htc-Px-fit-sho.pdf', format='pdf', dpi=100)
#plt.show()


plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close(fig) # Close a figure window

if 1:
	# new figure; fitting parameters (effective lambda, wv)
	fig2, grid2 = plt.subplots(nrows=1, ncols=2)
	sqr2 = 1.0/np.sqrt(2); # to get lam_eff from x_0
	par = par*sqr2; # to get wv_eff factor also affected
	
	#grid2[0].set_title(r'$\lambda^\prime_{\pm}=\mp3$')
	#grid2[1].set_title(r'$\nu=0.2$')

	grid2[0].set_yscale('log')
	#grid2[1].set_yscale('log')

	for im in range(nms):
		plot(im)








plt.tight_layout()

plt.savefig('effective-sho-parameters.pdf', format='pdf', dpi=100)

#plt.show()


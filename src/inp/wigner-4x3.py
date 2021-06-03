
# ---------------------------------------------
#  Usage: 
#  Steps:
# 1. First run wigner.py to calculate wigner functions from density matrices
# 2. Then, desired file names are to be placed in the list 'files' below in order to be appeared in 2x3 panel from top left to right and repeat from left to right for next lower line.
# 3. Then this script should be run.
# ---------------------------------------------
# Note: xlim and ylim in this script should be the same as used in wigner.py to calculate wigner functions.
# ---------------------------------------------


#  The cbar is also of different size


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.axes as axx
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.ticker import MaxNLocator,FixedLocator,IndexLocator
#matplotlib.use('MacOSX')

import matplotlib.patches as mpatch




# font size of labels: lambda_0, N
flam = 18;
fn = 18;
# font size of tics labels, x,y axes labels, cmap
ftics = 14;
fxl = 18;
fyl= fxl;
fcmp = 14
# title size: conditions, latex
fts = 18


files=["rho0_vib-N=2-M=50-wr=1.2.dat-lambda-2.0-plot-wf-0.pdf-w.txt","rho1_vib-N=2-M=50-wr=1.2.dat-lambda-2.0-plot-wf-1.pdf-w.txt","rho2_vib-N=2-M=50-wr=1.2.dat-lambda-2.0-plot-wf-2.pdf-w.txt","rho0_vib-N=2-M=50-wr=1.2.dat-lambda-2.5-plot-wf-0.pdf-w.txt","rho1_vib-N=2-M=50-wr=1.2.dat-lambda-2.5-plot-wf-1.pdf-w.txt","rho2_vib-N=2-M=50-wr=1.2.dat-lambda-2.5-plot-wf-2.pdf-w.txt"];

# the values used for the wigner
#xi=-4.5; xf=2
#yi=-2.2; yf=2.2;
# for fliplr(W) plot, here we dont give xvec,yvec so that we will just used reversed xi,xf etc.
xi=-2; xf=5
#yi=-2.2; yf=2.2;
yi=-1.75; yf=1.75;

#xtcs=range(1,-5,-1);
xtcs=[-1,0,1,2,3,4];#range(-2,5,1);
ytcs=range(2,-3,-1);
#ytcs = [1.5,1,0.5,0,-0.5,-1,-1.5]

# row 1, col1
#read first data set from file:
filename=files[0]; 
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
cdat=data[:,:];

# row 1, col2
# Loop over all data files:
filename=files[1];
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
# add to cdat 
cdat=np.dstack((cdat, data))

# row 1, col3
filename=files[2]; 
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
# add to cdat 
cdat=np.dstack((cdat, data))



# row 2, col1
# read all other data sets from files:
filename=files[3]; 
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
# add to cdat 
cdat=np.dstack((cdat, data))

# row 2, col2
filename=files[4];  
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
# add to cdat 
cdat=np.dstack((cdat, data))

# row 2, col3
filename=files[5]; 
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
# add to cdat 
cdat=np.dstack((cdat, data))








files=["rho0_vib-N=20-M=8-wr=1.0.dat-lambda-2.0-plot-wf-0.pdf-w.txt",
"rho1_vib-N=20-M=8-wr=1.0.dat-lambda-2.0-plot-wf-1.pdf-w.txt",
"rho2_vib-N=20-M=8-wr=1.0.dat-lambda-2.0-plot-wf-2.pdf-w.txt",
"rho0_vib-N=20-M=8-wr=1.0.dat-lambda-2.5-plot-wf-0.pdf-w.txt",
"rho1_vib-N=20-M=8-wr=1.0.dat-lambda-2.5-plot-wf-1.pdf-w.txt",
"rho2_vib-N=20-M=8-wr=1.0.dat-lambda-2.5-plot-wf-2.pdf-w.txt"];

# row 1, col1
filename=files[0];
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
# add to cdat 
cdat=np.dstack((cdat, data))

# row 1, col2
# Loop over all data files:
filename=files[1];
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
# add to cdat 
cdat=np.dstack((cdat, data))

# row 1, col3
filename=files[2]; 
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
# add to cdat 
cdat=np.dstack((cdat, data))



# row 2, col1
# read all other data sets from files:
filename=files[3]; 
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
# add to cdat 
cdat=np.dstack((cdat, data))

# row 2, col2
filename=files[4];  
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
# add to cdat 
cdat=np.dstack((cdat, data))

# row 2, col3
filename=files[5]; 
a=[];
fin = open(filename,'r')
for line in fin.readlines():
	a.append( [ float (x) for x in line.split(',') ] )
data=np.array(a[:][:])
# add to cdat 
cdat=np.dstack((cdat, data))





# min and max values for cbar
v1 = np.amin(cdat)
v2 = np.amax(cdat)
print(v1, v2)
print(np.amin(data), np.amax(data))

# Generate some data that where each slice has a different range
# (The overall range is from 0 to 2)
# data = np.random.random((6,10,10))
# data *= np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])[:,None,None]

    
fig = plt.figure(1, figsize=(7.45, 8), dpi=80) # 

grid = AxesGrid(fig, 111,  # similar to subplot(111)
                    nrows_ncols=(4, 3),
                    axes_pad=0.0,
                    share_all=True,
                    label_mode="L",
			  cbar_mode="single",
			 cbar_location="right"
)

print("cdat.shape",cdat.shape)
              
cmap = cm.get_cmap('Blues', 10)    #'PiYG', 11 discrete colors

#i=-1
# for dat in cdat:
for i in range(12):
	#i=i+1
	print("i = ",i)
	im = grid[i].imshow(np.fliplr(cdat[:,:,i]), vmin=v1, vmax=v2,origin="lower", extent=[xi,xf,yi,yf],aspect=(xf-xi)/(yf-yi),cmap=cmap)
	ax = im.get_axes( )
	ax.set_xlabel(r'$x$', fontsize=fxl, fontweight='bold')
#		ax.set_ylabel(r'$p$', fontsize=fyl, fontweight='bold')
#	ax.xticks(xvec, fontsize=9)
#	ax.yticks(yvec, fontsize=9)
	ax.xaxis.set_major_locator(FixedLocator(xtcs))
	ax.set_xticklabels(xtcs,fontsize=ftics)
	ax.yaxis.set_major_locator(FixedLocator(ytcs))
	ax.set_yticklabels(ytcs,fontsize=ftics)

	if i>2 and i<6:
		ax.axhline(y=yi,linewidth=6, color='k')
		
	if i==0:
		ax.set_title(r'$\left | P_{}^{}\right\rangle$', fontsize=fts, y=1.03)
		ax.set_ylabel('$\lambda_0 = 2.0$\n$p$', fontsize=flam, fontweight='bold')
	if i==1:
		ax.set_title(r'$\left|X_{}\right\rangle_i$', fontsize=fts, y=1.03)
	if i==2:
		ax.set_title(r'$\left|X_{}\right\rangle_{j \ne i}$', fontsize=fts,y=1.03)
	if i==3:
	 ax.set_ylabel('$\lambda_0 = 2.5$ \n $p$', fontsize=flam, fontweight='bold')

	if i==0+6:
	 ax.set_ylabel('$\lambda_0 = 2.0$ \n $p$', fontsize=flam, fontweight='bold')
	if i==3+6:
	 ax.set_ylabel('$\lambda_0 = 2.5$ \n $p$', fontsize=flam, fontweight='bold')



# for paper, no title
#fig.suptitle(r'$W(x,p): N=2,\,\Omega=0.2,\, \omega_R = 1.25$', fontsize=24, fontweight='bold',color='blue')


##cax = fig.add_axes([0.9, 0.1, 0.03, 0.5])
##fig.colorbar(im, cax=cax)
#ax.Axes.set_xlabel("x")
#ax.Axes.set_ylabel("y")
cbar=grid.cbar_axes[0].colorbar(im,format='%.2f')
cbar.ax.tick_params(labelsize=fcmp) 

#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
fig.text(0.0175, 0.72, '$N=2$',rotation='vertical', fontsize=fn, fontweight='bold')#, bbox=props)
fig.text(0.0175, 0.32, '$N=20$',rotation='vertical', fontsize=fn, fontweight='bold')#, bbox=props)

#plt.subplots_adjust(left=0.1, right=0.97, top=0.9, bottom=0.1)

#plt.tight_layout()

plt.savefig('fig-1.pdf', format='pdf', dpi=1000)

plt.show()

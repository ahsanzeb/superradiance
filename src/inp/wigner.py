
 
#  smaller xvec for larger N

# WARNING: RUN 'python3.4 wigner-function.py' 
# because 'python qutip' does not run quitp
import numpy as np
import glob
import os
from qutip import *
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
# import matplotlib.animation as anim
# import matplotlib.image as mpimg

import ntpath


nxmax= 500
nymax= 300
xi=-7; xf=2
yi=-4; yf=4;

# save W(x,p) at these lambda (if they exist in lambin, obviously)
lambsvW=[0.5,1.0,	1.5, 2.0]

xvec = np.linspace(xi,xf, nxmax)
yvec = np.linspace(yi,yf, nymax)
dp=yvec[1]-yvec[0]; # dp

#path='/Users/panda/Desktop/calc/M-sites/dm/nb-rho-data'
path=os.getcwd();

	
probx0=[];probx1=[];probx2=[];
probn0=[];probn1=[];probn2=[];
prob0=[]; prob1=[]; prob2=[]; 
prob=[]; 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Loops over files starts here
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
for filename in glob.glob(os.path.join(path, '*.dat')):
	fin = open(filename,'r')
	param=[];lambin=[];
	dm=[]; dmall=[];
	iline=0; irow=0;
	for line in fin.readlines():
		iline=iline+1;
		if iline == 1 : # param
			param.append( [ float (x) for x in line.split(',') ] )
			m=int(param[0][1]);
			#print("param = ",param[0][:])
			#print("m  = ",m)
		elif iline==2: # lambda values array
			lambin.append( [ float (x) for x in line.split(',') ] )
			lambin=np.array(lambin[0][:])
			# print("lambin = ",lambin)
			# print("lambin.shape = ",lambin.shape)
			zzz = lambin.shape;
			nlmax = zzz[0]
			
			# print("nlmax = ",nlmax)
		else : # density matrix 
			if irow <= m : # dm at a given lambda
				irow=irow+1;
				dm.append( [ float (x) for x in line.split(',') ] )
				#print("irow, m = ", irow, m)
				if irow > m: # dmall=dm_all, at all lambdas
					#print("----------")
					#print("dm  = ", dm)
					irow =0;
					dmall.append(dm);
					dm=[];
		



	il=-1
	for dm in dmall:
		il=il+1;
		dm=np.array(dm);
		# normalise
		tra=np.trace(dm);
		dm0=dm[:][:];# dm0 not normalised
		dm = dm0/tra;# dm normalised
	
		rho=Qobj(dm)
		W = wigner(rho, xvec, yvec)
	
		# calc prob
		px=[]; pn=[];
		px=np.sum(W,axis=0) # P(x)= integrate W(x,y) over y
		px = px*dp # normalise
		pn=np.diagonal(dm) # P(n)= diag elem of dm
	
		pp=0;
		pp0=np.diagonal(dm0) # diag elem of dm0
		pp=np.sum(pp0) # trace(dm)


		fig=plt.figure();
		plt2 = plt.contourf(xvec, yvec, W, 100)
		plt.colorbar(plt2)
		plt.title(" ")
		plt.xlabel(r'$x$', fontsize=20)
		plt.ylabel(r'$p$', fontsize=20)
		plt.tight_layout() #plt.show()
		ttl0=r"${N="+str(int(param[0][0]))+", M="+str(int(param[0][1]))+', \Omega='+str(param[0][2])+', \omega_R='+str('%.1f'%(param[0][3]))+', \lambda_0='+str('%.2f'%(lambin[il]))+", P="+str('%.3f'%(pp))+"}$"

		#ttl0=r"$N=$"+str(int(param[0][0]))+", M="+str(int(param[0][1]))+r', $\Omega=$'+str(param[0][2])+r', $\omega_R=$'+str('%.1f'%(param[0][3]))+r', $\lambda_0=$'+str('%.2f'%(lambin[il]))+", P="+str('%.3f'%(pp))

		filename0=ntpath.basename(filename)
		#print(filename0)
		if filename0.startswith('rho0'):
			#print('rho0 read correctly')
			fout = filename+'-lambda-'+str(lambin[il])+'-plot-wf'+'-0.pdf'
			ttlx=r'$\mid 1_c0_x \rangle :$'
			probx0.append(px)
			probn0.append(pn)
			prob0.append(pp)
		elif filename0.startswith('rho1'):
			fout = filename+'-lambda-'+str(lambin[il])+'-plot-wf'+'-1.pdf'
			ttlx=r'$|0_c1_x \rangle: $'
			probx1.append(px)
			probn1.append(pn)
			prob1.append(pp)
		elif filename0.startswith('rho2'):
			fout = filename+'-lambda-'+str(lambin[il])+'-plot-wf'+'-2.pdf'
			ttlx=r'$|0_c0_x \rangle :$'
			probx2.append(px)
			probn2.append(pn)
			prob2.append(pp)
		else: 
			print("error?: data files other than rho0,1,2... ")
			exit()

		ttl=ttlx+ttl0
		plt.title(ttl, fontsize=14, color='blue')
		plt.subplots_adjust(top=0.9) # use a lower number to make more vertical space

		plt.savefig(fout)
		plt.close(fig)
	# np.savetxt(fout+"-ww.txt", W, fmt='%10.6f',delimiter=',')

	# save W(x,p) for desired values of lambda
		ll = lambin[il]
		print("lamda = ", ll)
		for givenl in lambsvW:
			if abs(ll-givenl) < 0.001:
				np.savetxt(fout+"-w.txt", W, fmt='%10.6f',delimiter=',')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Loops over files completes here
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# save all data
probab=np.column_stack((lambin, prob0, prob1, prob2))
np.savetxt("prob.txt",probab, fmt='%10.6f', delimiter=',')

np.savetxt("probx-0.txt",probx0, fmt='%10.6f', delimiter=',') 
np.savetxt("probn-0.txt",probn0, fmt='%10.6f', delimiter=',')

np.savetxt("probx-1.txt",probx1, fmt='%10.6f', delimiter=',') 
np.savetxt("probn-1.txt",probn1, fmt='%10.6f', delimiter=',')

np.savetxt("probx-2.txt",probx2, fmt='%10.6f', delimiter=',') 
np.savetxt("probn-2.txt",probn2, fmt='%10.6f', delimiter=',')

 # plot the above saved data


fig=plt.figure()
 # plt.plot(lamb, prob0, marker='o', linestyle='--', color='r', label='P(1c)')
plt.plot(lambin, prob0, color='r', label=r'$P(1_c0_x)$')
plt.plot(lambin, prob1, color='b', label=r'$P(0_c1_x)$')
plt.plot(lambin, prob2,  color='g', label=r'$P(0_c0_x)$')

plt.xlabel(r'$\lambda_0$', fontsize=20)
plt.ylabel('Conditional Probabilities', fontsize=14)
plt.title('')
plt.legend()
 # plt.show()

ttl0=r"${N="+str(int(param[0][0]))+", M="+str(int(param[0][1]))+', \Omega='+str(param[0][2])+', \omega_R='+str('%.1f'%(param[0][3]))+"}$"


 #ttl0="N="+str(int(param[0][0]))+", M="+str(int(param[0][1]))+r', $\Omega$'+str(param[0][2])+r', $\omega_R$'+str('%.1f'%(param[0][3]))
plt.title(ttl0, fontsize=14, color='blue')
plt.subplots_adjust(top=0.9)
fout="prob-vs-lambda.pdf"
plt.savefig(fout)
plt.close(fig)




 # m = np.loadtxt("probx-0.txt", delimiter=" ")

for i in range(3):
	if i==0:
		ttlx = r'$P(x|1_c0_x): $'
		fig=plt.figure()
		probx0=np.array(probx0)
		plt2 = plt.contourf(lambin, xvec, probx0.transpose(), 100)
		fout="prob0x.pdf"
	if i==1:
		ttlx = r'$P(x|0_c1_x): $'
		fig=plt.figure()
		probx1=np.array(probx1)
		plt2 = plt.contourf(lambin, xvec, probx1.transpose(), 100)
		fout="prob1x.pdf"
	if i==2:
		ttlx = r'$P(x|0_c0_x): $'
		fig=plt.figure()
		probx2=np.array(probx2)
		plt2 = plt.contourf(lambin, xvec, probx2.transpose(), 100)
		fout="prob2x.pdf"

	plt.colorbar(plt2)
	plt.title(" ")
	plt.xlabel(r'$\lambda_0$', fontsize=20)
	plt.ylabel(r'$x$', fontsize=20)
	plt.tight_layout() #plt.show()

	ttl0=r"${N="+str(int(param[0][0]))+", M="+str(int(param[0][1]))+', \Omega='+str(param[0][2])+', \omega_R='+str('%.1f'%(param[0][3]))+"}$"

	#ttl0="N="+str(int(param[0][0]))+", M="+str(int(param[0][1]))+r', $\Omega$'+str(param[0][2])+r', $\omega_R$'+str('%.1f'%(param[0][3]))

	ttl=ttlx+ttl0
	plt.title(ttl, fontsize=14, color='blue')
	plt.subplots_adjust(top=0.9)

	plt.savefig(fout)
	plt.close(fig)




ndim=len(dm);
nvec=range(ndim);

for i in range(3):

	if i==0:
		ttlx = r'$P(n|1_c0_x): $'
		fig=plt.figure()
		probn0=np.array(probn0)
		plt2 = plt.contourf(lambin, nvec, probn0.transpose(), 100)
		fout="prob0n.pdf"
	if i==1:
		ttlx = r'$P(n|0_c1_x): $'
		fig=plt.figure()
		probn1=np.array(probn1)
		plt2 = plt.contourf(lambin, nvec, probn1.transpose(), 100)
		fout="prob1n.pdf"
	if i==2:
		ttlx = r'$P(n|0_c0_x): $'
		fig=plt.figure()
		probn2=np.array(probn2)
		plt2 = plt.contourf(lambin, nvec, probn2.transpose(), 100)
		fout="prob2n.pdf"

	plt.colorbar(plt2)
	plt.title(" ")
	plt.xlabel(r'$\lambda_0$', fontsize=20)
	plt.ylabel(r'$n$', fontsize=20)
	plt.tight_layout() #plt.show()
	ttl0=r"${N="+str(int(param[0][0]))+", M="+str(int(param[0][1]))+', \Omega='+str(param[0][2])+', \omega_R='+str('%.1f'%(param[0][3]))+"}$"


	#ttl0="N="+str(int(param[0][0]))+", M="+str(int(param[0][1]))+r', $\Omega$'+str(param[0][2])+r', $\omega_R$'+str('%.1f'%(param[0][3]))

	ttl=ttlx+ttl0
	plt.title(ttl, fontsize=14, color='blue')
	plt.subplots_adjust(top=0.9)

	plt.savefig(fout)
	plt.close(fig)









	module correlation
	use modmain, only: Hf, eig, kappa2,prntstep,fft,
     .							 task, eig0, basis, nsym, nph, ntotb, ntot
	implicit none

	public :: tcorr
	private :: rk4step, yprime
	private :: FourierT, writecorrt, writecorrw
		
	double precision :: dt, dth, dt6, twopi,kappa2r
	double complex :: iotam, phase
	double complex, dimension(:), allocatable:: corrt,corrw
	double precision, dimension(:), allocatable:: ws, ts,wss
	integer:: nt, nw
	character :: fname*100

	double complex,dimension(:), allocatable :: psi0, psit
	double precision, dimension(:), allocatable ::psitnorm
	integer :: ij
	
	contains
! Ahsan Zeb, 04 Nov, 2019, modified from poalriton code routine tcorr.f
! This subroutine takes the hamiltonian, initial state and other parameters, evolve the state and compute the two time correlation that is given as output.

	subroutine tcorr(dt0,w1,w2,nt0,nw0,ij0, task) 
	implicit none
	double precision, intent(in) :: dt0,w1,w2
	integer, intent(in) :: nt0,nw0, ij0
	integer, intent(in):: task
	! local
	integer :: i, nsize
	double precision:: dw
	double complex :: e0dt
	
	dt = dt0;	
	nt = nt0;
	if(fft) then
		nw = nt0; !FFTW: nw0;
	else
		nw = nw0;
	endif
	ij = ij0;

	! set variables for RKK4 integration
	dth = dt/2.0d0;
	dt6 = dt/6.0d0;
	iotam = (0.0d0,-1.0d0);

	kappa2r = kappa2 !/m;
	
	allocate(corrt(nt))
	allocate(psitnorm(nt))
	allocate(ws(nw))
	!allocate(wss(nw))
	allocate(ts(nt))
	
	! calc ws: frequencies....
	ws = 0.0d0;
	twopi= 2.0d0*3.14159265359d0;
	if(fft) then
		dw = twopi/(nt*dt); !(w2-w1)*1.0d0/nw;
	else
		dw = (w2-w1)*1.0d0/(nw-1);
	endif

	do i=1,nw
		!ws(i) = w1 + (i-1)*dw 
		ws(i) =  w1 + (i-1)*dw 
	end do
	! calc ts: times
	ts=0.0d0;
	do i=1,nt
		ts(i) = (i-1)*dt 
	end do

	! ntot given in input will be either global ntot or global ntotg
	! depending on the task.
	allocate(psi0(ntot))
	allocate(psit(ntot))

	!write(6,*)'tcorr: ntot, Hf%ntot = ',ntot, Hf%ntot
	!write(6,*)'Hf%Col(:) = ', Hf%Col(:)
	
	! construct init state for time evolution
	if(task == 101) then 	!absorption...
		call makepsi0ad(ij0)
	else
		stop "correlation: task =101 only"
	endif

	write(6,*)'norm of psi0 = ', sum(abs(psi0)**2)
	!write(6,*)'kappa2 = ',kappa2r
	
	! set psi0 and corr(t=0)
	psitnorm(1) = real(DOT_PRODUCT(psi0, psi0));
	psi0 = psi0/dsqrt(psitnorm(1));
	psitnorm(1)=1.0d0;
	psit(:) = psi0(:);
	corrt(1) = dcmplx(1.0, 0.0); !DOT_PRODUCT(psi0, psit);  

	e0dt = (-iotam*eig0%eval(1)-kappa2r)*dt;
	
	write(6,'(a,f15.10)')'ELP= ',eig(ij0)%eval(1);
	do i=2,nt
		call rk4step(psit)
		!phase = zexp((-iotam*eig(ij0)%eval(1)-kappa2r)*(i-1)*dt);
		phase = zexp((i-1)*e0dt);
		corrt(i) = DOT_PRODUCT(psi0, psit) * phase

		psitnorm(i) = real(DOT_PRODUCT(psit, psit))

		if (mod(i,prntstep) == 0) then
			write(6,'(a20,2x,i10,a10,i10)') ' time evolution step ',i,
     .  ' out of ',nt
			write(6,'(a,E15.10)')'|psit|/|psi| =', psitnorm(i)/psitnorm(1)
		endif

		if(psitnorm(i) > 10.0*psitnorm(1)) then
			write(*,*)'corr: RKK4 diverging... stopping..., decrease dt'
			write(*,*)'corr: psitnorm = ', psitnorm(i)
			stop
		endif
	end do

	!call writecorrt(task)
	call writeatnode(ij0,nt,ts, corrt,'temp-t')
	if(fft) then
		call FFourierT()
	else
		call FourierT()
	endif
	!call writecorrw(task)
	call writeatnode(ij0,nw,ws, corrw,'temp-w')
	
	deallocate(psi0,psit,corrt,corrw,psitnorm,ws,ts)

	write(6,'(a)')'tcorr: finished...!'
	return
	end subroutine tcorr
	!----------------------------------------
	subroutine yprime(vin,vout)
	implicit none
	double complex, intent(in) :: vin(ntot)
	double complex, intent(out) :: vout(ntot)
	integer :: i,j	, k

	! Hf upper or lower triangular only; diag elements halved.
	vout = 0.0d0
	do i=1,Hf%ntot !nrp-1
		do j=Hf%RowPntr(i),Hf%RowPntr(i+1)-1
			k = Hf%Col(j); 
			vout(i) = vout(i) + vin(k)*Hf%dat(j);
			vout(k) = vout(k) +  vin(i)*Hf%dat(j);
		end do
	end do
	vout = vout*iotam - kappa2r*vin
	return
	end subroutine yprime
	!----------------------------------------
	subroutine rk4step(y)
	implicit none
	double complex, intent(inout) :: y(ntot)
	! local
	double complex, dimension(ntot) :: k1,k2,k3
      call yprime(y,k1)
      call yprime(y+k1*dth,k2)
      k1 = k1 + 2*k2
      call yprime(y+k2*dth,k3)
      k1 = k1 + 2*k3
      call yprime(y+k3*dt,k2)
      k1 = k1 + k2
      y = y + k1*dt6
	return
	end subroutine rk4step
	!----------------------------------------

	! fix sqrt(n_photon) terms in makepsi0ad(), makepsi0a()
	
	! add a photon to the given state to make psi0 for evolution
	! fully inside the sym space... no projection on non-sym space... 
	subroutine makepsi0ad(ijob)
	implicit none
	integer, intent(in) :: ijob
	integer :: i1,i2,p, i3
	double precision :: fac

	psi0 = 0.0d0
	! psi0 is complex, set real part.
	! add/create a photon in init state
	i1 = 0;
	do p=0,nph-1; 
	 fac = dsqrt(dble(p+1)); ! sqrt(n_photon + 1) init state
	 i2 = i1 + ntotb;
	 i3 = i2 + ntotb; ! shift by ntotb, i.e., one block
	 psi0(i2+1:i3) = fac * eig0%evec(i1+1:i2,1); ! psi0 is complex, set real part.
	 i1 = i2;
	enddo

	!dcmplx(eig(ijob)%evec(:,1), 0.0d0); 
	return
	end 	subroutine makepsi0ad
	!----------------------------------------

	!----------------------------------------
	subroutine FourierT()
	!use fftpack5
	implicit none
	integer :: i

	allocate(corrw(nw))
	do i=1,nw
		corrw(i) = sum(zexp(-iotam*ts*ws(i))*corrt)
	end do
	corrw = corrw/dsqrt(twopi*nw*1.0d0);
	
	return
	end 	subroutine FourierT
	!----------------------------------------

	subroutine FFourierT()
	!use fftpack5
	implicit none
	integer :: i
	
	call cfftnd(1,nt,+1,corrt)

	! nw = nt ? cant we have more ws?
	allocate(corrw(nt))
	do i=1,nw/2 ! nw power of 2
		corrw(i) = corrt(nw/2+i)
		corrw(nw/2 + i) = corrt(i)
		ws = ws - ws(nw/2)
	end do
	corrw = corrw/dsqrt(twopi*nw*1.0d0);
	
	return
	end 	subroutine FFourierT
	!----------------------------------------
	subroutine writecorrt(task)
	implicit none
	integer, intent(in) :: task
	integer :: i
	open(10,file='correlation-vs-t.dat',
     .             action='write',position='append')

	if(task==101) then
		write(10,'(a)')'# Photon'
	elseif(task==105) then
		write(10,'(a)')'# Exciton |up><up|'
	elseif(task==106) then
		write(10,'(a)')'# Exciton |dn><dn|'
	endif
	do i=1,nt
		write(10,'(f20.10,3x,2E20.10,2x,f15.10)') ts(i), 
     .   real(corrt(i)), aimag(corrt(i)), psitnorm(i)
	end do
	write(10,*)
	write(10,*)
	close(10)
	return
	end 	subroutine writecorrt
	!----------------------------------------
	subroutine writecorrw(task)
	integer, intent(in) :: task
	integer :: i
	if(task==101) then
	 open(10,file='a-vs-w.dat',
     .             action='write',position='append')
		write(10,'(a)')'# Photon'
	elseif(task==102) then
	 open(10,file='pl-vs-w.dat',
     .             action='write',position='append')
		write(10,'(a)')'# Photon'
	elseif(task==105) then
	 open(10,file='dd-vs-w.dat',
     .             action='write',position='append')
		write(10,'(a)')'# Exciton |up><up|'
	elseif(task==106) then
	 open(10,file='dd-vs-w.dat',
     .             action='write',position='append')
		write(10,'(a)')'# Exciton |dn><dn|'
	endif
	do i=1,nw
		write(10,'(f20.10,3x,2E20.10)') ws(i),
     .             real(corrw(i)), aimag(corrw(i))
	end do
	write(10,*)
	write(10,*)
	close(10)
	
	return
	end 	subroutine writecorrw
	!----------------------------------------
	subroutine mkfname(n)
	implicit none
	integer, intent(in) :: n
	character :: dir*30
	logical :: ex
	
		if(n<10) then
			write(dir,'(I1)') n
		elseif(n<100) then
			write(dir,'(I2)') n
		elseif(n<1000) then
			write(dir,'(I3)') n
		elseif(n<10000) then
			write(dir,'(I4)') n
		else
			write(*,*)'ERROR: n too large!!!'
			stop
		endif

	dir = 'n-'//trim(dir);
	inquire (file=trim(dir), exist=ex)
	!if(.not. ex) stop 'data directory does not exist!'

	fname = trim(dir)//'/wfs.dat'
	write(6,'(a)') fname
	
	return
	end subroutine mkfname
	!-------------------------------------------
	! writes w or t 1-dimensional data: corrw, corrt
	subroutine writeatnode(ij1,nw,ws, corrw,filename)
	use modmain, only: node
	implicit none
	integer, intent(in) :: ij1, nw
	double complex, dimension(nw), intent(in) :: corrw
	double precision, dimension(nw), intent(in) :: ws
	character(len=*), intent(in) :: filename
	integer :: i
	character :: rank*30, fname*100
	
	write(rank,'(i6.6)') node
	fname = trim(filename)//'-'//trim(rank)
	! write unformatted file
	if(ij1==0) then
		open(1,file=trim(fname), form="formatted", action="write")
	else
		open(1,file=trim(fname), form="formatted", action="write",
     .                                      position="append")
	endif


	do i=1,nw
		write(1,'(f20.10,3x,2E20.10)') ws(i),
     .             real(corrw(i)), aimag(corrw(i))
	end do
	write(1,*)
	write(1,*)


	close(1)
	return
	end 	subroutine writeatnode
	!=============================================================
	subroutine rwallnodesx(filename,fout,nw)
	use modmain, only: node, jobs, njobs, num_procs
	implicit none

	integer, intent(in) :: nw
	character(len=*), intent(in) :: filename, fout
	! local
	double precision, dimension(nw,3) :: dat
	integer :: i, j
	character :: rank*30, fname*100

	dat = 0.0d0;
	
	! output file with data from all nodes
	fname = trim(fout)//'.dat'
	! read unformatted file
	open(2,file=trim(fname), form="formatted", action="write",
     .      position="append") ! remove append?

	do i=0,min(njobs,num_procs)-1
		write(rank,'(i6.6)') i
		fname = trim(filename)//'-'//trim(rank)
		! read node specific output file
		open(1,file=trim(fname), form="formatted", action="read")
		do ij=jobs(i)%i1,jobs(i)%i2!1,jobs(i)%njobs
			! read
			do j=1,nw
				read(1,*) dat(j,1:3)
			end do
			read(1,*) 
			read(1,*) 
			! write
			do j=1,nw
				write(2,'(f20.10,3x,2E20.10)') dat(j,1:3)
			end do
			write(2,*)
			write(2,*)
		end do
		close(1, status='delete')
		!shift = jobs(i)%njobs
	end do
	close(2)

	return
	end 	subroutine rwallnodesx
	!=============================================================


	end module correlation


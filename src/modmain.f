

! Copyright (C) 2018 M. Ahsan Zeb
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

	module modmain
	implicit none

	! to set kind of integers in ksub
	! kind 1 (8bits) can go upto 127, which should be enough 
	! for us if nsites <= 127.
	!	kind 2 (16 bits) max is 32,767
	! selected_int_kind(R) ==> smallest kind, to 10^R exclusive 
	integer, parameter :: isk=selected_int_kind(2);

	! tasks(maxtasks)
	!integer, parameter :: maxtasks = 100; ! some large numeber
	!integer, dimension(maxtasks) :: tasks
	integer :: task
	integer :: n, mv, nph
	double precision :: wr, delta, lambda, wv,wc ! set in main? fro ith job
	integer :: ntot, ntotb
	integer :: nact, nsym
	double precision :: uscfac ! Ultra-strong coupling factor;
	
	! time evolution to get correlation
	logical :: td, fft !, melem, hopresp, fixrhoex, amelem
	integer :: chi
	double precision :: dt,w1,w2, kappa2
	integer :: nt,nw, prntstep, nstepf

	logical :: debug

	logical :: ddiagOK ! direct diagonalisation OK? calc and saves wavefunctions... 
	integer :: node, num_procs, njobs ! number of processors, number of jobs
	integer :: mynode
	!integer(kind=MPI_OFFSET_KIND), dimension(:,:), allocatable :: iodisp 

	integer :: memlim ! memory limit per node, in bytes

	logical :: writewfs ! write wavefunctions
	logical :: writecdms	! calc and write cond density matrices
	logical :: writebmap ! wirte basis and map files and stop
	logical :: onlyenergy ! write only energy eigenvalues 
	
	type :: njblocksforcmds
		integer :: ij1, nj
	end type njblocksforcmds
	type(njblocksforcmds), dimension(:), allocatable :: njpntr


	
	! to manage independent parallel jobs on num_procs nodes.
	type :: cluster
		integer :: njobs, i1, i2
	end type cluster
	type(cluster), dimension(:), allocatable :: jobs

	type :: parameters
		double precision :: wr, del, lam, wv, j
		integer :: ntot, ntotdn, ntotb, ntotup, ntotg
	end type parameters
	type(parameters), dimension(:), allocatable :: param

	integer :: nwr, ndel, nlam, nwv
	double precision, dimension(:), allocatable::
     .             wrs, dels, lams, wvs
	double precision :: lamd, lamd2wv    




	type :: BSectors
		integer :: ntot ! will see later if these ntot are needed?
		integer(kind=isk), allocatable :: sets(:,:)
		double precision, allocatable :: Nv(:), P(:), r(:,:) ! actually integers, but define real to avoid upper limit.; 
		!r for ratios:
		!# avoids large factorials
		!# sec(p)%r[m*,j_{p-1}] = sqrt( P_{p-1}[j]/P_p[k] ); where k=map[m*,j];

		integer, allocatable :: f(:,:) ! to store freq of occupations for all perm symm basis; used in Hb calc.	
	end type BSectors

	type :: BasisSet
		integer :: ntot
		!logical :: xst
		! number of active sites
		!integer :: n
		integer :: maxk ! max k of k-sub it has.
		integer, allocatable :: pntr(:)
		type(BSectors), allocatable :: sec(:)
	end type BasisSet

	type :: Eigensystems
		integer :: ntot, n1,n2 ! size of hilbert space, =nrows, ncols=nevecs
		double precision, allocatable :: eval(:)
		double precision, allocatable :: evec(:,:)
		integer, allocatable:: par(:) ! parity: 0/1 for even/odd
	end type Eigensystems

	type :: HSectors
		integer :: ntot, nnu, nmu, nnz
		integer, allocatable :: row(:)
		integer, allocatable :: col(:)
		!integer, allocatable :: rowpntr(:)
		double precision, allocatable :: vdat(:)
		!double precision :: hg ! fixed number, pth-sec: Sqrt[(m-p)(p+1)(N-p)]
	end type HSectors


	type :: Ham
	!Hg%ntot,Hg%nnz, Hg%coodat, Hg%coo1, Hg%coo2, Hg%dat, Hg%col, Hg%rowpntr
		logical :: xst,dense
		integer :: n,m,m1
		integer :: nev, ncv
		integer :: ntot ! HilbertSpace dimension
		integer :: nnz  ! no of non-zero elements
		integer :: srptr ! size of row pointers, not eq to ntot
		! to hold data before CSR conversion
		integer, allocatable :: coo1(:), coo2(:)
		double precision, allocatable :: coodat(:)
		! final CSR format 
		integer, allocatable :: col(:)
		integer, allocatable :: rowpntr(:)
		double precision, allocatable :: dat(:)
		! direct diagon, dense form
		double precision, allocatable :: h(:,:)
		!integer, allocatable :: row(:)
		 ! pointers to ksub sectors in rowpntr
		integer, allocatable :: spntr(:)
		! to hold data for largest m1 possible for reuse for smaller m1 values
		type(HSectors), allocatable :: sec(:)
	end type Ham


	type(BasisSet):: basis
	! change name of HilbertSpace to something like eigensystem ??
	!	13 hamiltonians and eigensystems
	type(Eigensystems), dimension(:), allocatable :: eig

	type(Eigensystems) :: eig0

	double precision, dimension(:,:), allocatable :: eigp ! parity eigenstates

	
!	type(Ham), allocatable:: Hg
	type(Ham) :: Hhtc, Hg, Hb, Hf ! off-diagonal terms
	type(Ham) :: Hg1, Hb1, Hg1f, Hb1f ! off-diagonal terms

	!diagonal terms, simple forms, maybe taken as simple arrays... 
	double precision, dimension(:), allocatable :: Hc,Hs,Ht, Hdv
	double precision, dimension(:), allocatable :: Hc1f,Hs1f,Ht1f

	double precision, dimension(:,:,:,:), allocatable :: dmup, dmdn


	integer, dimension(:), allocatable :: origin
	integer :: nev ! number of eigenvectors to compute

	
	! indexes: map(0:mv, 1,ntotmaxsites), keep in mind when allocate it
	integer, dimension(:,:), allocatable :: map ! stores map of p-1 sites to p sites.

	! above this size, iterative arpack solver will be used
	!	to diagonalise Hamiltonians
	integer:: nmaxddiag ! max size for direct diagonalisation
	integer :: diagmaxitr! key to set it in input: DirectSolverSize
	
	contains
	!===================================================================
	integer (kind=16) function factorial(p) 
	! length of block that would be added
	! to the map if we add a site to p-1 sites
	implicit none
	integer, intent(in):: p
	integer :: i
	factorial = 1;
	do i=2,p
		factorial = factorial*i
	end do
	return
	end function factorial
	!===================================================================
	integer (kind=16) function binomial(n,m) 
	! length of block that would be added
	! to the map if we add a site to p-1 sites
	implicit none
	integer, intent(in):: n,m
	integer :: i, smaller, greater

	if(n<m) then
		binomial = 0;
	elseif(n==m) then
		binomial = 1;		
	else
		!binomial = factorial(n)/(factorial(m)*factorial(n-m))	
		! remove the largest factorial
		if(m .ge. n-m) then
			greater = m
			smaller = n-m
		else
			greater = n-m
			smaller = m
		endif

		binomial = 1
		do i=greater+1,n
				binomial = binomial*i
		end do
		binomial = binomial/factorial(smaller)
			
	endif
	return
	end function binomial
!--------------- timer -----------------
	double precision function clock()
		real*4:: etime, tm(2)
		clock = etime( tm )
	return
	end function clock
!---------------------------------------



	Subroutine timestamp(node)
!--------------------------------------------------------------------
! 	prints welcome message with date and time
!--------------------------------------------------------------------

	implicit none
	integer, intent(in):: node
	integer	:: dt(8), x
	character*10 :: bb(3)
	double precision :: dummy
	logical, save :: start = .true.;
	integer, save	:: dti(8)
!--------------------------------------------------------------------

	call date_and_time(bb(1), bb(2), bb(3), dt)
	! start random sequence with seed = milliseconds in current time
	call srand(dt(8)) !rand(dt(8))
	if(node .ne. 0 ) return


	if (start) then
		dti = dt;
		start = .false.;
		!call srand(0)
		!write(6,*) 'random = ', rand(0)		
		write(6,*)
		!write(6,*)"**************************"//
    ! .		"******************************"
		write(6,'(/,a)')"			TransOC Started! "
		!write(6,'(/,a)')" TransOC: Transport in Organic Cavities"
		!write(6,'(/,a)')
    ! .        '	        . . . Running v 0.0 . . . '
		write(6,'(/,a,i2,a,i2,2a,i2,a,i2.2,a,i4,a)')
     . 	"	   ",
     .                   dt(5)," Baj kar ",dt(6), " min",
     .               " aor tarekh hai ",dt(3),"/",dt(2),"/",dt(1)
		write(6,*)"         -----------------"//
     .		"----------------------------"
	else
		! assume prog will run on the same day
		x= dt(3)*86400+dt(5)*3600 + dt(6)*60 + dt(7) ! sec
		x = x - (dti(3)*86400+dti(5)*3600 + dti(6)*60 + dti(7));
		dt(3) = x/86400; ! days
		x = x - dt(3)*86400;
		dt(5) = x/3600; ! hours
		x = x - dt(5)*3600;
		dt(6) = x/60; ! min
		x = x - dt(6)*60;
		dt(7) = x; ! sec

		write(6,
     . '("Time taken: ",i2,":",i2.2,":",i2.2,":",i2.2,"(d:h:m:s)")')
     .   dt(3),dt(5:7)
	endif

	return
	END Subroutine timestamp
!========================================================	
	subroutine input(node)
	implicit none
	integer, intent(in):: node
	! local
	character(40) :: block
	integer :: iostat
	logical :: file_exists
	logical :: parameters
	double precision :: memory
	integer :: i, nevdm, nexdm, anev
	logical :: fixrhoexdm
	
	debug = .false.
	! set default
	nmaxddiag = 20
	diagmaxitr = 500
	n = 2
	mv = 2
	nph=10;
	wc = 2.0d0; ! singlet energy
	parameters = .false.
	memory = 1.00 ; ! 1 GB per node
	nev = 1;
	writewfs = .false.
	writecdms = .true.
	writebmap = .false.
	onlyenergy = .false.
	anev = 1;
	ddiagOK = .true.;

	td= .true.;
	chi= 1;
	fft=.true.
	nt = 1000;
	nw=500;
	w1=-0.25d0;
	w2=+2.5d0;
	dt = 0.01;
	kappa2=0.05;
	prntstep = 1000;
	uscfac = 1.0d0;
	nstepf = 10.0d0;
	!--------------------------!
	!     read from input.in   !
	!--------------------------!
	INQUIRE(FILE="input.in", EXIST=file_exists)
	if (.not. file_exists)	then
		write(6,*) 'Error: input.in does not exist!'
		stop
	endif

	open(50,file='input.in',action='READ',
     . 					status='OLD',form='FORMATTED',iostat=iostat)
	if (iostat.ne.0) then
		write(*,*)
		write(*,'("Error(readinput): error opening input.in")')
		write(*,*)
		stop
	end if

10		continue

	read(50,*,end=30) block
	! check for a comment
	if ((scan(trim(block),'!').eq.1).or.
     .  (scan(trim(block),'#').eq.1)) goto 10

	select case(trim(block))

	! task: 
	! Response functions: 100 series
	! 101 abs
	! Density Matrices: 300 series
	! 310 photon reduced dm

	case('task') 
	 read(50,*,err=20) task

	case('nstates') 
	 read(50,*,err=20) nev

	case('USCfac') 
	 read(50,*,err=20) uscfac ! logical

	case('TimeEvolutionParam')
		read(50,*,err=20) kappa2, dt, nt, nw, w1,w2, fft
		kappa2 = kappa2/2.0d0; ! just set it here....
	case('TimeEvolutionPrintStep')
		read(50,*,err=20) prntstep

	case('TimeEvolutionFissionStep')
		read(50,*,err=20) nstepf

	case('nph') 
	 read(50,*,err=20) nph

	case('w0') ! singlet energy, fixed, does not change with detuning.
	 read(50,*,err=20) wc

	case('N')
		read(50,*,err=20) nact !n_original
		if(nact .le. 1) then
			stop 'Error(input): Set N to at least two sites.... '
		endif
	case('debug')
		read(50,*,err=20) debug
	case('parameters')
		parameters = .true.
		read(50,*,err=20) nwr, ndel, nlam, nwv
		!write(*,*)"allocated(wrs) = ",allocated(wrs)
		allocate(wrs(nwr))
		allocate(dels(ndel))
		allocate(lams(nlam))
		allocate(wvs(nwv))
		read(50,*,err=20) wrs
		read(50,*,err=20) dels
		read(50,*,err=20) lams
		read(50,*,err=20) wvs
	case('DirectSolverSize')
		read(50,*,err=20) nmaxddiag
	case('IterSolverMaxIter')
		read(50,*,err=20) diagmaxitr
	case('')
		goto 10
	case default
	write(*,*)
	write(*,'("Error(readinput): invalid block: ",A)')trim(block)
	write(*,*)
	stop
	end select
	goto 10
20		continue
	write(*,*)
	write(*,'("Error(readinput): error reading from input.in")')
	write(*,'("Problem occurred in ''",A,"'' block")') trim(block)
	write(*,'("Check input convention in manual")')
	write(*,*)
	stop
30			continue
	close(50)


	 	nsym = nact;
	  n = nsym


	! parameters block must be given by user
	if(.not. parameters) then
	write(*,*)
		write(*,'("Error(readinput): parameters not given in input.in")')	
		stop
	endif

	return


99		write(*,*)
	write(*,'("Error(readinput): invalid block: ",A)')trim(block)
	write(*,*)
	stop		

	
	end subroutine input
!==================================================


	end module


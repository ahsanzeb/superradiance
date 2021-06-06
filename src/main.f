

! Ahsan Zeb, 20 oct 2018
! HTC model with N molecules, m exciton-photon excitations, M+1 vibrational levels.

! We will use permutation symmetry:
! blocks with respect to number of up spins:
!	pth block: p up spins.
!					p up sites are permutaiton symmetric; N-p down sites are perm symm. 
!					m1=Min(m,N)+1 total blocks; i.e., p=0,1,2,3,...,m1
! further details in doc.


	program super
	use modmain
	use maps, only: getmap, writemap
	use bases, only: PermSymBasis, writebasis, writebasisf
	use hamiltonian, only: MakeHhtc, HamParts
	use diag, only: diagonalise
	use mpi
	use dmat, only: rwallnodes, rdmmol, rdmf
	
	implicit none

	double precision:: tottime, dtau, ed,ew
	logical :: goforcdms
	! mpi
	integer:: ierr, nj
	integer(kind=MPI_OFFSET_KIND),
     .           dimension(:), allocatable :: disp 

	! 
	integer :: i, ijob, ij1, njl,k, j,ll,il,idark,nnzcoef
	integer :: ix,jux,ind 
!======================================================================
	call MPI_INIT(ierr)
	!find out MY process ID, and how many processes were started.
	call MPI_COMM_RANK (MPI_COMM_WORLD, node, ierr)
	call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

	if(node==0) 	call timestamp(node) ! time stamp
	! set mv=2 for 3 states of molecules: G,T,S
	mv = 2;

	!write(*,*) "Node = ",node," num_procs, ntp = ",num_procs, ntp 
	! readinput file
	call input(node)

	! just write basis and map?
	if(writebmap) then
		call MPI_FINALIZE(ierr) ! done only at master node
		if(node==0) then
			call getmap(n,mv)
			call PermSymBasis(n,mv)
			! write output files
			call writemap()
			call writebasis()
		endif
		stop "htc: basis and map written to files... done!"
	endif


	! things need to be done only once	
	! allocate space for basis, and calc maps, etc.
	!---------- called once for a given n,mv ---------
	call getmap(n,mv) !,ntot) ! ntot output
	call PermSymBasis(n,mv)

	call setparam()
	
	call setjobs(node)

	!total jobs on this node
	nj = jobs(node)%njobs

	!allocate eig
	allocate(eig(nj))

	open(114,file="energy.dat",action='write')

	ij1 = 0;
	!+++++++++++++++++++++START IJOB LOOP ++++++++++++++++++++++++++
	do i=1,nj,1

		!write(*,*)'jobs(node)%i1 , i2= ',jobs(node)%i1,jobs(node)%i2
		if(node==0)write(*,'(a,i5,a,i5)')'Node 0: job ',i,' out of ',nj
		ijob = jobs(node)%i1+i-1

	  call groundstate(i, ijob, ij1)


	 !call dmphot()
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if(1==1)  then
	! temporary, to write evals. run in serial
	write(6,'(a5x100000f25.15)') "EVALS: ",
     .  	eig(i)%eval(1:min(20,eig(i)%n2))
	!write(*,*)'eig(i)%n2 = ',eig(i)%n2

	write(114,'(100000f25.15)') eig(i)%eval
	!write(*,'(100000f10.3)') eig(i)%evec(:,2)

	endif

	end do ! i jobs
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	call MPI_FINALIZE(ierr)

	! read files	written by all nodes and write a single file.
	if(node==0) then
	 ! combine all dm files written by diff nodes

	 !call rwallnodes("dmmol",n,2)
	 !call rwallnodes("dmfield",n,nph)
	 close(114)

		write(*,*)"Super: everything done.... " 
		call timestamp(node)

	endif

!======================================================================
	contains
!======================================================================
	subroutine printnode(node)
	implicit none
	integer, intent(in) :: node
	if(node==0) then
		write(*,'("main: this is the master node ",I5)')node
	else
		write(*,'("main: this is the slave node ",I5)') node
	endif
	return
	end 	subroutine printnode
!=============================================	
	subroutine CheckError(ierr, node, oper)
	implicit none
	integer, intent(in):: ierr, node
	character(len=*):: oper
	if(ierr .ne. 0) then
		write(*,'("main: mpi operation ",a," at node ",I3," ier = ",I3)')
     .          oper, node, ierr
		stop
	endif
	
	! comment
	!write(*,'("main: mpi operation ",a," at node ",I3," ier = ",I3)')
  !   .          oper, node, ierr
	return
	end 	subroutine CheckError
!=============================================

!=============================================
	subroutine setparam()
	! sets param(:).
	! every ndoes sets its own copy, and uses its part of the job,
	! node=node will use param(jobs(node)%i1:jobs(node)%i2)
	implicit none
	integer :: i,i1,i2,i3,i4,i5

	! total number of jobs
	njobs = nwr * ndel * nlam * nwv;
	allocate(param(njobs))

	!if(node==0) write(*,*) 'njobs = ',njobs

	i = 0;
		do i2=1,nwr,1
			do i3=1,ndel,1
				do i4=1,nlam,1
					do i5=1,nwv,1
						i = i + 1;
						param(i)%wr = wrs(i2) ! omega_R
						param(i)%del = dels(i3) ! delta for omega_c, photon energy
						param(i)%lam = lams(i4) ! lam for soc
						param(i)%wv = wvs(i5) ! wv for omega_T triplet energy
						!param(i)%wc = wc
					end do
				end do
			end do
		end do

	return
	end 	subroutine setparam
!=============================================

!=============================================
	subroutine setjobs(node)
	implicit none
	integer, intent(in) :: node
	integer :: i,j,k

	allocate(jobs(0:num_procs-1))
	
	jobs(:)%njobs = 0; ! set on all processes,
									! independently, and used only for the relevant node
	jobs(:)%i1 = 0; 	jobs(:)%i2 = 0;
	
	! assign jobs one by one
	k =0; ! do i=1,large can be replaced by an open do loop, or a do while loop.
	do i=1, njobs ! i=1, a large dummy number(something equal or bigger than njobs)
		do j=0,num_procs-1
			k = k + 1;
			jobs(j)%njobs = jobs(j)%njobs + 1
			if(k == njobs) exit
		enddo
		if(k == njobs) exit
	end do

	!if(node==0) write(*,*) 'jobs(:)%njobs = ',jobs(:)%njobs

	! assign interval range, only jobs(node) will actually be used later by this node.
	jobs(0)%i1= 1;
	jobs(0)%i2 = jobs(0)%i1 + jobs(0)%njobs - 1; 
	do j=1,num_procs-1,1
		jobs(j)%i1 = jobs(j-1)%i2+1
		jobs(j)%i2 = jobs(j)%i1 + jobs(j)%njobs - 1
	enddo

	return
	end 	subroutine setjobs

!=============================================
	subroutine checkifgoforcdms(nj,i,goforcdms)
	implicit none
	integer, intent(in) :: nj,i
	logical, intent(out) :: goforcdms
	integer :: valueRSS, k

	goforcdms = .false.;
	! memory limit read from input file
	
	! estimate the size of memory allocated
	!call system_mem_usage(valueRSS)
	do k=1,nj
		valueRSS = valueRSS + size(eig(k)%evec)*8; ! 8 bytes for a double precision
	end do
	!write(*,*)'size of eig, memlim = ',valueRSS, memlim
	!if (valueRSS == -1 ) stop "main: valueRSS not fouond!"
	if ((i.eq.nj) .or. (valueRSS .ge. memlim)) then
		goforcdms = .true.
	endif

	return
	end 	subroutine checkifgoforcdms
!=============================================
	subroutine system_mem_usage(valueRSS)
	implicit none
	!use ifport !if on intel compiler
	integer, intent(out) :: valueRSS
	character(len=200):: filename=' '
	character(len=80) :: line
	character(len=8)  :: pid_char=' '
	integer :: pid
	logical :: ifxst	

	valueRSS=-1    ! return negative number if not found

	!--- get process ID
	pid=getpid()
	write(pid_char,'(I8)') pid
	filename='/proc/'//trim(adjustl(pid_char))//'/status'

	!--- read system file
	inquire (file=filename,exist=ifxst)
	if (.not.ifxst) then
		write (*,*) 'system file does not exist'
		return
	endif

	open(unit=100, file=filename, action='read')
	do
		read (100,'(a)',end=120) line
		if (line(1:6).eq.'VmRSS:') then
			read (line(7:),*) valueRSS
			exit
		endif
	enddo
120	continue
	close(100)
	return
	end subroutine system_mem_usage
!=============================================
!=============================================================
!=============================================
! write lowest nev eigenvectors, and eigvals etc
	subroutine iowfnode(i,ij) ! i,ij: local and global job indeces
	implicit none
	integer, intent(in) :: i,ij
	character :: rank*30, fname*100
	
	write(rank,'(i6.6)') node
	fname = 'wf'//'-'//trim(rank)
	! write unformatted file
	if(i==1) then
		open(1,file=trim(fname), form="unformatted", action="write")
	else
		open(1,file=trim(fname), form="unformatted", action="write",
     .                                      position="append")
	endif
	write(1)param(ij)%wr,
     .      param(ij)%del,param(ij)%lam, param(ij)%wv
	!...................... basis%pntr ............
	write(1) basis%pntr
	!...............................................
	write(1) eig(i)%ntot, eig(i)%n2
	write(1) eig(i)%eval 
	write(1) eig(i)%evec	
	close(1)
	return
	end 	subroutine iowfnode
	!=============================================================






	!=============================================================
	subroutine iowfcombine()
	implicit none
	integer :: i,ij,k,j, ntot, nev, nsize
	character :: rank*30, fname*100
	!local
	integer :: nnz
	integer, dimension(:), allocatable :: coo,pntr
	double precision, dimension(:), allocatable :: val,vec,vx
	double precision, dimension(4) :: rvals
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
	if(.not. ex) call system('mkdir '//trim(dir))

	! write output file with data from all nodes
	fname = trim(dir)//'/wfs.dat'
	open(10,file=trim(fname), form="unformatted", 
     .                          action="write",position='append')

	nsize = 0;
	allocate(val(1))
	allocate(vec(1))

	do i=0,min(njobs,num_procs)-1
		write(rank,'(i6.6)') i
		fname = 'wf-'//trim(rank)
		! read unformatted file
		open(1,file=trim(fname), form="unformatted", action="read")
		do ij=1,jobs(i)%njobs
			read(1) rvals
			write(10) rvals

			allocate(pntr(0:n+1))
			read(1) pntr
			write(10) pntr
			deallocate(pntr)
			
			read(1) ntot, nev
			write(10) ntot, nev

			if( nsize /= ntot*nev) then
				nsize=ntot*nev
				deallocate(val,vec)
				allocate(val(nev))
				allocate(vec(nsize))		
			endif
			read(1) val 
			read(1) vec	
			! write to final output file
			write(10) val 
			write(10) vec	
		end do
		close(1, status='delete')
	end do
	close(10)

	return
	end 	subroutine iowfcombine
	!=============================================================	
	subroutine checktrace(l,a)
	implicit none
	integer, intent(in) :: l
	double precision, dimension(l,l), intent(in) :: a
	double precision :: tr
	integer :: i
	
	tr =0.0d0
	do i=1,l
		tr = tr + a(i,i)
	end do

	write(*,*) 'tr = ',tr
	
	return
	end 	subroutine checktrace
	!=============================================================	
	!.....................................................
	! writes output: node=0 combines all output files
	!.....................................................
	subroutine writeout()
	implicit none
	 call rwallnodes('dmmol',nact,2)
	 call rwallnodes('dmfield',nact,nph)
	return
	end subroutine writeout
	!=============================================================	


	subroutine groundstate(i, ijob, ij1)
	implicit none
	integer, intent(in) :: i, ijob
	integer, intent(inout) :: ij1

	! first job? calculate Hamiltonian's parts. 
		if(i==1) then
			call HamParts(nsym,nph)
		endif
	! things need to be done for every job
	call MakeHhtc(nsym, ijob) ! ijob===> sets wr,delta,lambda,wv values
	! diagonalise
	call diagonalise(i)

	return
	 ! calc dms to free mem or wait for more jobs?
	 call checkifgoforcdms(nj,i,goforcdms)
	 if(goforcdms) then
		njl = i-ij1
		 !write(*,*) "ij1, njl = ",ij1, njl
		 call rdmmol(ij1, njl,n,nph,nev)
		 call rdmf(ij1, njl,n,nph,nev)
		! reset variables for next iteration
		goforcdms = .false.;
		ij1 = i;
	 end if


	return
	end 	subroutine groundstate
	!.....................................................



	end program super

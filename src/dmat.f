

	module dmat
	use modmain
	implicit none

	public:: cdms, rwallnodes, cdms1, dmphot
	private:: cdmup, cdmdn 

	contains

	!--------------------------------------------------------------------
	! p-th block cdm for spin up
	!--------------------------------------------------------------------
	subroutine cdmup(ij1,nj,n,p,mv,nev,dm)
	implicit none
	integer, intent(in) :: ij1,nj, n,p, mv,nev ! ij1: now nj can be smaller than actual total nj for a node
	double precision,dimension(mv+1,mv+1,nj,nev),intent(out) :: dm
	integer :: l1,l2, lout, nnp,ntotp
	!double precision, dimension(:,:), allocatable :: v
	integer :: ii,i1,j,j1,i,jj,i11,j11,ij, ijob, m1ij,k,xj
	double precision :: fac, Pii, Pi1,Pj1,fac2
	type(Eigensystems), dimension(:), allocatable :: v
	
	l1 = basis%pntr(p); l2 = basis%pntr(p+1)

	!write(*,*)'p, l1, l2 = ',p, l1,l2
	!lout = basis%sec(p)%ntot;
	nnp = basis%sec(n-p)%ntot;
	!allocate(v(lout))
	! tracemu(l,lout,nnp,v) ! traces over mu type
	!v = tracemu(l2-l1, lout, nnp, eig%evec(l1:l2-1,1)) 
	!dm = trmol(n,p,mv,lout,v) ! function to trace over p-1 other molecules' vibrational states

	ntotp = min(basis%sec(p-1)%ntot,basis%sec(ndummy-1)%ntot); ! states for which p block states are available to map onto... 

	dm =0.0d0

	allocate(v(nj))
	k = 0
	do i=ij1+1,ij1+nj
		k = k + 1
		ijob = jobs(node)%i1 + i -1
		m1ij = min(param(ijob)%m, n)
		if(p .le. m1ij) then
			allocate(v(k)%evec(l2-l1,nev))
			v(k)%evec = eig(i)%evec(l1:l2-1,:)
		endif
	enddo
	!allocate(v(l2-l1,nj)) ! lout*nnp = l2-l1
	!do ij=1,nj
	!	v(:,ij) = eig(ij)%evec(l1:l2-1)
	!end do
	!================================================
	do ii=1, ntotp !basis%sec(p-1)%ntot
		!Pii = basis%sec(p-1)%P(ii)
		do i=0,mv
			i1 = map(i,ii) ! {nu|mu}_(k-1) + i ==> i1-th basis state of {nu|mu}_k
			!Pi1 = basis%sec(p)%P(i1)
			!write(*,*)'i1 = ', i1
			!fac = Pii/dsqrt(Pi1); ! fac due to vib permutations
			fac = basis%sec(p-1)%r(i,ii) ! r(i,ii) = dsqrt(Pii/Pi1)
			i1 = (i1-1)*nnp;
		
			do j=i,mv ! upper triangular only
				j1 = map(j,ii) 
				!Pj1 = basis%sec(p)%P(j1)
				!fac2 = fac/dsqrt(Pj1); ! fac due to vib permutations
				fac2 = fac*basis%sec(p-1)%r(j,ii) ! r(j,ii) = dsqrt(Pii/Pj1)
				!write(*,*)'j1 = ', j1
				j1 = (j1-1)*nnp;
				k = 0
				do ij=ij1+1,ij1+nj,1
					k=k+1
					ijob = jobs(node)%i1 + ij - 1;
					if(p .le. min(param(ijob)%m,n)) then
						do xj=1,nev
						dm(i+1,j+1,k,xj) = dm(i+1,j+1,k,xj)
!     .         + fac2*sum(v(i1+1:i1+nnp,ij)*v(j1+1:j1+nnp,ij))
     . + fac2*sum(v(k)%evec(i1+1:i1+nnp,1)*v(k)%evec(j1+1:j1+nnp,xj))
     				enddo ! xj final state from LP
					endif
				end do
			enddo
		enddo
	end do ! ii
	! from perm sym of exciton-photon states
	dm = dm*p*1.0d0/n ! p out of n are up
	!================================================
	
	deallocate(v)
	return
	end subroutine cdmup


	!--------------------------------------------------------------------
	! full cdms combine all blocks
	!--------------------------------------------------------------------
	subroutine cdms(ij1, nj,n,m,mv,nev)
	implicit none
	integer, intent(in) :: ij1, nj,n,m,mv,nev
	double precision,dimension(mv+1,mv+1,nj,nev):: dmup,dmdn
	double precision,dimension(mv+1,mv+1,nj,nev) :: dm
	integer :: p,m1,ij,m1max, k
	double precision :: tr
	
	dmup = 0.0d0;
	do p=0,nph
			call cdmup(ij1,nj,n,p,mv,nev,dm)
			dmup = dmup + dm;
		endif	
	end do


	! write output files at each node
	call writeatnode(ij1,nj,dmup,'dmup')

	
	return
	end subroutine cdms
	!--------------------------------------------------------------------
















!======================================================================
	subroutine dmphot(ij1,nj,n,nph,nev,dm)
	implicit none
	integer, intent(in) :: ij1,nj, n,nph, nev 
	double precision,dimension(nph+1,nph+1,nj,nev),intent(out) :: dm
	integer :: ntotb,k1,k2,k1l,k2l,i,k1i,k2i,ij
	double precision :: 
	
	ntotb = basis%sec(n)%ntot; 
	! size of the molecular block for every photon state

	dm =0.0d0

	do k1=0,nph ! photon states
	 k1l = k1*ntotb
	 do k2=0,nph ! photon states
	  k2l = k2*ntotb
	  do i=1,ntotb ! mol blocks
	   k1i=k1l + i ! global index
	   k2i=k2l + i
	   do ij = ij1+1,ij1+nj ! jobs
	    ! (...,:) for nev lowest eigenstates
	    dm(k1,k2,ij,:) = dm(k1,k2,ij,:) + 
     .       eig(ij)%evec(k1i,:)*eig(ij)%evec(k2i,:)
	   end do! ij
	  end do !i
	 end do ! k2
	end do! k1

	return
	end subroutine dmphot
!======================================================================



































	subroutine writecdmsfiles(mv,nj,nev, dmup, dmdn)
	implicit none
	integer, intent(in) :: mv, nj,nev
	double precision,dimension(mv+1,mv+1,nj,nev),
     .                         intent(in)::dmup,dmdn
	integer :: i,j, k1, k2


	! write output files
	open(10,file='dmdn.dat',action='write',position='append')
	open(12,file='dmup.dat',action='write',position='append')

	do i=1, nj
	 do j=1, nev
	  ! down
	  do k1=1,mv+1
	   write(10,'(1000f25.15)') (dmdn(k1,k2,i,j), k2=1,mv+1)
	  end do
		! up
	  do k1=1,mv+1
	   write(12,'(1000f25.15)') (dmup(k1,k2,i,j), k2=1,mv+1)
	  end do
	 enddo
	enddo
	
	close(10)
	close(12)

	return
	end subroutine writecdmsfiles




	subroutine writeatnode(ij1,nj,dm,filename)
	implicit none
	integer, intent(in) :: ij1,nj
	double precision, dimension(mv+1,mv+1,nj,nev), intent(in) :: dm
	character(len=*), intent(in) :: filename
	integer, dimension(3) :: dsize
	integer :: thefile
	integer :: i,ij,j,xj
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
	do ij=1,nj
		do xj=1,nev
			do i=1,mv+1
				write(1,*) (dm(i,j,ij,xj), j=1,mv+1)
			end do
		enddo ! xj
	end do
	close(1)
	return
	end 	subroutine writeatnode
	!=============================================================
	subroutine rwallnodes(filename,n)
	implicit none
	character(len=*), intent(in) :: filename
	integer, intent(in) :: n
	integer :: i,ij,k,j,shift,xj
	character :: rank*30, fname*100
	!local
	double precision, dimension(mv+1,mv+1,njobs,nev) :: dm
	character :: dir*30
	logical :: ex


	! trim() does not work if you make dir by "write(dir,*) n"
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

	dm =0.0d0

	shift=0;
	do i=0,min(njobs,num_procs)-1
		write(rank,'(i6.6)') i
		fname = trim(filename)//'-'//trim(rank)
		! read unformatted file
		open(1,file=trim(fname), form="formatted", action="read")
		do ij=jobs(i)%i1,jobs(i)%i2!1,jobs(i)%njobs
			do xj=1,nev
				do k=1,mv+1
					read(1,*) (dm(k,j,ij,xj), j=1,mv+1)
				end do
			enddo
			!call checktrace(mv+1,dm(:,:,ij))
		end do
		!close(1)
		close(1, status='delete')
		!shift = jobs(i)%njobs
	end do

	! write output file with data from all nodes
	fname = trim(dir)//'/'//trim(filename)//'.dat'
	! read unformatted file
	open(1,file=trim(fname), form="formatted", action="write",
     .      position="append")
	do ij=1,njobs
		do xj=1,nev
			call symmetrise(mv+1,dm(:,:,ij,xj)) ! complete lower triangular
			do i=1,mv+1
				write(1,*) (dm(i,j,ij,xj), j=1,mv+1)
			end do
		end do
	end do
	close(1)

	return
	end 	subroutine rwallnodes
	!=============================================================

	
	subroutine symmetrise(l,a)
	implicit none
	integer, intent(in) :: l
	double precision, dimension(l,l), intent(inout) :: a
	integer :: i,j
	do i=1,l-1
		do j=i+1,l,1
			a(j,i) = a(i,j)
		end do
	end do
	return
	end 	subroutine symmetrise
	!=============================================================




	end 	module dmat

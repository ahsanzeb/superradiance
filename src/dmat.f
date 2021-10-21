

	module dmat
	use modmain
	implicit none

	public :: rdmmol, rdmf, rwallnodes, parity, parityeig, setparity
	public :: mixparity
	private :: writeatnode
	contains


!======================================================================
	!--------------------------------------------------------------------
	! reduced density matrix of a molecule
	!--------------------------------------------------------------------
	subroutine rdmmol(ij1, nj,n,nph,nev)
	implicit none
	integer, intent(in) :: ij1, nj,n,nph,nev
	double precision,dimension(0:2,0:2,nj,nev) :: dm
	integer :: jj,k1,k2,i1,i2,j1,j2,ij,p,ntotb, is
	double precision :: x1, x2 !, ns, nt
	
	dm = 0.0d0;
	ntotb = basis%sec(n)%ntot; ! size of mol block
	 
	 do jj=1,basis%sec(n-1)%ntot; ! N-1 mol state
	  do k1=0,2; ! target mol states
	   j1 = map(k1,jj);  ! N mol state
	   x1 = basis%sec(n-1)%r(k1,jj)
	   do k2=0,2; ! target mol states
	    j2 = map(k2,jj);
	    x2 = x1 * basis%sec(n-1)%r(k2,jj)
	    do p=0,nph; ! photon states
	     i1 = p*ntotb + j1; ! global index of the basis state with k1 
	     i2 = p*ntotb + j2; ! global index of the basis state with k2 
	     do ij=1,nj! ij1+1,ij1+nj; ! jobs
	      dm(k1,k2,ij,:) =  dm(k1,k2,ij,:) + 
     .          x2 *eig(ij1+ij)%evec(i1,:)*eig(ij1+ij)%evec(i2,:)
	     end do ! ij	
	    end do ! p 
	   end do ! k2
	  end do ! k1
	 end do ! jj
	
	! write output files at each node
	call writeatnode(ij1,nj,2,dm,'dmmol')


	! calculate <ns> and <nt>
	!ns = dm(2,2);
	!nt = dm(1,1);

	! write output - serial version at the moment....
	if(ij1==0) then
		open(13,file="mol-pops.dat", form="formatted", action="write")
	else
		open(13,file="mol-pops.dat", form="formatted", action="write",
     .                                      position="append")
	endif


	do ij = 1,nj !ij1+1,ij1+nj ! jobs
		            !   S            T
	  write(13,*) dm(2,2,ij,1), dm(1,1,ij,1), ! +ve side
     .          dm(2,2,ij,2), dm(1,1,ij,2) !  -ve side
	end do ! ij

	close(13)
	
	return
	end subroutine rdmmol
!======================================================================
	!--------------------------------------------------------------------
	! reduced density matrix of the field mode
	!--------------------------------------------------------------------
	subroutine rdmf(ij1,nj,n,nph,nev)
	implicit none
	integer, intent(in) :: ij1,nj, n,nph, nev 
	double precision,dimension(0:nph,0:nph,nj,nev) :: dm
	integer :: ntotb,k1,k2,k1l,k2l,i,k1i,k2i,ij,is

	double precision,dimension(2) :: a(2), ada(2) ! expectations for the lowest two states
	logical, dimension(nj) :: reorder

	reorder = .false.;
	
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
	   do ij = 1,nj !ij1+1,ij1+nj ! jobs
	    ! (...,:) for nev lowest eigenstates
	    dm(k1,k2,ij,:) = dm(k1,k2,ij,:) + 
     .       eig(ij1+ij)%evec(k1i,:)*eig(ij1+ij)%evec(k2i,:)
	   end do! ij
	  end do !i
	 end do ! k2
	end do! k1

	! calculate <a> and <a^+a>
	! write output - serial version at the moment....
	if(ij1==0) then
		open(12,file="order-param.dat", form="formatted", action="write")
	else
		open(12,file="order-param.dat", form="formatted", action="write",
     .                                      position="append")
	endif

	a = 0.0d0; 
	ada = 0.0d0;

	do ij = 1,nj !ij1+1,ij1+nj ! jobs
	 do is=1,2
	  do i=1,nph; ! i=0 has zero contribution
		 a(is) = a(is) + dsqrt(dble(i)) * dm(i,i-1,ij,is)
		 ada(is) = ada(is) + dble(i) * dm(i,i,ij,is)
	  end do
	 end do ! is

	 ! sort positive and negative <a>:
	 ! write positive <a> first and remember the order for mol pops and dmf.
	 if(a(1) >= a(2)) then
	 	reorder(ij) = .false.
	 	write(12,*) a(1), a(2), ada
	 else
		reorder(ij) = .true.
	 	write(12,*) a(2), a(1), ada(2), ada(1)
	 endif

	end do ! ij
	
	close(12)



	! write output files at each node, considering reorder
	call writeatnodeph(ij1,nj,nph,dm,'dmfield',reorder)




	return
	end subroutine rdmf
!------------------------------------------------------------------
	subroutine writeatnodeph(ij1,nj,mv,dm,filename, reorder)
	implicit none
	integer, intent(in) :: ij1,nj,mv
	double precision, dimension(mv+1,mv+1,nj,nev), intent(in) :: dm
	character(len=*), intent(in) :: filename
	logical,dimension(nj), intent(in) :: reorder ! determined from <a>
	integer, dimension(3) :: dsize
	integer :: thefile
	integer :: i,ij,j,xj, iu1,iu2
	character :: rank*30, fname*100
	double precision:: tr

	tr = 0.0d0	
	write(rank,'(i6.6)') node
	fname = trim(filename)//'-'//trim(rank)
	! write unformatted file
	if(ij1==0) then
		open(1,file=trim(fname//'-1'), form="formatted", action="write")
		open(2,file=trim(fname)//'-2', form="formatted", action="write")
	else
		open(1,file=trim(fname//'-1'), form="formatted", action="write",
     .                                      position="append")
		open(2,file=trim(fname//'-2'), form="formatted", action="write",
     .                                      position="append")
	endif


	
	do ij=1,nj
		if(reorder(ij)) then
	   iu1 = 2; iu2 = 1
	  else
	   iu1 = 1; iu2 = 2
	  endif
	  	  
		!do xj=1,nev
		  !tr = 0.0d0
			do i=1,mv+1
				write(iu1,*) (dm(i,j,ij,1), j=1,mv+1)
				write(iu2,*) (dm(i,j,ij,2), j=1,mv+1)
			end do
			!write(*,*)"ij, is, trace: ", ij, xj, tr
		!enddo ! xj
	end do
	close(1)

	return
	end 	subroutine writeatnodeph


!------------------------------------------------------------------
	subroutine writeatnode(ij1,nj,mv,dm,filename)
	implicit none
	integer, intent(in) :: ij1,nj,mv
	double precision, dimension(mv+1,mv+1,nj,nev), intent(in) :: dm
	character(len=*), intent(in) :: filename
	integer, dimension(3) :: dsize
	integer :: thefile
	integer :: i,ij,j,xj
	character :: rank*30, fname*100
	double precision:: tr

	tr = 0.0d0	
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
		  !tr = 0.0d0
			do i=1,mv+1
			  !tr = tr + dm(i,i,ij,xj)
				write(1,*) (dm(i,j,ij,xj), j=1,mv+1)
			end do
			!write(*,*)"ij, is, trace: ", ij, xj, tr
		enddo ! xj
	end do
	close(1)
	
	return
	end 	subroutine writeatnode

!------------------------------------------------------------------
	subroutine rwallnodes(filename,n,mv)
	implicit none
	character(len=*), intent(in) :: filename
	integer, intent(in) :: n,mv
	integer :: i,ij,k,j,xj
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

	do i=0,min(njobs,num_procs)-1
		write(rank,'(i6.6)') i
		fname = trim(filename)//'-'//trim(rank)
		! read unformatted file
		open(1,file=trim(fname), form="formatted", action="read")
		do ij=jobs(i)%i1,jobs(i)%i2
			do xj=1,nev
				do k=1,mv+1
					read(1,*) (dm(k,j,ij,xj), j=1,mv+1)
				end do
			enddo
		end do
		close(1, status='delete')
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
		!write(*,'(1000f10.5)') (dm(1,j,ij,1), j=1,mv+1)
	end do
	close(1)

	return
	end 	subroutine rwallnodes
!------------------------------------------------------------------

	
	subroutine symmetrise(l,a) ! set lower triangular equal to the upper triangular.
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
!------------------------------------------------------------------

	!--------------------------------------------------------------------
	! check the parity of the eigenstates: w.r.t total number of excitations
	!--------------------------------------------------------------------
	subroutine parity(ij1, nj,n,nph,nev)
	implicit none
	integer, intent(in) :: ij1, nj,n,nph,nev
	double precision, dimension(nj,nev) :: dm
	integer :: jj,k1,k2,i1,i2,j1,j2,ij,p,ntotb,i, par, nx
	double precision :: x1, x2
	
	dm = 0.0d0;
	ntotb = basis%sec(n)%ntot; ! size of mol block
	 
	 do i=1,basis%sec(n)%ntot; ! N-1 mol state
	  nx = basis%sec(n)%f(1,i) + basis%sec(n)%f(2,i); ! number of excitons, S+T
	  do p=0,nph; ! photon states
	   i1 = p*ntotb + i; ! global index of the basis state
		 par = mod(nx + p,2)	   
	   do ij=1,nj! ij1+1,ij1+nj; ! jobs
	     dm(ij,:) = dm(ij,:) + par*dabs(eig(ij1+ij)%evec(i1,:))**2
	   end do ! ij	
	  end do ! p 
	 end do ! jj

	
	write(*,'(a)') "Parities of eigenstates: "
	write(*,'(a)') "ijob, parities: "
	do ij=1,nj
	 write(*,'(i5,3x,1000f8.2)') ij, dm(ij,:)
	end do
	
	return
	end subroutine parity

!------------------------------------------------------------------

	!--------------------------------------------------------------------
	! build the parity eigenstates
	!--------------------------------------------------------------------
	subroutine parityeig(n,nph) ! eigp global in modmain.
	implicit none
	integer, intent(in) :: n,nph
	integer :: ntot, ntotb, i,nx,par,i1, p

	ntotb = basis%sec(n)%ntot; ! size of mol block
	ntot = ntotb * (nph + 1); ! size of our full space.
	if(allocated(eigp)) deallocate(eigp)
	allocate(eigp(ntot,2))
	eigp = 0.0d0

	! build parity eigenstates eigp:
	 do i=1,basis%sec(n)%ntot;! ntotb ! N-1 mol state
	  nx = basis%sec(n)%f(1,i) + basis%sec(n)%f(2,i); ! number of excitons, S+T
	  do p=0,nph; ! photon states
	   i1 = p*ntotb + i; ! global index of the basis state
		 par = mod(nx + p,2)	   
		 if(par == 0) then ! even
	     eigp(i1,2) = 1.0d0;
		 else ! odd
	     eigp(i1,1) = 1.0d0; 
		 endif
	  end do ! p 
	 end do ! i

	return
	end subroutine parityeig
!------------------------------------------------------------------


	!--------------------------------------------------------------------
	! fix the parity of the energy eigenstates by projecting onto 
	! even or odd parity sectors and normalising. 
	!--------------------------------------------------------------------
	subroutine setparity(ij1, nj, nev) ! eigp global in modmain.
	implicit none
	integer, intent(in) :: ij1, nj, nev
	integer :: ij, is
	double precision :: w1, w2, norm
	double precision, dimension(eig(1)%ntot,2) :: proj ! projected states, onto definite parity sectors
	double precision, dimension(eig(1)%ntot,nev) :: auxg, auxu !g,u german
	integer, dimension(nev) :: evenodd

	evenodd = -1;
	! set first to even and second to odd parity; consistent with the low light-matter coupling or normal phase. 
	! the higher eigenstates are sorted according to their larger parity component.
	do ij=1,nj
	 auxg = 0.0d0; auxu = 0.0d0;
	 do is=1,nev
	 	proj = 0.0d0;
	  proj(:,1) = eigp(:,1) * eig(ij1+ij)%evec(:,is);
	  proj(:,2) = eigp(:,2) * eig(ij1+ij)%evec(:,is);
		w1 = sum(dabs(proj(:,1))**2); ! weight on odd
		w2 = sum(dabs(proj(:,2))**2); ! weight on even
		write(*,*)"is, w1,w2 = ",is, w1,w2
	  if(w1 >= w2) then ! set to odd
	   norm = 1.0d0/dsqrt(w1);
	   auxu(:,is) = norm * proj(:,1);
	   evenodd(is) = 1
	   !write(*,*) "Odd: is=",is
	   if (w1==w2) then
	    write(*,*)"Warning(dmat): W1 = W2"
	   endif
	  elseif(w2 > w1)then !if(w2 > w1) then ! set to even
	   !write(*,*) "Even: is=",is
	   norm = 1.0d0/dsqrt(w2);
	   auxg(:,is) = norm * proj(:,2);
	   evenodd(is) = 0
	  endif
	 end do ! is

 
	! get the lowest energy even/odd states, not interested in finding/sorting higher even/odd states at the moment. might do later sometimes...
	! even
	eig(ij1+ij)%evec = 0.0d0
	do is=1,nev
	 if (evenodd(is)==0) then
	  eig(ij1+ij)%evec(:,1) = auxg(:,is);
	  exit
	 endif
	end do
	! odd
	if(nev>1) then
	do is=1,nev
	 if (evenodd(is)==1) then
	  eig(ij1+ij)%evec(:,2) = auxu(:,is);
	  exit
	 endif
	end do
	endif

	
	end do ! ij

	write(*,'(a,100i3)')"evenodd: ",evenodd
!"Warning(dmat): W1 = W2":
! two states, one this one and one another one later or before this, would both be assigned to the same parity if exactly equal wt (w1,w2: double precision, a one digit difference is enough to eliminate this issue) in either sector.
	    ! practically less likely but still there is a chance that this happens.

	return
	end subroutine setparity
!------------------------------------------------------------------






	!--------------------------------------------------------------------
	! makes +,- superpositions of even and odd parity eigenstates to make
	! states that have non-zero expectation of photon annihilation 
	! operator in the superradiance phase.
	! IMPORTANT: call after setparity() has created/resolved even and odd parity states
	!--------------------------------------------------------------------
	subroutine mixparity(ij1, nj, nev) ! eigp global in modmain.
	implicit none
	integer, intent(in) :: ij1, nj, nev
	integer :: ij, is
	double precision :: w1, w2, norm
	double precision, dimension(eig(1)%ntot,2) :: proj ! projected states, onto definite parity sectors

	do ij=1,nj
	 do is=1,nev-1,2
	 	proj = 0.0d0;
	  proj(:,1) = eig(ij1+ij)%evec(:,is);
	  proj(:,2) = eig(ij1+ij)%evec(:,is+1);
	  eig(ij1+ij)%evec(:,is) = proj(:,1) + proj(:,2); ! even + odd
	  eig(ij1+ij)%evec(:,is+1) = proj(:,1) - proj(:,2); ! even - odd
	 end do ! is

	end do ! ij

	return
	end subroutine mixparity
!------------------------------------------------------------------







	end 	module dmat

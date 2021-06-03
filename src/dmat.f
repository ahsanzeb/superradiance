

	module dmat
	use modmain
	implicit none

	public :: rdmmol, rdmf, rwallnodes
	private :: writeatnode
	contains


!======================================================================
	!--------------------------------------------------------------------
	! reduced density matrix of a molecule
	!--------------------------------------------------------------------
	subroutine rdmmol(ij1, nj,n,nph,nev)
	implicit none
	integer, intent(in) :: ij1, nj,n,nph,nev
	double precision,dimension(3,3,nj,nev) :: dm
	integer :: jj,k1,k2,i1,i2,j1,j2,ij
	double precision :: x1, x2
	
	dm = 0.0d0;

	 do jj=1,basis%sec(n-1)%ntot; ! N-1 mol state
	  do k1=1,3; ! target mol states
	   j1 = map(k1,jj);  ! N mol state
	   x1 = basis%sec(n-1)%r(k1,jj)
	   do k2=1,3; ! target mol states
	    j2 = map(k2,jj);
	    x2 = x1 * basis%sec(n-1)%r(k1,jj)
	    do p=0,nph; ! photon states
	     i1 = p*ntotb + j1; ! global index of the basis state with k1 
	     i2 = p*ntotb + j2; ! global index of the basis state with k2 
	     do ij=ij1+1,ij1+nj; ! jobs
	      dm(k1,k2,ij,:) =  dm(k1,k2,ij,:) + 
     .                   x2 *eig(ij)%evec(i1,:)*eig(ij)%evec(i2,:)
	     end do ! ij	
	    end do ! p 
	   end do ! k2
	  end do ! k1
	 end do ! jj

	! write output files at each node
	call writeatnode(ij1,nj,2,dm,'dmmol')

	
	return
	end subroutine rdmmol
!======================================================================
	!--------------------------------------------------------------------
	! reduced density matrix of the field mode
	!--------------------------------------------------------------------
	subroutine rdmf(ij1,nj,n,nph,nev)
	implicit none
	integer, intent(in) :: ij1,nj, n,nph, nev 
	double precision,dimension(nph+1,nph+1,nj,nev) :: dm
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

	! write output files at each node
	call writeatnode(ij1,nj,nph,dm,'dmfield')

	return
	end subroutine rdmf
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
	end do
	close(1)

	return
	end 	subroutine rwallnodes
!------------------------------------------------------------------

	
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
!------------------------------------------------------------------




	end 	module dmat

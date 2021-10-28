
	module hamiltonian
	use modmain	
	implicit none
	
	public :: MakeHhtc, HamParts
	private:: MakeHgbBlock,MakeHgb
	private:: MakeHd, MakeH1R
	public :: coocsr, coo2dense
	contains
 
! ===========================================
	subroutine HamParts(n,nph)
	implicit none
	integer, intent(in):: n,nph

	!write(*,*)'------------------1 '
	call MakeHgbBlock(n)
	!write(*,*)'------------------2 '
	call MakeHgb(n,nph)
	!write(*,*)'------------------3 '
	call MakeHd(n,nph)

	return
	end subroutine HamParts

	!----------------------------------------------------------
	subroutine MakeHgbBlock(n)
	implicit none
	integer, intent(in) :: n
	integer :: p,nnz,i,j,jj,k
	double precision :: f2j,f0i,f1k
	integer :: ntotp, ntotp1

	p = 1; ! for the block for n sites, only one block needed for superradiance code
	if(allocated(Hg%sec)) deallocate(Hg%sec)
	allocate(Hg%sec(p)) ! light-matter coupling: G and S
	allocate(Hb%sec(p)) ! SOC: S and T

		ntotp = basis%sec(n)%ntot
		ntotp1 = basis%sec(n-1)%ntot
		
		nnz = ntotp1
		Hg%sec(p)%nnz = nnz
		allocate(Hg%sec(p)%row(nnz))
		allocate(Hg%sec(p)%col(nnz))
		allocate(Hg%sec(p)%vdat(nnz))

		Hb%sec(p)%nnz = nnz
		allocate(Hb%sec(p)%row(nnz))
		allocate(Hb%sec(p)%col(nnz))
		allocate(Hb%sec(p)%vdat(nnz))

		do jj = 1, ntotp1
			i = map(0,jj); ! combine jj with a molecule in G
			j = map(2,jj); ! combine jj with a molecule in S
			k = map(1,jj); ! combine jj with a molecule in T

			
			Hg%sec(p)%row(jj) = i;
			Hg%sec(p)%col(jj) = j; ! sigma^+ [G to S]

			Hb%sec(p)%row(jj) = k;
			Hb%sec(p)%col(jj) = j; ! sigma^+ tau^- [T to S]

			
			f0i = basis%sec(n)%f(0,i); 
			f2j = basis%sec(n)%f(2,j);
			f1k = basis%sec(n)%f(1,k); 
			
			Hg%sec(p)%vdat(jj) = dsqrt(dble(f2j * f0i));
			Hb%sec(p)%vdat(jj) = dsqrt(dble(f2j * f1k));

		end do

!	write(*,'(1000i5)') Hb%sec(p)%row
!	write(*,'(1000i5)') Hb%sec(p)%col

!	write(*,'(1000i5)') Hg%sec(p)%row
!	write(*,'(1000i5)') Hg%sec(p)%col


	return
	end 	subroutine MakeHgbBlock
	!----------------------------------------------------------




	!----------------------------------------------------------
	subroutine MakeHgb(n,nph)
	! constructs sparse Hg using Hg%sec
	! saves the Hamiltonians in CSR format
	implicit none
	integer, intent(in) :: n,nph
	double precision:: x,y
	integer :: nnz, i,p,i1,i2,l,q,n0,n0dn,n0up,j1,j2


	l = basis%sec(n)%ntot; ! == number of MOLECULAR states in the block
	q = basis%sec(n-1)%ntot; ! == number of nnz of hgblock
	
	! total nnz
	nnz = 2*q * nph !* 2 for co+counter rotating term;
	! nph instead of nph+1: 
	!for missing q terms in co and missing q in counter in p=0 and p=nph blocks.
	
	Hg%nnz = nnz
	! aux arrays to hold data from all p-blocks in coo format
	if(allocated(Hg%coo1))deallocate(Hg%coo1)
	if(allocated(Hg%coo2))deallocate(Hg%coo2)
	if(allocated(Hg%coodat))deallocate(Hg%coodat)
	allocate(Hg%coo1(nnz))
	allocate(Hg%coo2(nnz))
	allocate(Hg%coodat(nnz))


	! blocks for soc term
	nnz = q * (nph + 1)
	Hb%nnz = nnz
	if(allocated(Hb%coo1))deallocate(Hb%coo1)
	if(allocated(Hb%coo2))deallocate(Hb%coo2)
	if(allocated(Hb%coodat))deallocate(Hb%coodat)
	allocate(Hb%coo1(nnz))
	allocate(Hb%coo2(nnz))
	allocate(Hb%coodat(nnz))


	ntot = l * (nph+1);
	Hg%ntot = ntot

	! coo format but with global indexing 

	i1 = 0;
	! p=0 and nph cases dealt with seperatly
	do p=0,nph ! photons
	 x = dsqrt(dble(p)); ! a^-
	 y = dsqrt(dble(p+1)) * uscfac; ! a^+

	 ! absolute basis index:
	 ! n0up, n0dn:   up/dn here refer to changes in photon number
	 n0 = p*l; 
	 n0up = (p+1)*l; ! p < nph; 
	 n0dn = (p-1)*l; ! p > 0

	 ! counter rotating term: a^+ [total q terms]
	 if(p<nph) then
		i2 = i1 + q;
		Hg%coo1(i1+1:i2) = n0 + Hg%sec(1)%row;
		Hg%coo2(i1+1:i2) = n0up + Hg%sec(1)%col;			
		Hg%coodat(i1+1:i2) = y * Hg%sec(1)%vdat;
		i1 = i2;
	 endif
		
	 ! co rotating term: a^- [total q terms]
	 if(p>0) then
	  i2 = i1 + q; 
		Hg%coo1(i1+1:i2) = n0 + Hg%sec(1)%row;
		Hg%coo2(i1+1:i2) = n0dn + Hg%sec(1)%col;	 ! n0dn: one less photon
		Hg%coodat(i1+1:i2) = x * Hg%sec(1)%vdat;
		i1 = i2;
	 endif
		
	 ! SOC: sigma^+ tau^-
	 j1 = p*q; j2 = j1 + q; 
	 Hb%coo1(j1+1:j2) = n0 + Hb%sec(1)%row;
	 Hb%coo2(j1+1:j2) = n0 + Hb%sec(1)%col;		! n0: same photon number
	 Hb%coodat(j1+1:j2) = Hb%sec(1)%vdat; 

				
	end do

!	write(*,'(1000i5)') Hb%coo1
!	write(*,'(1000i5)') Hb%coo2

	return
	end 	subroutine MakeHgb
	!----------------------------------------------------------













	!----------------------------------------------------------
	! two fission sites space 3x3= 9 states
	! 9x9 blocks made up in the space of the rest of the system (cavity + symmtric molecules)
	! local indexing for the H1r and Hb1 blocks.
	! constructs sparse H1r block for one of the fission sites
	! size of the block is equal to the size of the symmetric space & cavity.
	subroutine MakeH2(n,nph)
	implicit none
	integer, intent(in) :: n,nph
	double precision:: x,y
	integer :: nnz, i,p,i1,i2,l,q,n0,n0dn,n0up,j1,j2
	integer :: ntotsym
	integer, allocatable, dimension(:) :: lrange


	ntotsym = Hg%ntot
	l = basis%sec(n)%ntot; ! ntotb
	q = l; ! ntotb

	
	! indices of symmetric molecular block with ntotb states
	allocate(lrange(l))
	lrange = (/ (i, i= 1, l) /); ! ntotb_range

	
	! total nnz
	nnz = 2* l * nph !* 2 for co+counter rotating term;
	! nph instead of nph+1: 
	!for missing q terms in co and missing q in counter in p=0 and p=nph blocks.
	
	Hg1%nnz = nnz
	! aux arrays to hold data from all p-blocks in coo format
	if(allocated(Hg1%coo1))deallocate(Hg1%coo1)
	if(allocated(Hg1%coo2))deallocate(Hg1%coo2)
	if(allocated(Hg1%coodat))deallocate(Hg1%coodat)
	allocate(Hg1%coo1(nnz))
	allocate(Hg1%coo2(nnz))
	allocate(Hg1%coodat(nnz))


	! blocks for soc term
	Hb1%nnz = ntotsym
	if(allocated(Hb1%coo1))deallocate(Hb1%coo1)
	if(allocated(Hb1%coo2))deallocate(Hb1%coo2)
	if(allocated(Hb1%coodat))deallocate(Hb1%coodat)
	allocate(Hb1%coo1(nnz))
	allocate(Hb1%coo2(nnz))
	allocate(Hb1%coodat(nnz))

	! coo format, local indexing

	i1 = 0;
	! p=0 and nph cases dealt with seperatly
	do p=0,nph ! photons
	 x = dsqrt(dble(p)); ! a^-
	 y = dsqrt(dble(p+1)) * uscfac; ! a^+

	 ! absolute basis index:
	 ! n0up, n0dn:   up/dn here refer to changes in photon number
	 n0 = p*l; 
	 n0up = (p+1)*l; ! p < nph; 
	 n0dn = (p-1)*l; ! p > 0

	 ! counter rotating term: a^+ [total q terms]
	 if(p<nph) then
		i2 = i1 + q;
		Hg1%coo1(i1+1:i2) = n0 + lrange;
		Hg1%coo2(i1+1:i2) = n0up + lrange;			
		Hg1%coodat(i1+1:i2) = y;
		i1 = i2;
	 endif
		
	 ! co rotating term: a^- [total q terms]
	 if(p>0) then
	  i2 = i1 + q; 
		Hg1%coo1(i1+1:i2) = n0 + lrange;
		Hg1%coo2(i1+1:i2) = n0dn + lrange;	 ! n0dn: one less photon
		Hg1%coodat(i1+1:i2) = x;
		i1 = i2;
	 endif
		
	end do ! p

	! SOC: sigma^+ tau^-
	! diagonal block of size Hg%ntot
	Hb1%coo1(1:ntotsym) = (/ (i, i=1,ntotsym)/);
	Hb1%coo2(1:ntotsym) = (/ (i, i=1,ntotsym)/);
	Hb1%coodat(1:ntotsym) = 1.0d0; 




	! using Hg1 and Hb1 blocks, build full space Hg1f and Hb1f
	Hg1f%nnz = 6*Hg1%nnz; ! 3 blocks for first mol (each for a diagonal block of the second) and 3 blocks for second mol, etc.
	if(allocated(Hg1f%coo1))deallocate(Hg1f%coo1)
	if(allocated(Hg1f%coo2))deallocate(Hg1f%coo2)
	if(allocated(Hg1f%coodat))deallocate(Hg1f%coodat)
	allocate(Hg1f%coo1(Hg1f%nnz))
	allocate(Hg1f%coo2(Hg1f%nnz))
	allocate(Hg1f%coodat(Hg1f%nnz))

	Hb1f%nnz = 6*Hb1%nnz;
	if(allocated(Hb1f%coo1))deallocate(Hb1f%coo1)
	if(allocated(Hb1f%coo2))deallocate(Hb1f%coo2)
	if(allocated(Hb1f%coodat))deallocate(Hb1f%coodat)
	allocate(Hb1f%coo1(Hb1f%nnz))
	allocate(Hb1f%coo2(Hb1f%nnz))
	allocate(Hb1f%coodat(Hb1f%nnz))



	p1 = 0; q1=0;
	do i1=1,3
	 k1 = (i1-1)*ntotsym
	 !l1 = k1 + ntotsym
	 do j1=1,3
	  m1 = (j1-1)*ntotsym
	  !n1 = m1 + ntotsym
	  do i2=1,3
	   k2 = (i2-1)*ntotsym
	   !l2 = k2 + ntotsym
	   do j2=1,3
	    m2 = (j2-1)*ntotsym
	    !n2 = m2 + ntotsym

			! molecule 1: |S1><G1| & diagonal in mol 2
			if(i1==3 .and. j1==1 .and. i2==j2) then 
			 p2 = p1 + ntotsym;
	     Hg1f%coo1(p1+1:p2) = k1 + Hg1%coo1
	     Hg1f%coo2(p1+1:p2) = m1 + Hg1%coo2
	     Hg1f%coodat(p1+1:p2) = Hg1%coodat
	     p1 = p2;
			endif

			! molecule 2: |S2><G2| & diagonal in mol 1
			if(i2==3 .and. j2==1 .and. i1==j1) then 
			 p2 = p1 + ntotsym;
	     Hg1f%coo1(p1+1:p2) = k2 + Hg1%coo1
	     Hg1f%coo2(p1+1:p2) = m2 + Hg1%coo2
	     Hg1f%coodat(p1+1:p2) = Hg1%coodat
	     p1 = p2;
			endif


			! molecule 1: |S1><T1| & diagonal in mol 2
			if(i1==2 .and. j1==1 .and. i2==j2) then 
			 q2 = q1 + ntotsym;
	     Hg1f%coo1(q1+1:q2) = k1 + Hg1%coo1
	     Hg1f%coo2(q1+1:q2) = m1 + Hg1%coo2
	     Hg1f%coodat(q1+1:q2) = Hg1%coodat
	     q1 = q2;
			endif

			! molecule 2: |S2><T2| & diagonal in mol 1
			if(i2==2 .and. j2==1 .and. i1==j1) then 
			 q2 = q1 + ntotsym;
	     Hg1f%coo1(q1+1:q2) = k2 + Hg1%coo1
	     Hg1f%coo2(q1+1:q2) = m2 + Hg1%coo2
	     Hg1f%coodat(q1+1:q2) = Hg1%coodat
	     q1 = q2;
			endif



	   end do
	  end do
	 end do
	end do








	return
	end 	subroutine MakeH2
	!----------------------------------------------------------






	!----------------------------------------------------------
	! diagonal parts of cavity-molecule hamiltonian:
	! Hc for cavity, Hs for singlet molecules, Ht for triplet molecules.
	subroutine MakeHd(n,nph) 
	implicit none
	integer, intent(in) :: n, nph

	integer :: i,i1,i2,p !, ntot,ntotb

	!ntot = Hg%ntot; ! set Hg%ntot first
	! ntot already set in MakeHgb alongwith Hg%ntot

	if(allocated(Hc)) deallocate(Hc)
	allocate(Hc(ntot)) 

	if(allocated(Hs)) deallocate(Hs)
	allocate(Hs(ntot))  

	if(allocated(Ht)) deallocate(Ht)
	allocate(Ht(ntot)) 


	! molecular block's ntot [the same for all photon numbers]
	ntotb = basis%sec(n)%ntot; ! == number of MOLECULAR states in the block

	!write(*,'(a,2x,10000i5)') "f(2,:) = ", basis%sec(n)%f(2,:)
	!write(*,'(a,2x,10000i5)') "f(1,:) = ", basis%sec(n)%f(1,:)

	! coo format but with global indexing 
	i1 = 0
	do p=0,nph
	  i2 = ntotb * (p+1) 
	  
		Hc(i1+1:i2) = dble(p)
		Hs(i1+1:i2) = (/(basis%sec(n)%f(2,i), i=1,ntotb) /)
		Ht(i1+1:i2) = (/ (basis%sec(n)%f(1,i), i=1,ntotb) /)

		i1 = i2;
	end do

	return
	end 	subroutine MakeHd
	!----------------------------------------------------------







	subroutine MakeHhtc(n,ijob)!, wr,delta,lambda,wv)
	implicit none
	integer, intent(in) :: n,ijob
	!double precision, intent(in) :: wr, delta, lambda, wv

	! local
	double precision :: g, lamwv, w0, wcc, wv
	integer :: i,n1,n2,n3,nnz, Hbnnzm

	! parameters for this job
	wr = param(ijob)%wr; ! omega_R
	delta = param(ijob)%del; ! detuning = w0-w
	lambda = param(ijob)%lam; ! for lambda_SOC

	g = wr/dsqrt(dble(n));
	lamwv = lambda;

	! half values needed for diagonal terms:  
	! *0.5 here for efficieny
	! use wc instead of param(ijob)%wc
	wcc = (wc)*0.5d0; ! Cavity
	w0 = (wc+delta)*0.5d0; !Signlet
	wv = (wc-param(ijob)%j)*0.5d0; ! triplet

	!write(*,*)"wc, w0, wt = ", 2*wcc, 2*w0, 2*wv

	! Hg has upper (a^+, counter-rotating terms) triangular
	!      AND lower (a^-, co-rotating terms) triangular elements.
	! no issues with having mixed upper/lower elements as long as they appear once in either upper or lower triangular.
	! diagonal: Hdv are halved for matvec()

	!..............................................
	! combine all terms to make full Hamiltnian in coo
	!..............................................



	nnz = Hg%nnz + Hb%nnz + Hg%ntot; ! +Hg%ntot for size of diagonal term, Hv+Hd 


!	write(*,*) "Hg%nnz, Hb%nnz, Hg%ntot, Hhtc%nnz = ", 
!     . Hg%nnz, Hb%nnz, Hg%ntot, nnz
	
	Hhtc%nnz = nnz;
	if(allocated(Hhtc%coo1))deallocate(Hhtc%coo1)
	if(allocated(Hhtc%coo2))deallocate(Hhtc%coo2)
	if(allocated(Hhtc%coodat))deallocate(Hhtc%coodat)
	allocate(Hhtc%coo1(nnz))
	allocate(Hhtc%coo2(nnz))
	allocate(Hhtc%coodat(nnz))

	if(allocated(Hdv)) deallocate(Hdv)
	allocate(Hdv(Hg%ntot))  ! set Hg%ntot first

	Hdv = wcc*Hc + w0*Hs +  wv*Ht; ! summ arrays, values.
	! diagonal will be kept seperate

!	write(*,'(1000i5)') Hg%coo1


	! Hg first
	n1 = 0;
	n3 = Hg%nnz;
	n2 = n1 + n3;
!	write(*,*)"Hg: n1,n2,n3 = ", n1,n2,n3
	Hhtc%coo1(n1+1:n2) = Hg%coo1(1:n3)
	Hhtc%coo2(n1+1:n2) = Hg%coo2(1:n3)
	Hhtc%coodat(n1+1:n2) = g * Hg%coodat(1:n3)



	! Hb afterwards
	n1 = n2;
	n3 = Hb%nnz;
	n2 = n1 + n3;
!	write(*,*)"Hb: n1,n2,n3 = ", n1,n2,n3
	Hhtc%coo1(n1+1:n2) = Hb%coo1(1:n3)
	Hhtc%coo2(n1+1:n2) = Hb%coo2(1:n3)
	Hhtc%coodat(n1+1:n2) = lamwv * Hb%coodat(1:n3)


	! now diagonal terms
	n1 = n2;
	n3 = Hg%ntot;
	n2 = n1 + n3;
!	write(*,*)"Hd: n1,n2,n3 = ", n1,n2,n3
	Hhtc%coo1(n1+1:n2) = (/ (i,i=1,n3) /)
	Hhtc%coo2(n1+1:n2) = (/ (i,i=1,n3) /) 
	Hhtc%coodat(n1+1:n2) = Hdv(1:n3) ! already halved for matvec(), see wc,w0,wt above.

	!write(*,'(a,3x,10000f15.10)') 'Hdv = ',Hhtc%coodat(n1+1:n2)

	!write(*,*) 'Hdv: n1, n2 = ',n1, n2 


!	write(*,'(1000i5)') Hhtc%coo1
!	write(*,'(1000i5)') Hhtc%coo2


	!write(*,*)'============ htc =========='
	
	Hf%dense = .false.
	Hf%ntot = Hg%ntot ! dimension of full hilbert space
	Hf%nnz = Hhtc%nnz; 
	if (Hf%ntot .le. nmaxddiag .and. ddiagOK) then 
		! Hhtc in dense format and use direct diagonalisation
		Hf%dense = .true.
		n1 =Hf%ntot;
		if(allocated(Hf%h))deallocate(Hf%h)
		allocate(Hf%h(n1,n1))
		call coo2dense(n1, n2, Hhtc%coo1(1:n2),
     .      Hhtc%coo2(1:n2),Hhtc%coodat(1:n2), Hf%h)

	 !write(*,*)'************* H *******************'
	 !do i=1,n1
	 ! write(*,'(1000f6.2)') Hf%h(i,:)
	 !end do
	 !write(*,*)'************************************'

	else ! iterative diaognalisation
		!..............................................
		! change the format to CSR 
		!..............................................
		if(allocated(Hf%col))deallocate(Hf%col)
		if(allocated(Hf%dat))deallocate(Hf%dat)
		if(allocated(Hf%rowpntr))deallocate(Hf%rowpntr)
		allocate(Hf%col(nnz))
		allocate(Hf%dat(nnz))
		allocate(Hf%rowpntr(Hf%ntot + 1))
		call coocsr(Hf%ntot, Hf%nnz, 
     .  Hhtc%coodat, Hhtc%coo1, Hhtc%coo2,  
     .  Hf%dat, Hf%col, Hf%rowpntr)
	endif

	
!	endif


	return
	end subroutine MakeHhtc
	!--------------------------------------------------





	subroutine coo2dense(ntot,nnz,row,col,val,a)
	implicit none
	integer, intent(in) :: ntot,nnz
	integer, dimension(nnz), intent(in)  :: row, col
	double precision, dimension(nnz), intent(in)  :: val
	double precision, dimension(ntot,ntot), intent(out) :: a
	integer :: i,j,k
		a =0.0d0
		do k=1,nnz
			a(row(k),col(k)) = val(k)
		end do
		! symmetrise
		a = a + transpose(a); 
		! diagonal terms were halved so
		! no need to minus diagonalmatrix(diagonal(a))
	return
	end 	subroutine coo2dense

	!--------------------------------------------------
	subroutine coocsr( nrow, nnz, a, ir, jc, ao, jao, iao )
	!----------------------------------------------------------
	! MAZ Oct 2018, Modified from sparskit, @ Youcef Saad, 07 January 2004
!! COOCSR converts COO to CSR.
!    This routine converts a matrix that is stored in COO coordinate format
!    a, ir, jc into a CSR row general sparse ao, jao, iao format.
! IN:
!    Input, integer NROW, the row dimension of the matrix.
!    Input, integer NNZ, the number of nonzero elements.
!
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
!         the elements, ir(k) = its row number and jc(k) = its column
!        number. The order of the elements is arbitrary.
!
! OUT:
!    Output, real AO(*), JAO(*), IAO(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
	implicit none

	! in
	integer, intent(in):: nrow, nnz
	integer, dimension(nnz), intent(in):: ir, jc
	double precision, dimension(nnz), intent(in) :: a
	! out
	integer, dimension(nnz), intent(out):: jao
	integer, dimension(nrow+1), intent(out):: iao
	double precision, dimension(nnz), intent(out) :: ao

	! local
	integer:: i, iad,j,k,k0
	double precision:: x

	iao(1:nrow+1) = 0
!
!  Determine the row lengths.
!
	do k = 1, nnz
		iao(ir(k)) = iao(ir(k)) + 1
	end do
!
!  The starting position of each row.
!
	k = 1
	do j = 1, nrow+1
		k0 = iao(j)
		iao(j) = k
		k = k + k0
	end do
!
!  Go through the structure once more.  Fill in output matrix.
!
	do k = 1, nnz
		i = ir(k)
		j = jc(k)
		x = a(k)
		iad = iao(i)
		ao(iad) = x
		jao(iad) = j
		iao(i) = iad + 1
	end do
!
!  Shift back IAO.
!
	do j = nrow, 1, -1
		iao(j+1) = iao(j)
	end do
	iao(1) = 1
	
	return
	end subroutine coocsr 
	!----------------------------------------------------------






	subroutine MakeHhtc2(n,ijob)!, wr,delta,lambda,wv)
	implicit none
	integer, intent(in) :: n,ijob
	!double precision, intent(in) :: wr, delta, lambda, wv

	! local
	double precision :: g, lamwv, w0, wcc, wv
	integer :: i,n1,n2,n3,nnz, Hbnnzm

	! parameters for this job
	wr = param(ijob)%wr; ! omega_R
	delta = param(ijob)%del; ! detuning = w0-w
	lambda = param(ijob)%lam; ! for lambda_SOC

	g = wr/dsqrt(dble(n));
	lamwv = lambda;

	! half values needed for diagonal terms:  
	! *0.5 here for efficieny
	! use wc instead of param(ijob)%wc
	wcc = (wc)*0.5d0; ! Cavity
	w0 = (wc+delta)*0.5d0; !Signlet
	wv = (wc-param(ijob)%j)*0.5d0; ! triplet

	!write(*,*)"wc, w0, wt = ", 2*wcc, 2*w0, 2*wv

	! Hg has upper (a^+, counter-rotating terms) triangular
	!      AND lower (a^-, co-rotating terms) triangular elements.
	! no issues with having mixed upper/lower elements as long as they appear once in either upper or lower triangular.
	! diagonal: Hdv are halved for matvec()

	!..............................................
	! combine all terms to make full Hamiltnian in coo
	!..............................................



	nnz = Hg%nnz + Hb%nnz + Hg%ntot; ! +Hg%ntot for size of diagonal term, Hv+Hd 


!	write(*,*) "Hg%nnz, Hb%nnz, Hg%ntot, Hhtc%nnz = ", 
!     . Hg%nnz, Hb%nnz, Hg%ntot, nnz
	
	Hhtc%nnz = nnz;
	if(allocated(Hhtc%coo1))deallocate(Hhtc%coo1)
	if(allocated(Hhtc%coo2))deallocate(Hhtc%coo2)
	if(allocated(Hhtc%coodat))deallocate(Hhtc%coodat)
	allocate(Hhtc%coo1(nnz))
	allocate(Hhtc%coo2(nnz))
	allocate(Hhtc%coodat(nnz))

	if(allocated(Hdv)) deallocate(Hdv)
	allocate(Hdv(Hg%ntot))  ! set Hg%ntot first

	Hdv = wcc*Hc + w0*Hs +  wv*Ht; ! summ arrays, values.
	! diagonal will be kept seperate

!	write(*,'(1000i5)') Hg%coo1


	! Hg first
	n1 = 0;
	n3 = Hg%nnz;
	n2 = n1 + n3;
!	write(*,*)"Hg: n1,n2,n3 = ", n1,n2,n3
	Hhtc%coo1(n1+1:n2) = Hg%coo1(1:n3)
	Hhtc%coo2(n1+1:n2) = Hg%coo2(1:n3)
	Hhtc%coodat(n1+1:n2) = g * Hg%coodat(1:n3)



	! Hb afterwards
	n1 = n2;
	n3 = Hb%nnz;
	n2 = n1 + n3;
!	write(*,*)"Hb: n1,n2,n3 = ", n1,n2,n3
	Hhtc%coo1(n1+1:n2) = Hb%coo1(1:n3)
	Hhtc%coo2(n1+1:n2) = Hb%coo2(1:n3)
	Hhtc%coodat(n1+1:n2) = lamwv * Hb%coodat(1:n3)


	! now diagonal terms
	n1 = n2;
	n3 = Hg%ntot;
	n2 = n1 + n3;
!	write(*,*)"Hd: n1,n2,n3 = ", n1,n2,n3
	Hhtc%coo1(n1+1:n2) = (/ (i,i=1,n3) /)
	Hhtc%coo2(n1+1:n2) = (/ (i,i=1,n3) /) 
	Hhtc%coodat(n1+1:n2) = Hdv(1:n3) ! already halved for matvec(), see wc,w0,wt above.

	!write(*,'(a,3x,10000f15.10)') 'Hdv = ',Hhtc%coodat(n1+1:n2)

	!write(*,*) 'Hdv: n1, n2 = ',n1, n2 


!	write(*,'(1000i5)') Hhtc%coo1
!	write(*,'(1000i5)') Hhtc%coo2


	!write(*,*)'============ htc =========='
	
	Hf%dense = .false.
	Hf%ntot = Hg%ntot ! dimension of full hilbert space
	Hf%nnz = Hhtc%nnz; 
	if (Hf%ntot .le. nmaxddiag .and. ddiagOK) then 
		! Hhtc in dense format and use direct diagonalisation
		Hf%dense = .true.
		n1 =Hf%ntot;
		if(allocated(Hf%h))deallocate(Hf%h)
		allocate(Hf%h(n1,n1))
		call coo2dense(n1, n2, Hhtc%coo1(1:n2),
     .      Hhtc%coo2(1:n2),Hhtc%coodat(1:n2), Hf%h)

	 !write(*,*)'************* H *******************'
	 !do i=1,n1
	 ! write(*,'(1000f6.2)') Hf%h(i,:)
	 !end do
	 !write(*,*)'************************************'

	else ! iterative diaognalisation
		!..............................................
		! change the format to CSR 
		!..............................................
		if(allocated(Hf%col))deallocate(Hf%col)
		if(allocated(Hf%dat))deallocate(Hf%dat)
		if(allocated(Hf%rowpntr))deallocate(Hf%rowpntr)
		allocate(Hf%col(nnz))
		allocate(Hf%dat(nnz))
		allocate(Hf%rowpntr(Hf%ntot + 1))
		call coocsr(Hf%ntot, Hf%nnz, 
     .  Hhtc%coodat, Hhtc%coo1, Hhtc%coo2,  
     .  Hf%dat, Hf%col, Hf%rowpntr)
	endif

	
!	endif


	return
	end subroutine MakeHhtc2
	!--------------------------------------------------



	end !module


	


	module hamiltonian
	use modmain	
	implicit none
	
	public :: MakeHhtc, HamParts
	private:: MakeHgSec,MakeHg
	private:: MakeHd,MakeHv
	private:: MakeHbSec,MakeHb
	private:: sumdup, rowadddup
	public :: coocsr, coo2dense
	contains
 
! ===========================================
	subroutine HamParts(n) ! arg n is norig here, set in main.f
	implicit none
	integer, intent(in):: n

	!write(*,*)'------------------1 '
	call MakeHgBlock(n)
	!write(*,*)'------------------2 '
	if(.not. writewfs) then
		call MakeHg(n,m)
	else
		call MakeHgAndVx(n,m)
	endif
	!write(*,*)'------------------3 '
	call MakeHd(n,m)
	!write(*,*)'------------------4 '
	call MakeHv(n,m)
	!write(*,*)'------------------5 '
	if(mv > 0) then
		call MakeHbSec(n,m,mv)
		!write(*,*)'------------------6 '
		call MakeHb(n,m,mv)
	else
		Hb%nnz = 0; !use to calc size of sparse Hhtc in MakeHhtc
	endif
	!write(*,*)'------------------7 '
	return
	end subroutine HamParts
! ===========================================
!	subroutine MakeHam(n,m,mv)
!	implicit none
!	integer, intent(in):: n,m,mv
!
!	call MakeHhtc(n, wr,delta,lambda,wv)
!
!	subroutine
!	end 	subroutine MakeHam
! ===========================================

	!----------------------------------------------------------
	subroutine MakeHgbBlock(n)
	implicit none
	integer, intent(in) :: n
	integer :: p,nnz,i,j,jj
	double precision :: f2j,f0i,f1k
	integer :: ntotp, ntotp1

	p = 1; ! for the block for n sites
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

			
			Hg%sec(p)%row(jj) = i
			Hg%sec(p)%col(jj) = j; ! sigma^+ a^-

			Hb%sec(p)%row(jj) = k
			Hb%sec(p)%col(jj) = j; ! sigma^+ tau^-

			
			f0i = basis%sec(n)%f(0,i); 
			f2j = basis%sec(n)%f(2,j);
			f1k = basis%sec(n)%f(1,k); 
			
			Hg%sec(p)%vdat(jj) = dsqrt(dble(f2j * f0i));
			Hb%sec(p)%vdat(jj) = dsqrt(dble(f2j * f1k));

		end do

	return
	end 	subroutine MakeHgbBlock
	!----------------------------------------------------------




	!----------------------------------------------------------
	subroutine MakeHg(n,nph)
	! constructs sparse Hg using Hg%sec
	! saves the Hamiltonians in CSR format
	implicit none
	integer, intent(in) :: n,nph
	double precision:: x,y
	integer :: nnz,ntot,i,p,i1,i2,l,q,n0,n0dn,n0up


	l = basis%sec(n)%ntot; ! == number of states in the block
	q = basis%sec(n-1)%ntot; ! == number of nnz of hgblock
	
	! total nnz
	nnz = 2*q * nph !* 2 for co+counter rotating term
	Hg%nnz = nnz
	! aux arrays to hold data from all p-blocks in coo format
	if(allocated(Hg%coo1))deallocate(Hg%coo1)
	if(allocated(Hg%coo2))deallocate(Hg%coo2)
	if(allocated(Hg%coodat))deallocate(Hg%coodat)
	allocate(Hg%coo1(nnz))
	allocate(Hg%coo2(nnz))
	allocate(Hg%coodat(nnz))

	ntot = l * nph;
	Hg%ntot = ntot

	! coo format but with global indexing 

	! p=0 and nph cases dealt with seperatly
	do p=0,nph ! photons
		x = dsqrt(dble(p)); ! a^-
		y = dsqrt(dble(p+1)); ! a^+

		! absolute basis index:
		! n0up, n0dn:   up/dn here refer to changes in photon number
		n0 = p*l; 
		n0up = (p+1)*l; ! p < nph; 
	  n0dn = (p-1)*l; ! p > 0

		! counter rotating term: a^+
		i1 = p*2*q; i2 = i1 + q; ! *2*q because 2q terms in co+counter.
		Hg%coo1(i1+1:i2) = n0 + Hg%sec(p)%row;
		Hg%coo2(i1+1:i2) = n0up + Hg%sec(p)%col;			
		Hg%coodat(i1+1:i2) = y * Hg%sec(p)%vdat; 

		! co rotating term: a^-
		i1 = i2; i2 = i2 + q; 
		Hg%coo1(i1+1:i2) = n0 + Hg%sec(p)%row;
		Hg%coo2(i1+1:i2) = n0dn + Hg%sec(p)%col;	 ! n0dn: one less photon
		Hg%coodat(i1+1:i2) = x * Hg%sec(p)%vdat; 

		! SOC: sigma^+ tau^-
		Hb%coo1(i1+1:i2) = n0 + Hb%sec(p)%row;
		Hb%coo2(i1+1:i2) = n0 + Hb%sec(p)%col;		! n0: same photon number
		Hb%coodat(i1+1:i2) = Hb%sec(p)%vdat; 
				
	end do

	return
	end 	subroutine MakeHg
	!----------------------------------------------------------








	!----------------------------------------------------------
	subroutine MakeHd(n,m) ! diagonal part of cavity-exciton hamiltonian
	! constructs sparse Hd for m excitations using Hg%sec
	implicit none
	integer, intent(in) :: n,m
	integer :: i1,m1,ntot,p

	if(allocated(Hd)) deallocate(Hd)
	allocate(Hd(Hg%ntot))  ! set Hg%ntot first
	!write(*,*) 'Hg%ntot = ', Hg%ntot 
	m1=min(m,n);
	! coo format but with global indexing 
	i1 = 0;
	do p=0,m1
		ntot = basis%sec(p)%ntot * basis%sec(n-p)%ntot ! N_{nu_p} * N_{mu_p}
		!write(*,*) 'p, pth_ntot = ', p, ntot
		Hd(i1+1:i1+ntot) = (m-p)*1.0d0
		i1 = i1 + ntot;
	end do


	return
	end 	subroutine MakeHd
	!----------------------------------------------------------



	!----------------------------------------------------------
	subroutine MakeHv(n,m) ! diagonal part of vibrational hamiltonian
	! constructs sparse Hv using basis%sec
	implicit none
	integer, intent(in) :: n,m
	integer :: i1,m1,ntot,p, Nvi, i2,i

	if(allocated(Hv)) deallocate(Hv)
	allocate(Hv(Hg%ntot)) ! set Hg%ntot first

	m1 = min(n,m);
	! coo format but with global indexing 
	i1 = 0;
	do p=0,m1
		!ntot = basis%sec(p)%ntot * basis%sec(n-p)%ntot ! N_{nu_p} * N_{mu_p}
		! vib on p excited and n-p unexcited molecules
		do i=1,basis%sec(p)%ntot
			!write(*,*)'i,p, basis%sec(p)%ntot =',i,p, basis%sec(p)%ntot
			Nvi = basis%sec(p)%Nv(i)
			i2 = i1 + basis%sec(n-p)%ntot
			Hv(i1+1:i2) = Nvi + basis%sec(n-p)%Nv(:)
			i1 = i2; 
		end do
	end do

	return
	end 	subroutine MakeHv
	!----------------------------------------------------------

	



	!----------------------------------------------------------
	subroutine MakeHbSec(n,m,mv) ! m for mv==M
	implicit none
	integer, intent(in) ::n,m, mv
	integer :: fi, fi1
	integer :: p, nnz,nnp,ind,k,i,j,ii !,indi,indf 
	double precision :: Pnui,Pnuf
	integer:: ntotp
	
	! basis info: Normalisations of vibrational states
	! define indexes from 0 for sec of Ham and basis so that p sec is p-th elem.
	! similarly, map ind from 0,0
	! vibrational maps that connect perm symm states of p sites with those of p+1 sites... 

	!m1max set in main program

	if(allocated(Hb%sec)) deallocate(Hb%sec)
	allocate(Hb%sec(0:m1max))
	!allocate(Hg%nnzs(m1max))	! array for number of nnz in all blocks
	!write(6,*)'Hbsec: m1max=',m1max
	! loop over exciton-photon blocks.
	! p=0, no hb
	do p=1,m1max ! m1=Min(N,m); m1max= max if calc for multiple m's are desired.
		! nnz = N_{p-1}*M*N_{n-p} non-zero elements in pth block of Hb: some may have the same coordinates, will be combined before converting to csr.
		!nnz=(basis%sec(p-1)%ntot)*mv !*(basis%sec(n-p)%ntot) ! N_{n-p} states with down spins.; Hb is diagonal in these. treated later.
		ntotp = min(basis%sec(p-1)%ntot, basis%sec(ndummy-1)%ntot);
		nnz= ntotp * mv;
		Hb%sec(p)%nnz = nnz
		allocate(Hb%sec(p)%row(nnz))
		allocate(Hb%sec(p)%col(nnz))
		allocate(Hb%sec(p)%vdat(nnz))

		nnp = basis%sec(n-p)%ntot; ! number of perm symm basis for n-p sites [for mu, spin down n-p sies]

		Hb%sec(p)%ntot = basis%sec(p)%ntot ! nu_p basis only, mu_{n-p} treated later

		ind=0;
		do ii = 1, ntotp !basis%sec(p-1)%ntot
		do k=0,mv-1 
			ind = ind + 1; 
			i = map(k,ii);
			j = map(k+1,ii);	

			Hb%sec(p)%row(ind) = i !(l,l=(i-1)*nnp+1,i*nnp)
			Hb%sec(p)%col(ind) = j !(l,l=(j-1)*nnp+1,j*nnp)

			! number of total perm that basis states i, j stand for.
			!Pnui=basis%sec(p)%P(i);
			!Pnuf=basis%sec(p)%P(j);
			! frequency of k in i-th perm symm state, fi = ?
			fi = basis%sec(p)%f(k,i)	 ! from basis .... tabulated
			fi1 = basis%sec(p)%f(k+1,i)	

			!Hb%sec(p)%vdat(ind) = fi*dsqrt((k+1)*Pnui/Pnuf) ! ahsan, 14-01-2020
			!Hb%sec(p)%vdat(ind) = fi * dsqrt((k+1)*1.0d0/(fi+1))
			Hb%sec(p)%vdat(ind) = dsqrt((k+1)*1.0d0*fi*(fi1+1))

			! for details, see notepad starting 13 jan 2020 on chi2 code...
			! Pnui/Pnuf = 1/(fi+1)
		end do ! k 
		end do ! ii
	end do ! p

	return
	end 	subroutine MakeHbSec
	!----------------------------------------------------------







	!----------------------------------------------------------
	subroutine MakeHb(n,m,mv) !,spref)
	! constructs sparse Hb for m excitations using Hb%sec
	! NOTE: only first spref elements of Hb%coo1 , Hb%coo12, Hb%coodat arrays are significant.
	implicit none
	integer, intent(in) :: n,m,mv
	!integer, intent(out) :: spref
	integer :: spref
	integer, dimension(:), allocatable :: iro, jco
	double precision, dimension(:), allocatable :: aoo
	integer :: i1,i2,l1,l2,j1,j2,ii,k,i,j,m1,maxsize
	integer :: nnp,nnp1,nnz,nnzu,p,nrow,ref
	double precision :: val
	
	! find max size of sector to allocate aux arrays
	maxsize = 0; nnz = 0;
	m1=min(m,n);
	do p=1,m1 ! exclude p=0 for Hb 
		!write(*,*) 'p, Hb%sec(p)%nnz = ',p,Hb%sec(p)%nnz
		maxsize=max(maxsize,Hb%sec(p)%nnz) ! maxsize in nu_p space
		nnp = basis%sec(n-p)%ntot ! N_{N-p}; (for unexcited sites)
		nnz = nnz + Hb%sec(p)%nnz * nnp ! nnz elem in full [nu x mu] space
	enddo




	!write(*,*)'maxsize = ',maxsize
	
	! aux arrays:
	allocate(iro(maxsize))
	allocate(jco(maxsize))
	allocate(aoo(maxsize))
	
	! full Hb coo arrays: in FULL [nu x mu] space.
	if(allocated(Hb%coo1))deallocate(Hb%coo1)
	if(allocated(Hb%coo2))deallocate(Hb%coo2)
	if(allocated(Hb%coodat))deallocate(Hb%coodat)
	allocate(Hb%coo1(nnz))
	allocate(Hb%coo2(nnz))
	allocate(Hb%coodat(nnz))
	if(allocated(Hb%spntr))deallocate(Hb%spntr)
	allocate(Hb%spntr(m1+1))
	Hb%spntr(1) = 1;

	! total nnz
	!m1=min(m,n);
	ref = basis%sec(n-0)%ntot !N_{p=0}; all photons, no exciton block: binomial(n+m,m) states.
	spref = 0;
	do p=1,m1
		! sorted coo format & sum over duplicate entries.
		! subroutine sumdup( nrow, nnz, ir, jc, a, nnzu, iro, jo, aoo )
		nrow = basis%sec(p)%ntot ! not full Hg%sec(p)%ntot
		nnz = Hb%sec(p)%nnz;
		!write(*,*) 'p =', p, 'nrow = ',nrow, 'nnz = ', nnz
		call sumdup(nrow, nnz, 
     .          Hb%sec(p)%row, Hb%sec(p)%col, Hb%sec(p)%vdat,
     .          nnzu, iro(1:nnz), jco(1:nnz), aoo(1:nnz))
		!...........................................................
		! get the global indices by taking into account diagonal mu states...
		! and combine data from all sectors into a larger coo sparse...
		!...........................................................
		! local indexes in p-th block: diagonal in mu
		!	I_{inu,imu} = (inu-1)*N_{N-p} + imu
		!	J_{jnu,imu} = (jnu-1)*N_{N-p} + imu
		nnp = basis%sec(n-p)%ntot ! N_{N-p}; number of perm sym basis for n-p sites; (for unexcited sites)
		l1 = 0; i1=0; j1=0;
		do ii=1,nnzu ! only nnzu elements of iro,jco,aoo are significant
			! diagonal in mu_{n-p} basis; nnp total states of type mu			
			i = iro(ii);
			j = jco(ii);
			val = aoo(ii);
			!..........................................
			i1 = ref + (i-1)*nnp; ! row indices global
			i2 = i1 + nnp;
			!..........................................	
			j1 = ref + (j-1)*nnp;! col indices global
			j2 = j1 + nnp;
			!..........................................			
			l1 = spref + (ii-1)*nnp; ! sparse coo indices
			l2 = l1 + nnp;
			!..........................................			 
			Hb%coo1(l1+1:l2) = (/(k,k=i1+1,i2)/);
			Hb%coo2(l1+1:l2) = (/(k,k=j1+1,j2)/);				
			Hb%coodat(l1+1:l2) = val
		end do	
		ref = ref + (basis%sec(p)%ntot)*nnp !nnp=(basis%sec(n-p)%ntot) ! global position of p-th block
		spref = spref + nnzu*nnp; ! sparse coo reference for p-th block
		Hb%spntr(p+1) = spref + 1
	end do

	! use this to construct full hamiltonian, multiply wv, lambda etc there....
	! NOTE: only first spref elements of Hb%coo1 , Hb%coo12, Hb%coodat arrays are significant.
	Hb%nnz = spref ! only significant size of coo Hb.
	! Hb%nnz set to spref[= significant no of elements]
	! Do not reset it to sum of Hb%sec(p)%nnz
	return
	end 	subroutine MakeHb
	!----------------------------------------------------------

	subroutine MakeHhtc(n,ijob,mode)!, wr,delta,lambda,wv)
	implicit none
	integer, intent(in) :: n,ijob,mode
	!double precision, intent(in) :: wr, delta, lambda, wv

	! local
	double precision :: g, lamwv
	integer :: i,n1,n2,n3,nnz, Hbnnzm

	! parameters for this job
	wr = param(ijob)%wr;
	delta = param(ijob)%del;
	lambda = param(ijob)%lam;
	wv = param(ijob)%wv;

	! hg : upper triangular
	! Hb might have lower triangular elements?! 
	! no issues with having mixed upper/lower elements as long as they appear once in either upper or lower.
	! diagonal: Hdv are halved for matvec()

	!..............................................
	! combine all terms to make full Hamiltnian in coo
	!..............................................

	g = wr/dsqrt(nact*1.0d0);
	lamwv = -lambda*wv;

	nnz = Hg%nnz + Hb%nnz ; !+ Hg%ntot ! Hg%ntot for size of diagonal term, Hv+Hd
	if(mode==1) nnz = nnz + Hg%ntot; ! Hg%ntot for size of diagonal term, Hv+Hd 

	!write(*,*)'Hg%nnz, Hb%nnz, Hg%ntot =',Hg%nnz,Hb%nnz,Hg%ntot
	
	Hhtc%nnz = nnz;
	if(allocated(Hhtc%coo1))deallocate(Hhtc%coo1)
	if(allocated(Hhtc%coo2))deallocate(Hhtc%coo2)
	if(allocated(Hhtc%coodat))deallocate(Hhtc%coodat)
	allocate(Hhtc%coo1(nnz))
	allocate(Hhtc%coo2(nnz))
	allocate(Hhtc%coodat(nnz))

	if(allocated(Hdv)) deallocate(Hdv)
	allocate(Hdv(Hg%ntot))  ! set Hg%ntot first

	!write(*,*)'Hg%ntot = ',Hg%ntot

	Hdv = delta*Hd +  wv*Hv; ! summ arrays, values.
	! diagonal will be kept seperate

	n1 = 0;
	n3 = Hg%nnz;
	n2 = n3;

	!write(*,*) '-----   Hg%nnz = ',Hg%nnz

	! Hg first
	Hhtc%coo1(n1+1:n2) = Hg%coo1(1:n3)
	Hhtc%coo2(n1+1:n2) = Hg%coo2(1:n3)
	Hhtc%coodat(n1+1:n2) = g * Hg%coodat(1:n3)

	if (mv > 0) then
		!write(*,*) 'Hg: n1, n2 = ',n1, n2 
		n1 = n2;
		n3 = Hb%nnz; ! Hb%nnz set to spref[= significant no of elements]
								! Do not reset it to sum of Hb%sec(p)%nnz
		n2 = n1 + n3;
		! Hb afterwards
		Hhtc%coo1(n1+1:n2) = Hb%coo1(1:n3)
		Hhtc%coo2(n1+1:n2) = Hb%coo2(1:n3)
		Hhtc%coodat(n1+1:n2) = lamwv * Hb%coodat(1:n3)
	endif

	if(mode==1) then ! make Hf sparse here....

	n1 = n2;
	n3 = Hg%ntot;
	n2 = n1 + n3;
	! now diagonal terms
	Hhtc%coo1(n1+1:n2) = (/ (i,i=1,n3) /)
	Hhtc%coo2(n1+1:n2) = (/ (i,i=1,n3) /) 
	Hhtc%coodat(n1+1:n2) = Hdv(1:n3)*0.5d0 ! half it for matvec() using upper triangular only

	!write(*,'(a,3x,10000f15.10)') 'Hdv = ',Hhtc%coodat(n1+1:n2)

	!write(*,*) 'Hdv: n1, n2 = ',n1, n2 

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

	
	endif


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








	!--------------------------------------------------
	subroutine sumdup( nrow, nnz, ir, jc, a, nnzu, iro, jo, aoo )
	!----------------------------------------------------------
	! MAZ Oct 2018, sum duplicate entries in coo format.
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
	!    Output, iro, jo, aoo unique entries in  COO
	!
	implicit none

	! in
	integer, intent(in):: nrow, nnz
	integer, dimension(nnz), intent(in):: ir, jc
	double precision, dimension(nnz), intent(in) :: a
	! out
	integer, intent(out) :: nnzu ! nnz-unique
	integer, dimension(nnz), intent(out):: jo, iro ! col and row indexes; only first nnzout significant
	double precision, dimension(nnz), intent(out) :: aoo ! values; only first nnzout significant

	! local
	integer i, iad,j,k,k0, nnzout, i1,i2, ind1, ind2,nnzr
	double precision:: x
	integer, dimension(nrow+1):: iao	
	integer, dimension(nnz):: jao
	double precision, dimension(nnz) :: ao
	double precision, dimension(:), allocatable :: valout
	integer, dimension(:), allocatable :: colout


	!write(*,*)'0    ir = ',ir
	!write(*,*)'0    jc = ',jc
	! ........................................
	!			sort w.r.t row index
	! ........................................
	iao(1:nrow+1) = 0
	!  Determine the row lengths.
	do k = 1, nnz
		iao(ir(k)) = iao(ir(k)) + 1
	end do

	!write(*,*)'0    iao = ',iao


	!  The starting position of each row.
	k = 1
	do j = 1, nrow+1
		k0 = iao(j)
		iao(j) = k
		k = k + k0
	end do
	!  Go through the structure once more.  Fill in output matrix.
	do k = 1, nnz
		i = ir(k)
		j = jc(k)
		x = a(k)
		iad = iao(i)
		ao(iad) = x
		iro(iad) = i
		jao(iad) = j
		iao(i) = iad + 1
	end do
	!write(*,*)'rev iao = ',iao

	! iro, jao, ao: sorted w.r.t row, in coo format
	!  Shift back IAO.
	do j = nrow, 1, -1
		iao(j+1) = iao(j)
	end do
	iao(1) = 1
	!write(*,*)'iao = ',iao
	!write(*,*)'jao = ',jao
	
	! ........................................
	!			sum over duplicate terms.
	! ........................................
	ind1 = 0;
	do i=1, nrow
		! range of entries for this row
		nnzr = iao(i+1)-iao(i); ! number of elem in ith row
		!......................................................
		if (nnzr<1) cycle ! go to next iteration; here, its probably last iteration.
		!......................................................		
		i1 = iao(i); ! starting index
		i2 = iao(i+1) - 1; ! end index
		! rowadddup(nnz,col,val, nnzout, colout,valout)
		allocate(colout(nnzr))
		allocate(valout(nnzr))
		!write(*,*) '----- nnzr = ',nnzr

		call rowadddup(nnzr,jao(i1:i2),ao(i1:i2),nnzout,colout,valout)
		ind2 = ind1 + nnzout;
		iro(ind1+1:ind2) = i; ! row index final
		!iaoo(i+1) = iaoo(i) + nnzout; ! rowpntr final
		jo(ind1+1:ind2) = colout(1:nnzout) ! col index final
		aoo(ind1+1:ind2) = valout(1:nnzout) ! value final
		ind1 = ind2; ! advance index for next iteration
		deallocate(colout,valout)
	enddo
	nnzu = ind2; ! total number of unique nnz values

	return
	end subroutine sumdup
!----------------------------------------------------------------------


!----------------------------------------------------------------------
	subroutine rowadddup(nnz,col,val, nnzout, colout,valout)
	! adds duplicate terms for a given row.
	implicit none
	integer, intent(in):: nnz ! nnz in this row === nnzrow 
	integer, dimension(nnz), intent(in):: col ! col index, with possible repetion
	double precision, dimension(nnz), intent(in):: val ! values
	integer, intent(out):: nnzout ! number of elements with distinct columns indices
	integer, dimension(nnz), intent(out):: colout ! distict col indeces 
	double precision, dimension(nnz), intent(out):: valout ! sum of duplicate values
	! local
	integer:: l,nset,k,j
	logical:: found

	!write(*,*) '----- nnz = ',nnz
	colout=0;
	valout=0.0d0;

	nset=0;
	! set the first element
	nset = nset + 1; 	! set the number of distinct col values set so far.
	colout(nset) = col(nset);
	valout(nset) = val(nset);

	found = .false.
	do l=2,nnz ! 
		j = col(l);
		do k=1,nset
			if(j==colout(k)) then ! found
				valout(k) = valout(k) + val(l) ! add l-th value to this box
				found = .true. 
				exit ! exit the inner do loop, over k.
			endif
		enddo ! k		
		if (.not. found) then
			nset = nset + 1; ! advance index to make space for new distict j value
			colout(nset) = j;
			valout(nset) = val(l);				
		endif
		found = .false. ! reset found for next iteration
	end do ! l

	!set total non-zero elements in out arrays
	nnzout = nset
	return
	end subroutine rowadddup
!----------------------------------------------------------------------


	end !module


	

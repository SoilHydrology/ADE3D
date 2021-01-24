	program MAIN
! -----------------------------------------------------------------
!	3D FEM
! -----------------------------------------------------------------
	implicit none
	integer Ne, Nn, Nb, Nm, Nd, Ng, Ns, NnNd	! array parameters
	integer	Npol, Nf, Ndof				! array parameters
	integer iout, idbg, ipost, ires, ijmp
	integer lastA, lastL, lastB, lastQc1, lastQd1, lastF1, lastQc2, lastQd2, lastF2, &
		lastQc3, lastQd3, lastF3 
	integer lastY, lastZ
	integer Nt, n
	integer ldw, maxit
	real*8 Kappa, tc, aopt, por
	real*8 time, tmax, dt, dto, timeo
	real*8 theta, theta1, Ctol
	integer ipar(16), ipar0(16)			! bCGstab integer parameters array
	real*8 fpar(16), fpar0(16)			! bCGstab real    parameters array
	real*8, allocatable :: Ae(:,:), Le(:,:), Be(:,:)	! element arrays
	real*8, allocatable :: Qce1(:,:), Qde1(:,:)		! element arrays
	real*8, allocatable :: Qce2(:,:), Qde2(:,:)		! element arrays
	real*8, allocatable :: Qce3(:,:), Qde3(:,:)		! element arrays
	real*8 Ain(0:7, 0:7)				! Legendre polynomial Ain coefficients
	real*8, allocatable:: Pn(:,:), dPn(:,:)		! Legendre polynomial values and derivatives vectors
	integer, allocatable :: rA (:), rL (:), rB (:)	! global  arrays (compact rows)
	integer, allocatable :: cA (:), cL (:), cB (:)	! global  arrays (compact columns)
	integer, allocatable :: rY (:), cY (:)		! global  arrays
	integer, allocatable :: rZ (:), cZ (:)		! global  arrays
	integer, allocatable :: rQc1(:), rQd1(:), rF1 (:)	! global  arrays (compact rows)
	integer, allocatable :: rQc2(:), rQd2(:), rF2 (:)	! global  arrays (compact rows)
	integer, allocatable :: rQc3(:), rQd3(:), rF3 (:)	! global  arrays (compact rows)
	integer, allocatable :: cQc1(:), cQd1(:), cF1 (:)	! global  arrays (compact columns)
	integer, allocatable :: cQc2(:), cQd2(:), cF2 (:)	! global  arrays (compact columns)
	integer, allocatable :: cQc3(:), cQd3(:), cF3 (:)	! global  arrays (compact columns)
	integer, allocatable :: ie(:,:)			! global connectivity array
	integer, allocatable :: epol(:)			! global element polynomial array
	integer, allocatable :: nmat(:,:)		! global nodal materials array
	integer, allocatable :: BCe(:,:), BCi(:)	! BC element and local element
						 	! face numbers
	integer, allocatable :: face(:,:)		! element, local face, node 1 - node 4
	real*8, allocatable :: BCn(:,:)			! BC outwards normal
	real*8, allocatable :: BCvalue(:,:,:)		! BC value (j_bar, q_bar or c_bar)
	character*1, allocatable ::  BCtype(:)		! BC type ('R', 'N' or 'D')
	real*8, allocatable :: x(:,:)			! global coordinates array
	real*8, allocatable :: C   (:,:), T   (:,:)	! global  arrays
	real*8, allocatable :: Cold(:,:), Told(:,:) 	! global  arrays
	real*8, allocatable :: Cit (:,:)		! global  arrays
	real*8, allocatable :: Son (:,:), Soe(:,:,:)	! global  arrays
	real*8, allocatable :: Sn (:,:)			! global  arrays (nodal S at t^n+1)
	real*8, allocatable :: Sn0(:,:)			! global  arrays (nodal S at t^n)
	real*8, allocatable :: phi (:), phio(:)		! global  arrays
	real*8, allocatable :: vA (:), vL (:), vB (:)	! global  arrays (compact values)
	real*8, allocatable :: vY (:), vZ(:)		! global  arrays (compact values)
	real*8, allocatable :: vQc1(:), vQd1(:), vF1(:)	! global  arrays (compact values)
	real*8, allocatable :: vQc2(:), vQd2(:), vF2(:)	! global  arrays (compact values)
	real*8, allocatable :: vQc3(:), vQd3(:), vF3(:)	! global  arrays (compact values)
	real*8, allocatable :: V(:,:), D(:,:,:)		! global  arrays
	real*8, allocatable :: fn(:,:,:)		! element face unit normals
	real*8, allocatable :: Rn(:,:), w(:,:), work(:)	! work arrays
	real*8, allocatable :: xg(:), wg(:)		! Gauss abscissas [-1,+1] and weights
	real*8, allocatable :: J(:,:,:,:,:), Ji(:,:,:,:,:), Jac(:,:,:), Jaci(:,:,:)
							! geometric entities
	real*8, allocatable :: aii(:,:,:,:)		! the minors a11, a22, a33 of Jij
	real*8, allocatable :: Shp(:,:,:,:), dNdr(:,:,:,:,:), Wgt(:,:,:,:)
							! shape and weight functions
	real*8, allocatable :: Shpb(:,:,:,:)		! boundary shape functions
	real*8, allocatable :: Vm(:,:)			! nodal averaged array
	real	 time_begin, time_end

	integer Nr					! ### new parameters ###
	integer, allocatable :: order(:,:)		! ### new parameters ###
	real*8, allocatable ::  Rr(:,:)			! ### new parameters ###
	real*8, allocatable :: I_Omega(:,:)		! ### error parameters ###
	real*8, allocatable :: I_Gamma(:,:,:)		! ### error parameters ###
	real*8, allocatable :: I_Jump(:,:)		! ### error parameters ###
	real*8, allocatable :: R_Omega(:)		! ###     residual on element ###
	real*8, allocatable :: J_Gamma(:)		! ###     jump at faces ###

	integer e, Nc, i, s, it
	integer k, pol
	real*8 Cerr, Cerrmax
	real*8 Dn(0:7)					! Legendre polynomial denominators Dn

	data iout/3/, idbg/2/, ipost/4/, ires/12/, ijmp/13/
	data n/0/, time/0./, Nc/0/

	call CPU_TIME ( time_begin )

! open files
	open(iout, file='fout.txt', status='unknown')
	open(idbg, file='fdbg.txt', status='unknown')
	open(ipost,file='post.msh', status='unknown')
	open(ires, file='resi.txt', status='unknown')
	open(ijmp, file='jump.txt', status='unknown')

	write(iout,*) '**********************************'	! version ID
	write(iout,*) '3D p-FEM v02 of 21/12/20'		! version ID
	write(iout,*) '**********************************'	! version ID

! read parameters
	call READPARAM(iout, idbg, Ne, Nn, Nb, Nm, Nd, Ng, Ns, Npol, Nf)

! allocate arrays
	Ndof = Nn + Ne*(Npol-1)	! # of DOFs
	write(iout,*) 'Ndof = ', Ndof

	ldw  = 8*Ndof	! storage for SPARSKIT BCGSTAB
!	NnNd = (Npol+2)*Ndof	! Npol+2 entries for 1D elements (incl. off diagonals - mainly in A)
!	NnNd = (Npol+9)*Ndof	! Npol+9 entries for 3D elements (incl. off diagonals - mainly in A)
	NnNd = (Npol+29)*Ndof	!!!! Npol+29 entries for 3D elements (incl. off diagonals - mainly in A)
	write(iout,*) 'ldw = ', ldw
	allocate ( ie(Ne,9), x(Nn,3), C(Ndof,Ns), T(Ndof,Ns), Cold(Ndof,Ns), Told(Ndof,Ns) )
	allocate ( Soe(Ne,Ns,8), Son(Ndof,Ns), Sn(Ndof,Ns), Sn0(Ndof,Ns) )
	allocate ( Vm(Nn,3), Cit(Nn,Ns) )
	allocate ( vA (NnNd ), vL (NnNd ), vB(NnNd ), vZ(NnNd ) )
	allocate ( rA (Ndof+1) , rL (Ndof+1) , rB(Ndof+1) , rZ(Ndof+1) )
	allocate ( cA (NnNd) , cL (NnNd) , cB(NnNd) , cZ(NnNd) )
	allocate ( vQc1(NnNd), vQd1(NnNd), vF1(NnNd) )
	allocate ( vQc2(NnNd), vQd2(NnNd), vF2(NnNd) )
	allocate ( vQc3(NnNd), vQd3(NnNd), vF3(NnNd) )
	allocate ( rQc1(Ndof+1), rQd1(Ndof+1), rF1(Ndof+1) )
	allocate ( rQc2(Ndof+1), rQd2(Ndof+1), rF2(Ndof+1) )
	allocate ( rQc3(Ndof+1), rQd3(Ndof+1), rF3(Ndof+1) )
	allocate ( cQc1(NnNd), cQd1(NnNd), cF1(NnNd) )
	allocate ( cQc2(NnNd), cQd2(NnNd), cF2(NnNd) )
	allocate ( cQc3(NnNd), cQd3(NnNd), cF3(NnNd) )
	allocate ( V(Ne,3), D(Ne,3,3), nmat(Nn,0:Nd) )
	allocate ( Rn(Ndof,Ns), w(Ndof,Ns), work(ldw) )
	allocate ( xg(Ng), wg(Ng), J(3,3,0:Ng+1,0:Ng+1,0:Ng+1), Ji(3,3,0:Ng+1,0:Ng+1,0:Ng+1), Jac(0:Ng+1,0:Ng+1,0:Ng+1), &
			Jaci(0:Ng+1,0:Ng+1,0:Ng+1), aii(3,0:Ng+1,0:Ng+1,0:Ng+1) )
	allocate ( Shp(Npol+7,0:Ng+1,0:Ng+1,0:Ng+1), dNdr(Npol+7,3,0:Ng+1,0:Ng+1,0:Ng+1), Wgt(Npol+7,0:Ng+1,0:Ng+1,0:Ng+1) )
	allocate ( Pn(0:Npol,0:Ng+1), dPn(0:Npol,0:Ng+1) )
	allocate ( Ae(Npol+7,Npol+7), Le(Npol+7,Npol+7), Be(Npol+7,Npol+7) )
	allocate ( Qce1(Npol+7,Npol+7), Qde1(Npol+7,Npol+7) )
	allocate ( Qce2(Npol+7,Npol+7), Qde2(Npol+7,Npol+7) )
	allocate ( Qce3(Npol+7,Npol+7), Qde3(Npol+7,Npol+7) )
	allocate (                Shpb(Nb,4,Ng,Ng)                  )
	allocate ( BCe(Nb,5), BCi(Nb), BCn(Nb,3), BCvalue(Nb,Ns,4), BCtype(Nb) )
	allocate (order(Ndof,0:1), Rr(Ndof,Ns) )		! ### new parameters ###
	allocate ( I_Omega(Ne,Npol+7), I_Gamma(Ne,Npol+7,6), I_Jump(Nf,8) )	! ### error parameters ###
	allocate ( R_Omega (Ne), J_Gamma (Nf) )				! ###     residual on element ###
	allocate ( epol (Ne) )						! ### global element polynomial array ###
	allocate ( fn  (Ne,6,3) )					! element face unit normals
	allocate ( face(Nf,8) )						! element, local face, node 1 - node 4

! read input
	call READIN(iout, idbg, ipost, Ne, Nn, Nb, Nm, Ng, Ns, Npol, Nf, Ndof, NnNd, ldw, &
			Kappa, tc, aopt, tmax, dt, dto, Nt, theta, theta1, ie, face, &
			x, C, V, D, xg, wg, BCe, BCi, BCvalue, BCtype, Soe, por, &
			vA, vL, vB, vQc1, vQd1, vF1, vQc2, vQd2, vF2, vQc3, vQd3, vF3, ipar, fpar, maxit, Ctol, epol, &
			rA, rL, rB, rQc1, rQd1, rF1, rQc2, rQd2, rF2, rQc3, rQd3, rF3, &
			cA, cL, cB, cQc1, cQd1, cF1, cQc2, cQd2, cF2, cQc3, cQd3, cF3,  &
			lastA, lastL, lastB, lastQc1, lastQd1, lastF1, lastQc2, lastQd2, lastF2, lastQc3, lastQd3, lastF3)

! read parameters from legendre.txt
	open(9, file='legendre.txt', status='old')
	do s = 0,7
	  read (9,*) Dn(s), (Ain(s,i), i=0,7)
	enddo	! s
	close(9)

	 Pn = 0.d0	! initialize (vector form)
	dPn = 0.d0	! initialize (vector form)

	do pol = 0, Npol
	  Ain(pol,:) = Ain(pol,:) / Dn(pol)	! calculate Ain
! compute Legendre polynomial Pn(xi) and its derivative dPn(xi) at Gauss points and at the boundaries
	  call LEGENDRE(idbg, Ng, Npol, xg, Ain, Pn, dPn, pol)
	enddo ! pol

	if(Kappa .ne. 0.) then
! arrays allocation for reactive transpost
	  allocate ( phi(Ndof), vY (NnNd ), rY (Ndof+1), cY (NnNd) )
	endif

! Build a matrix of materials at each node
	call NODALMAT(iout, idbg, Ne, Nn, Nd, ie, nmat)
	ipar0 = ipar	! store the original ipar
	fpar0 = fpar	! store the original fpar

! loop on elements
	do e = 1, Ne

! calculate element shape functions
	  call SHAPE(iout, idbg, Ne, Nn, Ng, Npol, V, ie, x, xg, e, &
			  Shp, dNdr, Wgt, J, Jac, Ji, Jaci, aii, aopt, Pn, dPn)

! calculate element matrices
	  call SHAPEM(iout, idbg, Ne, Ng, Npol, epol, V, D, &
			 Ae, Le, Be, wg, e, Shp, dNdr, Wgt, Jac, Ji)

! calculate nodal flux matrices
	  call SHAPEQ(iout, idbg, Ne, Nn, Ng, Npol, V, D, dNdr, Ji, &
			 Qce1, Qde1, Qce2, Qde2, Qce3, Qde3, epol, e)

! assemble arrays
	  call ASSEMBLE(iout, idbg, Ne, Nn, Nd, Npol, Ndof, NnNd, &
			vA, vL, vB, vQc1, vQd1, vQc2, vQd2, vQc3, vQd3, &
			rA, rL, rB, rQc1, rQd1, rQc2, rQd2, rQc3, rQd3, &
			cA, cL, cB, cQc1, cQd1, cQc2, cQd2, cQc3, cQd3, &
			lastA, lastL, lastB, lastQc1, lastQd1, lastQc2, lastQd2, lastQc3, lastQd3, &
			Ae, Le, Be, Qce1, Qde1, Qce2, Qde2, Qce3, Qde3, ie, nmat, e)

	enddo	! e

! initialize
	call INIT(iout, idbg, Ne, Nn, Ndof, Nb, Nm, Nd, Ns, NnNd, dt, &
			theta, theta1, BCe, BCi, BCn, BCvalue, BCtype, ie, x, nmat, &
			C, T, Cold, Told, Kappa, tc, &
			vA, vL, vB, vQc1, vQd1, vF1, vQc2, vQd2, vF2, vQc3, vQd3, vF3, vY, vZ, &
			rA, rL, rB, rQc1, rQd1, rF1, rQc2, rQd2, rF2, rQc3, rQd3, rF3, rY, rZ, &
			cA, cL, cB, cQc1, cQd1, cF1, cQc2, cQd2, cF2, cQc3, cQd3, cF3, cY, cZ, &
			lastA, lastL, lastB, lastY, lastZ, lastQc1, lastQd1, lastF1, &
			lastQc2, lastQd2, lastF2, lastQc3, lastQd3, lastF3, V, Vm, &
			phi, fn, &
			Nr, order, Rr)		! ### new parameters ###

! calculate the independent source [A]{So} (time-independent)
	call SOURCE0(iout, idbg, Ne, Nn, Ns, Ndof, NnNd, ie, tc, Soe, &
			 vZ, rZ, cZ, lastZ, Son, Sn)
	deallocate (Soe, vZ, rZ, cZ)

! calculate element shape functions on boundary faces
	call SHAPEB(iout, idbg, Ne, Ng, Npol, Nb, BCe, Bci, ie, &
			  Shp, Shpb)

! element error initialization
	call ERROR_INIT(iout, idbg, Ne, Ng, Ns, Npol, Nf, epol, face, &
				V, D, fn, wg, Shp, dNdr, Jac, Ji, &
				aii, I_Omega, I_Gamma, I_Jump)

! write solver output, C, T and nodal fluxes
	call OUT(iout, idbg, ipost, Nn, Ns, Ndof, NnNd, ldw, &
			time, Nc, C, T, w, work, &
			lastQc1, lastQd1, lastF1, lastQc2, lastQd2, lastF2, lastQc3, lastQd3, lastF3, &
			vQc1, vQd1, vF1, rQc1, rQd1, rF1, cQc1, cQd1, cF1, &	
			vQc2, vQd2, vF2, rQc2, rQd2, rF2, cQc2, cQd2, cF2, &	
			vQc3, vQd3, vF3, rQc3, rQd3, rF3, cQc3, cQd3, cF3)	
	timeo = dto

	call CPU_TIME ( time_end )
	write (iout,*) 'Initialization time=', time_end - time_begin, ' seconds'

! time loop
	do n = 1, Nt

	  Sn0 = Sn	! update nodal S at t^n (use matrix form)
! advance in time
	  time = time + dt

! implicit BC and source iteration loop
	  do it = 1, maxit
	    Cit = C	! save former iteration C
! update BC
	    call BC(iout, idbg, Nn, Nb, Ng, Ns, Ndof, BCe, Bci, BCn, BCvalue, BCtype, &
			C, T, Son, Vm, wg, aii, Shpb)

	    if (Kappa .gt. 0.)	then
! update reaction source
	      call SOURCE(iout, idbg, Nn, Ns, Ndof, NnNd, &
			vY, rY, cY, lastY, Kappa, T, phi)
	    endif

! solve equations and update
	    call SOLVE(iout, idbg, Nn, Nb, Ns, Ndof, NnNd, ldw, &
			Kappa, theta, theta1, ipar, fpar, ipar0, fpar0, &
			BCe, BCtype, C, T, Cold, Told, Rn, &
			vA, vL, vB, rA, rL, rB, cA, cL, cB, lastA, lastL, lastB, &
			phi, &
			w, work, &
			Nr, order, Rr)		! ### new parameters ###

! check iteration convegence
	    Cerrmax = 0.	! initialize C error
	    do i = 1, Ndof
	      do s = 1, Ns
		Cerr = abs(C(i,s)-Cit(i,s))
		if(Cerr .gt. Cerrmax)	Cerrmax = Cerr
	      enddo	! s
	    enddo	! i
	    if(Cerrmax .le. Ctol)	exit
	  enddo		! it
	  if(it .gt. maxit)	then
	    write(iout,*) '*** ABORT: BC/source converge failure. it, Cerrmax = ', &
		 it, Cerrmax
	    stop
	  endif

! element error calculation
	call ERROR(iout, idbg, ires, ijmp, Ne, Nn, Nb, Ns, Npol, Nf, Ndof, &
			 time, dt, theta, theta1, ie, epol, face, C, Sn, Cold, Sn0, BCe, &
			 BCvalue, BCtype, I_Omega,  I_Gamma, I_Jump, R_Omega, J_Gamma)

	  if( abs(time - timeo) .le. 0.5*dt ) then
! write solver output, C, T and nodal fluxes
	    call OUT(iout, idbg, ipost, Nn, Ns, Ndof, NnNd, ldw, &
			time, Nc, C, T, w, work, &
			lastQc1, lastQd1, lastF1, lastQc2, lastQd2, lastF2, lastQc3, lastQd3, lastF3, &
			vQc1, vQd1, vF1, rQc1, rQd1, rF1, cQc1, cQd1, cF1, &	
			vQc2, vQd2, vF2, rQc2, rQd2, rF2, cQc2, cQd2, cF2, &	
			vQc3, vQd3, vF3, rQc3, rQd3, rF3, cQc3, cQd3, cF3)	
	    timeo = time + dto
	  endif

! update Cold and Told for the next step
	  Cold = C				! use matrix form
	  Told = T				! use matrix form
	enddo	! n

	call CPU_TIME ( time_end )
	write (iout,*) 'Total time=', time_end - time_begin, ' seconds'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	write(iout,*) '+++ TEMPORARY stop +++'
!	call FLUSH()    ! flush all units
!	STOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end

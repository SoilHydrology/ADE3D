	subroutine READIN(iout, idbg, ipost, Ne, Nn, Nb, Nm, Ng, Ns, Npol, Nf, Ndof, NnNd, ldw, &
			Kappa, tc, aopt, tmax, dt, dto, Nt, theta, theta1, ie, face, &
			x, C, V, D, xg, wg, BCe, BCi, BCvalue, BCtype, Soe, por, &
			vA, vL, vB, vQc1, vQd1, vF1, vQc2, vQd2, vF2, vQc3, vQd3, vF3, ipar, fpar, maxit, Ctol, epol, &
			rA, rL, rB, rQc1, rQd1, rF1, rQc2, rQd2, rF2, rQc3, rQd3, rF3, &
			cA, cL, cB, cQc1, cQd1, cF1, cQc2, cQd2, cF2, cQc3, cQd3, cF3,  &
			lastA, lastL, lastB, lastQc1, lastQd1, lastF1, lastQc2, lastQd2, lastF2, lastQc3, lastQd3, lastF3)
! read input

	implicit none
	integer iout, idbg, ipost
	integer Ne, Nn, Nb, Nm, Ng, Ns, NnNd, ldw	! array parameters
	integer	Npol, Nf, Ndof				! array parameters
	integer lastA, lastL, lastB, lastQc1, lastQd1, lastF1, lastQc2, lastQd2, lastF2, lastQc3, lastQd3, lastF3
	integer Nt				! number of time increments
	integer maxit
	real*8 Kappa, tc, aopt, por
	real*8 tmax, dt, dto
	real*8 theta, theta1, Ctol
	integer ipar(16)			! bCGstab integer parameters array
	real*8 fpar(16)				! bCGstab real    parameters array
	integer BCe(Nb,5), BCi(Nb)		! BC element and local element face numbers
	integer rA (Ndof+1), rL (Ndof+1), rB(Ndof+1)	! global  arrays (compact rows)
	integer cA (NnNd), cL (NnNd), cB(NnNd)	! global  arrays (compact columns)
	integer rQc1(Ndof+1), rQd1(Ndof+1), rF1(Ndof+1)	! global  arrays (compact rows)
	integer rQc2(Ndof+1), rQd2(Ndof+1), rF2(Ndof+1)	! global  arrays (compact rows)
	integer rQc3(Ndof+1), rQd3(Ndof+1), rF3(Ndof+1)	! global  arrays (compact rows)
	integer cQc1(NnNd), cQd1(NnNd), cF1(NnNd)	! global  arrays (compact columns)
	integer cQc2(NnNd), cQd2(NnNd), cF2(NnNd)	! global  arrays (compact columns)
	integer cQc3(NnNd), cQd3(NnNd), cF3(NnNd)	! global  arrays (compact columns)
	real*8 BCvalue(Nb,Ns,4)			! BC value (jx_bar, qx_bar or c_bar)
	character*1 BCtype(Nb)			! BC type ('R', 'N' or 'D')
	integer ie(Ne,9)			! global connectivity array
	integer epol(Ne)			! global element polynomial array
	real*8 x(Nn,3)				! global coordinates array
	real*8 Soe(Ne,Ns,8)			! element independent source nodal values
	real*8 V(Ne,3), D(Ne,3,3)		! global  arrays
	real*8 C   (Ndof,Ns)			! global  arrays
	real*8 vA (NnNd ), vL (NnNd ), vB(NnNd)	! global  arrays (compact values)
	real*8 vQc1(NnNd), vQd1(NnNd), vF1(NnNd)! global  arrays (compact values)
	real*8 vQc2(NnNd), vQd2(NnNd), vF2(NnNd)! global  arrays (compact values)
	real*8 vQc3(NnNd), vQd3(NnNd), vF3(NnNd)! global  arrays (compact values)
	real*8 xg(Ng), wg(Ng)			! Gauss abscissas [-1,+1] and weights
	integer face(Nf,8)			! element, local face, node 1 - node 4

	integer i, e, m, p, s, ierror, i1, i2, i3, i4, f
	real*8 Peclet, h, Vabs, Dabs	! cell Pe, size, |Vi|, |Dik|
	real*8 h17, h28, h35, h46, CFLc, CFLd, CFLcmax, CFLdmax, Pecletmax

	data ierror /0/, CFLcmax/-1.d12/, CFLdmax/-1.d12/, Pecletmax/-1.d12/
!	write(idbg,'(a)') ' --- READIN ---'	! ### TEMPORARY ###

! read input from inp.txt
	read (1,*) Kappa, tc, por		! properties
	read (1,*) tmax, dt, dto
	read (1,*) theta			! temporal scheme parameter
	read (1,*) ipar(2), ipar(3), ipar(6)	! arguments needed to biCGstab
	read (1,*) fpar(1), fpar(2), fpar(11)	! arguments needed to biCGstab
	read (1,*) maxit, Ctol			! implicit BC and source convergence parameters
	close(1)

	ipar(1) = 0				! initialize bCGstab
	ipar(4) = ldw				! work array dimension for bCGstab
	ipar(5) = 0				! unused for bCGstab
	Nt = nint( tmax / dt )
	theta1 = 1.d0 - theta

! reset all global arrays (ARRAY ASSIGNMENT)
	vA   = 0.
	vL   = 0.
	vB   = 0.
	vQc1 = 0.
	vQd1 = 0.
	vF1  = 0.
	vQc2 = 0.
	vQd2 = 0.
	vF2  = 0.
	vQc3 = 0.
	vQd3 = 0.
	vF3  = 0.
	rA   = 1
	rL   = 1
	rB   = 1
	rQc1 = 1
	rQd1 = 1
	rF1  = 1
	rQc2 = 1
	rQd2 = 1
	rF2  = 1
	rQc3 = 1
	rQd3 = 1
	rF3  = 1
	cA   = 0
	cL   = 0
	cB   = 0
	cQc1 = 0
	cQd1 = 0
	cF1  = 0
	cQc2 = 0
	cQd2 = 0
	cF2  = 0
	cQc3 = 0
	cQd3 = 0
	cF3  = 0
! reset last indices for the sparse arrays
	lastA   = 0
	lastL   = 0
	lastB   = 0
	lastQc1 = 0
	lastQd1 = 0
	lastF1  = 0
	lastQc2 = 0
	lastQd2 = 0
	lastF2  = 0
	lastQc3 = 0
	lastQd3 = 0
	lastF3  = 0

! read element connectivity from elements.txt
	open(1, file='elements.txt', status='old')
	do e = 1, Ne
	  read (1,*) i, (ie(e,m), m=1,9)	! element node 1 - node 8, mat
	enddo		! e
	close(1)

! read element faces from faces.txt
	open(1, file='faces.txt', status='old')
	do f = 1, Nf
	  read (1,*) i, (face(f,m), m=1,8)
			! element, local face, node 1 - node 4, neighbour element, neighbour face
	enddo		! f
	close(1)

! read element polynomial order from porder.txt
	open(1, file='porder.txt', status='old')
	do e = 1, Ne
	  read (1,*) i, epol(e)			! element order
	enddo		! e
	close(1)

! read Darcian flux data from V.txt
	open(1, file='V.txt', status='old')
	do e = 1, Ne
	  read (1,*) m, (V(e,i), i=1,3)		! Vx, Vy, Vz
	enddo		! e
	close(1)
	V = V/por				! change from Darcian flux to velocity

! read dispersion data from D.txt
	open(1, file='D.txt', status='old')
	do e = 1, Ne
	  read (1,*) m, ((D(e,i,p), i=1,3), p=1,3)	! Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz
	enddo		! e
	close(1)

! reset C global array (ARRAY ASSIGNMENT)
	C = 0.

! read nodes data from nodes.txt
	open(1, file='nodes.txt', status='old')
! nodal coordinates & IC
	do i = 1, Nn
	  read (1,*) m, x(i,1), x(i,2), x(i,3), (C(i,s), s=1, Ns)
	enddo		! i
! read nodeless IC
!	i = Nn
!	do p = 2, Npol
!	  do e = 1, Ne
!	    i = i+1
!	do i = Nn+1, Ndof
!	  read (1,*)            (C(i,s), s=1, Ns)
!	enddo		! i
!	  enddo		! e
!	enddo		! p
	close(1)

! read element face BC data from BCs.txt
	open(1, file='BCs.txt', status='old')
! BC element, BC face, BC type and 4 BC values at the BC element face nodes
	do i = 1, Nb
	  read (1,*) BCi(i), BCe(i,5), BCtype(i), ((BCvalue(i,s,m), m=1,4), s=1,Ns)
	enddo		! i
	close(1)

! read element independent source data from source.txt
	open(1, file='source.txt', status='old')
! 8 independent source values at the element nodes
	do e = 1, Ne
	  do s = 1, Ns
	    read (1,*) m, p, (Soe(e,s,i), i=1,8)
	  enddo		! s
	enddo		! e
	close(1)


! Gauss abscissas [-1,+1] and weights
	call LEGENDRE_SET ( Ng, xg, wg )
! --------------------------------- ECHO INPUT ---------------------------------
	write(iout,*) 'Kappa, tc, por = ', Kappa, tc, por
	write(iout,*) 'tmax, dt, dto, Nt'
	write(iout,*)  tmax, dt, dto, Nt
	write(iout,*) 'theta= ', theta
	write(iout,*) 'ipar = ', ipar
	write(iout,*) 'fpar = ', fpar
	write(iout,*) 'maxit, Ctol =', maxit, Ctol

	write(iout,*)  'i, xg(i), wg(i)'
	do i = 1, Ng
	  write(iout,*) i, xg(i), wg(i)
	enddo		! i

! the stability criteria for explicit scheme are:
! o v*dt/h   <1    (advection)
! o D*dt/h^2 < 0.5 (dispersion)
! assume the cell size, h, is the larger main diagonal
! assume V = |Vi|, D = |Dik|
	do e = 1, Ne

	  h17 = ( x(ie(e,7),1) - x(ie(e,1),1) )**2 + &
		( x(ie(e,7),2) - x(ie(e,1),2) )**2 + &
		( x(ie(e,7),3) - x(ie(e,1),3) )**2			! diagonal 1->7 ^2
	  h28 = ( x(ie(e,8),1) - x(ie(e,2),1) )**2 + &
		( x(ie(e,8),2) - x(ie(e,2),2) )**2 + &
		( x(ie(e,8),3) - x(ie(e,2),3) )**2			! diagonal 2->8 ^2
	  h35 = ( x(ie(e,5),1) - x(ie(e,3),1) )**2 + &
		( x(ie(e,5),2) - x(ie(e,3),2) )**2 + &
		( x(ie(e,5),3) - x(ie(e,3),3) )**2			! diagonal 3->5 ^2
	  h46 = ( x(ie(e,6),1) - x(ie(e,4),1) )**2 + &
		( x(ie(e,6),2) - x(ie(e,4),2) )**2 + &
		( x(ie(e,6),3) - x(ie(e,4),3) )**2			! diagonal 4->6 ^2

	  h = sqrt( max( h17, h28, h35, h46 ) )				! element size
	  Vabs = sqrt(  V(e,1)**2 + V(e,2)**2 + V(e,3)**2 )		! |Vi|
	  Dabs = sqrt(  D(e,1,1)**2 + D(e,1,2)**2 + D(e,1,3)**2 + &
			D(e,2,1)**2 + D(e,2,2)**2 + D(e,2,3)**2 + &
			D(e,3,1)**2 + D(e,3,2)**2 + D(e,3,3)**2 )	! sqrt(Dik*Dik)
	  Peclet = 0.5d0 * Vabs * h / Dabs	! Pe
	  CFLc = Vabs*dt/h			! advection  CFL
	  CFLd = Dabs*dt/h**2			! dispersion CFL
	  if(CFLc   .gt. CFLcmax)	CFLcmax = CFLc
	  if(CFLd   .gt. CFLdmax)	CFLdmax = CFLd
	  if(Peclet .gt. Pecletmax)	Pecletmax = Peclet
	enddo		! e

! calculate aopt = ctgh(|Pe|) - 1/|Pe| where Pe = (|Vi|*h) / (2*D) (half of my Pe_h !!!)
	Peclet = Pecletmax
	aopt = ( exp(Peclet) + exp(-Peclet) ) / ( exp(Peclet) - exp(-Peclet) ) - 1./Peclet
	if(Pecletmax .le. 1.  ) aopt = 0.				! SUPG not needed
	if(Pecletmax .ge. 700.) aopt = 1.				! saturation reached

	write(iout,*)  'maximum advection   CFL (Vabs * dt / h)       = ', CFLcmax
	write(iout,*)  'maximum dispersion  CFL (Dabs * dt / h**2)    = ', CFLdmax
	write(iout,*)  'maximum Peclet number (0.5 * Vabs * h / Dabs) = ', Pecletmax
	write(iout,*)  'SUPG parameter, aopt = ', aopt
! --------------------------------- ECHO INPUT ---------------------------------

! write post-processing Gmsh file
! header
	write(ipost,'(a)') '$MeshFormat'
	write(ipost,'(a)') '2.2 0 8	// version-number file-type data-size'
	write(ipost,'(a)') '$EndMeshFormat'
	write(ipost,*)
! nodes
	write(ipost,'(a)') '$Nodes'
	write(ipost,'(i10, a)') Nn, '	// number-of-nodes'
	do i = 1, Nn
	  write(ipost,'(i10,3g14.7)') i, x(i,1), x(i,2), x(i,3)
	enddo		! i
	write(ipost,'(a)')  '// node-number x-coord y-coord z-coord'
	write(ipost,'(a)') '$EndNodes'
	write(ipost,*)
! elements
	write(ipost,'(a)') '$Elements'
	write(ipost,'(i10, a)') Ne, '	// number-of-elements'
	do e = 1, Ne
	  write(ipost,'(99i10)') e, 5, 1, ie(e,9), (ie(e,i),i=1,8)
	enddo		! i
	write(ipost,'(a)') '// elm-number elm-type number-of-tags < tag > ... node-number-list'
	write(ipost,'(a)') '$EndElements'
	write(ipost,*)
! velocity vector
	write(ipost,'(a)') '$ElementData'
	write(ipost,'(i10, a)') 1, '	// number-of-string-tags'
	write(ipost,'(a)') '"velocity vector"	// < "string-tag" > ...'
	write(ipost,'(i10, a)') 1, '	// number-of-real-tags'
	write(ipost,'(g14.7, a)') 0.d0, '	// < real-tag > ... time'
	write(ipost,'(i10, a)') 3, '	// number-of-integer-tags'
	write(ipost,'(i10, a)') 0, '	// step #'
	write(ipost,'(i10, a)') 3, '	// number of components'
	write(ipost,'(i10, a)') Ne, '	// number of elements'
	do e = 1, Ne
	  write(ipost,'(i10, 3g14.7)') e, (V(e,i),i=1,3)
	enddo		! e
	write(ipost,'(a)') '// elm-number number-of-nodes-per-element value ... Vx, Vy, Vz'
	write(ipost,'(a)') '$EndElementData'
	write(ipost,*)
! diperssion tensor
	write(ipost,'(a)') '$ElementData'
	write(ipost,'(i10, a)') 1, '	// number-of-string-tags'
	write(ipost,'(a)') '"diperssion tensor"	// < "string-tag" > ...'
	write(ipost,'(i10, a)') 1, '	// number-of-real-tags'
	write(ipost,'(g14.7, a)') 0.d0, '	// < real-tag > ... time'
	write(ipost,'(i10, a)') 3, '	// number-of-integer-tags'
	write(ipost,'(i10, a)') 0, '	// step #'
	write(ipost,'(i10, a)') 9, '	// number of components'
	write(ipost,'(i10, a)') Ne, '	// number of elements'
	do e = 1, Ne
	  write(ipost,'(i10, 9g14.7)') e,   (D(e,1,i), i=1,3), &
				        (D(e,2,i), i=1,3), &
					(D(e,3,i), i=1,3)
	enddo		! e
	write(ipost,'(a)') '// elm-number number-of-nodes-per-element value ... Dxx, Dxy,.. Dzz'
	write(ipost,'(a)') '$EndElementData'
	write(ipost,*)
! element polynomial order
	write(ipost,'(a)') '$ElementData'
	write(ipost,'(i10, a)') 1, '	// number-of-string-tags'
	write(ipost,'(a)') '"element polynomial order"	// < "string-tag" > ...'
	write(ipost,'(i10, a)') 1, '	// number-of-real-tags'
	write(ipost,'(g14.7, a)') 0.d0, '	// < real-tag > ... time'
	write(ipost,'(i10, a)') 3, '	// number-of-integer-tags'
	write(ipost,'(i10, a)') 0, '	// step #'
	write(ipost,'(i10, a)') 1, '	// number of components'
	write(ipost,'(i10, a)') Ne, '	// number of elements'
	do e = 1, Ne
	  write(ipost,'(2i10)') e, (epol(e))
	enddo		! e
	write(ipost,'(a)') '// elm-number number-of-nodes-per-element value ... epol'
	write(ipost,'(a)') '$EndElementData'
	write(ipost,*)

! check general input
	if(por .le. 0. .or. por .gt. 1.) then
	  write(iout,*) '*** ABORT: porosity is not between 0 and 1'
	  ierror = ierror + 1
	endif
	if(Kappa .lt. 0.) then
	  write(iout,*) '*** ABORT: Kappa < 0'
	  ierror = ierror + 1
	endif
	if(Kappa .gt. 0. .and. Ns .ne. 3) then
	  write(iout,*) '*** ABORT: for Kappa>0 Ns must be 3'
	  ierror = ierror + 1
	endif
	if(dt .le. 0. .or. dto .le. 0.) then
	  write(iout,*) '*** ABORT: dt .le. 0 or dto .le. 0'
	  ierror = ierror + 1
	endif
	if(Nt .lt. 0) then
	  write(iout,*) '*** ABORT: Nt .lt. 0'
	  ierror = ierror + 1
	endif
	if(Ng .lt. 2 .or. Ng .gt. 8 ) then
	  write(iout,*) '*** ABORT: Ng between 2 and 8 only are currently allowed for'
	  ierror = ierror + 1
	endif
	if(theta .lt. 0. .or. theta .gt. 1.) then
	  write(iout,*) '*** ABORT: theta <0 or theta>1'
	  ierror = ierror + 1
	endif
	if(fpar(1) .le. 0. .or. fpar(2) .le. 0.) then
	  write(iout,*) '*** ABORT: fpar(1) .le. 0. .or. fpar(2) .le. 0.'
	  ierror = ierror + 1
	endif
	if(ipar(2) .lt. 0 .or. ipar(2) .gt. 3) then
	  write(iout,*) '*** ABORT: ipar(2) .lt. 0 .or. ipar(2) .gt. 3'
	  ierror = ierror + 1
	endif
	if(ipar(3) .lt. -2 .or. ipar(2) .gt. 2) then
	  write(iout,*) '*** ABORT: ipar(3) .lt. -2 .or. ipar(2) .gt. 2'
	  ierror = ierror + 1
	endif
! check element input
	do e = 1, Ne
	  if(ie(e,1) .lt. 1 .or. ie(e,1) .gt. Nn .or. &
	     ie(e,2) .lt. 1 .or. ie(e,2) .gt. Nn .or. &
	     ie(e,3) .lt. 1 .or. ie(e,3) .gt. Nn .or. &
	     ie(e,4) .lt. 1 .or. ie(e,4) .gt. Nn .or. &
	     ie(e,5) .lt. 1 .or. ie(e,5) .gt. Nn .or. &
	     ie(e,6) .lt. 1 .or. ie(e,6) .gt. Nn .or. &
	     ie(e,7) .lt. 1 .or. ie(e,7) .gt. Nn .or. &
	     ie(e,8) .lt. 1 .or. ie(e,8) .gt. Nn .or. &
	     ie(e,9) .lt. 1 .or. ie(e,9) .gt. Nm ) then
	    write(iout,*) '*** ABORT: illegal ie(e,*), e, ie(e,*) = ', &
						       e, (ie(e,m), m=1,9)
	    ierror = ierror + 1
	  endif
	enddo		! e
	do e = 1, Ne
	  if(epol(e) .lt. 1 .or. epol(e) .gt. Npol) then
	    write(iout,*) '*** ABORT: illegal epol(e), e, epol(e) = ', &
						       e, epol(e)
	    ierror = ierror + 1
	  endif
	enddo		! e
	do e = 1, Ne
	  if(   D(e,1,1).lt.0. .or. D(e,1,2).lt.0. .or. D(e,1,3).lt.0. .or. &
		D(e,2,1).lt.0. .or. D(e,2,2).lt.0. .or. D(e,2,3).lt.0. .or. &
		D(e,3,1).lt.0. .or. D(e,3,2).lt.0. .or. D(e,3,3).lt.0. ) then
	    write(iout,*) '*** ABORT: D(e,i,j)<0 for e = ', e
	    ierror = ierror + 1
	  endif
	enddo		! e
! check element face input
	do f = 1, Nf
	  if(face(f,1) .lt. 0 .or. face(f,1) .gt. Ne .or. &
	     face(f,2) .lt. 0 .or. face(f,2) .gt. 6  .or. &
	     face(f,3) .lt. 0 .or. face(f,3) .gt. Nn .or. &
	     face(f,4) .lt. 0 .or. face(f,4) .gt. Nn .or. &
	     face(f,5) .lt. 0 .or. face(f,5) .gt. Nn .or. &
	     face(f,6) .lt. 0 .or. face(f,6) .gt. Nn)	then
	    write(iout,*) '*** ABORT: illegal face(f,*), f, face(f,*) = ', &
						         f,(face(f,m), m=1,6)
	    ierror = ierror + 1
	  endif
	enddo		! f
! check BCs input
	do i = 1, Nb
	  if(BCe(i,5) .lt. 1 .or. BCe(i,5) .gt. Ne) then
	    write(iout,*) '*** ABORT: illegal BCe'
	    ierror = ierror + 1
	  endif
	  if(BCi(i) .lt. 1 .or. BCi(i) .gt. 6) then
	    write(iout,*) '*** ABORT: illegal BCi'
	    ierror = ierror + 1
	  endif
	  if( BCtype(i).ne.'R' .and. BCtype(i).ne.'N' .and. BCtype(i).ne.'D') then
	    write(iout,*) '*** ABORT: illegal i, BCtype(i) = ', i, BCtype(i)
	    ierror = ierror + 1
	  endif
	enddo		! i

	if(ierror .ne. 0) then
	  write(iout,*) '*** ABORT: ierror = ', ierror
	  stop
	endif
	
	return
	end

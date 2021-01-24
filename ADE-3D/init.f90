	subroutine INIT(iout, idbg, Ne, Nn, Ndof, Nb, Nm, Nd, Ns, NnNd, dt, &
			theta, theta1, BCe, BCi, BCn, BCvalue, BCtype, ie, x, nmat, &
			C, T, Cold, Told, Kappa, tc, &
			vA, vL, vB, vQc1, vQd1, vF1, vQc2, vQd2, vF2, vQc3, vQd3, vF3, vY, vZ, &
			rA, rL, rB, rQc1, rQd1, rF1, rQc2, rQd2, rF2, rQc3, rQd3, rF3, rY, rZ, &
			cA, cL, cB, cQc1, cQd1, cF1, cQc2, cQd2, cF2, cQc3, cQd3, cF3, cY, cZ, &
			lastA, lastL, lastB, lastY, lastZ, lastQc1, lastQd1, lastF1, &
			lastQc2, lastQd2, lastF2, lastQc3, lastQd3, lastF3, V, Vm, &
			phi, fn, &
			Nr, order, Rr)		! ### new parameters ###
! initialize

	implicit none
	integer iout, idbg
	integer Ne, Nn, Nb, Nm, Nd, Ns, NnNd	! array parameters
	integer	Ndof				! array parameters
	integer lastA, lastL, lastB, lastQc1, lastQd1, lastF1, lastQc2, lastQd2, lastF2, lastQc3, lastQd3, lastF3 
	integer lastY, lastZ
	real*8 Kappa, tc
	real*8 dt
	real*8 theta, theta1
	integer BCe(Nb,5), BCi(Nb)		! BC element and local element face numbers
	integer rA (Ndof+1), rL (Ndof+1), rB(Ndof+1)	! global  arrays (compact rows)
	integer cA (NnNd), cL (NnNd), cB(NnNd)	! global  arrays (compact columns)
	integer rY (Ndof+1), cY (NnNd)		! global  arrays
	integer rZ (Ndof+1), cZ (NnNd)		! global  arrays
	integer rQc1(Ndof+1), rQd1(Ndof+1), rF1(Ndof+1)	! global  arrays (compact rows)
	integer rQc2(Ndof+1), rQd2(Ndof+1), rF2(Ndof+1)	! global  arrays (compact rows)
	integer rQc3(Ndof+1), rQd3(Ndof+1), rF3(Ndof+1)	! global  arrays (compact rows)
	integer cQc1(NnNd), cQd1(NnNd), cF1(NnNd)	! global  arrays (compact columns)
	integer cQc2(NnNd), cQd2(NnNd), cF2(NnNd)	! global  arrays (compact columns)
	integer cQc3(NnNd), cQd3(NnNd), cF3(NnNd)	! global  arrays (compact columns)
	real*8 BCvalue(Nb,Ns,4)			! BC value (jx_bar, qx_bar or c_bar)
	real*8 BCn(Nb,3)			! BC outwards normal
	character*1 BCtype(Nb)			! BC type ('R', 'N' or 'D')
	integer ie(Ne,9)			! global connectivity array
	real*8 x(Nn,3)				! global coordinates array
	integer nmat(Nn,0:Nd)			! global nodal materials array
	real*8 V(Ne,3)				! global  arrays
	real*8 fn(Ne,6,3)			! element face unit normals
	real*8 C   (Ndof,Ns), T   (Ndof,Ns)	! global  arrays
	real*8 Cold(Ndof,Ns), Told(Ndof,Ns) 	! global  arrays
	real*8 phi(Ndof)			! global  arrays
	real*8 vA (NnNd ), vL (NnNd ), vB(NnNd)	! global  arrays (compact values)
	real*8 vY (NnNd ), vZ (NnNd )		! global  arrays (compact values)
	real*8 vQc1(NnNd), vQd1(NnNd), vF1(NnNd)! global  arrays (compact values)
	real*8 vQc2(NnNd), vQd2(NnNd), vF2(NnNd)! global  arrays (compact values)
	real*8 vQc3(NnNd), vQd3(NnNd), vF3(NnNd)! global  arrays (compact values)
	real*8 Vm(Nn,3)				! nodal averaged array

	integer N_D, Nr				! ### new parameters ###
	integer order(Ndof,0:1)			! ### new parameters ###
	real*8 Rr(Ndof,Ns)			! ### new parameters ###
	real*8, allocatable ::  Cr(:,:)		! ### new parameters ###

	integer, allocatable :: ec(:)			! element counter at each node

	integer, allocatable :: rBo(:)			! global  arrays (compact rows)
	integer, allocatable :: cBo(:)			! global  arrays (compact columns)
	real*8, allocatable :: vBo (:)			! global  arrays (compact values)
	integer i, j, k, m, n, p, e, i1, i2, i3, i4, s, f
	real*8 alpha, beta, sx1, sy1, sz1, sx2, sy2, sz2, nx, ny, nz, Lx

!	write(idbg,'(a)') ' --- INIT ---'	! ### TEMPORARY ###
    
! reset T to 0
	T = 0		! use matrix form

! store A in Z
	vZ = vA	! use matrix form
	rZ = rA	! use matrix form
	cZ = cA	! use matrix form
	lastZ = lastA

! average V on all elements interfacing node i
	allocate ( ec(Nn) )	! allocate ec
	ec = 0			! use matrix form
	Vm = 0.d0		! use matrix form
	do e = 1, Ne
	  do i = 1, 8
	    n = ie(e,i)		! node i of element e
	    ec(n) = ec(n) + 1	! sum number of elements at this node
	    do i1 = 1, 3
	      Vm(n,i1) = Vm(n,i1) + V(e,i1)
	    enddo	! i1
	  enddo		! i
	enddo		! e
! average nodal velocity by number of elements at each node
	do n = 1, Nn
	  do i1 = 1, 3
	    Vm(n,i1) = Vm(n,i1) / ec(n)
	  enddo		! i1
	enddo		! n
	deallocate ( ec )		! deallocate ec

! calculate B=P-L and store in B
	do i = 1 ,Ndof
	  do j = 1 ,Ndof
	    sx1 = 0.
! access rank 2 sparse arrays
	    call ACCESS2(iout, idbg, Ndof, NnNd, rL, cL, i, j, lastL, m)
	    call ACCESS2(iout, idbg, Ndof, NnNd, rB, cB, i, j, lastB, n)

	    if (m .ne. 0)	then
	      sx1 = vL(m)		! store L(i,j) in sx1
	    endif

	    if (n .ne. 0)	then
	      vB(n) = vB(n) - sx1	! store B(i,j)-L(i,j) in vB
	    endif

	  enddo	! j
	enddo	! i

! store A in Y
	  if(Kappa .ne. 0.) then
	    vY = vA	! use matrix form
	    rY = rA	! use matrix form
	    cY = cA	! use matrix form
	    lastY = lastA
	  endif

! calculate alpha = A/dt + B*theta  and store in A
! calculate beta  = A/dt - B*theta1 and store in L
	  do i = 1, Ndof

	    if(Kappa .ne. 0.) then
	      phi(i) = C(i,1)*C(i,2)	! initialize phi
	    endif
	    do j = 1 ,Ndof
	      sx1 = 0.
	      sy1 = 0.
! access rank 2 sparse arrays
	      call ACCESS2(iout, idbg, Ndof, NnNd, rA, cA, i, j, lastA, k)
	      call ACCESS2(iout, idbg, Ndof, NnNd, rL, cL, i, j, lastL, m)
	      call ACCESS2(iout, idbg, Ndof, NnNd, rB, cB, i, j, lastB, n)

	      if (k .ne. 0)	then
	        sx1 = vA(k)	! store A(i,j) in sx1
	      endif

	      if (n .ne. 0)	then
	        sy1 = vB(n)	! store B(i,j) in sy1
	      endif

	      if (k .ne. 0 .or. n .ne. 0)	then	! skip if both A and B are 0
	        alpha = sx1/dt + sy1*theta
	        beta =  sx1/dt - sy1*theta1
	      endif

	      if (k .ne. 0)	then
	        vA(k) = alpha
	      endif

	      if (m .ne. 0)	then
	        vL(m) = beta
	      endif
	    enddo	! j
	  enddo		! i

	allocate ( vBo(NnNd), cBo(NnNd), rBo(Ndof+1) )

! create a list of the original vs reordered nodes
	call DORDER(iout, idbg, Ne, Nn, Nb, Ndof, BCe, BCi, BCtype, ie, &
			N_D, order)	! ### new parameters ###

! permute A rows & columns into Bo if Dirichlet BC
	call DPERM (Ndof, vA, cA, rA, vBo, cBo, rBo, order(:,1), order(:,1), 1)

	Nr = Ndof - N_D		! reduced rows and columns
	allocate( Cr(Ndof,Ns) )	! constant part
	Cr = 0.			! initialize Cr

! calculate element faces unit normals
	do e = 1,Ne
	  do f = 1,6
	    if     (f .eq. 1) then
! face 1 - west
	      i1 = ie(e,1)
	      i2 = ie(e,4)
	      i3 = ie(e,8)
	      i4 = ie(e,5)
	    else if(f .eq. 2) then
! face 2 - south
	      i1 = ie(e,1)
	      i2 = ie(e,5)
	      i3 = ie(e,6)
	      i4 = ie(e,2)
	    else if(f .eq. 3) then
! face 3 - bottom
	      i1 = ie(e,1)
	      i2 = ie(e,2)
	      i3 = ie(e,3)
	      i4 = ie(e,4)
	    else if(f .eq. 4) then
! face 4 - east
	      i1 = ie(e,2)
	      i2 = ie(e,3)
	      i3 = ie(e,7)
	      i4 = ie(e,6)
	    else if(f .eq. 5) then
! face 5 - north
	      i1 = ie(e,4)
	      i2 = ie(e,8)
	      i3 = ie(e,7)
	      i4 = ie(e,3)
	    else if(f .eq. 6) then
! face 6 - top
	      i1 = ie(e,5)
	      i2 = ie(e,6)
	      i3 = ie(e,7)
	      i4 = ie(e,8)
	    endif
! i1->i3 vector
	    sx1 = x(i3,1)-x(i1,1)	! i1->i3 x component
	    sy1 = x(i3,2)-x(i1,2)	! i1->i3 y component
	    sz1 = x(i3,3)-x(i1,3)	! i1->i3 z component
! i2->i4 vector
	    sx2 = x(i4,1)-x(i2,1)	! i2->i4 x component
	    sy2 = x(i4,2)-x(i2,2)	! i2->i4 y component
	    sz2 = x(i4,3)-x(i2,3)	! i2->i4 z component
! normal vector
	    nx = sy1*sz2 - sz1*sy2
	    ny = sz1*sx2 - sx1*sz2
	    nz = sx1*sy2 - sy1*sx2

	    Lx = sqrt(nx**2 + ny**2 + nz**2)	! normal length
! switch normal direction for west, south and bottom faces
	    if(f.eq.1 .or. f.eq.2 .or. f.eq.3)	Lx= -Lx
! element e, face f unit normal componet i
	    fn(e,f,1) = nx/Lx		! unit normal x component
	    fn(e,f,2) = ny/Lx		! unit normal y component
	    fn(e,f,3) = nz/Lx		! unit normal z component
	  enddo	! f
	enddo	! e

! calculate boundary faces unit normals
	do n = 1, Nb

	  e = BCe(n,5)	! BC element number
	  i = BCi(n)	! BC local element face number
	  if     (i .eq. 1) then
! face 1 - west
	    i1 = ie(e,1)
	    i2 = ie(e,4)
	    i3 = ie(e,8)
	    i4 = ie(e,5)
	  else if(i .eq. 2) then
! face 2 - south
	    i1 = ie(e,1)
	    i2 = ie(e,5)
	    i3 = ie(e,6)
	    i4 = ie(e,2)
	  else if(i .eq. 3) then
! face 3 - bottom
	    i1 = ie(e,1)
	    i2 = ie(e,2)
	    i3 = ie(e,3)
	    i4 = ie(e,4)
	  else if(i .eq. 4) then
! face 4 - east
	    i1 = ie(e,2)
	    i2 = ie(e,3)
	    i3 = ie(e,7)
	    i4 = ie(e,6)
	  else if(i .eq. 5) then
! face 5 - north
	    i1 = ie(e,4)
	    i2 = ie(e,8)
	    i3 = ie(e,7)
	    i4 = ie(e,3)
	  else if(i .eq. 6) then
! face 6 - top
	    i1 = ie(e,5)
	    i2 = ie(e,6)
	    i3 = ie(e,7)
	    i4 = ie(e,8)
	  endif
! i1->i3 vector
	  sx1 = x(i3,1)-x(i1,1)		! i1->i3 x component
	  sy1 = x(i3,2)-x(i1,2)		! i1->i3 y component
	  sz1 = x(i3,3)-x(i1,3)		! i1->i3 z component
! i2->i4 vector
	  sx2 = x(i4,1)-x(i2,1)		! i2->i4 x component
	  sy2 = x(i4,2)-x(i2,2)		! i2->i4 y component
	  sz2 = x(i4,3)-x(i2,3)		! i2->i4 z component
! normal vector
	  nx = sy1*sz2 - sz1*sy2
	  ny = sz1*sx2 - sx1*sz2
	  nz = sx1*sy2 - sy1*sx2

	  Lx = sqrt(nx**2 + ny**2 + nz**2)	! normal length
! switch normal direction for west, south and bottom faces
	  if(i.eq.1 .or. i.eq.2 .or. i.eq.3)	Lx= -Lx

	  BCn(n,1) = nx/Lx		! unit normal x component
	  BCn(n,2) = ny/Lx		! unit normal y component
	  BCn(n,3) = nz/Lx		! unit normal z component

	  BCe(n,1) = i1			! store BC node 1 in BCe
	  BCe(n,2) = i2			! store BC node 2 in BCe
	  BCe(n,3) = i3			! store BC node 3 in BCe
	  BCe(n,4) = i4			! store BC node 4 in BCe

! replace C(i,s) by Cbar if Dirichlet BC
	  if(BCtype(n) .eq. 'D') then
	    do p = 1, 4
	      i = BCe(n,p)		! BC global node p number
	      do s = 1, Ns
	        Cr(i,s) = BCvalue(n,s,p)
	        C (i,s) = BCvalue(n,s,p)
	      enddo	! s
	    enddo	! p
	  endif
	enddo		! n

	do s = 1, Ns
! permute Cr vector (in-place)
	  call DVPERM (Ndof, Cr(:,s), order(:,1)) 

! calculate {Ra} = [Baa]{0} + [Bab]{Cb}
	  call AMUX(Ndof, Cr(:,s), Rr(:,s), vBo, cBo, rBo)
! make {Rb} = 0
	  do i = Nr+1, Ndof
	    Rr(i,s) = 0.
	  enddo	! i
	enddo	! s

! extract submatrix Baa into A
	call SUBMAT (Ndof ,1, 1, Nr, 1, Nr, vBo, cBo, rBo, Nr, Nr, vA, cA, rA)
	lastA = rA(Nr+1) - rA(1)	! update lastA
	deallocate ( Cr )		! deallocate Cr
	deallocate ( vBo, cBo, rBo )	! deallocate Bo

! reset C and T of the former time step, Cold and Told, to C and T, respectively
	Cold = C	! use matrix form
	Told = T	! use matrix form
	
!	F = Qc + Qd			! flux

! copy Qd to F
	rF1 = rQd1
	rF2 = rQd2
	rF3 = rQd3
	cF1 = cQd1
	cF2 = cQd2
	cF3 = cQd3
	vF1 = vQd1
	vF2 = vQd2
	vF3 = vQd3
	lastF1 = lastQd1
	lastF2 = lastQd2
	lastF3 = lastQd3

	do j = 1, Ndof
	  do k = 1, Ndof
	    sx1 = 0.
! access rank 2 sparse arrays
	    call ACCESS2(iout, idbg, Ndof, NnNd, rQc1, cQc1, j, k, lastQc1, m)
	    call ACCESS2(iout, idbg, Ndof, NnNd, rF1 , cF1 , j, k, lastF1 , n)

	    if (m .ne. 0)	then
	      sx1 = vQc1(m)		! store Qc1(j,k) in sx1
	    endif

	    if (n .ne. 0)	then
	      vF1(n) = vF1(n) + sx1	! store Qd1(j,k)+Qc1(j,k) in F1
	    endif
	    sy1 = 0.
! access rank 2 sparse arrays
	    call ACCESS2(iout, idbg, Ndof, NnNd, rQc2, cQc2, j, k, lastQc2, m)
	    call ACCESS2(iout, idbg, Ndof, NnNd, rF2 , cF2 , j, k, lastF2 , n)

	    if (m .ne. 0)	then
	      sy1 = vQc2(m)		! store Qc2(j,k) in sy1
	    endif

	    if (n .ne. 0)	then
	      vF2(n) = vF2(n) + sy1	! store Qd2(j,k)+Qc2(j,k) in F2
	    endif
	    sz1 = 0.
! access rank 2 sparse arrays
	    call ACCESS2(iout, idbg, Ndof, NnNd, rQc3, cQc3, j, k, lastQc3, m)
	    call ACCESS2(iout, idbg, Ndof, NnNd, rF3 , cF3 , j, k, lastF3 , n)

	    if (m .ne. 0)	then
	      sz1 = vQc3(m)		! store Qc3(j,k) in sz1
	    endif

	    if (n .ne. 0)	then
	      vF3(n) = vF3(n) + sz1	! store Qd3(j,k)+Qc3(j,k) in F3
	    endif

	  enddo	! k
	enddo	! j

	return
	end


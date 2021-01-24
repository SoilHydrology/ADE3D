	subroutine BC(iout, idbg, Nn, Nb, Ng, Ns, Ndof, BCe, Bci, BCn, BCvalue, BCtype, &
			C, T, Son, Vm, wg, aii, Shpb)
! update BC
! *** IMPORTANT: T is positive for outflow and negative for inflow ***
! NOTE: on the boundaries W = N by definition

	implicit none
	integer iout, idbg
	integer Nn, Nb, Ng, Ns			! array parameters
	integer	Ndof				! array parameters
	integer BCe(Nb,5), BCi(Nb)		! BC element and local element face numbers
	real*8 BCn(Nb,3)			! BC outwards normal
	real*8 BCvalue(Nb,Ns,4)			! BC value (jx_bar, qx_bar or c_bar)
	character*1 BCtype(Nb)			! BC type ('R', 'N' or 'D')
	real*8 C   (Ndof,Ns), T   (Ndof,Ns)	! global  arrays
	real*8 Son(Ndof,Ns)			! global  arrays
	real*8 wg(Ng)				! Gauss weights
	real*8 aii(3,0:Ng+1,0:Ng+1,0:Ng+1)	! the minors a11, a22, a33 of Jij
	real*8 Shpb(Nb,4,Ng,Ng)			! boundary shape functions
	real*8 Vm(Nn,3)				! nodal averaged array

	integer i, j, k, n, s, g1b, g2b, g3b, Lf, f	! local indices
	integer ii				! global index
	real*8 BCv(4), Nx, Ny, Nz, w, Te(4), Ij(4), Vn(4), Fe(4,4)

!	write(idbg,'(a)') ' --- BC ---'	! ### TEMPORARY ###

! reset T to Son in case that a node has more than one BC
	T = Son		! use matrix form

	do s = 1, Ns
	  do n = 1, Nb
	    Nx = BCn(n,1)		! BC Nx
	    Ny = BCn(n,2)		! BC Ny
	    Nz = BCn(n,3)		! BC Ny
	    BCv(:) = BCvalue(n,s,:)	! BC value

	    Lf = BCi(n)			! BC local face number
	    if     (Lf.eq.1) then
! west face
	      f = 1
	      g3b = 0
	    else if(Lf.eq.4) then
! east face
	      f = 1
	      g3b = Ng+1
	    else if(Lf.eq.2) then
! south face
	      f = 2
	      g3b = 0
	    else if(Lf.eq.5) then
! north face
	      f = 2
	      g3b = Ng+1
	    else if(Lf.eq.3) then
! bottom face
	      f = 3
	      g3b = 0
	    else if(Lf.eq.6) then
! top face
	      f = 3
	      g3b = Ng+1
	    endif

! reset Te, Fe to 0
	    Te = 0.
	    Fe = 0.

	    if      (BCtype(n) .eq. 'R') then
! Robin BC
!---------
	      do i = 1,4
	        ii = BCe(n,i)					! BC global node i number
	        Vn(i) = Vm(ii,1)*Nx + Vm(ii,2)*Ny + Vm(ii,3)*Nz	! normal velocity at local node i
	          Ij(i) = C(ii,s)			! I <- C
	      enddo	! i
	    
	      do g1b = 1, Ng
	        do g2b = 1, Ng
	          w = wg(g1b) * wg(g2b)
	          do i = 1,4
	            do j = 1,4
	              Te(i) = Te(i) + w * Shpb(n,i,g1b,g2b) * BCv(j) * &
					aii(f,g1b,g2b,g3b)
! calculate integral[(W_I*n_i*v_i*N_J)]dG
	              Fe(i,j) = Fe(i,j) + w * Shpb(n,i,g1b,g2b) * Vn(i) * Shpb(n,j,g1b,g2b) * &
					aii(f,g1b,g2b,g3b)
	            enddo	! j
	          enddo		! i
	        enddo		! g2b
	      enddo		! g1b

! calculate and store in T
	      do i = 1,4
	        ii= BCe(n,i)				! BC global node i number
	        do j = 1,4
		  T(ii,s) = T(ii,s) - Fe(i,j)*Ij(j)	! add first term
	        enddo	! j
	        T(ii,s) = T(ii,s) + Te(i)		! add second term
	      enddo	! i

!	    else if (BCtype(n) .eq. 'D') then
! Dirichlet BC is implemented in INIT and SOLVE
!-------------

	    else if (BCtype(n) .eq. 'N') then
! Neumann BC
!-----------
	      do g1b = 1, Ng
	        do g2b = 1, Ng
	          w = wg(g1b) * wg(g2b)
	          do i = 1,4
	            do j = 1,4
	              Te(i) = Te(i) + w * Shpb(n,i,g1b,g2b) * BCv(j) * &
					aii(f,g1b,g2b,g3b)
	            enddo	! j
	          enddo		! i
	        enddo		! g2b
	      enddo		! g1b

! calculate and store in T
	      do i = 1,4
	        ii= BCe(n,i)				! BC global node i number
	        T(ii,s) = T(ii,s) + Te(i)		! add second term
	      enddo	! i

	    endif
	  enddo	! n
	enddo	! s

	return
	end


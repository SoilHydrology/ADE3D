	subroutine ERROR(iout, idbg, ires, ijmp, Ne, Nn, Nb, Ns, Npol, Nf, Ndof, &
			 time, dt, theta, theta1, ie, epol, face, C, Sn, Cold, Sn0, BCe, &
			 BCvalue, BCtype, I_Omega,  I_Gamma, I_Jump, R_Omega, J_Gamma)
! element error calculation

	implicit none
	integer iout, idbg, ires, ijmp
	integer Ne, Nn, Nb, Ns				! array parameters
	integer	Npol, Nf, Ndof				! array parameters
	real*8 time, dt
	real*8 theta, theta1
	integer ie(Ne,3)				! global connectivity array
	integer epol(Ne)				! global element polynomial array
	integer face(Nf,8)				! element, local face, node 1 - node 4
	integer BCe(Nb,5)				! BC element numbers
	real*8 BCvalue(Nb,Ns,4)				! BC value (jx_bar, qx_bar or c_bar)
	character*1 BCtype(Nb)				! BC type ('R', 'N' or 'D')
	real*8 C   (Ndof,Ns)				! global  arrays
	real*8 Cold(Ndof,Ns)				! global  arrays
	real*8 Sn (Ndof,Ns)				! global  arrays (nodal S at t^n+1)
	real*8 Sn0(Ndof,Ns)				! global  arrays (nodal S at t^n  )
	real*8 I_Omega(Ne,Npol+7)			! ### error parameters ###
	real*8 I_Gamma(Ne,Npol+7,6)			! ### error parameters ###
	real*8 I_Jump(Nf,8)				! ### error parameters ###
	real*8 R_Omega(Ne)				! ###     residual on element ###
	real*8 J_Gamma(Nf)				! ###     jump at faces ###

	integer i, n, s, e, p, f, m(4)
	integer ii					! global index
!	real*8 Rmax, Jmax

!	write(idbg,'(a)') ' --- ERROR ---'	! ### TEMPORARY ###

	R_Omega  = 0.d0	! use matrix form
	J_Gamma  = 0.d0	! use matrix form
	do s = 1, Ns

! element residual (integrated along the element)
	  do e = 1, Ne
! add to R_Omega 
! NOTE that C_I = 0 at the nodes for higher-than-linear DOFs (I>2)
	    do p = 1, epol(e)+7
	      if(p .le. 8) then
! for p<=8 p is the local DOF
	        ii = ie(e,p)		! global node #
	      else
! for p> 8 p is the order
	        ii = (p-9)*Ne + e + Nn	! global element DOF # for higher order p>2
	      endif
	      R_Omega(e)  = R_Omega(e) + &
				( C   (ii,s)/dt - theta *Sn (ii,s) - &
				( Cold(ii,s)/dt + theta1*Sn0(ii,s)  ) ) * I_Omega(e,p)
! add the element face terms
	      do f = 1,6
	        R_omega(e) = R_Omega(e) + ( theta*C(ii,s) + theta1*Cold(ii,s) ) * I_Gamma(e,p,f) 
	      enddo	! f    
	    enddo	! p
	  enddo		! e
	enddo		! s

! nodal jump between adjacent elements
! NOTE that C_p = 0 at the nodes for higher-than-linear DOFs (p>8)
	do s = 1, Ns
	  do f = 1, Nf
	    e = face(f,1)		! element #
	    do i = 1,4
	      m(i) = face(f,i+2)		! node i of face f
	      J_Gamma(f) = J_Gamma(f) + ( theta*C( m(i), s) + theta1*Cold( m(i) ,s) ) * I_Jump( f, i+2 )
	    enddo	! i
	  enddo		! f

! at Robin BC faces, add the BC flux to J_Gamma at the appropriate face
	  do n = 1, Nb
	    if(BCtype(n) .eq. 'R') then
! for Robin BC
	      do i = 1,4
	        ii= BCe(n,i)		! BC global node i number
	        J_Gamma(ii) = J_Gamma(ii) + BCvalue(n,s,i)   	! add BC to J_Gamma
	      enddo	! i
	    endif
	  enddo		! n
	enddo		! s

! output R_Omega to file resi.txt
	  write(ires,*) 'time = ', time
	  write(ires,*) 'R_Omega(e)'
	do e = 1, Ne
	  write(ires,*)  R_Omega(e)
	enddo		! e
!	Rmax = maxval(abs(R_Omega))	! max{ |R_Omega| }
!	write(ires,*)  Rmax

! output J_Gamma to file jump.txt
	  write(ijmp,*) 'time = ', time
	  write(ijmp,*) 'J_Gamma(f)'
	do f = 1, Nf
	  write(ijmp,*)  J_Gamma(f)
	enddo		! f
!	Jmax = maxval(abs(J_Gamma))	! max{ |J_Gamma| }
!	write(ijmp,*)  Jmax

	return
	end


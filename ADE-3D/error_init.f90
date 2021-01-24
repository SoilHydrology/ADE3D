	subroutine ERROR_INIT(iout, idbg, Ne, Ng, Ns, Npol, Nf, epol, face, &
				V, D, fn, wg, Shp, dNdr, Jac, Ji, &
				aii, I_Omega, I_Gamma, I_Jump)
! element error initialization

	implicit none
	integer iout, idbg
	integer Ne, Ng, Ns				! array parameters
	integer	Npol, Nf				! array parameters
	integer epol(Ne)				! global element polynomial array
	integer face(Nf,8)				! element, local face, node 1 - node 4
	real*8 V(Ne,3), D(Ne,3,3)			! global arrays
	real*8 fn(Ne,6,3)				! element face unit normals
	real*8 wg(Ng)					! Gauss weights
	real*8 aii(3,0:Ng+1,0:Ng+1,0:Ng+1)		! the minors a11, a22, a33 of Jij
	real*8 Ji(3,3,0:Ng+1,0:Ng+1,0:Ng+1), Jac(0:Ng+1,0:Ng+1,0:Ng+1)			! geometric entities
	real*8 Shp(Npol+7,0:Ng+1,0:Ng+1,0:Ng+1), dNdr(Npol+7,3,0:Ng+1,0:Ng+1,0:Ng+1)	! shape and weight functions
	real*8 I_Omega(Ne,Npol+7)			! ### error parameters ###
	real*8 I_Gamma(Ne,Npol+7,6)			! ### error parameters ###
	real*8 I_Jump(Nf,8)				! ### error parameters ###

	integer e, s

	integer g1, g2, g3, g1b, g1t, g2b, g2t, g3b, g3t, i, p, q, jj, n, f, Lf, Gf
	integer e_neigh, f_neigh
	real*8 w, Nu(3), w1, w2, w3

!	write(idbg,'(a)') ' --- ERROR_INIT ---'	! ### TEMPORARY ###

	I_Omega  = 0.d0	! use matrix form
	I_Gamma  = 0.d0	! use matrix form
	I_Jump   = 0.d0	! use matrix form

	do e = 1, Ne

! calculate the element volume integral, I_Omega

	  do g1 = 1, Ng
	    do g2 = 1, Ng
	      do g3 = 1, Ng
	        w = wg(g1) * wg(g2) * wg(g3)
	          do jj = 1, epol(e)+7
	            I_Omega(e,jj) = I_Omega(e,jj) +      w * Shp(jj,g1,g2,g3) * Jac(g1,g2,g3)
	          enddo	! jj
	      enddo	! g3
	    enddo	! g2
	  enddo		! g1

! calculate the nodal flux jump integral, I_Jump

	  do Lf = 1,6
	    do i = 1,3
	      Nu(i) = fn(e,Lf,i)		! element e, face Lf unit normal componet i
	    enddo		! i

	    if     (Lf.eq.1 .or. Lf.eq.4) then
! west   or east  face
	      f = 1
	    else if(Lf.eq.2 .or. Lf.eq.5) then
! south  or north face
	      f = 2
	    else if(Lf.eq.3 .or. Lf.eq.6) then
! bottom or top   face
	      f = 3
	    endif

! calculate integral[n_i*(N_J*v_i - N_J,jDij)]dG
	    if      (Lf .eq. 1) then
	      g1b = 0
	      g1t = 0
	    else if(Lf .eq. 4) then
	      g1b = Ng+1
	      g1t = Ng+1
	    else
	      g1b = 1
	      g1t = Ng
	    endif
	    do g1 = g1b, g1t
	      if (f .eq. 1) then
	        w1 = 1.d0
	      else
	        w1 = wg(g1)
	      endif

	      if     (Lf .eq. 2) then
	        g2b = 0
	        g2t = 0
	      else if(Lf .eq. 5) then
	        g2b = Ng+1
	        g2t = Ng+1
	      else
	        g2b = 1
	        g2t = Ng
	      endif
	      do g2 = g2b, g2t
	        if (f .eq. 2) then
	          w2 = 1.d0
	        else
	          w2 = wg(g2)
	        endif

	        if     (Lf .eq. 3) then
	          g3b = 0
	          g3t = 0
	        else if(Lf .eq. 6) then
	          g3b = Ng+1
	          g3t = Ng+1
	        else
	          g3b = 1
	          g3t = Ng
	        endif
	        do g3 = g3b, g3t
	          if (f .eq. 3) then
	            w3 = 1.d0
	          else
	            w3 = wg(g3)
	          endif

	          w = w1 * w2 * w3

	          do jj = 1, epol(e)+7
	            do i = 1,3
	              I_Gamma(e,jj,Lf) = I_Gamma(e,jj,Lf) +				&
				V(e,i) * Shp(jj,g1,g2,g3) *				&
				w * Nu(i) * aii(f,g1,g2,g3)
	              do p = 1,3
	                do q = 1,3
	                  I_Gamma(e,jj,Lf) = I_Gamma(e,jj,Lf) - 			&
			        D(e,i,p) * dNdr(jj,q,g1,g2,g3) * Ji(q,p,g1,g2,g3) *	&
				w * Nu(i) * aii(f,g1,g2,g3)
	                enddo	! q
	              enddo	! p
	            enddo	! i
	          enddo		! jj
	        enddo		! g3
	      enddo		! g2
	    enddo		! g1
	  enddo			! Lf
	enddo			! e

! for jj<=8 calculate nodal jumps
	do Gf = 1, Nf
	  e	  = face(Gf,1)	! element #
	  Lf	  = face(Gf,2)	! local face #
	  e_neigh = face(Gf,7)	! neighbour element #
	  f_neigh = face(Gf,8)	! local neighbour face #
	  do jj = 1,8
	    I_Jump(Gf,jj) = I_Jump(Gf,jj) + I_Gamma(e,jj,Lf)
	    if(e_neigh.ge.1 .and. e_neigh.le. Ne)	then
! if the neighbour element exists, add its I_Gamma
	      I_Jump(Gf,jj) = I_Jump(Gf,jj) + I_Gamma(e_neigh,jj,f_neigh)
	    endif
	  enddo	! jj
	enddo	! Gf

	return
	end

	subroutine SHAPE(iout, idbg, Ne, Nn, Ng, Npol, V, ie, x, xg, e, &
			  Shp, dNdr, Wgt, J, Jac, Ji, Jaci, aii, aopt, Pn, dPn)
! calculate element shape functions

	implicit none
	integer iout, idbg
	integer Ne, Nn, Ng			! array parameters
	integer	Npol				! array parameters
	integer ie(Ne,9)			! global connectivity array
	real*8 x(Nn,3)				! global coordinates array
	real*8 V(Ne,3)				! global  arrays
	real*8 xg(Ng)				! Gauss abscissas [-1,+1]
	real*8 J(3,3,0:Ng+1,0:Ng+1,0:Ng+1), Ji(3,3,0:Ng+1,0:Ng+1,0:Ng+1), &
		Jac(0:Ng+1,0:Ng+1,0:Ng+1), Jaci(0:Ng+1,0:Ng+1,0:Ng+1)	! geometric entities
	real*8 aii(3,0:Ng+1,0:Ng+1,0:Ng+1)	! the minors a11, a22, a33 of Jij
	real*8 Shp(Npol+7,0:Ng+1,0:Ng+1,0:Ng+1), dNdr(Npol+7,3,0:Ng+1,0:Ng+1,0:Ng+1), &
		Wgt(Npol+7,0:Ng+1,0:Ng+1,0:Ng+1)	! shape and weight functions
	real*8 Pn(0:Npol,0:Ng+1), dPn(0:Npol,0:Ng+1)	! Legendre polynomial values and derivatives vectors (output)
	real*8 aopt				! SUPG alpha_opt
	integer e

	integer i1, i2, i3, i4, i5, i6, i7, i8, g1, g2, g3, ierror
	integer i, m, k, p, q, idx
	real*8 r, s, t, fact
	real*8 h17, h28, h35, h46, h, Vabs

	data ierror /0/

!	write(idbg,'(a)') ' --- SHAPE ---'	! ### TEMPORARY ###

	i1 = ie(e,1)			! 1st node
	i2 = ie(e,2)			! 2nd node
	i3 = ie(e,3)			! 3rd node
	i4 = ie(e,4)			! 4th node
	i5 = ie(e,5)			! 5th node
	i6 = ie(e,6)			! 6th node
	i7 = ie(e,7)			! 7th node
	i8 = ie(e,8)			! 8th node
! assume the cell size, h, is the larger main diagonal
	h17 =   ( x(ie(e,7),1) - x(ie(e,1),1) )**2 + &
		( x(ie(e,7),2) - x(ie(e,1),2) )**2 + &
		( x(ie(e,7),3) - x(ie(e,1),3) )**2			! diagonal 1->7 ^2
	h28 =   ( x(ie(e,8),1) - x(ie(e,2),1) )**2 + &
		( x(ie(e,8),2) - x(ie(e,2),2) )**2 + &
		( x(ie(e,8),3) - x(ie(e,2),3) )**2			! diagonal 2->8 ^2
	h35 =   ( x(ie(e,5),1) - x(ie(e,3),1) )**2 + &
		( x(ie(e,5),2) - x(ie(e,3),2) )**2 + &
		( x(ie(e,5),3) - x(ie(e,3),3) )**2			! diagonal 3->5 ^2
	h46 =   ( x(ie(e,6),1) - x(ie(e,4),1) )**2 + &
		( x(ie(e,6),2) - x(ie(e,4),2) )**2 + &
		( x(ie(e,6),3) - x(ie(e,4),3) )**2			! diagonal 4->6 ^2


	h = sqrt( max( h17, h28, h35, h46 ) )				! element size
	Vabs = sqrt(  V(e,1)**2 + V(e,2)**2 + V(e,3)**2 )		! |Vi|

	do g1 = 0, Ng+1
	  if      (g1 .eq. 0   ) then
	    r = -1.d0
	  else if (g1 .eq. Ng+1) then
	    r =  1.d0
	  else
	    r = xg(g1)
	  endif

	  do g2 = 0, Ng+1
	    if      (g2 .eq. 0   ) then
	      s = -1.d0
	    else if (g2 .eq. Ng+1) then
	      s =  1.d0
	    else
	      s = xg(g2)
	    endif

	    do g3 = 0, Ng+1
	      if      (g3 .eq. 0   ) then
	        t = -1.d0
	      else if (g3 .eq. Ng+1) then
	        t =  1.d0
	      else
	        t = xg(g3)
	      endif

! linear 3D shape functions
! -1 < r, s, t < +1
! Ni(r,s,t) = (1 +/- r)(1 +/- s)(1 +/- t) / 8
	      Shp(1,g1,g2,g3) = 0.125d0 * (1.-r) * (1.-s) * (1.-t)	! N1(r,s,t)
	      Shp(2,g1,g2,g3) = 0.125d0 * (1.+r) * (1.-s) * (1.-t)	! N2(r,s,t)
	      Shp(3,g1,g2,g3) = 0.125d0 * (1.+r) * (1.+s) * (1.-t)	! N3(r,s,t)
	      Shp(4,g1,g2,g3) = 0.125d0 * (1.-r) * (1.+s) * (1.-t)	! N4(r,s,t)
	      Shp(5,g1,g2,g3) = 0.125d0 * (1.-r) * (1.-s) * (1.+t)	! N5(r,s,t)
	      Shp(6,g1,g2,g3) = 0.125d0 * (1.+r) * (1.-s) * (1.+t)	! N6(r,s,t)
	      Shp(7,g1,g2,g3) = 0.125d0 * (1.+r) * (1.+s) * (1.+t)	! N7(r,s,t)
	      Shp(8,g1,g2,g3) = 0.125d0 * (1.-r) * (1.+s) * (1.+t)	! N8(r,s,t)

	      dNdr(1,1,g1,g2,g3) =-0.125d0 * (1.-s) * (1.-t)		! dN1(r,s,t)/dr
	      dNdr(2,1,g1,g2,g3) = 0.125d0 * (1.-s) * (1.-t)		! dN2(r,s,t)/dr
	      dNdr(3,1,g1,g2,g3) = 0.125d0 * (1.+s) * (1.-t)		! dN3(r,s,t)/dr
	      dNdr(4,1,g1,g2,g3) =-0.125d0 * (1.+s) * (1.-t)		! dN4(r,s,t)/dr
	      dNdr(5,1,g1,g2,g3) =-0.125d0 * (1.-s) * (1.+t)		! dN5(r,s,t)/dr
	      dNdr(6,1,g1,g2,g3) = 0.125d0 * (1.-s) * (1.+t)		! dN6(r,s,t)/dr
	      dNdr(7,1,g1,g2,g3) = 0.125d0 * (1.+s) * (1.+t)		! dN7(r,s,t)/dr
	      dNdr(8,1,g1,g2,g3) =-0.125d0 * (1.+s) * (1.+t)		! dN8(r,s,t)/dr

	      dNdr(1,2,g1,g2,g3) =-0.125d0 * (1.-r) * (1.-t)		! dN1(r,s,t)/ds
	      dNdr(2,2,g1,g2,g3) =-0.125d0 * (1.+r) * (1.-t)		! dN2(r,s,t)/ds
	      dNdr(3,2,g1,g2,g3) = 0.125d0 * (1.+r) * (1.-t)		! dN3(r,s,t)/ds
	      dNdr(4,2,g1,g2,g3) = 0.125d0 * (1.-r) * (1.-t)		! dN4(r,s,t)/ds
	      dNdr(5,2,g1,g2,g3) =-0.125d0 * (1.-r) * (1.+t)		! dN5(r,s,t)/ds
	      dNdr(6,2,g1,g2,g3) =-0.125d0 * (1.+r) * (1.+t)		! dN6(r,s,t)/ds
	      dNdr(7,2,g1,g2,g3) = 0.125d0 * (1.+r) * (1.+t)		! dN7(r,s,t)/ds
	      dNdr(8,2,g1,g2,g3) = 0.125d0 * (1.-r) * (1.+t)		! dN8(r,s,t)/ds

	      dNdr(1,3,g1,g2,g3) =-0.125d0 * (1.-r) * (1.-s)		! dN1(r,s,t)/dt
	      dNdr(2,3,g1,g2,g3) =-0.125d0 * (1.+r) * (1.-s)		! dN2(r,s,t)/dt
	      dNdr(3,3,g1,g2,g3) =-0.125d0 * (1.+r) * (1.+s)		! dN3(r,s,t)/dt
	      dNdr(4,3,g1,g2,g3) =-0.125d0 * (1.-r) * (1.+s)		! dN4(r,s,t)/dt
	      dNdr(5,3,g1,g2,g3) = 0.125d0 * (1.-r) * (1.-s)		! dN5(r,s,t)/dt
	      dNdr(6,3,g1,g2,g3) = 0.125d0 * (1.+r) * (1.-s)		! dN6(r,s,t)/dt
	      dNdr(7,3,g1,g2,g3) = 0.125d0 * (1.+r) * (1.+s)		! dN7(r,s,t)/dt
	      dNdr(8,3,g1,g2,g3) = 0.125d0 * (1.-r) * (1.+s)		! dN8(r,s,t)/dt

! Jij is the Jacobian matrix
	      J(1,1,g1,g2,g3) = dNdr(1,1,g1,g2,g3)*x(i1,1) + dNdr(2,1,g1,g2,g3)*x(i2,1) + &
			        dNdr(3,1,g1,g2,g3)*x(i3,1) + dNdr(4,1,g1,g2,g3)*x(i4,1) + &
			        dNdr(5,1,g1,g2,g3)*x(i5,1) + dNdr(6,1,g1,g2,g3)*x(i6,1) + &
			        dNdr(7,1,g1,g2,g3)*x(i7,1) + dNdr(8,1,g1,g2,g3)*x(i8,1)	! x,r
	      J(1,2,g1,g2,g3) = dNdr(1,2,g1,g2,g3)*x(i1,1) + dNdr(2,2,g1,g2,g3)*x(i2,1) + &
		     	        dNdr(3,2,g1,g2,g3)*x(i3,1) + dNdr(4,2,g1,g2,g3)*x(i4,1)	 + &
		     	        dNdr(5,2,g1,g2,g3)*x(i5,1) + dNdr(6,2,g1,g2,g3)*x(i6,1)	 + &
		     	        dNdr(7,2,g1,g2,g3)*x(i7,1) + dNdr(8,2,g1,g2,g3)*x(i8,1)	! x,s
	      J(1,3,g1,g2,g3) = dNdr(1,3,g1,g2,g3)*x(i1,1) + dNdr(2,3,g1,g2,g3)*x(i2,1) + &
		     	        dNdr(3,3,g1,g2,g3)*x(i3,1) + dNdr(4,3,g1,g2,g3)*x(i4,1)	 + &
		     	        dNdr(5,3,g1,g2,g3)*x(i5,1) + dNdr(6,3,g1,g2,g3)*x(i6,1)	 + &
		     	        dNdr(7,3,g1,g2,g3)*x(i7,1) + dNdr(8,3,g1,g2,g3)*x(i8,1)	! x,t

	      J(2,1,g1,g2,g3) = dNdr(1,1,g1,g2,g3)*x(i1,2) + dNdr(2,1,g1,g2,g3)*x(i2,2) + &
		     	        dNdr(3,1,g1,g2,g3)*x(i3,2) + dNdr(4,1,g1,g2,g3)*x(i4,2) + &
		     	        dNdr(5,1,g1,g2,g3)*x(i5,2) + dNdr(6,1,g1,g2,g3)*x(i6,2) + &
		     	        dNdr(7,1,g1,g2,g3)*x(i7,2) + dNdr(8,1,g1,g2,g3)*x(i8,2)	! y,r
	      J(2,2,g1,g2,g3) = dNdr(1,2,g1,g2,g3)*x(i1,2) + dNdr(2,2,g1,g2,g3)*x(i2,2) + &
		     	        dNdr(3,2,g1,g2,g3)*x(i3,2) + dNdr(4,2,g1,g2,g3)*x(i4,2) + &
		     	        dNdr(5,2,g1,g2,g3)*x(i5,2) + dNdr(6,2,g1,g2,g3)*x(i6,2) + &
		     	        dNdr(7,2,g1,g2,g3)*x(i7,2) + dNdr(8,2,g1,g2,g3)*x(i8,2)	! y,s
	      J(2,3,g1,g2,g3) = dNdr(1,3,g1,g2,g3)*x(i1,2) + dNdr(2,3,g1,g2,g3)*x(i2,2) + &
		     	        dNdr(3,3,g1,g2,g3)*x(i3,2) + dNdr(4,3,g1,g2,g3)*x(i4,2) + &
		     	        dNdr(5,3,g1,g2,g3)*x(i5,2) + dNdr(6,3,g1,g2,g3)*x(i6,2) + &
		     	        dNdr(7,3,g1,g2,g3)*x(i7,2) + dNdr(8,3,g1,g2,g3)*x(i8,2)	! y,t

	      J(3,1,g1,g2,g3) = dNdr(1,1,g1,g2,g3)*x(i1,3) + dNdr(2,1,g1,g2,g3)*x(i2,3) + &
			        dNdr(3,1,g1,g2,g3)*x(i3,3) + dNdr(4,1,g1,g2,g3)*x(i4,3) + &
			        dNdr(5,1,g1,g2,g3)*x(i5,3) + dNdr(6,1,g1,g2,g3)*x(i6,3) + &
			        dNdr(7,1,g1,g2,g3)*x(i7,3) + dNdr(8,1,g1,g2,g3)*x(i8,3)	! z,r
	      J(3,2,g1,g2,g3) = dNdr(1,2,g1,g2,g3)*x(i1,3) + dNdr(2,2,g1,g2,g3)*x(i2,3) + &
		     	        dNdr(3,2,g1,g2,g3)*x(i3,3) + dNdr(4,2,g1,g2,g3)*x(i4,3)	 + &
		     	        dNdr(5,2,g1,g2,g3)*x(i5,3) + dNdr(6,2,g1,g2,g3)*x(i6,3)	 + &
		     	        dNdr(7,2,g1,g2,g3)*x(i7,3) + dNdr(8,2,g1,g2,g3)*x(i8,3)	! z,s
	      J(3,3,g1,g2,g3) = dNdr(1,3,g1,g2,g3)*x(i1,3) + dNdr(2,3,g1,g2,g3)*x(i2,3) + &
		     	        dNdr(3,3,g1,g2,g3)*x(i3,3) + dNdr(4,3,g1,g2,g3)*x(i4,3)	 + &
		     	        dNdr(5,3,g1,g2,g3)*x(i5,3) + dNdr(6,3,g1,g2,g3)*x(i6,3)	 + &
		     	        dNdr(7,3,g1,g2,g3)*x(i7,3) + dNdr(8,3,g1,g2,g3)*x(i8,3)	! z,t

! aii are the minors a11, a22, a33 of Jij
	      aii(1,g1,g2,g3) = J(2,2,g1,g2,g3)*J(3,3,g1,g2,g3) - J(2,3,g1,g2,g3)*J(3,2,g1,g2,g3)	
	      aii(2,g1,g2,g3) = J(3,3,g1,g2,g3)*J(1,1,g1,g2,g3) - J(3,1,g1,g2,g3)*J(1,3,g1,g2,g3)	
	      aii(3,g1,g2,g3) = J(1,1,g1,g2,g3)*J(2,2,g1,g2,g3) - J(1,2,g1,g2,g3)*J(2,1,g1,g2,g3)	

	      Jac (g1,g2,g3) = J(1,1,g1,g2,g3) * ( J(2,2,g1,g2,g3)*J(3,3,g1,g2,g3)   - &
						   J(2,3,g1,g2,g3)*J(3,2,g1,g2,g3) ) - &
			       J(1,2,g1,g2,g3) * ( J(2,1,g1,g2,g3)*J(3,3,g1,g2,g3)   - &
						   J(2,3,g1,g2,g3)*J(3,1,g1,g2,g3) ) + &
			       J(1,3,g1,g2,g3) * ( J(2,1,g1,g2,g3)*J(3,2,g1,g2,g3)   - &
						   J(2,2,g1,g2,g3)*J(3,1,g1,g2,g3) ) ! Jacobian determinant, |J|

	      Jaci(g1,g2,g3) = 1. / Jac(g1,g2,g3)		! 1/|J|

! inv(Jij)
	      Ji(1,1,g1,g2,g3) = ( -J(2,3,g1,g2,g3)*J(3,2,g1,g2,g3) + &
				    J(2,2,g1,g2,g3)*J(3,3,g1,g2,g3) ) *Jaci(g1,g2,g3)
	      Ji(1,2,g1,g2,g3) = ( +J(1,3,g1,g2,g3)*J(3,2,g1,g2,g3) - &
				    J(1,2,g1,g2,g3)*J(3,3,g1,g2,g3) ) *Jaci(g1,g2,g3)
	      Ji(1,3,g1,g2,g3) = ( -J(1,3,g1,g2,g3)*J(2,2,g1,g2,g3) + &
				    J(1,2,g1,g2,g3)*J(2,3,g1,g2,g3) ) *Jaci(g1,g2,g3)

	      Ji(2,1,g1,g2,g3) = (  J(2,3,g1,g2,g3)*J(3,1,g1,g2,g3) - &
				    J(2,1,g1,g2,g3)*J(3,3,g1,g2,g3) ) *Jaci(g1,g2,g3)
	      Ji(2,2,g1,g2,g3) = ( -J(1,3,g1,g2,g3)*J(3,1,g1,g2,g3) + &
				    J(1,1,g1,g2,g3)*J(3,3,g1,g2,g3) ) *Jaci(g1,g2,g3)
	      Ji(2,3,g1,g2,g3) = ( +J(1,3,g1,g2,g3)*J(2,1,g1,g2,g3) - &
				    J(1,1,g1,g2,g3)*J(2,3,g1,g2,g3) ) *Jaci(g1,g2,g3)

	      Ji(3,1,g1,g2,g3) = ( -J(2,2,g1,g2,g3)*J(3,1,g1,g2,g3) + &
				    J(2,1,g1,g2,g3)*J(3,2,g1,g2,g3) ) *Jaci(g1,g2,g3)
	      Ji(3,2,g1,g2,g3) = ( +J(2,1,g1,g2,g3)*J(3,1,g1,g2,g3) - &
				    J(1,1,g1,g2,g3)*J(3,2,g1,g2,g3) ) *Jaci(g1,g2,g3)
	      Ji(3,3,g1,g2,g3) = ( -J(1,2,g1,g2,g3)*J(2,1,g1,g2,g3) + &
				    J(1,1,g1,g2,g3)*J(2,2,g1,g2,g3) ) *Jaci(g1,g2,g3)

! Wi(r,s,t)
	      Wgt(1,g1,g2,g3) =  Shp(1,g1,g2,g3) + 0.5d0 * aopt * h / Vabs * &
			    ( dNdr(1,1,g1,g2,g3)*V(e,1)*Ji(1,1,g1,g2,g3) + &
			      dNdr(1,2,g1,g2,g3)*V(e,1)*Ji(2,1,g1,g2,g3) + &
			      dNdr(1,3,g1,g2,g3)*V(e,1)*Ji(3,1,g1,g2,g3) + &
			      dNdr(1,1,g1,g2,g3)*V(e,2)*Ji(1,2,g1,g2,g3) + &
			      dNdr(1,2,g1,g2,g3)*V(e,2)*Ji(2,2,g1,g2,g3) + &
			      dNdr(1,3,g1,g2,g3)*V(e,2)*Ji(3,2,g1,g2,g3) + &
			      dNdr(1,1,g1,g2,g3)*V(e,3)*Ji(1,3,g1,g2,g3) + &
			      dNdr(1,2,g1,g2,g3)*V(e,3)*Ji(2,3,g1,g2,g3) + &
			      dNdr(1,3,g1,g2,g3)*V(e,3)*Ji(3,3,g1,g2,g3) )	! W1(r,s,t)


	      Wgt(2,g1,g2,g3) =  Shp(2,g1,g2,g3) + 0.5d0 * aopt * h / Vabs * &
			    ( dNdr(2,1,g1,g2,g3)*V(e,1)*Ji(1,1,g1,g2,g3) + &
			      dNdr(2,2,g1,g2,g3)*V(e,1)*Ji(2,1,g1,g2,g3) + &
			      dNdr(2,3,g1,g2,g3)*V(e,1)*Ji(3,1,g1,g2,g3) + &
			      dNdr(2,1,g1,g2,g3)*V(e,2)*Ji(1,2,g1,g2,g3) + &
			      dNdr(2,2,g1,g2,g3)*V(e,2)*Ji(2,2,g1,g2,g3) + &
			      dNdr(2,3,g1,g2,g3)*V(e,2)*Ji(3,2,g1,g2,g3) + &
			      dNdr(2,1,g1,g2,g3)*V(e,3)*Ji(1,3,g1,g2,g3) + &
			      dNdr(2,2,g1,g2,g3)*V(e,3)*Ji(2,3,g1,g2,g3) + &
			      dNdr(2,3,g1,g2,g3)*V(e,3)*Ji(3,3,g1,g2,g3) )	! W2(r,s,t)

	      Wgt(3,g1,g2,g3) =  Shp(3,g1,g2,g3) + 0.5d0 * aopt * h / Vabs * &
			    ( dNdr(3,1,g1,g2,g3)*V(e,1)*Ji(1,1,g1,g2,g3) + &
			      dNdr(3,2,g1,g2,g3)*V(e,1)*Ji(2,1,g1,g2,g3) + &
			      dNdr(3,3,g1,g2,g3)*V(e,1)*Ji(3,1,g1,g2,g3) + &
			      dNdr(3,1,g1,g2,g3)*V(e,2)*Ji(1,2,g1,g2,g3) + &
			      dNdr(3,2,g1,g2,g3)*V(e,2)*Ji(2,2,g1,g2,g3) + &
			      dNdr(3,3,g1,g2,g3)*V(e,2)*Ji(3,2,g1,g2,g3) + &
			      dNdr(3,1,g1,g2,g3)*V(e,3)*Ji(1,3,g1,g2,g3) + &
			      dNdr(3,2,g1,g2,g3)*V(e,3)*Ji(2,3,g1,g2,g3) + &
			      dNdr(3,3,g1,g2,g3)*V(e,3)*Ji(3,3,g1,g2,g3) )	! W3(r,s,t)

	      Wgt(4,g1,g2,g3) =  Shp(4,g1,g2,g3) + 0.5d0 * aopt * h / Vabs * &
			    ( dNdr(4,1,g1,g2,g3)*V(e,1)*Ji(1,1,g1,g2,g3) + &
			      dNdr(4,2,g1,g2,g3)*V(e,1)*Ji(2,1,g1,g2,g3) + &
			      dNdr(4,3,g1,g2,g3)*V(e,1)*Ji(3,1,g1,g2,g3) + &
			      dNdr(4,1,g1,g2,g3)*V(e,2)*Ji(1,2,g1,g2,g3) + &
			      dNdr(4,2,g1,g2,g3)*V(e,2)*Ji(2,2,g1,g2,g3) + &
			      dNdr(4,3,g1,g2,g3)*V(e,2)*Ji(3,2,g1,g2,g3) + &
			      dNdr(4,1,g1,g2,g3)*V(e,3)*Ji(1,3,g1,g2,g3) + &
			      dNdr(4,2,g1,g2,g3)*V(e,3)*Ji(2,3,g1,g2,g3) + &
			      dNdr(4,3,g1,g2,g3)*V(e,3)*Ji(3,3,g1,g2,g3) )	! W4(r,s,t)

	      Wgt(5,g1,g2,g3) =  Shp(5,g1,g2,g3) + 0.5d0 * aopt * h / Vabs * &
			    ( dNdr(5,1,g1,g2,g3)*V(e,1)*Ji(1,1,g1,g2,g3) + &
			      dNdr(5,2,g1,g2,g3)*V(e,1)*Ji(2,1,g1,g2,g3) + &
			      dNdr(5,3,g1,g2,g3)*V(e,1)*Ji(3,1,g1,g2,g3) + &
			      dNdr(5,1,g1,g2,g3)*V(e,2)*Ji(1,2,g1,g2,g3) + &
			      dNdr(5,2,g1,g2,g3)*V(e,2)*Ji(2,2,g1,g2,g3) + &
			      dNdr(5,3,g1,g2,g3)*V(e,2)*Ji(3,2,g1,g2,g3) + &
			      dNdr(5,1,g1,g2,g3)*V(e,3)*Ji(1,3,g1,g2,g3) + &
			      dNdr(5,2,g1,g2,g3)*V(e,3)*Ji(2,3,g1,g2,g3) + &
			      dNdr(5,3,g1,g2,g3)*V(e,3)*Ji(3,3,g1,g2,g3) )	! W5(r,s,t)


	      Wgt(6,g1,g2,g3) =  Shp(6,g1,g2,g3) + 0.5d0 * aopt * h / Vabs * &
			    ( dNdr(6,1,g1,g2,g3)*V(e,1)*Ji(1,1,g1,g2,g3) + &
			      dNdr(6,2,g1,g2,g3)*V(e,1)*Ji(2,1,g1,g2,g3) + &
			      dNdr(6,3,g1,g2,g3)*V(e,1)*Ji(3,1,g1,g2,g3) + &
			      dNdr(6,1,g1,g2,g3)*V(e,2)*Ji(1,2,g1,g2,g3) + &
			      dNdr(6,2,g1,g2,g3)*V(e,2)*Ji(2,2,g1,g2,g3) + &
			      dNdr(6,3,g1,g2,g3)*V(e,2)*Ji(3,2,g1,g2,g3) + &
			      dNdr(6,1,g1,g2,g3)*V(e,3)*Ji(1,3,g1,g2,g3) + &
			      dNdr(6,2,g1,g2,g3)*V(e,3)*Ji(2,3,g1,g2,g3) + &
			      dNdr(6,3,g1,g2,g3)*V(e,3)*Ji(3,3,g1,g2,g3) )	! W6(r,s,t)


	      Wgt(7,g1,g2,g3) =  Shp(7,g1,g2,g3) + 0.5d0 * aopt * h / Vabs * &
			    ( dNdr(7,1,g1,g2,g3)*V(e,1)*Ji(1,1,g1,g2,g3) + &
			      dNdr(7,2,g1,g2,g3)*V(e,1)*Ji(2,1,g1,g2,g3) + &
			      dNdr(7,3,g1,g2,g3)*V(e,1)*Ji(3,1,g1,g2,g3) + &
			      dNdr(7,1,g1,g2,g3)*V(e,2)*Ji(1,2,g1,g2,g3) + &
			      dNdr(7,2,g1,g2,g3)*V(e,2)*Ji(2,2,g1,g2,g3) + &
			      dNdr(7,3,g1,g2,g3)*V(e,2)*Ji(3,2,g1,g2,g3) + &
			      dNdr(7,1,g1,g2,g3)*V(e,3)*Ji(1,3,g1,g2,g3) + &
			      dNdr(7,2,g1,g2,g3)*V(e,3)*Ji(2,3,g1,g2,g3) + &
			      dNdr(7,3,g1,g2,g3)*V(e,3)*Ji(3,3,g1,g2,g3) )	! W7(r,s,t)


	      Wgt(8,g1,g2,g3) =  Shp(8,g1,g2,g3) + 0.5d0 * aopt * h / Vabs * &
			    ( dNdr(8,1,g1,g2,g3)*V(e,1)*Ji(1,1,g1,g2,g3) + &
			      dNdr(8,2,g1,g2,g3)*V(e,1)*Ji(2,1,g1,g2,g3) + &
			      dNdr(8,3,g1,g2,g3)*V(e,1)*Ji(3,1,g1,g2,g3) + &
			      dNdr(8,1,g1,g2,g3)*V(e,2)*Ji(1,2,g1,g2,g3) + &
			      dNdr(8,2,g1,g2,g3)*V(e,2)*Ji(2,2,g1,g2,g3) + &
			      dNdr(8,3,g1,g2,g3)*V(e,2)*Ji(3,2,g1,g2,g3) + &
			      dNdr(8,1,g1,g2,g3)*V(e,3)*Ji(1,3,g1,g2,g3) + &
			      dNdr(8,2,g1,g2,g3)*V(e,3)*Ji(2,3,g1,g2,g3) + &
			      dNdr(8,3,g1,g2,g3)*V(e,3)*Ji(3,3,g1,g2,g3) )	! W8(r,s,t)

! see Szabo & Babuska's book, 2011, Eq. (2.68)
	    fact = 1.d0								! reset    factor
	    idx = 8								! reset    index
	    do i = 3, Npol+1
	      fact     = ( 2 * (2*i-3) ) ** (-0.5)				! i facor
	      do m = 3, Npol+1
	        fact   = ( 2 * (2*m-3) ) ** (-0.5) * fact			! m facor * i facor
	        do k = 3, Npol+1
	          fact = ( 2 * (2*k-3) ) ** (-0.5) * fact			! k facor * m facor * i facor
		  idx = idx + 1							! increment index
	          Shp (idx,  g1,g2,g3) = ( Pn (i-1,g1) -  Pn(i-3,g1) ) * &
	          			 ( Pn (m-1,g2) -  Pn(m-3,g2) ) * &
	          			 ( Pn (k-1,g3) -  Pn(k-3,g3) ) * fact	!  N(r,s,t)
		  do p = 1,3							! p component
	            dNdr(idx,p,g1,g2,g3) = (dPn (i-1,g1) - dPn(i-3,g1) ) * &
	          			   (dPn (m-1,g2) - dPn(m-3,g2) ) * &
	          			   (dPn (k-1,g3) - dPn(k-3,g3) ) * fact	! dN(r,s,t)/dxi_p
		    do q = 1,3							! q component
	              Wgt (idx,  g1,g2,g3) =  Shp(idx,g1,g2,g3) + 0.5d0 * aopt * h / Vabs * &
			       ( dNdr(idx,q,g1,g2,g3)*V(e,p)*Ji(q,p,g1,g2,g3) )	! Wi(r,s,t)
		    enddo	! q
		  enddo		! p

	        enddo		! k
	      enddo		! m
	    enddo		! i

! check element geometry
	    if(Jac(g1,g2,g3) .le. 0.) then
	      write(iout,*) '*** ABORT: Jac(g1,g2,g3) <= 0, e, g1, g2, g3, Jac(g1,g2,g3) = ', &
							 e ,g1, g2 ,g3, Jac(g1,g2,g3)
	      ierror = ierror + 1
	    endif

	    enddo	! g3
	  enddo		! g2
	enddo		! g1

	if(ierror .ne. 0) then
	  write(iout,*) '*** ABORT: ierror = ', ierror
	  stop
	endif
	
	return
	end

	subroutine SHAPEM(iout, idbg, Ne, Ng, Npol, epol, V, D, &
			 Ae, Le, Be, wg, e, Shp, dNdr, Wgt, Jac, Ji)
! calculate element matrices

	implicit none
	integer iout, idbg
	integer Ne, Ng				! array parameters
	integer Npol				! array parameters
	real*8 Ae(Npol+7,Npol+7), Le(Npol+7,Npol+7), Be(Npol+7,Npol+7)	! element arrays
	integer epol(Ne)			! global element polynomial array
	real*8 V(Ne,3), D(Ne,3,3)		! global  arrays
	real*8 wg(Ng)				! Gauss weights
	real*8 Ji(3,3,0:Ng+1,0:Ng+1,0:Ng+1), Jac(0:Ng+1,0:Ng+1,0:Ng+1)			! geometric entities
	real*8 Shp(Npol+7,0:Ng+1,0:Ng+1,0:Ng+1), dNdr(Npol+7,3,0:Ng+1,0:Ng+1,0:Ng+1), Wgt(Npol+7,0:Ng+1,0:Ng+1,0:Ng+1)	! shape and weight functions
	integer e

	integer g1, g2, g3, i, p, k, q, ii, jj
	real*8 w

!	write(idbg,'(a)') ' --- SHAPEM ---'	! ### TEMPORARY ###

	Ae = 0.d0	! use matrix form
	Le = 0.d0	! use matrix form
	Be = 0.d0	! use matrix form

	do g1 = 1, Ng
	  do g2 = 1, Ng
	    do g3 = 1, Ng
	      w = wg(g1)*wg(g2)*wg(g3)

! calculate the mass matrix, Ae
! consistent mass matrix:
	      do ii = 1, epol(e)+7
	        do jj = 1, epol(e)+7
	          Ae(ii,jj) = Ae(ii,jj) + w * Wgt(ii,g1,g2,g3) * Shp(jj,g1,g2,g3) * Jac(g1,g2,g3)
	        enddo	! jj
	      enddo	! ii

! calculate the advection matrix, Le (Non-Symmetric)
	      do i = 1,3
	        do k = 1,3
	          do ii = 1, epol(e)+7
	            do jj = 1, epol(e)+7
		      Le(ii,jj) = Le(ii,jj) - w * Wgt(ii,g1,g2,g3) * V(e,i) * &
				  dNdr(jj,k,g1,g2,g3) * Ji(k,i,g1,g2,g3) * Jac(g1,g2,g3)
	            enddo	! jj
	          enddo		! ii
	        enddo		! k
	      enddo		! i

! calculate the conduction matrix, Be
	      do i = 1,3
	        do p = 1,3
	          do k = 1,3
	            do q = 1,3
	              do ii = 1, epol(e)+7
	                do jj = 1, epol(e)+7
!!! dN/dr should be replaced by dW/dr in the following line !!!
			  Be(ii,jj) = Be(ii,jj) + w * dNdr(ii,k,g1,g2,g3) * Ji(k,i,g1,g2,g3) * &
				      D(e,i,p) * dNdr(jj,q,g1,g2,g3) * Ji(q,p,g1,g2,g3) * Jac(g1,g2,g3)
	                enddo	! jj
	              enddo	! ii
	            enddo	! q
	          enddo		! k
	        enddo		! p
	      enddo		! i

	    enddo	! g3
	  enddo		! g2
	enddo		! g1
	
	return
	end

	subroutine SHAPEQ(iout, idbg, Ne, Nn, Ng, Npol, V, D, dNdr, Ji, &
			 Qce1, Qde1, Qce2, Qde2, Qce3, Qde3, epol, e)
! calculate nodal flux matrices

! NOTE that N here is actually the mapping functions Qk inSzabo's 2001 book,
! Eq. (2.60), which are the LINEAR shape functions N1, N2

	implicit none
	integer iout, idbg
	integer Ne, Nn, Ng			! array parameters
	integer Npol				! array parameters
	real*8 Ji(3,3,0:Ng+1,0:Ng+1,0:Ng+1)	! geometric entities
	real*8 dNdr(Npol+7,3,0:Ng+1,0:Ng+1,0:Ng+1)		! shape and weight functions
	real*8 Qce1(Npol+7,Npol+7), Qde1(Npol+7,Npol+7)		! element arrays
	real*8 Qce2(Npol+7,Npol+7), Qde2(Npol+7,Npol+7)		! element arrays
	real*8 Qce3(Npol+7,Npol+7), Qde3(Npol+7,Npol+7)		! element arrays
	integer epol(Ne)			! global element polynomial array
	real*8 V(Ne,3), D(Ne,3,3)		! global  arrays
	integer e

	integer i, j, k, m, p, ii, jj, g1, g2, g3

!	write(idbg,'(a)') ' --- SHAPEQ ---'	! ### TEMPORARY ###

! calculate the nodal advection flux matrix, Qce
	Qce1 = 0.		! reset Qce1
	Qce2 = 0.		! reset Qce2
	Qce3 = 0.		! reset Qce3
! the shape functions for i=9:epol(e) are 0 at the element boundaries
! but their derivaties are not

!	do ii = 1, epol(e)+7
	do ii = 1, 8
	  Qce1(ii,ii) = V(e,1)
	  Qce2(ii,ii) = V(e,2)
	  Qce3(ii,ii) = V(e,3)
	enddo			! ii	

! calculate the nodal conduction flux matrix, Qde
	Qde1 = 0.		! reset Qde1
	Qde2 = 0.		! reset Qde2
	Qde3 = 0.		! reset Qde3

!	do ii = 1, epol(e)+7
	do ii = 1, 8
!	  do jj = 1, epol(e)+7
	  do jj = 1, 8
	    do i = 1,2
! loop on endpoints (i=1,2 at r=-1, +1, or Gauss point 0, Ng+1, resp.)
	      if     (i .eq. 1) then
	        g1 = 0
	      else
	        g1 = Ng+1
	      endif
	      do j = 1,2
! loop on endpoints (j=1,2 at s=-1, +1, or Gauss point 0, Ng+1, resp.)
	        if     (jj .eq. 1) then
	          g2 = 0
	        else
	          g2 = Ng+1
	        endif
	        do k = 1,2
! loop on endpoints (k=1,2 at t=-1, +1, or Gauss point 0, Ng+1, resp.)
	          if     (k .eq. 1) then
	            g3 = 0
	          else
	            g3 = Ng+1
	          endif
	        enddo		! k
	      enddo		! j
	    enddo		! i
	    do m = 1,3		! direction loop
	      do p = 1,3	! direction loop
	        Qde1(ii,jj) = Qde1(ii,jj) - D(e,1,m) * dNdr(jj,1,g1,g2,g3) * Ji(m,p,g1,g2,g3)
	        Qde2(ii,jj) = Qde1(ii,jj) - D(e,2,m) * dNdr(jj,2,g1,g2,g3) * Ji(m,p,g1,g2,g3)
	        Qde3(ii,jj) = Qde1(ii,jj) - D(e,3,m) * dNdr(jj,3,g1,g2,g3) * Ji(m,p,g1,g2,g3)
	      enddo		! p
	    enddo		! m
	  enddo			! jj
	enddo			! ii

	return
	end

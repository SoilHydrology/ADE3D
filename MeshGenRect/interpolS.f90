	subroutine INTERPOLS(Ns, Nxy1, Nxy2, f, xy1, xy2, vl1, vl2, dcoor1, dcoor2, coor1, coor2, BCv, ierror)
! 2D interpolation at all nodes of an orthogonal quad face

! f is face
! s is species

	implicit none
	integer Ns, Nxy1(0:6), Nxy2(0:6), f, ierror
	real*8 xy1(6,Nxy1(0)), xy2(6,Nxy1(0)), vl1(6,Ns,Nxy1(0)), vl2(6,Ns,Nxy2(0))
	real*8 dcoor1, dcoor2, coor1(2), coor2(2), Bcv(6,Ns,4) 

	integer i, j, k, n, nnn
	real*8 tol

	tol = max(dcoor1, dcoor2) * 1.d-6	! geometric tolerance

! linear interpolation for the 2x2 nodes
	n = 0
	do i = 1,2
	  do j = 1,2
	    n = n+1

! switch n=3 and n=4
	    nnn = n
	    if(n .eq. 3)	nnn = 4
	    if(n .eq. 4)	nnn = 3

	    do k = 1,Nxy1(f)-1
	      if( coor1(i) .ge. xy1(f,k)-tol .and. coor1(i) .le. xy1(f,k+1)+tol .and.	&
	          coor2(j) .ge. xy2(f,k)-tol .and. coor2(j) .le. xy2(f,k+1)+tol ) then
! found - interpolate k interval
	        BCv(f,:,nnn) = ( vl1(f,:,k) + ( coor1(i) - xy1(f,k) ) *				&
			   ( vl1(f,:,k+1) - vl1(f,:,k) ) / ( xy1(f,k+1) - xy1(f,k) ) ) *	&
		               ( vl2(f,:,k) + ( coor2(j) - xy2(f,k) ) *				&
			   ( vl2(f,:,k+1) - vl2(f,:,k) ) / ( xy2(f,k+1) - xy2(f,k) ) )
	       endif
	    enddo	! k
	  enddo		! i
	enddo		! i

	return
	end

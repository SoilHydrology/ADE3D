	subroutine ASSEMBLE(iout, idbg, Ne, Nn, Nd, Npol, Ndof, NnNd, &
			vA, vL, vB, vQc1, vQd1, vQc2, vQd2, vQc3, vQd3, &
			rA, rL, rB, rQc1, rQd1, rQc2, rQd2, rQc3, rQd3, &
			cA, cL, cB, cQc1, cQd1, cQc2, cQd2, cQc3, cQd3, &
			lastA, lastL, lastB, lastQc1, lastQd1, lastQc2, lastQd2, lastQc3, lastQd3, &
			Ae, Le, Be, Qce1, Qde1, Qce2, Qde2, Qce3, Qde3, ie, nmat, e)
! assemble arrays

	implicit none
	integer iout, idbg
	integer Ne, Nn, Nd, NnNd		! array parameters
	integer Npol, Ndof			! array parameters
	integer lastA, lastL, lastB, lastQc1, lastQd1, lastQc2, lastQd2, lastQc3, lastQd3
	real*8 Ae(Npol+7,Npol+7), Le(Npol+7,Npol+7), Be(Npol+7,Npol+7)	! element arrays
	real*8 Qce1(Npol+7,Npol+7), Qde1(Npol+7,Npol+7)		! element arrays
	real*8 Qce2(Npol+7,Npol+7), Qde2(Npol+7,Npol+7)		! element arrays
	real*8 Qce3(Npol+7,Npol+7), Qde3(Npol+7,Npol+7)		! element arrays
	integer ie(Ne,9)			! global connectivity array
	integer nmat(Nn,0:Nd)			! global nodal materials array
	integer rA (Ndof+1), rL (Ndof+1), rB(Ndof+1)	! global  arrays (compact rows)
	integer cA (NnNd), cL (NnNd), cB(NnNd)! global  arrays (compact columns)
	integer rQc1(Ndof+1), rQd1(Ndof+1)		! global  arrays (compact rows)
	integer rQc2(Ndof+1), rQd2(Ndof+1)		! global  arrays (compact rows)
	integer rQc3(Ndof+1), rQd3(Ndof+1)		! global  arrays (compact rows)
	integer cQc1(NnNd), cQd1(NnNd)		! global  arrays (compact columns)
	integer cQc2(NnNd), cQd2(NnNd)		! global  arrays (compact columns)
	integer cQc3(NnNd), cQd3(NnNd)		! global  arrays (compact columns)
	real*8 vA (NnNd), vL (NnNd), vB(NnNd)	! global  arrays (compact values)
	real*8 vQc1(NnNd), vQd1(NnNd)		! global  arrays (compact values)
	real*8 vQc2(NnNd), vQd2(NnNd)		! global  arrays (compact values)
	real*8 vQc3(NnNd), vQd3(NnNd)		! global  arrays (compact values)
	integer e

	integer i, j, ii, jj

!	write(idbg,'(a)') ' --- ASSEMBLE ---'	! ### TEMPORARY ###

! assemble rank 2 sparse arrays
	do i = 1,Npol+7			! # of DOFs in the element
	  if(i .le. 8) then
	    ii = ie(e,i)		! global node #
	  else
	    ii = (i-9)*Ne + e + Nn	! global element DOF # for higher order i>2		
	  endif
	  do j = 1,Npol+7
	    if(j .le. 8) then
	      jj = ie(e,j)		! global node #
	    else
	      jj = (j-9)*Ne + e + Nn	! global element DOF # for higher order j>2
	    endif

! change roundoff entries to 0
	    if(abs( Ae(i,j) / Ae(1,1) ) .le. 1.d-10) Ae(i,j) = 0.d0
	    if(abs( Le(i,j) / Le(1,1) ) .le. 1.d-10) Le(i,j) = 0.d0
	    if(abs( Be(i,j) / Be(1,1) ) .le. 1.d-10) Be(i,j) = 0.d0

	    call ASSEMBLE2(iout, idbg, Ndof, NnNd, lastA, vA, cA, rA, &
			ii, jj, Ae(i,j), 1)

	    call ASSEMBLE2(iout, idbg, Ndof, NnNd, lastL, vL, cL, rL, &
			ii, jj, Le(i,j), 1)

	    call ASSEMBLE2(iout, idbg, Ndof, NnNd, lastB, vB, cB, rB, &
			ii, jj, Be(i,j), 1)

	  enddo		! j
	enddo		! i

! average rank 2 sparse arrays
!	do i = 1,8
	do i = 1,Npol+7			! # of DOFs in the element
	  if(i .le. 8) then
	    ii = ie(e,i)		! global node #
	  else
	    ii = (i-9)*Ne + e + Nn	! global element DOF # for higher order i>2		
	  endif
	  do j = 1,Npol+7
	    if(j .le. 8) then
	      jj = ie(e,j)		! global node #
	    else
	      jj = (j-9)*Ne + e + Nn	! global element DOF # for higher order j>2
	    endif

! change roundoff entries to 0
	    if(abs( Qce1(i,j) / Qce1(1,1) ) .le. 1.d-10) Qce1(i,j) = 0.d0
	    if(abs( Qde1(i,j) / Qde1(1,1) ) .le. 1.d-10) Qde1(i,j) = 0.d0
	    if(abs( Qce2(i,j) / Qce2(1,1) ) .le. 1.d-10) Qce2(i,j) = 0.d0
	    if(abs( Qde2(i,j) / Qde2(1,1) ) .le. 1.d-10) Qde2(i,j) = 0.d0
	    if(abs( Qce3(i,j) / Qce3(1,1) ) .le. 1.d-10) Qce3(i,j) = 0.d0
	    if(abs( Qde3(i,j) / Qde3(1,1) ) .le. 1.d-10) Qde3(i,j) = 0.d0

	    call ASSEMBLE2(iout, idbg, Ndof, NnNd, lastQc1, vQc1, cQc1, rQc1, &
		ii, jj, Qce1(i,j), nmat(ie(e,i),0))	! advection  comp. 1

	    call ASSEMBLE2(iout, idbg, Ndof, NnNd, lastQd1, vQd1, cQd1, rQd1, &
		ii, jj, Qde1(i,j), nmat(ie(e,i),0))	! dispersion comp. 1

	    call ASSEMBLE2(iout, idbg, Ndof, NnNd, lastQc2, vQc2, cQc2, rQc2, &
		ii, jj, Qce2(i,j), nmat(ie(e,i),0))	! advection  comp. 2

	    call ASSEMBLE2(iout, idbg, Ndof, NnNd, lastQd2, vQd2, cQd2, rQd2, &
		ii, jj, Qde2(i,j), nmat(ie(e,i),0))	! dispersion comp. 2

	    call ASSEMBLE2(iout, idbg, Ndof, NnNd, lastQc3, vQc3, cQc3, rQc3, &
		ii, jj, Qce3(i,j), nmat(ie(e,i),0))	! advection  comp. 3

	    call ASSEMBLE2(iout, idbg, Ndof, NnNd, lastQd3, vQd3, cQd3, rQd3, &
		ii, jj, Qde3(i,j), nmat(ie(e,i),0))	! dispersion comp. 3

	  enddo		! j
	enddo		! i
				
	return
	end


	subroutine SOURCE(iout, idbg, Nn, Ns, Ndof, NnNd, &
			vY, rY, cY, lastY, Kappa, T, phi)
! update reaction source
! NOTE: on the boundaries W = N by definition

	implicit none
	integer iout, idbg
	integer Nn, Ns, NnNd			! array parameters
	integer	Ndof				! array parameters
	integer lastY
	real*8 Kappa
	integer rY (Ndof+1)			! global  arrays (compact rows)
	integer cY (NnNd)			! global  arrays (compact columns)
	real*8 vY (NnNd )			! global  arrays (compact values)
	real*8 T  (Ndof,Ns)			! global  arrays
	real*8 phi(Ndof)			! global  arrays

	integer m, p				! local indices
	integer ii, jj				! global index
	real*8 aa, src

!	write(idbg,'(a)') ' --- SOURCE ---'	! ### TEMPORARY ###

! volumetric term
!----------------
	do ii = 1, Ndof
	  do jj = 1, Ndof
	    aa = 0.
! access rank 2 sparse arrays
	    call ACCESS2(iout, idbg, Ndof, NnNd, rY, cY, ii, jj, lastY, m)

	    if (m .ne. 0)	then
	      aa = vY(m)		! retrieve the Y (=A)
	    endif

	      src = -aa*Kappa*phi(jj)	! volumetric source term

! add volumetric term
	    T(ii,1) = T(ii,1) - src
	    T(ii,2) = T(ii,2) - src
	    T(ii,3) = T(ii,3) + src
	  enddo		! jj
	enddo		! ii

	return
	end

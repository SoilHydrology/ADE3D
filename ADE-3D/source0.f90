	subroutine SOURCE0(iout, idbg, Ne, Nn, Ns, Ndof, NnNd, ie, tc, Soe, &
			 vZ, rZ, cZ, lastZ, Son, Sn)
! calculate the independent source [A]{So} (time-independent)

	implicit none
	integer iout, idbg
	integer Ne, Nn, Ns, NnNd		! array parameters
	integer	Ndof				! array parameters
	integer lastZ
	real*8 tc				! global  arrays (compact values)
	integer ie(Ne,9)			! global connectivity array
	integer rZ(Ndof+1)			! global  arrays (compact rows)
	integer cZ(NnNd)			! global  arrays (compact columns)
	real*8 vZ(NnNd)				! global  arrays (compact values)
	real*8 Soe(Ne,Ns,8)			! element independent source nodal values
	real*8 Son(Ndof,Ns)			! global  arrays
	real*8 Sn (Ndof,Ns)			! global  arrays (nodal S at t^n+1)

	integer e, m, ii, jj, s
	real*8 aa, src

!	write(idbg,'(a)') ' --- SOURCE0 ---'	! ### TEMPORARY ###

! reset Sn and Son to 0
	Sn = 0.			! use matrix form
	Son = 0.		! use matrix form

! store the element So values in the nodal Sn array
	do e = 1,Ne
	  do m = 1,8
	    do s = 1,Ns
	      Sn(ie(e, m),s) = Soe(e,s,m)
	    enddo	! s
	  enddo		! m
	enddo		! e

	do jj = 1, Ndof
	  do ii = 1, Ndof
	    aa = 0.
! access rank 2 sparse arrays
	    call ACCESS2(iout, idbg, Ndof, NnNd, rZ, cZ, ii, jj, lastZ, m)

	    if (m .ne. 0)	then
	      aa = vZ(m)			! retrieve Z, the original A
	    endif

	    do s = 1,Ns
	      src = aa*Sn(jj,s)
	      Son(ii,s) = Son(ii,s) - src	! -[A]{S}
	    enddo	! s
	  enddo		! ii
	enddo		! jj

	return
	end

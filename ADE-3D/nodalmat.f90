	subroutine NODALMAT(iout, idbg, Ne, Nn, Nd, ie, nmat)
! Build a matrix of materials at each node

	implicit none
	integer iout, idbg
	integer Ne, Nn, Nd			! array parameters
	integer ie(Ne,9)			! global connectivity array
	integer nmat(Nn,0:Nd)			! global nodal materials array

	integer e, n, m, ierror, i, j

!	write(idbg,'(a)') ' --- NODALMAT ---'	! ### TEMPORARY ###

! initialize
	nmat = 0	! use matrix form
	ierror = 0

	do e = 1, Ne
	  m = ie(e,9)

	  do i = 1,8
! node i
	    n = ie(e,i)
	    j = nmat(n,0)	! nmat(n,0) is the node n counter
	    if(j .lt. Nd) then
	      nmat(n,0) = j+1	! increment counter
	      nmat(n,j+1) = m	! add material number
	    else
	      write(iout,*) '*** ABORT: nmat overflow at n=', n
	      ierror = ierror + 1
	    endif
	  enddo	! i

	enddo	! e

! error check
	if(ierror .ne. 0) then
	  write(iout,*) '*** ABORT: ierror = ', ierror
	  stop
	endif
	
	return
	end


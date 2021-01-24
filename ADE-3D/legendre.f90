	subroutine LEGENDRE(idbg, Ng, Npol, xg, Ain, Pn, dPn, pol)

! compute Legendre polynomial Pn(xi) and its derivative dPn(xi) at Gauss points and at the boundaries
!	  i=pol
! Pn(xi) = sum (Ain * xi^i)
!	  i=0

	implicit none
	integer idbg
	integer	Npol				! array parameters
	integer pol				! polynomial order (input)
	integer Ng				! vectors length (input)

	real*8 Ain(0:7, 0:7)			! Legendre polynomial Ain coefficients (input)
	real*8 xg (Ng)				! Legendre polynomial abcissas vector (input)
	real*8 Pn(0:Npol,0:Ng+1), dPn(0:Npol,0:Ng+1)	! Legendre polynomial values and derivatives vectors (output)

	integer i, k
	real*8 xi

!	write(idbg,'(a)') ' --- LEGENDRE ---'	! ### TEMPORARY ###

	do k = 0, Ng+1
	  if      (k .eq. 0   ) then
	    xi = -1.d0
	  else if (k .eq. Ng+1) then
	    xi =  1.d0
	  else
	    xi = xg(k)
	  endif

	  do i = 0, pol
	     Pn(pol,k) =  Pn(pol,k) + Ain(pol,i)     * xi** i		! value
	    if(i .ne. 0)					&
	    dPn(pol,k) = dPn(pol,k) + Ain(pol,i) * i * xi**(i-1)	! derivative
	  enddo	! i
	enddo	! k

	return
	end


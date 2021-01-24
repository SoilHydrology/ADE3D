c from SPARSKIT2\ITSOL\iters.f (see BCGSTAB.f)

c     (2) ALL iterative solvers require a user-supplied DOT-product
c     routine named DISTDOT. The prototype of DISTDOT is
c
c     real*8 function distdot(n,x,ix,y,iy)
c     integer n, ix, iy
c     real*8 x(1+(n-1)*ix), y(1+(n-1)*iy)
c
c     This interface of DISTDOT is exactly the same as that of
c     DDOT (or SDOT if real == real*8) from BLAS-1. It should have
c     same functionality as DDOT on a single processor machine. On a
c     parallel/distributed environment, each processor can perform
c     DDOT on the data it has, then perform a summation on all the
c     partial results.

c ### TEMPORARILY ### a call to DDOT is made 
c modified from SPARSKIT2\ITSOL\itaux.f

	function distdot(n,x,ix,y,iy)
	implicit none
	integer n, ix, iy
	real*8 distdot, x(1+(n-1)*ix), y(1+(n-1)*iy), ddot
	external ddot
	distdot = ddot(n,x,ix,y,iy)
	return
	end
c-----end-of-distdot


	subroutine SOLVE(iout, idbg, Nn, Nb, Ns, Ndof, NnNd, ldw, &
			Kappa, theta, theta1, ipar, fpar, ipar0, fpar0, &
			BCe, BCtype, C, T, Cold, Told, Rn, &
			vA, vL, vB, rA, rL, rB, cA, cL, cB, lastA, lastL, lastB, &
			phi, &
			w, work, &
			Nr, order, Rr)		! ### new parameters ###
! solve equations and update

	implicit none
	integer iout, idbg
	integer Nn, Nb, Ns, NnNd		! array parameters
	integer	Ndof				! array parameters
	integer lastA, lastL, lastB
	integer ldw
	real*8 Kappa
	integer ipar(16), ipar0(16)		! bCGstab integer parameters array
	real*8 fpar(16), fpar0(16)		! bCGstab real    parameters array
	real*8  theta, theta1
	integer rA (Ndof+1), rL (Ndof+1), rB(Ndof+1)	! global  arrays (compact rows)
	integer cA (NnNd), cL (NnNd), cB(NnNd)	! global  arrays (compact columns)
	integer BCe(Nb,5)			! BC element numbers
	character*1 BCtype(Nb)			! BC type ('R', 'N' or 'D')
	real*8 C   (Ndof,Ns), T   (Ndof,Ns)	! global  arrays
	real*8 Cold(Ndof,Ns), Told(Ndof,Ns) 	! global  arrays
	real*8 phi(Ndof)			! global  arrays
	real*8 Rn(Ndof,Ns), w(Ndof,Ns), work(ldw)	! work arrays
	real*8 vA (NnNd ), vL (NnNd ), vB(NnNd)	! global  arrays (compact values)

	integer Nr				! ### new parameters ###
	integer order(Ndof,0:1)			! ### new parameters ###
	real*8 Rr(Ndof,Ns)			! ### new parameters ###

	integer i, n, p, s

!	write(idbg,'(a)') ' --- SOLVE ---'	! ### TEMPORARY ###

! solve A dC/dt + B*I + T = 0

! A is alpha
! L is beta
! Cold and Told are C(n  ) and T(n ), respectively, where n is time step number
! C    and T    are C(n+1) and Tn+1), respectively, on exit

! calculate Rn+1 (the RHS at n+1) and store in Rn
! Rn+1 = beta*Cn - { theta*Tn+1 + (1-theta)*Tn }
!	 - B * sum on p#0 of { [ theta1 +  theta*Mp(dt)/Ap ] * Inp }
	Rn = theta*T + theta1*Told	! use matrix form (implicit)

! compute {Rn} = [L]{Cold} - {Rn} using AMUX from SPARSKIT2
	do s = 1, Ns
	  call AMUX(Ndof, Cold(1,s), w(1,s), vL, cL, rL)	! {w} = [L]{Cold}
	enddo		! s
	Rn = w - Rn

! modify Rn if Dirichlet BC
! order Rn
! permute Rn vector (in-place)
	do s = 1, Ns
	  call DVPERM(Ndof, Rn(:,s), order(:,1))
	enddo	! s
! subtract Rr
	Rn = Rn - Rr
! permute C vector (in-place)
	do s = 1, Ns
	  call DVPERM(Ndof, C (:,s), order(:,1))
	enddo	! s

! solve {[alpha]{Cn+1} = {Rn+1} using bCGstab
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do s = 1, Ns
	  ipar = ipar0	! restore the original ipar
	  fpar = fpar0	! restore the original fpar

10	  call bCGstab(Nr, Rn(1,s), C(1,s), ipar, fpar , work)
	  if (ipar(1) .eq. 1) then
	    call AMUX(Nr, work(ipar(8)), work(ipar(9)), vA, cA, rA)
	    go to 10
	  else if(ipar(1) .ne. 0)	then
	    write(idbg,'(a)') ' --- SOLVE ---'	! ### TEMPORARY ###
	    write(idbg,*) 'ipar = ', ipar
	    write(idbg,*) 'fpar = ', fpar
	    stop 'problem in bCGstab'
	  endif
	enddo	! s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! permute C vector (in-place) back to the original order
	do s = 1, Ns
	  call DVPERM(Ndof, C (:,s), order(:,0))
	enddo	! s

! force non-negative C
	  do i = 1, Nn
	    do s = 1, Ns
	      if(C(i,s) .lt. 0.)	C(i,s) = 0.
	    enddo		! s
	  enddo			! i

	  if(Kappa .ne. 0) then
	    do i = 1, Ndof
	      phi(i) = C(i,1)*C(i,2)	! store Ca*Cb in phi
	    enddo		! i
	  endif

	return
	end

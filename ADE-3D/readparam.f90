	subroutine READPARAM(iout, idbg, Ne, Nn, Nb, Nm, Nd, Ng, Ns, Npol, Nf)
! Read parameters

	implicit none
	integer iout, idbg
	integer Ne, Nn, Nb, Nm, Nd, Ng, Ns, Npol, Nf	! array parameters

!	write(idbg,'(a)') ' --- READPARAM ---'	! ### TEMPORARY ###

! read parameters from inp.txt
	open(1, file='finp.txt', status='old')

! parameters:
	read (1,*) Ne, Nn, Nb, Nm
	read (1,*) Nd, Ns, Npol, Nf

	Ng = Npol + 1	! Ng >= Npol + 1/2 for exact quadrature

	write(iout,*) 'Ne = ', Ne
	write(iout,*) 'Nn = ', Nn
	write(iout,*) 'Nb = ', Nb
	write(iout,*) 'Nm = ', Nm
	write(iout,*) 'Nd = ', Nd
	write(iout,*) 'Npol = ', Npol
	write(iout,*) 'Ns = ', Ns
	write(iout,*) 'Nf = ', Nf
! error check
	if(Ne .le. 0 .or. Nn .le. 0 .or. Nb .le. 0 .or. Nm .le. 0 .or. &
	   Nd .le. 0 .or. Ns .le. 0 .or. Ns .gt. 3 .or. Nf .le. 0 .or. &
           Npol .lt. 1 .or. Npol .gt. 7 ) then
	  write(iout,*) '*** ABORT: illegal Ne, Nn, Nb, Nm, Nd, Ng, Ns, Nf or Npol'
! NOTE that currently Ns should not exceed 3 for the reaction A+B->C
! NOTE that Npol, the polynomial order, should be 1<= Npol <= M, where currently M=7 is assumed
	  stop
	endif
	
	return
	end


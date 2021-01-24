	subroutine ASSEMBLE2(iout, idbg, Ndof, NnNd, last, vY, cY, rY, r, c, v, m)
! assemble (or average) rank 2 sparse arrays

	implicit none
	integer iout, idbg
	integer Ndof, NnNd
	integer last, r, c
	integer m				! elements at node
	real*8 v
	real*8  vY(NnNd)
	integer cY(NnNd), rY(Ndof+1)

	integer iadd, i, j, ibeg
	
!	write(idbg,'(a)') ' --- ASSEMBLE2 ---'	! ### TEMPORARY ###

! existing v(r,c) (UNSORTED matrix assumed)
	call ACCESS2(iout, idbg, Ndof, NnNd, rY, cY, r, c, last, iadd)

	if(iadd .ne. 0) then
! found nonzero. assemble (or average) and exit
	  vY(iadd) = v/m + vY(iadd)
	else
! found zero
! check for overflow
	  if(last+1 .gt. NnNd) then
	    write(iout,*) 'from ASSEMBLE2: last, NnNd = ', last, NnNd
	    write(iout,'(a)') 'ABORT: last+1 .gt. NnNd'
	    stop
	  endif

! move vY and cY from the beginning of row r to the last 1 position
	  ibeg = rY(r  )	! first nonzero column in row r
	  do j = last, ibeg, -1
	    vY(j+1) = vY(j)	! move vY's
	    cY(j+1) = cY(j)	! move cY's
	  enddo	! j
! insert c and v (or its average) at the beginning of row r
	  vY(ibeg) = v/m	! insert v (or its average) to vY(ibeg) 
	  cY(ibeg) = c		! insert c to cY(ibeg) 
! increment rY for rows r+1 onwards
	  do i = Ndof+1, r+1, -1
	    rY(i) = rY(i)+1	! increment rY's
	  enddo	! i

	  last = last+1	! increment last

	endif
	return
	end


	subroutine ACCESS2(iout, idbg, Ndof, NnNd, rY, cY, r, c, last, iadd)
! access rank 2 sparse arrays
! modified from SPARSKIT getelm.f

	implicit none
	integer iout, idbg
	integer r, c, last			! current row, column and index
	integer Ndof, NnNd			! array parameters
	integer rY(Ndof+1), cY(NnNd)		! row and column arrays
	integer iadd				! returned index

	integer ibeg, iend, k

!	write(idbg,'(a)') ' --- ACCESS2 ---'	! ### TEMPORARY ###

! search index iadd for r,c

! initialization 

	ibeg = rY(r)
	iend = rY(r+1)-1

! scan the row - exit as soon as Y(r,c) is found

	do k = ibeg, iend
	  if (cY(k) .eq.  c) then
! found. return iadd=k
	    iadd = k
	    return
	  endif
	enddo	! k
   
! not found. return iadd=0
	iadd = 0 
	return
	end 


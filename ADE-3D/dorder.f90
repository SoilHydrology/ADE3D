	subroutine DORDER(iout, idbg, Ne, Nn, Nb, Ndof, BCe, BCi, BCtype, ie, &
			N_D, order)	! ### new parameters ###
! create a list of the original vs reordered nodes

! ### new parameters ###
! N_D is the total number of D type nodes
!	(output)
! D_node(Ndof) is an integer array. on return it holds D type nodes
! order(Ndof,0:1) is an integer array. on return
!	its 0 column holds the original D type number and 
!	its 1 column holds the reodered D type number
!	to be allocated in MAIN and passed
!	(input and output)

	implicit none
	integer iout, idbg
	integer Ne, Nn, Nb			! array parameters
	integer	Ndof				! array parameters
	integer N_D				! ### new parameter ###
	integer BCe(Nb,5), BCi(Nb)		! BC element and local element face numbers
	integer ie(Ne,9)			! global connectivity array
	integer order(Ndof,0:1)			! ### new parameter ###
	character*1 BCtype(Nb)			! BC type ('R', 'N' or 'D')

	integer, allocatable :: D_node(:)	! ### new parameter ###
	integer i, e, i1, i2, i3, i4, n

!	write(idbg,'(a)') ' --- DORDER ---'	! ### TEMPORARY ###
	allocate (D_node(Ndof))

! reset D_node to 0
	D_node = 0		! use matrix form
	N_D = 0

	do n = 1, Nb
! find nodes with Dirichlet BC
	  if(BCtype(n) .eq. 'D') then
	    e = BCe(n,5)	! BC element number
	    i = BCi(n)		! BC local element face number
	    if     (i .eq. 1) then
! face 1 - west
	      i1 = ie(e,1)
	      i2 = ie(e,4)
	      i3 = ie(e,8)
	      i4 = ie(e,5)
	    else if(i .eq. 2) then
! face 2 - south
	      i1 = ie(e,1)
	      i2 = ie(e,5)
	      i3 = ie(e,6)
	      i4 = ie(e,2)
	    else if(i .eq. 3) then
! face 3 - bottom
	      i1 = ie(e,1)
	      i2 = ie(e,2)
	      i3 = ie(e,3)
	      i4 = ie(e,4)
	    else if(i .eq. 4) then
! face 4 - east
	      i1 = ie(e,2)
	      i2 = ie(e,3)
	      i3 = ie(e,7)
	      i4 = ie(e,6)
	    else if(i .eq. 5) then
! face 5 - north
	      i1 = ie(e,4)
	      i2 = ie(e,8)
	      i3 = ie(e,7)
	      i4 = ie(e,3)
	    else if(i .eq. 6) then
! face 6 - top
	      i1 = ie(e,5)
	      i2 = ie(e,6)
	      i3 = ie(e,7)
	      i4 = ie(e,8)
	    endif

! store the Dirichlet BC node numbers and count them
! (if they were not yet stored or counted)
	    if(D_node(i1).eq. 0) then
	      N_D = N_D +1
	      D_node(i1) = i1
	    endif
	    if(D_node(i2).eq. 0) then
	      N_D = N_D +1
	      D_node(i2) = i2
	    endif
	    if(D_node(i3).eq. 0) then
	      N_D = N_D +1
	      D_node(i3) = i3
	    endif
	    if(D_node(i4).eq. 0) then
	      N_D = N_D +1
	      D_node(i4) = i4
	    endif
	  endif
	enddo		! n

	i1 = 0		! reset the initial non-D type location
	i2 = Ndof-N_D	! reset the initial     D type location
	do n = 1, Ndof
	  if(D_node(n) .eq. 0) then
! order the non-D type nodes (reorder to the first Ndof-N_D locations)
	    i1 = i1 + 1
	    order(n,1) = i1		! reordered
	  else
! order the     D type nodes (reorder to the last     N_D locations)
	    i2 = i2 + 1
	    order(D_node(n),1) = i2	! reordered
	  endif
	enddo		! n

	do n = 1, Ndof
	  order(order(n,1),0) = n	! original
	enddo		! n

	deallocate (D_node)

	return
	end


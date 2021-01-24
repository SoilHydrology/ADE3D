	subroutine OUT(iout, idbg, ipost, Nn, Ns, Ndof, NnNd, ldw, &
			time, Nc, C, T, w, work, &
			lastQc1, lastQd1, lastF1, lastQc2, lastQd2, lastF2, lastQc3, lastQd3, lastF3, &
			vQc1, vQd1, vF1, rQc1, rQd1, rF1, cQc1, cQd1, cF1, &	
			vQc2, vQd2, vF2, rQc2, rQd2, rF2, cQc2, cQd2, cF2, &	
			vQc3, vQd3, vF3, rQc3, rQd3, rF3, cQc3, cQd3, cF3)	
! write solver output, C, T and nodal fluxes

	implicit none
	integer iout, idbg, ipost
	integer Nn, Ns, ldw, NnNd			! array parameters
	integer	 Ndof					! array parameters
	integer lastQc1, lastQd1, lastF1, lastQc2, lastQd2, lastF2, lastQc3, lastQd3, lastF3
	integer rQc1(Ndof+1), rQd1(Ndof+1), rF1(Ndof+1)	! global  arrays (compact rows)
	integer rQc2(Ndof+1), rQd2(Ndof+1), rF2(Ndof+1)	! global  arrays (compact rows)
	integer rQc3(Ndof+1), rQd3(Ndof+1), rF3(Ndof+1)	! global  arrays (compact rows)
	integer cQc1(NnNd), cQd1(NnNd), cF1(NnNd)	! global  arrays (compact columns)
	integer cQc2(NnNd), cQd2(NnNd), cF2(NnNd)	! global  arrays (compact columns)
	integer cQc3(NnNd), cQd3(NnNd), cF3(NnNd)	! global  arrays (compact columns)
	real*8 time
	real*8 C   (Ndof,Ns), T   (Ndof,Ns)		! global  arrays
	real*8 vQc1(NnNd), vQd1(NnNd), vF1(NnNd)	! global  arrays (compact values)
	real*8 vQc2(NnNd), vQd2(NnNd), vF2(NnNd)	! global  arrays (compact values)
	real*8 vQc3(NnNd), vQd3(NnNd), vF3(NnNd)	! global  arrays (compact values)
	real*8 w(Ndof,Ns), work(ldw)			! work arrays
	integer Nc

	integer i, p, ii, jj, s

	write(idbg,'(a)') ' --- OUT ---'	! ### TEMPORARY ###

	Nc = Nc + 1		! counter for Gmsh

! output C
	do s = 1, Ns
	  write(ipost,'(a)') '$NodeData'
	  write(ipost,'(i10, a)') 1, '	// number-of-string-tags'
	  write(ipost,'(a2,i1,a)') '"C', s , '"	// < "string-tag" > ...'
	  write(ipost,'(i10, a)') 1, '	// number-of-real-tags'
	  write(ipost,'(g14.7, a)') time, '	// < real-tag > ... time'
	  write(ipost,'(i10, a)') 3, '	// number-of-integer-tags'
	  write(ipost,'(i10, a)') Nc, '	// step #'
	  write(ipost,'(i10, a)') 1, '	// number of components'
	  write(ipost,'(i10, a)') Nn, '	// number of nodes'
	do i = 1, Nn
	    write(ipost,'(i10, 3g14.7)') i,  C(i,s)
	enddo	!i
	  write(ipost,'(a)') '// node-number value ...'
	  write(ipost,'(a)') '$EndNodeData'
	  write(ipost,*)
	enddo	! s
! output T
	do s = 1, Ns
	  write(ipost,'(a)') '$NodeData'
	  write(ipost,'(i10, a)') 1, '	// number-of-string-tags'
	  write(ipost,'(a2,i1,a)') '"T', s , '"	// < "string-tag" > ...'
	  write(ipost,'(i10, a)') 1, '	// number-of-real-tags'
	  write(ipost,'(g14.7, a)') time, '	// < real-tag > ... time'
	  write(ipost,'(i10, a)') 3, '	// number-of-integer-tags'
	  write(ipost,'(i10, a)') Nc, '	// step #'
	  write(ipost,'(i10, a)') 1, '	// number of components'
	  write(ipost,'(i10, a)') Nn, '	// number of nodes'
	  do i = 1, Nn
	    write(ipost,'(i10, 3g14.7)') i,  T(i,s)
	  enddo	! i
	  write(ipost,'(a)') '// node-number value ...'
	  write(ipost,'(a)') '$EndNodeData'
	  write(ipost,*)
	enddo	! s
! nodal advection flux, qc
! compute and write {qc} = [Qc]{C}
	do s = 1, Ns
	  write(ipost,'(a)') '$NodeData'
	  write(ipost,'(i10, a)') 1, '	// number-of-string-tags'
	  write(ipost,'(a3,i1,a)') '"qc', s , '" // < "string-tag" > ...'
	  write(ipost,'(i10, a)') 1, '	// number-of-real-tags'
	  write(ipost,'(g14.7, a)') time, '	// < real-tag > ... time'
	  write(ipost,'(i10, a)') 3, '	// number-of-integer-tags'
	  write(ipost,'(i10, a)') Nc, '	// step #'
	  write(ipost,'(i10, a)') 3, '	// number of components'
	  write(ipost,'(i10, a)') Nn, '	// number of nodes'
! compute {work(       1:  Ndof)} = [Qc1]{C} using AMUX from SPARSKIT2
! compute {work(  Ndof+1:2*Ndof)} = [Qc2]{C} using AMUX from SPARSKIT2
! compute {work(2*Ndof+1:3*Ndof)} = [Qc3]{C} using AMUX from SPARSKIT2
	  call AMUX(Ndof, C(1,s), work(       1), vQc1, cQc1, rQc1)	! {work(       1:  Ndof)} = [Qc1]{C}
	  call AMUX(Ndof, C(1,s), work(  Ndof+1), vQc2, cQc2, rQc2)	! {work(  Ndof+1:2*Ndof)} = [Qc2]{C}
	  call AMUX(Ndof, C(1,s), work(2*Ndof+1), vQc3, cQc3, rQc3)	! {work(2*Ndof+1:3*Ndof)} = [Qc3]{C}
	  do ii = 1, Nn
	    write(ipost,'(i10, 3g14.7)') ii, work(ii), work(Ndof+ii), work(2*Ndof+ii)
	  enddo	! ii
	  write(ipost,'(a)') '// node-number value ...'
	  write(ipost,'(a)') '$EndNodeData'
	  write(ipost,*)
	enddo	! s
! nodal dispersion flux, qd
! compute and write {qd} = [Qd]{C}
	do s = 1, Ns
	  write(ipost,'(a)') '$NodeData'
	  write(ipost,'(i10, a)') 1, '	// number-of-string-tags'
	  write(ipost,'(a3,i1,a)') '"qd', s , '" // < "string-tag" > ...'
	  write(ipost,'(i10, a)') 1, '	// number-of-real-tags'
	  write(ipost,'(g14.7, a)') time, '	// < real-tag > ... time'
	  write(ipost,'(i10, a)') 3, '	// number-of-integer-tags'
	  write(ipost,'(i10, a)') Nc, '	// step #'
	  write(ipost,'(i10, a)') 3, '	// number of components'
	  write(ipost,'(i10, a)') Nn, '	// number of nodes'
! compute {work(       1:  Ndof)} = [Qd1]{C} using AMUX from SPARSKIT2
! compute {work(  Ndof+1:2*Ndof)} = [Qd2]{C} using AMUX from SPARSKIT2
! compute {work(2*Ndof+1:3*Ndof)} = [Qd3]{C} using AMUX from SPARSKIT2
	  call AMUX(Ndof, C(1,s), work(       1), vQd1, cQd1, rQd1)	! {work(       1:  Ndof)} = [Qd1]{C}
	  call AMUX(Ndof, C(1,s), work(  Ndof+1), vQd2, cQd2, rQd2)	! {work(  Ndof+1:2*Ndof)} = [Qd2]{C}
	  call AMUX(Ndof, C(1,s), work(2*Ndof+1), vQd3, cQd3, rQd3)	! {work(2*Ndof+1:3*Ndof)} = [Qd3]{C}
	  do ii = 1, Nn
	    write(ipost,'(i10, 3g14.7)') ii, work(ii), work(Ndof+ii), work(2*Ndof+ii)
	  enddo	! ii
	  write(ipost,'(a)') '// node-number value ...'
	  write(ipost,'(a)') '$EndNodeData'
	  write(ipost,*)
	enddo	! s
! nodal flux, j
	do s = 1, Ns
	    w = C					! use matrix form
! compute and write {j} = [F]{w}
	  write(ipost,'(a)') '$NodeData'
	  write(ipost,'(i10, a)') 1, '	// number-of-string-tags'
	  write(ipost,'(a2,i1,a)') '"j', s , '" // < "string-tag" > ...'
	  write(ipost,'(i10, a)') 1, '	// number-of-real-tags'
	  write(ipost,'(g14.7, a)') time, '	// < real-tag > ... time'
	  write(ipost,'(i10, a)') 3, '	// number-of-integer-tags'
	  write(ipost,'(i10, a)') Nc, '	// step #'
	  write(ipost,'(i10, a)') 3, '	// number of components'
	  write(ipost,'(i10, a)') Nn, '	// number of nodes'
! compute {work(       1:  Ndof)} = [F1]{w} using AMUX from SPARSKIT2
! compute {work(  Ndof+1:2*Ndof)} = [F2]{w} using AMUX from SPARSKIT2
! compute {work(2*Ndof+1:3*Ndof)} = [F3]{w} using AMUX from SPARSKIT2
	  call AMUX(Ndof, w(1,s), work(       1), vF1, cF1, rF1)	! {work(       1:  Ndof)} = [F1]{w}
	  call AMUX(Ndof, w(1,s), work(  Ndof+1), vF2, cF2, rF2)	! {work(  Ndof+1:2*Ndof)} = [F2]{w}
	  call AMUX(Ndof, w(1,s), work(2*Ndof+1), vF3, cF3, rF3)	! {work(2*Ndof+1:3*Ndof)} = [F3]{w}
	  do ii = 1, Nn
	    write(ipost,'(i10, 3g14.7)') ii, work(ii), work(Ndof+ii), work(2*Ndof+ii)
	  enddo		! ii
	  write(ipost,'(a)') '// node-number value ...'
	  write(ipost,'(a)') '$EndNodeData'
	   write(ipost,*)
	enddo		! s

	return
	end


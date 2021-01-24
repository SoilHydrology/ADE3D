	subroutine SHAPEB(iout, idbg, Ne, Ng, Npol, Nb, BCe, Bci, ie, &
			  Shp, Shpb)
! calculate element shape functions on boundary faces
! note that on element boundaries W=N

	implicit none
	integer iout, idbg
	integer Ne, Nb, Ng			! array parameters
	integer	Npol				! array parameters
	integer BCe(Nb,5), BCi(Nb)		! BC element and local element face numbers
	integer ie(Ne,9)			! global connectivity array
	real*8 Shp(Npol+7,0:Ng+1,0:Ng+1,0:Ng+1)	! shape and weight functions
	real*8 Shpb(Nb,4,Ng,Ng)			! boundary shape functions

	integer e, Lf, n, g1, g2, g3

!	write(idbg,'(a)') ' --- SHAPEB ---'	! ### TEMPORARY ###

! linear 3D shape functions on a boundary face
! Ni(r,s,t) = (1 +/- r)(1 +/- s)(1 +/- t) / 8 ; -1 < r, s, t < +1

!================================================================================
! from Rami\FEM-3D\01-ADE:
! o face order in both MeshGEn-Rect-3D (v01) FEM-3D (v01)	of 16/09/19
! o correction: n4 = 8->6 in east				of 23/09/19
! ================================== face order ==================================
! face	name	fixed	coor1	coor2	n1	n2	n3	n4	normal
! 1	west	x=Xmin	y	z	1	4	8	5	inwards
! 2	south	y=Ymin	z	x	1	5	6	2	inwards
! 3	bottom	z=Zmin	x	y	1	2	3	4	inwards
! 4	east	x=Xmax	y	z	2	3	7	6	outwards
! 5	north	y=Ymax	z	x	4	8	7	3	outwards
! 6	top	z=Zmax	x	y	5	6	7	8	outwards
!================================== face order ==================================
!================================================================================

	do n = 1, Nb
	  e  = BCe(n,5)		! BC element number
	  Lf = BCi(n)		! BC local face number

	  if     (Lf .eq. 1) then		! face 1, r=-1
	    g1 = 0
	    do g2 = 1,Ng
	      do g3 = 1,Ng
	        Shpb(n,1,g2,g3) = Shp(1,g1,g2,g3)
	        Shpb(n,2,g2,g3) = Shp(4,g1,g2,g3)
	        Shpb(n,3,g2,g3) = Shp(8,g1,g2,g3)
	        Shpb(n,4,g2,g3) = Shp(5,g1,g2,g3)
	      enddo	! g3
	    enddo	! g2
	  else if(Lf .eq. 4) then		! face 4, r=+1
	    g1 = Ng+1
	    do g2 = 1,Ng
	      do g3 = 1,Ng
	        Shpb(n,1,g2,g3) = Shp(2,g1,g2,g3)
	        Shpb(n,2,g2,g3) = Shp(3,g1,g2,g3)
	        Shpb(n,3,g2,g3) = Shp(7,g1,g2,g3)
	        Shpb(n,4,g2,g3) = Shp(6,g1,g2,g3)
	      enddo	! g3
	    enddo	! g2
	  else if(Lf .eq. 2) then		! face 2, s=-1
	    g2 = 0
	    do g1 = 1,Ng
	      do g3 = 1,Ng
	        Shpb(n,1,g1,g3) = Shp(1,g1,g2,g3)
	        Shpb(n,2,g1,g3) = Shp(5,g1,g2,g3)
	        Shpb(n,3,g1,g3) = Shp(6,g1,g2,g3)
	        Shpb(n,4,g1,g3) = Shp(2,g1,g2,g3)
	      enddo	! g3
	    enddo	! g1
	  else if(Lf .eq. 5) then		! face 5, s=+1
	    g2 = Ng+1
	    do g1 = 1,Ng
	      do g3 = 1,Ng
	        Shpb(n,1,g1,g3) = Shp(4,g1,g2,g3)
	        Shpb(n,2,g1,g3) = Shp(8,g1,g2,g3)
	        Shpb(n,3,g1,g3) = Shp(7,g1,g2,g3)
	        Shpb(n,4,g1,g3) = Shp(3,g1,g2,g3)
	      enddo	! g3
	    enddo	! g1
	  else if(Lf .eq. 3) then		! face 3, t=-1
	    g3 = 0
	    do g1 = 1,Ng
	      do g2 = 1,Ng
	        Shpb(n,1,g1,g2) = Shp(1,g1,g2,g3)
	        Shpb(n,2,g1,g2) = Shp(2,g1,g2,g3)
	        Shpb(n,3,g1,g2) = Shp(3,g1,g2,g3)
	        Shpb(n,4,g1,g2) = Shp(4,g1,g2,g3)
	      enddo	! g2
	    enddo	! g1
	  else if(Lf .eq. 6) then		! face 6, t=+1
	    g3 = Ng+1
	    do g1 = 1,Ng
	      do g2 = 1,Ng
	        Shpb(n,1,g1,g2) = Shp(5,g1,g2,g3)
	        Shpb(n,2,g1,g2) = Shp(6,g1,g2,g3)
	        Shpb(n,3,g1,g2) = Shp(7,g1,g2,g3)
	        Shpb(n,4,g1,g2) = Shp(8,g1,g2,g3)
	      enddo	! g2
	    enddo	! g1
	  endif

	enddo	! n
	return
	end

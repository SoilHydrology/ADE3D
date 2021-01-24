	program MeshGenRect_3D
! Rectangular box mesh generator (3D)
! assumptions: (1) material #=1, (2) cons Vi, Dij

	implicit none
	integer Nx, Ny, Nz, Ne, Nn, Nb, Nf, BCe, ie(9), Ns, Nxy1(0:6), Nxy2(0:6), ierror
	integer Nxx, Nyy, Nzz
	integer e, n, f, ix, iy, iz, k, s, Lf, Gf
	real*8 Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Lx, Ly, Lz, dx, dy, dz, Vx, Vy, Vz, Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, x(3)
	real*8 coor1(2), coor2(2), dcoor1, dcoor2
	logical Darcy					! T Darcy-2D, F for FEM-2D
	logical Dedge
	logical delta					! T for delta(Xdelta, Ydelta, Zdelta)
	real*8 Xdelta, Ydelta, Zdelta
	real*8 Vf, Vn
	character*1 BCt(6)				! south, east, north, west, bottom, top

	real*8, allocatable :: C(:), BCv(:,:,:), xy1(:,:), vl1(:,:,:), xy2(:,:), vl2(:,:,:)
	real*8, allocatable :: xx(:), yy(:), zz(:), fx(:,:), gy(:,:), hz(:,:), Sxy(:,:)

	data n/0/, Gf/0/, ierror/0/

! write param file (part of finp.txt)
	open(3, file='param.txt', status='unknown')
	write(3,*) '**********************************'	! version ID
	write(3,*) 'MeshGenRect_3D v01-f-1 of 02/12/20'	! version ID
	write(3,*) '**********************************'	! version ID

! read MeshGenRect input file
	open(1, file='meshinp.txt', status='old')
	read(1,*) Darcy	
	read(1,*) Xmin, Xmax, Ymin, Ymax, Zmin, Zmax
	if(Darcy) then
	  read(1,*) Nx, Ny, Nz
	  Ns = 1	! Ns=1 used for Darcy
	else
	  read(1,*) Nx, Ny, Nz, Ns
	endif

	allocate ( C(Ns), BCv(6,Ns,4) )

	if(.not. Darcy) then
	  read(1,*) Vx, Vy, Vz
	endif
	read(1,*) Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz

! an independent source is assumed to be given as S=f(x)*g(y)*h(z) for each species
	read(1,*) Nxx, Nyy, Nzz			! No. of x, y,z data points (must be >1)

	if(Nxx.lt.2 .or. Nyy.lt.2 .or. Nzz.lt.2) then
	  write(3,*) '!!! ERROR !!! Nxx or Nyy or Nzz <2 !'
	  stop
	endif

	allocate ( xx(Nxx), fx(Ns,Nxx), yy(Nyy), gy(Ns,Nyy), zz(Nzz), hz(Ns,Nzz), Sxy(8,Ns) )

! f(x)
	read(1,*) (xx(k), k=1,Nxx)	! coordinate x
	do s = 1,Ns
	  read(1,*) (fx(s,k), k=1,Nxx)	! f(x) value
	enddo	! s
! g(y)
	read(1,*) (yy(k), k=1,Nyy)	! coordinate y
	do s = 1,Ns
	  read(1,*) (gy(s,k), k=1,Nyy)	! g(y) value
	enddo	! s
! h(z)
	read(1,*) (zz(k), k=1,Nzz)	! coordinate z
	do s = 1,Ns
	  read(1,*) (hz(s,k), k=1,Nzz)	! h(z) value
	enddo	! s

	if(xx(1).ne.Xmin .or. xx(Nxx).ne.Xmax .or. &
           yy(1).ne.Ymin .or. yy(Nyy).ne.Ymax .or. &
           zz(1).ne.Zmin .or. zz(Nzz).ne.Zmax) then
	  write(3,*) '!!! ERROR !!! xx(1).ne.Xmin .or. xx(Nxx).ne.Xmax .or.', &
	 	'yy(1).ne.Ymin .or. yy(Nyy).ne.Ymax .or.', &
	 	'zz(1).ne.Zmin .or. zz(Nzz).ne.Zmax !'
	  stop
	endif

	read(1,*) delta, Xdelta, Ydelta, Zdelta
	read(1,*) (C(s), s=1,Ns)
	read(1,*) (BCt(f), f=1,6)
	read(1,*) (Nxy1(f), f=1,6)			! 1st coor. No. of data points (must be >1)
	read(1,*) (Nxy2(f), f=1,6)			! 2nd coor. No. of data points (must be >1)

	if(Nxy1(1).lt.2 .or.Nxy1(2).lt.2 .or.Nxy1(3).lt.2 .or.Nxy1(4).lt.2 .or.Nxy1(5).lt.2 .or.Nxy1(6).lt.2 .or. &
	   Nxy2(1).lt.2 .or.Nxy2(2).lt.2 .or.Nxy2(3).lt.2 .or.Nxy2(4).lt.2 .or.Nxy2(5).lt.2 .or.Nxy2(6).lt.2) then
	  write(3,*) '!!! ERROR !!! Nxy*(f) <2 !'
	  stop
	endif

	Nxy1(0) = max(Nxy1(1),Nxy1(2),Nxy1(3),Nxy1(4),Nxy1(5),Nxy1(6))	! 1st coor. max. nodes in each face
	Nxy2(0) = max(Nxy2(1),Nxy2(2),Nxy2(3),Nxy2(4),Nxy2(5),Nxy2(6))	! 2nd coor. max. nodes in each face
	allocate ( xy1(6,Nxy1(0)), vl1(6,Ns,Nxy1(0)) )
	allocate ( xy2(6,Nxy2(0)), vl2(6,Ns,Nxy2(0)) )

	do f=1,6
	  read(1,*) (xy1(f,k), k=1,Nxy1(f))	! 1st coor. coordinate (x, y or z)
	  do s = 1,Ns
	    read(1,*) (vl1(f,s,k), k=1,Nxy1(f))	! 1st coord. value
	  enddo	! s
	  read(1,*) (xy2(f,k), k=1,Nxy2(f))	! 2nd coor. coordinate (x, y or z)
	  do s = 1,Ns
	    read(1,*) (vl2(f,s,k), k=1,Nxy2(f))	! 2nd coord. value
	  enddo	! s
	enddo	! f
	close(1)

	if(xy1(1,1).ne.Ymin .or. xy1(1,Nxy1(1)).ne.Ymax .or. &  ! west
	   xy1(2,1).ne.Zmin .or. xy1(2,Nxy1(2)).ne.Zmax .or. &  ! south
	   xy1(3,1).ne.Xmin .or. xy1(3,Nxy1(3)).ne.Xmax .or. &  ! bottom
	   xy1(4,1).ne.Ymin .or. xy1(4,Nxy1(4)).ne.Ymax .or. &  ! east 
	   xy1(5,1).ne.Zmin .or. xy1(5,Nxy1(5)).ne.Zmax .or. &  ! north
	   xy1(6,1).ne.Xmin .or. xy1(6,Nxy1(6)).ne.Xmax .or. &  ! top

	   xy2(1,1).ne.Zmin .or. xy2(1,Nxy2(1)).ne.Zmax .or. &  ! west
	   xy2(2,1).ne.Xmin .or. xy2(2,Nxy2(2)).ne.Xmax .or. &  ! south
	   xy2(3,1).ne.Ymin .or. xy2(3,Nxy2(3)).ne.Ymax .or. &  ! bottom
	   xy2(4,1).ne.Zmin .or. xy2(4,Nxy2(4)).ne.Zmax .or. &  ! east
	   xy2(5,1).ne.Xmin .or. xy2(5,Nxy2(5)).ne.Xmax .or. &  ! north
	   xy2(6,1).ne.Ymin .or. xy2(6,Nxy2(6)).ne.Ymax) then   ! top
	  write(3,*) '!!! ERROR !!! xy*(f,1).ne.Xmin or Ymin Zmin .or. xx*(f,Nxy*(f)).ne.Xmax or Ymax Zmax!'
	  stop
	endif

	Lx = Xmax-Xmin	
	Ly = Ymax-Ymin
	Lz = Zmax-Zmin	

	Ne = Nx*Ny*Nz					! # of elements
	Nn = (Nx+1)*(Ny+1)*(Nz+1)			! # of nodes
	Nf = (Nx+1)*Ny*Nz + Nx*(Ny+1)*Nz + Nx*Ny*(Nz+1)	! # of element faces
	Nb = 2*(Nx*Nz + Ny*Nz + Nx*Ny)			! # of BC faces
	dx = Lx/Nx
	dy = Ly/Ny
	dz = Lz/Nz

	write(3,'(6i10)')  Ne,  Nn,  Nb,   1, Ns, Nf
	write(3,'(a   )') 'Ne,  Nn,  Nb,  Nm, Ns, Nf'

! write nodes file
	open(1, file='nodes.txt', status='unknown')

	do iz = 1, Nz+1
	  do iy = 1, Ny+1
	    do ix = 1, Nx+1
	      n = n+1
	      x(1) = dx * (ix-1) + Xmin
	      x(2) = dy * (iy-1) + Ymin
	      x(3) = dz * (iz-1) + Zmin
! Co = delta(x,y)
	      if(delta)	then
	        if(abs(x(1)-Xdelta).lt.0.6*dx .and. abs(x(2)-Ydelta).lt.0.6*dy .and. abs(x(3)-Zdelta).lt.0.6*dz) then
		  write(1,'(i10, 2g14.6, 9g12.4)') n, x(1), x(2), x(3), (C(s), s = 1,Ns)
	        else
		  write(1,'(i10, 2g14.6, 9g12.4)') n, x(1), x(2), x(3), (0.  , s = 1,Ns)
	        endif
	      else
! Co = const
	        write(1,'(i10 ,2g14.6 ,9g12.4)') n, x(1), x(2), x(3), (C(s), s = 1,Ns)
	      endif
	    enddo	! ix
	  enddo		! iy
	enddo		! iz
	write(1,'(a)') 'n, x, y, z, c(x,y,z|0, 1:Ns)'
	close(1)

! write elements, V, D, independent source and faces files
	open(1, file='elements.txt', status='unknown')
	open(2, file=       'D.txt', status='unknown')
	open(4, file=  'source.txt', status='unknown')
	if(.not. Darcy)  open(5, file=       'V.txt', status='unknown')
	open(7, file=  'faces.txt' , status='unknown')

	do iz = 1, Nz
	  do iy = 1, Ny
	    do ix = 1, Nx
	      e     = ix + (iy-1) *  Nx    + (iz-1) *  Nx *     Ny
	      ie(1) = ix + (iy-1) * (Nx+1) + (iz-1) * (Nx+1) * (Ny+1)
	      ie(2) = ie(1) + 1
	      ie(3) = ie(2) + 1 + Nx
	      ie(4) = ie(3) - 1
	      ie(5) = ie(1) + (Ny+1) * (Nx+1)
	      ie(6) = ie(2) + (Ny+1) * (Nx+1)
	      ie(7) = ie(3) + (Ny+1) * (Nx+1)
	      ie(8) = ie(4) + (Ny+1) * (Nx+1)
	      ie(9) = 1
	      write(1,'(i10, 9i10)') e, (ie(k), k=1,9)
	      write(2,'(i10, 9g12.4)') e, Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz
	      if(.not. Darcy)  write(5,'(i10, 3g12.4)') e, Vx, Vy, Vz

! 3D interpolation at all nodes of an orthogonal hexa element (x1,x2)*(y1,y2)*(z1,z2)
	      call INTERPOLV(Ns, Xmin, Ymin, Zmin, Nxx, Nyy, Nzz, xx, yy, zz, &
				fx, gy, hz, dx, dy, dz, ix, iy, iz, Sxy, ierror)
	      do s = 1, Ns
	        write (4,'(2i10, 8g12.4)') e, s, (Sxy(k,s), k=1,8)
	      enddo	! s

	      Gf = Gf+1
	      write(7,'(9i10)') Gf, e, 1, ie(1), ie(4), ie(8), ie(5), &
				    e -  1   , 4			! Lf=1 (west ) face
	      Gf = Gf+1
	      write(7,'(9i10)') Gf, e, 2, ie(1), ie(5), ie(6), ie(2), &
				    e - Nx   , 5			! Lf=2 (south) face
	      Gf = Gf+1
	      write(7,'(9i10)') Gf, e, 3, ie(1), ie(2), ie(3), ie(4), &
				    e - Nx*Ny, 6			! Lf=3 (bot  ) face

! add to faces file
	      if     (ix.eq.Nx)	then
	        Gf = Gf+1
		write(7,'(9i10)') Gf, e, 4, ie(2), ie(3), ie(7), ie(6), &
				      e   + 1  , 1			! Lf=4 (east ) face
	        endif
	      if(iy.eq.Ny)	then
	        Gf = Gf+1
		write(7,'(9i10)') Gf, e, 5, ie(4), ie(8), ie(7), ie(3), &
				      e + Nx   , 2			! Lf=5 (north) face
	      endif
	      if(iz.eq.Nz)	then
	        Gf = Gf+1
		write(7,'(9i10)') Gf, e, 6, ie(5), ie(6), ie(7), ie(8), &
				      e + Nx*Ny, 3			! Lf=6 (top  ) face
	     endif

	    enddo	! ix
	  enddo		! iy
	enddo		! iz
	write(1,'(a)') 'e, n1, n2, n3, n4, n5, n6, n7, n8, mat'
	close(1)
	write(2,'(a)') 'e, Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz'
	close(2)
	if(.not. Darcy) then
	  write(5,'(a)') 'e, Vx, Vy, Vz'
	  close(5)
	endif

	write(7,'(a)') 'global_face, element,  local_face, node1, node2, node3, node4, neighbour_element, neighbour_face'
	close(7)

	write (4,'(a)') 'e, s, (Sxy(k,s), k=1,8)'
	close(4)

! write BCs file
	open(1, file='BCs.txt', status='unknown')

	f = 1	! west (Ymin < y < Ymax, Zmin < z < Zmax, x=Xmin)
	ix = 1
	do iy = 1, Ny
	  coor1(1) = dy * (iy-1) + Ymin	! y at node 1
	  coor1(2) = dy * (iy  ) + Ymin	! y at node 2
	  dcoor1 = dy
	  do iz = 1, Nz
	    coor2(1) = dz * (iz-1) + Zmin	! z at node 1
	    coor2(2) = dz * (iz  ) + Zmin	! z at node 2
	    dcoor2 = dz
! 2D interpolation at all nodes of an orthogonal quad face
	    call INTERPOLS(Ns, Nxy1, Nxy2, f, xy1, xy2, vl1, vl2, dcoor1, dcoor2, coor1, coor2, BCv, ierror)
	    BCe = ix + (iy-1)*Nx + (iz-1)*Nx*Ny
	    write(1,'(2i10, a5, 18g12.4)') f, BCe, BCt(f), &
		  (BCv(f,s,1), BCv(f,s,2), BCv(f,s,3), BCv(f,s,4), s = 1,Ns)
	  enddo	! iz
	enddo	! iy

	f = 2	! south (Zmin < z < Zmax, Xmin < x < Xmax, y=Ymin)
	iy = 1
	do iz = 1, Nz
	  coor1(1) = dz * (iz-1) + Zmin	! z at node 1
	  coor1(2) = dz * (iz  ) + Zmin	! z at node 2
	  dcoor1 = dz
	  do ix = 1, Nx
	    coor2(1) = dx * (ix-1) + Xmin	! x at node 1
	    coor2(2) = dx * (ix  ) + Xmin	! x at node 2
	    dcoor2 = dx
! 2D interpolation at all nodes of an orthogonal quad face
	    call INTERPOLS(Ns, Nxy1, Nxy2, f, xy1, xy2, vl1, vl2, dcoor1, dcoor2, coor1, coor2, BCv, ierror)
	    BCe = ix + (iy-1)*Nx + (iz-1)*Nx*Ny
	    write(1,'(2i10, a5, 18g12.4)') f, BCe, BCt(f), &
		  (BCv(f,s,1), BCv(f,s,2), BCv(f,s,3), BCv(f,s,4), s = 1,Ns)
	  enddo	! ix
	enddo	! iz

	f = 3	! bottom (Xmin < x < Xmax, Ymin < y < Ymax, z=Zmin)
	iz = 1
	do ix = 1, Nx
	  coor1(1) = dx * (ix-1) + Xmin	! x at node 1
	  coor1(2) = dx * (ix  ) + Xmin	! x at node 2
	  dcoor1 = dx
	  do iy = 1, Ny
	    coor2(1) = dy * (iy-1) + Ymin	! y at node 1
	    coor2(2) = dy * (iy  ) + Ymin	! y at node 2
	    dcoor2 = dy
! 2D interpolation at all nodes of an orthogonal quad face
	    call INTERPOLS(Ns, Nxy1, Nxy2, f, xy1, xy2, vl1, vl2, dcoor1, dcoor2, coor1, coor2, BCv, ierror)
	    BCe = ix + (iy-1)*Nx + (iz-1)*Nx*Ny
	    write(1,'(2i10, a5, 18g12.4)') f, BCe, BCt(f), &
		  (BCv(f,s,1), BCv(f,s,2), BCv(f,s,3), BCv(f,s,4), s = 1,Ns)
	  enddo	! iy
	enddo	! ix

	f = 4	! east (Ymin < y < Ymax, Zmin < z < Zmax, x=Xmax)
	ix = Nx
	do iy = 1, Ny
	  coor1(1) = dy * (iy-1) + Ymin	! y at node 1
	  coor1(2) = dy * (iy  ) + Ymin	! y at node 2
	  dcoor1 = dy
	  do iz = 1, Nz
	    coor2(1) = dz * (iz-1) + Zmin	! z at node 1
	    coor2(2) = dz * (iz  ) + Zmin	! z at node 2
	    dcoor2 = dz
! 2D interpolation at all nodes of an orthogonal quad face
	    call INTERPOLS(Ns, Nxy1, Nxy2, f, xy1, xy2, vl1, vl2, dcoor1, dcoor2, coor1, coor2, BCv, ierror)
	    BCe = ix + (iy-1)*Nx + (iz-1)*Nx*Ny
	    write(1,'(2i10, a5, 18g12.4)') f, BCe, BCt(f), &
		  (BCv(f,s,1), BCv(f,s,2), BCv(f,s,3), BCv(f,s,4), s = 1,Ns)
	  enddo	! iz
	enddo	! iy

	f = 5	! north (Zmin < z < Zmax, Xmin < x < Xmax, y=Ymax)
	iy = Ny
	do iz = 1, Nz
	  coor1(1) = dz * (iz-1) + Zmin	! z at node 1
	  coor1(2) = dz * (iz  ) + Zmin	! z at node 2
	  dcoor1 = dz
	  do ix = 1, Nx
	    coor2(1) = dx * (ix-1) + Xmin	! x at node 1
	    coor2(2) = dx * (ix  ) + Xmin	! x at node 2
	    dcoor2 = dx
! 2D interpolation at all nodes of an orthogonal quad face
	    call INTERPOLS(Ns, Nxy1, Nxy2, f, xy1, xy2, vl1, vl2, dcoor1, dcoor2, coor1, coor2, BCv, ierror)
	    BCe = ix + (iy-1)*Nx + (iz-1)*Nx*Ny
	    write(1,'(2i10, a5, 18g12.4)') f, BCe, BCt(f), &
		  (BCv(f,s,1), BCv(f,s,2), BCv(f,s,3), BCv(f,s,4), s = 1,Ns)
	  enddo	! ix
	enddo	! iz

	f = 6	! top (Xmin < x < Xmax, Ymin < y < Ymax, z=Zmax)
	iz = Nz
	do ix = 1, Nx
	  coor1(1) = dx * (ix-1) + Xmin	! x at node 1
	  coor1(2) = dx * (ix  ) + Xmin	! x at node 2
	  dcoor1 = dx
	  do iy = 1, Ny
	    coor2(1) = dy * (iy-1) + Ymin	! y at node 1
	    coor2(2) = dy * (iy  ) + Ymin	! y at node 2
	    dcoor2 = dy
! 2D interpolation at all nodes of an orthogonal quad face
	    call INTERPOLS(Ns, Nxy1, Nxy2, f, xy1, xy2, vl1, vl2, dcoor1, dcoor2, coor1, coor2, BCv, ierror)
	    BCe = ix + (iy-1)*Nx + (iz-1)*Nx*Ny
	    write(1,'(2i10, a5, 18g12.4)') f, BCe, BCt(f), &
		  (BCv(f,s,1), BCv(f,s,2), BCv(f,s,3), BCv(f,s,4), s = 1,Ns)
	  enddo	! iy
	enddo	! ix

	write(1,'(a)') 'f, BCe, BCtype, (BCv(f,s,1), BCv(f,s,2), BCv(f,s,3), BCv(f,s,4),s = 1,Ns)'
	close(1)

! check if Dirichlet BC edges have the same values
!!	do f = 1,12
!!	  n = f+1		! neighbour face
!!	  if(f .eq. 6)	n = 1	! allow for f periodity
!! 	  Dedge = BCt(f) .eq. 'D' .and. BCt(n) .eq. 'D'
!!	  do s = 1, Ns
!!	    if     (f .eq. 1) then
!!	      Vf = vl1(f,s,Nxy1(1) )	! faces 1 & 2
!!	      Vn = vl1(n,s,1      )
!!	    else if(f .eq. 2) then
!!	      Vf = vl1(f,s,Nxy1(2) )	! faces 2 & 3
!!	      Vn = vl1(n,s,Nxy1(3) )
!!	    else if(f .eq. 3) then
!!	      Vf = vl1(f,s,1      )	! faces 3 & 4
!!	      Vn = vl1(n,s,Nxy1(4) )
!!	    else
!!	      Vf = vl1(f,s,1      )	! faces 4 & 1
!!	      Vn = vl1(n,s,1      )
!!	    endif
!!	    if( Dedge .and. Vf.ne.Vn) then	! Dirichlet BC edge
!!	      write (3,*) '!!! error !!! inconsistent Dirichlet BC edge values'
!!	      write (3,*) ' species s, face f and neighbour face n: s, f, n = ', s, f, n
!!	      write (3,*) ' f- and n- face values, Vf, Vn  = ', Vf, Vn
!!	      ierror = ierror + 1
!!	    endif
!!	  enddo		! s
!!	enddo		! f

! check for errors
	if(ierror .ne. 0)	then
	  write (3,*) '!!! ABORTED !!!'
	  stop
	endif
	close(3)

	end

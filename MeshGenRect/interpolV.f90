	subroutine INTERPOLV(Ns, Xmin, Ymin, Zmin, Nxx, Nyy, Nzz, xx, yy, zz, &
				fx, gy, hz, dx, dy, dz, ix, iy, iz, Sxy, ierror)
! 3D interpolation at all nodes of an orthogonal hexa element (x1,x2)*(y1,y2)*(z1,z2)

! s is species

	implicit none
	integer Nx, Ny, Nz, Ns, Nxx, Nyy, Nzz, ix, iy, iz, ierror
	real*8 Xmin, Ymin, Zmin, dx, dy, dz
	real*8 xx(Nxx), fx(Ns,Nxx), yy(Nyy), gy(Ns,Nyy), zz(Nzz), hz(Ns,Nzz), Sxy(8,Ns)

	integer i, k, m, n, s
	real*8 x1, x2, y1, y2, z1, z2, x(3), tol, xn(8), yn(8), zn(8)
	real*8 v1(Ns), v2(Ns), v3(Ns), v4(Ns), v5(Ns), v6(Ns), v7(Ns), v8(Ns)
	real*8 rr, ss, tt	! normalized coordinates

	tol = max(dx,dy) * 1.d-6	! geometric tolerance

! x1, x2, y1, y2
	x1 = dx * (ix-1) + Xmin
	y1 = dy * (iy-1) + Ymin
	z1 = dz * (iz-1) + Zmin
	x2 = x1 + dx
	y2 = y1 + dy
	z2 = z1 + dz

! nodal coordinates
	xn(1) = x1
	xn(2) = x2
	xn(3) = x2
	xn(4) = x1
	xn(5) = x1
	xn(6) = x2
	xn(7) = x2
	xn(8) = x1

	yn(1) = y1
	yn(2) = y1
	yn(3) = y2
	yn(4) = y2
	yn(5) = y1
	yn(6) = y1
	yn(7) = y2
	yn(8) = y2

	zn(1) = z1
	zn(2) = z1
	zn(3) = z1
	zn(4) = z1
	zn(5) = z2
	zn(6) = z2
	zn(7) = z2
	zn(8) = z2

! linear interpolation for the 8 nodes
	do k = 1,Nxx-1
	  do m = 1,Nyy-1
	    do n = 1,Nzz-1
	      do i = 1,8
	        if( xn(i) .ge. xx(k)-tol .and. xn(i) .le. xx(k+1)+tol .and. &
	            yn(i) .ge. yy(m)-tol .and. yn(i) .le. yy(m+1)+tol .and. &
	            zn(i) .ge. zz(n)-tol .and. zn(i) .le. zz(n+1)+tol ) then
! found - interpolate (k,m,n) interval

!!!!!!!!!!!!!!!!!!!!! from Mathematica !!!!!!!!!!!!!!!!!!!!!
! r = 2/(Xk2 - Xk1)*(x - Xk1) - 1
! s = 2/(Ym2 - Ym1)*(y - Ym1) - 1
! t = 2/(Zm2 - Zm1)*(z - Zm1) - 1
! v = 1/8*(1-r)*(1-s)*(1-t)*v1 + 1/8*(1+r)*(1-s)*(1-t)*v2 + 1/8*(1+r)*(1+s)*(1-t)*v3 + 1/8*(1-r)*(1+s)*(1-t)*v4 + 
!     1/8*(1-r)*(1-s)*(1+t)*v5 + 1/8*(1+r)*(1-s)*(1+t)*v6 + 1/8*(1+r)*(1+s)*(1+t)*v7 + 1/8*(1-r)*(1+s)*(1+t)*v8

!FullSimplify[v]
!!((v3 x - v4 x - v3 Xk1 + v4 Xk2) (y - Ym1) - v2 (x - Xk1) (y - Ym2) + 
!! v1 (x - Xk2) (y - Ym2))/((Xk1 - Xk2) (Ym1 - Ym2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! normalized coordinates
		  rr = 2. / (x2-x1) * (xn(i)-x1) - 1.	! r
		  ss = 2. / (y2-y1) * (yn(i)-y1) - 1.	! s
		  tt = 2. / (z2-z1) * (zn(i)-z1) - 1.	! t

	          do s = 1, Ns
	            v1(s) = fx(s,k  ) * gy(s,m  ) * hz(s,n  )
	            v2(s) = fx(s,k+1) * gy(s,m  ) * hz(s,n  )
	            v3(s) = fx(s,k+1) * gy(s,m+1) * hz(s,n  )
	            v4(s) = fx(s,k  ) * gy(s,m+1) * hz(s,n  )
	            v5(s) = fx(s,k  ) * gy(s,m  ) * hz(s,n+1)
	            v6(s) = fx(s,k+1) * gy(s,m  ) * hz(s,n+1)
	            v7(s) = fx(s,k+1) * gy(s,m+1) * hz(s,n+1)
	            v8(s) = fx(s,k  ) * gy(s,m+1) * hz(s,n+1)

!!		    Sxy(i,s) = ( ( v3(s) * xn(i) - v4(s) * xn(i) - v3(s)        * &
!!			    xx(k  ) + v4(s) * xx(k+1) ) * ( yn(i) - yy(m  ) )   - &
!!			    v2(s) * ( xn(i) - xx(k  ) ) * ( yn(i) - yy(m+1) )   + &
!!			    v1(s) * ( xn(i) - xx(k+1) ) * ( yn(i) - yy(m+1) ) ) / &
!!			    ( ( xx(k  ) - xx(k+1) ) * ( yy(m  ) - yy(m+1) ) )
		    Sxy(i,s) = 0.125 * (				  &
				(1.-rr) * (1.-ss) * (1.-tt) * v1(s)	+ &
				(1.+rr) * (1.-ss) * (1.-tt) * v2(s)	+ &
				(1.+rr) * (1.+ss) * (1.-tt) * v3(s)	+ &
				(1.-rr) * (1.+ss) * (1.-tt) * v4(s)	+ &
				(1.-rr) * (1.-ss) * (1.+tt) * v5(s)	+ &
				(1.+rr) * (1.-ss) * (1.+tt) * v6(s)	+ &
				(1.+rr) * (1.+ss) * (1.+tt) * v7(s)	+ &
				(1.-rr) * (1.+ss) * (1.+tt) * v8(s)	)
	          enddo	! s
	        endif
	      enddo	! i
	    enddo	! n
	  enddo		! m
	enddo		! k

	return
	end

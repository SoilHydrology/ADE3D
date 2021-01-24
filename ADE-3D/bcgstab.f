c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c         Basic Iterative Solvers with Reverse Communication           c
c----------------------------------------------------------------------c
c     This file currently has several basic iterative linear system    c
c     solvers. They are:                                               c
c     CG       -- Conjugate Gradient Method                            c
c     CGNR     -- Conjugate Gradient Method (Normal Residual equation) c
c     BCG      -- Bi-Conjugate Gradient Method                         c
c     DBCG     -- BCG with partial pivoting                            c
c     BCGSTAB  -- BCG stabilized                                       c
c     TFQMR    -- Transpose-Free Quasi-Minimum Residual method         c
c     FOM      -- Full Orthogonalization Method                        c
c     GMRES    -- Generalized Minimum RESidual method                  c
c     FGMRES   -- Flexible version of Generalized Minimum              c
c                 RESidual method                                      c
c     DQGMRES  -- Direct versions of Quasi Generalize Minimum          c
c                 Residual method                                      c
c----------------------------------------------------------------------c
c     They all have the following calling sequence:
c      subroutine solver(n, rhs, sol, ipar, fpar, w)
c      integer n, ipar(16)
c      real*8 rhs(n), sol(n), fpar(16), w(*)
c     Where
c     (1) 'n' is the size of the linear system,
c     (2) 'rhs' is the right-hand side of the linear system,
c     (3) 'sol' is the solution to the linear system,
c     (4) 'ipar' is an integer parameter array for the reverse
c     communication protocol,
c     (5) 'fpar' is an floating-point parameter array storing
c     information to and from the iterative solvers.
c     (6) 'w' is the work space (size is specified in ipar)
c
c     They are preconditioned iterative solvers with reverse
c     communication. The preconditioners can be applied from either
c     from left or right or both (specified by ipar(2), see below).
c
c     Author: Kesheng John Wu (kewu@mail.cs.umn.edu) 1993
c
c     NOTES:
c
c     (1) Work space required by each of the iterative solver
c     routines is as follows:
c       CG      == 5 * n
c       CGNR    == 5 * n
c       BCG     == 7 * n
c       DBCG    == 11 * n
c       BCGSTAB == 8 * n
c       TFQMR   == 11 * n
c       FOM     == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
c       GMRES   == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
c       FGMRES  == 2*n*(m+1) + (m+1)*m/2 + 3*m + 2 (m = ipar(5),
c                  default m=15)
c       DQGMRES == n + lb * (2*n+4) (lb=ipar(5)+1, default lb = 16)
c
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
c
c     (3) To use this set of routines under SPMD/MIMD program paradigm,
c     several things are to be noted: (a) 'n' should be the number of
c     vector elements of 'rhs' that is present on the local processor.
c     (b) if RHS(i) is on processor j, it is expected that SOL(i)
c     will be on the same processor, i.e. the vectors are distributed
c     to each processor in the same way. (c) the preconditioning and
c     stopping criteria specifications have to be the same on all
c     processor involved, ipar and fpar have to be the same on each
c     processor. (d) DISTDOT should be replaced by a distributed
c     dot-product function.
c
c     ..................................................................
c     Reverse Communication Protocols
c
c     When a reverse-communication routine returns, it could be either
c     that the routine has terminated or it simply requires the caller
c     to perform one matrix-vector multiplication. The possible matrices
c     that involve in the matrix-vector multiplications are:
c     A       (the matrix of the linear system),
c     A^T     (A transposed),
c     Ml^{-1} (inverse of the left preconditioner),
c     Ml^{-T} (inverse of the left preconditioner transposed),
c     Mr^{-1} (inverse of the right preconditioner),
c     Mr^{-T} (inverse of the right preconditioner transposed).
c     For all the matrix vector multiplication, v = A u. The input and
c     output vectors are supposed to be part of the work space 'w', and
c     the starting positions of them are stored in ipar(8:9), see below.
c
c     The array 'ipar' is used to store the information about the solver.
c     Here is the list of what each element represents:
c
c     ipar(1) -- status of the call/return.
c     A call to the solver with ipar(1) == 0 will initialize the
c     iterative solver. On return from the iterative solver, ipar(1)
c     carries the status flag which indicates the condition of the
c     return. The status information is divided into two categories,
c     (1) a positive value indicates the solver requires a matrix-vector
c     multiplication,
c     (2) a non-positive value indicates termination of the solver.
c     Here is the current definition:
c       1 == request a matvec with A,
c       2 == request a matvec with A^T,
c       3 == request a left preconditioner solve (Ml^{-1}),
c       4 == request a left preconditioner transposed solve (Ml^{-T}),
c       5 == request a right preconditioner solve (Mr^{-1}),
c       6 == request a right preconditioner transposed solve (Mr^{-T}),
c      10 == request the caller to perform stopping test,
c       0 == normal termination of the solver, satisfied the stopping
c            criteria,
c      -1 == termination because iteration number is greater than the
c            preset limit,
c      -2 == return due to insufficient work space,
c      -3 == return due to anticipated break-down / divide by zero,
c            in the case where Arnoldi procedure is used, additional
c            error code can be found in ipar(12), where ipar(12) is
c            the error code of orthogonalization procedure MGSRO:
c               -1: zero input vector
c               -2: input vector contains abnormal numbers
c               -3: input vector is a linear combination of others
c               -4: trianguler system in GMRES/FOM/etc. has nul rank
c      -4 == the values of fpar(1) and fpar(2) are both <= 0, the valid
c            ranges are 0 <= fpar(1) < 1, 0 <= fpar(2), and they can
c            not be zero at the same time
c      -9 == while trying to detect a break-down, an abnormal number is
c            detected.
c     -10 == return due to some non-numerical reasons, e.g. invalid
c            floating-point numbers etc.
c
c     ipar(2) -- status of the preconditioning:
c       0 == no preconditioning
c       1 == left preconditioning only
c       2 == right preconditioning only
c       3 == both left and right preconditioning
c
c     ipar(3) -- stopping criteria (details of this will be
c     discussed later).
c
c     ipar(4) -- number of elements in the array 'w'. if this is less
c     than the desired size, it will be over-written with the minimum
c     requirement. In which case the status flag ipar(1) = -2.
c
c     ipar(5) -- size of the Krylov subspace (used by GMRES and its
c     variants), e.g. GMRES(ipar(5)), FGMRES(ipar(5)),
c     DQGMRES(ipar(5)).
c
c     ipar(6) -- maximum number of matrix-vector multiplies, if not a
c     positive number the iterative solver will run till convergence
c     test is satisfied.
c
c     ipar(7) -- current number of matrix-vector multiplies. It is
c     incremented after each matrix-vector multiplication. If there
c     is preconditioning, the counter is incremented after the
c     preconditioning associated with each matrix-vector multiplication.
c
c     ipar(8) -- pointer to the input vector to the requested matrix-
c     vector multiplication.
c
c     ipar(9) -- pointer to the output vector of the requested matrix-
c     vector multiplication.
c
c     To perform v = A * u, it is assumed that u is w(ipar(8):ipar(8)+n-1)
c     and v is stored as w(ipar(9):ipar(9)+n-1).
c
c     ipar(10) -- the return address (used to determine where to go to
c     inside the iterative solvers after the caller has performed the
c     requested services).
c
c     ipar(11) -- the result of the external convergence test
c     On final return from the iterative solvers, this value
c     will be reflected by ipar(1) = 0 (details discussed later)
c
c     ipar(12) -- error code of MGSRO, it is
c                  1 if the input vector to MGSRO is linear combination
c                    of others,
c                  0 if MGSRO was successful,
c                 -1 if the input vector to MGSRO is zero,
c                 -2 if the input vector contains invalid number.
c
c     ipar(13) -- number of initializations. During each initilization
c                 residual norm is computed directly from M_l(b - A x).
c
c     ipar(14) to ipar(16) are NOT defined, they are NOT USED by
c     any iterative solver at this time.
c
c     Information about the error and tolerance are stored in the array
c     FPAR. So are some internal variables that need to be saved from
c     one iteration to the next one. Since the internal variables are
c     not the same for each routine, we only define the common ones.
c
c     The first two are input parameters:
c     fpar(1) -- the relative tolerance,
c     fpar(2) -- the absolute tolerance (details discussed later),
c
c     When the iterative solver terminates,
c     fpar(3) -- initial residual/error norm,
c     fpar(4) -- target residual/error norm,
c     fpar(5) -- current residual norm (if available),
c     fpar(6) -- current residual/error norm,
c     fpar(7) -- convergence rate,
c
c     fpar(8:10) are used by some of the iterative solvers to save some
c     internal information.
c
c     fpar(11) -- number of floating-point operations. The iterative
c     solvers will add the number of FLOPS they used to this variable,
c     but they do NOT initialize it, nor add the number of FLOPS due to
c     matrix-vector multiplications (since matvec is outside of the
c     iterative solvers). To insure the correct FLOPS count, the
c     caller should set fpar(11) = 0 before invoking the iterative
c     solvers and account for the number of FLOPS from matrix-vector
c     multiplications and preconditioners.
c
c     fpar(12:16) are not used in current implementation.
c
c     Whether the content of fpar(3), fpar(4) and fpar(6) are residual
c     norms or error norms depends on ipar(3). If the requested
c     convergence test is based on the residual norm, they will be
c     residual norms. If the caller want to test convergence based the
c     error norms (estimated by the norm of the modifications applied
c     to the approximate solution), they will be error norms.
c     Convergence rate is defined by (Fortran 77 statement)
c     fpar(7) = log10(fpar(3) / fpar(6)) / (ipar(7)-ipar(13))
c     If fpar(7) = 0.5, it means that approximately every 2 (= 1/0.5)
c     steps the residual/error norm decrease by a factor of 10.
c
c     ..................................................................
c     Stopping criteria,
c
c     An iterative solver may be terminated due to (1) satisfying
c     convergence test; (2) exceeding iteration limit; (3) insufficient
c     work space; (4) break-down. Checking of the work space is
c     only done in the initialization stage, i.e. when it is called with
c     ipar(1) == 0. A complete convergence test is done after each
c     update of the solutions. Other conditions are monitored
c     continuously.
c
c     With regard to the number of iteration, when ipar(6) is positive,
c     the current iteration number will be checked against it. If
c     current iteration number is greater the ipar(6) than the solver
c     will return with status -1. If ipar(6) is not positive, the
c     iteration will continue until convergence test is satisfied.
c
c     Two things may be used in the convergence tests, one is the
c     residual 2-norm, the other one is 2-norm of the change in the
c     approximate solution. The residual and the change in approximate
c     solution are from the preconditioned system (if preconditioning
c     is applied). The DQGMRES and TFQMR use two estimates for the
c     residual norms. The estimates are not accurate, but they are
c     acceptable in most of the cases. Generally speaking, the error
c     of the TFQMR's estimate is less accurate.
c
c     The convergence test type is indicated by ipar(3). There are four
c     type convergence tests: (1) tests based on the residual norm;
c     (2) tests based on change in approximate solution; (3) caller
c     does not care, the solver choose one from above two on its own;
c     (4) caller will perform the test, the solver should simply continue.
c     Here is the complete definition:
c      -2 == || dx(i) || <= rtol * || rhs || + atol
c      -1 == || dx(i) || <= rtol * || dx(1) || + atol
c       0 == solver will choose test 1 (next)
c       1 == || residual || <= rtol * || initial residual || + atol
c       2 == || residual || <= rtol * || rhs || + atol
c     999 == caller will perform the test
c     where dx(i) denote the change in the solution at the ith update.
c     ||.|| denotes 2-norm. rtol = fpar(1) and atol = fpar(2).
c
c     If the caller is to perform the convergence test, the outcome
c     should be stored in ipar(11).
c     ipar(11) = 0 -- failed the convergence test, iterative solver
c     should continue
c     ipar(11) = 1 -- satisfied convergence test, iterative solver
c     should perform the clean up job and stop.
c
c     Upon return with ipar(1) = 10,
c     ipar(8)  points to the starting position of the change in
c              solution Sx, where the actual solution of the step is
c              x_j = x_0 + M_r^{-1} Sx.
c              Exception: ipar(8) < 0, Sx = 0. It is mostly used by
c              GMRES and variants to indicate (1) Sx was not necessary,
c              (2) intermediate result of Sx is not computed.
c     ipar(9)  points to the starting position of a work vector that
c              can be used by the caller.
c
c     NOTE: the caller should allow the iterative solver to perform
c     clean up job after the external convergence test is satisfied,
c     since some of the iterative solvers do not directly
c     update the 'sol' array. A typical clean-up stage includes
c     performing the final update of the approximate solution and
c     computing the convergence information (e.g. values of fpar(3:7)).
c
c     NOTE: fpar(4) and fpar(6) are not set by the accelerators (the
c     routines implemented here) if ipar(3) = 999.
c
c     ..................................................................
c     Usage:
c
c     To start solving a linear system, the user needs to specify
c     first 6 elements of the ipar, and first 2 elements of fpar.
c     The user may optionally set fpar(11) = 0 if one wants to count
c     the number of floating-point operations. (Note: the iterative
c     solvers will only add the floating-point operations inside
c     themselves, the caller will have to add the FLOPS from the
c     matrix-vector multiplication routines and the preconditioning
c     routines in order to account for all the arithmetic operations.)
c
c     Here is an example:
c     ipar(1) = 0	! always 0 to start an iterative solver
c     ipar(2) = 2	! right preconditioning
c     ipar(3) = 1	! use convergence test scheme 1
c     ipar(4) = 10000	! the 'w' has 10,000 elements
c     ipar(5) = 10	! use *GMRES(10) (e.g. FGMRES(10))
c     ipar(6) = 100	! use at most 100 matvec's
c     fpar(1) = 1.0E-6	! relative tolerance 1.0E-6
c     fpar(2) = 1.0E-10 ! absolute tolerance 1.0E-10
c     fpar(11) = 0.0	! clearing the FLOPS counter
c
c     After the above specifications, one can start to call an iterative
c     solver, say BCG. Here is a piece of pseudo-code showing how it can
c     be done,
c
c 10   call bcg(n,rhs,sol,ipar,fpar,w)
c      if (ipar(1).eq.1) then
c         call amux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
c         goto 10
c      else if (ipar(1).eq.2) then
c         call atmux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
c         goto 10
c      else if (ipar(1).eq.3) then
c         left preconditioner solver
c         goto 10
c      else if (ipar(1).eq.4) then
c         left preconditioner transposed solve
c         goto 10
c      else if (ipar(1).eq.5) then
c         right preconditioner solve
c         goto 10
c      else if (ipar(1).eq.6) then
c         right preconditioner transposed solve
c         goto 10
c      else if (ipar(1).eq.10) then
c         call my own stopping test routine
c         goto 10
c      else if (ipar(1).gt.0) then
c         ipar(1) is an unspecified code
c      else
c         the iterative solver terminated with code = ipar(1)
c      endif
c
c     This segment of pseudo-code assumes the matrix is in CSR format,
c     AMUX and ATMUX are two routines from the SPARSKIT MATVEC module.
c     They perform matrix-vector multiplications for CSR matrices,
c     where w(ipar(8)) is the first element of the input vectors to the
c     two routines, and w(ipar(9)) is the first element of the output
c     vectors from them. For simplicity, we did not show the name of
c     the routine that performs the preconditioning operations or the
c     convergence tests.
c-----------------------------------------------------------------------
      subroutine bcgstab(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,8)
c-----------------------------------------------------------------------
c     BCGSTAB --- Bi Conjugate Gradient stabilized (BCGSTAB)
c     This is an improved BCG routine. (1) no matrix transpose is
c     involved. (2) the convergence is smoother.
c
c
c     Algorithm:
c     Initialization - r = b - A x, r0 = r, p = r, rho = (r0, r),
c     Iterate -
c     (1) v = A p
c     (2) alpha = rho / (r0, v)
c     (3) s = r - alpha v
c     (4) t = A s
c     (5) omega = (t, s) / (t, t)
c     (6) x = x + alpha * p + omega * s
c     (7) r = s - omega * t
c     convergence test goes here
c     (8) beta = rho, rho = (r0, r), beta = rho * alpha / (beta * omega)
c         p = r + beta * (p - omega * v)
c
c     in this routine, before successful return, the fpar's are
c     fpar(3) == initial (preconditionied-)residual norm
c     fpar(4) == target (preconditionied-)residual norm
c     fpar(5) == current (preconditionied-)residual norm
c     fpar(6) == current residual norm or error
c     fpar(7) == current rho (rhok = <r, r0>)
c     fpar(8) == alpha
c     fpar(9) == omega
c
c     Usage of the work space W
c     w(:, 1) = r0, the initial residual vector
c     w(:, 2) = r, current residual vector
c     w(:, 3) = s
c     w(:, 4) = t
c     w(:, 5) = v
c     w(:, 6) = p
c     w(:, 7) = tmp, used in preconditioning, etc.
c     w(:, 8) = delta x, the correction to the answer is accumulated
c               here, so that the right-preconditioning may be applied
c               at the end
c-----------------------------------------------------------------------
c     external routines used
c
      real*8 distdot
      logical stopbis, brkdn
      external distdot, stopbis, brkdn
c
      real*8 one
      parameter(one=1.0D0)
c
c     local variables
c
      integer i
      real*8 alpha,beta,rho,omega
      logical lp, rp
      save lp, rp
c
c     where to go
c
      if (ipar(1).gt.0) then
         goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110) ipar(10)
      else if (ipar(1).lt.0) then
         goto 900
      endif
c
c     call the initialization routine
c
      call bisinit(ipar,fpar,8*n,1,lp,rp,w)
      if (ipar(1).lt.0) return
c
c     perform a matvec to compute the initial residual
c
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + n
      do i = 1, n
         w(i,1) = sol(i)
      enddo
      ipar(10) = 1
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      do i = 1, n
         w(i,1) = rhs(i) - w(i,2)
      enddo
      fpar(11) = fpar(11) + n
      if (lp) then
         ipar(1) = 3
         ipar(10) = 2
         return
      endif
c
 20   if (lp) then
         do i = 1, n
            w(i,1) = w(i,2)
            w(i,6) = w(i,2)
         enddo
      else
         do i = 1, n
            w(i,2) = w(i,1)
            w(i,6) = w(i,1)
         enddo
      endif
c
      fpar(7) = distdot(n,w,1,w,1)
      fpar(11) = fpar(11) + 2 * n
      fpar(5) = sqrt(fpar(7))
      fpar(3) = fpar(5)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2 * n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
      endif
      if (ipar(3).ge.0) fpar(6) = fpar(5)
      if (ipar(3).ge.0 .and. fpar(5).le.fpar(4) .and.
     +     ipar(3).ne.999) then
         goto 900
      endif
c
c     beginning of the iterations
c
c     Step (1), v = A p
 30   if (rp) then
         ipar(1) = 5
         ipar(8) = 5*n+1
         if (lp) then
            ipar(9) = 4*n + 1
         else
            ipar(9) = 6*n + 1
         endif
         ipar(10) = 3
         return
      endif
c
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 5*n+1
      endif
      if (lp) then
         ipar(9) = 6*n + 1
      else
         ipar(9) = 4*n + 1
      endif
      ipar(10) = 4
      return
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 4*n + 1
         ipar(10) = 5
         return
      endif
c
 60   ipar(7) = ipar(7) + 1
c
c     step (2)
      alpha = distdot(n,w(1,1),1,w(1,5),1)
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(alpha, ipar)) goto 900
      alpha = fpar(7) / alpha
      fpar(8) = alpha
c
c     step (3)
      do i = 1, n
         w(i,3) = w(i,2) - alpha * w(i,5)
      enddo
      fpar(11) = fpar(11) + 2 * n
c
c     Step (4): the second matvec -- t = A s
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = n+n+1
         if (lp) then
            ipar(9) = ipar(8)+n
         else
            ipar(9) = 6*n + 1
         endif
         ipar(10) = 6
         return
      endif
c
 70   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = n+n+1
      endif
      if (lp) then
         ipar(9) = 6*n + 1
      else
         ipar(9) = 3*n + 1
      endif
      ipar(10) = 7
      return
 80   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 3*n + 1
         ipar(10) = 8
         return
      endif
 90   ipar(7) = ipar(7) + 1
c
c     step (5)
      omega = distdot(n,w(1,4),1,w(1,4),1)
      fpar(11) = fpar(11) + n + n
      if (brkdn(omega,ipar)) goto 900
      omega = distdot(n,w(1,4),1,w(1,3),1) / omega
      fpar(11) = fpar(11) + n + n
      if (brkdn(omega,ipar)) goto 900
      fpar(9) = omega
      alpha = fpar(8)
c
c     step (6) and (7)
      do i = 1, n
         w(i,7) = alpha * w(i,6) + omega * w(i,3)
         w(i,8) = w(i,8) + w(i,7)
         w(i,2) = w(i,3) - omega * w(i,4)
      enddo
      fpar(11) = fpar(11) + 6 * n + 1
c
c     convergence test
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 7*n + 1
         ipar(9) = 6*n + 1
         ipar(10) = 9
         return
      endif
      if (stopbis(n,ipar,2,fpar,w(1,2),w(1,7),one))  goto 900
 100  if (ipar(3).eq.999.and.ipar(11).eq.1) goto 900
c
c     step (8): computing new p and rho
      rho = fpar(7)
      fpar(7) = distdot(n,w(1,2),1,w(1,1),1)
      omega = fpar(9)
      beta = fpar(7) * fpar(8) / (fpar(9) * rho)
      do i = 1, n
         w(i,6) = w(i,2) + beta * (w(i,6) - omega * w(i,5))
      enddo
      fpar(11) = fpar(11) + 6 * n + 3
      if (brkdn(fpar(7),ipar)) goto 900
c
c     end of an iteration
c
      goto 30
c
c     some clean up job to do
c
 900  if (rp) then
         if (ipar(1).lt.0) ipar(12) = ipar(1)
         ipar(1) = 5
         ipar(8) = 7*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      endif
 110  if (rp) then
         call tidycg(n,ipar,fpar,sol,w(1,7))
      else
         call tidycg(n,ipar,fpar,sol,w(1,8))
      endif
c
      return
c-----end-of-bcgstab
      end


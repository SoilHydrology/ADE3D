subroutine LEGENDRE_SET ( Ng, Xg, Wg )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! copied from example.zip, quadrule.f90
! with the following changes:
! norder	-> Ng
! xtab		-> Xg
! weight	-> Wg
! truncated after Ng = 8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!*******************************************************************************
!
!! LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.0D+00
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) dX
!
!  Approximate integral:
!
!    Sum ( 1 <= I <= Ng ) Wg(I) * F ( Xg(I) )
!
!  Precision:
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*Ng-1).
!
!  Note:
!
!    The abscissas of the rule are the zeroes of the Legendre polynomial
!    P(Ng)(X).
!
!    The integral produced by a Gauss-Legendre rule is equal to the
!    integral of the unique polynomial of degree Ng-1 which
!    agrees with the function at the Ng abscissas of the rule.
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    MacMillan, 1962.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    18 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer Ng, the order of the rule.
!    Ng must be between 1 and 8.
!
!    Output, double precision Xg(Ng), the abscissas of the rule.
!
!    Output, double precision Wg(Ng), the weights of the rule.
!    The weights are positive, symmetric and should sum to 2.
!
  implicit none
!
  integer Ng
!
  double precision Xg(Ng)
  double precision Wg(Ng)
!
!	write(idbg,'(a)') ' --- LEGENDRE ---'	! ### TEMPORARY ###
 if ( Ng == 1 ) then

    Xg(1) =   0.0D+00

    Wg(1) = 2.0D+00

  else if ( Ng == 2 ) then

    Xg(1) = - 0.577350269189625764509148780502D+00
    Xg(2) =   0.577350269189625764509148780502D+00

    Wg(1) = 1.0D+00
    Wg(2) = 1.0D+00

  else if ( Ng == 3 ) then

    Xg(1) = - 0.774596669241483377035853079956D+00
    Xg(2) =   0.0D+00
    Xg(3) =   0.774596669241483377035853079956D+00

    Wg(1) = 5.0D+00 / 9.0D+00
    Wg(2) = 8.0D+00 / 9.0D+00
    Wg(3) = 5.0D+00 / 9.0D+00

  else if ( Ng == 4 ) then

    Xg(1) = - 0.861136311594052575223946488893D+00
    Xg(2) = - 0.339981043584856264802665759103D+00
    Xg(3) =   0.339981043584856264802665759103D+00
    Xg(4) =   0.861136311594052575223946488893D+00

    Wg(1) = 0.347854845137453857373063949222D+00
    Wg(2) = 0.652145154862546142626936050778D+00
    Wg(3) = 0.652145154862546142626936050778D+00
    Wg(4) = 0.347854845137453857373063949222D+00

  else if ( Ng == 5 ) then

    Xg(1) = - 0.906179845938663992797626878299D+00
    Xg(2) = - 0.538469310105683091036314420700D+00
    Xg(3) =   0.0D+00
    Xg(4) =   0.538469310105683091036314420700D+00
    Xg(5) =   0.906179845938663992797626878299D+00

    Wg(1) = 0.236926885056189087514264040720D+00
    Wg(2) = 0.478628670499366468041291514836D+00
    Wg(3) = 0.568888888888888888888888888889D+00
    Wg(4) = 0.478628670499366468041291514836D+00
    Wg(5) = 0.236926885056189087514264040720D+00

  else if ( Ng == 6 ) then

    Xg(1) = - 0.932469514203152027812301554494D+00
    Xg(2) = - 0.661209386466264513661399595020D+00
    Xg(3) = - 0.238619186083196908630501721681D+00
    Xg(4) =   0.238619186083196908630501721681D+00
    Xg(5) =   0.661209386466264513661399595020D+00
    Xg(6) =   0.932469514203152027812301554494D+00

    Wg(1) = 0.171324492379170345040296142173D+00
    Wg(2) = 0.360761573048138607569833513838D+00
    Wg(3) = 0.467913934572691047389870343990D+00
    Wg(4) = 0.467913934572691047389870343990D+00
    Wg(5) = 0.360761573048138607569833513838D+00
    Wg(6) = 0.171324492379170345040296142173D+00

  else if ( Ng == 7 ) then

    Xg(1) = - 0.949107912342758524526189684048D+00
    Xg(2) = - 0.741531185599394439863864773281D+00
    Xg(3) = - 0.405845151377397166906606412077D+00
    Xg(4) =   0.0D+00
    Xg(5) =   0.405845151377397166906606412077D+00
    Xg(6) =   0.741531185599394439863864773281D+00
    Xg(7) =   0.949107912342758524526189684048D+00

    Wg(1) = 0.129484966168869693270611432679D+00
    Wg(2) = 0.279705391489276667901467771424D+00
    Wg(3) = 0.381830050505118944950369775489D+00
    Wg(4) = 0.417959183673469387755102040816D+00
    Wg(5) = 0.381830050505118944950369775489D+00
    Wg(6) = 0.279705391489276667901467771424D+00
    Wg(7) = 0.129484966168869693270611432679D+00

  else if ( Ng == 8 ) then

    Xg(1) = - 0.960289856497536231683560868569D+00
    Xg(2) = - 0.796666477413626739591553936476D+00
    Xg(3) = - 0.525532409916328985817739049189D+00
    Xg(4) = - 0.183434642495649804939476142360D+00
    Xg(5) =   0.183434642495649804939476142360D+00
    Xg(6) =   0.525532409916328985817739049189D+00
    Xg(7) =   0.796666477413626739591553936476D+00
    Xg(8) =   0.960289856497536231683560868569D+00

    Wg(1) = 0.101228536290376259152531354310D+00
    Wg(2) = 0.222381034453374470544355994426D+00
    Wg(3) = 0.313706645877887287337962201987D+00
    Wg(4) = 0.362683783378361982965150449277D+00
    Wg(5) = 0.362683783378361982965150449277D+00
    Wg(6) = 0.313706645877887287337962201987D+00
    Wg(7) = 0.222381034453374470544355994426D+00
    Wg(8) = 0.101228536290376259152531354310D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of Ng = ', Ng
    write ( *, '(a)' ) '  Legal values are 1 to 8.'
    stop

  end if

  return
end


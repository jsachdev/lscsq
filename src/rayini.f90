SUBROUTINE lscsq_RayIni(RayIniErr)
!
!     Integration is in coords  R     Z     Phi
!     Starting uses             Rad   Pol   Phi (no B in Rad direction)
!
!     Kr       1   **  Bz  Br  **  ** Krad **
!         =  ----  *            *  *        *
!     Kz     Bpol  ** -Br  Bz  **  ** Kpol **
!
!     Kpar = (KtorBtor + KpolBpol)/B
!     Krad2 = Kperp2 - (KpolBtor - KtorBpol)**2/B**2
!     Krad is positive for slow wave, negative for fast.
!     Bpol has the sign of (Bz(R-Rmaj) - BrZ) in order to keep the
!     sign of Kr, Kz, Krad straight.  This should work except in
!     quite pathological shapes.

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : pi, twopi, vc, deg2rad
  use lscsq_mod, only: y,Epar,Eper,d1,d2,d4,woc2,woc4
  use lscsq_mod, only: lfast, begin
  use lscsq_mod, only: pe2min, Exy, izind, rzind
  use lscsq_mod, only: enpar, enpol
  use lscsq_mod, only: thet0, enth, ngrps, fghz, omega
  use lscsq_mod, only: pofray, nparry, nperry, vpar, lstop
  use lscsq_mod, only: timery, neofry, rtPsRy, Bthray, Bphray
  use lscsq_mod, only: delpsi,Powrry,RofRay,zofray
  use lscsq_mod, only: Ezsq,distry,detrry
  use lscsq_mod, only: iray, nzones, npsij 
  use lscsq_mod, only: sleave, tleave, senter, tenter, psiary
  use lscsq_mod, only: lh_const !, lh_out
  use lscsq_mod, only: lh_inp

  implicit none

  integer, intent(out) :: RayIniErr
  integer :: izindold
  real(fp) :: rzindold

  real(fp) :: Kpar2, Kpol, Krad, Ktor
  complex(fp), dimension(3) :: zc
  real(fp), dimension(4) :: Azplr

  real(fp) :: r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
  integer ::  i
  integer :: dummyer=0
  real(fp) :: Bpol,Bpol2,Btot2, lscsq_disprela,                            &
              det, fast,plas,Qpar,Sig,slow,try,                          &
              zi1,zi2,zi3,zr1,zr2,zr3, discrimt
  real(fp) :: CosT,drad,rad,rstar,SinT,t 
  real(fp) :: SEARCHINCR = 5.0e-03_fp

  real(fp):: RzindNew

  real(fp) :: re41
  real(fp) :: Btot, Bphi 
  real(fp) :: Kpar, Kper, kper2 

  T = enth * deg2rad    
  
  RayIniErr = 0
  Rad = lh_inp%Rmax-lh_inp%Raxis
  If ( Cos(T) .LT. 0.0_fp ) Rad=lh_inp%Raxis-lh_inp%Rmin
  CosT = Abs(Cos(T))+1.0e-20_fp
  SinT = Abs(Sin(T))+1.0e-20_fp
  Rad  = min(Rad/CosT, lh_inp%Zmax/SinT)
  rstar= -SEARCHINCR
  drad = abs(rstar)

5     do I=1,10000
     Rad = Rad - drad
     r    = lh_inp%Raxis + Rad*Cos(T)
     z    =        Rad*Sin(T)
     Kpol = omega*enpol/vc
     Ktor = omega*enpar/vc
     if (enpar.eq.0.0_fp) then
        RayIniErr = 2
        CALL lscsq_LSCwarn(' enpar == 0 in RayIni ')
        return
     endif
     CALL lscsq_plasma2d(r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)

     If (pe2.LE.pe2min) cycle    

      Bpol2= Br*Br + Bz*Bz
      Btot2= Bpol2 + (RBphi/r)**2
      Bpol = sqrt(Bpol2)
   
      Sig  = (r-lh_inp%Raxis)*Bz - z*Br
      If (Sig.LT.0.0_fp) Bpol = -Bpol
      Kpar2 = (Kpol*Bpol+Ktor*RBphi/R)**2/Btot2
      call lscsq_Eps( r, z, Kpar2, 0.0_fp)

10    Azplr(4) = -(aio/fghz(iray)**4 + ael)
      Azplr(3) =  Eper
      Qpar = Kpar2 - woc2*Eper
      Azplr(2) = Qpar*(Epar+Eper) +  woc2*Exy**2
      Azplr(1) = Epar*( Qpar**2 - woc4*Exy**2 )
 
      ! Estimate the roots far from lmc
      discrimt = Azplr(2)**2-4.0_fp*Azplr(3)*Azplr(1)
      if ( discrimt .LE. 0.0_fp ) then
        discrimt = 0.0_fp
        RayIniErr= 1
        CALL lscsq_LSCwarn(' no accessibity found for this ray ')
        return
      endif
 
      CALL lscsq_ZPLRCnr(3,Azplr,ZC)
      ! write(0,*) 'Azplr ', Azplr
      ! write(0,*) 'zc ', zc
!     Finds Zeroes of a Polynomial with Laguerre's method if
!     Real Coefficients
!     using Numerical Recipes code so we dont depend on IMSL
      zr1 = REAL(zc(1),kind=fp)
      zr2 = REAL(zc(2),kind=fp)
      zr3 = REAL(zc(3),kind=fp)
      ! zi1 = AIMAG(zc(1))
      ! zi2 = AIMAG(zc(2))
      ! zi3 = AIMAG(zc(3))

!     The normal case is that all roots are real, and we see if we can start.
22    continue
      plas =  max (zr1,zr2,zr3)
      fast =  min (zr1,zr2,zr3)
      slow = -1.0e30_fp 
      if (zr1.LT.plas .and. zr1.GT.fast) slow = zr1
      if (zr2.LT.plas .and. zr2.GT.fast) slow = zr2
      if (zr3.LT.plas .and. zr3.GT.fast) slow = zr3
      if (slow .EQ. -1.0e+30_fp) cycle    
      if (lfast .EQ. 1) try = fast
      if (lfast .EQ. 0) try = slow
24    if (try .LT. 0.0_fp ) cycle    
      try = try - ( Kpol*RBphi/R - Ktor*Bpol )**2/Btot2
      if ( try .LT. 0.0_fp ) cycle    
      try = sqrt(try)

      krad = try
      if(lfast .EQ. 1) krad = - try

      go to 35
  enddo
  RayIniErr =2
  CALL lscsq_LSCwarn(' cant find a starting point for this ray')
  return
35    continue

  y(1)  =  r
  y(2)  =  z
  y(3)  =  0.0_fp
  y(4)  =  ( Krad*Bz + Kpol*Br ) / Bpol
  y(5)  =  (-Krad*Br + Kpol*Bz ) / Bpol
  y(6)  =  Ktor*r
  y(7)  =  0.0_fp
  ! to begin the time as fn of distance
  y(8) = begin
  det = lscsq_DispRela(y(1),y(2),y(4),y(5),y(6))

  call lscsq_ftion(DummyEr)
!     only to initialize dkper & wdddw for damping calculation

  Kper=sqrt(abs(Y(4)**2+Y(5)**2+(Y(6)/Y(1))**2-Kpar2))

  RzindOld = (Psi-lh_inp%plflx(1))/DelPsi + 1.5_fp
  RE41 = RzindOld
  IzindOld = int(RE41)
!  izindold = minloc(abs(psiary-psi),1)
 
  izind(1,iray)  = IzindOld
  rzind(1,iray)  = RzindOld
  PowrRy(1,iray) = 1.0_fp
  RofRay(1,iray) = y(1)
  ZofRay(1,iray) = y(2)
  PofRay(1,iray) = y(3)
  NparRy(1,iray) = sqrt(Kpar2)/sqrt(woc2)
  NperRy(1,iray) = Kper/sqrt(woc2)
  rtPsRy(1,iray) = sqrt((psi-lh_inp%plflx(1))/(lh_inp%plflx(npsij)-lh_inp%plflx(1)))
  TimeRy(1,iray) = y(7)
  DistRy(1,iray) = y(8)
  NeofRy(1,iray) = pe2/lh_const%Pe2Fac*1.0e+14_fp
  BthRay(1,iray) = sqrt(Br**2+Bz**2)
  BphRay(1,iray) = RBphi/RofRay(1,iray)
  DetrRy(1,iray) = det/max(abs(d1),abs(d2),abs(d4))
  tenter(1,iray) = y(7)
  senter(1,iray) = y(8)

end subroutine lscsq_rayini
!                                                                      |
SUBROUTINE lscsq_zplrcnr (degree, coefr, zeroc)

  use iso_c_binding, only : fp => c_double
  implicit none

!     finds Zeros of Polynomials with Laguerre's method assuming
!     Real Coefficients by calling a routine from
!     Numerical Recipes.
!     The calling convention is set up to look like IMSL CALL lscsq_ZPLRC.
!     f(x) = 0 = coefr(1) + coefr(2)*x + coefr(3)*x**2 ...etc
  integer :: degree, i, polish
  integer, parameter :: degmax=10
  real(fp), dimension(degree+1) :: coefr
  complex(fp), dimension(degmax+1) :: coefc
  complex(fp), dimension(degree) :: zeroc

  ! zroots needs complex coef's
  do i = 1, degree+1
     coefc(i) = cmplx ( coefr(i) , 0.0_fp , kind=fp)
  enddo
 
  polish = 1
  CALL lscsq_zrootsnr ( coefc, degree, zeroc, polish )
 
end subroutine lscsq_zplrcnr
!
!     ------------------------------------------------------------------
!
SUBROUTINE lscsq_laguernr ( coef, degree, x, epsilon, polish )

  use iso_c_binding, only : fp => c_double
  implicit none

!     Given the  DEGREE  and  DEGREE+1  complex COEF's of the polynomial
!     Sum COEF(i) X**(i-1) and given  epsilon the desired fractional
!     accuracy, and given a complex value  X , this routine improves  X
!     by Laguerre's method until it converges to a root of the given
!     polynomial.  For normal use  POLISH  should be input as 0 (false).
!     When  POLISH  is 1 (true) the routine ignores  EPSILON  and
!     instead attempts to improve  X  (assumed to be a good initial
!     guess) to the achievable roundoff limit.
!     Ref: Numerical Recipes in Fortran, page 264.
 
  integer :: degree, iter, j, polish
  integer :: maxit=100
  complex(fp) :: czero=(0.0_fp,0.0_fp)

  complex(fp) :: coef(degree+1)
  complex(fp) :: x, dx, x1,b,d,f,g,h,sq,gp,gm,g2
  real(fp) :: epss = 6.0e-7_fp

  real(fp) :: epsilon, err, abx, cdx

  do iter = 1, MAXIT
     b   = coef ( degree+1 )
     err =  abs ( b )
     d   = CZERO
     f   = CZERO
     abx =  abs ( x )
!                                       Efficient computation of the
!                                       polynomial and its first 2
!                                       derivatives
     do j = degree, 1, -1
        f   = x*f + d
        d   = x*d + b
        b   = x*b + coef ( j )
        err =  abs ( b ) + abx * err
     enddo
!                                       Estimate of roundoff error in
!                                       evaluating polynomial
     err = EPSS * err
     if (  abs ( b ) .LE. err ) then
        return
     ! The generic case: use Laguerre's formula
     else
        g  = d / b
        g2 = g * g
        h  = g2 - 2.0_fp * f / b
        sq =  sqrt ( real(degree-1,kind=fp) * (real(degree,kind=fp)*h - g2) )
        gp = g + sq
        gm = g - sq
        if (  abs (gp) .LT.  abs (gm) ) gp = gm
        dx = real(degree,kind=fp) / gp
     endif
     x1    = x - dx
     if ( x .EQ. x1 ) return
     x     = x1
     cdx   =  abs ( dx )
     if ( polish .EQ. 0 ) then
        if ( cdx .LE. epsilon* abs ( x ) ) return
     endif
  enddo     
 
  CALL lscsq_LSCstop( ' too many iterations in LAGUER root finder ')

end subroutine lscsq_laguernr 
!
!     ------------------------------------------------------------------
!
subroutine lscsq_zrootsnr ( coef, degree, roots, polish )

  use iso_c_binding, only : fp => c_double
  implicit none

!!!!!!     EXTERNAL laguernr
!     Given the DEGREE and the  DEGREE+1  complex COEF's of the polynomial
!     Sum COEF(i) x**(i-1) this routine successively calls LAGUER and
!     finds all DEGREE complex ROOTS.  The integer variable POLISH
!     should be input as 1 (true) if polishing is desired, of 0 (false)
!     if the roots will be subsequently polished by other means.
!     Ref: Numerical Recipes in Fortran, page 265.
!
      INTEGER degree, i, j, jj, polish
  integer, parameter :: maxdegre=101
  real(fp) :: epsilon = 1.0e-6_fp
  complex(fp) :: coef(degree+1), x, b, c, roots(degree), defl(MAXDEGRE)
!
!                                       Copy coef's for successive deflation

  defl(1:degree+1) = coef(1:degree+1)
!
  do j = degree, 1, -1
     ! Start at 0 to favor smallest remaining root
     x = CMPLX ( 0.0_fp, 0.0_fp, kind=fp )
     CALL lscsq_laguernr ( defl, j, x, EPSILON, 0 )
     if (abs(REAL(AIMAG(x),kind=fp)).LE.(2.0_fp*EPSILON**2*abs(REAL(x,kind=fp))) )      &
        x = cmplx(REAL(x,kind=fp),0.0_fp , kind=fp)
     roots( j ) = x
     b = defl( j+1 )
     ! Forward deflation
     do jj = j, 1, -1
        c = defl(jj)
        defl(jj) = b
        b = x*b + c
     enddo
  enddo
 
  if ( polish .EQ. 1 ) then
     ! Polish roots using undeflated coefs
     do j = 1, degree
        CALL lscsq_laguernr ( coef, degree, roots(j), EPSILON, 1 )
     enddo
  endif
 
  ! Sort roots by real part by straight insertion
  do j = 2, degree
     x = roots (j)
     do i = j-1, 1, -1
        if ( REAL( roots(i),kind=fp ) .LE. REAL( x,kind=fp ) ) go to 10
        roots ( i+1 ) = roots ( i)
     enddo
     i=0
10    roots ( i+1 ) = x
  enddo        

end subroutine lscsq_zrootsnr 
!                                                                      |

!

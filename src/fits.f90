!
!     Reference: ``Current in wave-driven plasmas,'' by
!     Charles F. F. Karney and Nathaniel J. Fisch,
!     Phys. Fluids {\bf 29} 180-192 (1986)
!
!     Fits.F written March 1991 by D. W. Ignat.  All DATA statements
!     and logic based on the reference.
!     Copyright D. W. Ignat, E. J. Valeo, C. F. F. Karney, N. J. Fisch,
!     Princeton University, Plasma Physics Laboratory
!
!     These are memonics for the numerical constants:
!     R         Run-away probability
!     W         Energy imparted to the electric field by stopped electrons
!     WP        dW/du /u
!     p         mu is plus:  +1
!     m         mu is minus: -1
!     Z         effective charge
!
!     These are other usages:
!     RunProb   Run-away probability, or the R(u,mu) function of reference.
!     WsloDwn   Energy imparted to electric field by an electron as it
!               slows down, or the W_s(u,mu) function of reference.
!     WsloPrm   \partial W_s/\partial u / u function of reference.
!               This is the ratio P_{\rm el} to P_{\rm in} , or ratio of
!               power coupled from the rf source into electromagnetic energy
!               to the rf power absobed by the electrons.
!               Very roughly, Pel/Pin rises from 0 to 1 as u_par goes from
!               0 to 5, with most rise as u_par goes from 0 to 3.
!               See Fig. 7a of the reference.
!
!     u         velocity normalized to run-away velocity (magnitude!)
!     Mu        direction cosine, v-par/v-total
!     Z         Z_effective
!     x         u*u, or (u-1) in RunProb
!
!
!     Review of definitions and conventions of the reference:
!     q         carries sign of electron charge == - e
!     E         is parallel to B, and in the positive direction
!     v_par     is positive in the direction of E and B
!     \Gamma == \ln \Lambda \frac{n q^4}{4\pi \epsilon_0^2 m^2}
!     v_r    == -sign(qE) \sqrt{m\Gamma / \abs{qE} } (a positive number)
!     Dreicer Velocity v_D = - \sqrt{2 + Z} v_r
!     u         v /(\abs{ v_r})
!     u_\par    v_\par / v_r
!     u_\perp   v_\perp / \abs{v_r}
!
!     The following numbered lines are a copy of the MACSYMA fit equation,
!     kept here for redundancy.
!c01  WsloDwn = Mu1 * U04 / (5.+Z)
!c02 ^   -    (2. + Z + 3.*Mu2)                                * U06 /
!c03 ^                                     ( 3.*(3.+Z)*(5.+Z) )
!c04 ^   + 2.*( (24.+19.*Z+3.*Z**2)*Mu1 + (9.+Z)*Mu3 )         * U08 /
!c05 ^                       ( (3.+Z)*(5.+Z)*(7.+3.*Z)*(9.+Z) )
!c06 ^  -(
!c07 ^    (1041.+1864.*Z + 1189.*Z**2 + 316.*Z**3 + 30.*Z**4 ) * Mu0
!c08 ^   +( 417.+ 497.*Z +  181.*Z**2 +  21.*Z**3) * 10.       * Mu2
!c09 ^   +  (9.+Z)*(13.+3.*Z)                      *  5.       * Mu4
!c10 ^   )                                                     * U10 /
!c11 ^      ( 5.*(2.+Z)*(3.+Z)*(5.+Z)*(7.+3.*Z)*(9.+Z)*(13.+3.*Z) )
!
!
!     -----------------------------------------------------------------|-------
!
      FUNCTION lscsq_RunProb ( uGiven , MuGiven , ZGiven )
use iso_c_binding, only : fp => c_double
implicit none
      REAL(fp) lscsq_RunProb
      REAL(fp)                                                              &
     &        u , Mu , Z , uGiven, MuGiven, ZGiven,                     &
     &        x
      real(fp)    Y
      real(fp)                                                              &
     &     RpZ01a0,  RpZ01a1,  RpZ01a2,  RpZ01a3,  RpZ01b2,  RpZ01b3 ,  &
     &     RpZ02a0,  RpZ02a1,  RpZ02a2,  RpZ02a3,  RpZ02b2,  RpZ02b3 ,  &
     &     RpZ05a0,  RpZ05a1,  RpZ05a2,  RpZ05a3,  RpZ05b2,  RpZ05b3 ,  &
     &     RpZ10a0,  RpZ10a1,  RpZ10a2,  RpZ10a3,  RpZ10b2,  RpZ10b3

      DATA                                                              &
     &     RpZ01a0,  RpZ01a1,  RpZ01a2,  RpZ01a3,  RpZ01b2,  RpZ01b3 /  &
     &    -3.68063,  4.23913, -4.55894, -0.39755, -1.22774,  1.41450 /
      DATA                                                              &
     &     RpZ02a0,  RpZ02a1,  RpZ02a2,  RpZ02a3,  RpZ02b2,  RpZ02b3 /  &
     &    -4.97636,-16.09015,  0.83188,  0.21737,  6.84615, -0.98649 /
      DATA                                                              &
     &     RpZ05a0,  RpZ05a1,  RpZ05a2,  RpZ05a3,  RpZ05b2,  RpZ05b3 /  &
     &    -4.27687, -4.33629,  0.30338,  0.05697,  3.21315, -0.47749 /
      DATA                                                              &
     &     RpZ10a0,  RpZ10a1,  RpZ10a2,  RpZ10a3,  RpZ10b2,  RpZ10b3 /  &
     &    -4.94597, -1.53482,  0.10112,  0.03087,  2.45288, -0.36896 /

!
      u    = abs(uGiven)
      Mu = + MuGiven
      Z    = ZGiven
!
!     .                                 Electrons of velocity less than
!     .                                 the runaway velocity simply do not
!     .                                 run away:
      if ( u*u .LE. 1.d0 ) then
        lscsq_RunProb = 0.d0
        return
      endif
!
!     Backward-running electrons run away easily once u > 1, and it
!     does not depend much on Z.  However, there is no fit given in the
!     paper, even though a graph is given.  This is MY crude fit to
!     the Fig. 2 of the reference.
      if ( Mu .LE. -0.90d0 ) then
        x = 0.85 * (10.d0/Z)**2 * ( (u*u-1.d0)/3. )**2
        if (x .GT. 0.85d0) x = 1.0d0 - 0.15*exp(-(x-0.85d0)**2)
        lscsq_RunProb = x
        return
      endif

!     No fit given, no graph given for -1 < Mu < 1 , so I set RunProb to 0.0.
      if ( Mu .GT. -0.90d0 .and. Mu .LT. 0.90d0 ) then
        lscsq_RunProb = 0.0d0
        return
      endif

!     So we are left with u > 1, and Mu = 1. (or at least Mu > 0.90)
!     This uses the fit data from Table I., page 190.
      x = (u - 1.d0)
      Y = x
      if        ( Z .LE. 1.5d0 ) then
        lscsq_RunProb = exp (                                             &
     &  ( RpZ01a0 + RpZ01a1*Y + RpZ01a2*Y**2 + RpZ01a3*Y**3 ) /         &
     &  (               1.0*Y + RpZ01b2*Y**2 + RpZ01b3*Y**3 )           &
     &                )
        else if ( Z .LE. 3.0d0 ) then
        lscsq_RunProb = exp (                                             &
     &  ( RpZ02a0 + RpZ02a1*Y + RpZ02a2*Y**2 + RpZ02a3*Y**3 ) /         &
     &  (               1.0*Y + RpZ02b2*Y**2 + RpZ02b3*Y**3 )           &
     &                )
        else if ( Z .LE. 7.0d0 ) then
        lscsq_RunProb = exp (                                             &
     &  ( RpZ05a0 + RpZ05a1*Y + RpZ05a2*Y**2 + RpZ05a3*Y**3 ) /         &
     &  (               1.0*Y + RpZ05b2*Y**2 + RpZ05b3*Y**3 )           &
     &                )
        else if ( Z .LE. 10.d0 ) then
        lscsq_RunProb = exp (                                             &
     &  ( RpZ10a0 + RpZ10a1*Y + RpZ10a2*Y**2 + RpZ10a3*Y**3 ) /         &
     &  (               1.0*Y + RpZ10b2*Y**2 + RpZ10b3*Y**3 )           &
     &                )
      endif
      return
      END
!
!     -----------------------------------------------------------------|-------
!
      FUNCTION lscsq_WsloDwn ( uGiven , MuGiven , ZGiven )
use iso_c_binding, only : fp => c_double
implicit none
!     WsloDwn   is the inergy in units of m v_runaway^2 imparted to the
!               electric field by an electron as it slows down.
      REAL(fp) lscsq_WsloDwn
      REAL(fp)                                                              &
     &        u , Mu , Z , uGiven, MuGiven, ZGiven,                     &
     &        x
      REAL(fp)                                                              &
     &      Mu0 , Mu1 , Mu2 , Mu3 , Mu4 , U02, U04, U06, U08, U10

      real(fp)    Y
      real(fp)                                                              &
     &     WpZ01a2,  WpZ01a3,  WpZ01a4,  WpZ01b1,  WpZ01b2,  WpZ01b3 ,  &
     &     WpZ02a2,  WpZ02a3,  WpZ02a4,  WpZ02b1,  WpZ02b2,  WpZ02b3 ,  &
     &     WpZ05a2,  WpZ05a3,  WpZ05a4,  WpZ05b1,  WpZ05b2,  WpZ05b3 ,  &
     &     WpZ10a2,  WpZ10a3,  WpZ10a4,  WpZ10b1,  WpZ10b2,  WpZ10b3
      DATA                                                              &
     &     WpZ01a2,  WpZ01a3,  WpZ01a4,  WpZ01b1,  WpZ01b2,  WpZ01b3 /  &
     &     0.16612, -0.01495,  0.00775,  0.37136,  0.02240,  0.01645 /
      DATA                                                              &
     &     WpZ02a2,  WpZ02a3,  WpZ02a4,  WpZ02b1,  WpZ02b2,  WpZ02b3 /  &
     &     0.14200, -0.04048,  0.01145,  0.12253,  0.00384,  0.02440 /
      DATA                                                              &
     &     WpZ05a2,  WpZ05a3,  WpZ05a4,  WpZ05b1,  WpZ05b2,  WpZ05b3 /  &
     &     0.09880, -0.05152,  0.01113, -0.19484,  0.00559,  0.02362 /
      DATA                                                              &
     &     WpZ10a2,  WpZ10a3,  WpZ10a4,  WpZ10b1,  WpZ10b2,  WpZ10b3 /  &
     &     0.06537, -0.03895,  0.00738, -0.32456,  0.02797,  0.01526 /


      real(fp)                                                              &
     &     WmZ01a2,  WmZ01a3,  WmZ01a4,  WmZ01a5 ,                      &
     &     WmZ02a2,  WmZ02a3,  WmZ02a4,  WmZ02a5 ,                      &
     &     WmZ05a2,  WmZ05a3,  WmZ05a4,  WmZ05a5 ,                      &
     &     WmZ10a2,  WmZ10a3,  WmZ10a4,  WmZ10a5
      DATA                                                              &
     &     WmZ01a2,  WmZ01a3,  WmZ01a4,  WmZ01a5 /                      &
     &    -0.16483, -0.13420,  0.15346, -0.24314 /
      DATA                                                              &
     &     WmZ02a2,  WmZ02a3,  WmZ02a4,  WmZ02a5 /                      &
     &    -0.14186, -0.09297,  0.06661, -0.12870 /
      DATA                                                              &
     &     WmZ05a2,  WmZ05a3,  WmZ05a4,  WmZ05a5 /                      &
     &    -0.09975, -0.04781,  0.00606, -0.03545 /
      DATA                                                              &
     &     WmZ10a2,  WmZ10a3,  WmZ10a4,  WmZ10a5 /                      &
     &    -0.06651, -0.02797, -0.00247, -0.00934 /
!
!
      u    = abs(uGiven)
      Mu = + MuGiven
      Z    = ZGiven
!
!
!     If Mu < 0.5, then use the MACSYMA-derived formula of page 191.
      U02 = u*u
      if ( U02 .LE. 0.25d0 ) then

        Mu0 = 1.0d0
        Mu1 = Mu
        Mu2 = Mu*Mu
        Mu3 = Mu*Mu2
        Mu4 = Mu2*Mu2

        U04 = U02*U02
        U06 = U02*U04
        U08 = U04*U04
        U10 = U04*U06
        Y   = Z
!       --------------------------------------------------------------|
        lscsq_WsloDwn = Mu1 * U04 / (5.+Y)                                &
     &   -    (2. + Y + 3.*Mu2)                                * U06 /  &
     &                                     ( 3.*(3.+Y)*(5.+Y) )         &
     &   + 2.*( (24.+19.*Y+3.*Y**2)*Mu1 + (9.+Y)*Mu3 )         * U08 /  &
     &                       ( (3.+Y)*(5.+Y)*(7.+3.*Y)*(9.+Y) )         &
     &  -(                                                              &
     &    (1041.+1864.*Y + 1189.*Y**2 + 316.*Y**3 + 30.*Y**4 ) * Mu0    &
     &   +( 417.+ 497.*Y +  181.*Y**2 +  21.*Y**3) * 10.       * Mu2    &
     &   +  (9.+Y)*(13.+3.*Y)                      *  5.       * Mu4    &
     &   )                                                     * U10 /  &
     &      ( 5.*(2.+Y)*(3.+Y)*(5.+Y)*(7.+3.*Y)*(9.+Y)*(13.+3.*Y) )

!       --------------------------------------------------------------|

        return
      endif

!     So if u > 0.5 but -1 < Mu < 1 , then set WsloDwn to 0.00
!     The paper gives no fits for this region, but Fig. 4 does graph it.
      if ( Mu .GT. -0.90d0 .and. Mu .LT. 0.90d0 ) then
        lscsq_WsloDwn = 0.00d0
        return
      endif


      x = u*u
      Y = x
!     Mu is +1, so use the Table II. on page 190.
!     This is valid for u < 5.
      if ( Mu .GE. 0.90d0 .and. u .LE. 5.0d0 ) then
             if ( Z  .LE. 1.5d0  )  then
                lscsq_WsloDwn =                                           &
     &                 ( WpZ01a2*Y**2 + WpZ01a3*Y**3 + Wpz01a4*Y**4 ) / &
     &            ( 1. + WpZ01b1*Y    + WpZ01b2*Y**2 + WpZ01b3*Y**3 )
        else if ( Z  .LE. 3.0d0  )  then
                lscsq_WsloDwn =                                           &
     &                 ( WpZ02a2*Y**2 + WpZ02a3*Y**3 + Wpz02a4*Y**4 ) / &
     &            ( 1. + WpZ02b1*Y    + WpZ02b2*Y**2 + WpZ02b3*Y**3 )

        else if ( Z  .LE. 7.0d0  )  then
                lscsq_WsloDwn =                                           &
     &                 ( WpZ05a2*Y**2 + WpZ05a3*Y**3 + Wpz05a4*Y**4 ) / &
     &            ( 1. + WpZ05b1*Y    + WpZ05b2*Y**2 + WpZ05b3*Y**3 )

        else
                lscsq_WsloDwn =                                           &
     &                 ( WpZ10a2*Y**2 + WpZ10a3*Y**3 + Wpz10a4*Y**4 ) / &
     &            ( 1. + WpZ10b1*Y    + WpZ10b2*Y**2 + WpZ10b3*Y**3 )
        endif
      return
      endif

!     If Mu = -1 and 0 < u < 1 then use the fit of Table III. on page 191.
      if ( Mu .LE.-0.90d0 .and. u*u .LE. 1.0d0 ) then
             if ( Z  .LE. 1.5d0  )  then
                lscsq_WsloDwn = WmZ01a2*Y**2                              &
     &                 + WmZ01a3*Y**3 + WmZ01a4*Y**4 + WmZ01a5*Y**5
        else if ( Z  .LE. 3.0d0  )  then
                lscsq_WsloDwn = WmZ02a2*Y**2                              &
     &                 + WmZ02a3*Y**3 + WmZ02a4*Y**4 + WmZ02a5*Y**5
        else if ( Z  .LE. 7.0d0  )  then
                lscsq_WsloDwn = WmZ05a2*Y**2                              &
     &                 + WmZ05a3*Y**3 + WmZ05a4*Y**4 + WmZ05a5*Y**5
        else
                lscsq_WsloDwn = WmZ10a2*Y**2                              &
     &                 + WmZ10a3*Y**3 + WmZ10a4*Y**4 + WmZ10a5*Y**5
        endif
      return
      endif

!     Unanticipated input parameter, no data; set return to 0.00.
      lscsq_WsloDwn = 0.0d0
      return
      END

!
!     -----------------------------------------------------------------|------
!

      SUBROUTINE lscsq_WsloPrm ( uGiven , MuGiven , ZGiven,               &
     &                     dWsduou, iWhichWay)

!     Provides for a linear interpolation of Zeff in the region 1.0 <Zeff< 10.
!     The interpolation is calculated with nodes at Z = 1, 2, 5, 10 as given
!     in the Karney and Fisch paper.

!     Written by D. Enright, June 1992.

use iso_c_binding, only : fp => c_double
implicit none
      INTEGER iWhichWay
      REAL(fp)                                                              &
     &        Z , uGiven, MuGiven, ZGiven,                              &
     &        dWsduou, dWsduou1, dWsduou2, dWsduou5, dWsduou10
      REAL(fp)       ONE   ,       TWO   ,       FIVE  ,       TEN
      DATA       ONE   ,       TWO   ,       FIVE  ,       TEN    /     &
     &     1.0d0, 2.0d0, 5.0d0, 10.0d0/
!
      Z    = ZGiven
!
!
      if (Z .LE. ONE) then
         CALL lscsq_WsloPrmZ(uGiven, MuGiven, ONE, dWsduou, iWhichWay)
      else if (Z .LE. TWO) then
         CALL lscsq_WsloPrmZ(uGiven, MuGiven, ONE, dWsduou1, iWhichWay)
         CALL lscsq_WsloPrmZ(uGiven, MuGiven, TWO, dWsduou2, iWhichWay)
         dWsduou = (dWsduou2 - dWsduou1)/(TWO-ONE)*(Z-ONE) + dWsduou1
      else if (Z .LE. FIVE) then
         CALL lscsq_WsloPrmZ(uGiven, MuGiven, TWO, dWsduou2, iWhichWay)
         CALL lscsq_WsloPrmZ(uGiven, MuGiven,FIVE, dWsduou5, iWhichWay)
         dWsduou = (dWsduou5 - dWsduou2)/(FIVE-TWO)*(Z-TWO)+ dWsduou2
      else if (Z .LE. TEN) then
         CALL lscsq_WsloPrmZ(uGiven, MuGiven,FIVE, dWsduou5, iWhichWay)
         CALL lscsq_WsloPrmZ(uGiven, MuGiven, TEN, dWsduou10, iWhichWay)
         dWsduou = (dWsduou10 - dWsduou5)/(TEN-FIVE)*(Z-FIVE)+ dWsduou5
      else
         CALL lscsq_WsloPrmZ(uGiven, MuGiven, TEN , dWsduou, iWhichWay)
      endif

      END

!
!     -----------------------------------------------------------------|-------
!
      SUBROUTINE lscsq_WsloPrmZ( uGiven , MuGiven , ZGiven,               &
     &                     dWsduou, iWhichWay)
!     Ratio of power coupled form the rf source into electromagnetic energy to
!     the rf power absorbed by electrons; P_el / P_in.
use iso_c_binding, only : fp => c_double
implicit none
      INTEGER iWhichWay
      REAL(fp)                                                              &
     &        u , Mu , Z , uGiven, MuGiven, ZGiven,                     &
     &        x, dWsduou
      real(fp)    Y
      REAL(fp)    e01, e02, e03, e04, e1, e2, A0,  A1,  A2
      real(fp)                                                              &
     &    WPpZ01a1, WPpZ01a2, WPpZ01a3, WPpZ01b1, WPpZ01b2, WPpZ01b3 ,  &
     &    WPpZ02a1, WPpZ02a2, WPpZ02a3, WPpZ02b1, WPpZ02b2, WPpZ02b3 ,  &
     &    WPpZ05a1, WPpZ05a2, WPpZ05a3, WPpZ05b1, WPpZ05b2, WPpZ05b3 ,  &
     &    WPpZ10a1, WPpZ10a2, WPpZ10a3, WPpZ10b1, WPpZ10b2, WPpZ10b3
      DATA                                                              &
     &    WPpZ01a1, WPpZ01a2, WPpZ01a3, WPpZ01b1, WPpZ01b2, WPpZ01b3 /  &
     &     0.66445, -0.36032,  0.07328,  0.17769, -0.25452,  0.07278 /
      DATA                                                              &
     &    WPpZ02a1, WPpZ02a2, WPpZ02a3, WPpZ02b1, WPpZ02b2, WPpZ02b3 /  &
     &     0.56760, -0.38984,  0.08634, -0.04019, -0.24673,  0.08508 /
      DATA                                                              &
     &    WPpZ05a1, WPpZ05a2, WPpZ05a3, WPpZ05b1, WPpZ05b2, WPpZ05b3 /  &
     &     0.39906, -0.32879,  0.07670, -0.28281, -0.16275,  0.07436 /
      DATA                                                              &
     &    WPpZ10a1, WPpZ10a2, WPpZ10a3, WPpZ10b1, WPpZ10b2, WPpZ10b3 /  &
     &     0.27028, -0.23261,  0.05272, -0.39140, -0.07526,  0.04981 /

      real(fp)                                                              &
     &    WPmZ01a1, WPmZ01a2, WPmZ01a3, WPmZ01a4 ,                      &
     &    WPmZ02a1, WPmZ02a2, WPmZ02a3, WPmZ02a4 ,                      &
     &    WPmZ05a1, WPmZ05a2, WPmZ05a3, WPmZ05a4 ,                      &
     &    WPmZ10a1, WPmZ10a2, WPmZ10a3, WPmZ10a4

      DATA                                                              &
     &    WPmZ01a1, WPmZ01a2, WPmZ01a3, WPmZ01a4 /                      &
     &    -0.63673, -1.39960,  3.37662, -4.23684 /
      DATA                                                              &
     &    WPmZ02a1, WPmZ02a2, WPmZ02a3, WPmZ02a4 /                      &
     &    -0.55777, -0.80763,  1.43144, -2.03866 /
      DATA                                                              &
     &    WPmZ05a1, WPmZ05a2, WPmZ05a3, WPmZ05a4 /                      &
     &    -0.39704, -0.33811,  0.23607, -0.51011 /
      DATA                                                              &
     &    WPmZ10a1, WPmZ10a2, WPmZ10a3, WPmZ10a4 /                      &
     &    -0.26600, -0.17342,  0.01896, -0.13349 /
!
!
      u    = abs(uGiven)
      Mu = + MuGiven
      Z    = ZGiven
      iWhichWay = 0
!
!     .....
!
!     For Mu = 1 and 0 < u < 5 use Table IV. on page 191.
      if ( Mu .GE. 0.90d0 ) then

             if ( u .GT. 5.0d0 ) then
!     .                                 If the velocity is too large, limit it
!     .                                 to a value covered by the table, and
!     .                                 report WhichWay the velocity is large.
               u  = 5.0d0
               iWhichWay = +1
             endif
!
             x = u*u
             Y = x
!
             if ( Z  .LE. 1.5d0  )  then
                dWsduou =                                               &
     &           (   WPpZ01a1*Y    + WPpZ01a2*Y**2 + WPpZ01a3*Y**3 ) /  &
     &       ( 1. +  WPpZ01b1*Y    + WPpZ01b2*Y**2 + WPpZ01b3*Y**3 )

        else if ( Z  .LE. 3.0d0  )  then
                dWsduou =                                               &
     &           (   WPpZ02a1*Y    + WPpZ02a2*Y**2 + WPpZ02a3*Y**3 ) /  &
     &       ( 1. +  WPpZ02b1*Y    + WPpZ02b2*Y**2 + WPpZ02b3*Y**3 )

        else if ( Z  .LE. 7.0d0  )  then
                dWsduou =                                               &
     &           (   WPpZ05a1*Y    + WPpZ05a2*Y**2 + WPpZ05a3*Y**3 ) /  &
     &       ( 1. +  WPpZ05b1*Y    + WPpZ05b2*Y**2 + WPpZ05b3*Y**3 )

        else
                dWsduou =                                               &
     &           (   WPpZ10a1*Y    + WPpZ10a2*Y**2 + WPpZ10a3*Y**3 ) /  &
     &       ( 1. +  WPpZ10b1*Y    + WPpZ10b2*Y**2 + WPpZ10b3*x**3 )

        endif
      return
      endif

!     For Mu = -1 and 0 < u < 1   use Table V. on page 192.
      if ( Mu .LE.-0.90d0 .and.  u .LE. 1.0d0 ) then
!
             x = u*u
             Y = x
!
             if ( Z  .LE. 1.5d0  )  then
                dWsduou = WPmZ01a1*Y                                    &
     &                 + WPmZ01a2*Y**2 + WPmZ01a3*Y**3 + WPmZ01a4*Y**4

        else if ( Z  .LE. 3.0d0  )  then
                dWsduou = WPmZ02a1*Y                                    &
     &                 + WPmZ02a2*Y**2 + WPmZ02a3*Y**3 + WPmZ02a4*Y**4

        else if ( Z  .LE. 7.0d0  )  then
                dWsduou = WPmZ05a1*Y                                    &
     &                 + WPmZ05a2*Y**2 + WPmZ05a3*Y**3 + WPmZ05a4*Y**4

        else
                dWsduou = WPmZ10a1*Y                                    &
     &                 + WPmZ10a2*Y**2 + WPmZ10a3*Y**3 + WPmZ10a4*Y**4

        endif
      return
      endif


!     For Mu = -1 and 0 < u < 1   use Table V. on page 192.
      if ( Mu .LE.-0.90d0 .and.  u .GT. 1.0d0 ) then
!     .                                 If the velocity is too large,
!     .                                 damp the results at high u in such
!     .                                 a way that the value and derivative
!     .                                 are continuous at u=e01 (1), but the
!     .                                 values returned at u > e01 (1) are not
!     .                                 too large.  The table is not valid.
!     .                                 The hope is to induce E field
!     .                                 in TSC.
!

!     e01 etc:   power 01 etc of the u^2 (E-like) at which we expand
!     e1, e2:    first and second power of the expansion parameter u^2-e01
!     A0,01,02:  expansion coeficients around u^2=e01
             x = u*u
             iWhichWay = -1
             e01 =1.0d0
             e02=e01*e01
             e03=e01*e02
             e04=e01*e03
             e1 =(x-e01)
             e2 =(x-e01)**2
!
             if ( Z  .LE. 1.5d0  )  then
                A0 =    WPmZ01a1*e01 +    WPmZ01a2*e02                  &
     &             +    WPmZ01a3*e03 +    WPmZ01a4*e04
                A1 =    WPmZ01a1     + 2.*WPmZ01a2*e01                  &
     &             + 3.*WPmZ01a3*e02 + 4.*WPmZ01a4*e03
                A2 =                      WPmZ01a2                      &
     &             + 3.*WPmZ01a3*e01 + 6.*WPmZ01a4*e02

        else if ( Z  .LE. 3.0d0  )  then
                A0 =    WPmZ02a1*e01 +    WPmZ02a2*e02                  &
     &             +    WPmZ02a3*e03 +    WPmZ02a4*e04
                A1 =    WPmZ02a1     + 2.*WPmZ02a2*e01                  &
     &             + 3.*WPmZ02a3*e02 + 4.*WPmZ02a4*e03
                A2 =                      WPmZ02a2                      &
     &             + 3.*WPmZ02a3*e01 + 6.*WPmZ02a4*e02

        else if ( Z  .LE. 7.0d0  )  then
                A0 =    WPmZ05a1*e01 +    WPmZ05a2*e02                  &
     &             +    WPmZ05a3*e03 +    WPmZ05a4*e04
                A1 =    WPmZ05a1     + 2.*WPmZ05a2*e01                  &
     &             + 3.*WPmZ05a3*e02 + 4.*WPmZ05a4*e03
                A2 =                      WPmZ05a2                      &
     &             + 3.*WPmZ05a3*e01 + 6.*WPmZ05a4*e02

        else
                A0 =    WPmZ10a1*e01 +    WPmZ10a2*e02                  &
     &             +    WPmZ10a3*e03 +    WPmZ10a4*e04
                A1 =    WPmZ10a1     + 2.*WPmZ10a2*e01                  &
     &             + 3.*WPmZ10a3*e02 + 4.*WPmZ10a4*e03
                A2 =                      WPmZ10a2                      &
     &             + 3.*WPmZ10a3*e01 + 6.*WPmZ10a4*e02


        endif
        dWsduou = A0 + A1*e1 + A2*e2
      return
      endif


!     Unanticipated parameter range, set return to 0.00.
      dWsduou = 0.0d0
      return
      END
!                                                                      |
!     Fits.F ends       -----------------------------------------------|

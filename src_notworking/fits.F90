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
!               0 to 5_fp, with most rise as u_par goes from 0 to 3.
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
!c01  WsloDwn = Mu1 * U04 / (5.0_fp+Z)
!c02 ^   -    (2. + Z + 3.0_fp*Mu2)                                * U06 /
!c03 ^                                     ( 3.0_fp*(3.0_fp+Z)*(5.0_fp+Z) )
!c04 ^   + 2.0_fp*( (24.0_fp+19.0_fp*Z+3.0_fp*Z**2)*Mu1 + (9.0_fp+Z)*Mu3 )         * U08 /
!c05 ^                       ( (3.0_fp+Z)*(5.0_fp+Z)*(7.0_fp+3.0_fp*Z)*(9.0_fp+Z) )
!c06 ^  -(
!c07 ^    (1041.0_fp+1864.0_fp*Z + 1189.0_fp*Z**2 + 316.0_fp*Z**3 + 30.0_fp*Z**4 ) * Mu0
!c08 ^   +( 417.0_fp+ 497.0_fp*Z +  181.0_fp*Z**2 +  21.0_fp*Z**3) * 10.       * Mu2
!c09 ^   +  (9.0_fp+Z)*(13.0_fp+3.0_fp*Z)                      *  5.       * Mu4
!c10 ^   )                                                     * U10 /
!c11 ^      ( 5.0_fp*(2.0_fp+Z)*(3.0_fp+Z)*(5.0_fp+Z)*(7.0_fp+3.0_fp*Z)*(9.0_fp+Z)*(13.0_fp+3.0_fp*Z) )
!
!
!     -----------------------------------------------------------------|-------
!
     ! REAL(fp) 
FUNCTION lscsq_RunProb ( uGiven , MuGiven , ZGiven )
use iso_c_binding, only : fp => c_double
implicit none
      REAL(fp)                                                              &
     &        u , Mu , Z , uGiven, MuGiven, ZGiven,                     &
     &        x
      real(fp)    Y, lscsq_RunProb
      real(fp)                                                              &
     &     RpZ01a0,  RpZ01a1,  RpZ01a2,  RpZ01a3,  RpZ01b2,  RpZ01b3 ,  &
     &     RpZ02a0,  RpZ02a1,  RpZ02a2,  RpZ02a3,  RpZ02b2,  RpZ02b3 ,  &
     &     RpZ05a0,  RpZ05a1,  RpZ05a2,  RpZ05a3,  RpZ05b2,  RpZ05b3 ,  &
     &     RpZ10a0,  RpZ10a1,  RpZ10a2,  RpZ10a3,  RpZ10b2,  RpZ10b3

      DATA                                                              &
     &     RpZ01a0,  RpZ01a1,  RpZ01a2,  RpZ01a3,  RpZ01b2,  RpZ01b3 /  &
     &    -3.68063_fp,  4.23913_fp, -4.55894_fp, -0.39755_fp, -1.22774_fp,  1.41450_fp /
      DATA                                                              &
     &     RpZ02a0,  RpZ02a1,  RpZ02a2,  RpZ02a3,  RpZ02b2,  RpZ02b3 /  &
     &    -4.97636_fp,-16.09015_fp,  0.83188_fp,  0.21737_fp,  6.84615_fp, -0.98649_fp /
      DATA                                                              &
     &     RpZ05a0,  RpZ05a1,  RpZ05a2,  RpZ05a3,  RpZ05b2,  RpZ05b3 /  &
     &    -4.27687_fp, -4.33629_fp,  0.30338_fp,  0.05697_fp,  3.21315_fp, -0.47749_fp /
      DATA                                                              &
     &     RpZ10a0,  RpZ10a1,  RpZ10a2,  RpZ10a3,  RpZ10b2,  RpZ10b3 /  &
     &    -4.94597_fp, -1.53482_fp,  0.10112_fp,  0.03087_fp,  2.45288_fp, -0.36896_fp /

      lscsq_RunProb = 0._fp
      u    = abs(uGiven)
      Mu = + MuGiven
      Z    = ZGiven
!
!     .                                 Electrons of velocity less than
!     .                                 the runaway velocity simply do not
!     .                                 run away:
      if ( u*u .LE. 1._fp ) then
        lscsq_RunProb = 0._fp
        return
      endif
!
!     Backward-running electrons run away easily once u > 1_fp, and it
!     does not depend much on Z.  However, there is no fit given in the
!     paper, even though a graph is given.  This is MY crude fit to
!     the Fig. 2 of the reference.
      if ( Mu .LE. -0.90_fp ) then
        x = 0.85 * (10._fp/Z)**2 * ( (u*u-1._fp)/3. )**2
        if (x .GT. 0.85_fp) x = 1.0_fp - 0.15*exp(-(x-0.85_fp)**2)
        lscsq_RunProb = x
        return
      endif

!     No fit given, no graph given for -1 < Mu < 1 , so I set RunProb to 0.0.
      if ( Mu .GT. -0.90_fp .and. Mu .LT. 0.90_fp ) then
        lscsq_RunProb = 0.0_fp
        return
      endif

!     So we are left with u > 1_fp, and Mu = 1. (or at least Mu > 0.90)
!     This uses the fit data from Table I., page 190.
      x = (u - 1._fp)
      Y = x
      if        ( Z .LE. 1.5_fp ) then
        lscsq_RunProb = exp (                                             &
     &  ( RpZ01a0 + RpZ01a1*Y + RpZ01a2*Y**2 + RpZ01a3*Y**3 ) /         &
     &  (               Y + RpZ01b2*Y**2 + RpZ01b3*Y**3 )           &
     &                )
        else if ( Z .LE. 3.0_fp ) then
        lscsq_RunProb = exp (                                             &
     &  ( RpZ02a0 + RpZ02a1*Y + RpZ02a2*Y**2 + RpZ02a3*Y**3 ) /         &
     &  (               Y + RpZ02b2*Y**2 + RpZ02b3*Y**3 )           &
     &                )
        else if ( Z .LE. 7.0_fp ) then
        lscsq_RunProb = exp (                                             &
     &  ( RpZ05a0 + RpZ05a1*Y + RpZ05a2*Y**2 + RpZ05a3*Y**3 ) /         &
     &  (               Y + RpZ05b2*Y**2 + RpZ05b3*Y**3 )           &
     &                )
        else if ( Z .LE. 10._fp ) then
        lscsq_RunProb = exp (                                             &
     &  ( RpZ10a0 + RpZ10a1*Y + RpZ10a2*Y**2 + RpZ10a3*Y**3 ) /         &
     &  (               Y + RpZ10b2*Y**2 + RpZ10b3*Y**3 )           &
     &                )
      endif
      return
      END
!
!     -----------------------------------------------------------------|-------
!
      !REAL(fp) 
      FUNCTION lscsq_WsloDwn ( uGiven , MuGiven , ZGiven )
use iso_c_binding, only : fp => c_double
implicit none
!     WsloDwn   is the inergy in units of m v_runaway^2 imparted to the
!               electric field by an electron as it slows down.
      REAL(fp)                                                              &
     &        u , Mu , Z , uGiven, MuGiven, ZGiven,                     &
     &        x
      REAL(fp)                                                              &
     &      Mu0 , Mu1 , Mu2 , Mu3 , Mu4 , U02, U04, U06, U08, U10

      real(fp)    Y, lscsq_WsloDwn
      real(fp)                                                              &
     &     WpZ01a2,  WpZ01a3,  WpZ01a4,  WpZ01b1,  WpZ01b2,  WpZ01b3 ,  &
     &     WpZ02a2,  WpZ02a3,  WpZ02a4,  WpZ02b1,  WpZ02b2,  WpZ02b3 ,  &
     &     WpZ05a2,  WpZ05a3,  WpZ05a4,  WpZ05b1,  WpZ05b2,  WpZ05b3 ,  &
     &     WpZ10a2,  WpZ10a3,  WpZ10a4,  WpZ10b1,  WpZ10b2,  WpZ10b3
      DATA                                                              &
     &     WpZ01a2,  WpZ01a3,  WpZ01a4,  WpZ01b1,  WpZ01b2,  WpZ01b3 /  &
     &     0.16612_fp, -0.01495_fp,  0.00775_fp,  0.37136_fp,  0.02240_fp,  0.01645_fp /
      DATA                                                              &
     &     WpZ02a2,  WpZ02a3,  WpZ02a4,  WpZ02b1,  WpZ02b2,  WpZ02b3 /  &
     &     0.14200_fp, -0.04048_fp,  0.01145_fp,  0.12253_fp,  0.00384_fp,  0.02440_fp /
      DATA                                                              &
     &     WpZ05a2,  WpZ05a3,  WpZ05a4,  WpZ05b1,  WpZ05b2,  WpZ05b3 /  &
     &     0.09880_fp, -0.05152_fp,  0.01113_fp, -0.19484_fp,  0.00559_fp,  0.02362_fp /
      DATA                                                              &
     &     WpZ10a2,  WpZ10a3,  WpZ10a4,  WpZ10b1,  WpZ10b2,  WpZ10b3 /  &
     &     0.06537_fp, -0.03895_fp,  0.00738_fp, -0.32456_fp,  0.02797_fp,  0.01526_fp /


      real(fp)                                                              &
     &     WmZ01a2,  WmZ01a3,  WmZ01a4,  WmZ01a5 ,                      &
     &     WmZ02a2,  WmZ02a3,  WmZ02a4,  WmZ02a5 ,                      &
     &     WmZ05a2,  WmZ05a3,  WmZ05a4,  WmZ05a5 ,                      &
     &     WmZ10a2,  WmZ10a3,  WmZ10a4,  WmZ10a5
      DATA                                                              &
     &     WmZ01a2,  WmZ01a3,  WmZ01a4,  WmZ01a5 /                      &
     &    -0.16483_fp, -0.13420_fp,  0.15346_fp, -0.24314_fp /
      DATA                                                              &
     &     WmZ02a2,  WmZ02a3,  WmZ02a4,  WmZ02a5 /                      &
     &    -0.14186_fp, -0.09297_fp,  0.06661_fp, -0.12870_fp /
      DATA                                                              &
     &     WmZ05a2,  WmZ05a3,  WmZ05a4,  WmZ05a5 /                      &
     &    -0.09975_fp, -0.04781_fp,  0.00606_fp, -0.03545_fp /
      DATA                                                              &
     &     WmZ10a2,  WmZ10a3,  WmZ10a4,  WmZ10a5 /                      &
     &    -0.06651_fp, -0.02797_fp, -0.00247_fp, -0.00934_fp /
!
!
      u    = abs(uGiven)
      Mu = + MuGiven
      Z    = ZGiven
!
!
!     If Mu < 0.5_fp, then use the MACSYMA-derived formula of page 191.
      U02 = u*u
      if ( U02 .LE. 0.25_fp ) then

        Mu0 = 1.0_fp
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
        lscsq_WsloDwn = Mu1 * U04 / (5.0_fp+Y)                                &
     &   -    (2._fp + Y + 3.0_fp*Mu2)                                * U06 /  &
     &                                     ( 3.0_fp*(3.0_fp+Y)*(5.0_fp+Y) )         &
     &   + 2.0_fp*( (24.0_fp+19.0_fp*Y+3.0_fp*Y**2)*Mu1 + (9.0_fp+Y)*Mu3 )         * U08 /  &
     &                       ( (3.0_fp+Y)*(5.0_fp+Y)*(7.0_fp+3.0_fp*Y)*(9.0_fp+Y) )         &
     &  -(                                                              &
     &    (1041.0_fp+1864.0_fp*Y + 1189.0_fp*Y**2 + 316.0_fp*Y**3 + 30.0_fp*Y**4 ) * Mu0    &
     &   +( 417.0_fp+ 497.0_fp*Y +  181.0_fp*Y**2 +  21.0_fp*Y**3) * 10.       * Mu2    &
     &   +  (9.0_fp+Y)*(13.0_fp+3.0_fp*Y)                      *  5.       * Mu4    &
     &   )                                                     * U10 /  &
     &      ( 5.0_fp*(2.0_fp+Y)*(3.0_fp+Y)*(5.0_fp+Y)*(7.0_fp+3.0_fp*Y)*(9.0_fp+Y)*(13.0_fp+3.0_fp*Y) )

!       --------------------------------------------------------------|

        return
      endif

!     So if u > 0.5 but -1 < Mu < 1 , then set WsloDwn to 0.00
!     The paper gives no fits for this region, but Fig. 4 does graph it.
      if ( Mu .GT. -0.90_fp .and. Mu .LT. 0.90_fp ) then
        lscsq_WsloDwn = 0.0_fp
        return
      endif


      x = u*u
      Y = x
!     Mu is +1_fp, so use the Table II. on page 190.
!     This is valid for u < 5.
      if ( Mu .GE. 0.90_fp .and. u .LE. 5.0_fp ) then
             if ( Z  .LE. 1.5_fp  )  then
                lscsq_WsloDwn =                                           &
     &                 ( WpZ01a2*Y**2 + WpZ01a3*Y**3 + Wpz01a4*Y**4 ) / &
     &            ( 1.0_fp + WpZ01b1*Y    + WpZ01b2*Y**2 + WpZ01b3*Y**3 )
        else if ( Z  .LE. 3.0_fp  )  then
                lscsq_WsloDwn =                                           &
     &                 ( WpZ02a2*Y**2 + WpZ02a3*Y**3 + Wpz02a4*Y**4 ) / &
     &            ( 1.0_fp + WpZ02b1*Y    + WpZ02b2*Y**2 + WpZ02b3*Y**3 )

        else if ( Z  .LE. 7.0_fp  )  then
                lscsq_WsloDwn =                                           &
     &                 ( WpZ05a2*Y**2 + WpZ05a3*Y**3 + Wpz05a4*Y**4 ) / &
     &            ( 1.0_fp + WpZ05b1*Y    + WpZ05b2*Y**2 + WpZ05b3*Y**3 )

        else
                lscsq_WsloDwn =                                           &
     &                 ( WpZ10a2*Y**2 + WpZ10a3*Y**3 + Wpz10a4*Y**4 ) / &
     &            ( 1.0_fp + WpZ10b1*Y    + WpZ10b2*Y**2 + WpZ10b3*Y**3 )
        endif
      return
      endif

!     If Mu = -1 and 0 < u < 1 then use the fit of Table III. on page 191.
      if ( Mu .LE.-0.9_fp .and. u*u .LE. 1.0_fp ) then
             if ( Z  .LE. 1.5_fp  )  then
                lscsq_WsloDwn = WmZ01a2*Y**2                              &
     &                 + WmZ01a3*Y**3 + WmZ01a4*Y**4 + WmZ01a5*Y**5
        else if ( Z  .LE. 3.0_fp  )  then
                lscsq_WsloDwn = WmZ02a2*Y**2                              &
     &                 + WmZ02a3*Y**3 + WmZ02a4*Y**4 + WmZ02a5*Y**5
        else if ( Z  .LE. 7.0_fp  )  then
                lscsq_WsloDwn = WmZ05a2*Y**2                              &
     &                 + WmZ05a3*Y**3 + WmZ05a4*Y**4 + WmZ05a5*Y**5
        else
                lscsq_WsloDwn = WmZ10a2*Y**2                              &
     &                 + WmZ10a3*Y**3 + WmZ10a4*Y**4 + WmZ10a5*Y**5
        endif
      return
      endif

!     Unanticipated input parameter, no data; set return to 0.00.
      lscsq_WsloDwn = 0.0_fp
      return
      END

!
!     -----------------------------------------------------------------|------
!

      SUBROUTINE lscsq_WsloPrm ( uGiven , MuGiven , ZGiven,               &
     &                     dWsduou, iWhichWay)

!     Provides for a linear interpolation of Zeff in the region 1.0 <Zeff< 10.
!     The interpolation is calculated with nodes at Z = 1_fp, 2_fp, 5_fp, 10 as given
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
     &     1.0_fp, 2.0_fp, 5.0_fp, 10.0_fp/
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
      REAL(fp)    e01, e02, e03, e04, e1, e2, a0,  a1,  A2
      real(fp)                                                              &
     &    WPpZ01a1, WPpZ01a2, WPpZ01a3, WPpZ01b1, WPpZ01b2, WPpZ01b3 ,  &
     &    WPpZ02a1, WPpZ02a2, WPpZ02a3, WPpZ02b1, WPpZ02b2, WPpZ02b3 ,  &
     &    WPpZ05a1, WPpZ05a2, WPpZ05a3, WPpZ05b1, WPpZ05b2, WPpZ05b3 ,  &
     &    WPpZ10a1, WPpZ10a2, WPpZ10a3, WPpZ10b1, WPpZ10b2, WPpZ10b3
      DATA                                                              &
     &    WPpZ01a1, WPpZ01a2, WPpZ01a3, WPpZ01b1, WPpZ01b2, WPpZ01b3 /  &
     &     0.66445_fp, -0.36032_fp,  0.07328_fp,  0.17769_fp, -0.25452_fp,  0.07278_fp /
      DATA                                                              &
     &    WPpZ02a1, WPpZ02a2, WPpZ02a3, WPpZ02b1, WPpZ02b2, WPpZ02b3 /  &
     &     0.56760_fp, -0.38984_fp,  0.08634_fp, -0.04019_fp, -0.24673_fp,  0.08508_fp /
      DATA                                                              &
     &    WPpZ05a1, WPpZ05a2, WPpZ05a3, WPpZ05b1, WPpZ05b2, WPpZ05b3 /  &
     &     0.39906_fp, -0.32879_fp,  0.07670_fp, -0.28281_fp, -0.16275_fp,  0.07436_fp /
      DATA                                                              &
     &    WPpZ10a1, WPpZ10a2, WPpZ10a3, WPpZ10b1, WPpZ10b2, WPpZ10b3 /  &
     &     0.27028_fp, -0.23261_fp,  0.05272_fp, -0.39140_fp, -0.07526_fp,  0.04981_fp /

      real(fp)                                                              &
     &    WPmZ01a1, WPmZ01a2, WPmZ01a3, WPmZ01a4 ,                      &
     &    WPmZ02a1, WPmZ02a2, WPmZ02a3, WPmZ02a4 ,                      &
     &    WPmZ05a1, WPmZ05a2, WPmZ05a3, WPmZ05a4 ,                      &
     &    WPmZ10a1, WPmZ10a2, WPmZ10a3, WPmZ10a4

      DATA                                                              &
     &    WPmZ01a1, WPmZ01a2, WPmZ01a3, WPmZ01a4 /                      &
     &    -0.63673_fp, -1.39960_fp,  3.37662_fp, -4.23684_fp /
      DATA                                                              &
     &    WPmZ02a1, WPmZ02a2, WPmZ02a3, WPmZ02a4 /                      &
     &    -0.55777_fp, -0.80763_fp,  1.43144_fp, -2.03866_fp /
      DATA                                                              &
     &    WPmZ05a1, WPmZ05a2, WPmZ05a3, WPmZ05a4 /                      &
     &    -0.39704_fp, -0.33811_fp,  0.23607_fp, -0.51011_fp /
      DATA                                                              &
     &    WPmZ10a1, WPmZ10a2, WPmZ10a3, WPmZ10a4 /                      &
     &    -0.26600_fp, -0.17342_fp,  0.01896_fp, -0.13349_fp /
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
      if ( Mu .GE. 0.9_fp ) then

             if ( u .GT. 5.0_fp ) then
!     .                                 If the velocity is too large, limit it
!     .                                 to a value covered by the table, and
!     .                                 report WhichWay the velocity is large.
               u  = 5.0_fp
               iWhichWay = +1
             endif
!
             x = u*u
             Y = x
!
             if ( Z  .LE. 1.5_fp  )  then
                dWsduou =                                               &
     &           (   WPpZ01a1*Y    + WPpZ01a2*Y**2 + WPpZ01a3*Y**3 ) /  &
     &       ( 1. +  WPpZ01b1*Y    + WPpZ01b2*Y**2 + WPpZ01b3*Y**3 )

        else if ( Z  .LE. 3.0_fp  )  then
                dWsduou =                                               &
     &           (   WPpZ02a1*Y    + WPpZ02a2*Y**2 + WPpZ02a3*Y**3 ) /  &
     &       ( 1. +  WPpZ02b1*Y    + WPpZ02b2*Y**2 + WPpZ02b3*Y**3 )

        else if ( Z  .LE. 7.0_fp  )  then
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
      if ( Mu .LE.-0.9_fp .and.  u .LE. 1.0_fp ) then
!
             x = u*u
             Y = x
!
             if ( Z  .LE. 1.5_fp  )  then
                dWsduou = WPmZ01a1*Y                                    &
     &                 + WPmZ01a2*Y**2 + WPmZ01a3*Y**3 + WPmZ01a4*Y**4

        else if ( Z  .LE. 3.0_fp  )  then
                dWsduou = WPmZ02a1*Y                                    &
     &                 + WPmZ02a2*Y**2 + WPmZ02a3*Y**3 + WPmZ02a4*Y**4

        else if ( Z  .LE. 7.0_fp  )  then
                dWsduou = WPmZ05a1*Y                                    &
     &                 + WPmZ05a2*Y**2 + WPmZ05a3*Y**3 + WPmZ05a4*Y**4

        else
                dWsduou = WPmZ10a1*Y                                    &
     &                 + WPmZ10a2*Y**2 + WPmZ10a3*Y**3 + WPmZ10a4*Y**4

        endif
      return
      endif


!     For Mu = -1 and 0 < u < 1   use Table V. on page 192.
      if ( Mu .LE.-0.9_fp .and.  u .GT. 1.0_fp ) then
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
!     e1_fp, e2:    first and second power of the expansion parameter u^2-e01
!     a0,01_fp,02:  expansion coeficients around u^2=e01
             x = u*u
             iWhichWay = -1
             e01 =1.0_fp
             e02=e01*e01
             e03=e01*e02
             e04=e01*e03
             e1 =(x-e01)
             e2 =(x-e01)**2
!
             if ( Z  .LE. 1.5_fp  )  then
                A0 =    WPmZ01a1*e01 +    WPmZ01a2*e02                  &
     &             +    WPmZ01a3*e03 +    WPmZ01a4*e04
                A1 =    WPmZ01a1     + 2.0_fp*WPmZ01a2*e01                  &
     &             + 3.0_fp*WPmZ01a3*e02 + 4.0_fp*WPmZ01a4*e03
                A2 =                      WPmZ01a2                      &
     &             + 3.0_fp*WPmZ01a3*e01 + 6.0_fp*WPmZ01a4*e02

        else if ( Z  .LE. 3.0_fp  )  then
                A0 =    WPmZ02a1*e01 +    WPmZ02a2*e02                  &
     &             +    WPmZ02a3*e03 +    WPmZ02a4*e04
                A1 =    WPmZ02a1     + 2.0_fp*WPmZ02a2*e01                  &
     &             + 3.0_fp*WPmZ02a3*e02 + 4.0_fp*WPmZ02a4*e03
                A2 =                      WPmZ02a2                      &
     &             + 3.0_fp*WPmZ02a3*e01 + 6.0_fp*WPmZ02a4*e02

        else if ( Z  .LE. 7.0_fp  )  then
                A0 =    WPmZ05a1*e01 +    WPmZ05a2*e02                  &
     &             +    WPmZ05a3*e03 +    WPmZ05a4*e04
                A1 =    WPmZ05a1     + 2.0_fp*WPmZ05a2*e01                  &
     &             + 3.0_fp*WPmZ05a3*e02 + 4.0_fp*WPmZ05a4*e03
                A2 =                      WPmZ05a2                      &
     &             + 3.0_fp*WPmZ05a3*e01 + 6.0_fp*WPmZ05a4*e02

        else
                A0 =    WPmZ10a1*e01 +    WPmZ10a2*e02                  &
     &             +    WPmZ10a3*e03 +    WPmZ10a4*e04
                A1 =    WPmZ10a1     + 2.0_fp*WPmZ10a2*e01                  &
     &             + 3.0_fp*WPmZ10a3*e02 + 4.0_fp*WPmZ10a4*e03
                A2 =                      WPmZ10a2                      &
     &             + 3.0_fp*WPmZ10a3*e01 + 6.0_fp*WPmZ10a4*e02


        endif
        dWsduou = A0 + A1*e1 + A2*e2
      return
      endif


!     Unanticipated parameter range, set return to 0.00.
      dWsduou = 0.0_fp
      return
      END
!                                                                      |
!     Fits.F ends       -----------------------------------------------|

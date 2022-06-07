!     Copyright 1991 by E. J. Valeo, C. F. F. Karney and N. J. Fisch,
!     Princeton University Plasma Physics Laboratory.
!     Ref: Current in wave-driven plasmas, by CFFK and NJF, Physics of
!     Fluids <29> pp 180-192 (1986).  See esp. eqn. 21b.
!     In constructing units, the rest of the code has velocities normalized
!     to c, and f contains a density in cm^{-3}.
!     Dql has a c^2 taken out of it.
!     E is in v/m we assume.
!     \Gamma   = n_e q^4     \ln \lambda / (4 \pi \epsilon_0^2 m^2)
!              = n_e q^4 c^4 \ln \lambda \mu_0/(4\pi) \mu_0 / m^2
!                 with units of velocity^3 / time
!     v_r      = -sign(qE) \sqrt{m\Gamma/\abs(qE)}
!                 where q has the sign of the electron charge = - 1.6e-19
!                 and the electrons run away in the negative direction if
!                 v_r is positive
!     v_{norm} = \abs(v_r) / c
!              = sqrt ( \Gamma / \abs(qE/m) ) / c
!              = sqrt ( n/E q^3/m c^2 lnLambda mu_0 (mu_0/4PI) )
!                 dimensionless, in units of c
!     gmrun    = n_e q^4 ln(Lambda) / ( 4 PI eps_0^2 m_e^2 )          /  c^3
!              = n_e q^4 c^1 / m_e^2 * lnLambda * (mu_0/(4PI))*mu_0
!                 with units of time^{-1}
!
!     nuRuna   = gmrun/vnorm^3
!                 with units of time^{-1}
!     NeAry in cm^{-3}
!     gmrun in s^{-1}
!     nuRuna in s^{-1}
!     js in ampere/m^2
!     dfdv in cm^{-3}
!     dlnJdlnE = dJ/dE * (E/J) is dimensionless
!     and can be shown for small E (E is the E_{dc}) to be approximately
!     eE_{dc}/(m Gamma) * v_phase^2 * (1/(2 \mu)) ( (2+Z+3\mu^2)/(3+Z) )
!     =~ eE_{dc} / lnLamba /n_\parallel^2 /
!                  [ e^2/(4\pi\epsilon_0 c^2/\omega_p^2)] * Order(1)
!
!     ------------------------------------------------------------------
!
!     Keeping the signs straight is a problem.  We say that positive
!     phase velocity tends to drive current in the supportive, or correct
!     direction.  Positive current drive.  Thus, presumably, the Edc
!     is also in the positive direction.  Electrons tend to run away
!     and be accelerated in the opposite direction.  Thus for positive
!     electron velocity which goes with positive phase velocity, the
!     normalized   u \equiv v/v_r is negative, and the mu for this positive
!     v is negative.
!
!     v  E      u    mu
!     +  +      -     -
!     -  +      +     +
!     +  -      +     +
!     -  -      -     -
!
!     Code is written with the presumtion of E positive and mu minus
!     for positive velocity and mu plus for negative velocity.  If
!     E is negative, the quantitiy muminus becomes +1; muplus becomes -1
!
!     ------------------------------------------------------------------
!
subroutine lscsq_VmaxCalc
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : ivind,izind,Vparmaxn,Vparmaxp, vpar, npsi
  use lscsq_mod, only : nzones, nrays, ok_ray
  implicit none

  integer :: ip,ir,iz,iv
  real(fp):: vl

  VparMaxN(1:npsi)= -0.05_fp
  VparMaxP(1:npsi)= +0.05_fp
 
  do ir=1,nrays
     if (ok_ray(ir).eq.1) then
        do iz=1,nzones
           iv = ivind(iz,ir)
           ip = izind(iz,ir)
           if (iv .EQ. 0 .or. ip .EQ. 0) exit
           vl = vpar(iv)
           if(VparMaxP(ip) .LT. vl) VparMaxP(ip) = vl
           if(VparMaxN(ip) .GT. vl) VparMaxN(ip) = vl
        enddo     
     endif
  enddo         

end subroutine lscsq_VmaxCalc
!
!     ------------------------------------------------------------------
!
SUBROUTINE lscsq_JdepCalc

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: jray, edcary,vnormPos,VparMaxP,vnormNeg,vnorm
  use lscsq_mod, only: vparmaxN,vrunidx,nrundot,jrundot,js
  use lscsq_mod, only: nv, npsi
  implicit none

  integer :: ip, ipLastRuna, iGotRuna, nGotRuna, iFillJray
  integer :: nrLines = 0
  character(len=70) :: ErrMsg
  real(fp), dimension(npsi) :: vnormnok

  ! compute rf driven current given power dissipation vs vpar, psi
  call lscsq_GetEdc(0.0_fp)
  ! Fill array EdcAry

  jray(1:nv,1:npsi) = 0.0
  CALL lscsq_VmaxCalc
  nGotRuna = 0
  vnormnOK(1:npsi) = 0.0

  do ip = 1, npsi
     CALL lscsq_jnorm(ip)
     !    set up normalization, tabulate u lookup values
     !    This is based on figure 2 of the Karney Fisch paper.  It shows that
     !    if you are trying to drive current in cooperation with the electric field
     !    at v/vrun you have had it.  If trying to drive against the electric field
     !    you can go up to v/vrun of 2 before runaways take over.
     if (EdcAry(ip) .GE. 0.0 ) then
        vnormPos(ip) = VparMaxP(ip)/vnorm 
        vnormNeg(ip) = 0.5_fp*VparMaxN(ip)/vnorm 
     else
        vnormPos(ip) = 0.5_fp*VparMaxP(ip)/vnorm 
        vnormNeg(ip) = VparMaxN(ip)/vnorm
     endif

     vnormNOK(ip) = max(int(vnormPos(ip)),abs(int(vnormNeg(ip))))

     ! compute jrf
     iFillJray = 1
     CALL lscsq_mkj(ip,js(ip), iGotRuna, iFillJray)
     ! Report if a runaway problem
     if (iGotRuna .GE. 1) then
        nGotRuna = nGotRuna + 1
        ipLastRuna = ip
     endif
     ! compute d/dt {n_runa; J-runa} as in eqn (20) and (21c) of
     ! Karney Fisch paper.  Aug 93
     CALL lscsq_ddtNrnJrn ( ip,nRunDot(ip), jRunDot(ip), vRunIdx(ip) )
  enddo
 
  if (nGotRuna .GE. 1) then
     write(ErrMsg,'(1x,i4,'' Rway shells; last at ip='',i3)') nGotRuna, ipLastRuna
     call lscsq_LSCwarn(ErrMsg)
  endif

end subroutine lscsq_Jdepcalc
!
!     ------------------------------------------------------------------
!
subroutine lscsq_JsplCalc
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : dEdcAmnt, jsp, npsi
  implicit none

  integer :: iFillJray=0
  integer :: ip, iGotRuna
  real(fp) :: jd  ! current driven with the E field incremented by dEdcAmnt

  !     compute rf driven current given power dissipation vs vpar, psi
  call lscsq_GetEdc(dEdcAmnt)
  ! Fill array EdcAry with dEdcAmnt over the actual Edc supplied from TSC
  do ip = 1, npsi
     CALL lscsq_jnorm(ip)
     ! set up normalization, tabulate u lookup values; These depend on Ez
     CALL lscsq_mkj(ip,jd, iGotRuna,iFillJray )
     ! compute jrf driven, jd
     jsp(ip) = jd
   enddo

   CALL lscsq_GetEdc(0.0_fp)

end subroutine lscsq_JsplCalc
!     ------------------------------------------------------------------
!
subroutine lscsq_jnorm(ip)

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : pi, twopi, vc, qe_eV, me_g
  use lscsq_mod, only: vnorm,muplus,muminus,gmrun,ugr,LnlAry,ugr
  use lscsq_mod, only: neary, Edcary, ivrun, vpar, vnmax, nuruna, nv
  implicit none

  integer, intent(in) :: ip ! flux surface index
  integer :: iv
  
  !     set up normalization, tabulate u lookup values
  gmrun = 1.0e-2_fp*NeAry(ip)*LnlAry(ip)*qe_eV**4/me_g**2*vc*2.0_fp*twopi 
  if(EdcAry(ip) .ne. 0.0) then
     vnorm = 1.0e-9_fp*sqrt(abs(NeAry(ip)/EdcAry(ip)))*    &
             sqrt(1e13_fp*qe_eV**3/me_g*vc**2 *2.0_fp*twopi*LnlAry(ip))
     if( EdcAry(ip) .GT. 0.) then
        muplus  = +1.0_fp
        muminus = -1.0_fp
     else
        muplus  = -1.0_fp
        muminus = +1.0_fp
     endif
  else
     vnorm = vnmax
     muplus = 1.0_fp
     muminus = -1.0_fp
  endif
  nuRuna = gmrun/vnorm**3
  ugr(1:nv) = Vpar(1:nv)/vnorm *muminus
 
  ivrun = 0
  if (muminus .EQ. -1.0_fp) then
     do iv = (nv+1)/2, nv-1
        if ( ugr(iv)  .LE. -1.0_fp .and. ugr(iv-1)  .GT. -1.0_fp) then
           ivrun = iv
           return
         endif
     enddo
  else
     do iv = 2,(nv+1)/2
        if ( ugr(iv)  .LE. 1.0_fp .and. ugr(iv-1)  .GT. 1.0_fp) then
           ivrun = iv
           return
        endif
     enddo
  endif

end subroutine lscsq_jnorm
!
!     ------------------------------------------------------------------
!
subroutine lscsq_mkj(ip, jd, iGotRuna, iFillJray)
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : pi, twopi, qe_eV, vc
  use lscsq_mod, only: Dql,vpar,dfdv,muplus,jray,vpar,dvplus
  use lscsq_mod, only: ismo,fe,ivzero,xmag,ugr,vtherm,muminus
  use lscsq_mod, only: Zbrary,ivlary, iitr, nuruna, nv
  implicit none

  integer, intent(in) :: ip   ! index of flux surface
  integer, intent(in) :: ifilljray
  integer, intent(out):: igotruna
  real(fp), intent(out) :: jd

  integer :: iv, iWhichWay, iSMOi
  real(fp) :: jdtemp, constfac, xr, eps, w
  real(fp) :: FeCutOff, TrapFac
  real(fp) :: FeCutFac=1.0e-5_fp
  real(fp) :: dWsdu, dWsduou
  real(fp) :: eps_p

!     compute jrf driven, jd, at flux surface ip
  jd = 0.0
  constfac = 1.0e6_fp*qe_eV*vc/nuruna   
  iGotRuna = 0
  iSMOi = mod(iSMO,2) + 1
  FeCutOff = FeCutFac*fe(IvZero,ip,iITR)
  eps = sqrt(iVlary(ip)/(2.0_fp*PI**2*xmag))/xmag
  eps_p = eps**0.77_fp  ! calculate this only once
  xr = 3.5_fp    ! fmp - what is this?

  ! first negative velocities
  do iv = 1, IvZero - 1
     ! Ignore jd where f_e is small
     if (fe(iv,ip,iITR).ge.FeCutOff ) then
        ! Muplus, because mu is plus for negative velocity if the Edc is >0
        ! See above for reversal if E<0
        CALL lscsq_WsloPrm(ugr(iv), muplus, ZbrAry(ip), dWsduou, iWhichWay)
        dWsdu = dWsduou*ugr(iv)
        ! Note that dWsdu*ugr is always positive
        if (iWhichWay.ne.0) iGotRuna = iGotRuna+1
        ! trapping effect from Ehst-Karney for LH limit of large vpar/vth
        w = abs(Vpar(iv))/Vtherm(ip)
        TrapFac = 1.0_fp-((eps_p*sqrt(xr**2+w**2))/(eps_p*xr+w))
        ! For negative velocity the current driven is negative because df/dv >0
        jdtemp = -Dql(iv,ip,iSMO)*dfdv(iv,ip,iSMOi)*TrapFac*dvplus(iv)*dWsdu*constfac
        jd = jd+jdtemp
        if(iFillJray.EQ.1) Jray(iv,ip) = jdtemp
!       if(iFillJray.EQ.1) Jray(iv,ip) = jd
     endif
  enddo
        
  ! then positive velocities
  do iv=IvZero+1, nv-1
     ! Ignore jd where f_e is small
     ! Muminus, because mu is neg for positive veloctiy if the Edc is >0
     ! See above for reversal if E<0
     if (fe(iv,ip,iITR).GE.FeCutOff) then
        CALL lscsq_WsloPrm(ugr(iv), muminus, ZbrAry(ip), dWsduou, iWhichWay)
        dWsdu = dWsduou  * ugr(iv)
        ! Note that dWsdu*ugr is always positive
        if (iWhichWay.ne.0) iGotRuna = iGotRuna+1
        ! trapping effect from Ehst-Karney for LH limit of large vpar/vth
        w = abs(Vpar(iv))/Vtherm(ip)
        TrapFac = 1.0_fp-((eps_p*sqrt(xr**2+w**2))/(eps_p*xr+w))
        ! For positive velocity the current driven is positive because df/dv<0
        jdtemp = -Dql(iv,ip,iSMO)*dfdv(iv,ip,iSMOi)*TrapFac*dvplus(iv)*dWsdu*constfac
        jd = jd + jdtemp
         if(iFillJray .EQ. 1) Jray(iv,ip) = jdtemp
!       if(iFillJray .EQ. 1) Jray(iv,ip) = jd
     endif
  enddo
 
end subroutine lscsq_mkj
 
!     ------------------------------------------------------------------
!
subroutine lscsq_GetEdc (DifAmt)
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : pi, twopi
  use lscsq_mod, only: Edcary, Edcvec, npsij, psivec, psiary, npsi
  implicit none

  real(fp), intent(in) :: DifAmt
  integer :: j
  real(fp):: psi, Edc
  real(fp) :: MinEdc = 1.0e-4_fp
  ! DifAmt is intended to be 0.0 when the E_dc is desired as given
  ! and intended to be something like 0.0001 when we are trying to form
  ! the derivative   d ln J / d ln E

  Edc = 0.0
  do j = 1, npsi
     psi = PsiAry(j)
     CALL lscsq_linr1d(NpsiJ, PsiVec, EdcVec, psi, Edc)
     Edc = Edc + DifAmt
     if (abs(Edc) .LE. MinEdc) then
        if (Edc .GE. 0.0 ) then
           Edc = +MinEdc
        else
           Edc = -MinEdc
        endif
     endif
     EdcAry(j) =  Edc
  enddo

end subroutine lscsq_getEdc
!
!     ------------------------------------------------------------------
subroutine lscsq_ddtNrnJrn (ip,nDot, jDot, vRunAwayIndex)

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : qe_eV, vc
  use lscsq_mod, only: dql, ivrun, ismo, dfdv, vpar, nv
  implicit none

  ! Added Aug 93 to help quantify the runaway situation for LSC report PPPL 2929
  integer, intent(in) :: ip      ! index of flux surface
  integer :: iSMOi
  real(fp):: nDot, jDot, vRunAwayIndex
 
  iSMOi = mod(iSMO,2) + 1

  nDot= 0.0
  JDot= 0.0
  vRunAwayIndex= 0.0
  if (ivrun .GT. 1 .and. ivrun .LT. nv ) then
     nDot =             Dql(ivrun, ip, iSMO) * dfdv(ivrun,ip,iSMOi)
     jDot = vpar(ivrun)*Dql(ivrun, ip, iSMO) * dfdv(ivrun,ip,iSMOi)
 
     nDot = abs ( nDot )
     jDot = abs(JDot)* 1.0e6_fp*qe_eV *vc
     vRunAwayIndex = REAL(ivrun,kind=fp)
  endif
 
end subroutine lscsq_ddtNrnJrn
 
!     ------------------------------------------------------------------
subroutine lscsq_SmoJandP(Radius)
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: neary, js,jsp, neary, dvol
  use lscsq_mod, only: Praysum, Praytot, LnlAry
  use lscsq_mod, only: npsi, diffujrf, prfspred
  use lscsq_mod, only : zcmtom
  use lscsq_mod, only : pi, twopi, me_g, qe_eV, vc
  implicit none

  real(fp), intent(in) :: Radius
  integer :: ip, idxFor0
  real(fp) :: fac

  real(fp), dimension(npsi)   :: Jdriven, Jdiffus
  real(fp), dimension(npsi) :: prforig, prfsmoo
  real(fp), dimension(npsi) :: nuslow, awk, bwk, cwk, jwk, xwk
  real(fp) :: PrfNorm
  real(fp), parameter :: nparcube=8.0_fp

  ! Do nothing if there is no diffusion; or no sensible radius
  if (Radius .LE. 0.0) return

!     Form the 'slowing down frequency'
!     \nu_slow = \ln \lambda n_e e^4 /(4 \pi \epsilon_0^2 m_e^2 c^3)
!                \times  n_\parallel^3
!
!     The assumed boudary conditions in SmoothJrf are that there is
!     a guard point on the low(central) side of the passed array:
!         > zero slope (in r-space) on inner side
!         > zero value on the outer side at idxFor0
!     However, idxFor0 must not be more than npsi+1

  fac = zcmtom*qe_eV**4/me_g**2*vc*4.0_fp*pi*nparcube
 
  do ip=1,npsi
     NuSlow(ip) =  fac * NeAry(ip) * LnlAry(ip) 
  enddo

  !     Smooth the current computed with the given field:
  idxFor0=0  ! DMC -- dangerous not to initialize this...
  do ip=1,npsi
     Jdriven(ip) = js(ip)
     if(Jdriven(ip) .ne. 0.0) idxFor0=ip
  enddo
  !---------------------------
  !  DMC bugfix, do nothing if J=0
  if(idxFor0.gt.0) then
     if(idxFor0 .LT. npsi-1) idxFor0 = (npsi + idxFor0)/2
     CALL lscsq_SmoothJrf(Radius, DiffuJrf, npsi, idxFor0,         &
                        NuSlow, Jdriven, Jdiffus,                    &
                        Awk, Bwk, Cwk, Jwk, Xwk,                     &
                        'PolFlux_Like', 'Zero2ndDeriv')
     js(1:npsi) = jdiffus(1:npsi)
 
     ! Smooth the current computed with the given field plus delta-field
     Jdriven(1:npsi) = jsp(1:npsi)
 
     CALL lscsq_SmoothJrf(Radius, DiffuJrf, npsi, idxFor0,         &
                        NuSlow, Jdriven, Jdiffus,                    &
                        Awk, Bwk, Cwk, Jwk, Xwk,                     &
                        'PolFlux_Like', 'Zero2ndDeriv')

     jsp(1:npsi) = Jdiffus(1:npsi)
 
     if(PrfSpred .LE. 0.0 .or. PrfSpred .GT. 1.0_fp) return
     !     Smooth Power according to
     !         Praytot(psi) = Prf-total n(psi) J-diffused(psi) dV(psi) /
     !                        [ \Sigma  n(psi) J-diffused(psi) dV(psi) ] *
     !                                                 PrfSpred     +
     !                  Prf-ray-undiffused(psi)* (1. - PrfSpred)
     !
     !     Calculate normalization factor
     PrfNorm = 0.00
     do ip=1,npsi
        PrfNorm = abs(js(ip))*NeAry(ip)*dVol(ip)+PrfNorm
     enddo
     PrfNorm = PraySum/PrfNorm
 
     !     Calculate raw smoothed power
     do ip=1,npsi
        PrfOrig(ip) = PRaytot(ip)
        PrfSmoo(ip) = abs(js(ip))*NeAry(ip)*dVol(ip) * PrfNorm
     enddo
 
     ! Fold raw smooth power with unsmoothed power
     do ip=1,npsi
        PrfSmoo(ip) = PrfSpred*PrfSmoo(ip)+(1.0_fp-PrfSpred)*PrfOrig(ip)
        PrayTot(ip) = PrfSmoo(ip)
     enddo
     !---------------------------
     !  end DMC bugfix
  endif

end subroutine lscsq_SmoJandP 
 
!     -----------------------------------------------------------------|
!
!     SUBROUTINE SmoothJrf(Radius, DiffCoef, NumJs, ....
!
!     Routine to Smooth Jrf
!     Copyright 1994 by
!     D. W. Ignat and S. C. Jardin
!     Plasma Physics Laboratory
!     Box 451
!     Princeton, New Jersey 08543
!
!     Smooth RF-driven current with a diffusion-like equation.
!     The idea is that the actual RF-driven current is modified from
!     the computed RF-driven current by slowing-down, characterized by
!     an inverse-time \nu_slow, and a cross-field diffusion,
!     characterized by DiffCoef.
!
!     Jdriven is given on points 1,2,...NumJs; regularly spaced
!     in poloidal flux. We add guard points at 0 and NumJs + 1.
!     We assume that poloidal flux is like radius squared, so that
!     the diffusion-like equation in radius
!
!     dJrf/dt = \nu_{slow} (Jdriven-Jrf) + 1/r d/dr( r DiffCoef dJrf/dr)
!
!     becomes
!
!     dJrf/dt = \nu_{slow} (Jdriven-Jrf) + d/dx(4 x DiffCoef dJrf/dx)/Radius^2
!
!     where  x  is normalized poloidal flux ranging from 0 to 1.
!
!     More correct use of generalised coordinates would follow this:
!         df/dt = 1/Jac  d/dp[ Jac (grad p)**2 D df/dp ]
!         where Jac is the Jacobian, or volume element factor.
!         Jac == (1/grad p) (\int^p dp/grad p)
!     In the case of p ~ r^2, grad p = 2 sqrt(p) and
!     \int^p dp/grad p== sqrt(p) so the J is 1/2, and cancels out.
!
!     We solve the equation for the steady state, from a difference
!     equation
!
!     Radius^2  *  \nu_slow
!     ----------------------  *  (Jrf - Jdriven) ==
!     4 DiffCoef (NumJs+1)^2
!
!    (xp+x)/2 * (Jrfp-Jrf) - (x+xm)/2 * (Jrf-Jrfm)
!
!     where xp is x-plus xm is x-minus etc and \delta x = 1/(NumJs+1)
!
!     The boundary condtions at index 0,1,2 are :
!        (1) that Jrf is constant (Jdiff_0 = Jdiff_1)              <or>
!        (2) that Jrf slope is constant (Jdiff_0 = 2 Jdiff_1 - Jdiff_2)
!
!     and at index idxFor0 that Jrf is zero.  Normally, one expects that
!     idxFor0 is NumJs+1, but maybe the current should not be zero at the edge-
!     it should be zero at idxFor0.
!
!     Boundary condition (1) goes naturally with radius-like coordinates;
!     Boundary condition (2) goes naturally with flux-like coordinates.
!
!     \nu_slow = \ln \lambda n_e e^4 /(4 \pi \epsilon_0^2 m_e^2 c^3)
!                \times  n_\parallel^3
!     where n_e is position dependent...but does not have to be.
!     Obviously, any fudge factor can be thrown in on the \nu_slow or on
!     the DiffCoef owing to the heuristic derivation of this treatment.
!     The extention to spatially varying DiffCoef is straightforward.
!
!     A modification for an input Jrf based on a radius-like grid
!     is very simple:
!
!     Radius^2  *  \nu_slow
!     ----------------------  *  r * (Jrf - Jdriven) ==
!     1 DiffCoef (NumJs+1)^2
!
!    (rp+r)/2 * (Jrfp-Jrf) - (r+rm)/2 * (Jrf-Jrfm)
!
!     where rp is r-plus rm is r-minus etc and \delta r = 1/(NumJs+1)
!
!     This spirit of this approach is similar to that of
!     V. Fuchs, I. P. Shkarofsky, R. A. Cairns, and P. T. Bonoli,
!     ``Simulations of lower hybrid current drive and ohmic transformer
!     recharge,'' Nuclear Fusion {\bf 29} 1479 1989.
!     It appears that the practical difference is that
!     in the Fuchs work a loss-rate at the edge is specified,
!     rather than the DiffCoef, and then the DiffCoef
!     is found by a shooting method.
!
!     In contrast, the present routine specifies the DiffCoef and finds
!     the modified (reduced) current.
!     (and the current loss by implication)
!
!     The following is a heuristic derivation of the formula,
!     loosely based on the Fuchs paper:                  (rectangular coords!)
!
!    df/dt = d/dv ( Dq df/dv )  + d/dv(Dc df/dv + \nu v f) + d/dr(D df/dr)
!            RF (q-l) souce term  Collisional terms        Spatial  Diffusion
!
!   multiply by e v and integrate over v,
!   integrating by parts in first 3 terms on RHS
!
!    dJ/dt = \int e(-Dq df/dv) dv + \int e(-Dc df/dv) dv - \int \nu e v f dv
!            RF q-l source term    Collisional diffusion  Dynamical friction
!
!           +                                                d/dr(D dJ/dr)
!                                                           Spatial Diffusion
!
!   In the dynamical friction term, pull \nu out of the integral
!   making it \nu_eff \cdot {fudge-factor}, but CALL lscsq_it \nu still.
!   Then you have - \nu J for that term
!
!   Tell yourself that if RF is there, the q-l source overwhelms the
!   collisional velocity diffusion, which we now ignore.  Then
!
!    dJ/dt = \int e(-Dq df/dv) dv  -\nu J + d/dr( D dJ/dr )
!
!   Tell yourself that the rf souce term can be written as
!   \nu^prime times the RF current found in the normal LSC calculation J_o;
!   and then tell yourself that we might as well treat \nu^prime == \nu
!
!   Then
!
!     dJ/dt = \nu (J_o - J) +  d/dr ( D dJ/dr )
!
!   This has some desirable properties, like if
!   D is small then the J moves to J_o on the \nu time scale
!   and if D is large then J and J_o can be quite different.
!   Also if \nu is space-dependent and D not then the J will
!   tend to move out (as Giruzzi finds) because there are
!   fewer collisions there.
!
!   It does not do anything like move the fast particles around
!   and change the damping because of moved fast particles.
!
!   The Edc treatment based on Karney Fisch is a confused
!   thing after this.
!
!
subroutine lscsq_SmoothJrf(Radius, DiffCoef, NumJs, idxFor0,        &
     &                     NuSlow, Jdriven, Jdiffus,                    &
     &                     Awk, Bwk, Cwk, Jwk, Xwk,                     &
     &                     IndeptVar, BoundCond)
  use iso_c_binding, only : fp => c_double
  implicit none
      INTEGER NumJs, idxFor0, idxFor0m1
      REAL(fp)    Radius, DiffCoef
      REAL(fp)    NuSlow(NumJs), Jdriven(NumJs), Jdiffus(NumJs)
      REAL(fp)    Awk(NumJs), Bwk(NumJs), Cwk(NumJs), Jwk(NumJs)
      REAL(fp)    Xwk(NumJs), xwk0
      REAL(fp)    TimeDum, Jp1
      INTEGER i, ierror
      CHARACTER*(*) IndeptVar, BoundCond
      REAL(fp)    BC1, BC2
 
      Jp1 = REAL(NumJs+1,kind=fp)

      if (idxFor0 .LT. NumJs) then
         idxFor0m1 = idxFor0 - 1
      else
         idxFor0m1 = NumJs
      endif
 
      !     Default Radial Independent Variable Type is:  PolFlux_Like
      !     This is Area_Like.
 
      if (IndeptVar .EQ. 'RadDist_Like') then
         TimeDum = Radius**2/(1.*DiffCoef) / Jp1**2
      else if (IndeptVar .EQ. 'PolFlux_Like')then
         TimeDum = Radius**2/(4.*DiffCoef) / Jp1**2
      else
         TimeDum = Radius**2/(4.*DiffCoef) / Jp1**2
      endif
 
      ! Zero slope boundary condition at i=0;i=1 or
      ! Zero second deriv boundary condition at i=0;i=1;i=2
      if(BoundCond .EQ. 'Zero1stDeriv') then
         BC1 = 1.0d0
         BC2 = 0.0d0
      else if(BoundCond .EQ. 'Zero2ndDeriv') then
         BC1 = 2.0d0
         BC2 =-1.0d0
      else
         BC1 = 1.0d0
         BC2 = 0.0d0
      endif
 
!     Xwk is filled with the average value of x(i) and x(i+1);
!     a more important purpose is workspace for the matrix inverter.
!     Xwk does not have to be filled above index idxFor0m1, but we
!     fill it anyway to emphasize that xwk goes with the geometry of
!     the problem, not with where we think the current must go to zero.
!
      do i=1,NumJs
         xwk(i) = (REAL(i,kind=fp)+0.5d0)/ Jp1
      enddo
      xwk0 = 0.5d0/Jp1
 
      Awk(1)     = 0.0d0
      Bwk(1)     = TimeDum*NuSlow(1)
      if(IndeptVar .EQ. 'RadDist_Like') Bwk(1) = Bwk(1)*1./Jp1
      Bwk(1)     = Bwk(1) + xwk0*(1.d0-BC1) + xwk(1)
      Cwk(1)     = - (xwk(1) + BC2*xwk0 )
!
!     Awk(1) and Cwk(NumJs) are never used;  setting to zero is cosmetic
!
      do i=2,idxFor0m1
         Awk(i) = - xwk(i-1)
         Bwk(i) = TimeDum*NuSlow(i)
         if(IndeptVar .EQ. 'RadDist_Like') Bwk(i) = Bwk(i)*REAL(i,kind=fp)/Jp1
         Bwk(i) = Bwk(i) + xwk(i-1) + xwk(i)
         Cwk(i) = - xwk(i)
      enddo
!
!     Zero value boundary condition at i=idxFor0
      Cwk(idxFor0m1) = 0.0d0
!
!     Awk(1) and Cwk(idxFor0-1) are never used;  setting to zero is cosmetic
!
      do i=1,idxFor0m1
         Jwk(i) = TimeDum*NuSlow(i)*Jdriven(i)
         if(IndeptVar .EQ. 'RadDist_Like') Jwk(i) = Jwk(i)*REAL(i,kind=fp)/Jp1
      enddo
 
!     Find the 'diffused' current;
!     Note that we are done with Xwk, and now pass it to tridiaNR for
!     it to use as workspace.
!
      CALL lscsq_tridiaNR(Awk,Bwk,Cwk,Jwk,Jdiffus,idxFor0m1,Xwk,ierror)
!     Fill remaining Jdiffused locations with zeros, if necessary.

      if (idxFor0 .LE. NumJs) then
         do i=idxFor0, NumJs
            Jdiffus(i) = 0.0d0
         enddo
      endif

      if (ierror .ne. 0) then
         do i=1,NumJs
            Jdiffus(i) = Jdriven(i)
         enddo
      endif
      return
      END
 
!     -----------------------------------------------------------------|
 
SUBROUTINE lscsq_tridiaNR (a,b,c, r,u,n,                            &
     &                              gam,ierror)
!     References:
!     Numerical Recipies in Fortran by William H. Press,
!     et al, Ch 2; p 40.
!
!     Solves for a vector  u  of length  n  the tridiagonal linear set:
!
!     b1   c1   0    0    0    0    0        u1          r1
!     a2   b2   c2   0    0    0    0        u2          r2
!     0    a3   b3   c3   0    0    0        u3          r3
!     0    0    a4   b4   c4   0    0        u4      =   r4
!     ..   ..   ..   ..   ..   ..   ..       ..          ..
!     ..   ..   ..   ..   aN-1 bN-1 cN-1     uN-1        rN-1
!     .              ..   0    aN   bN       uN          rN
!
!     a  b  c  r  are input vectors and are not modified.
!
  use iso_c_binding, only : fp => c_double
  implicit none
      INTEGER n, j, ierror
      REAL(fp)    a(n), b(n), c(n), r(n), u(n)
      REAL(fp)    bet
!     PARAMETER (NMAX = 200)  ! as given in the book
!     REAL    gam(NMAX)       ! as given in the book
      REAL(fp)    gam(n)
!
!     Trap error on ill-posed problem:
!     if (b(1) .eq. 0.) stop ' b1 is zero in tridiag' ! as given in the book
      ierror = 0
      if (b(1) .EQ. 0.d0) then
         ierror=1
         goto 1313
      endif
!
!     Normal beginning point:
      bet = b(1)
      u(1)= r(1)/bet
!
      do 11 j=2,n
        gam(j) = c(j-1)/bet
        bet    = b(j) - a(j)*gam(j)
!                if (bet .eq. 0.) stop ' tridiag fails' ! as given in the book
                 if (bet .EQ. 0.0d0) then
                    ierror=2
                    goto 1313
                 endif
        u(j)   = ( r(j) - a(j)*u(j-1))/bet
 11   continue
!
      do 12 j=n-1,1,-1
        u(j) = u(j)-gam(j+1)*u(j+1)
 12   continue
!
      return
!
!     Error condition...
 1313 continue
      return
!
      END
!

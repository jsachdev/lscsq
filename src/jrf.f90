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
!     E is negative, the quantity muminus becomes +1; muplus becomes -1
!
!     ------------------------------------------------------------------
!
subroutine lscsq_VmaxCalc
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : Vparmaxn,Vparmaxp, npsi
  use lscsq_mod, only : nzones, nrays
  use lscsq_mod, only : lh_out
  implicit none

  integer :: ip,ir,iz,iv
  real(fp):: vl

  VparMaxN(1:npsi)= -0.05_fp
  VparMaxP(1:npsi)= +0.05_fp
 
  do ir=1,nrays
     if (lh_out%ok_ray(ir).eq.1) then
        do iz=1,nzones
           iv = lh_out%ivind(iz,ir)
           ip = lh_out%izind(iz,ir)
           if (iv .EQ. 0 .or. ip .EQ. 0) exit
           vl = lh_out%vpar(iv)
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

  jray(1:nv,1:npsi) = 0.0_fp
  CALL lscsq_VmaxCalc
  nGotRuna = 0
  vnormnOK(1:npsi) = 0.0_fp

  do ip = 1, npsi
     CALL lscsq_jnorm(ip)
     !    set up normalization, tabulate u lookup values
     !    This is based on figure 2 of the Karney Fisch paper.  It shows that
     !    if you are trying to drive current in cooperation with the electric field
     !    at v/vrun you have had it.  If trying to drive against the electric field
     !    you can go up to v/vrun of 2 before runaways take over.
     if (EdcAry(ip) .GE. 0.0_fp ) then
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
 
#if DEBUG==2
  if (nGotRuna .GE. 1) then
     write(ErrMsg,'(1x,i4,'' Rway shells; last at ip='',i3)') nGotRuna, ipLastRuna
     call lscsq_LSCwarn(ErrMsg)
  endif
#endif

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
  use lscsq_mod, only: vnorm,muplus,muminus,gmrun,ugr,ugr
  use lscsq_mod, only: Edcary, ivrun, vnmax, nuruna, nv
  use lscsq_mod, only: lh_out
  implicit none

  integer, intent(in) :: ip ! flux surface index
  integer :: iv
  
  !     set up normalization, tabulate u lookup values
  gmrun = 1.0e-2_fp*lh_out%Ne(ip)*lh_out%logL(ip)*qe_eV**4/me_g**2*vc*2.0_fp*twopi 
  if(EdcAry(ip) .ne. 0.0) then
     vnorm = 1.0e-9_fp*sqrt(abs(lh_out%Ne(ip)/EdcAry(ip)))*    &
             sqrt(1e13_fp*qe_eV**3/me_g*vc**2 *2.0_fp*twopi*lh_out%logL(ip))
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
  ugr(1:nv) = lh_out%vpar(1:nv)/vnorm *muminus
 
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
  use lscsq_mod, only: muplus,jray,dvplus
  use lscsq_mod, only: ismo,ivzero,ugr,vtherm,muminus
  use lscsq_mod, only: ivlary, iitr, nuruna, nv
  use lscsq_mod, only: lh_inp, lh_out
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
  jd = 0.0_fp
  constfac = 1.0e6_fp*qe_eV*vc/nuruna   
  iGotRuna = 0
  iSMOi = mod(iSMO,2) + 1
  FeCutOff = FeCutFac*lh_out%fe(IvZero,ip,iITR)
  eps = sqrt(iVlary(ip)/(2.0_fp*PI**2*lh_inp%Raxis))/lh_inp%Raxis
  eps_p = eps**0.77_fp  ! calculate this only once
  xr = 3.5_fp    ! fmp - what is this?

  ! first negative velocities
  do iv = 1, IvZero - 1
     ! Ignore jd where f_e is small
     if (lh_out%fe(iv,ip,iITR).ge.FeCutOff ) then
        ! Muplus, because mu is plus for negative velocity if the Edc is >0
        ! See above for reversal if E<0
        CALL lscsq_WsloPrm(ugr(iv), muplus, lh_out%zbar(ip), dWsduou, iWhichWay)
        dWsdu = dWsduou*ugr(iv)
        ! Note that dWsdu*ugr is always positive
        if (iWhichWay.ne.0) iGotRuna = iGotRuna+1
        ! trapping effect from Ehst-Karney for LH limit of large vpar/vth
        w = abs(lh_out%vpar(iv))/Vtherm(ip)
        TrapFac = 1.0_fp-((eps_p*sqrt(xr**2+w**2))/(eps_p*xr+w))
        ! For negative velocity the current driven is negative because df/dv >0
        jdtemp = -lh_out%Dql(iv,ip,iSMO)*lh_out%dfdv(iv,ip,iSMOi)*TrapFac*dvplus(iv)*dWsdu*constfac
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
     if (lh_out%fe(iv,ip,iITR).GE.FeCutOff) then
        CALL lscsq_WsloPrm(ugr(iv), muminus, lh_out%zbar(ip), dWsduou, iWhichWay)
        dWsdu = dWsduou  * ugr(iv)
        ! Note that dWsdu*ugr is always positive
        if (iWhichWay.ne.0) iGotRuna = iGotRuna+1
        ! trapping effect from Ehst-Karney for LH limit of large vpar/vth
        w = abs(lh_out%vpar(iv))/Vtherm(ip)
        TrapFac = 1.0_fp-((eps_p*sqrt(xr**2+w**2))/(eps_p*xr+w))
        ! For positive velocity the current driven is positive because df/dv<0
        jdtemp = -lh_out%Dql(iv,ip,iSMO)*lh_out%dfdv(iv,ip,iSMOi)*TrapFac*dvplus(iv)*dWsdu*constfac
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
  use lscsq_mod, only: lh_inp, lh_out
  use lscsq_mod, only: Edcary, npsij, npsi
  implicit none

  real(fp), intent(in) :: DifAmt
  integer :: j
  real(fp):: Edc
  real(fp) :: MinEdc = 1.0e-4_fp
  ! DifAmt is intended to be 0.0 when the E_dc is desired as given
  ! and intended to be something like 0.0001 when we are trying to form
  ! the derivative   d ln J / d ln E

  do j = 1, npsi
     ! cannot use the map1d subroutine here because the DC electric field can
     ! have large gradients locally. This will be revised as the code is upgraded
     Edc = lh_out%Edc(j) + DifAmt
     if (abs(Edc) .LE. MinEdc) then
        if (Edc .GE. 0.0_fp ) then
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
  use lscsq_mod, only: ivrun, ismo, nv
  use lscsq_mod, only: lh_out
  implicit none

  ! Added Aug 93 to help quantify the runaway situation for LSC report PPPL 2929
  integer, intent(in) :: ip      ! index of flux surface
  integer :: iSMOi
  real(fp):: nDot, jDot, vRunAwayIndex
 
  iSMOi = mod(iSMO,2) + 1

  nDot= 0.0_fp
  JDot= 0.0_fp
  vRunAwayIndex= 0.0_fp
  if (ivrun .GT. 1 .and. ivrun .LT. nv ) then
     nDot = lh_out%Dql(ivrun,ip,iSMO) * lh_out%dfdv(ivrun,ip,iSMOi)
     jDot = lh_out%vpar(ivrun)*lh_out%Dql(ivrun,ip,iSMO) * lh_out%dfdv(ivrun,ip,iSMOi)
 
     nDot = abs ( nDot )
     jDot = abs(JDot)* 1.0e6_fp*qe_eV *vc
     vRunAwayIndex = REAL(ivrun,kind=fp)
  endif
 
end subroutine lscsq_ddtNrnJrn
 
 

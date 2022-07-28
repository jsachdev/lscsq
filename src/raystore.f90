SUBROUTINE lscsq_E2byPr(ninac1,accum1,rayqt1)

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : pi, twopi, eps0
  use lscsq_mod, only : qe_eV, me_g
  use lscsq_mod, only: Ezsq,distry,woc2,y,detrry, iray
  use lscsq_mod, only: izone, delpsi,Powrry,RofRay,zofray,rzind
  use lscsq_mod, only: pofray, nparry, nperry, lstop, npsij
  use lscsq_mod, only: timery, neofry, rtPsRy, Bthray, Bphray
  use lscsq_mod, only: d1, d2, d4, power, npar, iray, dlnPdsK, dlnPdsX
  use lscsq_mod, only: epsz, ecyc2, Epari, epql, epsl, dvol
  use lscsq_mod, only: omega, wdDdw, woc4, fghz
  use lscsq_mod, only: Epar, Eper, Exy, epsq
  use lscsq_mod, only: nzones, npsi
  use lscsq_mod, only: senter, tenter, sleave, tleave
  use lscsq_mod, only: lh_const, lh_out, lh_inp
  implicit none


!     E2byPr    Finds for all zones the ratio E_{z}^2/P;
!     damping constants  dlnPdsX  dlnPdsK;
!     n_{\parallel}; and polarization.
!     Written by D. W. Ignat, April 1991.
!     The time and path-position of entering and leaving
!     a zone is obtained by linear interpolation.  This is to allow
!     for best accuracy in obtaining  dt  and  ds .
!     This version does not anticipate crossing more than one zone
!     boundary in a step, but if this happens a warning is issued.
!
!     The approach is to average the quantity needed for each zone by
!     accumulating a total and dividing by the number of entries.  This
!     alone is sufficient for n_{\parallel} and polarization.
!     For E_z^2/P ... ezsq ... one needs to multiply by the time spent
!     in the zone, and divide by the volume of the zone.
!     For dlnPdsK  dlnPdsX we have chosen to multiply by   ds  at the end.
!     These variables are really, then,
!     dP/P -- Kernel to accept df/dv later; and dP/P -- assuming MaXwellian.
!
!     Zone boundaries are between the psi points, as follows.
!
!     +     .     +     .     +     .     +     .     +     .     +
! psimin                                                        psimax
!     1           2           3                     npsi-1        npsi
!  zone 1---><--zone 2--><--zone 3-->                        <--zone npsi
!
!     Ray power at the beginning is indexed 1.  As ray crosses into
!     the next zone, the value is indexed 2.  Power deposited in outer most
!     zone is therefore P_1 - P_2
!
!     sNew, sOld                New and old path length from start
!     tNew, tOld                and time duration from start.
!     RzindNew, RzindOld        The real and integer values of the zone
!     IzindNew, IzindOld        index.  These are related to the izind array.
!                               zind is small at the center and large at edge.
!                               izone=1 implies edge value of zind for
!                               the start of the ray.
!     Rzind* = (Psi - PsiMinx) / (PsiMaxx - PsiMinx) * (Npsi - 1) + 1.5
!     Rzind* = (Psi - PsiMinx) /  DelPsi                          + 1.5
!     Izind* = ifix(Rzind*)
!
!          Rzind*
!            ^
!     npsi   +
!            |    .                                  .
!     npsi-1 +
!            |         .    double jump              \
!     npsi-2 +            /                            double jump
!            |                                  .
!     .      +
!            |              .              .
!     .      +
!            |                   .    .
!     .      +
!            |
!     4      +
!   izind=3  |
!     3      +
!   izind=2  |
!     2      .____.____.____.____.____.____.____.____.____.
!                 0    1    2    3    4    5    6    7    8 > s (or) t
!
!
!     sSlope, tSlope    (New - Old)/(RzindNew - RzindOld)
!     sLeave(nzones), sEnter(nzones) are redundant: sLeave(i)=sEnter(i+1)
!     tLeave(nzones), tEnter(nzones) are redundant: tLeave(i)=tEnter(i+1)
!     accum(nRayQt)  an accumulator for averaging
!                    when the zone number is constant
!     RayQt(nRayQt)  the fresh parameter to be averaged
!     NinAc(nRayQt)  the number of entries now in the accumulator
!
!     RayQts are:
!     1 for ezsq    (epsz)^{-1} * 2.
!     2 for dlnPdsK dD/dK33 /(w/2)(dD/dw) / (dV/dwt)f * Im K33
!     3 for dlnPdsX but the _K means kernel for adding df/dv later
!                   and the _X means MaXwellian assumed.
!     4 npar
!     5 epolX
!     6 epolY
!
!     CLIGHT    speed of light, in 10^8 m/s
!     ee        electron epsilon == 1 + \omega_{pe}^2/\omega^2
!     eps0      (\mu_o c^2)^{-1} farads per meter; epsilon sub zero
!     epsz      epsilon sub z; converts E_z^2 to energy density
!     ex        E_x / E_z ; x-polarization
!     ey      i E_y / E_z ; y-polarization
!     PI        3.1415926
!     cEparIK   constant converts df_e/dv to Im{K_{33}}
!               ImEpar == - PI \omega_{pe}^2/k_{\parallel}^2 df/dv
!                      == - cEparIK df/dv / kpar2
!
!
!     [ K + n n - n^2 I ]  \cdot E = 0
!       =   - -       =          -
!
!     so if E = [ ex , i ey , 1 ] E_z , then
!     ex = - (Kzz - nperp^2)/(nperp npar)
!     ey = ex Kxy / (Kyy - n^2)
!
!
!     (1/P) (dP/dt) = - 2 Im(K_par) dD/dK_par / (dD/d\omega)
!     (1/P) (dP/ds) = - 2 Im(K_par) dD/dK_par / abs(dD/dk)
!           (dP/dV) = - 2 Im(K_par) dD/dK_par / (dD/d\omega) \times
!                         (eps_z  /2) E_z^2
!
!     P =  U  \vec{v_g} \cdot \vec{A} where A is area of ray
!     P = <U>  (dV/dt)  on flux surface average
!
!     U = (1/2) (eps_0 /2) E^* \cdot [K + d/d\omega(\omega K)] \cdot E
!                          -          =                    =         -
!          ^        ^                 ^           ^
!          |                      magnetic     electric and particle
!      avg amplitudes
!
!       = (E_z^2/2) (eps_0 /2) [ ex, -i ey, 1] \cdot
!   -                      -
!   |  2 ee    i Kxy    0  |   ex
!   | -i Kxy   2 ee     0  | i ey
!   |    0       0      2  |   1
!   -                      -
!
!     where   ee = (1 + \omega_{pe}^2 / \omega_{ce}^2)
!
!     U = E_z^2  (eps_0 /2) [ 1 + ee (ex^2 + ey^2) - Kxy ex ey ]
!       = E_z^2  (eps_z /2)
!
!     E_z^2 /P  = [(eps_z /2) dV/dt]^{-1}        == ezsq
!     (note that dt is normalized with \omega....to get MKS put back.)
!
!     dVol(NPSIDIM) is an array such that dVol(j) is the volume
!     centered on psi_j, etc.  This is used in calculating
!     dV/dt.  Once a ray enters a zone, it is assumed to be spread out
!     over all the volume of that zone.
!
!     Manifestly evident consistency with the QL form is consistent
!     with
!     U = E_z^2 (eps_L /2)
!     where
!                eps_L/eps0 = \omega/2 dD/d\omega / (dD/dEpar)
!
!
!     For some special studies, the components along field,
!     perp to field and along gradient, perp to field and along
!     flux surface are important.  Call the unit vectors \
!     e_11, e_rr, e_tt respectively.
!     Related unit vectors are toroidal, along gradient, poloidal,
!     e_ph, e_rr, e_pl
!
!     The code is written in unit vectors
!     e_ph, e_R , e_Z, usually referred to as
!     (e_R , e_Z, e_ph) , (y(1), y(2), y(3)) , (y(4), y(5), y(6))
!
!     e_rr =       Bz/Bp e_R     -    Br/Bp e_Z +     0 e_ph
!     e_tt = Br/Bp Bph/B e_R  + Bz/Bp Bph/B e_Z + -Bp/B e_ph
!     e_11 =       Br /B e_R  +        Bz/B e_Z + Bph/B e_ph
!
!     e_pl =       Br/Bp e_R  +       Bz/Bp e_z +     0 e_ph
!
!     B    =  sqrt ( Br^2 + Bz^2 + Bph^2 )
!     Bp   =  sqrt ( Br^2 + Bz^2 )

      INTEGER  lscsq_ivtabl

  character(len=70) :: ErrMsg

  real(fp), dimension(6), intent(inout) :: accum1, rayqt1
  integer, intent(inout) :: ninac1
  integer :: ifirst=1
  integer :: i, jzn, jry, IzindNew, IzindJmp
  real(fp) :: sOld, tOld 
  real(fp):: RzindNew, RzindCrs, sNew, tNew, sSlope, tSlope
  real(fp) :: rzindold
  integer :: izindold
  integer :: iv

  real(fp), dimension(4) :: bcompv
  real(fp) :: woc       
  real(fp) :: Btot, Bpol, Bphi, dDdEpar, dDdEparOld, ee, ex, ey, &
              Kpar, Kper, Kpar2, Kper2, psie2, qpar, veow2
  real(fp) :: psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
  real(fp) :: lscsq_DispRela, det, dum
  real(fp) :: NparDead = 9.8_fp 
  real(fp) :: NperDead = 147.0_fp
  real(fp) :: MaxwDead = 1.0e-09_fp !1.0e-12_fp
  real(fp) :: MxdPdz   = 0.90_fp
  real(fp) :: SMALL    = 1.0e-30_fp

  real(fp) :: re41 
!---------------------------------

!  DMC: always update these; frequency can change with ray index now
!
!  (omega was computed: setfreq SR, Rayini.F   omega   = fghz_now*TWO*PI*TENp9)
!
  woc = sqrt(woc2)
 
  !    Find location and various parameters.
  !    This is the NORMAL BEGINNING point, in that all initializations have been made
  !    and the ray is progressing through zones.

  CALL lscsq_plasma2d (y(1),y(2), psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
  Bpol  = sqrt(Br**2 + Bz**2)
  Bphi  = RBphi/y(1)
  Btot  = sqrt(Bpol**2 + Bphi**2)
  Kpar  = ( y(4)*Br + y(5)*Bz + y(6)/y(1)*Bphi)/Btot

  Kpar2 = Kpar**2
  Kper2 = Y(4)**2 + Y(5)**2 + (Y(6)/Y(1))**2 - Kpar2
  Kper  = sqrt(abs(Kper2))
  Qpar = Kpar2 - woc2*Eper
  dDdEparOld = (Qpar*(Qpar+Kper2)-woc4*Exy**2)
  dDdEpar    = (Qpar+kper2)*kpar2*kper2 / (kper2-woc2*Epar)

  EparI = 0.0_fp
  veow2 = 0.0445e-04_fp*tee/fghz(iray)**2
  psie2 = 2.0_fp*veow2*kpar2
  ! Trying to avoid overflows here.
  if ( psie2 .gt. 0.02_fp ) then
     psie2 = 1.0_fp/psie2
     EparI = 2.0_fp*sqrt(pi)*epsq*psie2*sqrt(psie2)*exp(-psie2)
  endif

  !     dDdEparOld = (Qpar*(Qpar + Kper2) - woc4*Exy**2)  is the term
  !     multiplying the Epar or K_{33} or K_{zz} in D.
  !     This is not manifestly positive, but it is believed to be positive.
  !     Now called Old.  We have had trouble with this evaluating wrong. 24apr92
  !
  !     Using the dispersion relation, dDdEpar can be rewritten
  !     dDdEpar = (Qpar + kper2) kpar2 kper2 / ( kper2 - woc2 Epar )
  !
  !     dDdkABS is partial of D wrt vector k, abs value
  !     The - is because of damping; the 2. is because power is field ^2

  ex    = -(woc2*Epar - Kper2)/(Kper*Kpar)
  ey    =(ex*(Eper - Kpar2/woc2) + Kper*Kpar/woc2) / Exy
  ee    = 1.0_fp + epsq/ecyc2
  epsz  = eps0*(1.0_fp + ee*(ex**2 + ey**2)-Exy*ex**2)
  epQL  = 0.5_fp * wdDdw / dDdEpar
  epsL  = eps0* epQL

  RayQt1(1) = 2.0_fp/epsL
  RayQt1(2) = +1.0_fp/epQL
  RayQt1(3) = -EparI/epQL
  RayQt1(4) = Kpar/woc
  RayQt1(5) = ex
  RayQt1(6) = ey

  tNew = y(7)
  sNew = y(8)

  if (izone.eq.1) then
     told = TimeRy(1,iray)
     sold = DistRy(1,iray)
     tenter(izone,iray) = told                 
     senter(izone,iray) = sold                
  else
     told = timeRy(izone-1,iray)
     sold = DistRy(izone-1,iray)
     tenter(izone,iray) = tleave(izone-1,iray)
     senter(izone,iray) = sleave(izone-1,iray)
  endif
  izindold = lh_out%izind(izone,iray)
  RzindOld = rzind(izone,iray)

!  if (izone.lt.10) write(*,*) 'izone=',izone
!  if (izone.lt.10) write(*,*) 'y:',y


  RofRay(izone,iray)   = y(1)
  ZofRay(izone,iray)   = y(2)
  PofRay(izone,iray)   = y(3)
  NparRy(izone,iray)   = Kpar/woc
  NperRy(izone,iray)   = Kper/woc
  rtPsRy(izone,iray)   = sqrt( (psi-lh_inp%plflx(1))/(lh_inp%plflx(npsij)-lh_inp%plflx(1)) )
  TimeRy(izone,iray)   = y(7)
  DistRy(izone,iray)   = y(8)
  NeofRy(izone,iray)   = pe2/lh_const%Pe2Fac * 1.0e+14_fp
  BthRay(izone,iray)   = Bpol
  BphRay(izone,iray)   = Bphi
  det = lscsq_DispRela ( y(1) , y(2) , y(4) , y(5) , y(6)  )
  d1 = abs(d1)
  d2 = abs(d2)
  d4 = abs(d4)
  det = det/max(d1, d2, d4)

!  if (izone.lt.10) write(*,*) 'det=',det

  DetrRy(izone,iray) = det

  RzindNew = (Psi-lh_inp%plflx(1))/DelPsi + 1.5_fp
  RE41     =  RzindNew
  IzindNew = int(RE41)
!  izindnew = minloc(abs(psiary-psi),1)
  ! START major IF/ELSE/ENDIF branch.
  ! If we are in the same zone as before, put RayQt in accumulator.
  
  if (IzindNew.EQ.IzindOld) then
     do i=1,6
        accum1(i)=accum1(i)+RayQt1(i)
     enddo
     NinAc1 = NinAc1 + 1
     ! We are in a new zone.
     ! Divide the accumulator by # entries but take care if NinAc = 0.
     ! Interpolate to find crossings.
  else
     do i=1,6
        if (NinAc1 .EQ. 0) then
           accum1(i)=RayQt1(i)
        else
           accum1(i)=accum1(i)/REAL(NinAc1,kind=fp)
        endif
     enddo
     NinAc1 = 0

     sSlope   = (sNew-sOld)/(RzindNew-RzindOld)
     tSlope   = (tNew-tOld)/(RzindNew-RzindOld)
!     RE41     = RzindNew
!     IzindNew = int(RE41)
     ! The crossing point is different depending on whether the ray is going
     ! outward:   IzindNew > IzindOld
     ! or inward: IzindNew < IzindOld
     if (IzindNew.GT.IzindOld) then
        RzindCrs = REAL(IzindOld,kind=fp) + 1.0_fp
     else
        RzindCrs = REAL(IzindOld,kind=fp)
     endif
     sLeave(izone,iray) = sOld+sSlope*(RzindCrs-RzindOld)
     tLeave(izone,iray) = tOld+tSlope*(RzindCrs-RzindOld)

     ! Compute the desired parameters for the zone we just left.
     ! npar according to the velocity table 15Jan93
     npar(izone,iray)  = accum1(4)
     lh_out%ivind(izone,iray) = lscsq_ivtabl(npar(izone,iray))
     iv = lh_out%ivind(izone,iray)
     npar(izone,iray)  = 1.0_fp/lh_out%vpar(iv) ! fmp-why this replacement?
     lh_out%izind(izone,iray) = IzindOld

     ezsq(izone,iray) = accum1(1)*(tLeave(izone,iray)-tEnter(izone,iray))/omega/dVol(IzindOld)
     dlnPdsK(izone,iray) = accum1(2)*lh_const%cEparik/(npar(izone,iray)*woc)**2 * (tLeave(izone,iray)-tEnter(izone,iray))
     dlnPdsX(izone,iray) = accum1(3)*(tLeave(izone,iray)-tEnter(izone,iray))

     if (ezsq(izone,iray) .LT. 0.0) then
        write(ErrMsg,'('' Ez2<0!zn ind iry epsL epsZ:'', i6,i6,i6,1pe10.2,1x,1pe10.2)')   &
                     izone, IzindNew, iray, epsL, epsZ
        CALL lscsq_LSCwarn( ErrMsg)
        ezsq(izone,iray)   = 0.0_fp
        dlnPdsK(izone,iray)= 0.0_fp
        dlnPdsX(izone,iray)= 0.0_fp
     endif
     ! Dont forget to save the integers.

     ! We are done with the information accumulated for the old zone, so
     ! set the accumulator to the value of the parameters we just calculated. Reset counter.
     accum1=RayQt1
     NinAc1 = 1

     ! If we have more zones to go, then INCREMENT THE ZONE COUNTER, and
     ! fill the arrays with the same values we just found, in case we run out of time steps (eg).
     if (izone .LT. nzones) then
        izone = izone + 1
        ! copy everything to the arrays for the next step. Values withh be
        ! updated/replaced by this subroutine as needed
        RofRay(izone,iray)        = RofRay(izone-1,iray)
        ZofRay(izone,iray)        = ZofRay(izone-1,iray)
        PofRay(izone,iray)        = PofRay(izone-1,iray)
        NparRy(izone,iray)        = NparRy(izone-1,iray)
        NperRy(izone,iray)        = NperRy(izone-1,iray)
        rtPsRy(izone,iray)        = rtPsRy(izone-1,iray)
        NeofRy(izone,iray)        = NeofRy(izone-1,iray)
        BthRay(izone,iray)        = BthRay(izone-1,iray)
        BphRay(izone,iray)        = BphRay(izone-1,iray)

        ezsq   (izone,iray)  = ezsq   (izone-1,iray)
        dlnPdsK(izone,iray)  = dlnPdsK(izone-1,iray)
        dlnPdsX(izone,iray)  = dlnPdsX(izone-1,iray)
        npar   (izone,iray)  = npar   (izone-1,iray)
        lh_out%ivind  (izone,iray)  = lh_out%ivind(izone-1,iray)
        dum = 1.0_fp+dlnPdsX(izone-1,iray)
        dum =  max (MxdPdZ,dum)
        dum =  min (dum,1.0_fp)

        PowrRy (izone,iray)       = PowrRy(izone-1,iray)*dum
        lh_out%izind  (izone,iray)  = IzindNew
     else
        ! But, if we used up all the zones, quit as fast as you can.
        izone = nzones
        Lstop = 1
     endif
 
     ! If Npar is so large that damping would be total, stop the calculation now.
     if ( abs(npar(izone,iray)) .GT. NparDead ) Lstop = 1
     ! If Nper is really large, stop now. No matter what.
     if ( abs(nperRy(izone,iray)).GT.NperDead ) Lstop = 1
     ! If all power is gone in Maxwellian stop the calculation now.
     if ( PowrRy(izone,iray) .LT. MaxwDead ) then
!     if ( (PowrRy(izone,iray) .LT. MaxwDead).and. &
!        (abs(PowrRy(izone,iray)-PowrRy(izone-1,iray)).lt.1e-6 )) then
!     if (abs(PowrRy(izone,iray)-PowrRy(izone-1,iray))/PowrRy(izone,iray).lt.1.0e-03) then
        Lstop = 1
     endif
  endif
  ! See if we jumped more than one zone.  Report if so.
  IzindJmp = IzindNew - IzindOld
  if (IzindJmp*IzindJmp .GT. 1) then
     if(izone .GE. nzones/10) CALL lscsq_LSCwarn( ' More than 1 zone jumped and izone >> 1 ')
  endif
  ! Set the  _Old parameters to the existing _New ones.
  sOld  = sNew
  tOld  = tNew
  IzindOld = IzindNew
  RzindOld = RzindNew

end subroutine lscsq_E2byPr
!
!                                                                      |
subroutine lscsq_RyZnInit

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: nzones, nrays
  use lscsq_mod, only: lh_out
  use lscsq_mod, only: rofray, zofray 
  use lscsq_mod, only: distry, detrry, ezsq, timery, powrry, npar, ntor
  use lscsq_mod, only: dlnPds, rtPsRy, nperry, pofray, dlnPdsK,dlnPdsX, NparRy
  implicit none

  !     Initialize  ray zones    -----------------------|
  RofRay (1:nzones,1:nrays) = 0.0_fp
  ZofRay (1:nzones,1:nrays) = 0.0_fp
  PofRay (1:nzones,1:nrays) = 0.0_fp
  NperRy (1:nzones,1:nrays) = 0.0_fp
  NparRy (1:nzones,1:nrays) = 0.0_fp
  rtPsRy (1:nzones,1:nrays) = 0.0_fp
  PowrRy (1:nzones,1:nrays) = 0.0_fp
  TimeRy (1:nzones,1:nrays) = 0.0_fp
  DistRy (1:nzones,1:nrays) = 0.0_fp
  DetrRy (1:nzones,1:nrays) = 0.0_fp

  ezsq (1:nzones,1:nrays) = 0.0_fp
  dlnPds (1:nzones,1:nrays) = 0.0_fp
  dlnPdsK(1:nzones,1:nrays) = 0.0_fp
  dlnPdsX(1:nzones,1:nrays) = 0.0_fp

  lh_out%izind(1:nzones,1:nrays) = 0
  lh_out%ivind(1:nzones,1:nrays) = 0
  npar(2:nzones,1:nrays) = 0.0_fp   ! But do npar separately since the

  NparRy(1,1:nrays) = npar(1,1:nrays)

end subroutine lscsq_ryzninit


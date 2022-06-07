module rayStore
   contains 
   SUBROUTINE lscsq_E2byPr(y1, y2, y3, y4, y5, y6, y7, y8, izone, lstop, arrys, &
                           p_data, a_data, o_data)!, local_prmtrs)!, all_vecs)

    
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
   !     tleave(nzones), arrys(izone-1)%tenter(nzones) are redundant: tleave(i)=arrys(izone-1)%tenter(i+1)
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
   !     (e_R , e_Z, e_ph) , (y1, y2, y3) , (y4, y5, y6)
   !
   !     e_rr =       Bz/Bp e_R     -    Br/Bp e_Z +     0 e_ph
   !     e_tt = Br/Bp Bph/B e_R  + Bz/Bp Bph/B e_Z + -Bp/B e_ph
   !     e_11 =       Br /B e_R  +        Bz/B e_Z + Bph/B e_ph
   !
   !     e_pl =       Br/Bp e_R  +       Bz/Bp e_z +     0 e_ph
   !
   !     B    =  sqrt ( Br^2 + Bz^2 + Bph^2 )
   !     Bp   =  sqrt ( Br^2 + Bz^2 )
                           
                           
   use pl_types 
   use rTrace
   use Tracing
   use iso_c_binding, only : fp => c_double
   use lscsq_mod, only : pi, eps0
   !use lscsq_mod, only : qe_eV, me_g
   use lscsq_mod, only: cEparIK, delpsi
   use lscsq_mod, only: psimin, psilim
   use lscsq_mod, only: pe2fac
   use lscsq_mod, only: dvol, vpar
   !use lscsq_mod, only: nzones!, nrays!, nrayqt
   !use lscsq_mod, only: NparDead, NperDead, MaxwDead, MxdPdz, SMALL

   
   implicit none
  ! integer, intent(IN) :: iray
   integer, intent(INOUT) :: lstop, izone
   !type(parameters), intent(inout) :: local_prmtrs

   type(pl_data), intent(inout) :: p_data
   type(accum_data),  intent(inout) :: a_data
   type(old_data), intent(inout) :: o_data
   type(storeAry), dimension(:), intent(inout) :: arrys
   Real(fp) :: y1, y2, y3, y4, y5, y6, y7, y8
   Real(fp) :: max_d, Epari, r, z
   INTEGER  lscsq_ivtabl, nzones
   integer :: i, IzindNew, IzindJmp
   Real(fp):: RzindNew, RzindCrs, sNew, tNew, tSlope!, sSlope
   Real(fp) :: Btot, Bpol, Bphi, dDdEpar, ee, ex, ey, &
               Kpar, Kper, Kpar2, Kper2, psie2, qpar, veow2
   Real(fp) :: Br,Bz,RBphi,omc,Tee,pi2,aio,ael, psi, pe2
   Real(fp) :: det, dum, npar, aion, ipsq
   Real(fp) :: epQL, epsL, epsz!, psimin, psilim

   Real(fp) :: NparDead = 9.8_fp 
   Real(fp) :: NperDead = 147.0_fp
   Real(fp) :: MaxwDead = 1.0e-12_fp
   Real(fp) :: MxdPdz   = 0.90_fp
   Real(fp) :: SMALL    = 1.0e-30_fp
   
   !type(vecs), dimension(npsij), intent(inout):: all_vecs

   !$acc routine seq
   !$acc routine(lscsq_plasma2d)
   !$acc routine(lscsq_LSCwarn)
   !$acc routine(lscsq_disprela)
   !$acc routine(lscsq_ivtabl)


   !---------------------------------

   !  DMC: always update these; frequency can change with ray index now
   !
   !  (omega was computed: setfreq SR, Rayini.F   omega   = fghz_now*TWO*PI*TENp9)
   !
   
   !    Find location and various parameters.
   !    This is the NORMAL BEGINNING point, in that all initializations have been made
   !    and the ray is progressing through zones.
   
   r   = y1
   z   = y2

   nzones = size(arrys)
   !write(*,*) 'zones', nzones

   call lscsq_plasma2d(o_data, p_data, r,z,psi,Br,&
                        Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)!,all_vecs)
   
   Bpol  = sqrt(Br**2 + Bz**2)
   Bphi  = RBphi/y1
   Btot  = sqrt(Bpol**2 + Bphi**2)
   Kpar  = ( y4*Br + y5*Bz + y6/y1*Bphi)/Btot

   Kpar2 = Kpar**2
   Kper2 = y4**2 + y5**2 + (y6/y1)**2 - Kpar2
   Kper  = sqrt(abs(Kper2))
   Qpar = Kpar2 - p_data%woc2*p_data%Eper
   !dDdEparOld = (Qpar*(Qpar+Kper2)-p_data%woc4*p_data%Exy**2)
   dDdEpar    = (Qpar+kper2)*kpar2*kper2 / (kper2-p_data%woc2*p_data%Epar)

   EparI = 0.0_fp
   veow2 = 0.0445e-04_fp*tee/p_data%fghz**2
   psie2 = 2.0_fp*veow2*kpar2
   ! Trying to avoid overflows here.
   
   if ( psie2 .gt. 0.02_fp ) then
      psie2 = 1.0_fp/psie2
      EparI = 2.0_fp*sqrt(pi)*p_data%epsq*psie2*sqrt(psie2)*exp(-psie2)
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

   ex    = - (p_data%woc2*p_data%Epar - Kper2)/(Kper*Kpar)
   ey    =(ex*(p_data%Eper - Kpar2/p_data%woc2) + Kper*Kpar/p_data%woc2)/p_data%Exy
   ee    = 1.0_fp + p_data%epsq/p_data%ecyc2
   epsz  = eps0*(1.0_fp + ee*(ex**2 + ey**2)-p_data%Exy*ex**2)
   epQL  = 0.5_fp * p_data%wdDdw / dDdEpar
   epsL  = eps0* epQL
   a_data%RayQt1 = 2.0_fp/epsL
   a_data%RayQt2 = 1.0_fp/epQL
   a_data%RayQt3 = -EparI/epQL
   a_data%RayQt4 = Kpar/p_data%woc
   ! a_data(5)%RayQt = ex
   ! a_data(6)%RayQt = ey

   tNew = y7
   sNew = y8

   if (izone.eq.1) then
      o_data%told = arrys(1)%TimeRy
      o_data%tenter = o_data%told
   !else
      !o_data%told = arrys(izone-1)%timeRy
   endif

   arrys(izone)%RofRay   = y1
   arrys(izone)%ZofRay   = y2
   arrys(izone)%PofRay   = y3
   arrys(izone)%Nparry   = Kpar/p_data%woc
   arrys(izone)%NperRy   = Kper/p_data%woc
   arrys(izone)%rtPsRy   = sqrt( (psi-psimin)/(psilim-psimin) )
   arrys(izone)%TimeRy   = y7
   arrys(izone)%DistRy   = y8
   arrys(izone)%NeofRy   = pe2/Pe2Fac * 1.0e+14_fp
   arrys(izone)%BthRay   = Bpol
   det = lscsq_DispRela (y1 , y2 , y4 , y5 , y6, o_data, p_data, max_d  )
   det = det/ max_d 

   arrys(izone)%DetrRy = det

   RzindNew = (Psi-PsiMin)/DelPsi + 1.5_fp
   IzindNew = int(RzindNew)

  ! write(*,*) 'Rz: ', RzindNew, IzindNew

   !  izindnew = minloc(abs(psiary-psi),1)
   ! START major IF/ELSE/ENDIF branch.
   ! If we are in the same zone as before, put RayQt in accumulator.
   if (IzindNew.EQ.o_data%IzindOld) then
      !do i=1,nRayQt
         a_data%accum1=a_data%accum1+a_data%RayQt1
         a_data%accum2=a_data%accum2+a_data%RayQt2
         a_data%accum3=a_data%accum3+a_data%RayQt3
         a_data%accum4=a_data%accum4+a_data%RayQt4

      !enddo
      o_data%NinAc = o_data%NinAc + 1
      ! We are in a new zone.
      ! Divide the accumulator by # entries but take care if NinAc = 0.
      ! Interpolate to find crossings.
   else
    !  do i=1,nRayQt
         if (o_data%NinAc .EQ. 0) then
            a_data%accum1=a_data%RayQt1
            a_data%accum2=a_data%RayQt2
            a_data%accum3=a_data%RayQt3
            a_data%accum4=a_data%RayQt4
         else
            a_data%accum1=a_data%accum1/REAL(o_data%NinAc,kind=fp)
            a_data%accum2=a_data%accum2/REAL(o_data%NinAc,kind=fp)
            a_data%accum3=a_data%accum3/REAL(o_data%NinAc,kind=fp)
            a_data%accum4=a_data%accum4/REAL(o_data%NinAc,kind=fp)
         endif
      !enddo
      o_data%NinAc = 0

      tSlope   = (tNew-o_data%tOld)/(RzindNew-o_data%RzindOld)

      ! The crossing point is different depending on whether the ray is going
      ! outward:   IzindNew > IzindOld
      ! or inward: IzindNew < IzindOld
      if (IzindNew.GT.o_data%IzindOld) then
         RzindCrs = REAL(o_data%IzindOld,kind=fp) + 1.0_fp
      else
         RzindCrs = REAL(o_data%IzindOld,kind=fp)
      endif

      o_data%tleave = o_data%tOld+tSlope*(RzindCrs-o_data%RzindOld)

      ! Compute the desired parameters for the zone we just left.
      ! npar according to the velocity table 15Jan93
      npar  = a_data%accum4
      arrys(izone)%ivind = lscsq_ivtabl(npar)
      npar = 1.0_fp / vpar(arrys(izone)%ivind) ! fmp-why this replacement?
      arrys(izone)%izind = o_data%IzindOld

      arrys(izone)%ezsq = a_data%accum1*(o_data%tleave-&
      o_data%tenter)/p_data%omega/&
                                 dVol(o_data%IzindOld)

      arrys(izone)%dlnPdsK = a_data%accum2*&
      cEparik/(npar*&
                                    p_data%woc)**2 * &
                                    (o_data%tleave-&
                                    o_data%tenter)
      arrys(izone)%dlnPdsX = a_data%accum3*(o_data%tleave-&
      o_data%tenter)/(abs(dVol(o_data%izindOld))+SMALL)

     ! write(*,*) 'save: ',  tNew,o_data%tOld, (RzindNew-o_data%RzindOld), RzindNew, o_data%RzindOld, Psi
      ! Dont forget to save the integers.
      ! We are done with the information accumulated for the old zone, so
      ! set the accumulator to the value of the parameters we just calculated. Reset counter.
      a_data%accum1=a_data%RayQt1
      a_data%accum2=a_data%RayQt2
      a_data%accum3=a_data%RayQt3
      a_data%accum4=a_data%RayQt4
      o_data%NinAc = 1

      ! If we have more zones to go, then INCREMENT THE ZONE COUNTER, and
      ! fill the arrays with the same values we just found, in case we run out of time steps (eg).
      if (izone .LT. nzones) then
         
         o_data%tenter= o_data%tleave
         o_data%told = arrys(izone)%timeRy
         izone = izone + 1

         ! copy everything to the arrays for the next step. Values withh be
         ! updated/replaced by this subroutine as needed
         arrys(izone)%RofRay       = arrys(izone-1)%RofRay
         arrys(izone)%ZofRay       = arrys(izone-1)%ZofRay
         arrys(izone)%PofRay       = arrys(izone-1)%PofRay
         arrys(izone)%Nparry       = arrys(izone-1)%NparRy
         arrys(izone)%NperRy        = arrys(izone-1)%NperRy
         arrys(izone)%rtPsRy        = arrys(izone-1)%rtPsRy
         arrys(izone)%NeofRy        = arrys(izone-1)%NeofRy
         arrys(izone)%BthRay        = arrys(izone-1)%BthRay
         arrys(izone)%ezsq     = arrys(izone-1)%ezsq  
         arrys(izone)%dlnPdsK  = arrys(izone-1)%dlnPdsK
         arrys(izone)%dlnPdsX  = arrys(izone-1)%dlnPdsX
         arrys(izone)%ivind  = arrys(izone-1)%ivind  
         dum = 1.0_fp+arrys(izone-1)%dlnPdsX
         dum =  max (MxdPdZ,dum)
         dum =  min (dum,1.0_fp)

         arrys(izone)%PowrRy = arrys(izone-1)%PowrRy*dum
         arrys(izone)%izind   = IzindNew
      else
         ! But, if we used up all the zones, quit as fast as you can.
         izone = nzones
         Lstop = 1
      endif
   
      ! If Npar is so large that damping would be total, stop the calculation now.
      if ( abs(npar) .GT. NparDead ) then
         Lstop = 1
      ! write(*,*) 'stop#1 ', npar , NparDead
      endif
      ! If Nper is really large, stop now. No matter what.
      if ( abs(arrys(izone)%nperRy).GT.NperDead ) then
         Lstop = 1
      ! write(*,*) 'stop#2 ', arrys(izone)%nperRy,NperDead
      endif
      ! If all power is gone in Maxwellian stop the calculation now.
      if ( arrys(izone)%PowrRy .LT. MaxwDead ) then
         Lstop = 1
          ! write(*,*) 'stop#2 ', arrys(izone)%PowrRy , MaxwDead
      
         ! need to transform these flags into a coding for stop reason in the CDF output
      endif
   endif
   ! See if we jumped more than one zone.  Report if so.
   IzindJmp = IzindNew - o_data%IzindOld
   if (abs(IzindJmp) .GT. 1) then
      write(*,*) 'More: ', IzindNew , o_data%IzindOld, r,z
      !if(izone .GE. nzones/10) CALL lscsq_LSCwarn( ' More than 1 zone jumped and izone >> 1 ')
   endif
   ! Set the  _Old parameters to the existing _New ones.
   o_data%IzindOld = IzindNew
   o_data%RzindOld = RzindNew

   end subroutine lscsq_E2byPr
end module rayStore 
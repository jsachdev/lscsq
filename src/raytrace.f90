subroutine lscsq_DoRay
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : twopi, vc
  use lscsq_mod, only : enpar, enpol, ntor, npol
  use lscsq_mod, only : iray, izone, ierror, ok_ray, fghz
  use lscsq_mod, only: npols, ntors, nrays, Ezsq, ind_ray
  use lscsq_mod, only: thgrid, thet0, enth, NparRy, nz_ind
  use lscsq_mod, only: ezsq, ivind, izind, dlnPdsK, dlnPdsX, npar
  use lscsq_mod, only: omega, woc2, woc4
  implicit none

  integer :: ity, ipy, ith, RayIniErr, iFillPsi, iErrCount
  integer :: ind

  integer, dimension(:), allocatable :: ninac
  real(fp), dimension(:,:), allocatable :: accum
  real(fp), dimension(:,:), allocatable :: rayqt

  if(.not.allocated(ninac)) allocate(ninac(nrays))
  if(.not.allocated(accum)) allocate(accum(6,nrays))
  if(.not.allocated(rayqt)) allocate(rayqt(6,nrays))

  iFillPsi =1

  iErrCount = 1

  ! The quantities to be collected (RayQt) are averaged over the
  ! zone by summing into accum, dividing by NinAc. This clears.
  accum = 0.0_fp
  rayqt = 0.0_fp
  NinAc = 0

  ! initialize all arrays for ray tracing
  CALL lscsq_RyZnInit

  do iray=1,nrays
     ity=ind_ray(1,iray)
     ipy=ind_ray(2,iray)
     ith=ind_ray(3,iray)
    
     enpar = ntor(ity)
     enpol = npol(ipy)
     enth  = thgrid(ith)
     izone = 1
     omega=  1.0e09_fp * twopi*fghz(iray)
     woc2 = (omega/vc)**2
     woc4 =  woc2**2

     CALL lscsq_RayIni(RayIniErr)
     if (RayIniErr .GE. 1) then
        ok_ray(iray) = 0
        iErrCount = iErrCount+1
        cycle   
     else
        call lscsq_predcLSC(ninac(iray),accum(:,iray),rayqt(:,iray))
     endif
  enddo

  if (iErrCount .GE. nrays) then
     iError = 1
     CALL lscsq_LSCtrace('RayIni')
  endif

end subroutine lscsq_doray
!---------------------------------------
subroutine lscsq_PredcLSC(ninac1,accum1,rayqt1)

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : pi, deg2rad
  use lscsq_mod, only : neqsp1, neqs, ierror, begin, lstop
  use lscsq_mod, only : f1, f2, f3, f
  use lscsq_mod, only : y1, y2, y3, y
  use lscsq_mod, only : hstplh, nstep 
  implicit none

!!!!!!     EXTERNAL RungeLSC, ftion, prtout
!     These integer variables allow stopping the
!     integration when some condition is met.
!     Lstop    is meant for serious conditions
!              such as going out of bounds.
!     The include statements are to catch: begin,h,nstep,lstop,
!     and especially NEQS and NEQSP1

  integer :: i,j, jstart, BoundsEr 
  integer :: iBndsErr = 0
  integer :: nBndsErr = 25 
  real(fp), dimension(6), intent(inout) :: accum1, rayqt1
  integer, intent(inout) :: ninac1
  real(fp):: yok(NEQSP1), thi, tho
  Real(fp):: pc(NEQS), c(NEQS),p(NEQS),dumy(NEQS)
  Real(fp) :: k1 = 2.65277776_fp 
  Real(fp) :: k2 = -1.4861111_fp 
  Real(fp) :: k3 = 1.5138888_fp
  Real(fp) :: k4 = -0.34722222_fp
  Real(fp) :: k5 = 0.94266666_fp
  Real(fp) :: k6 = 0.34722222_fp
  Real(fp) :: k7 = 1.2638888_fp
  Real(fp) :: k8 = 0.59722222_fp
  Real(fp) :: k9 = 0.125_fp
  Real(fp) :: k10= 0.05733333_fp


  jstart = 4
 10   continue

  pc(1:neqs) = 0.0_fp

  BoundsEr = 0
  lstop = 0
  CALL lscsq_RungeLSC(ninac1,accum1,rayqt1,BoundsEr)
  if (BoundsEr .ne. 0 ) then
     CALL lscsq_LSCendr(' Cant recover in PredcLSC')
     return
  endif

  yok(1:neqsp1) = y(1:neqsp1)


  do j = jstart, nstep
     yok(1:neqsp1) = y(1:neqsp1)
     if (iError.GE.1) return
     if (lstop.EQ.1) then
        call lscsq_E2byPr(ninac1,accum1,rayqt1)
        return
     endif
     y(NEQSP1) = begin+REAL(j,kind=fp)*HstpLH
     do i=1,NEQS
        dumy(i) = (y3(i)+y3(i)+y2(i))/3.0_fp
        p(i)=dumy(i)+HstpLH*(k1*f(i)+k2*f3(i)+k3*f2(i)+k4*f1(i))
        y1(i)=y2(i)
        y2(i)=y3(i)
        y3(i)=y(i)
        ! modify
        y(i)=p(i)-k5*pc(i)
        f1(i) = f2(i)
        f2(i) = f3(i)
        f3(i) = f(i)
     enddo
     CALL lscsq_ftion(BoundsEr)
     if(BoundsEr .ne. 0) go to 100
     do i = 1,NEQS
        ! correct
        c(i)= dumy(i)+HstpLH*(k6*f(i)+k7*f3(i)+k8*f2(i)+k9*f1(i))
        pc(i) = p(i)-c(i)
        ! final
        y(i) = c(i)+k10*pc(i)
     enddo

     call lscsq_E2byPr(ninac1,accum1,rayqt1) ! increment the zone counter
     call lscsq_ftion(BoundsEr)
     if(BoundsEr .ne. 0) go to 100
  enddo 
  return
  ! Attempt to recover from ray error, most likely running out of bounds,
  ! by specular reflection off the last good set of values.
100  continue

  jstart = j
!  if (ScatKdeg.gt.1.0e-3_fp) then
!     CALL lscsq_BounShft(yok, y, NEQSP1, ScatKdeg, thi, tho)
!     inciThet(iscatplt) = thi * deg2rad
!     scatThet(iscatplt) = tho * deg2rad
!     iscatplt=iscatplt+1
!  else
     CALL lscsq_BounceIt(yok, y, NEQSP1)
!  endif

  iBndsErr = iBndsErr + 1
  if (iBndsErr .GE. nBndsErr) then
     iBndsErr = 0
     write(*,'('' Rays out of bounds; recovered'',i4,'' times'')') nBndsErr
  endif
  goto 10

end subroutine lscsq_PredcLSC
!
!     -----------------------------------------------------------------
!
SUBROUTINE lscsq_RungeLSC(ninac1,accum1,rayqt1,BoundsEr)

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : neqsp1, neqs, hstplh
  use lscsq_mod, only : f1, f2, f3, f
  use lscsq_mod, only : y1, y2, y3, y
  
  implicit none

  integer :: ia,ib,ic,  i,k
  integer :: BoundsEr
  real(fp), dimension(6), intent(inout) :: accum1, rayqt1
  integer, intent(inout) :: ninac1
  real(fp), dimension(4) :: a = [0.5_fp, 0.29289322_fp, 1.70710675_fp, 0.1666667_fp]
  real(fp), dimension(4) :: b = [2.0_fp, 1.0_fp       , 1.0_fp       , 2.0_fp] 
  real(fp), dimension(4) :: c = [0.5_fp, 0.29289322_fp, 1.70710675_fp, 0.5_fp]
  real(fp):: q(NEQSP1)
  real(fp):: dum

  f(NEQSP1) = 1.0_fp
  q(NEQSP1) = 0.0_fp
  CALL lscsq_ftion(BoundsEr)
  if(BoundsEr .ne. 0) go to 100
  y1(1:neqs) = y(1:neqs)
  f1(1:neqs) = f(1:neqs)
  q(1:neqs) = 0.0_fp

  do ia=2,4
     do ib=1,4
        call lscsq_ftion(BoundsEr)
        if(BoundsEr .ne. 0) go to 100
        do ic=1,NEQSP1
           dum = a(ib)*(f(ic)-b(ib)*q(ic))
           y(ic) = y(ic) + HstpLH*dum
           q(ic) = q(ic) + 3.0_fp*dum-c(ib)*f(ic)
        enddo
     enddo
     k = ia-1
     call lscsq_E2byPr(ninac1,accum1,rayqt1)
     if (k.eq.1) then
        f2(1:neqs) = f(1:neqs)
        y2(1:neqs) = y(1:neqs)
     endif
     if (k.eq.2) then
        f3(1:neqs) = f(1:neqs)
        y3(1:neqs) = y(1:neqs)
     endif
     if (k.eq.3) cycle
  enddo

  return

100  continue

  BoundsEr = 1
  call lscsq_LSCwarn(' Boundary error in RungeLSC')

end subroutine lscsq_RungeLSC
!
!     -----------------------------------------------------------------
!
SUBROUTINE lscsq_BounceIt (yok, y, n)

  use iso_c_binding, only : fp => c_double
  implicit none

  integer, intent(in) :: n
  integer :: i

  real(fp), dimension(n), intent(inout) :: yok
  real(fp), dimension(n), intent(out):: y
  real(fp) ::    r,z,   psi, Br, Bz, RBphi, omc, Tee, pe2, pi2, aio, ael
  real(fp) ::    B2, crr, crz, czr, czz, KrNew, KzNew
!     Go back to the last ok solution point, and re-arrange kr and kz
!     such that k \dot B is unchanged but k \cross B changes sign.  Then
!     put this new information into the starting condition for y.
!     y(1)  y(2)  y(3)  y(4)  y(5)  y(6)
!      r     z    phi   k_r   k_z    n

  r  = yok(1)
  z  = yok(2)
  CALL lscsq_plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)

  B2  =  Br*Br + Bz*Bz
  crr = (Br*Br - Bz*Bz)/B2
  czz = - crr
  crz = 2.0_fp*  Br*Bz/B2
  czr = crz

  KrNew = crr * yok(4) + crz * yok(5)
  KzNew = czr * yok(4) + czz * yok(5)

  yok(4) = KrNew
  yok(5) = KzNew

  y(1:n) = yok(1:n)

end subroutine lscsq_bounceit
!
!     -----------------------------------------------------------------
!
FUNCTION lscsq_ran3(idum)
!     Transcribed 1993 by D. W. Ignat
!     References:
!     W. H. Press et al, Numerical Recipes in Fortran, p199
!     D. Knuth, Seminumerical Algorithms, Vol 2 of The Art of Computer
!          Programming
!     Returns a uniform random deviate between 0.0 and 1.0.
!     Set idum to any negative value to initialize or reinitialize the
!     sequence.
!     Substitute the CCCommented lines for the ones following to
!     render the routine entirely floating point.
  use iso_c_binding, only : fp => c_double
  implicit none
  real(fp) lscsq_ran3
  INTEGER :: idum, iff, mj,mk
  INTEGER :: i,ii,k, inext, inextp
!     .                                 55 dimension is special; no changes!
  integer :: ma(55)
  integer, parameter :: mbig=1000000000
  integer, parameter :: mseed=161803398
  integer, parameter :: mz=0
  real(fp), parameter :: fac=1.0_fp/real(mbig,kind=fp)
      DATA iff /0/
!                                       Initialize
      if(idum .LT. 0 .or. iff .EQ. 0) then
        iff=1
!     .                                 Initialize ma(55) using the seed idum
!     .                                 and the large number MSEED
        mj=MSEED-abs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
!     .                                 Now initialize the rest of the table
!     .                                 in a slightly random order
!     .                                 with numbers that are not especially
!     .                                 random
        do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk .LT. MZ) mk=mk+MBIG
          mj=ma(ii)
        enddo
! 11     continue
!     .                                 We randomize them by
!     .                                 'warming up the generator'
        do k=1,4
          do i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i) .LT. MZ) ma(i)=ma(i)+MBIG
          enddo
        enddo
!     .                                 Prepare indices for our first
!     .                                 generated number.   The constant 31
!     .                                 is special...see Knuth
        inext=0
        inextp=31
        idum=1
      endif
!
!     .                                 Here is where we start usually
!     .                                 Increment inext, wrap around 56 to 1
      inext=inext+1
      if(inext .EQ. 56)inext=1
      inextp=inext+1
!     .                                 Ditto for inextp
      if(inextp .EQ. 56)inextp=1
!     .                                 Now generate a new random number
      mj=ma(inext)-ma(inextp)
!     .                                 Be sure it is in the range
      if(mj .LT. MZ)mj=mj+MBIG
!     .                                 Store it, and output ran3
      ma(inext)=mj
      lscsq_ran3=mj*FAC
      return
END
!
!     -----------------------------------------------------------------
!
!
!
!                                                                      |
!                                                                      |
!     PredcLSC ends                         ---------------------------|
!
!     eps -- dielectric tensor              ---------------------------|
!                                                                      |
!                                                                      |
SUBROUTINE lscsq_eps ( r, z, kpar2, kper2 )

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : D11er, D33er, D12er
  use lscsq_mod, only : D11ar, D33ar, D12ar, D11w0, D33w0, D12w0
  use lscsq_mod, only : Eper , Epar , Exy  , Aion , Aelc 
  use lscsq_mod, only : ecyc , ecyc2, epsq , ipsq, fghz, iray
  implicit none

  real(fp) :: r,z,kpar2,kper2, veow2,elam,elamc,emli0,emli1, &
              exy0,lscsq_bsi0,lscsq_bsi1, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael

  real(fp) :: large=1.0e30_fp
  real(fp) :: Te2Ve=0.0445e-4_fp
      
  call lscsq_plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
  If (psi.GE.large) return
      ecyc  = omc/fghz(iray)
      ecyc2 = ecyc*ecyc
      epsq  = pe2/fghz(iray)**2
      ipsq  = pi2/fghz(iray)**2
      veow2 = Te2Ve*tee/fghz(iray)**2
      elamc = veow2/ecyc2
      elam  = elamc*Kper2
!     This is the electron Lambda - - (Kper RhoE)^2
!     Now form the quantities exp(-elam) I0,1 (elam)
      emli0= lscsq_bsi0(elam)
      emli1= lscsq_bsi1(elam)
      epar = 1.0_fp - epsq
      aion = aio/fghz(iray)**4
      aelc = ael
      eper = 1.0_fp+epsq/ecyc2-ipsq
      eper = eper-kper2*(aion+aelc)
      exy0 = epsq/ecyc
      exy  = exy0
!      return

!      ENTRY lscsq_epsdr(r, z, Kpar2, Kper2)
!     derivatives of the dielectric tensor elements
!     by (k-perp)**2 ; (k-par)**2 ; (omega)**2
!     this last is multiplied by (omega**2)

      d11er = -(aion + aelc)
      d11ar = 0.0_fp
      d11w0 = ipsq+2.0_fp*aion*kper2

      d12er = 0.0_fp
      d12ar = 0.0_fp
      d12w0 = -0.5_fp*exy

      d33er = 0.0_fp
      d33ar = 0.0_fp
      d33w0 = epsq
!      return

end subroutine lscsq_eps

!     FTION  function for integration       ---------------------------|
subroutine lscsq_ftion(BoundsEr)
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : y, f, pe2min
  use lscsq_mod, only : eper, Epar, woc2, woc4, dDdkABS, Exy
  use lscsq_mod, only : D11ar, D11er, D11w0, D12er, D12ar, D12w0
  use lscsq_mod, only : D33ar, D33er, D33w0
  use lscsq_mod, only : dkpar, dkper, wdDdw, denom
  implicit none

!!!!!!     EXTERNAL plasma2d, eps, epsdr
!     KdB   k \cdot B = ( k_r B_r + k_z B_z + (n/r) B_{\phi} )
!     KdK   k \cdot k = ( k_r^2   + k_z^2   + (n/r)^2 )
!     \p = \partial , the Greek for partial derivative
!     \R = {\bf r}  , the vector r
!     \K = {\bf k}  , the vector k
!     \abs{#1} = \mid #1 \mid , the absolute value
!     s  =  s       , the abs value of \R; ds = path length
!     D  = D(\R,\K,\omega) ; D = 0 is the dispersion relation
!     D  = Kper4 * Eper
!        + Kper2 * [ Qpar*(Epar+Eper) + woc2*Exy**2 ]  ! Kper2 * bb
!        +  Epar * [(Qpar )**2 - (woc2*Exy)**2 ]       ! Epar  * cc
!     bkb(3) = 0.5 \p k_{\parallel}^2 / \p (k_r,k_z,n)
!     per(3) = 0.5 \p k_{\perp}^2     / \p (k_r,k_z,n)
!     dRdwt  = -\frac{\p D/\p\K}{\omega \p D/\p\omega}
!     dKdwt  = +\frac{\p D/\p\R}{\omega \p D/\p\omega}
!     dsdwt  = \abs{dRdwt}
!     dRds   = dRdwt/dsdwt      == f(1,2,3)
!     dKds   = dKdwt/dsdwt      == f(4,5,6)
!     dwtds  = 1./dsdwt
!     wdDdw  = \omega \p D/\p \omega
!     DKpar ==  \p D/\p k_{\parallel}^2
!     DKper ==  \p D/\p k_{\perp}^2
!     Qpar = Kpar2 - woc2*Eper
!     QparE= \p Qpar/\p Kper2
!     QparA= \p Qpar/\p Kpar2
!     QparW= \p Qpar/\p \omega^2  \cdot \omega^2

  integer, intent(out) :: BoundsEr
  real(fp) ::     lscsq_DispRela
  real(fp) ::     KdB,KdK,  Kper2,Kper4,Kpar2
  real(fp) ::     Btot2, Qpar, QparE, QparA, QparW, bb, cc
  real(fp) ::     bkb(3), per(3), dRdwt(3), dKdwt(3),  dsdwt
  real(fp) ::     r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael

  real(fp) :: dr = 1.0e-4_fp      
  real(fp) :: dz = 1.0e-4_fp      

  r = Y(1)
  z = Y(2)
  call lscsq_plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
  BoundsEr = 0
  If (pe2.le.pe2min) then
     BoundsEr = 1
     return
  endif

  Btot2 =  Br*Br + Bz*Bz + (RBphi/r)**2
  KdK = Y(4)**2 + Y(5)**2 + (Y(6)/r)**2

  KdB = Y(4)*Br + Y(5)*Bz + Y(6)/r**2 *RBphi
  Kpar2 = KdB**2/Btot2
  Kper2 = KdK - Kpar2
  CALL lscsq_eps(r,z,Kpar2,Kper2)
!  CALL lscsq_epsdr ( r, z, Kpar2, Kper2 )
  bkb(1) = Br*KdB/Btot2
  bkb(2) = Bz*KdB/Btot2
  bkb(3) = rBphi*KdB/Btot2 /r**2 

  per(1) = Y(4) - bkb(1)
  per(2) = Y(5) - bkb(2)
  per(3) = Y(6)/r**2 - bkb(3)
10    continue
  Kper4 = Kper2**2
  Qpar = Kpar2 - woc2*Eper
  QparE= -woc2*D11er
  QparA= 1.0_fp -woc2*D11ar
  QparW= -woc2*(Eper+D11w0)

  bb = Qpar*(Epar+Eper) + woc2*Exy*Exy
  cc = Qpar*Qpar - woc4*Exy*Exy

  denom = Kper4*D11w0                                               &
           + Kper2*(QparW*(Epar+Eper) + Qpar*(D33w0+D11w0)          &
           +        woc2*Exy*Exy + 2.0_fp*woc2*EXY*D12w0      )        &
           + cc*D33w0                                               &
           + Epar*(2.0_fp*Qpar*QparW -2.0_fp*woc4*Exy*Exy                  &
                  -2.0_fp*woc4*Exy*D12w0                 )
  denom = 2.0_fp*denom
  wdDdw = denom

  DKper = 2.0_fp*Eper*Kper2 + bb                                       &
            +  D11er *Kper4                                         &
            + Kper2*(QparE*(Epar+Eper) + Qpar*(D33er+D11er)         &
                   + woc2*2.0_fp*Exy*D12er                      )      &
             + Epar*(Qpar*QparE - woc4*Exy*D12er)*2.0_fp               &
             +D33er*cc

  DKpar = Kper4*D11ar                                               &
           + Kper2*(QparA*(Epar+Eper) + Qpar*(D33ar+D11ar)          &
           +       2.0_fp*woc2*Exy*D12ar   )                           &
           + D33ar*cc                                               &
           + Epar*2.0_fp*(Qpar*QparA - woc4*Exy*D12ar )

  dRdwt(1) = -2.0_fp* ( DKper*per(1) + DKpar*bkb(1) )/denom
  dRdwt(2) = -2.0_fp* ( DKper*per(2) + DKpar*bkb(2) )/denom
  dRdwt(3) = -2.0_fp* ( DKper*per(3) + DKpar*bkb(3) )/denom

  dsdwt    = sqrt (dRdwt(1)**2 + dRdwt(2)**2 + r*r*dRdwt(3)**2)

  dKdwt(1) =(lscsq_DispRela(Y(1)+DR,Y(2), Y(4),Y(5),Y(6) ) -      &
             lscsq_DispRela(Y(1)-DR,Y(2), Y(4),Y(5),Y(6) ) ) /(2.0_fp*DR)/denom
  dKdwt(2) =(lscsq_DispRela(Y(1),Y(2)+DZ, Y(4),Y(5),Y(6) ) -      &
             lscsq_DispRela(Y(1),Y(2)-DZ, Y(4),Y(5),Y(6) ) ) /(2.0_fp*DZ)/denom
  dKdwt(3) = 0.0_fp
  f(1)     = dRdwt(1)/dsdwt
  f(2)     = dRdwt(2)/dsdwt
  f(3)     = dRdwt(3)/dsdwt
  f(4)     = dKdwt(1)/dsdwt
  f(5)     = dKdwt(2)/dsdwt
  f(6)     = 0.0_fp
  f(7)     = 1.0_fp     /dsdwt
  dDdkABS  = wdDdw  * dsdwt

end subroutine lscsq_ftion
!                                                                      |
!
subroutine lscsq_WhchRoot (r, z, Kr, Kz, Kphi, NperFs, NperSl)

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : woc2, woc4, Epar, Eper, Exy
  implicit none

   real(fp) :: NperFs, NperSl, AAD1, BBD2, CCD4, B2m4AC
   real(fp) :: psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
   real(fp) :: Kpar2, Kper2, KdK
   real(fp) :: Btot2, Qpar
   real(fp), intent(in) :: r, z, Kr, Kz, Kphi

 
   KdK  = Kr**2 + (Kz)**2 + (KPhi/r)**2
   CALL lscsq_plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
   Btot2 =  Br**2 + Bz**2 + (RBphi/r)**2
   Kpar2 = ( Kr*Br + Kz*Bz + Kphi*rBphi/r/r )**2/Btot2
   Kper2 = KdK - Kpar2

   CALL lscsq_eps (r, z, Kpar2, Kper2)
   Qpar = Kpar2 - woc2*Eper
   AAD1 = Eper
   BBD2 = ( (Epar+Eper)*Qpar + woc2*Exy**2 )/woc2
   CCD4 = Epar*( Qpar**2 - woc4*Exy**2 )/woc4
   B2m4AC = BBD2**2 - 4.0_fp*AAD1*CCD4
   NperFs = - BBD2 / (2.0_fp*AAD1)
   NperSl = NperFs
   if (B2m4AC .GT. 0.0_fp) then
      B2m4AC = sqrt(B2m4AC)/(2.0_fp*AAD1)
      NperFs = NperFs - B2m4AC
      NperSl = NperSl + B2m4AC
   endif

   NperFs = max (NperFs, 0.0_fp)
   NperSl = max(NperSl, 0.0_fp)
   NperFs = sqrt(NperFs)
   NperSl = sqrt(NperSl)

end subroutine lscsq_whchroot

!     FUNCTION DispRela                     ---------------------------|
!                                                                      |
FUNCTION lscsq_DispRela (r, z, Kr, Kz, Kphi)

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : woc2, woc4, Epar, Eper, Exy, D1, D2, D4
  implicit none

  real(fp) :: lscsq_DispRela
  real(fp), intent(in) :: r, z, Kr, Kz, Kphi
  real(fp) ::  psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
  real(fp) ::  Kpar2, Kper2, KdK
  real(fp) ::  Btot2, Qpar
!
  KdK  = Kr**2 + (Kz)**2 + (KPhi/r)**2
  CALL lscsq_plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
  Btot2 =  Br**2 + Bz**2 + (RBphi/r)**2
  Kpar2 = ( Kr*Br + Kz*Bz + Kphi*rBphi/r/r )**2/Btot2
  Kper2 = KdK - Kpar2

  CALL lscsq_eps ( r, z, Kpar2, Kper2 )
  Qpar = Kpar2 - woc2*Eper
  D1 = Kper2**2*Eper
  D2 = Kper2*( (Epar+Eper)*Qpar + woc2*Exy**2 )
  D4 = Epar*( Qpar**2 - woc4*Exy**2 )
  lscsq_DispRela  = D1 + D2 + D4

END function lscsq_DispRela
!                                                                      |
!     DispRela ends                         ---------------------------|
!                                                                      |
FUNCTION lscsq_bsi0(x)

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : lstop
  implicit none

  real(fp) :: lscsq_bsi0
  real(fp), intent(in) :: x
  real(fp) :: xx
 
  real(fp), parameter :: prettysmall=0.62e-4_fp
  real(fp), parameter :: large=20.0_fp
  real(fp), parameter :: small=1.0e-5_fp

!     BsI0(x) = e^{-x} * I_0 (x)

  if ( x .LT. large ) then
     lscsq_bsi0 = 1.0_fp
     if ( x .LT. small ) return
     xx = 0.5_fp*x
     xx = xx*xx
     lscsq_bsi0 = exp(-x)*(1.0_fp+xx+0.25_fp*xx**2+xx**3/36.0_fp)
     return
  else
     lscsq_bsi0 = PRETTYSMALL
     Lstop = 1
     CALL lscsq_LSCwarn (' lambda too large in function lscsq_bsi0')
  endif

end function lscsq_bsi0
!                                                                      
!                                                                      
FUNCTION lscsq_bsi1(x)
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : lstop
  implicit none
!
!     BsI1(x) = e^{-x} * I_1 (x)
!
  real(fp) :: lscsq_bsi1
  real(fp), intent(in) :: x
  real(fp) :: xx

  real(fp), parameter :: prettysmall=0.16e-3_fp
  real(fp), parameter :: large=20.0_fp
  real(fp), parameter :: small=1.0e-5_fp

  if ( x .LT. large ) then
     xx = 0.5_fp*x
     lscsq_bsi1 = exp(-x)*abs(xx)
     if ( x .LT. small) then
        return
     else
        xx = xx*xx
        lscsq_bsi1=lscsq_bsi1*(1.0_fp+0.5_fp*xx+xx**2/12.0_fp+xx**3/144.0_fp)
        return
     endif
  else
     lscsq_bsi1 = PRETTYSMALL
     Lstop = 1
     CALL lscsq_LSCwarn (' lambda too large in function lscsq_bsi1')
  endif

end function lscsq_bsi1
!

module RayIni
contains

  SUBROUTINE lscsq_RayIni(RayIniErr, y1, y2, y3, y4, y5, y6, y7, y8,&
    iray, arrys, o_data, p_data)!, local_prmtrs)!, all_vecs)
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
    use pl_types 
    use Tracing
    use rTrace
    use iso_c_binding, only : fp => c_double
    use lscsq_mod, only : twopi, vc, deg2rad
    use lscsq_mod, only: ntor, npol, thgrid, ind_ray
    use lscsq_mod, only: lfast, begin
    use lscsq_mod, only: rmax,rmin, zmax, rmaj
    use lscsq_mod, only: psimin,psilim
    !use lscsq_mod, only : npsij
    use lscsq_mod, only: pe2fac, SEARCHINCR
    use lscsq_mod, only: delpsi, nzones, nrays

    implicit none
    type(old_data), intent(inout) :: o_data
    type(pl_data), intent(inout) :: p_data
    type(storeAry), dimension(nzones), intent(inout) :: arrys
    !type(parameters), intent(inout) :: local_prmtrs

    Real(fp), intent(inout):: y1, y2, y3, y4, y5, y6, y7, y8
    !Real(fp), dimension(neqs) ::f_dumy

    integer, intent(inout) :: RayIniErr
    integer, intent(in) :: iray
    Real(fp) :: max_d, enpar, enpol, enth
    Real(fp) :: Kpar2, Kpol, Krad, Ktor
    complex(fp), dimension(3) :: zc
    Real(fp), dimension(4) :: Azplr
    Real(fp) :: r,z,Br,Bz,RBphi,omc,Tee,pi2,&
                aio,aion,ael, ipsq, psi, pe2
    integer ::  i, e_dumy
    Real(fp) :: Bpol,Bpol2,Btot2, &
                det, fast,plas,Qpar,Sig,slow,try, &
                zi1,zi2,zi3,zr1,zr2,zr3, discrimt, Kper
    Real(fp) :: CosT,drad,rad,SinT,t 
    character(len=70) :: ErrMsg

    !type(vecs),dimension(npsij),  intent(inout):: all_vecs

    !acc routine
    !acc routine(lscsq_plasma2d)
    !$cc routine(lscsq_LSCwarn)
    !acc routine(lscsq_disprela)
    !acc routine(lscsq_ZPLRCnr)
    !acc routine(lscsq_lscstop)

    

    enpar = ntor(ind_ray(1, iray))
    if (enpar==0.0_fp) then
      RayIniErr = 2
      write(*,*) ' enpar == 0 in RayIni for Ray#', iray
      !CALL lscsq_LSCwarn(' enpar == 0 in RayIni ')
      return
    endif
    
    enpol = npol(ind_ray(2,iray))
    enth  = thgrid(ind_ray(3,iray))

    Kpol = 1.0e9_fp*enpol*twopi*p_data%fghz/vc
    Ktor = 1.0e9_fp*enpar*twopi*p_data%fghz/vc
    T = enth * deg2rad  

    Rad = Rmax-Rmaj
    If ( Cos(T) < 0.0_fp ) Rad=Rmaj-Rmin
    CosT = Abs(Cos(T))+1.0e-20_fp
    SinT = Abs(Sin(T))+1.0e-20_fp
    Rad  = min(Rad/CosT, Zmax/SinT)
    drad = SEARCHINCR
  
    
    do I=1,10000
      
      Rad = Rad - drad
      r    = Rmaj + Rad*Cos(T)
      z    =        Rad*Sin(T)

      call lscsq_plasma2d(o_data, p_data, r,z,psi,Br,&
                          Bz,RBphi,omc,&
                          Tee,pe2,pi2,aio,ael)!,all_vecs)


      If (pe2 <= p_data%pe2min) cycle   

      Bpol2= Br*Br + Bz*Bz
      Btot2= Bpol2 + (RBphi/r)**2
      Bpol = sqrt(Bpol2)
      Sig  = (r-Rmaj)*Bz - z*Br
      If (Sig < 0.0_fp) Bpol = -Bpol
      Kpar2 = (Kpol*Bpol+Ktor*RBphi/R)**2/Btot2
      
      call lscsq_Eps(0.0_fp, psi, omc, pe2, pi2, tee, &
                    aio, aion, ael, ipsq, p_data)

      Azplr(4) = -(aio/p_data%fghz**4 + ael)
      Azplr(3) =  p_data%Eper
      Qpar = Kpar2 - p_data%woc2*p_data%Eper
      Azplr(2) = Qpar*(p_data%Epar+p_data%Eper) + p_data%woc2*p_data%Exy**2
      Azplr(1) = p_data%Epar*( Qpar**2 - p_data%woc4*p_data%Exy**2 )

      ! Estimate the roots far from lmc
      discrimt = Azplr(2)**2-4.0_fp*Azplr(3)*Azplr(1)
      if ( discrimt <= 0.0_fp ) then
          discrimt = 0.0_fp
          RayIniErr= 2
          !write(*,*) ' no accessibity found for this ray '!, iray
          !ErrMsg= ' no accessibity found for this ray ', iray
          !CALL lscsq_LSCwarn(ErrMsg)
          return
      endif

      CALL lscsq_ZPLRCnr(RayIniErr, 3,Azplr,ZC)
      if (RayIniErr >= 1) then
        write(0,*) ' too many iterations in LAGUER root finder for Ray# '!, iray
        return
      endif
  !     Finds Zeroes of a Polynomial with Laguerre's method if
  !     Real Coefficients
  !     using Numerical Recipes code so we dont depend on IMSL
        zr1 = REAL(zc(1),kind=fp)
        zr2 = REAL(zc(2),kind=fp)
        zr3 = REAL(zc(3),kind=fp)
        zi1 = Real(AIMAG(zc(1)), fp)
        zi2 = Real(AIMAG(zc(2)), fp)
        zi3 = Real(AIMAG(zc(3)), fp)

  !     The normal case is that all roots are real, and we see if we can start.

        plas =  max (zr1,zr2,zr3)
        fast =  min (zr1,zr2,zr3)
        slow = -1.0e30_fp 

        if ((zr1 < plas) .and. (zr1 > fast)) slow = zr1
        if ((zr2 < plas)  .and. (zr2 > fast)) slow = zr2
        if ((zr3 < plas) .and. (zr3 > fast)) slow = zr3

        if (slow == -1.0e+30_fp) cycle    
        if (lfast == 1) then
          try = fast
        !if (lfast == 0) 
        else
          try = slow
        endif
        if (try < 0.0_fp ) cycle    
        try = try - ( Kpol*RBphi/R - Ktor*Bpol )**2/Btot2
        if ( try < 0.0_fp ) cycle    
        try = sqrt(try)

        krad = try
        if(lfast == 1) krad = - try

        go to 35
    enddo
    RayIniErr =2
    write(0,*)  ' cant find a starting point for this ray'!, iray
    !CALL lscsq_LSCwarn(' cant find a starting point for this ray')
    return
  35  continue

    y1  =  r
    y2  =  z
    y3  =  0.0_fp
    y4  =  ( Krad*Bz + Kpol*Br ) / Bpol
    y5  =  (-Krad*Br + Kpol*Bz ) / Bpol
    y6  =  Ktor*r
    y7  =  0.0_fp
    ! to begin the time as fn of distance
    y8 = begin
    det =  lscsq_DispRela(y1,y2,y4,y5,y6, o_data, &
                          p_data, max_d)

    e_dumy = 0
    !call lscsq_ftion(e_dumy, y, f_dumy, o_data, p_data)
    
  !     only to initialize dkper & wdDdw for damping calculation

    Kper=sqrt(abs(y4**2+y5**2+(y6/y1)**2-Kpar2))

    o_data%RzindOld = (Psi-PsiMin)/DelPsi + 1.5_fp
    o_data%IzindOld = int(o_data%RzindOld)

    arrys(1)%izind = o_data%IzindOld
    arrys(1)%PowrRy   = 1.0_fp
    arrys(1)%RofRay   = y1
    arrys(1)%ZofRay   = y2
    arrys(1)%PofRay  = y3
    arrys(1)%NparRy   = sqrt(Kpar2)/sqrt(p_data%woc2)
    arrys(1)%NperRy   = Kper/sqrt(p_data%woc2)
    arrys(1)%rtPsRy   = sqrt((psi-psimin)/(psilim-psimin))
    arrys(1)%TimeRy   = y7
    arrys(1)%DistRy  = y8
    arrys(1)%NeofRy   = pe2/Pe2Fac*1.0e+14_fp
    arrys(1)%BthRay   = sqrt(Br**2+Bz**2)
    arrys(1)%DetrRy   = det/max_d 

    ! The quantities to be collected (RayQt) are averaged over the
    ! zone by summing into accum, dividing by NinAc. This clears.

  end subroutine lscsq_rayini
  !
  !     ------------------------------------------------------------------
  !
  SUBROUTINE lscsq_zplrcnr (RayIniErr, degree, coefr, zeroc)

    use iso_c_binding, only : fp => c_double
    implicit none

  !     finds Zeros of Polynomials with Laguerre's method assuming
  !     Real Coefficients by calling a routine from
  !     Numerical Recipes.
  !     The calling convention is set up to look like IMSL CALL lscsq_ZPLRC.
  !     f(x) = 0 = coefr(1) + coefr(2)*x + coefr(3)*x**2 ...etc
    integer, intent(in) :: degree
    integer :: i, polish
    integer, parameter :: degmax=10
    integer, intent(inout) :: RayIniErr
    Real(fp), intent(inout), dimension(degree+1) :: coefr
    complex(fp), dimension(degmax+1) :: coefc
    complex(fp), intent(inout), dimension(degree) :: zeroc

    !acc routine
    !acc routine(lscsq_zrootsnr)

      ! zroots needs complex coef's
    do i = 1, degree+1
      coefc(i) = cmplx ( coefr(i) , 0.0_fp, kind=fp)
   enddo
  
   polish = 1
   CALL lscsq_zrootsnr ( RayIniErr, coefc, degree, zeroc, polish )
  
 end subroutine lscsq_zplrcnr

  !
  !     ------------------------------------------------------------------
  !

 SUBROUTINE lscsq_laguernr ( RayIniErr, coef, degree, x, epsilon, polish )
  use Tracing
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
  integer, intent(inout) :: RayIniErr
  integer, intent(in) :: degree
  integer :: iter, j, polish
  integer :: maxit!=100
  complex(fp) :: czero!=(0.0_fp,0.0_fp)

  complex(fp), intent(in) :: coef(degree+1)
  complex(fp), intent(inout) :: x
  complex(fp) :: dx, x1,b,d,f,g,h,sq,gp,gm,g2
  Real(fp) :: epss! = 6.0e-7_fp
  Real(fp), intent(in) :: epsilon
  Real(fp) :: err, abx, cdx
  !character(len=70) :: ErrMsg

  !acc routine
  !acc routine(lscsq_lscstop)

  maxit=100
  epss = 6.0e-7_fp
  czero=(0.0_fp,0.0_fp)

  do iter = 1, MAXIT
    b   = coef ( degree+1 )
    err =  abs ( b )
    d   = CZERO
    f   = CZERO
    abx =  abs ( x )
    ! Efficient computation of the
    ! polynomial and its first 2
    ! derivatives
    do j = degree, 1, -1
        f   = x*f + d
        d   = x*d + b
        b   = x*b + coef ( j )
        err =  abs ( b ) + abx * err
    enddo
    ! Estimate of roundoff error in
    ! evaluating polynomial
    err = EPSS * err
    if (  abs ( b ) <= err ) then
        return
    ! The generic case: use Laguerre's formula
    else
        g  = d / b
        g2 = g * g
        h  = g2 - 2.0_fp * f / b
        sq =  sqrt ( (degree-1) * (degree*h - g2) )
        gp = g + sq
        gm = g - sq
        if (  abs (gp) <  abs (gm) ) gp = gm
        dx = degree / gp
    endif
    x1    = x - dx
    if ( x == x1 ) return
    x     = x1
    cdx   =  abs ( dx )
    if ( polish == 0 ) then
        if ( cdx <= epsilon* abs ( x ) ) return
    endif
  enddo     

  !CALL lscsq_LSCstop( ' too many iterations in LAGUER root finder ')
  !ErrMsg= ' too many iterations in LAGUER root finder '
  !CALL lscsq_LSCwarn(ErrMsg)
  RayIniErr = 2
end subroutine lscsq_laguernr 
  !
  !     ------------------------------------------------------------------
  !

  subroutine lscsq_zrootsnr ( RayIniErr, coef, degree, roots, polish )
    use Tracing
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
    integer, intent(inout) :: RayIniErr
    INTEGER degree, i, j, jj, polish
    integer, parameter :: maxdegre=101
    Real(fp) :: epsilon! = 1.0e-6_fp
    complex(fp) :: coef(degree+1), x, b, c, roots(degree), defl(MAXDEGRE)

  ! Copy coef's for successive deflation
    !acc routine
    !acc routine(lscsq_laguernr)

    epsilon = 1.0e-6_fp


    defl(1:degree+1) = coef(1:degree+1)
!
  do j = degree, 1, -1
     ! Start at 0 to favor smallest remaining root
     x = cmplx ( 0.0_fp, 0.0_fp, kind=fp)
     CALL lscsq_laguernr ( RayIniErr, defl, j, x, EPSILON, 0 )
     if (abs(Real(AIMAG(x),fp)).LE.2.*EPSILON**2*abs(REAL(x,kind=fp)))       &
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
        CALL lscsq_laguernr ( RayIniErr, coef, degree, roots(j), EPSILON, 1 )
     enddo
  endif
 
  ! Sort roots by real part by straight insertion
  do j = 2, degree
     x = roots (j)
     do i = j-1, 1, -1
        if ( REAL( roots(i),kind=fp ) <= REAL( x,kind=fp ) ) go to 10
        roots ( i+1 ) = roots ( i)
     enddo
     i=0
10    roots ( i+1 ) = x
  enddo        

end subroutine lscsq_zrootsnr     


  SUBROUTINE lscsq_setfreq(indx, p_data)
    !  (DMC Mar 2011: reset frequency; can be different for each ray, now)
      use pl_types
      use iso_c_binding, only : fp => c_double
      use lscsq_mod, only : pi, twopi, vc
    ! use lscsq_mod, only : fghz_now, omega, woc2, woc4
      use lscsq_mod, only: ngrps, fghz, powers
      implicit none

      type(pl_data), intent(inout) :: p_data
      integer, intent(in) :: indx  ! frequency index [1:ngrps] or 0 or -1
    !
    !        if indx>0 set the frequency for the indicated group
    !        if indx==0 set frequency for group with maximum power
    !        if indx<0 set the frequency to zero
    !
    !        also set associated variables
    !
    !------------------------
      real(fp) :: pmax
      integer :: indx_maxp,igrp
    !------------------------
    
      IF((indx<0).OR.(indx>nGrps)) then
        p_data%fghz = 0.0
      ELSE IF(indx==0) then
        pmax=1
        indx_maxp=1
        DO igrp=2,nGrps
            if(powers(igrp)>pmax) then
              pmax=powers(igrp)
              indx_maxp=igrp
            endif
        ENDDO
        p_data%fghz = fghz(indx_maxp)
      ELSE
        p_data%fghz = fghz(indx)
      ENDIF

      p_data%omega =  1.0e09_fp * twopi*p_data%fghz 
      p_data%woc2 = (p_data%omega/vc)**2
      p_data%woc4 =  p_data%woc2**2

    
    end subroutine lscsq_setfreq


      subroutine lscsq_RayDataInit(arrys)
        use pl_types 
        use iso_c_binding, only : fp => c_double
        use lscsq_mod, only: npar, nzones, nrays

        implicit none
       
        integer :: iray, j
        type(storeAry), dimension(nzones, nrays), intent(inout) :: arrys

        !acc routine seq

        !acc parallel loop collapse(2)
        do iray = 1, nrays
          do j=1,nzones
            !write(*,*) j
            arrys(j,iray)%RofRay = 0.0_fp
            arrys(j,iray)%ZofRay = 0.0_fp
            arrys(j,iray)%PofRay = 0.0_fp
            arrys(j,iray)%nparry = 0.0_fp
            arrys(j,iray)%nperry = 0.0_fp
            arrys(j,iray)%rtpsry = 0.0_fp
            arrys(j,iray)%powrry = 0.0_fp
            arrys(j,iray)%timery = 0.0_fp
            arrys(j,iray)%distry = 0.0_fp
            arrys(j,iray)%detrry = 0.0_fp
            arrys(j,iray)%ivind = 0
            arrys(j,iray)%ezsq = 0.0_fp
            arrys(j,iray)%dlnPdsk  = 0.0_fp
            arrys(j,iray)%dlnPdsx  = 0.0_fp
          enddo
        enddo

        !acc parallel loop 
        do iray = 1, nrays
          arrys(1,iray)%NparRy = npar(1,iray)
        enddo
       
      end subroutine lscsq_RayDataInit


      subroutine lscsq_PlDataInit(iray, &
        p_data, o_data, a_data)
        use pl_types 
        use iso_c_binding, only : fp => c_double
        use lscsq_mod, only: npar
 
        use lscsq_mod, only : twopi, vc!, nzones, nrays 
        use lscsq_mod, only : fghz, pe2Vec, NpsiJ
        
        implicit none
  
        integer, intent(in) :: iray
        type(pl_data), intent(inout) :: p_data
        type(accum_data), intent(inout) :: a_data
        type(old_data), intent(inout) :: o_data
        
        !acc routine seq

        
        p_data%pe2min = pe2Vec(NpsiJ)
        p_data%fghz = fghz(iray)
        p_data%omega =  1.0e09_fp * twopi*fghz(iray) 
        p_data%woc2 = (p_data%omega/vc)**2
        p_data%woc4 =  p_data%woc2**2
        p_data%woc = sqrt(p_data%woc2)
        p_data%dtdv = 0.0_fp
        p_data%epsq = 0.0_fp
        p_data%ecyc2 = 0.0_fp
        !p_data%ecyc = 0.0_fp
    
        o_data%first_call = .true.
        o_data%iread = 0
        o_data%NinAc = 0
    
        a_data%accum1 = 0.0_fp
        a_data%rayqt1 = 0.0_fp
        a_data%accum2 = 0.0_fp
        a_data%rayqt2 = 0.0_fp
        a_data%accum2 = 0.0_fp
        a_data%rayqt3 = 0.0_fp
        a_data%accum4 = 0.0_fp
        a_data%rayqt4 = 0.0_fp

      end subroutine lscsq_PlDataInit

      
end module rayIni

  
  !
  !     ------------------------------------------------------------------
  !

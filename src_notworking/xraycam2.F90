! -*-f90-*-
!23456789-123456789-123456789-123456789-123456789-123456789-123456789-23
! import_py.inf: TRANSP import script
!     ------------------------------------------------------------------
!     Authors:
!     E. J. Valeo and  D. P. Enright, August 1992
!     Modified for loading with LSC by D. W. Ignat, October 1992
!     Modified August 24 1993 by Ignat to get w/cm^2; also fix some
!     indexing errors.
!     Copyright:
!     Princeton University, Plasma Physics Laboratory, 1992, 1993, by
!     E. J. Valeo, D. W. Ignat, S. von Goeler, S. C. Jardin,
!     P. G. Roney, J. E. Stevens, and D. P. Enright
!
!     Units: Energy in MeV mostly; kev in absorber routines;
!     mc^2 in the evaluation of Koch and Motz formula.
!     Bremsstrahlung photon radiation in w/cm^2 on camera; w/cm^3 in plasma
!
!     References:
!     H. W. Koch and J. W. Motz,
!     "Bremsstrahlung Cross-Section Formulas and Related Data,"
!     \RMP{31}{920}{59}.
!     Formula 2BN on p 924
!     '2' denotes differential in photon energy and angle
!     'B' denotes Born approximation
!     'N' denotes no screening
!     and
!     X ray crosss sections based on W. H. McMaster, el al, UCRL 50174
!     Sec. II rev. 1,  May 1969.
!
!     .         --------------------------------------------------------
!
#ifndef USE_X_RAY_CAMERA
      SUBROUTINE lsc_XrayCam2
      return
      END
      SUBROUTINE lsc_XRopen(i,a)
      INTEGER i
      CHARACTER*24 a
      return
      ENTRY lsc_XRclos(i)
      return
      END
      SUBROUTINE lsc_XRdscw(x, y, n, a)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER n
      REAL    x(n), y(n)
      CHARACTER*(*) a
      return
      END
!      SUBROUTINE XRdscr
!      return
!      END
#else
      SUBROUTINE lsc_XrayCam2
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_xparams.h"
      CALL lsc_xRayInp
      CALL lsc_xBremSetup(ZbrAry(1))
      CALL lsc_xRaySetup
      CALL lsc_xIntens(NVELDIM,NPSIDIM,iITR,fe, nv,                     &
     &             npsi, NeAry, ZbrAry, Vpar)
      CALL lsc_xEmit
!
      return
      END
!
!     ------------------------------------------------------------------
      SUBROUTINE lsc_xRayInp
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none
!                                                                      |
      NAMELIST /inpxry/                                                 &
     &        dFoilTCM, FoilCode, iAbsXray,                             &
     &        nEbins, nMUbins, E_max, E_min, nr_source, nz_source,      &
     &        Z_bound_min, Z_bound_max, R_bound_max, R_bound_min,       &
     &        n_pixel_x, n_pixel_y,                                     &
     &        E_ph_min, dE_ph,                                          &
     &        Rpinhole, Zpinhole, phi_pinhole, pinhole_size,            &
     &        focal_length, screen_d, R_tangent, z_tangent,             &
     &        dlxray, npoints_max,                                      &
     &        PusherMajor, PusherMinor
!                                                                      |
!                                                                      |
!     NameList ends                         ---------------------------|
!
!#include "emitter.inc"
!
      Z_plasm_max =  max (ZlcfsMax,-ZlcfsMin)
      Z_plasm_min = -Z_plasm_max
      R_plasm_max = RlcfsMax
      R_plasm_min = RlcfsMin

      Z_bound_max = Z_plasm_max*1.005d0
      Z_bound_min =-Z_plasm_max*1.005d0
      R_bound_max = RlcfsMax*1.005
      R_bound_min =  min (RlcfsMin/1.005d0,R_bound_min)
!
!#include "camera.inc"
!
!
!#include "numerics.inc"
      open (nTSCunus,status='old',file = 'input.xry', err=1300 )
      read (nTSCunus, inpxry)
      close (nTSCunus)
!
      goto 1301
 1300 CALL lsc_LSCwarn (' cant find input.xry; use defaults! ')
 1301 continue
      if (n_pixel_x .GT. NPIXDIM) then
        CALL lsc_LSCwarn('n_pixel_x reset to NPIXDIM')
        n_pixel_x = NPIXDIM
        endif

      if (n_pixel_y .GT. NPIXDIM) then
        CALL lsc_LSCwarn('n_pixel_y reset to NPIXDIM')
        n_pixel_y = NPIXDIM
        endif

      if (nr_source .GT. NRDIM) then
        CALL lsc_LSCwarn('nr_source reset to NRDIM')
        nr_source = NRDIM
        endif

      if (nz_source .GT. NZDIM) then
        CALL lsc_LSCwarn('nz_source reset to NZDIM')
        nz_source = NZDIM
        endif

      if (nMUbins .GT. NMUDIM) then
        CALL lsc_LSCwarn('nMUbins reset to NMUDIM')
        nMUbins = NMUDIM
        endif

      if (2*(nMUbins/2) .EQ. nMUbins) then
        CALL lsc_LSCwarn('nMUbins must be odd; subtracting 1')
        nMUbins = nMUbins-1
        endif

      if (nEbins .GT. NENDIM) then
        CALL lsc_LSCwarn('nEbins reset to NENDIM')
        nEbins = NENDIM
        endif

      if (E_max .GT. 0.255) then
        CALL lsc_LSCwarn('E_max reset to 0.255 MeV')
        E_max = 0.255
        endif
!
      pix_size_x = screen_d/n_pixel_x
      pix_size_y = screen_d/n_pixel_y
!
      CALL lsc_GrafFixP('BEGIN PLOT', 0    )
      CALL lsc_GrafFixP(' NRDIM    ', NRDIM)
      CALL lsc_GrafFixP(' NZDIM    ', NZDIM)
      CALL lsc_GrafFixP(' NRZMXDIM ', NRZMXDIM)
      CALL lsc_GrafFixP(' NPIXDIM  ', NPIXDIM)
      CALL lsc_GrafFixP(' NRCHORDIM', NCHORDIM)
      CALL lsc_GrafFixP(' MAXPOINTS', MAXPOINTS)
      CALL lsc_GrafFixP(' SPACEDIM ', SPACEDIM)
      CALL lsc_GrafFixP(' NMUDIM   ', NMUDIM)
      CALL lsc_GrafFixP(' NENDIM   ', NENDIM)
      CALL lsc_GrafFixP(' iAbsXray ', iAbsXray)
      CALL lsc_GrafChrS(' FoilCode ', FoilCode)
      CALL lsc_GrafFltP(' dFoilTCM ', dFoilTCM)
      CALL lsc_GrafFixP(' nMUbins  ', nMUbins)
      CALL lsc_GrafFixP(' nEbins   ', nEbins)
      CALL lsc_GrafFltP(' mu_min   ', mu_min)
      CALL lsc_GrafFltP(' mu_max   ', mu_max)
      CALL lsc_GrafFltP(' E_min    ', E_min)
      CALL lsc_GrafFltP(' E_max    ', E_max)
      CALL lsc_GrafFltP(' E_ph_min ', E_ph_min)
      CALL lsc_GrafFltP(' dE_ph    ', dE_ph)
      CALL lsc_GrafFixP(' nr_source', nr_source)
      CALL lsc_GrafFixP(' nz_source', nz_source)
      CALL lsc_GrafFltP(' mu_0     ', mu_0)
      CALL lsc_GrafFltP(' mu_width ', mu_width)
      CALL lsc_GrafFltP(' ZboundMin', Z_bound_min)
      CALL lsc_GrafFltP(' ZboundMax', Z_bound_max)
      CALL lsc_GrafFltP(' RboundMin', R_bound_min)
      CALL lsc_GrafFltP(' RboundMax', R_bound_max)
      CALL lsc_GrafFltP(' ZplasmMin', Z_plasm_min)
      CALL lsc_GrafFltP(' ZplasmMax', Z_plasm_max)
      CALL lsc_GrafFltP(' RplasmMin', R_plasm_min)
      CALL lsc_GrafFltP(' RplasmMax', R_plasm_max)
      CALL lsc_GrafFltP(' ZlcfsMin ', ZlcfsMin)
      CALL lsc_GrafFltP(' ZlcfsMax ', ZlcfsMax)
      CALL lsc_GrafFltP(' RlcfsMin ', RlcfsMin)
      CALL lsc_GrafFltP(' RlcfsMax ', RlcfsMax)
      CALL lsc_GrafFltP(' PusherMaj', PusherMajor)
      CALL lsc_GrafFltP(' PusherMin', PusherMinor)
      CALL lsc_GrafFixP(' n_pixel_x', n_pixel_x)
      CALL lsc_GrafFixP(' n_pixel_y', n_pixel_y)
      CALL lsc_GrafFltP(' focal_len', focal_length)
      CALL lsc_GrafFltP(' screen_d ', screen_d)
      CALL lsc_GrafFltP(' pinhole_s', pinhole_size)
      CALL lsc_GrafFltP(' R_tangent', R_tangent)
      CALL lsc_GrafFltP(' Z_tangent', Z_tangent)
      CALL lsc_GrafFltP(' RcntChrd ', RcntChrd)
      CALL lsc_GrafFltP(' Rpinhole ', Rpinhole)
      CALL lsc_GrafFltP(' Zpinhole ', Zpinhole)
      CALL lsc_GrafFltP(' pix_sizex', pix_size_x)
      CALL lsc_GrafFltP(' pix_sizey', pix_size_y)
      CALL lsc_GrafFixP(' npointsMx', npoints_max)
      CALL lsc_GrafFltP(' dlxray   ', dlxray)
!     CALL lsc_GrafFixP(' .', .)
!     CALL lsc_GrafFltP(' .', .)
!     CALL lsc_GrafChrS(' .', .)
      CALL lsc_MkGrfLst(' List X ray Camera parameters ')
      CALL lsc_GrafFixP('END PLOT  ', 0   )
!
      return
      END
!     ..................................................................
      BLOCK DATA xRayBk
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_xparams.h"
!#include "lsc_emparams.h"
!#include "lsc_xray.h"

!#include "lsc_emitter.h"

!#include "lsc_camera.h"

!#include "lsc_numerics.h"

      DATA mu_min, mu_max / -1.d0, 1.d0 /
      DATA nMUbins / NMUDIM /
      DATA nEbins  / NENDIM /
      DATA E_ph_min / 1.0d-6 /
      DATA dE_ph / .01d0 /
      DATA E_min, E_max / 0.01d0 , 0.2d0 /
      DATA   Z_bound_min ,  Z_bound_max  /                              &
     &     -1.00d0, 1.00d0 /
      DATA  R_bound_min ,  R_bound_max  /                               &
     &     1.10d0, 2.00d0 /
      DATA mu_0 / 1.d0 / mu_width / .2d0 /
      DATA nr_source / NRDIM / nz_source / NZDIM /

      DATA Rpinhole / 2.661d0 /
      DATA Zpinhole /   0.d0 /
      DATA phi_pinhole / 0.d0 /
      DATA pinhole_size / 0.00625d0 /
      DATA focal_length, screen_d  /                                    &
     &          0.4905d0, 0.215d0     /
      DATA n_pixel_x, n_pixel_y  / NPIXDIM, NPIXDIM /
      DATA R_tangent / 1.524d0 /
      DATA RcntChrd  / 1.65d0  /
      DATA z_tangent / 0.d0 /
      DATA pix_size_x, pix_size_y / 0.01d0, 0.01d0 /

      DATA npoints_max / MAXPOINTS /
!     DATA eps / ACCEPTABLE_ERROR /
      DATA eps / 1.0d-3    /
      DATA dlxray / 0.005d0 /
      DATA  PusherMajor, PusherMinor / 1.152d0, 0.190d0 /
      DATA FoilCode / 'AG' /
      DATA dFoilTCM / 0.000d0 /
      DATA iAbsXray / 1 /
!
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_xRaySetup
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_xparams.h"
!#include "lsc_xray.h"
!#include "lsc_emparams.h"
!#include "lsc_emitter.h"
!     generate grids, cross-section database
!     and adjust the sigtot by the transmission and illumination factor
!
      CALL lsc_ugrid(muvec, nMUbins, mu_min, mu_max)
      dmu_inv = AREAL(nMUbins - 1) / (mu_max - mu_min)
!
      CALL lsc_ugrid(E_incvec, nEbins, E_min, E_max)
      CALL lsc_xSetupEindex
!     set up indices for table lookup into sigma
      CALL lsc_xGenEofv(nv, Eofv, dEofv, Vpar)
      if (iAbsXray .GT. 0) then
        CALL lsc_GetXmnFac(nEbins, E_incvec, XmnFac, dFoilTCM, FoilCode)
!        CALL lsc_lscpause
!        CALL lsc_vwrite(nEbins, E_incvec, 'E_incvec')
!        CALL lsc_vwrite(nEbins, XmnFac, 'XmnFac')
!        CALL lsc_lscpause
      endif
!     compute vector XmnFac of absorption coefficients
      CALL lsc_xSigma

      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_xGenEofv(n, E, dE, v)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER n, i
      real   signorm, Mevpermcsq
      COMMON / sigbrcom0 / signorm, Mevpermcsq
      REAL   E(n), dE(n), v(n)
      do i = 1, n
         E(i) = 0.5 * Mevpermcsq * v(i) * v(i)
      enddo
      do i = 2, n - 1
         dE(i) = 0.5 * abs(E(i + 1) - E(i - 1))
!     accounts for negative velocity
      enddo
      dE(1) = 0.5 * abs(E(2) - E(1))
      dE(n) = 0.5 * abs(E(n) - E(n - 1))
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_xSetupEindex
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none


!#include "lsc_xparams.h"
!#include "lsc_xray.h"
      INTEGER i, j
      REAL   Einc, MevPerMcsq
      PARAMETER ( MevPerMcsq = 0.511d0 )

      do i = 1, nv
         Einc = 0.5 * Mevpermcsq * Vpar(i) * Vpar(i)
         j = 1
         do while((E_incvec(j) .LE. Einc) .and. (j .LT. nEbins))
            j = j + 1
         enddo
         j = j - 1
         Eint(i) = j
         if(j .GE. 1 .and. j .LT. nEbins) then
           Efrac(i) = (Einc - E_incvec(j)) / (E_incvec(j + 1) -         &
     &        E_incvec(j))
         else
           Efrac(i) = 0.00d0
         endif
      enddo
      return
      END
!     .         --------------------------------------------------------
      REAL FUNCTION lsc_Xraytransmission(EmeV)
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none


!#include "lsc_xparams.h"
!#include "lsc_xray.h"
      INTEGER ie
      REAL frac, EmeV
      if(iAbsXray .GT. 0)then
         CALL lsc_Lookupindex(EmeV, E_incvec, nEbins, ie, frac)
         lsc_Xraytransmission = XmnFac(ie) * (1. - frac) +              &
     &        XmnFac(ie + 1) * frac
      else
         lsc_Xraytransmission = 1.d0
      endif
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_Lookupindex(c, vec, n, i, frac)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER n, i
      REAL                                                              &
     &        c, vec(n), frac, ri
      real    RE41
      ri = (c - vec(1)) / (vec(2) - vec(1))
!     i = ifix(ri)
      RE41 = ri
      i = int(RE41)
      frac = ri - AREAL(i)
      if(i .LT. 1)then
         i = 1
         frac = 0.
      else if(i .GE. n)then
         i = n - 1
         frac = 1.d0
      endif
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_xSigma
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none


!#include "lsc_xparams.h"
!#include "lsc_xray.h"
      REAL    sigg, E_photon, lsc_Xraytransmission
      INTEGER  jj, kk
      do 10 kk = 1, nEbins

         do 20 jj = 1, nMUbins
            sigtot(jj,kk) = 0.d0
            E_photon = E_ph_min + 0.5 * dE_ph
            do while(E_photon .LT. E_incvec(kk))
               CALL lsc_xSigGen(E_incvec(kk), E_photon, muvec(jj),      &
     &              sigg)
               sigtot(jj,kk) = sigtot(jj,kk) + sigg * E_photon * dE_ph  &
     &              * lsc_Xraytransmission(E_photon       )
!     weighted to yield rate of energy emission.
!
               E_photon = E_photon + dE_ph
            enddo
 20      continue
 10   continue

!cc    Ernies contour plot of sigtot contours 28Mar93 begin
!     CALL lsc_X_pl_intensity_cntrs(sigtot, NMUDIM, nMUbins, NEbins)
!     CALL lsc_EZfini(0,0)
!cc    Ernies contour plot of sigtot contours 28Mar93 begin

      return
      END

!     +---------------------------------------------------------------+

      SUBROUTINE lsc_xIntens(NVDIM,NPDIM,iIterate,FeOfv, nv, np,        &
     &     NeAry, ZbrAry, Vpar)
!     intensity subroutine calculates the value of the photon intensity
!     over a range of scattering angles mu and for each psi flux surface
!     Actually, the values given are I(mu,psi)/Clight.

use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none


!#include "lsc_xparams.h"
!#include "lsc_xray.h"
      INTEGER NPDIM, np, NVDIM, nv, iv, j, ip, iIterate
      INTEGER iv0
      INTEGER ith
      REAL    FeOfv(NVDIM, NPDIM, 2), sig_th(NMUDIM), Ifactor
      REAL    NeAry(NPDIM), ZbrAry(NPDIM), Vpar(NVDIM)
!     PARAMETER (MevPerMcsq = 0.511) delete this; add lines below
      real   signorm, Mevpermcsq
      COMMON / sigbrcom0 / signorm, Mevpermcsq
      REAL                                                              &
     &       ONE
      DATA   ONE/                                                       &
     &       1.0d0/
      iv0 = (nv + 1) / 2
!     velocity index for which Vpar(iv0) = 0
!     mu is measured relative to positive Vpar axis.
!     sigtot   INPUT  differential cross section times photon energy
!                     times transmission factor (.le. 1.0)
!                     w.r.t. solid angle(str) in
!                     MeV - cm^2 per atom per incident electron
!     sig_th   OUTPUT differential cross section times photon energy
!                     times transmission factor (.le. 1.0)
!                     w.r.t. solid angle(str) interpolated velocity in
!                     MeV - cm^2 per atom per incident electron
!     inten    OUTPUT joules/sec per cm^3 per steradian

      inten(:,:) = zero
!      CALL lsc_BrodCast(NMUDIM * NPDIM, inten, 0.0d0)
      do 100 iv = 1, nv
!.
        do 10 ith = 1, nMUbins
          if(Eint(iv) .GE. 1 .and. Eint(iv) .LT. nEbins) then
            sig_th(ith) = sigtot(ith, Eint(iv)) *                       &
     &                    (1.d0 - Efrac(iv)) +                   &
     &                    Efrac(iv) *sigtot(ith, Eint(iv) + 1)
          else
            sig_th(ith) = 0.00d0
          endif
 10     continue
!.
        if(iv .LE. iv0) then
          CALL lsc_xVecRefl(nMUbins, sig_th)
          Ifactor = 0.5 * abs((Vpar(iv + 1) + Vpar(iv)) *               &
     &         (Vpar(iv + 1) - Vpar(iv)))
       else
          Ifactor = 0.5 * abs((Vpar(iv) + Vpar(iv - 1)) *               &
     &         (Vpar(iv) - Vpar(iv - 1)))
        endif
!     note: dEofv ALWAYS .gt. 0, and Vpar(iv) increases from (-1., 1.),
!
        do 50 ip = 1, np
          do 40 j = 1, nMUbins
            inten(j, ip) = inten(j, ip) +                               &
     &                     FeOfv(iv, ip, iIterate) * sig_th(j) *        &
     &                     Ifactor
 40      continue
 50     continue
 100  continue
!
      do 150 ip = 1, np
        do 140 j = 1, nMUbins
          inten(j, ip) = inten(j, ip) * NeAry(ip)                       &
     &                   /  max (ONE ,ZbrAry(ip)) * 4.8e-3
!    ^                   / amax1(1.00,ZbrAry(ip)) * 4.8e-3

 140    continue
 150  continue
      return
      END

!     +---------------------------------------------------------------+

      SUBROUTINE lsc_xVecRefl(n, vec)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER i, n, i_comp
      REAL    vec(n), vec_tmp
      do i = 1, (n + 1) / 2
         i_comp = n + 1 - i
         vec_tmp = vec(i_comp)
         vec(i_comp) = vec(i)
         vec(i) = vec_tmp
      enddo
      return
      END

!     +---------------------------------------------------------------+

      REAL   FUNCTION lsc_delta_em(rr, zz, mu)
use iso_c_binding, only : fp => c_double
implicit none


      INTEGER ip, imu
      REAL(fp)    rr, zz, mu, ip_frac, imu_frac
      CALL lsc_xGetPsi(rr, zz, mu, ip, ip_frac, imu, imu_frac)
      CALL lsc_xPsiIntp(ip, ip_frac, imu, imu_frac, lsc_delta_em)
      lsc_delta_em = lsc_delta_em * Rmaj / rr
!     to account for convergence of flux tubes
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_xGetPsi(r, z, mu, ip, ip_frac, imu, imu_frac)
!     find psi of r,z
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none


!#include "lsc_emparams.h"
!#include "lsc_xparams.h"
!#include "lsc_xray.h"
      INTEGER ip, imu, jlo, jhi
      REAL    mu, ip_frac, imu_frac, rmu
      REAL    r,z,psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
      DATA    jlo / -1 /
      CALL lsc_plasma2d(r,z,psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      if (psi .GE. PsiAry(npsi) .or. pe2 .LE. 0.00d0 .or.        &
     &      r .GE. RlcfsMax     .or.   r .LE. RlcfsMin .or.             &
     &      z .LE. ZlcfsMin     .or.   z .GE. ZlcfsMax ) then
        ip      = 0
        ip_frac = 0.00d0
        imu     = 1
        imu_frac= 0.00d0
        return
      endif
      jlo = minloc(abs(psiary-psi),1)
      if (psi.lt.psiary(jlo)) jlo=jlo-1
!      CALL lsc_huntnr(PsiAry, npsi, psi, jlo)
      ip = jlo
      jhi= jlo+1
      if (jhi .GT. npsi) then
        ip_frac = 0.00d0
      else
        ip_frac = (psi - PsiAry(jlo))/(PsiAry(jlo+1) - PsiAry(jlo))
      endif
      rmu = (mu - mu_min) * dmu_inv + 1.00d0
      imu = int(rmu)
      imu_frac = rmu - AREAL(imu)
      return
      END

      SUBROUTINE lsc_xPsiIntp(ip, ip_frac, imu, imu_frac, lsc_delta_em)
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none


!#include "lsc_xparams.h"
!#include "lsc_xray.h"
      INTEGER ip, imu, ip_plus, imu_plus
      REAL    ip_frac, imu_frac, lsc_delta_em
!     linear intepolation
      imu_plus = min(imu + 1, nMUbins)
      ip_plus = min(ip + 1, npsi)
      if (ip .LE. 0 ) then
        lsc_delta_em = 0.00d0
      else
        lsc_delta_em = (1.d0 - ip_frac - imu_frac) * inten(imu, ip) &
     &     + ip_frac * inten(imu, ip_plus)                              &
     &     + imu_frac *inten(imu_plus, ip)
      endif
      return
      END
!
!     emit.F          --------------------------------------------------
!
!     produce bremstrahlung emission pattern from prescribed source
!     as observed by X-ray pinhole camera

      SUBROUTINE lsc_xEmit
use iso_c_binding, only : fp => c_double
implicit none

!cNCAR      CALL lsc_opngks
           CALL lsc_EZinit
      CALL lsc_X_setup_camera
      CALL lsc_X_compute_pattern
      CALL lsc_X_diagnose_trajectories
      CALL lsc_X_plot_pixel_data
      CALL lsc_X_plot_source_profile
           CALL lsc_MkGrfLst('X ray contour plots')
           CALL lsc_EZfini(0,0)
!!SLICES HERE?
           CALL lsc_EZinit
      CALL lsc_xSlices
      CALL lsc_eSlices
!     CALL lsc_XmnGraf

!     CALL lsc_X_wrres
!cNCAR      CALL lsc_clsgks
           CALL lsc_EZfini(0,0)

      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_rzphi_to_xyz(xyz, rzphi)
use iso_c_binding, only : fp => c_double
implicit none

      REAL   xyz(*), rzphi(*)
!     transform from cylindrical to cartesian coordinates
      xyz(3) = rzphi(2)
      xyz(1) = rzphi(1) * cos(rzphi(3))
      xyz(2) = rzphi(1) * sin(rzphi(3))
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_xyz_norm(norm, xyz)
use iso_c_binding, only : fp => c_double
implicit none

      REAL   norm, xyz(*)
      norm = sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2) + xyz(3) *          &
     &     xyz(3))
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_rzphi_norm(norm, rzphi)
use iso_c_binding, only : fp => c_double
implicit none

      REAL   norm, rzphi(*)
      norm = sqrt(rzphi(1) * rzphi(1) + rzphi(2) * rzphi(2))
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_det_33(det, mat)
use iso_c_binding, only : fp => c_double
implicit none

      REAL   det, mat(3, *)
      det = mat(1, 1) * mat(2, 2) * mat(3, 3) +                         &
     &     mat(1, 2) * mat(2, 3) * mat(3, 1) + mat(1, 3) *              &
     &     mat(2, 1) * mat(3, 2) - mat(1, 1) * mat(3, 2) *              &
     &     mat(2, 3) - mat(2, 1) * mat(1, 2) * mat(3, 3) -              &
     &     mat(3, 1) * mat(2, 2) * mat(1, 3)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_xyz_to_rzphi(rzphi, xyz)
use iso_c_binding, only : fp => c_double
use physconst_mod, only : pi
implicit none

!#include "lsc_numerics.h"
      REAL   xyz(*), rzphi(*)
!     transform from cartesian to cylindrical coordinates
      rzphi(2) = xyz(3)
      rzphi(1) = sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2))
      if(rzphi(1) .EQ. 0.d0)then
         rzphi(3) = 0.d0
      else if(xyz(1) .GT. 0.d0)then
         rzphi(3) = asin(xyz(2) / rzphi(1))
      else
         rzphi(3) = pi - asin(xyz(2) / rzphi(1))
      endif
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_setup_camera
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
!#include "lsc_emitter.h"
!#include "lsc_camera.h"
!#include "lsc_numerics.h"

      INTEGER ix, iy, i, j
      REAL   chi, mag_xy, del_R, del_z, rpy_sq(NPIXDIM),                &
     &     FocLen2, axis_orientation(3, 3),                             &
     &     chord_orientation(3), x_pix_min, x_pix_max, y_pix_min,       &
     &     y_pix_max, rpx_sq, rp_sq, phi_pix, mu_pix, xyz_pinhole(3),   &
     &     distance_pix_to_pinhole, det, norm, y_hat_xyz_pinhole(3)
      REAL    rt_rp_sq, asinarg, SMALL
      DATA    SMALL /1.0d-10/

      PI = 4. * atan(1.d0)
      R_bound_min_sq = R_bound_min * R_bound_min
      R_bound_max_sq = R_bound_max * R_bound_max

      CALL lsc_X_rzphi_to_xyz(xyz_pinhole, pinhole_loc)
!
!     begin by:
!     determining orientation of camera symmetry axis relative to
!     orientation of unit vector from pinhole (R, z, phi) to
!     tangency point prescribed by R_tangent, z_tangent
!
!     set up camera axis orientation -- spherical angles
      CALL lsc_X_normvec_xyz(3, y_hat_xyz_pinhole, xyz_pinhole)
      CALL lsc_X_xyz_norm(norm, y_hat_xyz_pinhole)
      chi = phi_pinhole + acos(R_tangent / Rpinhole)
      del_R = Rpinhole * cos(pi / 2. - (chi - phi_pinhole))
      del_z = (z_tangent - Zpinhole)
      mag_xy = sqrt(del_R * del_R + del_z * del_z)
      phi_axis = pi / 2. + chi
      mu_axis = acos(del_z / mag_xy)
      axis_orientation(1, 3) = sin(mu_axis) * cos(phi_axis)
      axis_orientation(2, 3) = sin(mu_axis) * sin(phi_axis)
      axis_orientation(3, 3) = cos(mu_axis)
      axis_orientation(1, 2) = sin(phi_axis)
      axis_orientation(2, 2) = - cos(phi_axis)
      axis_orientation(3, 2) = 0.d0
      axis_orientation(1, 1) = - cos(phi_axis) * cos(mu_axis)
      axis_orientation(2, 1) = - sin(phi_axis) * cos(mu_axis)
      axis_orientation(3, 1) = sin(mu_axis)
      CALL lsc_X_det_33(det, axis_orientation)
!     gives the (x, y, z) projections of the three camera axes
!
!     next:
!     set up to compute orientation (cartesian coordinates) of chords
!     originating from each element
!
!     first step:
!     determine angular offset, relative to camera axis, of chords
!     eminating from each CCD element
!     orientation of axes is that positive x axis lies in vertical plane
!     formed by camera symmetry axis and tokamak centerline (z axis)
!     in camera coordinates then:
!     symmetry axis is polar axis, x axis is origin for phi and
!     y = z cross vx
!     the pixels are assumed located on the x-y plane, a distance
!     focal_length from the pinhole, with dimensions
!     pix_sixe_x and pix_size_y
      FocLen2 = focal_length * focal_length
      x_pix_min = - 0.5 * AREAL(n_pixel_x) * pix_size_x
      x_pix_max = - x_pix_min
      y_pix_min = - 0.5 * AREAL(n_pixel_y) * pix_size_y
      y_pix_max = - y_pix_min
      CALL lsc_ugrid(x_pix_vec, n_pixel_x, x_pix_min, x_pix_max)
      CALL lsc_ugrid(y_pix_vec, n_pixel_y, y_pix_min, y_pix_max)
      do iy = 1, n_pixel_y
         rpy_sq(iy) = y_pix_vec(iy) * y_pix_vec(iy)
      enddo
      do ix = 1, n_pixel_x
         rpx_sq = x_pix_vec(ix) * x_pix_vec(ix)
         do iy = 1, n_pixel_y
            rp_sq = rpx_sq + rpy_sq(iy)
            rt_rp_sq = sqrt(rp_sq) * ( 1.0d0 + SMALL )
!           mu_pix = asin(sqrt(rp_sq / (rp_sq + FocLen2)))
            asinarg = sqrt( rp_sq / (rp_sq + FocLen2) )
            mu_pix = asin(asinarg)
            if(rp_sq .EQ. 0.)then
               phi_pix = 0.d0
            else if(x_pix_vec(ix) .GE. 0.d0)then
               asinarg = y_pix_vec(iy) / rt_rp_sq
!              CALL lsc_asinDebg(asinarg,ix,iy,n_pixel_x,n_pixel_y,2)
!              phi_pix = asin(y_pix_vec(iy) / sqrt(rp_sq))
               phi_pix = asin(asinarg)
            else
               asinarg = y_pix_vec(iy) / rt_rp_sq
!              CALL lsc_asinDebg(asinarg,ix,iy,n_pixel_x,n_pixel_y,3)
!              phi_pix = pi - asin(y_pix_vec(iy) / sqrt(rp_sq))
               phi_pix = pi - asin(asinarg)
            endif
!     generate the projection of chord orientation onto camera
!     axes
            phi_pix = pi + phi_pix
!     chord originates at pixel plane
            chord_orientation(1) = sin(mu_pix) * cos(phi_pix)
            chord_orientation(2) = sin(mu_pix) * sin(phi_pix)
            chord_orientation(3) = cos(mu_pix)
            CALL lsc_X_xyz_norm(norm, chord_orientation)
            distance_pix_to_pinhole = sqrt(rp_sq + FocLen2)
            do i = 1, 3
               y_hat_camera(i, ix, iy) = 0.d0
               do j = 1, 3
                  y_hat_camera(i, ix, iy) = y_hat_camera(i, ix, iy) +   &
     &                 axis_orientation(i, j) * chord_orientation(j)
!     in cartesian coordinates
               enddo
            enddo
            CALL lsc_X_xyz_norm(norm, y_hat_camera(1, ix, iy))
            do i = 1, 3
               chord_origin(i, ix, iy) = xyz_pinhole(i) -               &
     &              distance_pix_to_pinhole * y_hat_camera(i, ix, iy)
!     in cartesian coordinates
            enddo
         enddo
      enddo
      return
      END
!     .         --------------------------------------------------------
!      SUBROUTINE asinDebg(asinarg,ix,iy,n_pixel_x,n_pixel_y,icall)
!      if (asinarg .gt. 1.0 ) then
!        write(6,'('' asinarg,ix,iy,npixx,npixy,icall:''
!     ^  e20.13,5i4)')asinarg,ix,iy,n_pixel_x,n_pixel_y,icall
!        asinarg=+1.00
!      endif
!      if (asinarg .lt.-1.0 ) then
!        write(6,'('' asinarg,ix,iy,npixx,npixy,icall:''
!     ^  e20.13,5i4)')asinarg,ix,iy,n_pixel_x,n_pixel_y,icall
!        asinarg=-1.00
!      endif
!      return
!      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_normvec_xyz(n, y_hat, y)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER n
      REAL   norm, y_hat(n), y(n), ynorminv
      CALL lsc_X_xyz_norm(norm, y)
      ynorminv = 1.d0 / norm
      CALL lsc_X_vsmult(n, y_hat, y, ynorminv)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_vsmult(n, vout, vin, c)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER n, i
      REAL   vout(n), vin(n), c
      do i = 1, n
         vout(i) = c * vin(i)
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_diagnose_trajectories
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
!#include "lsc_camera.h"
      INTEGER ix, iy
      INTEGER pix_count(NPIXDIM)
      INTEGER iDBG
      DATA iDBG/0/
      if (iDBG .EQ. 0 ) return
      do ix = 1, NPIXDIM
         pix_count(ix) = ix
      enddo
      write(4, 101)(pix_count(iy), pix_count(iy), iy = 1, n_pixel_y)
 101  format(10(i5, 1x, i5, 2x))
      write(4, 102)
 102  format(/)
      do ix = 1, n_pixel_x
            write(4, 100)(r_tangent_actual(ix, iy),                     &
     &           z_tangent_actual(ix, iy), iy = 1, n_pixel_y)
 100        format(10(f5.1, 1x, f5.1, 2x))
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_compute_pattern
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
!#include "lsc_camera.h"
!#include "lsc_numerics.h"
      INTEGER ix, iy
      REAL   y_beg_rzphi(3), y_hat_xyz(3), y_end(3), SolAng
      solAng = 3.14159*(pinhole_size/focal_length)**2 * 0.25d+2
      do ix = 1, n_pixel_x
         do iy = 1, n_pixel_y
            CALL lsc_X_init_chord(ix, iy, y_beg_rzphi, y_hat_xyz)
            CALL lsc_X_c_phot(y_beg_rzphi, y_hat_xyz, dlxray, y_end, 3, &
     &           photon_count(ix, iy), npoints(ix, iy),                 &
     &           chord(1, 1, ix, iy), r_tangent_actual(ix, iy),         &
     &           z_tangent_actual(ix, iy))
            photon_count(ix, iy) =  photon_count(ix, iy)*SolAng
         enddo
      enddo
!      CALL lsc_LSCpause
!      do ix = 1, n_pixel_x
!         do iy = 1, n_pixel_y
!
!        if (pc(ix,iy) .gt. 1.e+20 .or.
!     ^      pc(ix,iy) .lt.-1.e+20  ) then
!            write(6,'('' pc is bad at i,j;'', e10.3, i4, i4)')
!     ^      pc(ix,iy), ix, iy
!            pc(ix,iy) = 0.
!         endif
!         enddo
!      enddo
!
!      CALL lsc_LSCpause

      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_init_chord(ix, iy, y_beg_rzphi, y_hat_xyz)
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
!#include "lsc_camera.h"
      INTEGER ix, iy
      REAL    y_beg_rzphi(*), y_hat_xyz(*)
      CALL lsc_X_xyz_to_rzphi(y_beg_rzphi, chord_origin(1, ix, iy))
      CALL lsc_xVecCopy(3, y_hat_xyz, y_hat_camera(1, ix, iy))
      return
      END
!     .         --------------------------------------------------------
      INTEGER FUNCTION lsc_inbound(rsq, rr, zz)
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
!#include "lsc_emitter.h"
      REAL   rsq, rr, zz
      if((rsq .GE. R_bound_min_sq) .and. (rsq .LE. R_bound_max_sq)      &
     &     .and.                                                        &
     &     (zz .GE. Z_bound_min) .and. (zz .LE. Z_bound_max)            &
     &     .and.                                                        &
     &     ( (rr-PusherMajor)**2 + zz**2 .GE. PusherMinor**2) ) then
         lsc_inbound = 1
      else
         lsc_inbound = 0
      endif
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_plot_source_profile
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
!#include "lsc_emitter.h"
!#include "lsc_camera.h"
      INTEGER iz, ir
!      REAL   lsc_delta_em, rr, zz
      CALL lsc_ugrid(r_source, nr_source, R_plasm_min, R_plasm_max)
      CALL lsc_ugrid(z_source, nz_source, Z_plasm_min, Z_plasm_max)
      do iz = 1, nz_source
         zz = z_source(iz)
         do ir = 1, nr_source
            rr = r_source(ir)
            source_profile(ir, iz) = lsc_delta_em(rr, zz, mu_0)
         enddo
      enddo
      CALL lsc_X_pl_emisivity_cntrs(source_profile, NRDIM, nr_source,   &
     &     nz_source,r_source,z_source)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_plot_pixel_data
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
!#include "lsc_emitter.h"
!#include "lsc_camera.h"
      INTEGER nr, DOCHORDPLOT
      REAL    Radii(2), Z_height
      DATA DOCHORDPLOT / 0 /
!     plot chords in (x, y) and (r, z) views
      nr = 2
      Radii(1) = R_bound_min
      Radii(2) = R_bound_max
      Z_height = Z_bound_max
      if(DOCHORDPLOT .EQ. 1 ) then
      CALL lsc_X_plan_view(nr, Radii, npoints, n_pixel_x, n_pixel_y, chord, &
     &     MAXPOINTS, NPIXDIM)
      CALL lsc_X_elevation_view(nr, Radii, - Z_height, Z_height, npoints, &
     &     n_pixel_x, n_pixel_y, chord, MAXPOINTS, NPIXDIM)
      endif
      CALL lsc_X_pl_intensity_cntrs(photon_count, NPIXDIM, n_pixel_x,   &
     &     n_pixel_y)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_xSlices
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
!#include "lsc_emitter.h"
!#include "lsc_camera.h"
!     DATA RcntChrd  / 1.65 / ! the nominal radius of the usual central
!                             ! chord plotted in the time development of the
!                             ! xray signal
      CHARACTER*70 MyString
      INTEGER i,j,icenter,jcenter, ioffset, joffset
      REAL    Tan2Pin, SHAVE
      REAL    slice(NPIXDIM,2),                                         &
     &        signal(NPIXDIM, NCHORDIM), sigMax, xmaxy
      DATA    SHAVE / 0.985d0 /
      REAL    ZERO, ONE
      DATA          ZERO   ,       ONE    /                             &
     &        0.0d0 , 1.0d0 /
      REAL TENm30
      DATA TENm30 / 1.0d-30 /

      Tan2Pin = sqrt ( Rpinhole**2 - R_tangent**2 )
!c23456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
      joffset = n_pixel_y *                                             &
     &  abs(RcntChrd-R_tangent) / Tan2Pin * focal_length / (screen_d)
      ioffset = n_pixel_x *                                             &
     &  abs(RcntChrd-R_tangent) / Tan2Pin * focal_length / (screen_d)
      icenter = (n_pixel_x + 1)/2
      jcenter = (n_pixel_y + 1)/2
!
!
      sigMax=0.00
      do 10 i=1, n_pixel_x
        slice(i,1) = AREAL(i)
        signal(i,1) = photon_count(i,(jcenter))
        signal(i,2) = photon_count(i,(jcenter-joffset))
        signal(i,3) = photon_count(i,(jcenter-joffset/2))
        signal(i,4) = photon_count(i,(jcenter+joffset))
        signal(i,5) = photon_count(i,(jcenter+joffset/2))
        sigMax=  max (sigMax, signal(i,1) )
        sigMax=  max (sigMax, signal(i,2) )
        sigMax=  max (sigMax, signal(i,3) )
        sigMax=  max (sigMax, signal(i,4) )
        sigMax=  max(sigMax, signal(i,5) )
 10   continue
      do 11 j=1, n_pixel_y
        slice(j,2) = AREAL(j)
        signal(j,6) = photon_count(icenter,j)
        sigMax=  max (sigMax, signal(j,6))
 11   continue
      CALL lsc_XRopen(8,'XRinput')
      CALL lsc_XRdscw(slice(1,1), signal(1,2), n_pixel_x,               &
     &            'Vertical x ray slice')
      CALL lsc_XRclos(8)

      do 15 i=1, n_pixel_x
        signal(i,1) = signal(i,1)/(sigMax+TENm30)*SHAVE
        signal(i,2) = signal(i,2)/(sigMax+TENm30)*SHAVE
        signal(i,3) = signal(i,3)/(sigMax+TENm30)*SHAVE
        signal(i,4) = signal(i,4)/(sigMax+TENm30)*SHAVE
        signal(i,5) = signal(i,5)/(sigMax+TENm30)*SHAVE
 15   continue
      do 16 j=1, n_pixel_y
        signal(j,6) = signal(j,6)/(sigMax+TENm30)*SHAVE
 16   continue
!
!     ULH corner
      CALL lsc_EZrnd2 (AREAL(n_pixel_x),xmaxy,i,j)
      CALL lsc_EZsets(100,500,450,700,                                  &
     & ZERO,xmaxy,ZERO, ONE, 1)
!    ^ 0.,xmaxy,0.00, 1.00, 1)
      CALL lsc_EZaxes(1,2,1,5)
      CALL lsc_EZcurv(slice(1,1),signal(1,1),n_pixel_x)
      write(MyString,'(''Vert slice /'',1pe9.2,                         &
     &    '' at tangent$'')') sigMax
      CALL lsc_EZwrit(100,725,MyString,0,0)
!     URH corner
      CALL lsc_EZsets(600,1000,450,700,                                 &
     & ZERO,xmaxy,ZERO, ONE, 1)
!    ^ 0.,xmaxy,0.00, 1.00, 1)
      CALL lsc_EZaxes(1,2,1,5)
      CALL lsc_EZcurv(slice(1,1),signal(1,2),n_pixel_x)
      CALL lsc_EZcros(slice(1,1),signal(1,3),n_pixel_x)

      write(MyString, '(''at chord'',i3,'' and'',i3,'' (dots)$'')')     &
     &                   (jcenter-joffset), (jcenter-joffset/2)
      CALL lsc_EZwrit(600,725,MyString,0,0)
!
!     LLH corner
      CALL lsc_EZsets(100,500,100,350,                                  &
     & ZERO,xmaxy,ZERO, ONE, 1)
!    ^ 0.,xmaxy,0.00, 1.00, 1)
      CALL lsc_EZaxes(1,2,1,5)
      CALL lsc_EZcurv(slice(1,1),signal(1,4),n_pixel_x)
      CALL lsc_EZcros(slice(1,1),signal(1,5),n_pixel_x)
      write(MyString, '(''at chord'',i3,'' and'',i3,'' (dots)$'')')     &
     &                   (jcenter+joffset), (jcenter+joffset/2)
      CALL lsc_EZwrit(100,375,MyString,0,0)
!
!     LRH corner
      CALL lsc_EZrnd2 (AREAL(n_pixel_y),xmaxy,i,j)
      CALL lsc_EZsets(600,1000,100,350,                                 &
     & ZERO,xmaxy,ZERO, ONE, 1)
!    ^ 0.,xmaxy,0.00, 1.00, 1)
      CALL lsc_EZaxes(1,2,1,5)
      CALL lsc_EZcurv(slice(1,2),signal(1,6),n_pixel_y)
      write(MyString,'(''Horiz slice '',''at mid plane$'')')
      CALL lsc_EZwrit(600,375,MyString,0,0)
!
           CALL lsc_MkGrfLst('X ray slices; 4 graphs')
           CALL lsc_EZfini(0,0)
!
!
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_eSlices
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
!#include "lsc_emitter.h"
!#include "camera.inc"
      CHARACTER*70 MyString
      INTEGER i,j,icenter,jcenter, ioffset, joffset
      REAL    SHAVE
      REAL    signal(NRZMXDIM, NCHORDIM), sigMax
      REAL    rminy, rmaxy, rdumy, zmaxy
      DATA    SHAVE / 0.985d0 /
      REAL          ZERO,          ONE
      DATA          ZERO   ,       ONE    /                             &
     &        0.0d0 , 1.0d0 /
      REAL TENm30
      DATA TENm30 / 1.0d-30 /

      ioffset = nr_source / 4
      joffset = nz_source / 4
      icenter = (nr_source + 1)/2
      jcenter = (nz_source + 1)/2
!
!
      sigMax=ZERO
      do 10 j=1, nz_source
        signal(j,1) = source_profile((icenter          ),j)
        signal(j,2) = source_profile((icenter+ioffset  ),j)
        signal(j,3) = source_profile((icenter+ioffset/2),j)
        signal(j,4) = source_profile((icenter-ioffset  ),j)
        signal(j,5) = source_profile((icenter-ioffset/2),j)
        sigMax=  max (sigMax, signal(j,1) )
        sigMax=  max (sigMax, signal(j,2) )
        sigMax=  max (sigMax, signal(j,3) )
        sigMax=  max (sigMax, signal(j,4) )
        sigMax=  max (sigMax, signal(j,5) )
 10   continue
      do 11 i=1, nr_source
        signal(i,6) = source_profile(i,jcenter)
        sigMax=  max (sigMax, signal(i,6))
 11   continue

      do 15 j=1, nz_source
        signal(j,1) = signal(j,1)/(sigMax+TENm30)*SHAVE
        signal(j,2) = signal(j,2)/(sigMax+TENm30)*SHAVE
        signal(j,3) = signal(j,3)/(sigMax+TENm30)*SHAVE
        signal(j,4) = signal(j,4)/(sigMax+TENm30)*SHAVE
        signal(j,5) = signal(j,5)/(sigMax+TENm30)*SHAVE
 15   continue
      do 16 i=1, nr_source
        signal(i,6) = signal(i,6)/(sigMax+TENm30)*SHAVE
 16   continue
!
!     ULH corner
      CALL lsc_EZrnd2( z_source(nz_source),zmaxy,i,j )
      CALL lsc_EZsets(100,500,450,700,                                  &
     & -zmaxy,zmaxy,ZERO, ONE, 1)
!    ^ -zmaxy,zmaxy,0.00, 1.00, 1)
      CALL lsc_EZaxes(i,j,1,5)
      CALL lsc_EZcurv(z_source,signal(1,1),nz_source)
      write(MyString,'(''Vert slice /'',1pe9.2,                         &
     &    '' at center$'')') sigMax
      CALL lsc_EZwrit(100,725,MyString,0,0)
!
!     URH corner
      CALL lsc_EZsets(600,1000,450,700,                                 &
     & -zmaxy,zmaxy,ZERO,  ONE, 1)
!    ^ -zmaxy,zmaxy,0.00, 1.00, 1)
      CALL lsc_EZaxes(i,j,1,5)
      CALL lsc_EZcurv(z_source,signal(1,2),nz_source)
      CALL lsc_EZcros(z_source,signal(1,3),nz_source)

      write(MyString, '(''at r_loc'',i3,'' and'',i3,'' (dots)$'')')     &
     &                   (icenter+ioffset), (icenter+ioffset/2)
      CALL lsc_EZwrit(600,725,MyString,0,0)
!
!     LLH corner
      CALL lsc_EZsets(100,500,100,350,                                  &
     & -zmaxy,zmaxy,ZERO,  ONE, 1)
!    ^ -zmaxy,zmaxy,0.00, 1.00, 1)
      CALL lsc_EZaxes(i,j,1,5)
      CALL lsc_EZcurv(z_source,signal(1,4),nz_source)
      CALL lsc_EZcros(z_source,signal(1,5),nz_source)
      write(MyString, '(''at r_loc'',i3,'' and'',i3,'' (dots)$'')')     &
     &                   (icenter-ioffset), (icenter-ioffset/2)
      CALL lsc_EZwrit(100,375,MyString,0,0)
!
!     LRH corner
      CALL lsc_EZrnd2 (r_source(nr_source), rmaxy, i,j)
      CALL lsc_EZrnd2 (rmaxy-r_source(1),rdumy, i,j)
      rminy = rmaxy - rdumy
      CALL lsc_EZsets(600,1000,100,350,                                 &
     & rminy,rmaxy,ZERO,  ONE, 1)
!    ^ rminy,rmaxy,0.00, 1.00, 1)
      CALL lsc_EZaxes(i,j,1,5)
      CALL lsc_EZcurv(r_source,signal(1,6),nr_source)
      write(MyString,'(''Horiz slice '',''at mid plane$'')')
      CALL lsc_EZwrit(600,375,MyString,0,0)
!
           CALL lsc_MkGrfLst('Emissivity slices; 4 graphs')
           CALL lsc_EZfini(0,0)
!
      return
      END
!     .         --------------------------------------------------------
!      SUBROUTINE XmnGraf
!#include "Implic.inc"
!#include "xparams.inc"
!#include "Xray.inc"

!...
!      return
!      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_wrres
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none


!#include "lsc_emparams.h"
!#include "lsc_xparams.h"
!#include "lsc_xray.h"
!#include "lsc_emitter.h"
!#include "lsc_camera.h"
!     write results to disk file

      INTEGER ix, iy
      INTEGER iDBG
      DATA iDBG/1/
      if (iDBG .EQ. 0 ) return
!.
      write(4, 102) 'photon_c'
 102  format(/,a10)
      do ix = 1, NPIXDIM
      write(4, 101)(photon_count(ix,iy), iy = 1, NPIXDIM)
 101  format(8(1pe10.1))
      enddo
      write(4, 102)
!.
      write(4, 102) 'source_p'
      do ix = 1, NRDIM
      write(4, 101)(source_profile(ix,iy), iy = 1, NZDIM)
      enddo
      write(4, 102)
!.
      write(4, 102) 'sigtot'
      do ix = 1, NMUDIM
      write(4, 101)(sigtot(ix,iy), iy = 1, NENDIM)
      enddo
      write(4, 102)
!.
      write(4, 102) 'inten'
      do ix = 1, NMUDIM
      write(4, 101)(inten(ix,iy), iy = 1, NPSIDIM)
      enddo
      write(4, 102)
!.

      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_c_phot(y_beg, y_hat, dlchord, y_end, n_space, pc, &
     &     npoints, chord, r_tangent_actual, z_tangent_actual)
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
!#include "lsc_numerics.h"
      INTEGER n_space, lsc_inbound, i, count, state,                    &
     &     current_status, previous_status
!     INTEGER*4 npoints
      INTEGER   npoints
      REAL                                                              &
     &     y_beg(n_space), dlchord, y(SPACEDIM), y_hat(SPACEDIM),       &
     &     lsc_delta_em, y_end(n_space), lsc_mu_yb, dy_xyz(3),          &
     &     r_tangent_actual, z_tangent_actual, radius_sq,               &
     &     radius_sq_min, rr, zz
!     REAL*4 chord(3, *), pc
      REAL   chord(3, *), pc
      CALL lsc_X_rzphi_to_xyz(y, y_beg)
      do i = 1, 3
         dy_xyz(i) = y_hat(i) * dlchord
      enddo
      pc = 0.d0
      count = 0
      state = 1
      radius_sq = y(1) * y(1) + y(2) * y(2)
      rr = sqrt (radius_sq)
      radius_sq_min = radius_sq
      previous_status = lsc_inbound(radius_sq, rr , y(3))
      current_status = previous_status
      do while((state .LT. 3) .and. (count .LT. npoints_max))
         do i = 1, 3
            chord(i, count + 1) = y(i)
            y(i) = y(i) + dy_xyz(i)
         enddo
         radius_sq = y(1) * y(1) + y(2) * y(2)
         rr = sqrt (radius_sq)
         if(radius_sq_min .GT. radius_sq)then
            radius_sq_min = radius_sq
            z_tangent_actual = y(3)
         endif
         count = count + 1
         if(state .EQ. 2)then
!           rr = sqrt(radius_sq) ! we already did this
            zz = y(3)
            pc = pc + lsc_delta_em(rr, zz, lsc_mu_yb(y, y_hat)) * dlchord
         endif
         previous_status = current_status
         current_status = lsc_inbound(radius_sq, rr, y(3))
         if(current_status .ne. previous_status)  then
            state = state + 1
         endif
      enddo
      r_tangent_actual = sqrt(radius_sq_min)
      npoints = count
      CALL lsc_X_xyz_to_rzphi(y_end, y)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_advance_y(y, dy_xyz)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER i
      REAL   y(*), dy_xyz(*), y_old(3)
      CALL lsc_X_rzphi_to_xyz(y_old, y)
!     CALL lsc_X_project_dy(y, dy_xyz, dy_rzphi)
      do i = 1, 3
         y(i) = y_old(i) + dy_xyz(i)
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_project_dy(y_rzphi, dy_xyz, dy_rzphi)
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
!#include "lsc_numerics.h"
!     given y (R, z, phi) and dy (dx, dy, dz), compute (dR, dz, dphi)
      INTEGER i, j
      REAL   y_rzphi(*), dy_xyz(*), dy_rzphi(*), T(3, 3), dRdx, dRdy,   &
     &     dRdz, dzdx, dzdy, dzdz, dphidx, dphidy, dphidz,              &
     &     xx, yy, Rsq_inv
      EQUIVALENCE                                                       &
     &     (T(1, 1), dRdx), (T(1, 2), dRdy), (T(1, 3), dRdz),           &
     &     (T(2, 1), dzdx), (T(2, 2), dzdy), (T(2, 3), dzdz),           &
     &     (T(3, 1), dphidx), (T(3, 2), dphidy), (T(3, 3), dphidz)
      DATA dRdz / 0.d0 /
      DATA        dzdx  ,       dzdy  ,       dzdz  /                   &
     &      0.0d0, 0.0d0, 1.d0 /
      DATA dphidz / 0.d0 /
      xx = y_rzphi(1) * cos(y_rzphi(3))
      yy = y_rzphi(1) * sin(y_rzphi(3))
      dRdx = xx / y_rzphi(1)
      dRdy = yy / y_rzphi(1)
      Rsq_inv = 1.0d0 / (y_rzphi(1) * y_rzphi(1))
      dphidx = - yy * Rsq_inv
      dphidy = xx * Rsq_inv
      do i = 1, 3
         dy_rzphi(i) = 0.d0
         do j = 1, 3
            dy_rzphi(i) = dy_rzphi(i) + T(i, j) * dy_xyz(j)
         enddo
      enddo
      return
      END
!     .         --------------------------------------------------------
      REAL   FUNCTION lsc_mu_yb(y, y_hat)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER i
      REAL   y(*), y_hat(*), inner_prod, B_vector(3), rr
!     B(1) = B_x, B(2) = B_y, B(3) = B_z
!     construct cos angle between magnetic field (assumed in phi direction)
!     and y_hat at point y (in r, Z, phi, coordinates) and
!     unit vector y_hat (in cartesian x, y, z, coordinates)
      rr = sqrt(y(1) * y(1) + y(2) * y(2))
      B_vector(1) = y(2) / rr
      B_vector(2) = - y(1) / rr
      B_vector(3) = 0.d0
      inner_prod = 0.d0
      do i = 1, 3
         inner_prod = inner_prod + B_vector(i) * y_hat(i)
      enddo
      lsc_mu_yb = inner_prod
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_xVecCopy(n, tovec, frvec)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER n, i
      REAL    tovec(n), frvec(n)
      do i = 1, n
         tovec(i) = frvec(i)
      enddo
      return
      END
!     .         --------------------------------------------------------
!
!     graphs    --------------------------------------------------------
      SUBROUTINE lsc_X_plan_view(nr, Radii, npoints,nx,ny, chord4, n1,n34)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER n34, nx, ny, nr, n1
!     INTEGER*4 npoints(n34, n34)
      INTEGER   npoints(n34, n34)
      REAL    R_min, R_max, Radii(nr)
!     REAL*4  chord4(3, n1, n34, n34)
      REAL    chord4(3, n1, n34, n34)
!      CALL lsc_VecMnMx( Radii, nr, R_min, R_max)
      R_min = minval(Radii)
      R_max = maxval(Radii)
      CALL lsc_xPlanSet (R_min, R_max)
      CALL lsc_PlPlBnds(nr, Radii)
      CALL lsc_PlPlChrds (npoints, nx, ny, chord4, n1, n34)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_xPlanSet (R_min, R_max)
use iso_c_binding, only : fp => c_double
implicit none

!     INTEGER*4 I1, I2
!     DATA      I1, I2 / 1, 2 /
      REAL   R_min, R_max
!     REAL*4 xdum4(2), ydum4(2), twormax
      REAL   twormax


      twormax = R_max + R_max
!cNCAR      CALL lsc_agsetf('X/MINIMUM.', - twormax)
!cNCAR      CALL lsc_agsetf('X/MAXIMUM.', twormax)
!cNCAR      CALL lsc_agsetf('Y/MINIMUM.', - twormax)
!cNCAR      CALL lsc_agsetf('Y/MAXIMUM.', twormax)
!cNCAR      CALL lsc_agstup(xdum4, 1, 1, 2, 1, ydum4, 1, 1, 2, 1)
!cNCAR      CALL lsc_agback
           CALL lsc_EZsets (50,550,  50,550,                            &
     &        -twormax, twormax, -twormax, twormax, 1)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_PlPlBnds(nr, Radii)
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_xgraphs.h"
!     INTEGER*4 nthet4, I1
      INTEGER   nthet4, I1
      INTEGER   i, nr, ir
      REAL      Radii(nr)
!     REAL*4  xva4(NWKDIX), yva4(NWKDIX),
!    ^        costhet4(NWKDIX), sinthet4(NWKDIX)
      REAL    xva4(NWKDIX), yva4(NWKDIX),                               &
     &        costhet4(NWKDIX), sinthet4(NWKDIX)
!     REAL*4    dthet4, pi, rim1, thet4
      REAL      dthet4, pi, rim1, thet4
      EQUIVALENCE (xva4(1), wkarx(1, 1)), (yva4(1), wkarx(1, 2)),       &
     &     (costhet4(1), wkarx(1, 3)), (sinthet4(1), wkarx(1, 4))
      DATA I1/ 1 /
      nthet4 = NWKDIx
      pi = 4. * atan(1.d0)
      dthet4 = (2. * pi) / AREAL(nthet4 - 1)
      do i = 1, nthet4
         rim1 = AREAL(i - 1)
         thet4 = dthet4 * rim1
         costhet4(i) = cos(thet4)
         sinthet4(i) = sin(thet4)
      enddo
      do ir = 1, nr
         do i = 1, nthet4
            xva4(i) = Radii(ir) * costhet4(i)
            yva4(i) = Radii(ir) * sinthet4(i)
         enddo
!cNCAR         CALL lsc_agcurv(xva4, 1, yva4, 1, nthet4, 1)
              CALL lsc_EZcurv(xva4, yva4, nthet4)
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_PlPlChrds(npoints, nx, ny, chord4, n1, n34)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER   n1, n34
      INTEGER   npoints(n34, n34), nx, ny, ix, iy, np
      REAL      chord4(3, n1, n34, n34)

      do ix = 1, nx, 2
        do iy = 1, ny, 2
!cncar     CALL lsc_agcurv(chord4(1, 1, ix, iy), 3, chord4(2, 1, ix, iy),3,
!cncar^    npoints(ix, iy), 1)
!         CALL lsc_EZcurvmd(chord4(1, 1, ix, iy), 3, chord4(2, 1, ix, iy),3,
!    ^    npoints(ix, iy))

          np = npoints(ix, iy) - 1

          CALL lsc_EZdra(chord4(1, 1, ix, iy), chord4(2, 1, ix, iy), 0 )
          CALL lsc_EZdra(chord4(1, 1, ix, iy), chord4(2, 1, ix, iy), 1 )
          CALL lsc_EZdra(chord4(1,np, ix, iy), chord4(2,np, ix, iy), 1 )

        enddo
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_elevation_view(nr, Radii, Z_min, Z_max, npoints, nx, &
     &     ny, chord4, n1, n34)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER   nr, n34, nx, ny, n1
      INTEGER*4 npoints(n34, n34)
      REAL   R_min, R_max, Radii(nr), Z_min, Z_max
!     REAL*4 chord4(3, n1, n34, n34)
      REAL   chord4
      REAL    ZERO
      DATA    ZERO         /                                            &
     &        0.0d0 /
!     CALL lsc_VecMnMx ( Radii, nr,  R_min, R_max)
      R_min = minval(Radii)
      R_max = maxval(Radii)
!     CALL lsc_xElevSet (0., R_max, Z_min, Z_max)
      CALL lsc_xElevSet (ZERO, R_max, Z_min, Z_max)
      CALL lsc_PlElBnds (R_min, R_max, Z_min, Z_max)
      CALL lsc_PlElChrds (npoints, nx, ny, chord4, n1, n34)
           CALL lsc_MkGrfLst('Xray plan, elevation chords')
           CALL lsc_EZfini(0,0)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_PlElBnds (R_min, R_max, Z_min, Z_max)
use iso_c_binding, only : fp => c_double
implicit none

!     INTEGER*4 i, I1, I2, I3, I4, I5
      INTEGER   i, I1, I2, I3, I4, I5
      REAL    R_min, R_max, Z_min, Z_max, ris
!     REAL*4  xva4 (5, 2), yva4 (5, 2)
      REAL    xva4 (5, 2), yva4 (5, 2)
      DATA I1, I2, I3, I4, I5 / 1, 2, 3, 4, 5 /
      REAL    ZERO
      DATA    ZERO /                                                    &
     &        0.0d0  /

!     i = 1: right half
!     i = 2: left half
      do i = I1, I2
         if(i .EQ. I1) then
            ris = 1.d0
         else
            ris = -1.d0
         endif
         xva4(I1, i) = ris * R_min
         yva4(I1, i) = Z_min
         xva4(I2, i) = ris * R_max
         yva4(I2, i) = Z_min
         xva4(I3, i) = ris * R_max
         yva4(I3, i) = Z_max
         xva4(I4, i) = ris * R_min
         yva4(I4, i) = Z_max
         xva4(I5, i) = xva4(I1, i)
         yva4(I5, i) = yva4(I1, i)

!cncar      CALL lsc_agcurv  (xva4(I1, i),I1, yva4(I1, i),I1,I5,I1)
           CALL lsc_EZcurvmd(xva4(I1, i),I1, yva4(I1, i),I1,I5)

      enddo
!          CALL lsc_EZdra(0.00, Z_min    , 0)
!          CALL lsc_EZdra(0.00, Z_min    , 1)
!          CALL lsc_EZdra(0.00, Z_min/10., 1)
!          CALL lsc_EZdra(0.00, 0.00     , 0)
!          CALL lsc_EZdra(0.00, 0.00     , 1)
!          CALL lsc_EZdra(0.00, Z_max/10., 0)
!          CALL lsc_EZdra(0.00, Z_max/10., 1)
!          CALL lsc_EZdra(0.00, Z_max    , 1)

           CALL lsc_EZdra(ZERO, Z_min    , 0)
           CALL lsc_EZdra(ZERO, Z_min    , 1)
           CALL lsc_EZdra(ZERO, Z_min/10., 1)
           CALL lsc_EZdra(ZERO, ZERO     , 0)
           CALL lsc_EZdra(ZERO, ZERO     , 1)
           CALL lsc_EZdra(ZERO, Z_max/10., 0)
           CALL lsc_EZdra(ZERO, Z_max/10., 1)
           CALL lsc_EZdra(ZERO, Z_max    , 1)

      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_xElevSet (R_min_pl, R_max, Z_min, Z_max)
use iso_c_binding, only : fp => c_double
implicit none

!     INTEGER*4 I1, I2
      INTEGER   I1, I2
      REAL   R_min_pl, R_max, Z_min, Z_max, Z_mn, Z_mx
!     REAL*4 xdum4(2), ydum4(2), twormax, twozmin, twozmax
!     REAL   xdum4(2), ydum4(2)
      REAL   twormax, twozmin, twozmax
!     REAL*4 R_min_pl4
      REAL   R_min_pl4
      DATA I1, I2 / 1,2 /

      Z_mn = Z_min - 0.5 * (Z_max - Z_min)
      Z_mx = Z_max + 0.5 * (Z_max - Z_min)

      R_min_pl4 = R_min_pl

      twormax = 2. * R_max
      twozmin = 2. * Z_mn
      twozmax = 2. * Z_mx


!cNCAR      CALL lsc_agsetf('X/MINIMUM.', R_min_pl4)
!cNCAR      CALL lsc_agsetf('X/MAXIMUM.', twormax)
!cNCAR      CALL lsc_agsetf('Y/MINIMUM.', twozmin)
!cNCAR      CALL lsc_agsetf('Y/MAXIMUM.', twozmax)
!cncar      CALL lsc_agstup(xdum4,I1,I1,I2,I1, ydum4,I1,I1,I2,I1)
!cNCAR      CALL lsc_agback
           CALL lsc_EZsets (700,1000, 150,450,                          &
     &        R_min_pl4, twormax, twozmin, twozmax, 1)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_PlElChrds (npoints, nx, ny, chord4, n1, n34)
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_xgraphs.h"
      INTEGER   n34, nx, ny, ix, iy, n1, ip, np
!     INTEGER*4 npoints(n34, n34), I1
      INTEGER   npoints(n34, n34), I1
!     REAL*4 chord4(3, n1, n34, n34)
      REAL   chord4(3, n1, n34, n34)
!     REAL*4 rr(1), zz(1)
      REAL   rr(NWKDIx), zz(NWKDIx)
      EQUIVALENCE (rr(1), wkarx(1, 1)), (zz(1), wkarx(1, 2))
      DATA I1 / 1 /
      do ix = 1, nx, 4
         do iy = 1, ny, 4
            do ip = 1, npoints(ix, iy)
               CALL lsc_xyz2rz4(chord4(1, ip, ix, iy),I1, rr(ip), zz(ip))
            enddo
!cncar            CALL lsc_agcurv  (rr,I1, zz,I1, npoints(ix, iy),I1)
!                CALL lsc_EZcurvmd(rr,I1, zz,I1, npoints(ix, iy))

            np = npoints(ix, iy) - 1

            CALL lsc_EZdra(rr( 1), zz( 1), 0 )
            do ip=1,np,30
              CALL lsc_EZdra(rr(ip), zz(ip), 1 )
            enddo
         enddo
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_xyz2rz4(xyz, n, rr, zz)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER n, i
      REAL    xyz(3, n)
!     REAL*4  rr(n), zz(n)
      REAL    rr(n), zz(n)
      do i = 1, n
         rr(i) = sqrt(xyz(1, i) * xyz(1, i) + xyz(2, i) * xyz(2, i))
         zz(i) = xyz(3, i)
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE lsc_X_pl_intensity_cntrs(pc, nxdim, nx, ny)
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_emparams.h"
      CHARACTER*5 whichside
      CHARACTER*30 ChStr
      INTEGER   nxdim, nx, ny, iside, ix1(2), ix2(2), iy1, iy2
!     REAL*4 pc (nxdim, ny)
      REAL   pc (nxdim, ny)
!.
      REAL    pcmin, pcmax, pcave
      REAL    xminy, xmaxy, xdumy, yminy, ymaxy
      REAL    xAry(NPIXDIM), yAry(NPIXDIM), clevelin(20)
      REAL    xGiv(*)      , yGiv(*)
      INTEGER i,j, kclev1, kclev2
      INTEGER ixmax,ixmin,ixstep
      INTEGER jymax,jymin,jystep
      DATA    ix1(1), ix2(1), ix1(2), ix2(2), iy1, iy2 /                &
     &           100,  500,    600,    1000 , 200, 600 /
      whichside = 'LEFT'
      iside = 1
      do 10 i=1,nx
        xAry(i) = AREAL(i)
 10   continue
      do 11 j=1,ny
        yAry(j) = AREAL(j)
 11   continue
      CALL lsc_EZrnd2 (xAry(nx)  , xmaxy, i,j)
      CALL lsc_EZrnd2 (yAry(ny)  , ymaxy, i,j)
      xminy = 0.00d0
      yminy = 0.00d0

      kclev1 = -5
      kclev2 = 2

!
      goto 15
      ENTRY lsc_X_pl_emisivity_cntrs(pc, nxdim, nx, ny,                 &
     &           xGiv, yGiv)
      whichside = 'RIGHT'
      iside = 2

!  round up from max found--> nice value; #major divns; #minor divns
      CALL lsc_EZrnd2 (xGiv(nx)  , xmaxy, i,j)
      CALL lsc_EZrnd2 (xmaxy-xGiv(1), xdumy, i,j)
      CALL lsc_EZrnd2 (yGiv(ny)  , ymaxy, i,j)
      if (2.*ymaxy .GT. xdumy) then
        xminy = xmaxy - 2.0*ymaxy
        yminy = - ymaxy
      else
        xminy = xmaxy - xdumy
        ymaxy = xdumy / 2.
        yminy =-ymaxy
      endif

      kclev1 = -5
      kclev2 = 2
!
!     ...................
!
 15   continue

      pcmin = pc(1,1)
      pcmax = pc(1,1)
      pcave = 0.00d0
      do 20 i=1,nx
      do 19 j=1,ny
        pcmin =  min (pcmin, pc(i,j))
        pcmax =  max (pcmax, pc(i,j))
        pcave = pcave +      pc(i,j)
 19   continue
 20   continue
      pcave = pcave / AREAL(nx *  ny)

!
!     NCAR Code STARTS:
!     Force PLOTCHAR to use characters of the lowest quality.
!cNCAR      CALL lsc_pcseti ('QU - QUALITY FLAG',2)
!     initialize the drawing of the contour plot.
!cNCAR      CALL lsc_cprect (pc, nxdim, nx, ny, wkarx, NWKDIX * 4,
!cNCAR     ^     iwkarx, NWKDIX)
!     Draw the default background.
!cNCAR      CALL lsc_cpback (pc, wkarx, iwkarx)
!     Draw contour lines and labels.
!cNCAR      CALL lsc_cpcldr (pc, wkarx, iwkarx)
!     Add the informational label and the high/low labels.
!cNCAR      CALL lsc_cplbdr (pc, wkarx, iwkarx)
!     Compute and print statistics for the plot, label it, and put a
!     boundary line around the edge of the plotter frame.
!     cALL CAPSAP ('counts', TIME, IAMA, 0)
!     cALL LABTOP ('EXAMPLE 1-1',.017)
!     cALL BNDARY
!     NCAR Code ENDS:
!
!     ------------------------------------------------------------------
!     kclev1: >0 --> clevelin contains levels to be used;
!                    dots for index less than kclev2;
!                    solid for index greater or equal kclev2
!     kclev1: =0 --> first contour at clevelin(1) with
!                    next one up by clevelin(2), and so on and so on
!     kclev1: <0 --> rcontr to choose -kclev1 equally spaced values between
!                    clevelin(1) and clevelin(2)
!     clevelin:      array of contour levels; this is output if kclev1<0
!     kclev2:        separates dots from solid lines
!     ------------------------------------------------------------------
!      previous to Mar 93 usage:
!      kclev1 =-10
!      kclev2 = 2
!      pcmin = pcmin + (pcmax-pcmin) / 10.
!      clevelin(1) = pcmin
!      clevelin(2) = pcmax-pcmin


      pcmin = pcmax / 100.
      clevelin(1) = pcmin
      clevelin(2) = pcmax*0.91

      ixmin = 1
      ixmax = nx
      jymin = 1
      jymax = ny
      ixstep= 1
      jystep= 1

      if (iside .EQ. 1) then
      CALL lsc_EZrcon(ix1(iside),ix2(iside),     iy1,        iy2,       &
     &                 xminy,     xmaxy,   yminy,      ymaxy,           &
     &  kclev1,clevelin,kclev2,                                         &
     &  pc,                                                             &
     &  nxdim,                                                          &
     &  xAry, ixmin, ixmax, ixstep,                                     &
     &  yAry, jymin, jymax, jystep )
      else if (iside .EQ. 2) then
      CALL lsc_EZrcon(ix1(iside),ix2(iside),     iy1,        iy2,       &
     &                 xminy,     xmaxy,   yminy,      ymaxy,           &
     &  kclev1,clevelin,kclev2,                                         &
     &  pc,                                                             &
     &  nxdim,                                                          &
     &  xGiv, ixmin, ixmax, ixstep,                                     &
     &  yGiv, jymin, jymax, jystep )
      else
         continue
      endif
!
      if (iside .EQ. 1)                                                 &
     &      write(ChStr,'(''min: '',1pe9.2,'' w/cm^2 ...$'')')pcmin
      if (iside .EQ. 2)                                                 &
     &      write(ChStr,'(''min: '',1pe9.2,'' w/cm^3 ...$'')')pcmin
      CALL lsc_EZwrit(ix1(iside), 750, ChStr, 0,0)
!
      write(ChStr,'(''max: '',1pe9.2,''  $'')')pcmax
      CALL lsc_EZwrit(ix1(iside), 725, ChStr, 0,0)
!
      write(ChStr,'(''ave: '', 1pe9.2,'' $'')')pcave
      CALL lsc_EZwrit(ix1(iside), 700, ChStr, 0,0)

      if (iside .EQ. 1) then
      write(ChStr,'(''Xmx/mn:'',2(1x,1pe9.2),''$'')')                   &
     &               xAry(ixmax),xary(1)
      CALL lsc_EZwrit(ix1(iside), 675, ChStr, 0,0)
      write(ChStr,'(''Ymx/mn:'',2(1x,1pe9.2),''$'')')                   &
     &               yAry(jymax),yary(1)
      CALL lsc_EZwrit(ix1(iside), 650, ChStr, 0,0)
      else if (iside .EQ. 2 ) then
      write(ChStr,'(''Xmx/mn:'',2(1x,1pe9.2),''$'')')                   &
     &               xGiv(ixmax),xGiv(1)
      CALL lsc_EZwrit(ix1(iside), 675, ChStr, 0,0)
      write(ChStr,'(''Ymx/mn:'',2(1x,1pe9.2),''$'')')                   &
     &               yGiv(jymax),yGiv(1)
      CALL lsc_EZwrit(ix1(iside), 650, ChStr, 0,0)

      endif

!
      return
      END
!
!     +---------------------------------------------------------------+
!
      SUBROUTINE lsc_xSigGen(Emev, Ephoton, mu, sigma)
use iso_c_binding, only : fp => c_double
implicit none

      REAL   EMev, Ephoton, mu, sigma, T0, k
      real   signorm, Mevpermcsq
      COMMON / sigbrcom0 / signorm, Mevpermcsq
!
!     converts to / from Mev to Koch, Motz units
!     input energies in Mev, angle in radians
!     return sigma in cm**2 / sR / MeV
!
      T0 = Emev / Mevpermcsq
      k = Ephoton / Mevpermcsq
      CALL lsc_xSigBrem(T0, k, mu, sigma)
      sigma = sigma / Mevpermcsq
      return
      END

!     +---------------------------------------------------------------+

      SUBROUTINE lsc_xBremSetup(Z)
use iso_c_binding, only : fp => c_double
use physconst_mod, only : pi, zmtocm, zm2tocm2
use physconst_mod, only : zel  ! electron charge in statcoulomb
use physconst_mod, only : vc   ! speed light in m/s
use physconst_mod, only : me_Kg
implicit none

      REAL   Z
      real   pie_local, r0, ECHARGE, EMASS, CLIGHT, HBAR, alpha
      PARAMETER(ECHARGE = 4.8e-10, EMASS = 9.1e-28,                     &
     &     CLIGHT = 3.e10, HBAR = 1.05e-27)
      real   signorm, Mevpermcsq
      COMMON / sigbrcom0 / signorm, Mevpermcsq
!
      MevPerMcsq = .511
!      pie_local = 4. * atan(1.)    ! fmp - this is pi
      alpha = zel * zel / (HBAR * vc * zmtocm)
!     r0 = zel * zel / (EMASS * vc*vc * zm2tocm2)
      r0 = zel * zel / (me_Kg * vc*vc * zm2tocm2)
!     signorm = alpha * r0 * r0 * Z * Z / ( 8. * pie_local )
      signorm = alpha * r0 * r0 * Z * Z / ( 8.0_fp * pi )
      return
      END

!     +-----------------------------------------------------------------+

      SUBROUTINE lsc_xSigBrem(T0, k, mu, sigma)
use iso_c_binding, only : fp => c_double
implicit none

!     Reference:
!     H. W. Koch and J. W. Motz,
!     "Bremsstrahlung Cross-Section Formulas and Related Data,"
!     \RMP{31}{920}{59}.
!     Formula 2BN on p 924
!     '2' denotes differential in photon energy and angle
!     'B' denotes Born approximation
!     'N' denotes no screening
!     T0       INPUT  initial k-e of electron in a collision; m_0 c^2 units
!     k        INPUT  energy of emitted photon;               m_0 c^2 units
!     mu       INPUT  cos theta_0; angle between initial e-momentum and emitted
!                     photon
!     sigma    OUTPUT differential cross section
!                     w.r.t. photon energy(k) & solid angle(\Omega_k) in
!                     cm^2 per atom per incident electron
!     signorm  Z^2 r_0^ /(8 \pi 137); units of cm^2

      REAL   sigma, k, mu, T0
      REAL   ksq, L, Delta0, p0, p0sq, eps, epsQ, Q, p, E, E0, T,       &
     &     Qsq, beta, beta0, Cth0, Sth0, EE0, term(13), E0sq,           &
     &     Esq, Delta0sq, Delta0_4, p0sqDsq, p0sqD4, Sth0sq
      real   signorm, Mevpermcsq
      COMMON / sigbrcom0 / signorm, Mevpermcsq
!
!     computation of much used factors
      ksq = k * k
      p0 = sqrt(T0 * (T0 + 2.d0))
      p0sq = p0 * p0
      E0 = 1.d0 + T0
      E0sq = E0 * E0
      E = E0 - k
      Esq = E * E
      EE0 = E * E0
      T = E - 1.d0
      p =  sqrt(T * (T + 2.))
      beta0 = sqrt(1.d0 - 1.d0 / (E0 * E0))
      beta = sqrt(1.d0 - 1.d0 / (E * E))
      Cth0 = mu
      Qsq = p0sq + ksq - 2. * p0 * k * Cth0
      Q = sqrt(Qsq)
      Sth0 = sqrt(1.d0 - mu * mu)
      Sth0sq = Sth0 * Sth0
      Delta0 = E0 - p0 * cth0
      Delta0sq = Delta0 * Delta0
      Delta0_4 = Delta0sq * Delta0sq
      p0sqDsq = p0sq * Delta0sq
      p0sqD4 = p0sq * Delta0_4
      eps = log((E + p) / (E - p))
      epsQ = log((Q + p) / (Q - p))
      L = log((EE0 - 1.d0 + p * p0) /                            &
     & (EE0 - 1.d0 - p * p0))
      term(1) = 8. * Sth0sq * ( 2. * E0sq + 1.d0) / p0sqD4
      term(2) = - 2. * (5. * E0sq + 2. * EE0 + 3.d0 ) / p0sqDsq
      term(3) = - 2. * (p0sq - ksq) / (Qsq * Delta0sq)
      term(4) = 4. * E / (p0sq * Delta0)
      term(5) = 4. * E0 * Sth0sq * (3. * k - p0sq * E) / p0sqD4
      term(6) = 4. * E0sq * (E0sq + Esq) / p0sqDsq
      term(7) = 2. * (1.d0 - (7. * E0sq - 3.* EE0 + Esq)) /      &
     &          p0sqDsq
      term(8) = 2. * k * (E0sq + EE0 - 1.d0) / (p0sq * Delta0)
      term(9) = - 4. * eps / (p * Delta0)
      term(10) = epsQ / (p * Q)
      term(11) = 4. / Delta0sq - 6. * k / Delta0 -                      &
     &      2. * k * (p0sq - ksq) / (Qsq * Delta0)
      term(12) = L / (p * p0)
      term(13) = signorm * p / (k * p0)
      sigma = term(13) * (                                              &
     &     term(1) + term(2) + term(3) + term(4)                        &
     &     + term(12) * ( term(5) + term(6) + term(7) + term(8) )       &
     &     + term(9) + term(10) * term(11)                              &
     &     )
      return
      END
!
!     graphs          --------------------------------------------------
!     End Xray Camera --------------------------------------------------
!
!
!     . ----------------------------------------------------------------
!     . ----------------------------------------------------------------
!     . Absorbers and Scintilator Treated Here -------------------------
!     . ----------------------------------------------------------------
!     . ----------------------------------------------------------------
!     .
!     Copyright D. W. Ignat, S. von Goeler, J. E. Stevens,
!     P. G. Roney, E. J. Valeo
!     Princeton University, Plasma Physics Laboratory, 1992, 1993
!
!     ------------------------------------------------------------------
!
      SUBROUTINE lsc_GetXmnFac(nEs, EbyMeV, XmnFac, dcu, FC )
use iso_c_binding, only : fp => c_double
use lsc_xraymod
implicit none

!#include "lsc_xparams.h"
      INTEGER nEs
      REAL    EbyMeV(nEs), XmnFac(nEs)
      REAL DCU
      CHARACTER*2 FC
      CHARACTER*32 ErrorTxt
!                       ! DCU: thickness of foil (cm)
!                       ! FC:  FoilCode: CU,AG,TA,MO, or 00 for no foil
      INTEGER  nds, ndsDIM
      PARAMETER ( NDS = 12, ndsDIM = 12)
      REAL    dcm(nds)
      INTEGER i
      INTEGER iVACUUM,  iCARBON,iALUMINM,   iIRON, iNICKEL,             &
     &                  iCOPPER, iIODINE, iCESIUM, iCSI,                &
     &                  iTANTALUM, iMOLYBDENUM, iSILVER
      PARAMETER (iVACUUM = 1, iCARBON = 2, iALUMINM = 3)
      PARAMETER (iIRON = 4, iNICKEL = 5)
      PARAMETER (iCOPPER = 6, iIODINE = 7 , iCESIUM = 8, iCSI = 9)
      PARAMETER (iTANTALUM = 10, iMOLYBDENUM = 11, iSILVER = 12)
      INTEGER nSigInt
      PARAMETER (nSigInt = NENDIM)
      REAL    Ekv(nSigInt)
      real    rhoCsCSI, rhoICsI, rhoCarbon
      DATA    rhoCsCSI, rhoICsI, rhoCarbon /                            &
     &           2.307,   2.203,       1.0 /
!                                            !must be looked up
!     k Edges at  33.169,       35.985,  kV
!     These Data Give the Shield Always in Place in the Actual Experiment:
!     DATA    dcm(iVACUUM), dcm(iCARBON), dcm(iALUMINM) /
!    ^                0.0 ,         0.2 ,     0.37489   /
!     DATA    dcm(iIRON)  , dcm(iNICKEL), dcm(iCOPPER)  /
!    ^            0.0031  ,       0.0137,         0.0   /
!     DATA    dcm(iIODINE), dcm(iCESIUM), dcm(iCSI)     /
!    ^              0.0   ,         0.0 ,     0.0297    /
!     DATA  dcm(iTANTALUM),dcm(iMOLYBDENUM),dcm(iSILVER)/
!    ^              0.0   ,         0.0 ,         0.0   /
!
!     If we wanted to have phony foils to look at Photon Temperature
!     at low energy then SvG recomments the following:
!        10 kv 0.014 cm Al; 20 kv 0.109 cm Al;
!        30 kv 0.335 cm Al; 40 kv 0.650 cm Al.
!

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!     put in thicknesses of magnetic screen over the xray tube --
!     identical to data statements above
      dcm(iVACUUM)   = 0.0
      dcm(iCARBON)   = 0.2
      dcm(iALUMINM)  = 0.37489
      dcm(iIRON)     = 0.0031
      dcm(iNICKEL)   = 0.0137
      dcm(iMOLYBDENUM)= 0.0006
!
!     put in active material on the screen
      dcm(iCSI)      = 0.0297
!
!     zero out thicknesses for optional materials -
      dcm(iCOPPER)    = 0.0
      dcm(iIODINE)    = 0.0
      dcm(iCESIUM)    = 0.0
      dcm(iTANTALUM)  = 0.0
      dcm(iSILVER)    = 0.0
!     .                                 determine subscript into dcm
!     .                                 corresponding to foil code, and
!     .                                 put foil thickness into place;
!     .
!     .                                 If PhonyFoil, ie  FC='PF', then
!     .                                 zero out the default thickness for
!     .                                 everything, and put in an Aluminum foil
!     .                                 of the given thickness (SvGoeler Dec93)
      if (FC.EQ.'CU') then
        dcm(iCOPPER) = dcu
      else if (FC.EQ.'TA') then
        dcm(iTANTALUM) = dcu
      else if (FC.EQ.'MO') then
        dcm(iMOLYBDENUM) = dcu + dcm(iMOLYBDENUM)
      else if (FC.EQ.'AG') then
        dcm(iSILVER) = dcu
!     .                                Phony Foil Here
      else if (FC.EQ.'PF') then
        do i=1,NDS
           dcm(i)=0.00
        enddo
!     .                                put back in thickness of the
!     .                                aluminum phony foil 'PF'
           dcm(iALUMINM) = dcu
!     .                                put back in active screen material
           dcm(iCSI)      = 0.0297
!     .
      else if (FC.EQ.'00') then
        continue
      else
        write(ErrorTxt,                                                 &
     &      '(a2,'' !! Bad foil code in GetXmnFac'')') FC
        CALL lsc_LSCwarn (ErrorTxt)
      endif


!     .                                 Convert units of mc^2? or MeV to keV
      do 10 i = 1,nEs
        Ekv(i)  =  EbyMeV(i)*1000.
 10   continue

        CALL lsc_GetXmiss (nEs, Ekv, XmnFac,  nds, dcm )

!
!cDBG        write(4,'(('' Ekv, Xmn: '', 2(1pe10.3,1x),/))')
!cDBG     ^  (Ekv(i), XmnFac(i), i=1,nEs )
      return
      END
!
!     ------------------------------------------------------------------
!
      SUBROUTINE lsc_Absorb12(Energy, iNAME,                            &
     &                   BarnPerA, MuByRho, nSigma, iError, ErrorTxt )
!     Reference:
!     X ray crosss sections based on W. H. McMaster, el al, UCRL 50174
!     Sec. II rev. 1
!     May 1969
!     \bibitem{mcmaster}
!     W. H. McMaster, \etal , UCRL 50174,
!     Sec. II rev. 1,  May 1969.

!     Copyright D. W. Ignat, S. von Goeler, J. E. Stevens,
!     P. G. Roney, E. J. Valeo
!     Princeton University, Plasma Physics Laboratory, 1992, 1993
!
use iso_c_binding, only : fp => c_double
use physconst_mod, only : zero, one
implicit none

      INTEGER iNAME, iError
      CHARACTER*32  BLANK
      CHARACTER*(*) ErrorTxt
      REAL    Energy, BarnPerA, MuByRho, nSigma
      INTEGER iVACUUM,  iCARBON,iALUMINM,   iIRON, iNICKEL,             &
     &                  iCOPPER, iIODINE, iCESIUM, iCSI,                &
     &                  iTANTALUM, iMOLYBDENUM, iSILVER
      PARAMETER (iVACUUM = 1, iCARBON = 2, iALUMINM = 3)
      PARAMETER (iIRON = 4, iNICKEL = 5)
      PARAMETER (iCOPPER = 6, iIODINE = 7 , iCESIUM = 8, iCSI = 9)
      PARAMETER (iTANTALUM = 10, iMOLYBDENUM = 11, iSILVER = 12)
      INTEGER iLshell, iKshell,  iCoher ,  iInCoh
      INTEGER nLshells, nLatoms
      PARAMETER (iLshell = 1, iKshell = 2,  iCoher = 3, iInCoh = 4)
      PARAMETER (nLshells = 4, nLatoms = 12)
!     nLshells is a misnomer: 1 refers to L interactions below k edge
!     .                       2 refers to K interactions above k edge
!     .                       3 refers to coherent scattering
!     .                       4 refers to incoherent scattering
!
!     If M shell interactions would be considered, then
!     nMshells would be 5: M, L, K, coherent, incoherent
!     Note, we are ignoring M shells in Cs and I
!
!
      real    kEdge(nLatoms), grPcc(nLatoms), AtmWt(nLatoms)
      REAL    Avogadro

      INTEGER iShell, i,j
      real    lnEnergy, lnSigma, epower(0:3)

      real    aL(4,nLshells,nLatoms)
      DATA    BLANK /'                                '/
!.    .              0123456789 123456789 123456789 12
      DATA    Avogadro / 6.02252d-1 /
!     DATA    Avogadro / 6.02252e-01 /  =  6.02 10^23 x 10^-24 for barns
      DATA     kEdge(iALUMINM), kEdge(iIRON), kEdge(iCOPPER)/           &
     &                   1.560,        7.112,          8.979/
      DATA     kEdge(iIODINE ), kEdge(iCESIUM),kEdge(iNICKEL)/          &
     &                  33.169,       35.985,          8.333 /
      DATA kEdge(iTANTALUM),kEdge(iMOLYBDENUM),kEdge(iSILVER)/          &
     &                  67.414,       19.999,         25.514 /

      DATA     grPcc(iALUMINM), grPcc(iIRON), grPcc(iCOPPER)/           &
     &                 2.702  ,       7.860 ,         8.940 /
      DATA     grPcc(iIODINE ), grPcc(iCESIUM),grPcc(iNICKEL)/          &
     &                 4.94   ,      1.873   ,        8.90   /
      DATA grPcc(iTANTALUM),grPcc(iMOLYBDENUM),grPcc(iSILVER) /         &
     &                16.600,       10.220,          10.500 /

      DATA     AtmWt(iALUMINM), AtmWt(iIRON), AtmWt(iCOPPER)/           &
     &                26.970  ,      55.850 ,        63.540 /
      DATA     AtmWt(iIODINE ), AtmWt(iCESIUM),AtmWt(iNICKEL)/          &
     &                126.90  ,      132.910,        58.690  /
      DATA AtmWt(iTANTALUM),AtmWt(iMOLYBDENUM),AtmWt(iSILVER)/          &
     &                180.950 ,       95.950,       107.880  /

      DATA     kEdge(iCARBON), grPcc(iCARBON) , AtmWt(iCARBON)  /       &
     &                  0.00 ,         1.580  ,      12.010     /

      DATA   ((aL(i,j      ,iVACUUM),i = 1,4),j = 1,nLshells)  /        &
     &        16*0.0                                       /

      DATA   ((aL(i,j,iCARBON ),i = 1,4),j = 1,nLshells)          /     &
     &         0.0,        0.0,        0.0,        0.0,                 &
     &         1.06879e1, -2.71400,   -2.00530e-1, 2.07248e-2,          &
     &         3.10861,   -2.60580e-1,-2.71974e-1, 1.35181e-2,          &
     &        -9.82878e-1, 1.46693,   -2.93743e-1, 1.56005e-2 /

      DATA   ((aL(i,j,iALUMINM),i = 1,4),j = 1,nLshells)          /     &
     &         1.08711E1, -2.77860,    1.75853e-1, 0.0,                 &
     &         1.31738e1, -2.18203,   -2.5896e-1,  2.2284e-2,           &
     &         4.51995  ,  1.40549e-1,-3.5244e-1,  1.9369e-2,           &
     &        -4.39322e-1, 1.30867 ,  -2.22648e-1, 7.54210e-3 /

      DATA   ((aL(i,j,iIRON   ),i = 1,4),j = 1,nLshells)          /     &
     &         1.36696E1, -2.39195,   -1.37648E-1, 0.,                  &
     &         1.43456E1, -1.23491,   -4.18785E-1, 3.21662E-2,          &
     &         5.93292,    2.25038E-1,-3.61748E-1, 1.93024E-2,          &
     &        -3.42379E-1 ,1.57245,   -2.53198E-1, 9.85822E-3 /

      DATA   ((aL(i,j,iNICKEL ),i = 1,4),j = 1,nLshells)          /     &
     &         1.39848e1, -2.48080,   -8.88115e-2, 0.0,                 &
     &         1.42388e1, -9.67736e-1,-4.78070e-1, 3.66138e-2,          &
     &         6.09204  , 2.52277e-1, -3.66568e-1, 1.96586e-2,          &
     &        -5.0436e-1, 1.70040,    -2.76443e-1, 1.12628e-2 /

      DATA   ((aL(i,j,iCOPPER ),i = 1,4),j = 1,nLshells)           /    &
     &        14.2439,    -2.58677,   -6.67398E-2, 0.    ,              &
     &        14.5808,    -1.18375,   -0.41385,    3.12088E-2 ,         &
     &         6.17739,    2.73123E-1,-3.7236E-1,  2.01638E-2 ,         &
     &        -5.7021E-1,  1.75042,   -2.84555E-1, 1.1693E-2  /

      DATA   ((aL(i,j,iIODINE ),i = 1,4),j = 1,nLshells)           /    &
     &         1.64086e1, -2.48214,   -5.07179e-2, 0.0,                 &
     &         1.21075e1,  1.43635,   -8.82038e-1, 6.03575e-2,          &
     &         7.27415,    3.77223e-1,-3.69728e-1, 1.86280e-2,          &
     &        -4.04420e-2, 1.65596,   -2.510676,   9.04874e-3   /

      DATA   ((aL(i,j,iCESIUM ),i = 1,4),j = 1,nLshells)           /    &
     &         1.65418e1, -2.46363,   -5.42849e-2, 0.0,                 &
     &         1.13757e1,  1.94161,   -9.83232e-1, 6.71986e-2,          &
     &         7.33490,    3.76825e-1,-3.65713e-1, 1.81843e-2,          &
     &         1.84861e-1, 1.50030,   -2.13333e-1, 6.24264e-3  /

      DATA  ((aL(i,j,iTANTALUM),i = 1,4),j = 1,nLshells)           /    &
     &         1.72410E1, -2.30313,   -5.91006E-2, 0.0,                 &
     &         8.65271,    3.73117,   -1.26359,    8.23539E-2,          &
     &         7.94534,    3.87299E-1,-3.47926E-1, 1.63299E-2,          &
     &         1.96871E-1, 1.50623,   -1.91396E-1, 3.70889E-3 /

      DATA  ((aL(i,j,iSILVER),i = 1,4),j = 1,nLshells)             /    &
     &         1.56869E+1,-2.22636,   -1.12223E-1, 0.0,                 &
     &         1.33926E+1, 4.41380E-1,-6.93711E-1, 4.82085E-2,          &
     &         7.06446,    3.63456E-1,-3.73597E-1, 1.92478E-2,          &
     &        -1.66475E-1, 1.65794,   -2.48740E-1, 8.66218E-3 /

      DATA  ((aL(i,j,iMOLYBDENUM),i = 1,4),j = 1,nLshells)         /    &
     &         1.53494E+1,-2.26646,   -1.16881E-1, 0.0,                 &
     &         1.39853E+1,-1.17426E-1,-5.91094E-1, 4.17843E-2,          &
     &         6.84600,    3.02797E-1,-3.51131E-1, 1.74403E-2,          &
     &        -5.62860E-2, 1.55778,   -2.33341E-1, 7.85506E-3 /
!
!     ............................................. Execution begins here
!
      REAL             TENp3
      DATA             TENp3 /                    &
     &     1.0d3/

      iError  =  0
      ErrorTxt= BLANK

      if (iNAME .EQ. 1) then
!     .                                 The practice of including 'data' for
!     .                                 vacuum came from Jim Stevens, and we
!     .                                 are not sure why.  Here it is bypassed
!%0.
        BarnPerA  =  ZERO
        MuByRho   =  ZERO
        nSigma    =  ZERO
        return
      endif

      if (Energy .LT. ONE ) then
        write(ErrorTxt,'('' invalid low energy: '',1pe10.3)')  Energy
        BarnPerA  =  ZERO
        MuByRho   =  ZERO
        nSigma    =  ZERO
        iError  =  1
        return
      endif

      if (Energy .GT. TENp3 ) then
        write(ErrorTxt,'('' invalid high energy: '',1pe10.3)') Energy
        BarnPerA  =  ZERO
        MuByRho   =  ZERO
        nSigma    =  ZERO
        iError  =  2
        return
      endif

      if (iNAME .GT. 12) then
        write(ErrorTxt,'('' invalid atom: '', i4)')           iNAME
        BarnPerA  =  ZERO
        MuByRho   =  ZERO
        nSigma    =  ZERO
        iError  =  3
        return
      endif


!     !                                 Major branch based on number of shells
!     !                                 (branch not used; M shell ignored now)

      IF (iNAME .LE. nLatoms) then
!     .                                 Use this code for L-shell atoms
        if (Energy .LT. kEdge(iNAME)) then
          iShell  =  1
        else
          iShell  =  2
        endif
!
        lnEnergy  =   log(Energy)
        epower(0)  =  ONE
        do 5 i = 1,3
          epower(i) = epower(i-1)*lnEnergy
 5      continue

          lnSigma  =  ZERO
        do 10 i = 1,4
          lnSigma  =  lnSigma + aL(i,iShell,iNAME)*epower(i-1)
  10    continue

        BarnPerA  =  exp(lnSigma)

        do 40 iShell = nLshells-1, nLshells
            lnSigma  =  ZERO
          do 30 i = 1,4
            lnSigma  =  lnSigma + aL(i,iShell,iNAME)*epower(i-1)
  30      continue
        BarnPerA  =  BarnPerA + exp(lnSigma)
  40    continue

!     .                                 Use this code for M-shell atoms
!     .                                 if and when this is added
      ELSE
        write(ErrorTxt,'('' invalid atom: '', i4)') iNAME
        BarnPerA  =  ZERO
        MuByRho   =  ZERO
        nSigma    =  ZERO
        iError  =  4
        return
      ENDIF
!     .                                 Finally, get N rho / AtWt * Sigma



      MuByRho   =  BarnPerA * Avogadro / AtmWt(iNAME)
      nSigma    =  MuByRho  * grPcc(iNAME)

      return
      END
!
!     ------------------------------------------------------------------
!
      SUBROUTINE lsc_GetXmiss (nSigInt, Eaxis, Xmiss,                   &
     &                     nds, dcm )
!     Copyright D. W. Ignat, S. von Goeler, J. E. Stevens,
!     P. G. Roney, E. J. Valeo
!     Princeton University, Plasma Physics Laboratory, 1992, 1993
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER nSigInt
      INTEGER nds
      CHARACTER*32 ErrorTxt
      REAL    dcm(nds)
      REAL    Eaxis(nSigInt), Xmiss(nSigInt)
      INTEGER iVACUUM,  iCARBON,iALUMINM,   iIRON, iNICKEL,             &
     &                  iCOPPER, iIODINE, iCESIUM, iCSI,                &
     &                  iTANTALUM, iMOLYBDENUM, iSILVER
      PARAMETER (iVACUUM = 1, iCARBON = 2, iALUMINM = 3)
      PARAMETER (iIRON = 4, iNICKEL = 5)
      PARAMETER (iCOPPER = 6, iIODINE = 7 , iCESIUM = 8, iCSI = 9)
      PARAMETER (iTANTALUM = 10, iMOLYBDENUM = 11, iSILVER = 12)

      INTEGER i, j, iError
      REAL    BarnPerA, MuByRho, nSigma
      REAL    MuRhoD, MuRhoCsI
      real    rhoCsCSI, rhoICsI, rhoCarbon
      DATA    rhoCsCSI, rhoICsI, rhoCarbon /                            &
     &           2.307,   2.203,       1.0 /
!                                            !must be looked up
!     k Edges at  33.169,       35.985,  kV
!
      REAL       ZERO,         ONE,          TENp3
      DATA       ZERO  ,       ONE   ,       TENp3 /                    &
     &     0.0d0, 1.0d0, 1.0d3/

      do 70 i = 1,nSigInt
!
      MuRhoD    =  ZERO
      MuRhoCsI  =  ZERO

      if (Eaxis(i) .GT. TENp3 ) then
        Xmiss(i) = ONE
      else if(Eaxis(i) .LT. ONE ) then
        Xmiss(i) = ZERO
      else
!
        do 50 j = 1,nds
          if (dcm(j) .GT. ZERO .and. j .ne. iCsI) then
            CALL lsc_Absorb12(Eaxis(i), j,                              &
     &                   BarnPerA, MuByRho, nSigma, iError, ErrorTxt)
!
            if (j .EQ. iCARBON) then
              MuRhoD  =  MuRhoD + MuByRho*rhoCarbon*dcm(j)
            else
              MuRhoD  =  MuRhoD + nSigma*dcm(j)
            endif
            if (iError .GE. 1) CALL lsc_LSCwarn(ErrorTxt)
          endif
!
          if (dcm(j) .GT. ZERO .and. j .EQ. iCSI) then
            CALL lsc_Absorb12(Eaxis(i), iIODINE,                        &
     &                   BarnPerA, MuByRho, nSigma, iError, ErrorTxt)
!
            MuRhoCsI  =  MuRhoCsI + MuByRho*rhoICsI*dcm(j)
!
            CALL lsc_Absorb12(Eaxis(i), iCESIUM,                        &
     &                   BarnPerA, MuByRho, nSigma, iError, ErrorTxt)
!
            MuRhoCsI  =  MuRhoCsI + MuByRho*rhoCsCsI*dcm(j)
!
            if (iError .GE. 1) CALL lsc_LSCwarn(ErrorTxt)
          endif
!
 50     continue
!
!
        if (MuRhoCsI .EQ. ZERO) then
          CALL lsc_LSCwarn(' no phosphor screen found; 1e-5 assumed')
          MuRhoCsI  =  1.0d-5
        endif

        Xmiss(i)  =  exp ( - MuRhoD )*(ONE - exp(-MuRhoCsI) )
      endif
 70   continue
!
      return
      END

!     ------------------------------------------------------------------
!     . ----------------------------------------------------------------
!     . Absorbers and Scintilators End Here ----------------------------
!     . ----------------------------------------------------------------



!                                                                      |
!     -----------------------------------------------------------------|
!                                                                      |
      SUBROUTINE lsc_XRopen(GivnUnit, FileFrstNa)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER GivnUnit, SomeUnit, linepointer
      CHARACTER*1 PERIOD, SPACE
      CHARACTER*8  FileFrstNa, FileScndNa
      CHARACTER*8  mmddyy, hhmmss
      CHARACTER*17 FileCompNa
      CHARACTER*75 SomeString

      COMMON / XRrwcom / linepointer,SomeUnit
      DATA        PERIOD, SPACE / '.' , ' ' /

      CHARACTER*8  zzdate
      CHARACTER*10 zztime
      CHARACTER*5  zzzone
!
!     zzdate  ccyymmdd    cc=century yy=year mm=month dd=day
!     zztime  hhmmss.ttt  hh=hour mm=minute ss=second ttt=millisecond
!     zzzone  Shhmm       S=sign (-=west, +=east) hh=hour mm=minute
!
      INTEGER ivals(8)
!                             EXAMPLE
!     ivals ( 1 ) year        "2000"
!     ivals ( 2 ) month       "10"
!     ivals ( 3 ) day         "12"
!     ivals ( 4 ) t - UTC     "-0400"
!     ivals ( 5 ) hour        "8"
!     ivals ( 6 ) minute      "21"
!     ivals ( 7 ) second      "23"
!     ivals ( 8 ) millisecond "621"
!
      call date_and_time(zzdate,zztime,zzzone,ivals)
      write(mmddyy,'(a)') zzdate(5:8)//zzdate(3:4)
      write(hhmmss,'(a)') zztime(1:6)

      SomeUnit = GivnUnit

      if (FileFrstNa .EQ. 'XRinput') then
        open (SomeUnit,file='XRinput',                                  &
     &        status='unknown')
        linepointer=0
      else
        FileScndNa = hhmmss
        open (SomeUnit,file=FileFrstNa//PERIOD//FileScndNa,             &
     &        status='unknown',err=1300 )
        linepointer=0
        SomeString=FileFrstNa//SPACE//mmddyy//SPACE//hhmmss
        write(SomeUnit,'(a75)') SomeString
        linepointer=1
      endif
 1300 continue
!c1300 continue need something real here
      return

      ENTRY lsc_XRclos(GivnUnit)
      close(GivnUnit)
      return
      END

      SUBROUTINE lsc_XRdscw(xarray, yarray, npoints, SomeString)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER SomeUnit, linepointer
      INTEGER i,npoints
      REAL    xarray(npoints), yarray(npoints)
      CHARACTER*(*) SomeString
      COMMON / XRrwcom / linepointer,SomeUnit
      write(SomeUnit,'(a75)') SomeString
        linepointer=linepointer+1
      do 10 i=1,npoints
        write(SomeUnit,'(1x,1pe10.3,1x,1pe10.3)')xarray(i),yarray(i)
        linepointer=linepointer+1
 10   continue
      return
      END
!
!      SUBROUTINE XRdscr(xarray, yarray, npoints, SomeString)
!      IMPLICIT NONE
!      INTEGER SomeUnit, linepointer
!      INTEGER i,npoints
!      REAL    xarray(npoints), yarray(npoints)
!      CHARACTER*75 SomeString
!      COMMON / XRrwcom / linepointer,SomeUnit
!
!      if(linepointer .eq. 0) then
!        read(SomeUnit,'(a75)') SomeString
!        linepointer=linepointer+1
!      endif
!
!      read(SomeUnit,'(a75)') SomeString
!        linepointer=linepointer+1
!      do 10 i=1,npoints
!        read (SomeUnit,'(1x,1pe10.3,1x,1pe10.3)')xarray(i),yarray(i)
!        linepointer=linepointer+1
! 10   continue
!
!      return
!      END
#endif

!subroutine lsc_BrambJES (fghz, phaseDeg, nrays, nparalls, amplitds,  &
!     &                     nparaMin, nslices, coupler, turnnegs)
subroutine lsc_BrambJES (fghz, nparalls, amplitds, nparaMin, coupler, turnnegs)
!     J. E. Stevens' Brambilla coupling code
!
!     Copyright 10/22/80     by      J. E. Stevens
!	edit:19-feb-87 [Phyllis Roney] convert graphics to sg for vax
!       (and other fortran 77 modifications).
!
!       Compute spectrum of lower hybrid antenna,
!       power flow, and reflection in guides
!
!     Copyright November 1991 by
!     J. Stevens, S. Bernabei, D. Ignat, E. Valeo, S. Jardin
!     with modifications by Ignat to be called by LSC/TSC.
! Modified by Budny, McCune, and Gao (May-1998) to include the HL1M antenna
! The height of the antenna had been set to a = 5.8cm previously.
!
! Modified by S. Ding (July 2009) to include the EAST antenna
!
! Modified by S. Ding (December 2010) to include the new 4.6GHz antenna for EAST (option 10)
!     -
!     arguments developed for LSC/TSC
!           phaseDeg    phase in degrees between waveguides               INPUT
!           nrays       number of nparallels to return to make spectrum   INPUT
!           nparalls    array of n_par values from low to high           OUTPUT
!           amplitds    array of amplitudes corresponding to those n_par OUTPUT
!           nparaMin    minimum value of abs(n_par) of interest, say 1.5  INPUT
!           nslices     number of slices taken on spectrum in -8 to 8     INPUT
!           turnnegs    negative going spectral elements are turned to    INPUT
!                       positive if this is 1 and if phaseDeg is different
!                       from 180 or 0 by more than 20 degrees
!
!
!     input variables for original code (JES)
!
!           nguide      no. of waveguides
!                       .lt. 0 => phase of total b field fixed
!           nd          total no. of integration steps
!           delphi      phase diff. in degrees to step
!                       .lt. 14  => use phase in input for guide
!           f0          rf frequency in MHz
!           epsr        relative dielectric constant of material in guides
!                       .lt. 0 => don't include guide impedance in calc.
!           x0          antenna to plasma spacing in cm
!           gradn       density gradient at plasma edge in cm-4
!                       .lt. 0 => compute for +4 orders of magnitude
!           a(i)        waveguide height in cm
!                       This had been set to 5.8cm for all couplers.
!                       Changed to 7.6cm for ToreS and to 10.8cm for HL1M
!                       in May, 1998
!                       10.92cm for EAST 2.45GHz antenna in July, 2009
!			5.00cm for 4.6GHz antenna in December, 2010
!           anmax       max. nz to use in integration
!           anplt       max. nz to use in plot
!           xprobe      rf probe position(cm) at which
!                       total e & b fields are calculated
!           edged       plasma density at the edge in cm-3
!                       .lt. 0 => compute +3 orders of magnatude
!           p0(i)       incident power in i'th guide in kw
!                       .le. 0. => passive guide
!           phi(i)      phase of electric field in i'th guide in degrees
!           az(i)       z distance to  one  edge of i'th guide in cm
!           bz(i)       z distance to other edge of i'th guide in cm
!
!     if not enough lines of az,bz,phi,p0 are specified, then
!     the program assumes az(n,1) is the septum width, bz(1)
!     is the guide to guide spacing, phi(1) is the relative
!     phase, and p0(1) is the power for the remaining guides.
!
!     The NumRay largest components are selected for later
!     ray tracing.

! The JET_LHCD launcher consists of 6 rows, by 8 columns of multijunctions.
! Each multijunction, consists at the interface with the plasma of 2 rows by 4
! collumbs of waveguides. The 4 waveguides in a row are phased, by fixed phase
! shifters in the waveguide, giving a 90 deg phase difference between the
! forward travelling wave in adjacent waveguides.
!
! In total the plasma launcher interface consists of 12 rows with
! 4X8 = 32 active waveguides in each. On either side of each row of waveguides,
! are 2 passive waveguides.
!
! The dimension of the wavguides are 72.1mm high by 9mm wide. The
! walls between the waveguides are 2mm wide. The depths of the passive
! waveguides is 3.6 cm for the two on the left of the row, and 5.9cm for the
! two on the right, as seen from the plasma towards the launcher. The fixed
! phasing of a multijunction is such that the phases as seen from the plasma
! towards the launcher is from left to right 0 -90 -180 90(-270) .
!
! The vertical position of the middle of the rows of waveguides, with
! respect to the midplane of the machine (z=0) i as follows:
!   row 1   407mm  row 2  333mm  row 3  259mm  row 4  185mm
!   row 5   111mm  row 6   37mm  row 7  -37mm  row 8  -111mm
!   row 9  -185mm  row10 -259mm  row11 -333mm  row12  -407mm
! The 90 deg phase difference between the wave travelling towards the
! plasma in adjacent waveguides only holds when the match to the plasma is
! good. Otherwise the properties of the multijunction has to be taken into
! account. Assuming 90 degrees is a good approximation in most relevant
! situations though, as we normally need to operate with a good match.
!
! The laucher is shaped at the plasma interface to fit the plasma, this
! means that the waveguides are cut at an agle different to 90 deg at the
! plasma interface. What influence this has on the launched spectra I don't
! know - we havent computed it, but I dont think it alters the spectra much as
! long as we have a good match.
!
! When we talk about a phasing of the launcher we refer to the phase
! difference between the leftmost waveguides in adjacent multijunctions. Most
! of our experiments have been carried out at 0 deg phasing which gives a
! n_\parallel peak at 1.84.  --- Morten Lennholm;  ml@jet.uk

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: ncpldim
  implicit none

!     Argument list first:
  character(len=8) ::  coupler
  integer :: turnnegs
      REAL(fp)                                                              &
     &        fghz,                                                     &
     &        phaseDeg, nparalls(nrays), amplitds(nrays),               &
     &                  nparaMin
!     Argument list ends.
!
!     These following are pretty much the original variables.
      INTEGER    DIM,      NGPBX,      NPDIM
      PARAMETER (DIM = 60, NGPBX = 32, NPDIM = 601)
      INTEGER    NGTdeV,   NGTors, NGSslow, NGJET, NGTFTR, NGHL1M
      PARAMETER (NGTdeV = 32, NGTors = 34, NGSslow = 16)
      PARAMETER (NGJET  = 32, NGTFTR = 32, NGHL1M = 12)

!     Three "user" couplers are put here as place holders.
!     The original entry for all was identical to TFTR.
!     The intent is that these entries can be modified
!     by the user.
!     Key variables are NGUSR1,   NGUSR2,   NGUSR3
!     Key names are       USR1SPEC, USR2SPEC, USR3SPEC
!
!cc      INTEGER    NGUSR1,      NGUSR2,      NGUSR3
!**************************************************************************
      INTEGER   NGEAST,      NGEAS2
!cc      PARAMETER (NGUSR1 = 32, NGUSR2 = 32, NGUSR3 =32)
      PARAMETER (NGEAST=37, NGEAS2=57)
!**************************************************************************
!     CHARACTER*20 filename
      CHARACTER*75 pline
      COMPLEX                                                           &
     &        e0(DIM),      r(DIM),    bt(DIM), et(DIM), bb(DIM),       &
     &        aa(DIM,DIM), ai(DIM,DIM),                                 &
     &        c0, c1, c2, c3,                                           &
     &        c5, c6, c7, c8, c9,                                       &
     &        xp, rho, enz, anx, anzi, denom
      INTEGER nplot, nguide, nd, nphase, nedgd, ngradn, nph, ngn
      INTEGER nloop, ned, nl, nmod, ngset, n, n2
      INTEGER idim, isw, iepsw, ifast, i, iskip, istop, idel, ix1, ix2
      INTEGER l, k, ll, ix, iy, ngpr1, ngpr2, npmin, npmax, j, j2
      INTEGER lsmx, ndif
      INTEGER NumRay, NumBeg
      INTEGER ii
      INTEGER G
      REAL(fp)                                                              &
     &        a1, a2, PI, delt, delphi, f0, epsr, x0, gradn, gradn0,    &
     &  a(NCPLDIM), anmax, anplt, anz, anzac, anx3, air, bir, aip, bip, &
     &        al0, w0, height,                                          &
     &        xprobe, edged, pmax, xmin, xmax,                          &
     &        edged0, dencr, pave, zmin, zmax,                          &
     &        twidth, wdc, w, beta, beta1, ratio,                       &
     &        alp, alph3, ph0, ph1, ph2, del, c10
      REAL(fp)                                                              &
     &        dnz, pplas, pint, area, rr, ri, etr, eti, btr, bti,       &
     &        arg, teff, pdave, epsrp, sign,                            &
     &        sum, smx, t1, danplt, asum1, asum2, asum3, fracc, t2, t
      real(fp)    az(NCPLDIM,DIM),  bz(NCPLDIM,DIM)
      REAL(fp)                                                              &
     &                                  b(DIM),                         &
     &        p0(DIM),  phi(DIM), phi0(DIM),  ei(DIM),                  &
     &        rm(DIM),  rp(DIM),  p0r(DIM),                             &
     &        etph(DIM),sp(NPDIM), s(NPDIM),    sn(NPDIM),              &
     &        etm(DIM), btm(DIM), btph(DIM),  phold(DIM)
      REAL(fp)                                                              &
     &        RelPha, MinNpar, TotPwr, MaxPwr
      DATA        a1,                a2,                 PI        /    &
     &    0.355028d0, 0.2588194d0, 3.1415926d0/
      DATA  idim, nplot, nd   /                                         &
     &        60,   301, 4000 /
      DATA      delphi,       f0      ,       epsr  /                   &
     &     0.d0, 4600.d0, 1.d0 /

      DATA       x0   ,       gradn   /                                 &
     &     0.d0, 1.0d12 /
!cc      DATA         a   ,       anmax ,       anplt  /
!cc     ^     5.8d0, 30.d0, 10.d0 /
!******************************************************************************
! make coupler heights an array to accomodate TORESUPRA & HL1M (Budny)
      DATA ( a(i), i = 1, NCPLDIM ) /                                   &
     & 5.8d0, 5.8d0, 5.8d0,                        &
     & 7.6d0, 5.8d0, 5.8d0, 5.8d0,          &
     &  10.80d0, 10.92d0, 5.00d0 /
!******************************************************************************
      DATA         anmax ,       anplt  /                               &
     &     30.d0, 10.d0 /

      REAL(fp)      R90     ,      R180     ,      R270
      DATA      R90     ,      R180     ,      R270      /              &
     &     90.0d0, 180.0d0, 270.0d0 /

      DATA       xprobe,       edged   /                                &
     &     0.0d0, 1.5d12 /
      DATA ( az(1,i), i = 1, NGPBX  )                               /   &
     &     00.00, 00.71,  1.42,  2.13,  3.15,  3.86,  4.57,  5.28,      &
     &      6.30,  7.01,  7.72,  8.43,  9.45, 10.16, 10.87, 11.58,      &
     &     12.60, 13.31, 14.02, 14.73, 15.57, 16.46, 17.17, 17.88,      &
     &     18.90, 19.61, 20.32, 21.03, 22.05, 22.76, 23.47, 24.18 /
      DATA ( bz(1,i), i = 1, NGPBX  )                               /   &
     &     00.50,  1.21,  1.92,  2.63,  3.65,  4.36,  5.07,  5.78,      &
     &      6.80,  7.51,  8.22,  8.93,  9.95, 10.66, 11.37, 12.08,      &
     &     13.10, 13.81, 14.52, 15.23, 16.25, 19.96, 17.67, 18.38,      &
     &     19.40, 20.11, 20.82, 21.53, 22.55, 23.26, 23.97, 24.68 /
      DATA ( az(2,i), i = 1, NGPBX  )                               /   &
     &      0.00,  0.61,  1.22,  1.83,  2.44,  3.05,  3.66,  4.27,      &
     &      4.88,  5.49,  6.10,  6.71,  7.32,  7.92,  8.53,  9.14,      &
     &      9.75, 10.36, 10.97, 11.58, 12.19, 12.80, 13.41, 14.02,      &
     &     14.63, 15.24, 15.85, 16.46, 17.07, 17.68, 18.29, 18.90 /
      DATA ( bz(2,i), i = 1, NGPBX  )                               /   &
     &      0.51,  1.12,  1.73,  2.34,  2.95,  3.56,  4.17,  4.78,      &
     &      5.38,  5.99,  6.60,  7.21,  7.82,  8.43,  9.04,  9.65,      &
     &     10.26, 10.87, 11.48, 12.09, 12.70, 13.31, 13.92, 14.53,      &
     &     15.14, 15.75, 16.36, 16.97, 17.58, 18.19, 18.80, 19.41 /

      DATA ( az(3,i), i = 1, NGTdeV  )                               /   &
     &     00.00, 00.70,  1.40,  2.10,  2.95,  3.65,  4.35,  5.05,      &
     &      5.90,  6.60,  7.30,  8.00,  8.85,  9.55, 10.25, 10.95,      &
     &     11.80, 12.50, 13.20, 13.90, 14.75, 15.45, 16.15, 16.85,      &
     &     17.70, 18.40, 19.10, 19.80, 20.65, 21.35, 22.05, 22.75 /
      DATA ( bz(3,i), i = 1, NGTdeV  )                               /  &
     &     00.55,  1.25,  1.95,  2.65,  3.50,  4.20,  4.90,  5.60,      &
     &      6.45,  7.15,  7.85,  8.55,  9.40, 10.10, 10.80, 11.50,      &
     &     12.35, 13.05, 13.75, 14.45, 15.30, 16.00, 16.70, 17.40,      &
     &     18.25, 18.95, 19.65, 20.35, 21.20, 21.90, 22.60, 23.30 /
!     DATA ( az(4,i), i = 1, NGTors  )                               /
!    ^      0.00,  0.95,  1.90,  2.85,  3.80,  4.85,  5.80,  6.75,
!    ^      7.70,  8.75,  9.70, 10.65, 11.60, 12.65, 13.60, 14.55,
!    ^     15.50, 16.55, 17.50, 18.45, 19.40, 20.45, 21.40, 22.35,
!    ^     23.30, 24.35, 25.30, 26.25, 27.20, 28.25, 29.20, 30.15,
!    ^     31.10, 32.05                                           /
!     DATA ( bz(4,i), i = 1, NGTors  )                               /
!    ^      0.85,  1.80,  2.75,  3.70,  4.65,  5.70,  6.65,  7.60,
!    ^      8.55,  9.60, 10.55, 11.50, 12.45, 13.50, 14.45, 15.40,
!    ^     16.35, 17.40, 18.35, 19.30, 20.25, 21.30, 22.25, 23.20,
!    ^     24.15, 25.20, 26.15, 27.10, 28.05, 29.10, 30.05, 31.00,
!     ^     31.95, 32.90                                           /
      DATA ( az(4,i), i = 1, NGTors  )                               /   &
     &      0.00,  1.05,  2.10,  3.15,  4.20,  5.60,  6.65,  7.70,      &
     &      8.75, 10.15, 11.20, 12.25, 13.30, 14.70, 15.75, 16.80,      &
     &     17.85, 19.25, 20.30, 21.35, 22.40, 23.80, 24.85, 25.90,      &
     &     26.95, 28.35, 29.40, 30.45, 31.50, 32.90, 33.95, 35.00,      &
     &     36.05, 37.10                                           /

      DATA ( bz(4,i), i = 1, NGTors  )                               /   &
     &      0.85,  1.90,  2.95,  4.00,  5.05,  6.45,  7.50,  8.55,      &
     &      9.60, 11.00, 12.05, 13.10, 14.15, 15.55, 16.60, 17.65,      &
     &     18.70, 20.10, 21.15, 22.20, 23.25, 24.65, 25.70, 26.75,      &
     &     27.80, 29.20, 30.25, 31.30, 32.35, 33.75, 34.80, 35.85,      &
     &     36.90, 37.95                                           /

      DATA ( az(5,i), i = 1, NGSslow  )                               / &
     &      0.00,  0.30,  0.60,  0.90,  1.20,  1.50,  1.80,  2.10,      &
     &      2.40,  2.70,  3.00,  3.30,  3.60,  3.90,  4.20,  4.50 /
      DATA ( bz(5,i), i = 1, NGSslow  )                               / &
     &      0.25,  0.55,  0.85,  1.15,  1.45,  1.75,  2.05,  2.35,      &
     &      2.65,  2.95,  3.25,  3.55,  3.85,  4.15,  4.45,  4.75 /

      DATA ( az(6,i), i = 1, NGJET  )                               /   &
     &      0.00,  1.10,  2.20,  3.30,  4.40,  5.50,  6.60,  7.70,      &
     &      8.80,  9.90, 11.00, 12.10, 13.20, 14.30, 15.40, 16.50,      &
     &     17.60, 18.70, 19.80, 20.90, 22.00, 23.10, 24.20, 25.30,      &
     &     26.40, 27.50, 28.60, 29.70, 30.80, 31.90, 33.00, 34.10 /

      DATA ( bz(6,i), i = 1, NGJET  )                               /   &
     &      0.90,  2.00,  3.10,  4.20,  5.30,  6.40,  7.50,  8.60,      &
     &      9.70, 10.80, 11.90, 13.00, 14.10, 15.20, 16.30, 17.40,      &
     &     18.50, 19.60, 20.70, 21.80, 22.90, 24.00, 25.10, 26.20,      &
     &     27.30, 28.40, 29.50, 30.60, 31.70, 32.80, 33.10, 34.20 /

      DATA ( az(7,i), i = 1, NGTFTR )                               /   &
     &      0.00,  0.70,  1.40,  2.10,  2.80,  3.50,  4.20,  4.90,      &
     &      5.60,  6.30,  7.00,  7.70,  8.40,  9.10,  9.80, 10.50,      &
     &     11.20, 11.90, 12.60, 13.30, 14.00, 14.70, 15.40, 16.10,      &
     &     16.80, 17.50, 18.20, 18.90, 19.60, 20.30, 21.00, 21.70 /

      DATA ( bz(7,i), i = 1, NGTFTR )                               /   &
     &      0.55,  1.25,  1.95,  2.65,  3.35,  4.05,  4.75,  5.45,      &
     &      6.15,  6.85,  7.55,  8.25,  8.95,  9.65, 10.35, 11.05,      &
     &     11.75, 12.45, 13.15, 13.85, 14.55, 15.25, 15.95, 16.65,      &
     &     17.35, 18.05, 18.75, 19.45, 20.15, 20.85, 21.55, 22.25 /

      DATA ( az(8,i), i = 1, NGHL1M )                               /   &
     &      0.00,  1.20,  2.40,  3.60,  4.80,  6.00,  7.20,  8.40,      &
     &      9.60, 10.80, 12.00, 13.20 /
      DATA ( bz(8,i), i = 1, NGHL1M  )                              /   &
     &      0.90,  2.10,  3.30,  4.50,  5.70,  6.90,  8.10,  9.30,      &
     &     10.50, 11.70, 12.90, 14.10 /

      DATA ( az(9,i), i = 1, NGEAST )                               /   &
     &      0.00,  0.85,  2.15,  3.45,  4.75,  6.05,  7.35,  8.65,      &
     &      9.95, 11.30, 12.65, 13.95, 15.25, 16.55, 17.85, 19.15,      &
     &     20.45, 21.75, 23.10, 24.45, 25.75, 27.05, 28.35, 29.65,      &
     &     30.95, 32.25, 33.55, 34.90, 36.25, 37.55, 38.85, 40.15,      &
     &     41.45, 42.75, 44.05, 45.35, 46.70  /
      DATA ( bz(9,i), i = 1, NGEAST )                               /   &
     &      0.50,  1.85,  3.15,  4.45,  5.75,  7.05,  8.35,  9.65,      &
     &     10.95, 12.30, 13.65, 14.95, 16.25, 17.55, 18.85, 20.15,      &
     &     21.45, 22.75, 24.10, 25.45, 26.75, 28.05, 29.35, 30.65,      &
     &     31.95, 33.25, 34.55, 35.90, 37.25, 38.55, 39.85, 41.15,      &
     &     42.45, 43.75, 45.05, 46.35, 47.20  /
!******************************************************************************
      DATA ( az(10,i), i = 1, NGEAS2 )                              /   &
     &      0.00,  0.85,  1.65,  2.45,  3.25,  4.05,  4.85,  5.70,      &
     &      6.55,  7.35,  8.15,  8.95,  9.75, 10.55, 11.40, 12.25,      &
     &     13.05, 13.85, 14.65, 15.45, 16.25, 17.10, 17.95, 18.75,      &
     &     19.55, 20.35, 21.15, 21.95, 22.80, 23.65, 24.45, 25.25,      &
     &     26.05, 26.85, 27.65, 28.50, 29.35, 30.15, 30.95, 31.75,      &
     &     32.55, 33.35, 34.20, 35.05, 35.85, 36.65, 37.45, 38.25,      &
     &     39.05, 39.90, 40.75, 41.55, 42.35, 43.15, 43.95, 44.75,      &
     &     45.60  /
      DATA ( bz(10,i), i = 1, NGEAS2 )                              /   &
     &      0.60,  1.45,  2.25,  3.05,  3.85,  4.65,  5.45,  6.30,      &
     &      7.15,  7.95,  8.75,  9.55, 10.35, 11.15, 12.00, 12.85,      &
     &     13.65, 14.45, 15.25, 16.05, 16.85, 17.70, 18.55, 19.35,      &
     &     20.15, 20.95, 21.75, 22.55, 23.40, 24.25, 25.05, 25.85,      &
     &     26.65, 27.45, 28.25, 29.10, 29.95, 30.75, 31.55, 32.35,      &
     &     33.15, 33.95, 34.80, 35.65, 36.45, 37.25, 38.05, 38.85,      &
     &     39.65, 40.50, 41.35, 42.15, 42.95, 43.75, 44.55, 45.35,      &
     &     46.20  /
!******************************************************************************
!
      DATA       pmax  ,       xmin  ,        xmax,         anzac  /    &
     &     1.0d0,-10.d0, +10.d0, 1.2d0 /
      DATA RelPha, MinNpar / 120.d0, 1.50d0 /
!
      INTEGER Curious
      DATA    Curious/0/
      REAL(fp)    RE81, RE82, RE8360
      RE8360 = 360.d0
!-------------------
! dmc -- changes to avoid use of un-initialized local variables

      phi = 0.0d0

!-------------------
!
!     Copy argument list into the local variables:
!
 1    continue
      f0 = fghz * 1000.
      nplot  = nslices
      NumRay = nrays
      MinNpar= abs(nparaMin)
      RE81   = phaseDeg
      RelPha =  mod(RE81,    RE8360)
!     RelPha = amod(phaseDeg,360.)
      do 2 i = 1, DIM
        p0(i)= ONE
        RE81    = REAL(i,kind=fp)*RelPha
        phi0(i) =  mod( RE81          , RE8360)
!       phi0(i) = amod(float(i)*RelPha, 360.)
 2    continue

      if      ( coupler .EQ. 'PBXMFAST' )  then
        G = 1
        nguide = NGPBX
      else if ( coupler .EQ. 'PBXMSLOW' )  then
        G = 2
        nguide = NGPBX
      else if ( coupler .EQ. 'TOKDEVAR' )  then
        G = 3
        nguide = NGTdeV
      else if ( coupler .EQ. 'TORSUPRA' )  then
        G = 4
        nguide = NGTors
        p0(1) =     - ONE
        p0(NGTors)= - ONE

        phi0(1) =    ZERO
        phi0(NGTors)=ZERO
        phi0(2) =    ZERO
        phi0(3) =   R90
        phi0(4) =  R180
        phi0(5) =  R270
        do i = 6, NGTors-1
          RE81    = phi0(i-4)+RelPha
          phi0(i) =  mod( RE81,              RE8360)
!         phi0(i) = amod((phi0(i-4)+RelPha), 360.)
        enddo
      else if ( coupler .EQ. 'SLOWSLOW' ) then
        G = 5
        nguide = NGSslow
      else if ( coupler .EQ. 'JET_LHCD' ) then
        G = 6
        nguide = NGJET
        phi0(1) =  ZERO
        phi0(2) =   R90
        phi0(3) =  R180
        phi0(4) =  R270
        do i = 5, NGTors
          RE81    = phi0(i-4)+RelPha
          phi0(i) =  mod( RE81,              RE8360)
!         phi0(i) = amod((phi0(i-4)+RelPha), 360.)
        enddo
      else if ( coupler .EQ. 'TFTRLHCD' ) then
        G = 7
        nguide = NGTFTR
      else if ( coupler .EQ. 'HL1MLHCD') then
        G = 8
        nguide = NGHL1M
      else if ( coupler .EQ. 'EASTLHCD') then
        G = 9
        nguide = NGEAST
    p0(1) =	- ONE
		p0(10)=	- ONE
		p0(19)= - ONE
		p0(28)= - ONE
		p0(NGEAST)= - ONE
		phi0(1)=	ZERO
		phi0(2)=	ZERO
		phi0(3)=	R90
		phi0(4)=	R180
		phi0(5)=	R270
		phi0(6)=	ZERO
		phi0(7)=	R90
		phi0(8)=	R180
		phi0(9)=	R270
		phi0(10)=	ZERO
		j2 = 10
		do i = 1, 3
		  do j = 1, 8
		    j2 = j2+1
		    RE81 = phi0(j2-9)+RelPha
			  phi0(j2) = mod(RE81, RE8360)
		  enddo
		  j2 = j2+1
		  phi0(j2) = ZERO
		enddo

      else if ( coupler .EQ. 'EASTLHC2') then
        G = 10
        nguide = NGEAS2
                p0(1) =	- ONE
		p0(8)=	- ONE
		p0(15)= - ONE
		p0(22)= - ONE
		p0(29)= - ONE
		p0(36)= - ONE
		p0(43)= - ONE
		p0(50)= - ONE
		p0(NGEAS2)= - ONE
		phi0(1)=	ZERO
		phi0(2)=	ZERO
		phi0(3)=	R90
		phi0(4)=	R180
		phi0(5)=	R270
		phi0(6)=	ZERO
		phi0(7)=	R90
		phi0(8) = ZERO
		j2 = 8
		do i = 1, 7
		  do j = 1, 6
		    j2 = j2+1
		    RE81 = phi0(j2-7)+RelPha
			  phi0(j2) = mod(RE81, RE8360)
		  enddo
		  j2 = j2+1
		  phi0(j2) = ZERO
		enddo

      else
        G = 1
        nguide = NGPBX
      endif

      if (nplot .GT. NPDIM) nplot = NPDIM

      isw=0
      if(nd.GT.5001) nd=5001
      if(nguide.LT.0) isw=2
      if(nguide.LT.0) nguide=-nguide
      if(nguide.GT.idim) go to 999
!         iterate phase
      nphase=1
      if(delphi .GT.14.d0)nphase = (R180/delphi)+1.1d0
!     if(delphi .gt. 14.) nphase = (180./delphi) + 1.1
!         include waveguide impedance?(iepsw=1)
      iepsw=1
      if(epsr.LT.ZERO) iepsw=0
      epsr=abs(epsr)
      if(epsr.LT.ZERO) epsr=ONE
      dencr=f0**2/0.009d0**2
!         iterate on edge density
      nedgd=1
      if(edged.LT.ZERO) nedgd=7
      edged0=abs(edged)
      if(nedgd .GT. ONE .and. (edged0/dencr .GT. 111.d0          &
     &                  .or.   edged0/dencr .LT. ONE))                  &
     & edged0=dencr
!         iterate density gradient
      ngradn=1
      if(gradn.LT.ZERO) ngradn=4
      gradn0=abs(gradn)
	iskip=0
	if(nphase .GT. 1 .or. ngradn .GT. 1 .or. nedgd .GT. 1)iskip=1
!
!         read in data for each waveguide
      pave=ZERO
      ifast=1
      zmin= 9999.d0
      zmax=-9999.d0
!      do 20 i=1,nguide
!      read(21,2020,END=18,err=999) p0(i),phi0(i),az(i),bz(i)
! 2020 format(6g10.0)
!	go to 19
!cc	assign values based on wg no. 1
! 18   p0(i)=p0(i-1)
!      phi0(i)=phi0(i-1)+phi0(1)
!      az(i)=bz(i-1)+az(1)
!      bz(i)=bz(i-1)+bz(1)
! 19   continue
!cc         check for passive guides
!      if(p0(i).lt.0.0) isw=1
!      if(p0(i).gt.0.0)pave=pave+p0(i)
!      if(az(i).lt.zmin) zmin=az(i)
!      if(bz(i).gt.zmax) zmax=bz(i)
!cc         check for symmetry of the waveguides(ifast=1)
!      if(i.gt.1.and.abs((az(i)-az(i-1))-(bz(i)-bz(i-1))).gt.0.01)ifast=0
! 20   continue
!      pave=pave/nguide
!      twidth=abs(zmax-zmin)
      pave = ONE
      twidth = abs ( bz(G,nguide) - az(G,1) )
!
!         calculate constants
      nloop=1
      if(isw.GT.0) nloop=5
      c0=(0.0d0,0.0d0)
      c1=(0.0d0,1.0d0)
      w=2.*PI*f0*1.d6
      wdc=w/2.997d10
      al0=2.997d4/f0 ! speed-of-light / freq = wavelength (cm) in vacuum
      height=a(G)     ! height of waveguide
      beta=sqrt(epsr-(al0/2./height)**2)
      beta1=beta
      if(iepsw.EQ.0) beta1=sqrt(epsr)
      xp=c1*wdc*xprobe*beta
!
!         begin main loop
!
!         iterate on phase
      do 800 nph=1,nphase
      if(delphi .LT. 14.d0) go to 26
      do 25 i=1,nguide
   25 phi0(i)=(nph-ONE)*(i-ONE)*delphi
   26 continue
!
!         iterate density gradient
      do 800 ngn=1,ngradn
      gradn=gradn0*10.**(ngn-1)
!
!         iterate edge density
      do 800 ned=1,nedgd
      edged=edged0*3.1622777**(ned-1.)
      ratio=edged/dencr-ONE
      if(ratio.LT.-0.999d0) edged=dencr
      if(ratio.LT.-0.999d0) ratio=ZERO
!          iterate for passive guides or phase feedback of density gradient
      do 788 nl=1,nloop
!         incident electric field in i'th guide in v/cm
      istop=1
      if(nl.EQ.1) istop=0
      do 30 i=1,nguide
      b(i)=abs(az(G,i)-bz(G,i))
      if(nl.ne.1.and.p0(i).LT.ZERO.and.p0(i).LE.p0r(i))                 &
     &                                           alp=0.9d0
      if(nl.ne.1.and.p0(i).LT.ZERO.and.p0(i).GT.p0r(i))                 &
     &                                           alp=1.5d0
      if(nl.ne.1.and.p0(i).LT.ZERO)                                     &
     &                            p0(i)=alp*p0r(i)+(ONE-alp)*p0(i)
      if(nl.EQ.1.and.p0(i).LT.ZERO) p0(i)=-pave/10.
      if(nl.ne.1.and.p0(i).LT.ZERO.and.abs(rm(i)-ONE) .GT.              &
     &   0.01d0)                           istop=0
      arg=4.d3*abs(p0(i))*120.*PI/height/b(i)/beta
      phold(i)=phi(i)
      phi(i)=phi0(i)
      if(nl.EQ.1.or.i.EQ.1.or.isw.ne.2) go to 29
!     if(nl.eq.1.or.i.eq.1.or.isw.ne.2) go to 30
      ph0=phi0(i)-phi0(i-1)
      if(phi0(i).LT.phi0(i-1)) ph0=ph0+RE8360
      ph1=btph(i)-btph(i-1)
      if(btph(i).LT.btph(i-1)) ph1=ph1+RE8360
      ph2=phold(i)-phold(i-1)
      if(phold(i).LT.phold(i-1)) ph2=ph2+RE8360
      del=ph0-ph1
      idel=del/R180
      del=del-idel*R180
      if(abs(del).GT.ONE)istop=0
      phi(i)=phi(i-1)+ph2+del*0.5
      nmod=phi(i)/RE8360
      phi(i)=phi(i)-nmod*RE8360
   29 RE81 = PI*phi(i)/R180
   30 e0(i)=sqrt(arg)*exp(RE81*c1)
!  30 e0(i)=sqrt(arg)*exp(PI*c1*phi(i)/R180)
      if(istop.EQ.1) go to 799
      alph3=(gradn/(f0/.009d0)**2/wdc)**0.333333
      RE81 = PI/(3.d0)
      c3=exp(RE81*c1)*a2/a1*alph3
!     c3=cexp(PI*c1/3.)*a2/a1*alph3
      RE81 = PI/2.d0
      c2=exp(RE81*c1)*a2/a1*alph3
!     c2=cexp(PI*c1/2.)*a2/a1*alph3
      c9=2.*c1*x0*wdc
!
!         zero matricies
      do 51 l=1,nguide
      r(l)=c0
!
        do 50 k=1,nguide
        ai(l,k)=c0
        aa(l,k)=c0
 50     continue
 51   continue
!
!         perform nz integration from -anmax to anmax by simpson's rule
      dnz=anmax/nd
      anz=0.0d0
!         if guides are symmetrical,
!         then don't need to compute all the matrix elements
      ngset=nguide
      if(ifast.EQ.1)ngset=1
      do 100 n=1,nd,2
      do 100 n2=1,2
      anzi=anz*c1*wdc
!         abs(nz) > 1
      if(abs(anz).LE.ONE) go to 60
      anx=c1*sqrt(anz*anz-1.)
      anx3=(anz**2-ONE)**(0.1666667)
      if(anx3.LT.0.01d0) anx3=0.01d0
      c5=c3/anx3
      if(ratio.LE.0.1d0) go to 70
!         density step present
      w0=-anx3**2*ratio/alph3**2
      CALL lsc_airy(air,bir,aip,bip,w0)
      c5=alph3*(aip+c1*bip)/(air+c1*bir)/anx3
      go to 70
!         abs(nz) < 1
   60 anx=sqrt(1.-anz**2)
      anx3=(ONE-anz**2)**(0.166667)
      if(anx3.LT.0.01d0) anx3=0.01d0
      c5=c2/anx3
      if(ratio.LE.0.1d0) go to 70
!         density step present
      w0=anx3**2*ratio/alph3**2
      CALL lsc_airy(air,bir,aip,bip,w0)
      if(abs(air).LT.0.001d0) air=0.001d0
      c5=-c1*alph3*aip/air/anx3
   70 denom=anx*anz*anz
      if( abs(denom).LT.0.001d0) denom=0.001d0
      rho=exp(c9*anx)*(ONE-c5)/(ONE+c5)
      c6=n2*(ONE-rho)/(ONE+rho)/denom
      do 80 l=1,ngset
      RE81 = bz(G,l)
      RE82 = az(G,l)
      c7=exp(anzi*RE81)   -exp(anzi*RE82)
!     c7=cexp(anzi*bz(G,l))-cexp(anzi*az(G,l))
      do 80 k=1,nguide
      RE81 = bz(G,k)
      RE82 = az(G,k)
      c8=c7*(exp(-anzi*RE81)    -exp(-anzi*RE82  ))
!     c8=c7*(cexp(-anzi*bz(G,k))-cexp(-anzi*az(G,k)))
   80 ai(l,k)=ai(l,k)+c6*(c8+conjg(c8))
  100 anz=anz+dnz
!
!         guides are symmetrical - use symmetry to compute matrix elements ai
      if(ifast.EQ.0) go to 200
      do 150 l=2,nguide
      do 130 k=l,nguide
  130 ai(l,k)=ai(l-1,k-1)
      do 150 k=1,l-1
  150 ai(l,k)=ai(k,l)
!         multiply 'ai' by constant, set up matrix equation
  200 c10=dnz/3./PI/wdc/beta1
      do 300 l=1,nguide
      bb(l)=e0(l)
      aa(l,l)=e0(l)
      do 300 k=1,nguide
      ai(l,k)=ai(l,k)/b(l)*c10
      bb(l)=bb(l)-e0(k)*ai(l,k)
  300 aa(l,k)=aa(l,k)+e0(k)*ai(l,k)
!       print 8000,((l,k,aa(l,k),bb(l),ai(l,k),k=1,nguide),l=1,nguide)
!  8000 format(2i5,6e12.4)
!
!         invert matrix 'aa'  -  put in array 'ai'
      CALL lsc_invert(idim,nguide,aa,ai)
!
!         multiply aa(inverse) * bb  =  r
      do 400 l=1,nguide
      do 400 k=1,nguide
!       print 8000,l,k,ai(l,k),aa(l,k)
  400 r(l)=r(l)+ai(l,k)*bb(k)
!
!         compute parameters in each guide
      pplas=ZERO
      pint=ZERO
      area=ZERO
      do 500 k=1,nguide
      if(p0(k).GT.ZERO)area=area+height*b(k)
      rm(k)= abs(r(k))
      rr=REAL(r(k),kind=fp)
      ri=AIMAG(r(k))
      rp(k)=atan2(ri,rr)*R180/PI
      if(p0(k).GT.ZERO)pint=pint+p0(k)
      p0r(k)=p0(k)*rm(k)**2
      if(p0(k).GT.ZERO)pplas=pplas+p0(k)-p0r(k)
      et(k)=e0(k)*(exp(xp)+r(k)*exp(-xp))
      ei(k)= abs(e0(k))
      etr=REAL(et(k),kind=fp)
      eti=AIMAG(et(k))
      etm(k)= abs(et(k))
      etph(k)=atan2(eti,etr)*R180/PI
      bt(k)=-e0(k)*(exp(xp)-r(k)*exp(-xp))
      btm(k)= abs(bt(k))
      btr=REAL(bt(k),kind=fp)
      bti=AIMAG(bt(k))
  500 btph(k)=atan2(bti,btr)*R180/PI
      teff=100.*pplas/PInt
      pdave=pint/area
!
!         print general information for run
      epsrp=epsr
      if(iepsw.EQ.0)epsrp=-epsr
!	loop=1+(nguide+15)/16
!	do 777 lll=1,loop

!	ix=1
!	iy=765
!                                       Take out the write statements, normally
                                        if (Curious .EQ. 1) then
!	write(nTSCscrn,2100) pint,pplas
!	write(nTSCscrn,2101) f0,nguide
!	write(nTSCscrn,2102) gradn,x0
!	write(nTSCscrn,2103) edged,twidth
!	write(nTSCscrn,2104) pdave,teff
!	write(nTSCscrn,2105) height,epsrp
!	write(nTSCscrn,2106) dnz,anmax
!	write(nTSCscrn,2107) xprobe,nl
 2100 format(' incident rf power(kW) =',f9.0,3x,                        &
     & ' rf power to plasma(kw)=' ,f9.0)
 2101 format(' frequency(MHz)        =' ,f9.0,3x,                       &
     & ' no. of waveguides     =' ,i9)
 2102 format(' density gradient(cm-4)=' ,1pe9.2,3x,                     &
     & ' plasma-guide gap(cm)  =' ,0pf9.3)
 2103 format(' edge density(cm-3)    =' ,1pe9.2,3x,                     &
     & ' total grill width(cm) =' ,0pf9.3)
 2104 format(' ave.pow.den.(kW/cm**2)=' ,f9.3,3x,                       &
     & ' % efficiency          =' ,f9.0)
 2105 format(' waveguide height(cm)  =' ,f9.3,3x,                       &
     & ' dielectric constant   =' ,f9.3)
 2106 format(' nz integration step   =' ,f9.3,3x,                       &
     & ' maximum nz            =' ,f9.3)
 2107 format(' rf probe position(cm) =' ,f9.3,3x,                       &
     & ' iteration             =' ,i9)
!
                                        endif
!
!         print results for each guide
!	ngpr1= 1+16*(lll-1)
!	ngpr2=16+16*(lll-1)
        ngpr1 = 1
        ngpr2 = nguide
	if(ngpr2.GT.nguide) ngpr2=nguide
	if(ngpr1.GT.nguide) go to 543
!                                       Take out more wri
!
                                        if (Curious .EQ. 1) then
!	write (nTSCscrn,2109)
!	write (nTSCscrn,2110)
 2109	format(1x)
 2110    format(' guide  p0 p0ref mag(r) ph(r)',                        &
     &          '   e0  phe0  etot  phet  btot  phbt b(cm)')
        do 520 i=ngpr1,ngpr2
!        write(nTSCscrn,2120)                                            &
!     &                i,p0(i),p0r(i),rm(i),rp(i),ei(i),phi(i),          &
!     &                etm(i),etph(i),btm(i),btph(i),b(i)
 2120     format(i4,1x,f6.0,f5.0,f7.3,2f6.0,f5.0,4f6.0,f6.2)
 520	continue
                                        endif
!
 543	continue
!	if(lll.ne.loop) go to 776
!
!         plot power vs. nz
      smx=ZERO
      sum=ZERO
      t1=height*beta/beta1*1.e-4/24./wdc/PI**2
      c2=c1*t1*alph3
      c9=c1*x0*wdc
	danplt=anplt*2./nplot
	npmin=1
	npmax=1
!               integrate 1/nz and 1/nz**2 over spectrum
	asum1=ZERO
	asum2=ZERO
	asum3=ZERO
	fracc=ZERO

      do 600 i=1,nplot
      anz=-anplt+(i-1)*anplt/(nplot*0.5-1.)
      sn(i)=anz
        if(anz.LT.xmin+danplt)npmin=i
        if(anz.LT.xmax)npmax=i
      c8=c0
      if(abs(anz).LT.1.1) go to 580
      anx=sqrt(anz**2-1.)*c1
      t2=ONE/anz**2
      anzi=anz*c1*wdc
      c6=c0
      c7=c0

      do 550 k=1,nguide
      RE81 = az(G,k)
      RE82 = bz(G,k)
      c5=exp(-anzi*RE81)   -exp(-anzi*RE82)
!     c5=cexp(-anzi*az(G,k))-cexp(-anzi*bz(G,k))
      c7=c7+anx*beta1*e0(k)*(1.-r(k))*c5
!     c7=c7+anx*beta1*e0(k)*(1.-r(k))*c5
  550 c6=c6+e0(k)*(1.+r(k))*c5

      enz=(c6+c7)/4./c1
      rho=(c6-c7)/(c6+c7)
      c8=enz*(exp(c9*anx)+rho*exp(-c9*anx))
      if(ratio.LE.0.1) go to 560
      w0=-(anz**2-1.)**0.333333*ratio/alph3**2
      CALL lsc_airy(air,bir,aip,bip,w0)
      c5=c2*(aip-c1*bip)/(air-c1*bir)
      go to 570
  560 RE81 = PI/3.
      c5=c2*a2/a1*exp(-RE81*c1)
! 560 c5=c2*a2/a1*cexp(-PI*c1/3.)
  570 c8=c8*conjg(c8)*t2*c5/(anz**2-1.)**0.666667
  580 sp(i)=REAL(c8,kind=fp)
      if(sp(i).GT.smx) smx=sp(i)
      sum=sum+sp(i)
      s(i)=sum
        sign=anz/abs(anz)
        delt=sign*sp(i)/anz**2
        asum1=asum1+delt
      if(abs(anz).LT.anzac) go to 600
        asum2=asum2+delt
        asum3=asum3+delt
        fracc=fracc+sp(i)
600   continue
!
!	ave. nz over all power nz>1
                asum1=abs(asum1/sum)**0.5
!	ave. nz over accessible power |nz|>nza
                asum2=abs(asum2/fracc)**0.5
!	ave. nz over accessible power, normalize to total power
                asum3=abs(asum3/sum)**0.5
                fracc=fracc/sum
!
!         plot nz spectrum in the plasma
      RE81= smx
      lsmx= log10(RE81)
!     lsmx=alog10(smx)
      t=10.**lsmx
      ll=smx/t+1
      smx=ll*t
	if(pmax.GT.ONE) smx=pmax
!
!
!                CALL lsc_initia
!
!       CALL lsc_ezsets(120,920,50,480,xmin,xmax,0.,smx,1)
	ndif=npmax-npmin+1
!       CALL lsc_ezcurv(sn(npmin),sp(npmin),ndif)

      do 700 i=npmin,npmax
        s(i)=s(i)*smx/sum
 700  continue
!       CALL lsc_ezpnts(sn(npmin),s(npmin),ndif)
	ix1=xmax-xmin+.01
        ix1=ix1/2
	ix2=2
	if(ix1.LT.5)ix2=10
!       CALL lsc_ezaxes(ix1,ix2,5,2)
!       CALL lsc_ezwrit(520,8,'nz$',1,0)
!       CALL lsc_ezwrit(5,225,'power(kw/nz)$',1,1)
!
		ix=650
		iy=485
		write(pline,2221) asum1
!                CALL lsc_EZwrit(ix, iy, pline, 0, 0)
                  iy = iy - 25
 2221	format('<|nz|**-2>t=',f6.3)
		write(pline,2222) asum2
!                CALL lsc_EZwrit(ix, iy, pline, 0, 0)
                  iy = iy - 25
 2222	format('<|nz|**-2>a=',f6.3)
		write(pline,2223) asum3
!                CALL lsc_EZwrit(ix, iy, pline, 0, 0)
                  iy = iy - 25
 2223	format('<|nz|**-2>at=',f6.3)
		write(pline,2224) anzac
!                CALL lsc_EZwrit(ix, iy, pline, 0, 0)
                  iy = iy - 25
 2224	format('  nz(acc)=',f7.3)
		write(pline,2225) fracc
!                CALL lsc_EZwrit(ix, iy, pline, 0, 0)
                  iy = iy - 25
 2225	format('fract(acc)=',f6.3)
!
 776	if(nloop.EQ.1.and.iskip.EQ.0) then
!         CALL lsc_tinput(idum)
        endif
! 777    continue
  788 continue

!         make hard copy
 799    continue
!       CALL lsc_finish
!
!     Plot the spectrum as a serites of bars
!     CALL lsc_ezsets(120,920,100,600,xmin,xmax,0.,smx,1)
!     CALL lsc_ezbars(sn(npmin),sp(npmin),ndif, 'y')
!     CALL lsc_ezaxes(ix1,ix2,5,2)
!     CALL lsc_ezwrit(520,25,'nz$',1,0)
!     CALL lsc_ezwrit(25,350,'rel power$',1,1)
!     write(pline,'(''phase, nrays, nslic: '',f8.4,1x,i3,1x,i3,''$'')')
!    ^             phaseDeg, nrays, nslic
!     -----------------------------------------------------------------|
!     CALL lsc_ezwrit(120,700, pline, 0, 0)

!     CALL lsc_finish

  800 continue
!
!     Find the NumRay largest peaks
!     But first eliminate the low nparallels
      do 820 i= npmin, npmax
        if (abs(sn(i)) .LE. MinNpar) sp(i) = ZERO
 820  continue
!     NumRay from argument list; NumBeg labels the first of NumRays in stack
      NumBeg = npmax - (NumRay-1)
!
!     Put the highest sp's at top of the array
      CALL lsc_PikSr2NR ( ndif  , sp(npmin) , sn(npmin)  )
!
!     Take the NumRay highest sp's and put them in order
!     accorting to the sn's
      CALL lsc_PikSr2NR ( NumRay, sn(NumBeg), sp(NumBeg) )
!
!     Normalize so the NumRay highest sp's add to 1.00,
!     and find out what is the largest value
      TotPwr = ZERO
      MaxPwr = ZERO
!cxxxx do 830 i = npmin + ndif - NumRay, npmin + ndif !Terpstra found this wrong
      do 830 i = NumBeg, npmax
        TotPwr = TotPwr + sp(i)
 830  continue
        ii = 1
!cxxxx do 840 i = npmin + ndif - NumRay, npmin + ndif ! This one too!
      do 840 i = NumBeg, npmax
        sp(i) = sp(i)/TotPwr
        nparalls(ii) = sn(i)
        amplitds(ii) = sp(i)
        if(sp(i) .GT. MaxPwr) MaxPwr = sp(i)
        ii = ii+1
 840  continue

      if ( abs(RelPha - R180) .GT. 10.d0   .and.                 &
     &     abs(RelPha - ZERO) .GT. 10.d0   .and.                 &
     &     turnnegs             .EQ. 1           ) then
        do 900 i = 1, NumRay
          nparalls(i) = abs(nparalls(i))
 900    continue
      endif
!
!
      return
  999 CALL lsc_LSCstop(' error in the Brambilla code')
      return
      END
!
!----------------------------------------------------------/invert
!         invert 'aa' - apply same operations to unit matrix as to 'aa'
      SUBROUTINE lsc_invert(idim,nguide,aa,ai)
use iso_c_binding, only : fp => c_double
implicit none

      INTEGER idim, nguide, l,k,l1
!     COMPLEX f,f0,aa(idim,1),ai(idim,1) ! Sun f90 wants proper dimensioning
      COMPLEX f,f0,aa(idim,idim),ai(idim,idim)
!         set 'ai' equal to unit matrix
      do 100 l=1,nguide
      do 100 k=1,nguide
      ai(l,k)=(0.0d0,0.0d0)
  100 ai(l,l)=(1.0d0,0.0d0)
!
      do 500 l=1,nguide
!
      if(l.EQ.nguide) go to 180
      do 150 l1=l+1,nguide
      f=aa(l,l)+aa(l1,l)
      f0=(1.0d0,0.0d0)
      if( abs(f).LT. abs(aa(l,l))) f0=(-1.0d0,0.0d0)
      do 120 k=1,nguide
      aa(l,k)=aa(l,k)+aa(l1,k)*f0
  120 ai(l,k)=ai(l,k)+ai(l1,k)*f0
  150 continue
!         set aa(l,l) = 1
  180 f=aa(l,l)
      do 200 k=1,nguide
      ai(l,k)=ai(l,k)/f
  200 aa(l,k)=aa(l,k)/f
!         set aa(l1,l) = 0  ,  l1 /= l
      do 400 l1=1,nguide
      if(l1.EQ.l) go to 400
      f=aa(l1,l)
      do 300 k=1,nguide
      ai(l1,k)=ai(l1,k)-f*ai(l,k)
      aa(l1,k)=aa(l1,k)-f*aa(l,k)
  300 continue
  400 continue
  500 continue
      return
      END
!cc
!         compute airy functions
      SUBROUTINE lsc_airy(ai,bi,aip,bip,x)
!         x < 0
use iso_c_binding, only : fp => c_double
implicit none
      INTEGER nm, m, n
      REAL(fp)                                                              &
     &        ai, bi, aip, bip, x, x3, xx, f, ft, g, gt
      real(fp)    a1, a2
      REAL(fp)                                                              &
     &        c, cn, sn, c0, c1, c2, f1, f2
      REAL(fp)    aa(2),bb(2)
      DATA    a1/0.355028/,a2/0.2588194/
      if(x.GT.7.d0) x=7.d0
      if(abs(x).GT.7.d0) go to 100
!         small argument expansion
      nm=3+2.5*abs(x)
      do 20 m=1,2
      xx=x+0.01d0*(m-1.5)
      x3=xx**3
      f=1.d0
      ft=f
      g=xx
      gt=g
      do 10 n=1,nm
      f=f*x3/(3.*n*(3.*n-1.))
      g=g*x3/(3.*n*(3.*n+1.))
      ft=ft+f
   10 gt=gt+g
      aa(m)=a1*ft-a2*gt
   20 bb(m)=sqrt(3.)*(a1*ft+a2*gt)
      go to 500
!         asymptotic expansion
  100 do 120 m=1,2
      xx=-x-0.01d0*(m-1.5)
      c=2./3.*xx**1.5
      cn=cos(c)
      sn=sin(c)
      c0=0.3989423*xx**(-0.25)
      c1=0.0694444/c
      c2=0.0371335/c**2
      f1=1.d0-c1-c2
      f2=1.d0+c1-c2
      aa(m)=c0*(cn*f1+sn*f2)
  120 bb(m)=c0*(cn*f2-sn*f1)
  500 ai=0.5*(aa(1)+aa(2))
      bi=0.5*(bb(1)+bb(2))
      aip=100.d0*(aa(2)-aa(1))
      bip=100.d0*(bb(2)-bb(1))
      return
      END
!

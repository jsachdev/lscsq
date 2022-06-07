!subroutine lscsq_main(LhPwrMWx, iRayTrsx, iErrorx, powtsc, curtsc, dJdE,dlJdlE)
!     LSC  main control module              ---------------------------|
!                                                                      |
!     LSC  -- Lower hybrid Simulation Code
!     AcDc -- A Current Drive Code
!     Copyright D. W. Ignat, E. J. Valeo,
!               N. J. Fisch, C. F. F. Karney,
!               S. C. Jardin, S. Bernabei, D. P. Enright.
!     Princeton University
!     1992, 1993, 1994, 1995, 1996, 2000, 2021
!
!     References:
!       D. W. Ignat, Phys. Fluids, <24>, 1110 (81);
!       E. J. Valeo and D. C. Eder, J. Comp. Physics, <69>, 341 (87).
!       C. F. F. Karney and N. J. Fisch, Phys. Fluids <29>, 180 (86).
!       D. W. Ignat, E. J. Valeo, and S. C. Jardin,
!       "Dynamic Modeling of Lower Hybrid Current Drive,"
!       Nuclear Fusion <34>, 837-852 (94)
!       D. W. Ignat,  B. P. LeBlanc, C. K. Phillips,
!       J. R. Wilson, and R. V. Budny,
!       "Computational Model for Fast Wave Current Drive,"
!       Proc. 11th Conf. on Radio Frequency Power in Plasmas,
!       (Palm Springs, CA, May 1995)
!       D. W. Ignat, R. Kaita, S. C. Jardin, and M. Okabayashi
!       "Spreading of Lower Hybrid Wave Driven Currents in PBX-M,"
!       Nucl. Fusion <36> 1733-1742 (96)
!       Myunghee Ju, Jinyong Kim, KSTAR Team,
!       "A predictive study of non-inductive current driven
!       KSTAR tokamak discharge modes using a new transport/heating
!       simulation package,"
!       Nucl. Fusion <40> 1859-1866 (2000)
!
!            Warning: Local variables are assumed to be SAVEd by the
!            compiler.  While this is a standard implementation, it is
!            not ANSI standard.
!
!
!    2021 - September  -- F. M. Poli
!         rewritten the code in modern Fortran
!         fixed bugs in generation of velocity grid, interface with temeprature
!         profiles, psi grid was mixed between TSC/LSC representation
!         all edge values of profiles were wrongly set
!         the extension to multiple antennae was wrong
!         indexing rays for parallelization
!         removed plotting options
!         created standalone version that uses standard eqdsk file
!         substantially simplified code
!         the code now uses distribution function from previous time step
!         the ray initialization is now done only in one place
!         the lscsq_ramppwr subroutine has been modified so that the power is not
!         ramped from zero and built from a Maxwellian distribution at each time
!         step (it should not) 
!         extended distribution of rays to a poloidal extent (antenna)
!
!     LSC_23Mar11 -- D. McCune @PPPL
!        I am implementing a request from Siye Ding (EAST tokamak, China)
!        to allow in LSC multiple antenna/couplers each with a different
!        frequency assigned.  The prior version has only allowed a single
!        frequency for all couplers or "groups"...
!           Regression testing protocol: original NTCC circular test case
!           + 2 cases from TRANSP (one from CMOD tokamak, one from EAST)
!           were tested and found identical in the original single frequency
!           and the modified multi frequency versions of LSC
!
!        1. RayWrk.inc modifications:
!           a) "fghz" namelist scalar -> fghz(NGRPDIM); code will use
!              namelist values fghz(1:nGrps).  Default is as before:
!              fghz(1)=4.6; fghz(2:NGRPDIM)=0.0.  Logic is added to
!              the code so that for igrp=1:Ngrps, if fghz(igrp)=0.0
!              then assign igrp the last non-zero frequency in the
!              array.  So, old namelists having just a single frequency
!              value will cause simulations to run as before, with the
!              same frequency assigned to all couplers or groups.
!           b) added "fghz_now" scalar meaning "frequency of ray now
!              being followed".
!           c) integer array indx_freq(NTORDIM) indicates a frequency
!              index associated with each toroidal wave number ntor(1:ntors);
!              this is used to recover the frequency for each launched ray.
!        2. Changes in AcDc.F:
!           a) logic to handle fghz(...) namelist input with backwards
!              compatibility;
!           b) removed frequency based initializations; see new subroutine
!              "setfreq" in Rayini.
!           c) graph labeling shows all frequencies instead of a single
!              frequency.
!           d) set up indx_freq(...) with toroidal spectrum...
!              1) if DoBram=0, count the number N of actual distinct
!                 frequencies; round up ntors to be an exact multiple of N;
!                 create ntor(...) spectrum with N duplicate values of each
!                 toroidal wave number; assign spectrum magnitude to
!                 (ntor,frequency) pairs instead of just (ntor) values.
!              2) if DoBram=1, pass to brambJES the separate frequency for
!                 each coupler; maintain association of frequency index;
!                 use PikSr3NRI instead of PikSr2NR for sorting (see
!                 grids.F changes).
!        3. Changes in Rayini.F:
!           use indx_freq(...) to set the frequency for each ray; call
!           setfreq to set fghz_now and derived quantities.
!        4. Changes in ql.F:
!           a) take 1/frequency out of DqlNorm; apply this instead
!              individually to each term in the Dql summation.
!           b) use indx_freq(...) to extract appropriate frequency for
!              each ray in the loop over rays for the Dql summation.
!        5. Changes in grids.F:
!           c) create new routine PikSr3NRI based on PikSr2NR; routine
!              does sort while carrying along a 3rd array: data for
!              indx_freq(...).
!
!
!     LSC_03Dec00
!        1. Portability enhancements.  Use the "const" macro (in upper case)
!           to do things like DATA /ONE/1.0D0/ | DATA /ONE/1.0E0/ for
!           double | single precision.
!        2. Remove intrinsics in favor of generics,
!           dsqrt --> sqrt;  alog --> log ; dexp | cexp --> exp ; amin1 --> min
!        3. Change Inpval.inc comments to show the more portable use of NAMELIST
!           inputs.  OLD plflag(1) = 1, 0, 0, .. ; NEW plflag = 1, 0, 0, ..
!        4. Change dimensioning of aa, ai in brambJES to be portable to Sun f90
!           OLD  COMPLEX f,f0,aa(idim,1),ai(idim,1)
!           NEW  COMPLEX f,f0,aa(idim,idim),ai(idim,idim)
!
!     LSC_16Oct00
!        1. Make changes found necessary for running with TSC
!           set nLSCcom2 = nTSCscrn
!           recycle ns30c1 to be nTSCgraf, the unit for gnuplot stuff
!           abandon writes of GrfTitle to ns30c1, or whatever
!        2. Use "call date_and_time(zzdate,zztime,zzzone,ivals)" in gnup.F
!           instead of "CALL lsc_date", "call itime" (or "call clock" for Cray)
!           "date_and_time" is a standard intrinsic in Fortran 90.   But,
!           under Fujitsu Linux an -X9 compiler switch is needed.
!           g77 knows about "date_and_time" .
!        3. Remove spurious "fdate" in XrayCam2.F; add "date_and_time".
!           Remove spurious mention of ezlilo, not used any more.
!     LSC_03Aug00
!        1. Correct errors in Jray _graphs_ which eliminated negative J-s
!        2. Remove Rayout, mktrajb, mktrajc.  These are ray plots in
!           the vpar-psi plane.  (I did not grasp the difference, surely
!           not obvious on looking at the output.)  These did not use EZ
!           calls, so were getting hard to support.
!        3. Adjust placement of n_\parallel, n_\perp vs rt psi;
!           3d isometric of ray, poloidal cross section of ray;
!           Adjust placement of n_\parallel enhancement and accessibility
!        4. Rename iFastPiece to be iOrigMaxwl in Giruzzi of pbxio.F
!           to make code more understandable.
!        5. Recast and simplify all Plot Print flags; see PlPr.inc
!        6. Add GNplot for rays, Jrf, Prf.  Revise SG-only wave
!           scattering plots to go in the middle of ray plots
!        7. Insert three "dummy" couplers USR1 USR2 USR3 and set them
!           equal to TFTR.  The intent is for these to be modified in the
!           source.
!        8. Fiddle to get gnuplot filename to take date/time from
!           both Fujitsu f90 and G77.
!        9. Many small changes aimed at getting things to work in
!           double precision without "-r8" switches.
!     LSC_27Jun00
!        0. Work on writes to the main log file, improve readability.
!           Small non-important bug in forgetting path length in ray log,
!           but have not found it yet.  (I remember once understanding the issue.)
!     LSC_21Jun00
!        0. Take out unit 6 writes; DoBern, DoComm,  DoEVsm, DoAdFs
!           DoEVsm is still available through iSMO, see below, but recompile.
!        1. Take out iRememFe nRampOld
!     LSC_12Jun00
!        0. Write Pql Pray and integrals for gnuplot
!     LSC_05Jun00
!        0. Remove phony symmetric hi n_parallel spectrum if DoBram = -1
!           (See 30Apr93)
!        1. Write spectrum to data file for gnuplot (SM ?) post processing
!     LSC_19May00
!        0. Take out DoRFHe, DoRFDa, DoPaus, DoScrn, DoDeBg, DoGrzi, DoHarv
!        1. Take out epol Err Ert Err2 Ert2 (NJ Fisch used these once)
!     LSC_17May00
!        0. Change setpwrl to do by default linear ramp up assuming about
!           100 ramp steps, with flat at the top, but keep in the geometric
!           ramp code, both original (2 * per) and trial (1.eps * per).
!           Invent parameter for NFLATDEF; use NRAMPDIM for default nrampup
!        1. Take out graph of Er Ntht Nrho Eth Err
!        2. Fix bug that diffused power can be negative where current is
!           negative.  Just take abs(Jrf) in the weighting formula.
!     LSC_15May00  -- simplify input and remove features not, or seldom, used
!        0. Split NAMELIST files in two pieces,
!           input.lhh containing  "inpvalue"  --- the most  basic  parameters
!           input.xpt containing  "inpexprt"  --- the experts-only parameters
!        1. Take out the stuff put in for Norton L Bretz ( NLB , bret )
!        2. Take out the stuff put in for Nathaniel J Fisch  (Drt DoFish) .
!        3. Take out the stuff put in for Bob Harvey, but save the code
!        4. Take out the stuff put in for Franco Paoletti, but save the code
!        5. Take out non uniform velocity grids
!     LSC_21Apr00  -- interim fixes, not a stable version
!        0.  Redefine rFudgDmp to be the fractional power lost in a zone
!            in the ray picture, and use that fraction constructing Dql
!            Abandon iFudgDmp, since rFudgDmp needed for all zones, rays.
!        1.  Greatly increase the allowed number of iterations on f_e to
!            200, and set up nflat with kludge to select a build up method:
!              nflat = 1,2 ... nramup-1  ==> old method
!              nflat = nrampup           ==> 1.15 increase each cycle
!              nflat > nrampup           ==> 1/nrampup increase each cycle
!        2.  Change the test for stopping the ray calculation so that
!            the calculation runs longer, and therefore can accomodate
!            large amounts of QL burn through.  It used to be that if
!            Maxwellian damping was "strong" then the ray stopped, but
!            now we say that Maxwellian damping cannot be more than 10%
!            per region, so it takes hundreds of strong damping zone crossings
!            to absorb the ray.
!     LSC_16Jan96
!        0.  Take out references to dpsi.  It was not actually
!        1.  Install the EhstKarney formula.
!     LSC_20Sep95
!        0.  Allow diffusing the rf power with PrfSpred if current has
!            been diffused with DiffuJrf such that
!            Praytot(psi) = Prf-total n(psi) J-diffused(psi) dV(psi) /
!                           [ \Sigma  n(psi) J-diffused(psi) dV(psi) ] *
!                                                 PrfSpred     +
!                  Prf-ray-undiffused(psi)* (1. - PrfSpred)
!        1.  Add phony listing of power to ions in transfer to calling
!            code; this to help with compatibility with TSC/FCD/TRANSP.
!     LSC_14Jul95
!        0.  Warning: Local variables are assumed to be SAVEd by the
!            compiler.  While this is a standard implementation, it is
!            not ANSI standard.  In particular, the NERSC A machine
!            compiler must be forced to do this with an option switch.
!        1.  CALL lsc_EZfini(0,0) consistently; remove 4-transition-point
!            grid stretcher with calls mkgrid eegrd vecadd and vars
!            frv1minus, frv2minus, v1minus, v2minus,
!            frv1plus,  frv2plus,  v1plus,  v2plus,  epsvgr
!
!     LSC_10Feb95
!        0.  Add coupler 'TFTRLHCD' which is 32 0.55cm wgs, 0.15cm septum
!            Move EXCLUDED to 1.2 from 1.5; better for TFTR
!        1.  Repair bug in the error trap on nrays, ntors, npols;
!            add trap of error nparmin>nparmax
!        2.  Replaced "100." with REXPMIN in a couple of  places in FeMaxw;
!            this avoids compiler problem under osf/1 v3.0, f77 v3.6:
!            evaluation of exp( -100. ) gave an arithmetic error.
!            In the TRANSP environment, however, it should be noted that
!            exp( -REXPMIN ) evaluates by silent underflow to zero. --DMC
!        3.  Add time of equilibrium to certain graphs; see EquTitle(81:96)
!     LSC_03Dec94
!        0.  Diffuse Jrf accorting to a method like the Fuchs method.
!            See extensive comments under SmoothJrf.
!            DiffuJrf in the namelist governs this; REPLACES JrfEnhc
!        1.  Add graphs of Jrf by velocity and psi for Bernabei 94 APS meeting
!            see JrfDiagn
!        2.  Move SUBROUTINE LSC to top of file...ahead of comments
!     LSC_30Aug94
!        0.  Enter option for JET_LHCD coupler
!        1.  Initialize iLastRy in Block Data as Fortran requires
!        2.  Remove infinite loop potential in RayIni....
!        3.  Add nslices to parameter output graph
!        4.  Trap read errors in RdTSC
!        5.  Trap powers(i) = 0.0; reset to 0.001
!     LSC_23May94
!        0.  Install Do0Edc so the Edc in lhcdoua/jardin.d etc can be ignored
!            and treated as zero.  This to investigate possible stability
!            problems with T/LSC combined on some shots like 313258
!        1.  Add graph of n-par enhancement on midplane to pitch(r),
!            Ne(r), Pray(r) attached to PlFlg(26)
!     LSC_28Mar94
!        0.  Change the Pitch of the field line graph to degrees; range: +/- 6
!        1.  Add plot of Pray(v) for various psi. Requested by Stefano for the
!            Boulder meeting....to show that PBX is relevant to TPX. Flag 19
!        2.  Also, add plot of Pray(iv) summed over all positions and
!            Plaunch at the edge in v, nparallel, and 1/nparallel(?) space.
!     LSC_10Dec93
!        0.  Allow the use of phony foils to see the Photon Temperature
!            as determined at lower energy might be higher.  Foil code
!            'PF' sets up only Aluminum with the thickness given
!        1.  Allow the introduction of a fast particle population such that
!            T/Tfast = epsT, nfast/n = epsN,  nfast Tfast /(nT) = epsP
!            using variables Tail?eps where ? = N, T, P
!        2.  Apply trapezoid rule for the integration over
!            \nuColl / (Dc + Dql) dv uniformly for +- v and starting and
!            subsequent iterations.
!        3.  Allow for operation on non-symmetric equilibria
!     LSC_18Oct93
!        0.  Allow re-calcluations of rays to take place 1 ray at a time
!            with Do1Rpr=1 (Do 1 Ray per retrace call)  See iLastRy, iDoRay
!        1.  Added entry points to help with TSC restarts:
!            ENTRY PutLSCrestart
!            ENTRY GetLSCrestart
!        2.  Add fast electron population with Doflag DoAdFs.
!            The npsi*nv numbers come from another code GGfe.F
!     LSC_03Sep93
!        0.  Give TSC the power and current !!centered!! on the psi values.
!            Due to a previous misunderstanding this had been given on
!            volume shells ending at the psi values passed.
!        1.  Allow TSC to input directly the new loop voltage and call
!            directly for the new current, and dJ/dE
!              SUBROUTINE JdepFromTSC(curtsc,djdets2,npsiTake) in jrf.F
!              ENTRY PutVinLSC(voltPut,npsiPut)                in grapgrd.F
!              note: location(#)1 is ignored---first data at second location
!              in LSC  input arrays are 1-->NpsiJ;   in TSC 1 -->npsiPut
!              in LSC output arrays are 1-->NpsiJ-1; in TSC 1 -->npsiTake
!     LSC_31Aug93
!        0.  Add EZdscw EZdscr EZopen EZclos to help capture numbers
!            for later manipulation; employ on vertical xray slice
!        1.  Fix error in paramater-list-graph for nrays
!        2.  work xray units: inten in w/cm^2 on camera; w/cm^3 in plasma;
!            including solid angle effect, see "solAng", and absorbers;
!            make frac = ri - i inst of ri - 1 in interpolator; move Cth0 = mu;
!            fix E-indexing; fix multiplication of normalizations; clamp Emax.
!        3.  Add nRunDot, jRunDot == d/dt {n_e_runaway,j_runaway_from_n_dot}
!            = Dql df/dv * {1; q_e v_parallel} | cooperative runaway velocity
!     LSC_19Jun93
!        0.  Add ability to launch rays in the middle, bypass the ql
!            calculation; for TFTR support as discussed with Norton L. Bretz.
!            DoFlag: DoBret; Starting points rstNLB, zstNLB, npaNLB, degNLB
!            ALL THIS UNDER "0." TAKEN OUT MAY 2000; Sorry, Dr Bretz (NLB)
!        1.  Add ability to write out the psi, r, z, n-par, n_rho, etc for
!            Paoletti and Bernabei who are exploring with MatLab the
!            'Volcano Limit' to n-par.
!        2.  Also for Paoletti, add r,z,rtpsi,ompe/omce,Bth,Bph
! !      3.  Fixed bugs in Wr2TSC.F and Rayio.F where Npsi was given instead
!            of NpsiJ as size of normalized flux array, for graphs of
!            approximate damping and dlnJdlnE
!     LSC_24May93
!        0.  Repair work on code that writes data for Bob Harvey.
!            sene was wrong; exde etc in doubt; wdnpar erratic
!        1.  Properly put DOLLAR at end of pline with CALL lsc_to DolAtEnd(pline)
!        2.  Call FastFrac every time at request of TRANSP
!        3.  Add plot of n-par enhancement factor in AccesiSB; puncture plot of
!            actual enhancement
!     LSC_30Apr93
!        0.  Allow for a phony symmetric hi n_parallel spectrum if DoBram = -1
!            This is for an IBW synergy simulation, but not with IBW disp rel.
!        1.  Short circuit warnings if iRayTrsi=0
!     LSC_20Apr93
!        0.  Density fluctuation scattering at bounces at edge; graphs w/ rays
!        1.  Eliminate NparEnhc in namelist; add scatKdeg;
!            add GetKscat and ran3 subroutines; modify BounShft extensively
! !!!    2.  Split the npar array into ntor(NTORDIM) npol(NPOLDIM)
!            These form components of npar(NRAYDIM); add nrays=ntors*npols;
!            ntor replaces npargrid; add npolmax/min;
!        3.  Stopped a ray with more than 4 zones of heavy damping
!            according to a Maxwellian calculation: MxdlnPdz = 0.001
!            This means no penalty for large number of steps.
!        4.  Moved NZONDIM to 2000;
!            removed restriction that NVELDIM should be greater than
!            NPSIDIM + NZONDIM by making WKVDIM = NVELDIM+NPSIDIM+NZONDIM;
!            made plot workspace large enough with MxPLTZON=NPLTDIM+NZONDIM &
!            moving IDM from 200 to NZONDIM.  These requested by Takahashi
!        5.  Changed LSCwarn to be left justified
!     LSC_29Mar93
!        0.  E. Valeo corrected the integral over absorber function,
!            and made other changes to the x ray camera. Constants still to be
!            made correct.  Forced nMUbins to be odd, but not sure this is
!            necessary.  Took out an extra density in the sigma.
!        1.  Small adjustments to absorber code reflecting refinements from
!            S. vonGoeler.  Dimensioning error in GetXmiss was fixed.
!     LSC_15Mar93
!        0.  Improved, added graphs having to do with Xray signal and emission
!        1.  Set ISIZElcfs = 800 because 100 was causing errors
!        2.  Made aspect ratio for coutours of flux and accessibility ok
!     LSC_23Feb93
!        0.  Adjusted Rayio.F on tracing flux surface at 0.3*PsiLim
!     LSC_16Feb93
!        0.  Allowed edge density to be reduced with DoBern=1.  See grapgrd.F
!        1.  Made namelist Inpval an include file at request of TRANSP
!     LSC_06Feb93
!        0.  Improved labeling of launched spectrum graph
!        1.  Trapped and reported possible errors on coupler names
!        2.  Added a very slow coupler, called SLOWSLOW with 16 wgs
!            This not accessible from TRANSP namelist which is integer-driven
! !!!    3.  Corrected error that fghz was not being reported to brambJES !
!     LSC_03Feb93
!        0.  Added plot of power in ray vs root psi, with the first 4 passes
!            at no damping followed in the graph
!        1.  Adjusted step in psi surface plotter to handle D-III-D extremes.
!            Also adjusted start search point to be 0.9 of Rmag
!        2.  Took out DoRFDa and DoRFHe as functioning; left in namelist
!            but took out of include file Doflags.inc
!     LSC_28Jan
!        0.  Made SMALL in fe.F 2.5e-05 for TRANSP use.
!     LSC_22Jan93
!        0.  Allowed Ip as well as ne Te Bo to be changed on DoBern Switch
! !!!    1.  Made adjustment to Dql in zones of high damping so Pql Pray
!            should be pretty much equal.  Key variables: iFudgDmp rFudgDmp
!        2.  Added TORSUPRA coupler option; phaseDeg follows CEN definition
!     LSC_15Jan93
!        0.  Made normalization of parametrized spectrum (DoBram=0) correct.
!        1.  Allowed for Bo to be changed on the DoBern switch
!        2.  Attempted to fix Pray/Pql; by using only discrete n_parallel
!     LSC_08Jan93 :
!        0.  Tokamak de Varennes and PBXM slow couplers selected by character
!            variable couplers(nGrups)
!        1.  dtdV substituted for dVdt in Raystore cause dt can be 0 in TRANSP
!        2.  Ray data written for RW Harvey & CQL3D  with switch DoHarv
!        3.  Grid for reporting Prf and Jrf shifted back to be consistent
!            with comments.  Thanks to Ted.
!        4.  Camera graphs buffer flushing worked on.  Seems ok now.
!     LSC_13Nov92 :
!        0.  copy of 7Nov TSC code for J/E iteration in jrf.F as comments
!     LSC_07Nov92 :
!        0.  SMALL set to 1.5e-5 in fe.F for TRANSP
!        1.  x ray camera w absorbers in XrayCam2 by Valeo, Enright, Ignat
!        2.  last closed flux surface R_, Z_" lcfs " min/max known for graphs;
!        3.  dJ/dE passed to TSC/TRANSP after E/J dJ/dE ;
!        4.  Effective date written with parameters and inputs
!
!subroutine lscsq_main(LhPwrMWx, iRayTrsx, iErrorx)
subroutine lscsq_main(iRayTrsx, iErrorx, arrys)
   use Tracing
   use pl_types
   use iso_c_binding, only : fp => c_double
   use lscsq_mod, only: rmax,rmaj,nrays, nzones,  ind_ray
   use lscsq_mod, only: npols, ntors, npar, ntor,  nth
   use lscsq_mod, only: iendry,iraytrsi, i1stcall
   use lscsq_mod, only: nant, power_inp
   use lscsq_mod, only: ierror!, pe2Vec, NpsiJ
   use lscsq_mod, only: totpwr, diffujrf
   use lscsq_mod, only: iraytrsi
   use lscsq_mod, only:  nz_ind
   !use lscsq_mod, only: lh_out
   use lscsq_mod, only: dql, fe, dfdv
   use lscsq_mod, only: powers_ant, centers_ant, widths_ant, Spec_ant, ntor_ant
   use lscsq_mod, only: spec, powers, centers, widths, nant, npeaks

   
 !  use lscsq_mod, only: powtsc, curtsc, dJdE,dlJdlE
   implicit none

   type(storeAry), dimension(nzones, nrays), intent(inout) :: arrys
   character(len=70) :: ErrMsg
 
   integer :: iry, ity,  i1, i2
   integer :: nzthisray, i, j
 
 !  real(fp), intent(in) :: LhPwrMWx   ! Lower Hybrid power in MW passed from TRANSP
 
   integer, intent(in) :: iRayTrsx    ! 0 use old ray data, and old f_e(v); use new E_dc for the current
                                      ! 1 calculate new rays, and f_e(v), from new equilibrium
                                      ! 2 use old ray data, but calculate new f_e(v)
                                      !   taking account of new n_e and T_e, use new E_dc for the current
  
   integer, intent(out) :: iErrorx    ! 0 LSC finishes without errors; 1 (or more) errors found
                                      !-1 LSC found an error; calling program can keep going
  
   real(fp)  :: PwrFnd, JrfFnd, PqlFnd
   integer :: nwarning

 
 !     Confusion possible over the suffix Vec vec and Ary ary; explained below:
 !     Vec as in PsiVec ... is on the TSC grid, not on the LSC grid: that is Ary
 !     Ary as in PsiAry ... in on the LSC grid, not on the TSC grid: that is Vec
 !
 !     -                                 iLastRy is the index of the last ray
 !     -                                 traced; or -1 if no rays have been
 !     -                                 traced.  iLastRy controls iDoRay
 !                                       iDoRay is a flag in DoRay
 !                                       if 1 thru nrays then trace that ray
 !                                       if otherwise, then trace all rays
 
 !  TotPwr   =  1.0e6_fp*LhPwrMWx 
   
   iRayTrsi =  iRayTrsx
   iError   =  0
   iErrorx  =  0
 
   !if (allocated(lh_out%fe)) then
   !   dql = lh_out%dql
   !   fe  = lh_out%fe
   !   dfdv= lh_out%dfdv
   !endif

 ! Clear warning counter; iEdc flag
   nWarning = 0
   iEndRy = 0
 
   ! the following calls should be made only for the first iteration, I think.
   ! the namelist values do not change, so do not the input profiles
   ! However, the Vloop and electric field do change ...
   ! probably modify the subroutine nml2lsc so that not all arrays are allocated
   ! and initialized in that call.
   ! NOTE: the namelist can change between TRANSP time steps, but cannot change
   ! when LSC is called within an iterative loop during a single LH time step.
 !  if (i1stcall.eq.1) then
 !     call nml2lsc
 !     call lscsq_RdVars
 !  endif
 
   if (iError .GT. 0) then
      CALL lscsq_LSCtrace ( ' RdVars')
      iErrorx = iError
      return
   endif
 
   if(iRayTrsi .EQ. 0 .and. i1stcall .EQ. 1) iRayTrsi =  2
   ! if iRayTrsi=0 and firstcall, then make iRayTrsi=2 so we find f_e etc
 
   if(iRayTrsi .ne. 0 .or. i1stcall .EQ. 1)  then
 !     ifirst = 0
      CALL lscsq_ValInit    ! initialization routines for quasilinear calculations
                           ! -- needed after each TSC read
      ! fmp - not sure. LSC should not re-initialize a Maxwellian distribution
      ! function each time the profiles change, it should use the distribution
      ! function from the previous LH time step
      if (iError .GT. 0) then
         CALL lscsq_LSCtrace ( ' ValInit')
         iErrorx = iError
         return
      endif
   endif
  
   if(iRayTrsi.eq.1 .or. iRayTrsi.eq.2)  then
      do i=1,nant
        TotPwr = power_inp(i)
        
        do j=1,npeaks(i)
           powers(j) = powers_ant(j,i)  
           centers(j)= centers_ant(j,i)  
           widths(j) = widths_ant(j,i)  
           call lscsq_SpecInit
           ! generate launch spectrum with
           ! nGrps humps at centers(nGrps) with widths(nGrps) & relative powers(nGrps).
           ! But if nGrps = 0, then get a Brambilla spectrum using the JE Stevens code.
         enddo
         ntor_ant(:,i) = ntor
         Spec_ant(:,i) = Spec
      enddo
      if (iError .GT. 0) then
         call lscsq_LSCtrace('SpecInit')
         iErrorx = iError
         return
      endif
   endif
   
 
   if (iRayTrsi.ne.0) then
      do i=1,nant
         ntor = ntor_ant(:,i)
         i1 = 1+(i-1)*ntors*npols*nth
         i2 = i*ntors*npols*nth
         do iry=i1,i2  
            ity = ind_ray(1,iry)
            npar(1,iry) = ntor(ity)
         enddo
         call lscsq_rspwr(1.0_fp)
      enddo
      ! apply the launch spectrum to the total power.
      if (iError .GT. 0) then
         CALL lscsq_LSCtrace ( ' PowrInit')
         iErrorx = iError
         return
      endif
   endif
 
   if(iRayTrsi.eq.1) then
      call lscsq_DoRay(arrys)
      if (iError .GT. 0) then
         CALL lscsq_LSCtrace(' DoRay')
         iErrorx = iError
         return
      endif

      do j=1,nrays
 !        if (ok_ray(j).eq.1) then
            nzThisRay=nzones
            do i=1,nzones
               if(arrys(i,j)%ivind .ne. 0) cycle  
               nzThisRay=i
               exit
            enddo
            nz_ind(j) = nzthisray
 !        endif
      enddo
   endif
 
   ! this part needs to be modified, now all rays are in memory
  ! if(iRayTrsi .ne. 1)         then
 !     write(*,*) 'debug rays at entry, print value sof R and Z at zone 10'
 !     write(*,*) 'Rofray:',Rofray(10,:)
 !     write(*,*) 'Zofray:',Zofray(10,:)
  ! endif


   if (iRayTrsi.ne.0) then
      call lscsq_RampPwr(arrys)
      
      if (iError .GT. 0) then
         call lscsq_LSCtrace (' RampPwr')
         iErrorx = iError
         return
      endif
      
      call lscsq_PdepCalc(arrys)    ! compute Prf
   endif
 
   call lscsq_JdepCalc(arrys)
 
   if (iError .GT. 0) then
      CALL lscsq_LSCtrace(' JdepCalc')
      iErrorx = iError
      return
   endif
 
   CALL lscsq_JsplCalc   ! compute Js at E_cd pl(us) dE_dc
 
   if (iError .GT. 0) then
      CALL lscsq_LSCtrace(' JsplCalc')
      iErrorx = iError
      return
   endif
 
   if (diffujrf.gt.0.0) CALL lscsq_SmoJandP(Rmax-Rmaj)
 
   call lscsq_output(PwrFnd, JrfFnd, Pqlfnd)
   write(*,'('' TotPwr PwrFnd JrfFnd PqlFnd: '',4(1pe10.2) )') sum(power_inp),PwrFnd,JrfFnd,PqlFnd
   
   if (abs(sum(power_inp)-PwrFnd).GE.0.2*sum(power_inp).or.abs(PqlFnd-PwrFnd).GE.0.2*sum(power_inp)) then
      write(ErrMsg,'('' TotPwr PwrFnd JrfFnd PqlFnd: '',4(1pe10.2) )') sum(power_inp),PwrFnd,JrfFnd,PqlFnd
      call lscsq_LSCwarn(ErrMsg)
   endif
   ! write a warning if power is not well absorbed or if ray/ql answers are not too close; give current too
  
   if (iError .GT. 0) then
      CALL lscsq_LSCtrace(' Wr2TSC')
      iErrorx = iError
      return
   endif
 
   if (iEndRy .GT. 0) then
      CALL lscsq_LSCwarn ( ' ray stopped early; see neg rc')
      iErrorx =-iEndRy
      return
   endif

   
 
 end subroutine lscsq_main
 
 
 subroutine lscsq_RdVars
   use Tracing
   use iso_c_binding, only : fp => c_double
   use lscsq_mod, only : one, pi, twopi, zel, me_g, vc, zmtocm
   use lscsq_mod, only : erg_eV, zev2keV 
   use lscsq_mod, only : ngrps, powers, begin
   implicit none
 
 
   integer :: i!, ios
   real(fp)::  powsum
   !real(fp) :: restenergy!, keVperERG
 
   !INTEGER :: j
   !INTEGER :: inum
  ! character(len=80) :: message
  ! LOGICAL icheck
 
   begin = 0.
 
   ! fmp - the section below must be removed, but make sure first that everything
   ! is defined
   ! these variables are never used
 !  if(i1stcall .EQ. 1) then
 !     !  (set later)         ifirst= 0
 !     dompesqdn = 4.0_fp*pi*zel**2/me_g
 !     restenergy = me_g*(vc*zmtocm)**2
 !     keVperERG = zev2keV / erg_eV
 !  endif
 
     continue
 
   powsum = 0.00

   do i=1,nGrps
      if(powers(i) .LE. 0.0_fp) then
         powers(i) = 0.001_fp
         CALL lscsq_LSCwarn( ' powers(i) being set to 0.001 ')
      endif
      powsum = powsum + powers(i)
   enddo
 
   !  (DMC) assure normalization of powers array;
   !  total power in MW is passed as argument to main routine LSC
   powers(1:nGrps) = powers(1:nGrps)/powsum
 
 !  if (DoBram .EQ. 1 ) then
 !     do i=1,nGrps
 !        icheck = .FALSE.
 !        do j=1,NCPLDIM
 !           if(couplers(i).eq.cplTyp(j)) icheck = .TRUE.
 !        enddo
 !        if(.not. icheck) then
 !           couplers(i)    = 'PBXMSLOW'
 !           message = ' '
 !           write(message,'(a,i2,a,a)') ' couplers(',i,' being set to ', couplers(i)
 !           CALL lscsq_LSCwarn (message)
 !        endif
 !     enddo
 !  endif
 
 end subroutine lscsq_rdvars
 !
 !     SpecInit begins -------------------------------------------------|
 subroutine lscsq_SpecInit
   use Tracing
   use iso_c_binding, only : fp => c_double
   use lscsq_mod, only : zero, one
   use lscsq_mod, only: ngrpdim, couplers
   use lscsq_mod, only: ntor, npol, spec, ntors, npols,ngrps 
   use lscsq_mod, only: ngrps
   use lscsq_mod, only: powers, centers, widths, fghz, phasedeg
   !use lscsq_mod, only:  nzones
   use lscsq_mod, only: nparmin, nparmax, npolmin, npolmax, dobram, iraytrsi
   use lscsq_mod, only:  nslices, i1stcall
   implicit none
 
   integer, save :: oldNtors, oldNpols, oldnGrps, oldNslic, oldDoBram
 ! real(fp), save :: oldPowrs(ngrps),oldCents(ngrps),oldWidts(ngrps)
 ! real(fp), save :: oldPhsDg(ngrps),oldFreqs(ngrps)
   real(fp) :: oldPowrs(ngrps),oldCents(ngrps),oldWidts(ngrps)
   real(fp) :: oldPhsDg(ngrps),oldFreqs(ngrps)
 !  integer, save :: infreq_diff    ! #of distinct frequencies
 
 ! character(len=8), SAVE :: oldCoupl(ngrps)
   character(len=8) :: oldCoupl(ngrps)
  ! character(len=80) :: pline
   character(len=40) :: MyString
   integer :: igrup, i, ity!, it0!, ipy!, imaxp
   !integer :: ilim1,ilim2
 !  integer :: nRyPrGr
 !  integer :: ifirst =1
   real(fp) :: RayGroup
   real(fp), dimension(ngrps) :: GrpNorm!,FrqPows
 !
 !  DMC Mar 2011: for multi-frequency code
 !  fmp - this part needs to be redesigned. Multi-frequency applies only to
 !  multiple antennae not to multiple groups
  ! integer, dimension(ngrps) :: inrank
  ! integer :: itmp!, iant
   !real(fp), dimension(ngrps) :: npll_min!, npll_max
   !real(fp), dimension(ntors,ngrps) ::  sofns
   !real(fp), dimension(ngrps*ntors) :: sofns1
  ! integer, dimension(ntors,ngrps) :: ifrqs
  ! integer, dimension(ntors*ngrps) :: ifrqs1
 
   real(fp) :: lscsq_SpShape
   !real(fp) :: TstPwr, tmpi
   real(fp) :: excluded = 1.2_fp
 
   ! If this is not the first call, then test to see if any of the parameters
   ! of the spectrum have changed.
   ! If stable, then return to avoid another calculation, which in brambJES is time consuming.
 
   if (i1stcall .EQ. 0 ) then
      i = 0 
      do igrup=1,nGrps
         if( oldFreqs(igrup) .ne. fghz    (igrup)) i=i+1
         if( oldPowrs(igrup) .ne. powers  (igrup)) i=i+1
         if( oldCents(igrup) .ne. centers (igrup)) i=i+1
         if( oldWidts(igrup) .ne. widths  (igrup)) i=i+1
         if( oldPhsDg(igrup) .ne. phaseDeg(igrup)) i=i+1
         if( oldCoupl(igrup) .ne. couplers(igrup)) i=i+1
      enddo
      if(oldNtors.ne.ntors .or. oldnGrps.ne.nGrps .or. oldNpols.ne.npols &
         .or. oldNslic.ne.nslices .or. oldDoBram.ne.DoBram) i=i+1
  
      if (i .EQ. 0) return
  
      iRayTrsi = 1
   else
 !     ifirst = 0
      oldFreqs(1:ngrps) = -99
      oldPowrs(1:ngrps) = -99
      oldCents(1:ngrps) = -99
      oldWidts(1:ngrps) = -99
      oldPhsDg(1:ngrps) = -99
      oldCoupl(1:ngrps) = ' '
      oldNtors = -99
      oldNpols = -99
      oldnGrps = -99
      oldNslic = -99
      oldDoBram= -99
  
      !  DMC Mar 2011: count # of distinct frequencies
      !  if they are nearly the same they count as the same
 !     CALL lscsq_freq_count(infreq_diff)
   endif
 
   ! Make the spectrum in poloidal number
   npol(1) = 0.5_fp*(npolmax+npolmin)
   if (npols.GT.1) call lscsq_ugrid(npol,npols,npolmin,npolmax)

   !     ........................................................................
   ! If not doing a Brambilla calculation, make a model spectrum using gaussian form.
   IF(DoBram .EQ. 0) then
        CALL lscsq_ugridEXC (ntor, ntors, nparmin, nparmax, EXCLUDED)
        !  this is the traditional (single frequency) LSC code.  It can create
        !  some nll values with very low associated power, which likely results
        !  in ray computations that do not affect results...  But, overlapping
        !  nll contributions from different couplers are nicely combined.
        
        do igrup = 1, nGrps
           GrpNorm(igrup) = zero
           do i = 1, ntors
              GrpNorm(igrup) = GrpNorm(igrup) + lscsq_SpShape(ntor(i), centers(igrup), widths(igrup))
           enddo
           GrpNorm(igrup) = one / GrpNorm(igrup)
        enddo
 
        Spec(1:ntors) = zero
        do ity = 1, ntors
           do igrup = 1, nGrps
              RayGroup = GrpNorm(igrup) * powers(igrup)*lscsq_SpShape(ntor(ity), centers(igrup), widths(igrup))
              Spec(ity) = Spec(ity) + RayGroup
           enddo
        enddo
     spec = spec/sum(spec)
 
   elseif (DoBram .EQ. 1) then
 !    ! Parcel out a ray count by groups. Over distribute if a it doesnt come out even
 !     nRyPrGr = ntors / nGrps
 !     if (nGrps * nRyPrGr .LT. ntors) nRyPrgr = nRyPrGr + 1
 !     if (nRyPrGr .GT. NTORDIM) nRyPrGr = NTORDIM
 !     sofns(:,:) = zero
 !     ens(:,:) = zero
 !     do i=1,nGrps
 !        CALL lscsq_BrambJES(fghz(i), phaseDeg(i),                         &
 !                            nRyPrGr, eNs(1,i), SofNs(1,i),              &
 !                            EXCLUDED, nslices, couplers(i), TurnNegs )
 !        ifrqs(1:NTORDIM,i) = i
 !     enddo
 !     do igrup=1,nGrps
 !        do ity=1,nRyPrGr
 !           SofNs(ity,igrup) = powers(igrup)*SofNs(ity,igrup)
 !        enddo
 !     enddo
 !     CALL lscsq_PikSr3NRi(ntors*ngrps, SofNs1, eNs1, ifrqs1)
 !     do i = 1, ntors
 !        ntor(i) = eNs1  (i + NTORDIM*NGRPDIM-ntors )
 !        Spec(i)     = SofNs1(i + NTORDIM*NGRPDIM-ntors )
 !        indx_freq(i) = ifrqs1(i + NTORDIM*NGRPDIM-ntors )
 !     enddo
 !
 !     CALL lscsq_PikSr3NRi(ntors, ntor, Spec, indx_freq)
 !     spec = spec/sum(spec)
   endif
   ! It can happen that there are two rays at the same n_-par. Check and report   
   do i = 2, ntors
      if ( ntor(i) .ne. ntor(i-1) ) cycle    
      write(MyString,*) 'identical nll at ', ntor(i)
      CALL lscsq_LSCwarn(MyString)
   enddo
  
   ! Prepare to return by saving old values 
   do igrup = 1,nGrps
      oldFreqs(igrup) = fghz(igrup)
      oldPowrs(igrup) = powers(igrup)
      oldCents(igrup) = centers(igrup)
      oldWidts(igrup) = widths(igrup)
      oldPhsDg(igrup) = phaseDeg(igrup)
      oldCoupl(igrup) = couplers(igrup)
   enddo
   oldNtors = ntors
   oldNpols = npols
   oldnGrps = nGrps
   oldNslic = Nslices
   oldDoBram= DoBram
 !
  !format(1,'  Ray Index     n-tor  Power(W)  indx_freq')
  !format(1x,     8x,i2,1x,1pe9.2,1x,1pe9.2,4x,i2  )
  
 end subroutine lscsq_SpecInit
 
 !------------------------------------------------
 !
 function lscsq_SpShape(ParValue, center, width)
 
   use iso_c_binding, only : fp => c_double
   implicit none
       
   real(fp) :: lscsq_SpShape
   real(fp), intent(in) :: ParValue, center, width
   real(fp) :: relnpar, exparg
   
   relnpar = (ParValue - center) / width
   exparg = 0.5_fp * relnpar * relnpar
   lscsq_SpShape = exp(-exparg)
 
 end function lscsq_SpShape
 !     -----------------------------------------------------------------|
 subroutine lscsq_ValInit
   use pl_types
   use iso_c_binding, only : fp => c_double
   use lscsq_mod, only : dql, nv, npsi
   implicit none
 
   CALL lscsq_MiscInit    !  doesn't fit elsewhere
   CALL lscsq_MkGrids     !  set up velocity, psi grids needed for quasilinear calculations
   CALL lscsq_ProfInit   !  interpolate Ne, Ni, Te, Ti to the PsiAry grid
   CALL lscsq_FeInit      !  compute constants needed for solution of Fe ....
   CALL lscsq_DqlInit     !  and for solution of Dql
   dql(1:nv,1:npsi,1:2) = 0.0
 
 end subroutine lscsq_ValInit
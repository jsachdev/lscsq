module lscsq_mod

  use iso_c_binding, only: fp => c_double, c_int
  implicit none

  public 

  ! physical constants - this section will be removed once a public library is
  ! defined
  !-----------------------------------------
  real(fp), parameter :: zero = 0.0_fp
  real(fp), parameter :: half = 0.5_fp
  real(fp), parameter :: one  = 1.0_fp
  real(fp), parameter :: two  = 2.0_fp
  !-----------------------------------------

  real(fp), parameter :: pi   = atan2(0.0_fp,-1.0_fp)
  real(fp), parameter :: twopi= PI+PI
  real(fp), parameter :: sqpi = sqrt(pi)
  real(fp), parameter :: kb_K = 1.38064852e-23_fp    !  J/K
  real(fp), parameter :: kb_eV = 1.60217653e-19_fp   !  J/eV
  real(fp), parameter :: erg_eV = 1.60217653e-12_fp  !  erg/eV
  real(fp), parameter :: kb_keV = 1.60217653e-16_fp  !  J/keV
  real(fp), parameter :: qe_eV  = 1.60217653e-19_fp  !  qe (eV) 
  real(fp), parameter :: me_Kg  = 9.10938215e-31_fp  ! electron mass (kg)
  real(fp), parameter :: me_g   = 9.10938215e-28_fp  ! electron mass (g)
  real(fp), parameter :: mp_Kg  = 1.672621637e-27_fp ! proton mass (kg)
  real(fp), parameter :: mp_g   = 1.672621637e-24_fp ! proton mass (g)
  real(fp), parameter :: amu_g  = 1.660539066e-24_fp ! AMU_g 
  !-----------------------------------------
  ! Physics constants, based on CODATA 2006
  real(fp), parameter :: ZEL  = 4.803206799125e-10_fp !ELECTRON CHARGE (STATCOULOMBS)
  real(fp), parameter :: VC   = 2.99792458e8_fp    ! speed of light
  real(fp), parameter :: RMU0 = (4.0e-7_fp)*PI     ! permeability
  real(fp), parameter :: EPS0 = ONE/(VC*VC*RMU0)   ! permittivity
  real(fp), parameter :: usdp = rmu0
  real(fp), parameter :: aee_amp = 9.578833391e+7  ! electron_charge/proton_mass (C*kg^-1)
  real(fp), parameter :: evptocm_sec = 1.384112291e+6_fp !sqrt(2*kb_keV*10^7/mp_g) 
                         !cm/sec for 1eV proton, note J=10^7 erg
  !Define conversion factors
  real(fp), parameter :: zcmtom = 1.0e-2_fp  ! Centimeter to meter
  real(fp), parameter :: zmtocm = 1.0e+2_fp  ! meters to centimeters
  real(fp), parameter :: zcm2tom2 = 1.0e-4_fp  ! cm^2 to m^2
  real(fp), parameter :: zm2tocm2 = 1.0e+4_fp  ! m^2 to cm^2
  real(fp), parameter :: zcm3tom3 = 1.0e-6_fp  ! cm^3 to m^3 | 1/m^3 to 1/cm^3
  real(fp), parameter :: zm3tocm3 = 1.0e+6_fp  ! m^3 to cm^3 | 1/cm^3 to 1/m^3
  real(fp), parameter :: zwtomw = 1.0e-6_fp  ! Watts to Mega-watts
  real(fp), parameter :: zatoma = 1.0e-6_fp  ! Amperes to Mega-amperes
  real(fp), parameter :: zev2kev = 1.0e-3_fp ! eV to keV
  real(fp), parameter :: zkev2ev = 1.0e+3_fp ! keV to eV
  real(fp), parameter :: cgs2mks = 1.0E-7_fp ! CGS TO MKS
  real(fp), parameter :: t2gauss = 1.0E+4_fp ! Tesla to Gauss
  real(fp), parameter :: gauss2t = 1.0E-4_fp ! Gauss to Tesla

  real(fp), parameter :: rad2deg = 180.0_fp/pi
  real(fp), parameter :: deg2rad = pi/180.0_fp

  real(fp), parameter :: epslon = 1.0e-34_fp    ! small number
  real(fp), parameter :: epsinv = 1.0e+34_fp    ! large number

  integer, parameter :: neqs=7
  integer, parameter :: neqsp1=8
  integer, parameter :: NCPLDIM=10      ! max number of coupler *types*
  integer, parameter :: NGRPDIM = 5     ! max number of groups in the spectrum

  integer :: nant  = 1
  integer :: ngrps = 3
  integer :: nstep = 2000 ! max steps in following each ray (500) 
  integer :: npsi = 40
  integer :: nzones = 2000
  integer :: nslices=301 ! num of n_par slices used in Brambilla calc (301) 
  integer :: nv = 199
  integer :: nth= 19
  integer :: nsmoo=9
  integer :: nsmw=3
  integer :: nrampup= 100 ! number of steps to ramp up power 
  integer :: nflat = 10 
  integer :: dobram = 0   ! 1 computes spectrum from JEStevens Brambilla code
                          ! 0 makes a spectrum out of arbitrary Gaussians
  integer :: doxcam = 0      ! 1 give pictures and plots like the 2d x ray camera
  integer :: ntors = 60
  integer :: npols = 1
  integer :: nrays = -1
  integer :: turnnegs = 0

  integer :: i1stcall
 
  integer, dimension(:,:), allocatable :: ind_ray
  real(fp) :: nparmin = 2.5_fp
  real(fp) :: nparmax = 5.5_fp
  real(fp) :: npolmin = -1.0_fp
  real(fp) :: npolmax = 1.0_fp
  real(fp) :: hstplh = 0.0
  real(fp) :: weghtitr=0.2_fp
  real(fp) :: thet0   = 0.0  ! angle of launch 0=> outside midplane, .25=> top
  real(fp) :: dthet = 0.0
  real(fp), dimension(:), allocatable :: thgrid

  real(fp) :: diffujrf = 0.0
  real(fp) :: prfspred = 0.0
  real(fp):: tailteps = 0.0
  real(fp):: tailpeps = 0.0
  real(fp):: tailneps = 0.0
  real(fp):: scatkdeg = 0.0

  real(fp) :: dx_grid=0.0
  real(fp) :: dz_grid=0.0

  real(fp) :: praysum = 0.0
  real(fp) :: pqlsum  = 0.0
  real(fp) :: ppwrsum = 0.0

  real(fp) :: fghz_now = 0.0
  real(fp) :: TotPwr = 0.0
  real(fp) :: enpar= 0.0   ! n_{\parallel} launched for ray being worked at the moment
  real(fp) :: enpol=0.0   ! n_{poloidal}  launched for all rays(zero is good enough)
  real(fp) :: enth =0.0
  real(fp) :: omega=0.0   ! RF frequency (radians/sec)

  real(fp) :: begin=0.0    ! value of path length to begin ray (0 at start)
  real(fp) :: btesl=0.0    ! B_{T}  field at nominal major radius
  real(fp) :: capo2=0.0    ! \omega_{ce}\omega_{ci} / \omega^2
                       !  no LH resonance at any density if .lt. 1
  real(fp) :: ipsq=0.0     ! \sum_i \omega_{pi}^2 / \omega^2
  real(fp) :: ecyc=0.0    ! \omega_{ce} / \omega
  real(fp) :: icyc=0.0    ! \omega_{pi} / \omega
  real(fp) :: epsq=0.0    ! \omega_{pe}^2 / \omega^2
  real(fp) :: ecyc2=0.0   ! ecyc^2
  real(fp) :: psilim=1.0  ! \psi_{lim}  Flux in webers per radian
  real(fp) :: psimin=0.0  ! \psi_{min}
  real(fp) :: kappa=1.0   ! elongation: \kappa
  real(fp) :: qlim=1.0    ! q_{lim}
  real(fp) :: q0=1.0      ! q_{0}
  real(fp) :: rmaj=0.0
  real(fp) :: Rmag=0.0    !  magnetic axis....same as xmag in other commons
  real(fp) :: rmax=0.0    !  outer radius of flux grid
  real(fp) :: rmin=0.0    !  inner radius of flux grid
  real(fp) :: wcei2=0.0   ! \omega_{ce} \omega_{cH} \sum_i ( n_i Z_i^2 / n_e) m_H/m_i
  real(fp) :: woc2=0.0    ! \omega^2/c^2
  real(fp) :: woc4=0.0    ! woc2^2
  real(fp) :: zmin=0.0    ! lower extent of flux grid
  real(fp) :: zmax=0.0    ! upper extent of flux grid

  integer :: lfast = 0
  integer :: lstop = 0
  integer :: iscatplt=0  !index of the scatter event; used for filling inciThet, scatThet
  real(fp) :: delpsi=0.0
  real(fp) :: Te=0.0
  real(fp) :: Ne=0.0
  real(fp) :: pe2min = 0.0

! from lscsq_febins
  integer :: ivZero = 1
  integer :: iITR = 1

  integer :: ismo=1

  real(fp) :: vmin = -1.0_fp
  real(fp) :: vmax = +1.0_fp

  real(fp) :: VthNorm = 0.0
  real(fp) :: fe0=0.0
  real(fp) :: nu0=0.0
  real(fp) :: TailVtrn=0.0

  character(len=8), dimension(ngrpdim) :: couplers
  character(len=8) :: cplTyp(NCPLDIM)

  real(fp), dimension(neqs) ::  f1, f2, f3, y1, y2, y3
  real(fp), dimension(neqsp1) :: f, y 
 
  integer :: iEdc
  real(fp) :: Edcinp

  real(fp) :: d1,d2,d4
  real(fp) :: denom,DKpar,DKper,wdDdw,dDdkABS,epQL, epsL, epsz
 
  real(fp) :: Eper, Epar, Exy, Aion, Aelc, OmEC2, EparI, EparIK, cEparIK
 
  real(fp) :: D11er, D33er, D12er, D11ar, D33ar, D12ar, D11w0, D33w0, D12w0

  integer :: DqlBox(4), nsmsym
  real(fp) :: DcollNorm, nuNorm, DqlHite, Pwrnorm, DqlNorm

  real(fp), dimension(:), allocatable :: powtsc, curtsc, dlJdlE,dJdE
  real(fp), dimension(:,:), allocatable :: sleave, tleave, senter, tenter

  real(fp), dimension(:,:), allocatable :: Dcoll, nucoll
  real(fp), dimension(:), allocatable :: qlsm
  real(fp), dimension(:,:,:), allocatable :: Dql

  real(fp), dimension(:), allocatable :: fghz, powers, centers, widths, phasedeg
  real(fp), dimension(:), allocatable :: fghz_ant, power_inp
  real(fp), dimension(:,:), allocatable :: powers_ant, centers_ant, widths_ant
  real(fp), dimension(:), allocatable :: pwrlevel, FeCvgAry
  real(fp), dimension(:), allocatable :: printgrl, pqintgrl
  real(fp), dimension(:), allocatable :: praytot, pqltot
  real(fp), dimension(:,:), allocatable :: pray, pql

  real(fp), dimension(:), allocatable :: wkzr, wkv, wkpsi, wkzx

  real(fp), dimension(:), allocatable :: ntor, Spec
  real(fp), dimension(:,:), allocatable :: ntor_ant, Spec_ant, npol_ant
  real(fp), dimension(:), allocatable :: npol, scatthet, incithet 

  real(fp), dimension(:), allocatable :: vpar, vtherm, fenorm, VperpSq, &
                                         dvplus,   &
                                         nu0psi, FstFracN, FstFracE

  real(fp), dimension(:,:,:), allocatable :: fe, dfdv


  real(fp), allocatable, dimension(:) :: NeAry, PsiAry, TeAry, ZbrAry, &
                                  iVlAry, dVol, EdcAry, LnlAry, &
                                  BetZAry, MidAry 

  real(fp), allocatable, dimension(:) :: js, jsp, nRunDot, jRunDot, ugr, &
                                         vRunIdx, IrIntgrl, IpIntgrl,    &
                                         vnormPos, vnormNeg, VparMaxP, VparMaxN

  real(fp), allocatable, dimension(:,:) :: Jray 
  
  integer, dimension(:), allocatable :: vnormOK
  integer, dimension(:), allocatable :: npeaks

  real(fp) :: vnorm, nuRuna, gmrun, muminus, muplus
  real(fp) :: dEdcAmnt=1.0e-4
  real(fp) :: vnmax = 0.99 
  integer :: ivrun

! from lscsq_raybins
  integer :: ipsi, iray, iznew, izold, izone   
 
  real(fp) ::  dnpar = 0.0
  real(fp) ::  dtdV = 0.0
  
  integer, dimension(:,:), allocatable :: izind, ivind ! , izcind             
  real(fp), dimension(:,:), allocatable :: rofray, zofray, Pofray,        &
                                         nparry, nperry, rtpsry, powrry,&
                                         TimeRy, DistRy, DetrRy,        &
                                         NeofRy, BthRay, BphRay, rzind

  real(fp), dimension(:,:), allocatable :: ezsq, npar, power, rFudgDmp, &
                                           dlnPds, dlnPdsK, dlnPdsX    

   integer :: nTSCscrn = 6

     integer         iRayTrsi, iXraysi, iError, iEndRy

  integer :: nx = 125
  integer :: nz = 159
  integer :: iplim = 1
  integer :: isym = 0
  integer :: npsij 

  real(fp) :: RlcfsMax, RlcfsMin, ZlcfsMin, ZlcfsMax

  real(fp), dimension(:,:), allocatable :: psigrd  

  REAL(fp) :: RBphi0,  pe2Fac, pi2Fac, pe2Fac14, AioFac, AelFac, ceiFac, OmcFac
  real(fp), allocatable, dimension(:) :: rho, Tekev, pary, ppary, gary, gpary,voltlp,vptemp


  real(fp) :: bgzero, rgzero, xmag, zmag              

  integer(c_int), dimension(:), allocatable :: nz_ind, ok_ray

  real(fp), allocatable, dimension(:) :: dVlVec, iVlVec, pi2Vec, AioVec, AelVec,           &
                          pe2Vec, RBpVec, VprVec, PsiVec, EdcVec, MidVec, TeeVec
      REAL(fp), dimension(4) :: pi2Coefs, AioCoefs, AelCoefs, pe2Coefs, RBphiCoefs, TeeCoefs

  type lscsq_set
       integer :: nant
       integer :: ngrps
       integer :: npeaks
       integer :: nstep
       integer :: nfreq
       integer :: npsi
       integer :: nzones
       integer :: nv
       integer :: nsmoo
       integer :: nsmw
       integer :: nrampup
       integer :: nflat
       integer :: dobram
       integer :: doxcam
       integer :: do1rpr
       integer :: ntors
       integer :: npols
       integer :: nrays
       integer :: turnnegs
       real(fp):: nparmin
       real(fp):: nparmax
       real(fp):: npolmin
       real(fp):: npolmax
       real(fp):: hstplh
       real(fp):: weghtitr
       real(fp):: the0
       real(fp):: diffujrf
       real(fp):: prfspred
       real(fp):: tailteps
       real(fp):: tailpeps
       real(fp):: tailneps
       real(fp):: scatkdeg
       real(fp), allocatable, dimension(:,:) :: powers
       real(fp), allocatable, dimension(:,:) :: centers
       real(fp), allocatable, dimension(:,:) :: widths
       real(fp), allocatable, dimension(:,:) :: phasedeg
       real(fp), allocatable, dimension(:,:) :: fghz  
       character(len=8), allocatable, dimension(:) :: couplers
  end type lscsq_set


  type(lscsq_set), allocatable, dimension(:) :: lscsq_par



  type lh_rays       
    real(fp), dimension(:), allocatable :: freq_ray    
    real(fp), dimension(:,:), allocatable :: dlnPdsK 
    real(fp), dimension(:,:), allocatable :: dlnPdsX
    real(fp), dimension(:,:), allocatable :: ezsq
    real(fp), dimension(:,:,:), allocatable :: fe   
    real(fp), dimension(:,:,:), allocatable :: dfdv 
    real(fp), dimension(:,:,:), allocatable :: dql  
  end type lh_rays   

  type lh_plasma
    real(fp) :: mass
    real(fp) :: chrg
    real(fp) :: Rminlcfs
    real(fp) :: Rmaxlcfs
    real(fp) :: Zminlcfs
    real(fp) :: Zmaxlcfs
    real(fp) :: Rmin
    real(fp) :: Rmax
    real(fp) :: Zmin
    real(fp) :: Zmax
    real(fp) :: dx_grid
    real(fp) :: dz_grid
    real(fp) :: B_axis  
    real(fp) :: Raxis  
    real(fp) :: Zaxis  
    real(fp), dimension(:), allocatable :: powerlh
    real(fp), dimension(:), allocatable :: freqlh
    real(fp), dimension(:), allocatable :: ne
    real(fp), dimension(:), allocatable :: te  
    real(fp), dimension(:), allocatable :: ni
    real(fp), dimension(:), allocatable :: ti  
    real(fp), dimension(:), allocatable :: Vloop
    real(fp), dimension(:), allocatable :: plflx
    real(fp), dimension(:), allocatable :: dvol  
    real(fp), dimension(:), allocatable :: vol  
    real(fp), dimension(:), allocatable :: g_eq  
    real(fp), dimension(:,:), allocatable :: psirz  
  end type lh_plasma


  type(lh_plasma) :: lh_inp
  type(lh_rays) :: lh_out


contains

subroutine lscsq_allocrays

  if(.not.allocated(senter)) allocate(senter(nzones+1,nrays))
  if(.not.allocated(tenter)) allocate(tenter(nzones+1,nrays))
  if(.not.allocated(sleave)) allocate(sleave(nzones+1,nrays))
  if(.not.allocated(tleave)) allocate(tleave(nzones+1,nrays))

  if (.not.allocated(nz_ind)) allocate(nz_ind(nrays))
  nz_ind(1:nrays) = 1
  if (.not.allocated(ok_ray)) allocate(ok_ray(nrays))
  ok_ray(1:nrays) = 1

  if (.not.allocated(thgrid)) allocate(thgrid(nth))
  thgrid(1:nth) = 0.0
  if (.not.allocated(Rofray)) allocate(RofRay(nzones,nrays))
  if (.not.allocated(Zofray)) allocate(ZofRay(nzones,nrays))
  if (.not.allocated(Pofray)) allocate(PofRay(nzones,nrays))
  if (.not.allocated(nparry)) allocate(nparry(nzones,nrays))
  if (.not.allocated(nperry)) allocate(nperry(nzones,nrays))
  if (.not.allocated(rtpsry)) allocate(rtpsry(nzones,nrays))
  if (.not.allocated(powrry)) allocate(powrry(nzones,nrays))
  if (.not.allocated(timery)) allocate(timery(nzones,nrays))
  if (.not.allocated(distry)) allocate(distry(nzones,nrays))
  if (.not.allocated(detrry)) allocate(detrry(nzones,nrays))
  if (.not.allocated(neofry)) allocate(neofry(nzones,nrays))
  if (.not.allocated(bthray)) allocate(bthray(nzones,nrays))
  if (.not.allocated(bphray)) allocate(bphray(nzones,nrays))
  if (.not.allocated(izind)) allocate(izind(nzones,nrays))
  if (.not.allocated(rzind)) allocate(rzind(nzones,nrays))
  rzind(1:nzones,1:nrays) = 0.0
  if (.not.allocated(ivind)) allocate(ivind(nzones,nrays))
  if (.not.allocated(ezsq)) allocate(ezsq(nzones,nrays))
  if (.not.allocated(npar)) allocate(npar(nzones,nrays))
  if (.not.allocated(power)) allocate(power(nzones,nrays))
  if (.not.allocated(dlnPds)) allocate(dlnPds(nzones,nrays))
  dlnPds(1:nzones,1:nrays) = 0.0
  if (.not.allocated(dlnPdsk)) allocate(dlnPdsk(nzones,nrays))
!  dlnPdsk(1:nzones,1:nrays) = 0.0
  if (.not.allocated(dlnPdsx)) allocate(dlnPdsx(nzones,nrays))
!  dlnPdsx(1:nzones,1:nrays) = 0.0
  if (.not.allocated(rfudgdmp)) allocate(rfudgdmp(nzones,nrays))
  rfudgdmp(1:nzones,1:nrays) = 0.0


  if(.not.allocated(ntor_ant)) allocate(ntor_ant(ntors,nant))
  ntor_ant(1:ntors,1:nant) = 0.0
  if(.not.allocated(spec_ant)) allocate(spec_ant(ntors,nant))
  spec_ant(1:ntors,1:nant) = 0.0
  if(.not.allocated(npol_ant)) allocate(npol_ant(npols,nant))
  npol_ant(1:npols,1:nant) = 0.0
  if(.not.allocated(ntor)) allocate(ntor(ntors))
  ntor(1:ntors) = 0.0
  if(.not.allocated(spec)) allocate(spec(ntors))
  spec(1:ntors) = 0.0
  if(.not.allocated(npol)) allocate(npol(npols))
  npol(1:npols) = 0.0
  if(.not.allocated(scatthet)) allocate(scatthet(npols))
  scatthet(1:npols) = 0.0
  if(.not.allocated(incithet)) allocate(incithet(npols))
  incithet(1:npols) = 0.0

  if(.not.allocated(lh_out%dlnPdsK)) allocate(lh_out%dlnPdsK(nzones,nrays))
  if(.not.allocated(lh_out%dlnPdsX)) allocate(lh_out%dlnPdsX(nzones,nrays))
  if(.not.allocated(lh_out%ezsq)) allocate(lh_out%ezsq(nzones,nrays))

  if(.not.allocated(lh_out%fe)) allocate(lh_out%fe(nv,npsi,2))
  if(.not.allocated(lh_out%dfdv)) allocate(lh_out%dfdv(nv,npsi,2))
  if(.not.allocated(lh_out%dql)) allocate(lh_out%dql(nv,npsi,2))

end subroutine lscsq_allocrays

subroutine lscsq_alloc

  if(.not.allocated(neary)) allocate(neary(npsi))
  neary(1:npsi) = 0.0
  if(.not.allocated(Teary)) allocate(Teary(npsi))
  Teary(1:npsi) = 0.0
  if(.not.allocated(psiary)) allocate(psiary(npsi))
  psiary(1:npsi) = 0.0
  if(.not.allocated(midary)) allocate(midary(npsi))
  midary(1:npsi) = 0.0
  if(.not.allocated(zbrary)) allocate(zbrary(npsi))
  zbrary(1:npsi) = 0.0
  if(.not.allocated(Edcary)) allocate(Edcary(npsi))
  Edcary(1:npsi) = 0.0
  if(.not.allocated(Lnlary)) allocate(Lnlary(npsi))
  Lnlary(1:npsi) = 0.0
  if(.not.allocated(iVlary)) allocate(iVlary(npsi))
  iVlary(1:npsi) = 0.0
  if(.not.allocated(dVol)) allocate(dVol(npsi))
  dVol(1:npsi) = 0.0
  if(.not.allocated(betZary)) allocate(betZary(npsi))
  betZary(1:npsi) = 0.0

  if(.not.allocated(qlsm)) allocate(qlsm(nv))
  if(.not.allocated(dql)) allocate(dql(nv,npsi,2))
  if(.not.allocated(dcoll)) allocate(dcoll(nv,npsi))
  if(.not.allocated(nucoll)) allocate(nucoll(nv,npsi))


  if(.not.allocated(vpar)) allocate(vpar(nv))
  vpar(1:nv) = 0.0  
  if(.not.allocated(dvplus)) allocate(dvplus(nv))
  if(.not.allocated(vtherm)) allocate(vtherm(npsi))
  vtherm(1:npsi) = 0.0
  if(.not.allocated(Fenorm)) allocate(Fenorm(npsi))
  if(.not.allocated(nu0psi)) allocate(nu0psi(npsi))
  nu0psi(1:npsi) = 0.0  
  if(.not.allocated(FstFracn)) allocate(FstFracn(npsi))
  fstfracn(1:npsi) = 0.0  
  if(.not.allocated(FstFracE)) allocate(FstFracE(npsi))
  fstfrace(1:npsi) = 0.0  
  if(.not.allocated(vperpsq)) allocate(vperpsq(npsi))
  vperpsq(1:npsi) = 0.0  
  if(.not.allocated(fe)) allocate(fe(nv,npsi,2))
  if(.not.allocated(dfdv)) allocate(dfdv(nv,npsi,2))


  if(.not.allocated(fghz_ant)) allocate(fghz_ant(nant))
  if(.not.allocated(power_inp)) allocate(power_inp(nant))
  if(.not.allocated(npeaks)) allocate(npeaks(nant))
  if(.not.allocated(powers_ant)) allocate(powers_ant(ngrps,nant))
  if(.not.allocated(centers_ant)) allocate(centers_ant(ngrps,nant))
  if(.not.allocated(widths_ant)) allocate(widths_ant(ngrps,nant))

  if(.not.allocated(powers)) allocate(powers(ngrps))
  if(.not.allocated(centers)) allocate(centers(ngrps))
  if(.not.allocated(widths)) allocate(widths(ngrps))
  if(.not.allocated(phasedeg)) allocate(phasedeg(ngrps))

  if(.not.allocated(pray)) allocate(pray(nv,npsi))
  if(.not.allocated(pql)) allocate(pql(nv,npsi))
  if(.not.allocated(praytot)) allocate(praytot(npsi))
  if(.not.allocated(pqltot)) allocate(pqltot(npsi))
  if(.not.allocated(pqintgrl)) allocate(pqintgrl(npsi))
  if(.not.allocated(printgrl)) allocate(printgrl(npsi))
  if(.not.allocated(vnormOK)) allocate(vnormOK(npsi))
  if(.not.allocated(js)) allocate(js(npsi))
  if(.not.allocated(jsp)) allocate(jsp(npsi))
  if(.not.allocated(nrundot)) allocate(nrundot(npsi))
  if(.not.allocated(Jrundot)) allocate(Jrundot(npsi))
  if(.not.allocated(vrunidx)) allocate(vrunidx(npsi))
  if(.not.allocated(irintgrl)) allocate(irintgrl(npsi))
  irintgrl(1:npsi) = 0.0 
  if(.not.allocated(ipintgrl)) allocate(ipintgrl(npsi))
  ipintgrl(1:npsi) = 0.0 
  if(.not.allocated(vnormpos)) allocate(vnormpos(npsi))
  vnormpos(1:npsi) = 0.0 
  if(.not.allocated(vnormneg)) allocate(vnormneg(npsi))
  vnormneg(1:npsi) = 0.0 
  if(.not.allocated(vparmaxp)) allocate(vparmaxp(npsi))
  vparmaxp(1:npsi) = 0.0 
  if(.not.allocated(vparmaxn)) allocate(vparmaxn(npsi))
  vparmaxn(1:npsi) = 0.0 
  if(.not.allocated(jray)) allocate(jray(nv,npsi))
  if(.not.allocated(ugr)) allocate(ugr(nv))
  ugr(1:nv) = 0.0 

  if(.not.allocated(pwrlevel)) allocate(pwrlevel(nrampup))
  if(.not.allocated(fecvgary)) allocate(fecvgary(nrampup))

  if(.not.allocated(wkv)) allocate(wkv(nv))

end subroutine lscsq_alloc


subroutine lscsq_constants

  implicit none

  !-----------------------------------------
  real(fp), parameter :: zero = 0.0_fp
  real(fp), parameter :: half = 0.5_fp
  real(fp), parameter :: one  = 1.0_fp
  real(fp), parameter :: two  = 2.0_fp
  !-----------------------------------------
  real(fp), parameter :: pi   = atan2(0.0_fp,-1.0_fp)
  real(fp), parameter :: twopi= PI+PI
  real(fp), parameter :: sqpi = sqrt(pi)
  real(fp), parameter :: me_Kg  = 9.10938215e-31_fp  ! electron mass (kg)
  real(fp), parameter :: me_g   = 9.10938215e-28_fp  ! electron mass (g)
  real(fp), parameter :: mp_Kg  = 1.672621637e-27_fp ! proton mass (kg)
  real(fp), parameter :: mp_g   = 1.672621637e-24_fp ! proton mass (g)
  real(fp), parameter :: qe_eV  = 1.60217653e-19_fp  !  qe (eV) 
  real(fp), parameter :: vc     = 2.99792458e8_fp    ! speed of light
  real(fp), parameter :: zcmtom = 1.0e-2_fp  ! Centimeter to meter
  real(fp), parameter :: zmtocm = 1.0e+2_fp  ! meters to centimeters
  real(fp), parameter :: zcm2tom2 = 1.0e-4_fp  ! cm^2 to m^2
  real(fp), parameter :: zm2tocm2 = 1.0e+4_fp  ! m^2 to cm^2
  real(fp), parameter :: zcm3tom3 = 1.0e-6_fp  ! cm^3 to m^3 | 1/m^3 to 1/cm^3
  real(fp), parameter :: zm3tocm3 = 1.0e+6_fp  ! m^3 to cm^3 | 1/cm^3 to 1/m^3
  real(fp), parameter :: zev2kev = 1.0e-3_fp ! eV to keV
  real(fp), parameter :: zkev2ev = 1.0e+3_fp ! keV to eV
  real(fp), parameter :: cgs2mks = 1.0E-7_fp ! CGS TO MKS
  real(fp), parameter :: t2gauss = 1.0E+4_fp ! Tesla to Gauss
  real(fp), parameter :: gauss2t = 1.0E-4_fp ! Gauss to Tesla

  real(fp), parameter :: rad2deg = 180.0_fp/pi
  real(fp), parameter :: deg2rad = pi/180.0_fp

end subroutine lscsq_constants

subroutine alloc_profs

  implicit none

  if(.not.allocated(dvlvec)) allocate(dvlvec(npsij))
  if(.not.allocated(ivlvec)) allocate(ivlvec(npsij))
  if(.not.allocated(pi2vec)) allocate(pi2vec(npsij))
  if(.not.allocated(aiovec)) allocate(aiovec(npsij))
  if(.not.allocated(aelvec)) allocate(aelvec(npsij))
  if(.not.allocated(pe2vec)) allocate(pe2vec(npsij))
  if(.not.allocated(edcvec)) allocate(edcvec(npsij))
  if(.not.allocated(midvec)) allocate(midvec(npsij))
  if(.not.allocated(voltlp)) allocate(voltlp(npsij))
  if(.not.allocated(psivec)) allocate(psivec(npsij))
  if(.not.allocated(vprvec)) allocate(vprvec(npsij))
  if(.not.allocated(rbpvec)) allocate(rbpvec(npsij))

  if(.not.allocated(TekeV)) allocate(TekeV(npsij))
  if(.not.allocated(powtsc)) allocate(powtsc(npsij))
  if(.not.allocated(curtsc)) allocate(curtsc(npsij))
  if(.not.allocated(dJdE)) allocate(dJdE(npsij))
  if(.not.allocated(dlJdlE)) allocate(dlJdlE(npsij))

end subroutine alloc_profs

!--------------------------------------------------------------


subroutine dealloc_profs

implicit none

  deallocate(dvlvec,ivlvec,pi2vec,aiovec,aelvec,pe2vec)
  deallocate(edcvec,midvec,voltlp,psivec,vprvec,rbpvec)
  deallocate(Tekev,powtsc,curtsc,dJdE,dlJdlE)


end subroutine dealloc_profs



end module lscsq_mod
!-----------------------





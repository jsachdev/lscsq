  module lscsq_mod

  use iso_c_binding, only: fp => c_double, c_int
  implicit none

  public 

  ! physical constants - this section will be removed once a public library is
  ! defined
  !-----------------------------------------
  Real(fp), parameter :: zero = 0.0_fp
  Real(fp), parameter :: half = 0.5_fp
  Real(fp), parameter :: one  = 1.0_fp
  Real(fp), parameter :: two  = 2.0_fp
  !-----------------------------------------

  Real(fp), parameter :: pi   = atan2(0.0_fp,-1.0_fp)
  Real(fp), parameter :: twopi= PI+PI
  Real(fp), parameter :: sqpi = sqrt(pi)
  Real(fp), parameter :: kb_K = 1.38064852e-23_fp    !  J/K
  Real(fp), parameter :: kb_eV = 1.60217653e-19_fp   !  J/eV
  Real(fp), parameter :: erg_eV = 1.60217653e-12_fp  !  erg/eV
  Real(fp), parameter :: kb_keV = 1.60217653e-16_fp  !  J/keV
  Real(fp), parameter :: qe_eV  = 1.60217653e-19_fp  !  qe (eV) 
  Real(fp), parameter :: me_Kg  = 9.10938215e-31_fp  ! electron mass (kg)
  Real(fp), parameter :: me_g   = 9.10938215e-28_fp  ! electron mass (g)
  Real(fp), parameter :: mp_Kg  = 1.672621637e-27_fp ! proton mass (kg)
  Real(fp), parameter :: mp_g   = 1.672621637e-24_fp ! proton mass (g)
  Real(fp), parameter :: amu_g  = 1.660539066e-24_fp ! AMU_g 
  !-----------------------------------------
  ! Physics constants, based on CODATA 2006
  Real(fp), parameter :: ZEL  = 4.803206799125e-10_fp !ELECTRON CHARGE (STATCOULOMBS)
  Real(fp), parameter :: VC   = 2.99792458e8_fp    ! speed of light
  Real(fp), parameter :: RMU0 = (4.0e-7_fp)*PI     ! permeability
  Real(fp), parameter :: EPS0 = ONE/(VC*VC*RMU0)   ! permittivity
  Real(fp), parameter :: usdp = rmu0
  Real(fp), parameter :: aee_amp = 9.578833391e+7_fp  ! electron_charge/proton_mass (C*kg^-1)
  Real(fp), parameter :: evptocm_sec = 1.384112291e+6_fp !sqrt(2*kb_keV*10^7/mp_g) 
                         !cm/sec for 1eV proton, note J=10^7 erg
  !Define conversion factors
  Real(fp), parameter :: zcmtom = 1.0e-2_fp  ! Centimeter to meter
  Real(fp), parameter :: zmtocm = 1.0e+2_fp  ! meters to centimeters
  Real(fp), parameter :: zcm2tom2 = 1.0e-4_fp  ! cm^2 to m^2
  Real(fp), parameter :: zm2tocm2 = 1.0e+4_fp  ! m^2 to cm^2
  Real(fp), parameter :: zcm3tom3 = 1.0e-6_fp  ! cm^3 to m^3 | 1/m^3 to 1/cm^3
  Real(fp), parameter :: zm3tocm3 = 1.0e+6_fp  ! m^3 to cm^3 | 1/cm^3 to 1/m^3
  Real(fp), parameter :: zwtomw = 1.0e-6_fp  ! Watts to Mega-watts
  Real(fp), parameter :: zatoma = 1.0e-6_fp  ! Amperes to Mega-amperes
  Real(fp), parameter :: zev2kev = 1.0e-3_fp ! eV to keV
  Real(fp), parameter :: zkev2ev = 1.0e+3_fp ! keV to eV
  Real(fp), parameter :: cgs2mks = 1.0e-7_fp ! CGS TO MKS
  Real(fp), parameter :: t2gauss = 1.0e+4_fp ! Tesla to Gauss
  Real(fp), parameter :: gauss2t = 1.0e-4_fp ! Gauss to Tesla

  Real(fp), parameter :: rad2deg = 180.0_fp/pi
  Real(fp), parameter :: deg2rad = pi/180.0_fp

  Real(fp), parameter :: epslon = 1.0e-34_fp    ! small number
  Real(fp), parameter :: epsinv = 1.0e+34_fp    ! large number
  integer, parameter :: nRayQt = 6
  integer, parameter :: neqs=7
  integer, parameter :: neqsp1=8
  Real(fp), parameter :: corrector_tol=0.01
  integer, parameter :: nCoeff = 5
  integer, parameter :: nCorrectIter = 5 ! number of iteration for corrector in 
	! Predictor-Corrector integrator
  !character(40), parameter :: predictorName = 'Adams 5'
  !character(40), parameter :: correctorName = 'Adams 6'

  !Real(fp), parameter :: aCorrect(1) = 1.0_fp
  Real(fp), parameter :: bCorrect_1 = 475.0_fp / 1440.0_fp
  Real(fp), parameter :: bCorrect_2 = 1427.0_fp / 1440.0_fp
  Real(fp), parameter :: bCorrect_3 = -798.0_fp / 1440.0_fp
  Real(fp), parameter :: bCorrect_4 = 482.0_fp / 1440.0_fp
  Real(fp), parameter :: bCorrect_5 = -173.0_fp / 1440.0_fp
  Real(fp), parameter :: bCorrect_6 = 27.0_fp / 1440.0_fp

  !aPredict(1) = 1.0_fp
  Real(fp), parameter :: bPredict_1 = 1901.0_fp / 720.0_fp
  Real(fp), parameter :: bPredict_2 = -2774.0_fp  / 720.0_fp
  Real(fp), parameter :: bPredict_3 = 2616.0_fp  / 720.0_fp
  Real(fp), parameter :: bPredict_4 = -1274.0_fp  / 720.0_fp
  Real(fp), parameter :: bPredict_5 = 251.0_fp  / 720.0_fp

  !Parameters used in rayStore.F90/lscsq_E2byPr
  ! Real(fp), parameter :: NparDead = 9.8_fp 
  ! Real(fp), parameter :: NperDead = 147.0_fp
  ! Real(fp), parameter :: MaxwDead = 1.0e-12_fp
  ! Real(fp), parameter :: MxdPdz   = 0.90_fp
  ! Real(fp), parameter :: SMALL    = 1.0e-30_fp

  !Parameters used in rayIni.F90/lscsq_RayIni
  Real(fp), parameter :: SEARCHINCR = 5.0e-03_fp 

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
  Real(fp) :: nparmin = 2.5_fp
  Real(fp) :: nparmax = 5.5_fp
  Real(fp) :: npolmin = -1.0_fp
  Real(fp) :: npolmax = 1.0_fp
  Real(fp) :: hstplh = 0.0_fp
  Real(fp) :: weghtitr=0.2_fp
  Real(fp) :: thet0   = 0.0_fp  ! angle of launch 0=> outside midplane, .25=> top
  Real(fp) :: dthet = 0.0_fp
  Real(fp), dimension(:), allocatable :: thgrid

  Real(fp) :: diffujrf = 0.0_fp
  Real(fp) :: prfspred = 0.0_fp
  Real(fp):: tailteps = 0.0_fp
  Real(fp):: tailpeps = 0.0_fp
  Real(fp):: tailneps = 0.0_fp
  Real(fp):: scatkdeg = 0.0_fp

  Real(fp) :: dx_grid=0.0_fp
  Real(fp) :: dz_grid=0.0_fp

  Real(fp) :: praysum = 0.0_fp
  Real(fp) :: pqlsum  = 0.0_fp
  Real(fp) :: ppwrsum = 0.0_fp

  Real(fp) :: TotPwr = 0.0_fp

  Real(fp) :: begin=0.0_fp    ! value of path length to begin ray (0 at start)
  Real(fp) :: btesl=0.0_fp    ! B_{T}  field at nominal major radius
  Real(fp) :: capo2=0.0_fp    ! \omega_{ce}\omega_{ci} / \omega^2
                       !  no LH resonance at any density if .lt. 1
  Real(fp) :: psilim=1.0_fp  ! \psi_{lim}  Flux in webers per radian
  Real(fp) :: psimin=0.0_fp  ! \psi_{min}
  Real(fp) :: rmaj=0.0_fp
  Real(fp) :: Rmag=0.0_fp    !  magnetic axis....same as xmag in other commons
  Real(fp) :: rmax=0.0_fp    !  outer radius of flux grid
  Real(fp) :: rmin=0.0_fp    !  inner radius of flux grid
  Real(fp) :: wcei2=0.0_fp   ! \omega_{ce} \omega_{cH} \sum_i ( n_i Z_i^2 / n_e) m_H/m_i
  Real(fp) :: zmin=0.0_fp    ! lower extent of flux grid
  Real(fp) :: zmax=0.0_fp    ! upper extent of flux grid

  integer :: lfast = 0
  integer :: iscatplt=0  !index of the scatter event; used for filling inciThet, scatThet
  Real(fp) :: delpsi=0.0_fp
  Real(fp) :: Te=0.0_fp
  Real(fp) :: Ne=0.0_fp

  integer :: ivZero = 1
  integer :: iITR = 1
  integer :: ismo=1

  Real(fp) :: vmin = -1.0_fp
  Real(fp) :: vmax = +1.0_fp

  Real(fp) :: VthNorm = 0.0_fp
  Real(fp) :: fe0=0.0_fp
  Real(fp) :: nu0=0.0_fp
  Real(fp) :: TailVtrn=0.0_fp

  character(len=8), dimension(ngrpdim) :: couplers
  character(len=8) :: cplTyp(NCPLDIM)
 
  integer :: iEdc
  Real(fp) :: Edcinp
 
  Real(fp) :: OmEC2, EparIK, cEparIK

  integer :: DqlBox(4), nsmsym
  Real(fp) :: DcollNorm, nuNorm, DqlHite, Pwrnorm, DqlNorm

  Real(fp), dimension(:), allocatable :: powtsc, curtsc, dlJdlE,dJdE
 
  Real(fp), dimension(:,:), allocatable :: Dcoll, nucoll
  Real(fp), dimension(:), allocatable :: qlsm
  Real(fp), dimension(:,:,:), allocatable :: Dql

  Real(fp), dimension(:), allocatable :: fghz, powers, centers, widths, phasedeg
  Real(fp), dimension(:), allocatable :: fghz_ant, power_inp
  Real(fp), dimension(:,:), allocatable :: powers_ant, centers_ant, widths_ant
  Real(fp), dimension(:), allocatable :: pwrlevel, FeCvgAry
  Real(fp), dimension(:), allocatable :: printgrl, pqintgrl
  Real(fp), dimension(:), allocatable :: praytot, pqltot
  Real(fp), dimension(:,:), allocatable :: pray, pql

  Real(fp), dimension(:), allocatable :: wkzr, wkv, wkpsi, wkzx

  Real(fp), dimension(:), allocatable :: ntor, Spec
  Real(fp), dimension(:,:), allocatable :: ntor_ant, Spec_ant, npol_ant
  Real(fp), dimension(:), allocatable :: npol, scatthet, incithet 

  Real(fp), dimension(:), allocatable :: vpar, vtherm, fenorm, VperpSq, &
                                         dvplus,   &
                                         nu0psi, FstFracN, FstFracE

  Real(fp), dimension(:,:,:), allocatable :: fe, dfdv


  Real(fp), allocatable, dimension(:) :: NeAry, PsiAry, TeAry, ZbrAry, &
                                  iVlAry, dVol, EdcAry, LnlAry, &
                                  BetZAry, MidAry 

  Real(fp), allocatable, dimension(:) :: js, jsp, nRunDot, jRunDot, ugr, &
                                         vRunIdx, IrIntgrl, IpIntgrl,    &
                                         vnormPos, vnormNeg, VparMaxP, VparMaxN

  Real(fp), allocatable, dimension(:,:) :: Jray 
  
  integer, dimension(:), allocatable :: vnormOK
  integer, dimension(:), allocatable :: npeaks

  Real(fp) :: vnorm, nuRuna, gmrun, muminus, muplus
  Real(fp) :: dEdcAmnt=1.0e-4_fp
  Real(fp) :: vnmax = 0.99_fp 
  integer :: ivrun

  integer :: ipsi, iznew, izold 
 
  Real(fp) ::  dnpar = 0.0_fp
  Real(fp) ::  dtdV = 0.0_fp

  Real(fp), dimension(:,:), allocatable :: npar, power, rFudgDmp!, &
  
  integer :: nTSCscrn = 6

  integer :: iRayTrsi, iXraysi, iError, iEndRy

  integer :: nx = 125
  integer :: nz = 159
  integer :: iplim = 1
  integer :: isym = 0
  integer :: npsij 

  Real(fp) :: RlcfsMax, RlcfsMin, ZlcfsMin, ZlcfsMax

  Real(fp), dimension(:,:), allocatable :: psigrd  

  Real(fp) :: RBphi0  , pi2Fac, pe2Fac14, AioFac, AelFac, ceiFac
  Real(fp) :: OmcFac, pe2Fac 
  Real(fp), allocatable, dimension(:) :: rho, Tekev, pary, ppary, gary, gpary,voltlp,vptemp


  Real(fp) :: bgzero, rgzero, xmag, zmag              

  integer(c_int), dimension(:), allocatable :: nz_ind, ok_ray

  Real(fp), allocatable, dimension(:) :: dVlVec, iVlVec, pi2Vec, AioVec, AelVec,           &
                          pe2Vec, RBpVec, VprVec, PsiVec, EdcVec, MidVec, TeeVec
  
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
       Real(fp):: nparmin
       Real(fp):: nparmax
       Real(fp):: npolmin
       Real(fp):: npolmax
       Real(fp):: hstplh
       Real(fp):: weghtitr
       Real(fp):: the0
       Real(fp):: diffujrf
       Real(fp):: prfspred
       Real(fp):: tailteps
       Real(fp):: tailpeps
       Real(fp):: tailneps
       Real(fp):: scatkdeg
       Real(fp), allocatable, dimension(:,:) :: powers
       Real(fp), allocatable, dimension(:,:) :: centers
       Real(fp), allocatable, dimension(:,:) :: widths
       Real(fp), allocatable, dimension(:,:) :: phasedeg
       Real(fp), allocatable, dimension(:,:) :: fghz  
       character(len=8), allocatable, dimension(:) :: couplers
  end type lscsq_set


  type(lscsq_set), allocatable, dimension(:) :: lscsq_par

  type lh_rays       
    Real(fp), dimension(:), allocatable :: freq_ray    
    Real(fp), dimension(:,:), allocatable :: dlnPdsK 
    Real(fp), dimension(:,:), allocatable :: dlnPdsX
    Real(fp), dimension(:,:), allocatable :: ezsq
    Real(fp), dimension(:,:,:), allocatable :: fe   
    Real(fp), dimension(:,:,:), allocatable :: dfdv 
    Real(fp), dimension(:,:,:), allocatable :: dql  
  end type lh_rays   

  type lh_plasma
    Real(fp) :: mass
    Real(fp) :: chrg
    Real(fp) :: Rminlcfs
    Real(fp) :: Rmaxlcfs
    Real(fp) :: Zminlcfs
    Real(fp) :: Zmaxlcfs
    Real(fp) :: Rmin
    Real(fp) :: Rmax
    Real(fp) :: Zmin
    Real(fp) :: Zmax
    Real(fp) :: dx_grid
    Real(fp) :: dz_grid
    Real(fp) :: B_axis  
    Real(fp) :: Raxis  
    Real(fp) :: Zaxis  
    Real(fp), dimension(:), allocatable :: powerlh
    Real(fp), dimension(:), allocatable :: freqlh 
    Real(fp), dimension(:), allocatable :: ne
    Real(fp), dimension(:), allocatable :: te  
    Real(fp), dimension(:), allocatable :: ni
    Real(fp), dimension(:), allocatable :: ti  
    Real(fp), dimension(:), allocatable :: Vloop
    Real(fp), dimension(:), allocatable :: plflx
    Real(fp), dimension(:), allocatable :: dvol  
    Real(fp), dimension(:), allocatable :: vol  
    Real(fp), dimension(:), allocatable :: g_eq  
    Real(fp), dimension(:,:), allocatable :: psirz  
  end type lh_plasma

  type(lh_plasma) :: lh_inp
  !type(lh_rays) :: lh_out


  !$acc declare create(psimin, psilim, hstplh, nzones, nrays)
  !$acc declare create(OmcFac, zmin, rmin, ivzero, pe2Fac , cEparIK, delpsi)
  !$acc declare create(psigrd, pe2vec, pi2vec, teKev, aiovec, aelvec, RBpvec, psivec)
  !$acc declare create(dvol, vpar, psigrd)
  !!$acc declare create(twopi, vc, deg2rad,ntor, npol, thgrid, ind_ray, lfast)
  !!$acc declare create(begin, rmax, zmax, rmaj, SEARCHINCR, fghz, npsij, ok_ray, npar)

contains

subroutine lscsq_allocrays
  use iso_c_binding, only : fp => c_double
  implicit none
  integer :: i, j
  if (.not.allocated(nz_ind)) allocate(nz_ind(nrays))
  if (.not.allocated(ok_ray)) allocate(ok_ray(nrays))
 
  do i=1,nrays
    nz_ind(i) = 1
    ok_ray(i) = 1
  enddo

  if (.not.allocated(thgrid)) allocate(thgrid(nth))

  do i=1,nth
    thgrid(i) = 0.0_fp
  enddo

  if (.not.allocated(npar)) allocate(npar(nzones,nrays))
  if (.not.allocated(power)) allocate(power(nzones,nrays))
  if (.not.allocated(rfudgdmp)) allocate(rfudgdmp(nzones,nrays))
  rfudgdmp(1:nzones,1:nrays) = 0.0_fp

  if(.not.allocated(ntor_ant)) allocate(ntor_ant(ntors,nant))
  if(.not.allocated(spec_ant)) allocate(spec_ant(ntors,nant))
  if(.not.allocated(npol_ant)) allocate(npol_ant(npols,nant))

  do i=1,ntors
    do J=1,nant
      ntor_ant(I,J) = 0.0_fp
      spec_ant(I,J) = 0.0_fp
    enddo
  enddo

  do i=1,npols
    do J=1,nant
      npol_ant(I,J) = 0.0_fp
    enddo
  enddo


  do i=1,nth
    thgrid(i) = 0.0_fp
  enddo
  if(.not.allocated(ntor)) allocate(ntor(ntors))
  if(.not.allocated(spec)) allocate(spec(ntors))
  
  if(.not.allocated(npol)) allocate(npol(npols))
  if(.not.allocated(scatthet)) allocate(scatthet(npols))
  if(.not.allocated(incithet)) allocate(incithet(npols))

  do i=1,ntors
    ntor(I) = 0.0_fp
    spec(I) = 0.0_fp
  enddo

  do i=1,npols
    npol(I) = 0.0_fp
    scatthet(I) = 0.0_fp
    incithet(I) = 0.0_fp
  enddo

  !if(.not.allocated(lh_out%dlnPdsK)) allocate(lh_out%dlnPdsK(nzones,nrays))
  !if(.not.allocated(lh_out%dlnPdsX)) allocate(lh_out%dlnPdsX(nzones,nrays))
  !if(.not.allocated(lh_out%ezsq)) allocate(lh_out%ezsq(nzones,nrays))

  !if(.not.allocated(lh_out%fe)) allocate(lh_out%fe(nv,npsi,2))
  !if(.not.allocated(lh_out%dfdv)) allocate(lh_out%dfdv(nv,npsi,2))
  !if(.not.allocated(lh_out%dql)) allocate(lh_out%dql(nv,npsi,2))

end subroutine lscsq_allocrays

subroutine lscsq_alloc
  use iso_c_binding, only : fp => c_double
  implicit none
  integer :: i
  
  if(.not.allocated(neary)) allocate(neary(npsi))
  if(.not.allocated(Teary)) allocate(Teary(npsi))
  if(.not.allocated(psiary)) allocate(psiary(npsi))
  if(.not.allocated(midary)) allocate(midary(npsi))
  if(.not.allocated(zbrary)) allocate(zbrary(npsi))
  if(.not.allocated(Edcary)) allocate(Edcary(npsi))
  if(.not.allocated(Lnlary)) allocate(Lnlary(npsi))
  if(.not.allocated(iVlary)) allocate(iVlary(npsi))
  if(.not.allocated(dVol)) allocate(dVol(npsi))
  if(.not.allocated(betZary)) allocate(betZary(npsi))



  if(.not.allocated(qlsm)) allocate(qlsm(nv))
 
  do i=1,NV
    qlsm(I) = 0.0_fp
  enddo

  if(.not.allocated(dql)) allocate(dql(nv,npsi,2))
  if(.not.allocated(dcoll)) allocate(dcoll(nv,npsi))
  if(.not.allocated(nucoll)) allocate(nucoll(nv,npsi))


  if(.not.allocated(vpar)) allocate(vpar(nv))
  if(.not.allocated(dvplus)) allocate(dvplus(nv))
  if(.not.allocated(vtherm)) allocate(vtherm(npsi))
  if(.not.allocated(Fenorm)) allocate(Fenorm(npsi))
  if(.not.allocated(nu0psi)) allocate(nu0psi(npsi))
  if(.not.allocated(FstFracn)) allocate(FstFracn(npsi)) 
  if(.not.allocated(FstFracE)) allocate(FstFracE(npsi)) 
  if(.not.allocated(vperpsq)) allocate(vperpsq(npsi))
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
  if(.not.allocated(ipintgrl)) allocate(ipintgrl(npsi))
  if(.not.allocated(vnormpos)) allocate(vnormpos(npsi))
  if(.not.allocated(vnormneg)) allocate(vnormneg(npsi))
  if(.not.allocated(vparmaxp)) allocate(vparmaxp(npsi))
  if(.not.allocated(vparmaxn)) allocate(vparmaxn(npsi))
  if(.not.allocated(jray)) allocate(jray(nv,npsi))
  if(.not.allocated(ugr)) allocate(ugr(nv))

  do i=1,npsi

    ugr(i) = 0.0_fp
    vparmaxn(i) = 0.0_fp
    vparmaxp(i) = 0.0_fp
    vnormneg(i) = 0.0_fp
    vnormpos(i) = 0.0_fp
    ipintgrl(i) = 0.0_fp
    vtherm(i) = 0.0_fp
    nu0psi(i) = 0.0_fp
    fstfracn(i) = 0.0_fp
    fstfrace(i) = 0.0_fp
    vperpsq(i) = 0.0_fp
    betZary(i) = 0.0_fp
    dVol(i) = 0.0_fp
    iVlary(i) = 0.0_fp
    Lnlary(i) = 0.0_fp
    Edcary(i) = 0.0_fp
    zbrary(i) = 0.0_fp
    midary(i) = 0.0_fp
    psiary(i) = 0.0_fp
    Teary(i) = 0.0_fp
    neary(i) = 0.0_fp
  enddo

  if(.not.allocated(pwrlevel)) allocate(pwrlevel(nrampup))
  if(.not.allocated(fecvgary)) allocate(fecvgary(nrampup))
  if(.not.allocated(wkv)) allocate(wkv(nv))

end subroutine lscsq_alloc


subroutine lscsq_constants

  implicit none

  !-----------------------------------------
  Real(fp), parameter :: zero = 0.0_fp
  Real(fp), parameter :: half = 0.5_fp
  Real(fp), parameter :: one  = 1.0_fp
  Real(fp), parameter :: two  = 2.0_fp
  !-----------------------------------------
  Real(fp), parameter :: pi   = atan2(0.0_fp,-1.0_fp)
  Real(fp), parameter :: twopi= PI+PI
  Real(fp), parameter :: sqpi = sqrt(pi)
  Real(fp), parameter :: me_Kg  = 9.10938215e-31_fp  ! electron mass (kg)
  Real(fp), parameter :: me_g   = 9.10938215e-28_fp  ! electron mass (g)
  Real(fp), parameter :: mp_Kg  = 1.672621637e-27_fp ! proton mass (kg)
  Real(fp), parameter :: mp_g   = 1.672621637e-24_fp ! proton mass (g)
  Real(fp), parameter :: qe_eV  = 1.60217653e-19_fp  !  qe (eV) 
  Real(fp), parameter :: vc     = 2.99792458e8_fp    ! speed of light
  Real(fp), parameter :: zcmtom = 1.0e-2_fp  ! Centimeter to meter
  Real(fp), parameter :: zmtocm = 1.0e+2_fp  ! meters to centimeters
  Real(fp), parameter :: zcm2tom2 = 1.0e-4_fp  ! cm^2 to m^2
  Real(fp), parameter :: zm2tocm2 = 1.0e+4_fp  ! m^2 to cm^2
  Real(fp), parameter :: zcm3tom3 = 1.0e-6_fp  ! cm^3 to m^3 | 1/m^3 to 1/cm^3
  Real(fp), parameter :: zm3tocm3 = 1.0e+6_fp  ! m^3 to cm^3 | 1/cm^3 to 1/m^3
  Real(fp), parameter :: zev2kev = 1.0e-3_fp ! eV to keV
  Real(fp), parameter :: zkev2ev = 1.0e+3_fp ! keV to eV
  Real(fp), parameter :: cgs2mks = 1.0e-7_fp ! CGS TO MKS
  Real(fp), parameter :: t2gauss = 1.0e+4_fp ! Tesla to Gauss
  Real(fp), parameter :: gauss2t = 1.0e-4_fp ! Gauss to Tesla

  Real(fp), parameter :: rad2deg = 180.0_fp/pi
  Real(fp), parameter :: deg2rad = pi/180.0_fp

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





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
  real(fp), parameter :: aee_amp = 9.578833391e+7_fp  ! electron_charge/proton_mass (C*kg^-1)
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

  integer :: nant  = 1
  integer :: ngrps = 3
  integer :: nstep = 2000 ! max steps in following each ray (500) 
  integer :: npsi = 40
  integer :: nzones = 2000
  integer :: nv = 199
  integer :: nth= 19
  integer :: nsmoo=9
  integer :: nsmw=3
  integer :: nrampup= 100 ! number of steps to ramp up power 
  integer :: nflat = 10 
                          ! 0 makes a spectrum out of arbitrary Gaussians
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
  real(fp) :: hstplh = 0.0_fp
  real(fp) :: weghtitr=0.2_fp
  real(fp) :: thet0   = 0.0_fp  ! angle of launch 0=> outside midplane, .25=> top
  real(fp) :: dthet = 0.0_fp
  real(fp), dimension(:), allocatable :: thgrid

  real(fp) :: diffujrf = 0.0_fp
  real(fp) :: prfspred = 0.0_fp
  real(fp):: tailteps = 0.0_fp
  real(fp):: tailpeps = 0.0_fp
  real(fp):: tailneps = 0.0_fp

  real(fp) :: praysum = 0.0_fp
  real(fp) :: pqlsum  = 0.0_fp
  real(fp) :: ppwrsum = 0.0_fp

  real(fp) :: TotPwr = 0.0_fp
  real(fp) :: enpar= 0.0_fp   ! n_{\parallel} launched for ray being worked at the moment
  real(fp) :: enpol=0.0_fp   ! n_{poloidal}  launched for all rays(zero is good enough)
  real(fp) :: enth =0.0_fp
  real(fp) :: omega=0.0_fp   ! RF frequency (radians/sec)

  real(fp) :: begin=0.0_fp    ! value of path length to begin ray (0 at start)
  real(fp) :: ipsq=0.0_fp     ! \sum_i \omega_{pi}^2 / \omega^2
  real(fp) :: ecyc=0.0_fp    ! \omega_{ce} / \omega
  real(fp) :: epsq=0.0_fp    ! \omega_{pe}^2 / \omega^2
  real(fp) :: ecyc2=0.0_fp   ! ecyc^2
  real(fp) :: wcei2=0.0_fp   ! \omega_{ce} \omega_{cH} \sum_i ( n_i Z_i^2 / n_e) m_H/m_i
  real(fp) :: woc2=0.0_fp    ! \omega^2/c^2
  real(fp) :: woc4=0.0_fp    ! woc2^2

  integer :: lfast = 0
  integer :: lstop = 0
  real(fp) :: delpsi=0.0_fp
  real(fp) :: Te=0.0_fp
  real(fp) :: Ne=0.0_fp
  real(fp) :: pe2min = 0.0_fp

! from lscsq_febins
  integer :: ivZero = 1
  integer :: iITR = 1

  integer :: ismo=1

  real(fp) :: vmin = -1.0_fp
  real(fp) :: vmax = +1.0_fp

  real(fp) :: TailVtrn=0.0_fp

  real(fp), dimension(neqs) ::  f1=0.0_fp
  real(fp), dimension(neqs) ::  f2=0.0_fp
  real(fp), dimension(neqs) ::  f3=0.0_fp
  real(fp), dimension(neqs) ::  y1=0.0_fp
  real(fp), dimension(neqs) ::  y2=0.0_fp
  real(fp), dimension(neqs) ::  y3=0.0_fp
  real(fp), dimension(neqsp1) :: f=0.0_fp 
  real(fp), dimension(neqsp1) :: y=0.0_fp 
 
  real(fp) :: d1=0.0_fp
  real(fp) :: d2=0.0_fp
  real(fp) :: d4=0.0_fp
  real(fp) :: denom=0.0_fp
  real(fp) :: DKpar=0.0_fp
  real(fp) :: DKper=0.0_fp
  real(fp) :: wdDdw=0.0_fp
  real(fp) :: dDdkABS=0.0_fp
  real(fp) :: epQL=0.0_fp
  real(fp) :: epsL=0.0_fp
  real(fp) :: epsz=0.0_fp
 
  real(fp) :: Eper=0.0_fp
  real(fp) :: Epar=0.0_fp
  real(fp) :: Exy=0.0_fp
  real(fp) :: Aion=0.0_fp
  real(fp) :: Aelc=0.0_fp
  real(fp) :: OmEC2=0.0_fp
  real(fp) :: EparI=0.0_fp
  real(fp) :: EparIK=0.0_fp
 
  real(fp) :: D11er=0.0_fp
  real(fp) :: D33er=0.0_fp
  real(fp) :: D12er=0.0_fp
  real(fp) :: D11ar=0.0_fp
  real(fp) :: D33ar=0.0_fp
  real(fp) :: D12ar=0.0_fp
  real(fp) :: D11w0=0.0_fp
  real(fp) :: D33w0=0.0_fp
  real(fp) :: D12w0=0.0_fp

  integer :: DqlBox(4)= 0
  integer :: nsmsym = 0
  real(fp) :: DcollNorm=0.0_fp
  real(fp) :: nuNorm=0.0_fp
  real(fp) :: DqlHite=0.0_fp
  real(fp) :: Pwrnorm=0.0_fp
  real(fp) :: DqlNorm=0.0_fp

  real(fp), dimension(:), allocatable :: powtsc, curtsc, dlJdlE,dJdE
  real(fp), dimension(:,:), allocatable :: sleave, tleave, senter, tenter

  real(fp), dimension(:,:), allocatable :: Dcoll, nucoll
  real(fp), dimension(:), allocatable :: qlsm
  real(fp), dimension(:,:,:), allocatable :: Dql

  real(fp), dimension(:), allocatable :: fghz, powers, centers, widths
  real(fp), dimension(:), allocatable :: fghz_ant, power_inp
  real(fp), dimension(:,:), allocatable :: powers_ant, centers_ant, widths_ant
  real(fp), dimension(:), allocatable :: pwrlevel, FeCvgAry
  real(fp), dimension(:), allocatable :: printgrl, pqintgrl
  real(fp), dimension(:), allocatable :: praytot, pqltot
  real(fp), dimension(:,:), allocatable :: pray, pql

  real(fp), dimension(:), allocatable :: wkzr, wkv, wkpsi, wkzx

  real(fp), dimension(:), allocatable :: ntor, Spec
  real(fp), dimension(:,:), allocatable :: ntor_ant, Spec_ant, npol_ant
  real(fp), dimension(:), allocatable :: npol  

  real(fp), dimension(:), allocatable :: vpar, vtherm, fenorm, VperpSq, &
                                         dvplus,   &
                                         nu0psi, FstFracN, FstFracE

  real(fp), dimension(:,:,:), allocatable :: fe, dfdv


  real(fp), allocatable, dimension(:) :: PsiAry,  &
                                  iVlAry, dVol, EdcAry,  &
                                  MidAry 

  real(fp), allocatable, dimension(:) :: js, jsp, nRunDot, jRunDot, ugr, &
                                         vRunIdx, IrIntgrl, IpIntgrl,    &
                                         vnormPos, vnormNeg, VparMaxP, VparMaxN

  real(fp), allocatable, dimension(:,:) :: Jray 
  
  integer, dimension(:), allocatable :: vnormOK
  integer, dimension(:), allocatable :: npeaks

  real(fp) :: vnorm = 0.0_fp
  real(fp) :: nuRuna = 0.0_fp
  real(fp) :: gmrun = 0.0_fp 
  real(fp) :: muminus = 0.0_fp
  real(fp) :: muplus = 0.0_fp
  real(fp) :: dEdcAmnt=1.0e-4_fp
  real(fp) :: vnmax = 0.99_fp 
  integer :: ivrun = 0

! from lscsq_raybins
  integer :: ips=0
  integer :: iray=0
  integer :: izone=0   
 
  integer, dimension(:,:), allocatable :: izind, ivind ! , izcind             
  real(fp), dimension(:,:), allocatable :: rofray, zofray, Pofray,        &
                                         nparry, nperry, rtpsry, powrry,&
                                         TimeRy, DistRy, DetrRy,        &
                                         NeofRy, BthRay, BphRay, rzind

  real(fp), dimension(:,:), allocatable :: ezsq, npar, power, rFudgDmp, &
                                           dlnPds, dlnPdsK, dlnPdsX    

  integer :: iRayTrsi=1
  integer :: iError=0
  integer :: iEndRy=0

  integer :: nx = 125
  integer :: nz = 159
  integer :: npsij 

  real(fp), allocatable, dimension(:) :: rho, voltlp, psivec

  integer(c_int), dimension(:), allocatable :: nz_ind, ok_ray

  real(fp), allocatable, dimension(:) :: pi2Vec, AioVec, AelVec,           &
                          pe2Vec, EdcVec, TeeVec
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
       real(fp), allocatable, dimension(:,:) :: powers
       real(fp), allocatable, dimension(:,:) :: centers
       real(fp), allocatable, dimension(:,:) :: widths
       real(fp), allocatable, dimension(:,:) :: fghz  
  end type lscsq_set


  type(lscsq_set), allocatable, dimension(:) :: lscsq_par

  type lh_rays       
    integer :: npsi
    integer :: nv
    integer :: nrays
    real(fp), dimension(:), allocatable :: psi
    real(fp), dimension(:), allocatable :: ne
    real(fp), dimension(:), allocatable :: Te
    real(fp), dimension(:), allocatable :: Edc
    real(fp), dimension(:), allocatable :: zbar
    real(fp), dimension(:), allocatable :: logL
    real(fp), dimension(:), allocatable :: betZ
    integer, dimension(:,:), allocatable :: izind
    integer, dimension(:,:), allocatable :: ivind
    integer, dimension(:), allocatable :: ok_ray
    real(fp), dimension(:), allocatable :: vpar
    real(fp), dimension(:), allocatable :: freq_ray    
    real(fp), dimension(:,:), allocatable :: dlnPdsK 
    real(fp), dimension(:,:), allocatable :: dlnPdsX
    real(fp), dimension(:,:), allocatable :: ezsq
    real(fp), dimension(:,:,:), allocatable :: fe   
    real(fp), dimension(:,:,:), allocatable :: dfdv 
    real(fp), dimension(:,:,:), allocatable :: dql  
  end type lh_rays   

  ! this type stores the plasma parameters from the cur_state file 
  type lh_plasma
    real(fp) :: mass
    real(fp) :: chrg
    real(fp) :: Rmin
    real(fp) :: Rmax
    real(fp) :: Zmin
    real(fp) :: Zmax
    real(fp) :: dx_grid
    real(fp) :: dz_grid
    real(fp) :: B_axis  
    real(fp) :: Raxis  
    real(fp) :: Zaxis  
    real(fp) :: Rminlcfs
    real(fp) :: Rmaxlcfs
    real(fp) :: Zminlcfs
    real(fp) :: Zmaxlcfs
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
    real(fp), dimension(:), allocatable :: rgrid 
    real(fp), dimension(:), allocatable :: zgrid 
    real(fp), dimension(:,:), allocatable :: Bphirz 
    real(fp), dimension(:,:), allocatable :: BrRZ   
    real(fp), dimension(:,:), allocatable :: BzRZ  
  end type lh_plasma

  type lh_param
     real(fp) :: nu0
     real(fp) :: fe0
     real(fp) :: pe2fac
     real(fp) :: pe2fac14
     real(fp) :: pi2fac
     real(fp) :: aiofac
     real(fp) :: aelfac
     real(fp) :: omcfac
     real(fp) :: ceifac
     real(fp) :: cEparIK
     real(fp) :: dqlnorm
     real(fp) :: vthnorm
     real(fp) :: pwrnorm
     real(fp) :: dcollnorm
     real(fp) :: nunorm
  end type lh_param

  type lh_map2D
     real(fp), dimension(:,:,:), allocatable :: Br
     real(fp), dimension(:,:,:), allocatable :: Bz
     real(fp), dimension(:,:,:), allocatable :: Bp
     real(fp), dimension(:,:,:), allocatable :: psi
     real(fp), dimension(:,:), allocatable :: rho
     real(fp), dimension(:,:), allocatable :: pe2
     real(fp), dimension(:,:), allocatable :: pi2
     real(fp), dimension(:,:), allocatable :: ael
     real(fp), dimension(:,:), allocatable :: aio
     real(fp), dimension(:,:), allocatable :: ne 
     real(fp), dimension(:,:), allocatable :: ni 
     real(fp), dimension(:,:), allocatable :: Te 
     real(fp), dimension(:,:), allocatable :: Ti 
  end type lh_map2D

  type(lh_map2D)  :: lh_coeff
  type(lh_plasma) :: lh_inp
  type(lh_rays)   :: lh_out
  type(lh_param)  :: lh_const


contains

subroutine lscsq_allocrays

  if(.not.allocated(senter)) allocate(senter(nzones+1,nrays))
  senter = 0.0_fp
  if(.not.allocated(tenter)) allocate(tenter(nzones+1,nrays))
  tenter = 0.0_fp
  if(.not.allocated(sleave)) allocate(sleave(nzones+1,nrays))
  sleave = 0.0_fp
  if(.not.allocated(tleave)) allocate(tleave(nzones+1,nrays))
  tleave = 0.0_fp
  if (.not.allocated(nz_ind)) allocate(nz_ind(nrays))
  nz_ind(1:nrays) = 1
  if (.not.allocated(ok_ray)) allocate(ok_ray(nrays))
  ok_ray(1:nrays) = 1

  if (.not.allocated(thgrid)) allocate(thgrid(nth))
  thgrid(1:nth) = 0.0_fp
  if (.not.allocated(Rofray)) allocate(RofRay(nzones,nrays))
  Rofray = 0.0_fp
  if (.not.allocated(Zofray)) allocate(ZofRay(nzones,nrays))
  Zofray = 0.0_fp
  if (.not.allocated(Pofray)) allocate(PofRay(nzones,nrays))
  Pofray = 0.0_fp
  if (.not.allocated(nparry)) allocate(nparry(nzones,nrays))
  nparry = 0.0_fp
  if (.not.allocated(nperry)) allocate(nperry(nzones,nrays))
  nperry = 0.0_fp
  if (.not.allocated(rtpsry)) allocate(rtpsry(nzones,nrays))
  rtpsry = 0.0_fp
  if (.not.allocated(powrry)) allocate(powrry(nzones,nrays))
  powrry = 0.0_fp
  if (.not.allocated(timery)) allocate(timery(nzones,nrays))
  timery = 0.0_fp
  if (.not.allocated(distry)) allocate(distry(nzones,nrays))
  distry = 0.0_fp
  if (.not.allocated(detrry)) allocate(detrry(nzones,nrays))
  detrry = 0.0_fp
  if (.not.allocated(neofry)) allocate(neofry(nzones,nrays))
  neofry = 0.0_fp
  if (.not.allocated(bthray)) allocate(bthray(nzones,nrays))
  bthray = 0.0_fp
  if (.not.allocated(bphray)) allocate(bphray(nzones,nrays))
  bphray = 0.0_fp
  if (.not.allocated(izind)) allocate(izind(nzones,nrays))
  izind = 0
  if (.not.allocated(rzind)) allocate(rzind(nzones,nrays))
  rzind(1:nzones,1:nrays) = 0.0_fp
  if (.not.allocated(ivind)) allocate(ivind(nzones,nrays))
  ivind = 0
  if (.not.allocated(ezsq)) allocate(ezsq(nzones,nrays))
  ezsq = 0.0_fp
  if (.not.allocated(npar)) allocate(npar(nzones,nrays))
  npar = 0.0_fp
  if (.not.allocated(power)) allocate(power(nzones,nrays))
  power = 0.0_fp
  if (.not.allocated(dlnPds)) allocate(dlnPds(nzones,nrays))
  dlnPds(1:nzones,1:nrays) = 0.0_fp
  if (.not.allocated(dlnPdsk)) allocate(dlnPdsk(nzones,nrays))
!  dlnPdsk(1:nzones,1:nrays) = 0.0_fp
  if (.not.allocated(dlnPdsx)) allocate(dlnPdsx(nzones,nrays))
!  dlnPdsx(1:nzones,1:nrays) = 0.0_fp
  if (.not.allocated(rfudgdmp)) allocate(rfudgdmp(nzones,nrays))
  rfudgdmp(1:nzones,1:nrays) = 0.0_fp


  if(.not.allocated(ntor_ant)) allocate(ntor_ant(ntors,nant))
  ntor_ant(1:ntors,1:nant) = 0.0_fp
  if(.not.allocated(spec_ant)) allocate(spec_ant(ntors,nant))
  spec_ant(1:ntors,1:nant) = 0.0_fp
  if(.not.allocated(npol_ant)) allocate(npol_ant(npols,nant))
  npol_ant(1:npols,1:nant) = 0.0_fp
  if(.not.allocated(ntor)) allocate(ntor(ntors))
  ntor(1:ntors) = 0.0_fp
  if(.not.allocated(spec)) allocate(spec(ntors))
  spec(1:ntors) = 0.0_fp
  if(.not.allocated(npol)) allocate(npol(npols))
  npol(1:npols) = 0.0_fp


  if(.not.allocated(lh_out%psi)) allocate(lh_out%psi(npsi))
  if(.not.allocated(lh_out%ne)) allocate(lh_out%ne(npsi))
  if(.not.allocated(lh_out%Te)) allocate(lh_out%Te(npsi))
  if(.not.allocated(lh_out%zbar)) allocate(lh_out%zbar(npsi))
  if(.not.allocated(lh_out%logL)) allocate(lh_out%logL(npsi))
  if(.not.allocated(lh_out%betZ)) allocate(lh_out%betZ(npsi))
  if(.not.allocated(lh_out%dlnPdsK)) allocate(lh_out%dlnPdsK(nzones,nrays))
  if(.not.allocated(lh_out%dlnPdsX)) allocate(lh_out%dlnPdsX(nzones,nrays))
  if(.not.allocated(lh_out%ezsq)) allocate(lh_out%ezsq(nzones,nrays))
  if(.not.allocated(lh_out%izind)) allocate(lh_out%izind(nzones,nrays))
  if(.not.allocated(lh_out%ivind)) allocate(lh_out%ivind(nzones,nrays))
  if(.not.allocated(lh_out%ok_ray)) allocate(lh_out%ok_ray(nrays))

  if(.not.allocated(lh_out%fe)) allocate(lh_out%fe(nv,npsi,2))
  if(.not.allocated(lh_out%dfdv)) allocate(lh_out%dfdv(nv,npsi,2))
  if(.not.allocated(lh_out%dql)) allocate(lh_out%dql(nv,npsi,2))
  if(.not.allocated(lh_out%vpar)) allocate(lh_out%vpar(nv))

end subroutine lscsq_allocrays

subroutine lscsq_alloc

  if(.not.allocated(psiary)) allocate(psiary(npsi))
  psiary(1:npsi) = 0.0_fp
  if(.not.allocated(midary)) allocate(midary(npsi))
  midary(1:npsi) = 0.0_fp
  if(.not.allocated(Edcary)) allocate(Edcary(npsi))
  Edcary(1:npsi) = 0.0_fp
  if(.not.allocated(iVlary)) allocate(iVlary(npsi))
  iVlary(1:npsi) = 0.0_fp
  if(.not.allocated(dVol)) allocate(dVol(npsi))
  dVol(1:npsi) = 0.0_fp

  if(.not.allocated(qlsm)) allocate(qlsm(nv))
  qlsm(1:nv) = 0.0_fp
  if(.not.allocated(dql)) allocate(dql(nv,npsi,2))
  if(.not.allocated(dcoll)) allocate(dcoll(nv,npsi))
  if(.not.allocated(nucoll)) allocate(nucoll(nv,npsi))


  if(.not.allocated(vpar)) allocate(vpar(nv))
  vpar(1:nv) = 0.0_fp  
  if(.not.allocated(dvplus)) allocate(dvplus(nv))
  if(.not.allocated(vtherm)) allocate(vtherm(npsi))
  vtherm(1:npsi) = 0.0_fp
  if(.not.allocated(Fenorm)) allocate(Fenorm(npsi))
  if(.not.allocated(nu0psi)) allocate(nu0psi(npsi))
  nu0psi(1:npsi) = 0.0_fp  
  if(.not.allocated(FstFracn)) allocate(FstFracn(npsi))
  fstfracn(1:npsi) = 0.0_fp  
  if(.not.allocated(FstFracE)) allocate(FstFracE(npsi))
  fstfrace(1:npsi) = 0.0_fp  
  if(.not.allocated(vperpsq)) allocate(vperpsq(npsi))
  vperpsq(1:npsi) = 0.0_fp  
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
  irintgrl(1:npsi) = 0.0_fp 
  if(.not.allocated(ipintgrl)) allocate(ipintgrl(npsi))
  ipintgrl(1:npsi) = 0.0_fp 
  if(.not.allocated(vnormpos)) allocate(vnormpos(npsi))
  vnormpos(1:npsi) = 0.0_fp 
  if(.not.allocated(vnormneg)) allocate(vnormneg(npsi))
  vnormneg(1:npsi) = 0.0_fp 
  if(.not.allocated(vparmaxp)) allocate(vparmaxp(npsi))
  vparmaxp(1:npsi) = 0.0_fp 
  if(.not.allocated(vparmaxn)) allocate(vparmaxn(npsi))
  vparmaxn(1:npsi) = 0.0_fp 
  if(.not.allocated(jray)) allocate(jray(nv,npsi))
  if(.not.allocated(ugr)) allocate(ugr(nv))
  ugr(1:nv) = 0.0_fp 

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

  if(.not.allocated(pi2vec)) allocate(pi2vec(npsij))
  if(.not.allocated(aiovec)) allocate(aiovec(npsij))
  if(.not.allocated(aelvec)) allocate(aelvec(npsij))
  if(.not.allocated(pe2vec)) allocate(pe2vec(npsij))
  if(.not.allocated(edcvec)) allocate(edcvec(npsij))
  if(.not.allocated(voltlp)) allocate(voltlp(npsij))

  if(.not.allocated(powtsc)) allocate(powtsc(npsij))
  if(.not.allocated(curtsc)) allocate(curtsc(npsij))
  if(.not.allocated(dJdE)) allocate(dJdE(npsij))
  if(.not.allocated(dlJdlE)) allocate(dlJdlE(npsij))

end subroutine alloc_profs

!--------------------------------------------------------------


subroutine dealloc_profs

implicit none

  deallocate(pi2vec,aiovec,aelvec,pe2vec)
  deallocate(edcvec,voltlp)
  deallocate(powtsc,curtsc,dJdE,dlJdlE)


end subroutine dealloc_profs



end module lscsq_mod
!-----------------------





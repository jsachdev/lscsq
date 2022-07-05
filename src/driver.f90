program driver

!  use iso_c_binding, only: fp => c_double, c_int
  use lscsq_mod, only: npsij, lh_inp, lh_coeff, lh_const
  use lscsq_mod, only: voltlp, Edcvec 
  use lscsq_mod, only: nx, nz, pe2min, wcei2, psivec
  use lscsq_mod, only: zmtocm, zkev2ev, zel, zcm3tom3, me_g, twopi
  use lscsq_mod, only: pe2vec, pi2vec, aelvec, aiovec
  use lscsq_mod, only: me_kg, mp_kg, pi, vc, qe_ev
  use lscsq_mod, only: i1stcall, ierror
  use lscsq_mod, only: alloc_profs, dealloc_profs

  use netcdf ! to read from CDF file
  use EZspline_obj      
  use EZspline
  
  implicit none

  character(len=80) :: cdffile
  !============
  integer :: ncid,vid,ierr, ier, istat
  integer :: nt, nr 
  integer :: j0, j, i, k, kk, idum, n

  real(fp) :: zchar, zhydmass
  real(fp) :: dum

  real(fp), dimension(:), allocatable :: zpowtsc, zcurtsc, ZdlJdlE, ZdJdE
  real(fp), dimension(:,:,:), allocatable :: f2d
  real(fp), dimension(:,:), allocatable :: f1d
  integer :: iraytr=1

  type(EZspline1) :: spline_o ! 1D interpolation
  integer :: sp_n1,sp_ier
  integer, dimension(2) :: sp_BCS1=0
  real(fp), dimension(:), allocatable :: sp_x1,sp_f ! function values
  real(fp) :: sp_fp  
  

  cdffile = 'cur_state.cdf'
  call lscsq_readpstate(cdffile)

  call alloc_profs

  ! all arrays with dimension set by the external code (e.g. TRANSP) 
  ! are allocated here

  voltlp = lh_inp%Vloop

  !	--------------------------------------------------------------

  ! here all operations
!     pe2Fac pi2Fac convert density in ^14 cm^-3 or ^20 m^-3
!     and frequency in GHz into

  lh_const%pe2fac = 1.0e-5_fp*(qe_eV*vc)**2/(pi*me_Kg) 
  lh_const%pe2Fac14= 1.0e-19_fp*(qe_eV*vc )**2/(PI*me_Kg)
  lh_const%pi2fac  = lh_const%pe2fac *(me_Kg/mp_Kg) 
  lh_const%AioFac  = 3.0e-15_fp * lh_const%pi2Fac * qe_eV/mp_Kg/TWOPI**2
  lh_const%AelFac  = 0.75e-15_fp * lh_const%pe2Fac * qe_eV/me_Kg/TWOPI**2 
  lh_const%OmcFac  = 1.0e-9_fp*qe_eV/me_Kg/TWOPI 
  lh_const%ceifac = 1.0e-18_fp*(qe_eV/twopi)**2/me_Kg/mp_Kg
  lh_const%cEparIK = 1.0e2_fp*twopi*twopi*qe_eV**2/me_g
  lh_const%DqlNorm = 2.5e-10_fp*(qe_eV/me_Kg)**2 / vc**2
  lh_const%VthNorm = sqrt(1.6e-12_fp*zkev2eV/me_g)/(vc*zmtocm)       
  lh_const%fe0 = 1.0_fp/sqrt(twopi)
  lh_const%nu0 = 4.0_fp*pi*zel**4/me_g**2
  lh_const%PwrNorm = 1.0e6_fp*me_Kg*vc**2
  lh_const%nuNorm =    zcm3tom3*lh_const%nu0/vc**3
  lh_const%DcollNorm = lh_const%nuNorm

  EdcVec(1:npsij) = voltlp(1:npsij)/(twopi*lh_inp%Raxis)

  pe2vec = 1.0e-20_fp*lh_const%pe2fac*lh_inp%ne
  AelVec(1:npsij) = 1.0e-20_fp*lh_const%AelFac*lh_inp%ne(1:npsij)*lh_inp%Te(1:npsij)
  ! we should sum over ion species, including impurities. However, LH does not
  ! heat on ions. Take for now only the background species (1 species for test)
  pi2Vec(1:npsij) = 1.0e-20_fp*lh_const%pi2Fac*lh_inp%ni(1:npsij)*(lh_inp%chrg/qe_eV)**2*mp_Kg/lh_inp%mass
  AioVec(1:npsij) = 1.0e-20_fp*lh_const%AioFac*lh_inp%ni(1:npsij)*(lh_inp%chrg/qe_eV)**2*(mp_Kg/lh_inp%mass)**2*lh_inp%Ti(1:npsij)

  wcei2 = lh_const%ceifac*lh_inp%ni(1)/lh_inp%ne(1)*(lh_inp%chrg/qe_eV)**2*mp_kg/lh_inp%mass
  wcei2 = wcei2*lh_inp%B_axis**2

  pe2min = pe2Vec(NpsiJ)
!  pe2min = 1.0e-20_fp*lh_const%pe2fac*lh_inp%ne(npsij)

  call nml2lsc

  ! here goes the call to the 2D pspline, to interpolate the following profiles
  ! that are used for the calculation of the dispersion relation:
  ! electron density
  ! electron and ion temperature
  ! components of the magnetic field
  ! psi(R,z)  
  ! the cyclotron frequency will be calculated locally
  ! all plasma and fundamental frequencies will be calculated locally from the
  ! plasma parameters interpolated along the ray trajectory. 
  ! the pspline coefficients need to be stored in a separate type that is passed
  ! to the ray tracing. We define the types here, they do not need to be defined
  ! as a global parameter, they can just be passed as an argument.   

  ! start with magnetic field only - it will be passed directly to avoid
  ! calculating derivatives of psi. This will eliminate the grap.f90 and
  ! simplify grapgrd.f90

  ! the spline coefficients will be stored in types, so that they can be passed
  ! as arguments => good for GPU
  allocate(lh_coeff%Bp(nx,nz,9))
  allocate(lh_coeff%Br(nx,nz,9))
  allocate(lh_coeff%Bz(nx,nz,9))
  allocate(lh_coeff%psi(nx,nz,9))
  allocate(f2d(nx,nz,9))

  f2d = 0.0_fp
  f2d(:,:,1) = lh_inp%BphiRZ
  call splnrz(f2d)
  lh_coeff%Bp(:,:,1) = lh_inp%BphiRZ
  lh_coeff%Bp(:,:,2:9) = f2d(:,:,2:9)

  f2d = 0.0_fp
  f2d(:,:,1) = lh_inp%BrRZ
  call splnrz(f2d)
  lh_coeff%Br(:,:,1) = lh_inp%BrRZ
  lh_coeff%Br(:,:,2:9) = f2d(:,:,2:9)

  f2d = 0.0_fp
  f2d(:,:,1) = lh_inp%BzRZ
  call splnrz(f2d)
  lh_coeff%Bz(:,:,1) = lh_inp%BzRZ
  lh_coeff%Bz(:,:,2:9) = f2d(:,:,2:9)

  f2d = 0.0_fp
  f2d(:,:,1) = lh_inp%psiRZ
  call splnrz(f2d)
  lh_coeff%psi(:,:,1) = lh_inp%psiRZ
  lh_coeff%psi(:,:,2:9) = f2d(:,:,2:9)
 
  ! derive profiles of ne, Te, Ti on a regular psi grid
  allocate(psivec(npsij))
  ! generate regular psi grid
  CALL lscsq_ugrid(Psivec, npsij, lh_inp%plflx(1), lh_inp%plflx(npsij)) 

  ! use EZspline to interpolate ne, Te, Ti profiles to regular grid
  allocate(sp_f(npsij))
  allocate(sp_x1(npsij))
  !     impose boundary conditions on first derivative/node (spline only)
      
  !     initialize call for spline interpolation 
  call EZspline_init(spline_o,npsij,sp_BCS1,sp_ier) ! spline
  call EZspline_error(sp_ier)
      
  !     set boundary conditions for 1st derivative
  spline_o%bcval1min = 0.0_fp
  spline_o%bcval1max = 0.0_fp
      
  ! radial basis for interpolation from transp.dat
  spline_o%x1(:) = psivec                
      
  ! ne profile      
  sp_f(:) = lh_inp%ne      
  ! setup interpolating function
  call EZspline_setup(spline_o,sp_f,sp_ier)
  call EZspline_error(sp_ier)

  allocate(lh_coeff%ne(npsij,3))
  ! interpolate ne           
  do i = 1,npsij 
     call EZspline_interp(spline_o,lh_inp%plflx(i),sp_fp,sp_ier)
     call EZspline_error(sp_ier)
     lh_coeff%ne(i,1) = sp_fp
  enddo
   
  allocate(f1d(npsij,3))
  f1d = 0.0_fp  
  f1d(:,1) = lh_coeff%ne(:,1)
  call spln1d(npsij,psivec,f1d)
  lh_coeff%ne(:,2:3) = f1d(:,2:3)

  ! ----------------------------------
  ! Te profile      
  sp_f(:) = lh_inp%Te      
  ! setup interpolating function
  call EZspline_setup(spline_o,sp_f,sp_ier)
  call EZspline_error(sp_ier)

  allocate(lh_coeff%te(npsij,3))
  do i = 1,npsij 
     call EZspline_interp(spline_o,lh_inp%plflx(i),sp_fp,sp_ier)
     call EZspline_error(sp_ier)
     lh_coeff%Te(i,1) = sp_fp
  enddo
  f1d = 0.0_fp  
  f1d(:,1) = lh_coeff%Te(:,1)
  call spln1d(npsij,psivec,f1d)
  lh_coeff%Te(:,2:3) = f1d(:,2:3)

  ! ----------------------------------
  ! Ti profile      
  sp_f(:) = lh_inp%Ti      
  ! setup interpolating function
  call EZspline_setup(spline_o,sp_f,sp_ier)
  call EZspline_error(sp_ier)

  allocate(lh_coeff%Ti(npsij,3))
  do i = 1,npsij 
     call EZspline_interp(spline_o,lh_inp%plflx(i),sp_fp,sp_ier)
     call EZspline_error(sp_ier)
     lh_coeff%Ti(i,1) = sp_fp
  enddo
  f1d = 0.0_fp  
  f1d(:,1) = lh_coeff%Ti(:,1)
  call spln1d(npsij,psivec,f1d)
  lh_coeff%Ti(:,2:3) = f1d(:,2:3)

  ! ----------------------------------
  ! ni profile      
  sp_f(:) = lh_inp%ni      
  ! setup interpolating function
  call EZspline_setup(spline_o,sp_f,sp_ier)
  call EZspline_error(sp_ier)

  allocate(lh_coeff%ni(npsij,3))
  do i = 1,npsij 
     call EZspline_interp(spline_o,lh_inp%plflx(i),sp_fp,sp_ier)
     call EZspline_error(sp_ier)
     lh_coeff%ni(i,1) = sp_fp
  enddo
  f1d = 0.0_fp  
  f1d(:,1) = lh_coeff%ni(:,1)
  call spln1d(npsij,psivec,f1d)
  lh_coeff%ni(:,2:3) = f1d(:,2:3)


  i1stcall = 1
  iraytr = 1  
  call lscsq_main(iraytr,iError)

  call lscsq_writecdf

  call dealloc_profs

contains


subroutine lscsq_readpstate(pstate_file)
  use netcdf
  use lscsq_mod 
  implicit none

  character(len=80), intent(in) :: pstate_file
  integer :: ncid, status
  integer :: rhodimid, reqdimid, dim_nrho, vid, vidimid, nsadimid, lhdim
  integer :: dm1_nrho, rh1dimid, nsdimid, dp1_nspec_th, dim_nrho_eq
  integer :: dp1_nspec_tha, dim_nr, dim_nz, rdimid, zdimid, dim_nlhrf_src
  real(fp) :: dum0
  real(fp), dimension(:), allocatable :: dums, dumsa, dum_r, dum_z, dumlh
  real(fp), dimension(:), allocatable :: dum, dumZ, rho_eq, dum1, dum1p
  real(fp), dimension(:,:), allocatable :: dum2d, dumrz

  status = nf90_open(pstate_file,NF90_NOWRITE,ncid)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_dimid(ncid, "dim_nrho", rhodimid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inquire_dimension(ncid, rhodimid, len = dim_nrho)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_dimid(ncid, "dim_nlhrf_src", lhdim)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inquire_dimension(ncid, lhdim, len = dim_nlhrf_src)
  if (status /= nf90_noerr) call handle_err(status)
 
  if(.not.allocated(dumlh)) allocate(dumlh(dim_nlhrf_src))
  if(.not.allocated(lh_inp%powerlh)) allocate(lh_inp%powerlh(dim_nlhrf_src))
  if(.not.allocated(lh_inp%freqlh))  allocate(lh_inp%freqlh(dim_nlhrf_src))
  ierr = nf90_inq_varid(ncid,'power_lh',vid)
  status = nf90_get_var(ncid,vid,dumlh)
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%powerlh = dumlh
  ierr = nf90_inq_varid(ncid,'freq_lh',vid)
  status = nf90_get_var(ncid,vid,dumlh)
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%freqlh = dumlh

  npsij = dim_nrho
  allocate(lh_inp%Ti(npsij))
  allocate(lh_inp%Te(npsij))
  allocate(lh_inp%ne(npsij))
  allocate(lh_inp%ni(npsij))
  allocate(lh_inp%Vloop(npsij))
  allocate(lh_inp%plflx(npsij))
  allocate(lh_inp%dvol(npsij))
  allocate(lh_inp%vol(npsij))
  allocate(lh_inp%g_eq(npsij))

  status = nf90_inq_dimid(ncid, "dim_nr", rdimid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inquire_dimension(ncid, rdimid, len = dim_nr)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_dimid(ncid, "dim_nz", zdimid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inquire_dimension(ncid, zdimid, len = dim_nz)
  if (status /= nf90_noerr) call handle_err(status)

  allocate(lh_inp%psirz(dim_nr,dim_nz))
  allocate(lh_inp%BphiRZ(dim_nr,dim_nz))
  allocate(lh_inp%BrRZ(dim_nr,dim_nz))
  allocate(lh_inp%BzRZ(dim_nr,dim_nz))
  allocate(lh_inp%rgrid(dim_nr))
  allocate(lh_inp%zgrid(dim_nz))
  nx = dim_nr
  nz = dim_nz

  status = nf90_inq_dimid(ncid, "dm1_nrho", rh1dimid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inquire_dimension(ncid, rh1dimid, len = dm1_nrho)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_dimid(ncid, "dp1_nspec_th", nsdimid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inquire_dimension(ncid, nsdimid, len = dp1_nspec_th)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_dimid(ncid, "dp1_nspec_th", nsadimid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inquire_dimension(ncid, nsadimid, len = dp1_nspec_tha)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_dimid(ncid, "dim_nrho_eq", reqdimid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inquire_dimension(ncid, reqdimid, len = dim_nrho_eq)
  if (status /= nf90_noerr) call handle_err(status)

  if(.not.allocated(rho_eq)) allocate(rho_eq(dim_nrho_eq))
  ierr = nf90_inq_varid(ncid,'rho_eq',vid)
  ierr = nf90_get_var(ncid,vid,rho_eq) 

  if(.not.allocated(rho)) allocate(rho(dim_nrho))
  ierr = nf90_inq_varid(ncid,'rho',vid)
  ierr = nf90_get_var(ncid,vid,rho) 

  if(.not.allocated(dumz)) allocate(dumz(dp1_nspec_th))
  if(.not.allocated(dums)) allocate(dums(dp1_nspec_th))
  if(.not.allocated(dumsa)) allocate(dumsa(dp1_nspec_tha))
  if(.not.allocated(dum1)) allocate(dum1(dm1_nrho))
  if(.not.allocated(dum1p)) allocate(dum1p(dim_nrho))
  if(.not.allocated(dum2d)) allocate(dum2d(dm1_nrho,dp1_nspec_th))

!        double vol(dim_nrho_eq) ;
  ierr = nf90_inq_varid(ncid,'vol',vid)
  status = nf90_get_var(ncid,vid,dum1p) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%vol = dum1p
  lh_inp%dvol(1) = 0.0_fp
  do i=2,dim_nrho_eq
     lh_inp%dvol(i) = lh_inp%vol(i)-lh_inp%vol(i-1) 
  enddo

  allocate(dumrz(dim_nr,dim_nz))
  ierr = nf90_inq_varid(ncid,'PsiRZ',vid)
  status = nf90_get_var(ncid,vid,dumrz) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%psirz = dumrz
  
  ierr = nf90_inq_varid(ncid,'BphiRZ',vid)
  status = nf90_get_var(ncid,vid,dumrz) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%BphiRZ = dumrz

  ierr = nf90_inq_varid(ncid,'BRRZ',vid)
  status = nf90_get_var(ncid,vid,dumrz) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%BrRZ = dumrz

  ierr = nf90_inq_varid(ncid,'BZRZ',vid)
  status = nf90_get_var(ncid,vid,dumrz) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%BzRZ = dumrz

!        double psipol(dim_nrho_eq) ;
  ierr = nf90_inq_varid(ncid,'psipol',vid)
  status = nf90_get_var(ncid,vid,dum1p) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%plflx = dum1p
  
  ierr = nf90_inq_varid(ncid,'Ti_bdy',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Ti(npsij) = dum0 

  ierr = nf90_inq_varid(ncid,'Te_bdy',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Te(npsij) = dum0 

  ierr = nf90_inq_varid(ncid,'Ts',vid)
  status = nf90_get_var(ncid,vid,dum2d) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Te(1) = dum2d(1,1)+0.5_fp*(dum2d(1,1)-dum2d(2,1)) 
  do i=1,dm1_nrho-1
     lh_inp%Te(i+1) = 0.5_fp*(dum2d(i,1)+dum2d(i+1,1)) 
  enddo 

  ierr = nf90_inq_varid(ncid,'Ti',vid)
  status = nf90_get_var(ncid,vid,dum1) 
  if (status /= nf90_noerr) call handle_err(status)
!  from zone center (?) to zone boundary
  lh_inp%Ti(1) = dum1(1)+0.5_fp*(dum1(1)-dum1(2)) 
  do i=1,dm1_nrho-1
     lh_inp%Ti(i+1) = 0.5_fp*(dum1(i)+dum1(i+1)) 
  enddo 

  ierr = nf90_inq_varid(ncid,'ns_bdy',vid)
  status = nf90_get_var(ncid,vid,dumZ) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%ne(npsij) = dumz(1)
  lh_inp%ni(npsij) = sum(dumz(2:dp1_nspec_th))

  ierr = nf90_inq_varid(ncid,'ns',vid)
  status = nf90_get_var(ncid,vid,dum2d) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%ne(1) = dum2d(1,1)+0.5_fp*(dum2d(1,1)-dum2d(2,1)) 
  do i=1,dm1_nrho-1
     lh_inp%ne(i+1) = 0.5_fp*(dum2d(i,1)+dum2d(i+1,1)) 
  enddo 

  ierr = nf90_inq_varid(ncid,'ni',vid)
  status = nf90_get_var(ncid,vid,dum1) 
  if (status /= nf90_noerr) call handle_err(status)
!  from zone center (?) to zone boundary
  lh_inp%ni(1) = dum1(1)+0.5_fp*(dum1(1)-dum1(2)) 
  do i=1,dm1_nrho-1
     lh_inp%ni(i+1) = 0.5_fp*(dum1(i)+dum1(i+1)) 
  enddo 
  ierr = nf90_inq_varid(ncid,'ni_bdy',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)


!        double V_loop(dim_nrho) ;
  ierr = nf90_inq_varid(ncid,'V_loop',vid)
  status = nf90_get_var(ncid,vid,dum1p) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Vloop = dum1p

  
  ierr = nf90_inq_varid(ncid,'q_S',vid)
  status = nf90_get_var(ncid,vid,dums) 
  if (status /= nf90_noerr) call handle_err(status)
!  write(*,*) 'q_S:',dums/1.602e-19
  lh_inp%chrg = dums(2)  ! take only D

  ierr = nf90_inq_varid(ncid,'q_Sa',vid)
  status = nf90_get_var(ncid,vid,dumsa) 
  if (status /= nf90_noerr) call handle_err(status)

  ierr = nf90_inq_varid(ncid,'m_S',vid)
  status = nf90_get_var(ncid,vid,dums) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%mass = dums(2)  ! take only D

  ierr = nf90_inq_varid(ncid,'R_min_box',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Rmin = dum0

  ierr = nf90_inq_varid(ncid,'R_max_box',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Rmax = dum0

  ierr = nf90_inq_varid(ncid,'R_min_lcfs',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Rminlcfs = dum0

  ierr = nf90_inq_varid(ncid,'R_max_lcfs',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Rmaxlcfs = dum0

  
  ierr = nf90_inq_varid(ncid,'Z_min_box',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Zmin = dum0

  ierr = nf90_inq_varid(ncid,'Z_max_box',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Zmax = dum0


  ierr = nf90_inq_varid(ncid,'Z_min_lcfs',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Zminlcfs = dum0

  ierr = nf90_inq_varid(ncid,'Z_max_lcfs',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Zmaxlcfs = dum0

  ierr = nf90_inq_varid(ncid,'g_eq',vid)
  status = nf90_get_var(ncid,vid,dum1p) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%g_eq = dum1p

  ierr = nf90_inq_varid(ncid,'R_axis',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%Raxis = dum0

  ierr = nf90_inq_varid(ncid,'Z_axis',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%zaxis = dum0

  ierr = nf90_inq_varid(ncid,'B_axis',vid)
  status = nf90_get_var(ncid,vid,dum0) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%B_axis = dum0
  
  allocate(dum_z(dim_nz)) 
  ierr = nf90_inq_varid(ncid,'Z_grid',vid)
  status = nf90_get_var(ncid,vid,dum_z) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%dz_grid=(dum_z(dim_nz)-dum_z(1))/(dim_nz-1)
  lh_inp%zgrid = dum_z 

  allocate(dum_r(dim_nr)) 
  ierr = nf90_inq_varid(ncid,'R_grid',vid)
  status = nf90_get_var(ncid,vid,dum_r) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%dx_grid = (dum_r(dim_nr)-dum_r(1))/(dim_nr-1)
  lh_inp%rgrid = dum_r 


  close(ncid)


end subroutine lscsq_readpstate



subroutine lscsq_writecdf

  use netcdf
  use lscsq_mod
  implicit none

  character(len=20) :: cdf_out = 'lsc_out.cdf'
  character(len=*), parameter :: UNITS = "units"

  integer :: xdimid, rdimid,xxdimid,zdimid, nzonid, nvid, nrayid
  integer :: sigdimid,nbdryid
  integer, dimension(3) :: dimindf
  integer, dimension(2) :: dimindr, dimvx
  integer :: ierr, ncid, i
  integer :: vid_psi, vid_r, vid_z, vid_zon, vid_vel, vid_ray
  integer :: vid_rray, vid_zray, vid_npar, vid_nper
  integer :: vid_pwry, vid_tmry, vid_dsry, vid_detry, vid_nery, vid_btry
  integer :: vid_rtpy, vid_pray,vid_edcin,vid_edcou,vid_near,vid_tear
  integer :: vid_psiou, vid_midou 
  integer :: vid_js, vid_jsp, vid_fe, vid_dfdv
  integer :: vid_vth, vid_neou, vid_teou, vid_nuc, vid_dc, vid_dql
  integer :: vid_pow, vid_pql, vid_nzon
  integer :: nx_tr, nxdim, vid_x, vid_powtr, vid_curtr
  integer :: vid_bdry, vid_rbdry, vid_zbdry



  ! -----------------------------------------------------
  ! Create netCDF file
  ierr=nf90_create(trim(cdf_out),nf90_clobber,ncid)
  ! --------------------------------------
  ! dimensions 
  ierr=nf90_def_dim(ncid,'npsi',npsi,xdimid)
  ierr=nf90_def_dim(ncid,'nr',nx,rdimid)
  ierr=nf90_def_dim(ncid,'nz',nz,zdimid)
  ierr=nf90_def_dim(ncid,'nzon',nzones,nzonid)
  ierr=nf90_def_dim(ncid,'nv',nv,nvid)
  ierr=nf90_def_dim(ncid,'nrays',nrays,nrayid)
  ierr=nf90_def_dim(ncid,'sigma',2,sigdimid) 

  dimindf=(/ nvid,xdimid,sigdimid /)
  dimindr=(/ nzonid,nrayid /)
  dimvx = (/ nvid,xdimid /)
  ! --------------------------------------
  ! definitions
  ierr=nf90_def_var(ncid,'psi',nf90_real,xdimid,vid_psi)
  ierr=nf90_def_var(ncid,'R',nf90_real,rdimid,vid_r)
  ierr=nf90_def_var(ncid,'Z',nf90_real,zdimid,vid_z)
  ierr=nf90_def_var(ncid,'vpar',nf90_real,nvid,vid_vel)
  
  ierr=nf90_def_var(ncid,'Rbdry',nf90_real,nbdryid,vid_rbdry)
  ierr=nf90_def_var(ncid,'Zbdry',nf90_real,nbdryid,vid_zbdry)
  ierr=nf90_def_var(ncid,'Rofray',nf90_real,dimindr,vid_rray)
  ierr=nf90_def_var(ncid,'Zofray',nf90_real,dimindr,vid_zray)
  ierr=nf90_def_var(ncid,'Pofray',nf90_real,dimindr,vid_pray)
  ierr=nf90_def_var(ncid,'Pwrray',nf90_real,dimindr,vid_pwry)
  ierr=nf90_def_var(ncid,'timray',nf90_real,dimindr,vid_tmry)
  ierr=nf90_def_var(ncid,'disray',nf90_real,dimindr,vid_dsry)
  ierr=nf90_def_var(ncid,'detray',nf90_real,dimindr,vid_detry)
  ierr=nf90_def_var(ncid,'neray',nf90_real,dimindr,vid_nery)
  ierr=nf90_def_var(ncid,'Btray',nf90_real,dimindr,vid_Btry)
  ierr=nf90_def_var(ncid,'rtpsray',nf90_real,dimindr,vid_rtpy)
  ierr=nf90_def_var(ncid,'npar',nf90_real,dimindr,vid_npar)
  ierr=nf90_def_var(ncid,'nper',nf90_real,dimindr,vid_nper)
  ierr=nf90_def_var(ncid,'Te_ou',nf90_real,xdimid,vid_Teou)
  ierr=nf90_def_var(ncid,'ne_ou',nf90_real,xdimid,vid_neou)
  ierr=nf90_def_var(ncid,'Edc_ou',nf90_real,xdimid,vid_edcou)
  ierr=nf90_def_var(ncid,'psi_ou',nf90_real,xdimid,vid_psiou)
  ierr=nf90_def_var(ncid,'mid_ou',nf90_real,xdimid,vid_midou)
  ierr=nf90_def_var(ncid,'Js',nf90_real,xdimid,vid_js)
  ierr=nf90_def_var(ncid,'Jsp',nf90_real,xdimid,vid_jsp)
  ierr=nf90_def_var(ncid,'Power',nf90_real,xdimid,vid_pow)
  ierr=nf90_def_var(ncid,'Power_ql',nf90_real,xdimid,vid_pql)
  ierr=nf90_def_var(ncid,'fe',nf90_real,dimindf,vid_fe)
  ierr=nf90_def_var(ncid,'Dql',nf90_real,dimindf,vid_dql)
  ierr=nf90_def_var(ncid,'vth',nf90_real,xdimid,vid_vth)
  ierr=nf90_def_var(ncid,'Dcoll',nf90_real,dimvx,vid_dc)
  ierr=nf90_def_var(ncid,'nucoll',nf90_real,dimvx,vid_nuc)
  ierr=nf90_def_var(ncid,'nzon',nf90_real,nrayid,vid_nzon)


  ierr=nf90_enddef(ncid)
 
  ! --------------------------------------
  ! Write variables to CDF file
  ierr=nf90_put_var(ncid,vid_psi,midary)
  if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr=nf90_put_var(ncid,vid_rray,rofray)
  ierr=nf90_put_var(ncid,vid_zray,zofray)
  ierr=nf90_put_var(ncid,vid_npar,nparry)
  ierr=nf90_put_var(ncid,vid_nper,nperry)
  ierr=nf90_put_var(ncid,vid_nper,nperry)
  ierr=nf90_put_var(ncid,vid_pray,Pofray)
  ierr=nf90_put_var(ncid,vid_pwry,Powrry)
  ierr=nf90_put_var(ncid,vid_tmry,Timery)
  ierr=nf90_put_var(ncid,vid_dsry,distry)
  ierr=nf90_put_var(ncid,vid_detry,detrry)
  ierr=nf90_put_var(ncid,vid_nery,Neofry)
  ierr=nf90_put_var(ncid,vid_Btry,Bthray)
  ierr=nf90_put_var(ncid,vid_rtpy,rtpsry)
  
  ierr=nf90_put_var(ncid,vid_neou,lh_out%ne)
  ierr=nf90_put_var(ncid,vid_teou,lh_out%te)
  ierr=nf90_put_var(ncid,vid_psiou,psiary)
  ierr=nf90_put_var(ncid,vid_Edcou,Edcary)
  ierr=nf90_put_var(ncid,vid_Edcin,Edcvec)
  ierr=nf90_put_var(ncid,vid_midou,midary)

  ierr=nf90_put_var(ncid,vid_js,js)
  ierr=nf90_put_var(ncid,vid_jsp,jsp)
  ierr=nf90_put_var(ncid,vid_pow,praytot)
  ierr=nf90_put_var(ncid,vid_pql,pqltot)
  ierr=nf90_put_var(ncid,vid_fe,fe)
  ierr=nf90_put_var(ncid,vid_vth,vtherm)
  ierr=nf90_put_var(ncid,vid_vel,vpar)
  ierr=nf90_put_var(ncid,vid_dc,dcoll)
  ierr=nf90_put_var(ncid,vid_nuc,nucoll)
  ierr=nf90_put_var(ncid,vid_dql,Dql)
  ierr=nf90_put_var(ncid,vid_nzon,nz_ind)

  ierr=nf90_close(ncid)

end subroutine lscsq_writecdf

subroutine handle_err(status)
      
  USE netcdf ! to read from CDF file
  IMPLICIT NONE 
      
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then
     write(0,*) trim(nf90_strerror(status))
     write(0,*) ''
     write(0,*) '*** Program stopped.'
     write(0,*) ''
     stop 
end if

end subroutine handle_err

end program driver

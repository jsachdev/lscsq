program driver

  use pl_types
  use RayIni
  use iso_c_binding, only: fp => c_double, c_int
  use lscsq_mod
  use netcdf ! to read from CDF file
  implicit none

  character(len=80) :: cdffile
  !============
  integer :: ierr, i

  integer :: iraytr=1
  type(storeAry), dimension(:, :), allocatable :: arrys
  




  cdffile = 'cur_state.cdf'
  call lscsq_readpstate(cdffile)

  call alloc_profs

  rmin = lh_inp%Rmin
  rmax = lh_inp%Rmax
  zmin = lh_inp%Zmin
  zmax = lh_inp%Zmax

  dx_grid = lh_inp%dx_grid
  dz_grid = lh_inp%dz_grid

  xmag = lh_inp%Raxis
  zmag = lh_inp%Zaxis
  
  allocate(psigrd(nx,nz))
  psigrd = lh_inp%psirz

  ! all arrays with dimension set by the external code (e.g. TRANSP) 
  ! are allocated here

  TekeV = lh_inp%Te
  voltlp = lh_inp%Vloop
  rbpvec = lh_inp%g_eq

  dVlvec = lh_inp%dvol
  iVlvec = lh_inp%vol

  psivec = lh_inp%plflx

  psimin = psivec(1)
  psilim = psivec(npsij)

  !	--------------------------------------------------------------
!  btesl = bcentr
!  rmaj = rcentr 
  rmaj = lh_inp%Raxis

  btesl = lh_inp%B_axis

  ! here all operations
!     pe2Fac pi2Fac convert density in ^14 cm^-3 or ^20 m^-3
!     and frequency in GHz into
!
  pe2fac = 1.0e-5_fp*(qe_eV*vc)**2/(pi*me_Kg) !/ 100000
  pe2Fac14= 1.0e-19_fp*(qe_eV*vc)**2/(PI*me_Kg)
  pi2fac  = pe2fac *(me_Kg/mp_Kg) 
  AioFac  = 3.0e-15_fp * pi2Fac * qe_eV/mp_Kg/TWOPI**2
  AelFac  = 0.75e-15_fp * pe2Fac * qe_eV/me_Kg/TWOPI**2 
  OmcFac  = 1.0e-9_fp*qe_eV/me_Kg/TWOPI 
  ceifac = 1.0e-18_fp*(qe_eV/twopi)**2/me_Kg/mp_Kg
  cEparIK = 1.0e2_fp*twopi*twopi*qe_eV**2/me_g
  midvec = psivec


  EdcVec(1:npsij) = voltlp(1:npsij)/(twopi*rmaj)

  pe2vec = 1.0e-20_fp*pe2fac*lh_inp%ne
  AelVec(1:npsij) = 1.0e-20_fp*AelFac*lh_inp%ne(1:npsij)*lh_inp%Te(1:npsij)
  ! we should sum over ion species, including impurities. However, LH does not
  ! heat on ions. Take for now only the background species (1 species for test)
  pi2Vec(1:npsij) = 1.0e-20_fp*pi2Fac*lh_inp%ni(1:npsij)*(lh_inp%chrg/qe_eV)**2*mp_Kg/lh_inp%mass
  AioVec(1:npsij) = 1.0e-20_fp*AioFac*lh_inp%ni(1:npsij)*(lh_inp%chrg/qe_eV)**2*(mp_Kg/lh_inp%mass)**2*lh_inp%Ti(1:npsij)

  !wcei2 = ceifac*lh_inp%ni(1)/lh_inp%ne(1)*(lh_inp%chrg/qe_eV)**2*mp_kg/lh_inp%mass
  !wcei2 = wcei2*Btesl**2

  RlcfsMax = lh_inp%Rmaxlcfs
  RlcfsMin = lh_inp%Rminlcfs
  zlcfsMax = lh_inp%Zmaxlcfs
  zlcfsMin = lh_inp%Zminlcfs

  call nml2lsc
  

  i1stcall = 1
  iraytr = 1  
!  call lscsq_main(LhPwrMW, iraytr,iError)

  allocate(arrys(nzones, nrays))
  call lscsq_main(iraytr,iError, arrys)
  call lscsq_writecdf(arrys)

  call dealloc_profs
  
contains


subroutine lscsq_readpstate(pstate_file)
  use netcdf
  use lscsq_mod 
  use Tracing
  implicit none

  character(len=80), intent(in) :: pstate_file
  integer :: ncid, status
  integer :: rhodimid, reqdimid, dim_nrho, vid, nsadimid, lhdim
  integer :: dm1_nrho, rh1dimid, nsdimid, dp1_nspec_th, dim_nrho_eq
  integer :: dp1_nspec_tha, dim_nr, dim_nz, rdimid, zdimid, dim_nlhrf_src
  real(fp) :: dum0
  real(fp), dimension(:), allocatable :: dums, dumsa, dum_r, dum_z,dumlh
  real(fp), dimension(:), allocatable :: dumZ, rho_eq, dum1, dum1p
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
  allocate(lh_inp%powerlh(dim_nlhrf_src))
  allocate(lh_inp%freqlh(dim_nlhrf_src))
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
  lh_inp%dvol(1) = 0.0
  do i=2,dim_nrho_eq
     lh_inp%dvol(i) = lh_inp%vol(i)-lh_inp%vol(i-1) 
  enddo

  allocate(dumrz(dim_nr,dim_nz))
  ierr = nf90_inq_varid(ncid,'PsiRZ',vid)
  status = nf90_get_var(ncid,vid,dumrz) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%psirz = dumrz


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
  lh_inp%Te(1) = dum2d(1,1)+0.5*(dum2d(1,1)-dum2d(2,1)) 
  do i=1,dm1_nrho-1
     lh_inp%Te(i+1) = 0.5*(dum2d(i,1)+dum2d(i+1,1)) 
  enddo 

  ierr = nf90_inq_varid(ncid,'Ti',vid)
  status = nf90_get_var(ncid,vid,dum1) 
  if (status /= nf90_noerr) call handle_err(status)
!  from zone center (?) to zone boundary
  lh_inp%Ti(1) = dum1(1)+0.5*(dum1(1)-dum1(2)) 
  do i=1,dm1_nrho-1
     lh_inp%Ti(i+1) = 0.5*(dum1(i)+dum1(i+1)) 
  enddo 

  ierr = nf90_inq_varid(ncid,'ns_bdy',vid)
  status = nf90_get_var(ncid,vid,dumZ) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%ne(npsij) = dumz(1)
  lh_inp%ni(npsij) = sum(dumz(2:dp1_nspec_th))

  ierr = nf90_inq_varid(ncid,'ns',vid)
  status = nf90_get_var(ncid,vid,dum2d) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%ne(1) = dum2d(1,1)+0.5*(dum2d(1,1)-dum2d(2,1)) 
  do i=1,dm1_nrho-1
     lh_inp%ne(i+1) = 0.5*(dum2d(i,1)+dum2d(i+1,1)) 
  enddo 

  ierr = nf90_inq_varid(ncid,'ni',vid)
  status = nf90_get_var(ncid,vid,dum1) 
  if (status /= nf90_noerr) call handle_err(status)
!  from zone center (?) to zone boundary
  lh_inp%ni(1) = dum1(1)+0.5*(dum1(1)-dum1(2)) 
  do i=1,dm1_nrho-1
     lh_inp%ni(i+1) = 0.5*(dum1(i)+dum1(i+1)) 
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
!  write(*,*) 'q_Sa:',dumsa/1.602e-19

  ierr = nf90_inq_varid(ncid,'m_S',vid)
  status = nf90_get_var(ncid,vid,dums) 
  if (status /= nf90_noerr) call handle_err(status)
!  write(*,*) 'q_S:',dums     
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

  allocate(dum_r(dim_nr)) 
  ierr = nf90_inq_varid(ncid,'R_grid',vid)
  status = nf90_get_var(ncid,vid,dum_r) 
  if (status /= nf90_noerr) call handle_err(status)
  lh_inp%dx_grid = (dum_r(dim_nr)-dum_r(1))/(dim_nr-1)

  close(ncid)

end subroutine lscsq_readpstate



subroutine lscsq_writecdf(arrys)
  use pl_types 
  use netcdf
  use lscsq_mod
  implicit none
  type(storeAry), dimension(nzones, nrays), intent(inout) :: arrys
  character(len=20) :: cdf_out = 'lsc_out.cdf'
  character(len=*), parameter :: UNITS = "units"

  integer :: xdimid, rdimid,zdimid, nzonid, nvid, nrayid
  integer :: sigdimid,nbdryid
  integer, dimension(3) :: dimindf
  integer, dimension(2) :: dimindr, dimvx
  integer :: ierr, ncid!, i
  integer :: vid_psi, vid_r, vid_z,  vid_vel
  integer :: vid_rray, vid_zray, vid_npar, vid_nper
  integer :: vid_pwry, vid_tmry, vid_dsry, vid_detry, vid_nery, vid_btry
  integer :: vid_rtpy, vid_pray,vid_edcin,vid_edcou
  integer :: vid_psiou, vid_midou, vid_psiin, vid_midin
  integer :: vid_js, vid_jsp, vid_fe
  integer :: vid_vth, vid_neou, vid_teou, vid_nuc, vid_dc, vid_dql
  integer :: vid_pow, vid_pql, vid_nzon
  integer ::  vid_rbdry, vid_zbdry

  ! -----------------------------------------------------
  ! Create netCDF file
  ierr=nf90_create(trim(cdf_out),nf90_clobber,ncid)
  ! --------------------------------------
  ! dimensions 
  ierr=nf90_def_dim(ncid,'npsi',npsi,xdimid)
!  ierr=nf90_def_dim(ncid,'npsij',npsij,xxdimid)
  ierr=nf90_def_dim(ncid,'nr',nx,rdimid)
  ierr=nf90_def_dim(ncid,'nz',nz,zdimid)
  ierr=nf90_def_dim(ncid,'nzon',nzones,nzonid)
  ierr=nf90_def_dim(ncid,'nv',nv,nvid)
!  ierr=nf90_def_dim(ncid,'nx',npsij-1,nxdim)
  ierr=nf90_def_dim(ncid,'nrays',nrays,nrayid)
  ierr=nf90_def_dim(ncid,'sigma',2,sigdimid) 
!  ierr=nf90_def_dim(ncid,'nbdry',nlcfs,nbdryid) 

  dimindf=(/ nvid,xdimid,sigdimid /)
  dimindr=(/ nzonid,nrayid /)
  dimvx = (/ nvid,xdimid /)
  ! --------------------------------------
  ! definitions
  ierr=nf90_def_var(ncid,'psi',nf90_real,xdimid,vid_psi)
  ierr=nf90_def_var(ncid,'R',nf90_real,rdimid,vid_r)
  ierr=nf90_def_var(ncid,'Z',nf90_real,zdimid,vid_z)
  ierr=nf90_def_var(ncid,'vpar',nf90_real,nvid,vid_vel)
!  ierr=nf90_def_var(ncid,'x',nf90_real,nxdim,vid_x)
 
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
!  ierr=nf90_def_var(ncid,'Edc_in',nf90_real,xxdimid,vid_edcin)
  ierr=nf90_def_var(ncid,'psi_ou',nf90_real,xdimid,vid_psiou)
!  ierr=nf90_def_var(ncid,'psi_in',nf90_real,xxdimid,vid_psiin)
  ierr=nf90_def_var(ncid,'mid_ou',nf90_real,xdimid,vid_midou)
!  ierr=nf90_def_var(ncid,'mid_in',nf90_real,xxdimid,vid_midin)
  ierr=nf90_def_var(ncid,'Js',nf90_real,xdimid,vid_js)
  ierr=nf90_def_var(ncid,'Jsp',nf90_real,xdimid,vid_jsp)
  ierr=nf90_def_var(ncid,'Power',nf90_real,xdimid,vid_pow)
  ierr=nf90_def_var(ncid,'Power_ql',nf90_real,xdimid,vid_pql)
  ierr=nf90_def_var(ncid,'fe',nf90_real,dimindf,vid_fe)
  ierr=nf90_def_var(ncid,'Dql',nf90_real,dimindf,vid_dql)
  ierr=nf90_def_var(ncid,'vth',nf90_real,xdimid,vid_vth)
  ierr=nf90_def_var(ncid,'Dcoll',nf90_real,dimvx,vid_dc)
  ierr=nf90_def_var(ncid,'nucoll',nf90_real,dimvx,vid_nuc)
!  ierr=nf90_def_var(ncid,'powtr',nf90_real,xxdimid,vid_powtr)
!  ierr=nf90_def_var(ncid,'curtr',nf90_real,xxdimid,vid_curtr)
  ierr=nf90_def_var(ncid,'nzon',nf90_real,nrayid,vid_nzon)


  ierr=nf90_enddef(ncid)
 
  ! --------------------------------------
  ! Write variables to CDF file
!  ierr=nf90_put_var(ncid,vid_r,xary)
!  ierr=nf90_put_var(ncid,vid_z,zary)
  ierr=nf90_put_var(ncid,vid_psi,midary)
  if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr=nf90_put_var(ncid,vid_rray,arrys(1:nzones,1:nrays)%rofray)
  ierr=nf90_put_var(ncid,vid_zray,arrys(1:nzones,1:nrays)%zofray)
  ierr=nf90_put_var(ncid,vid_npar,arrys(1:nzones,1:nrays)%nparry)
  ierr=nf90_put_var(ncid,vid_nper,arrys(1:nzones,1:nrays)%nperry)
  !ierr=nf90_put_var(ncid,vid_nper,arrys%nperry)
  ierr=nf90_put_var(ncid,vid_pray,arrys(1:nzones,1:nrays)%Pofray)
  ierr=nf90_put_var(ncid,vid_pwry,arrys(1:nzones,1:nrays)%Powrry)
  ierr=nf90_put_var(ncid,vid_tmry,arrys(1:nzones,1:nrays)%Timery)
  ierr=nf90_put_var(ncid,vid_dsry,arrys(1:nzones,1:nrays)%distry)
  ierr=nf90_put_var(ncid,vid_detry,arrys(1:nzones,1:nrays)%detrry)
  ierr=nf90_put_var(ncid,vid_nery,arrys(1:nzones,1:nrays)%Neofry)
  ierr=nf90_put_var(ncid,vid_Btry,arrys(1:nzones,1:nrays)%Bthray)
  ierr=nf90_put_var(ncid,vid_rtpy,arrys(1:nzones,1:nrays)%rtpsry)
  
  ierr=nf90_put_var(ncid,vid_neou,neary)
  ierr=nf90_put_var(ncid,vid_teou,teary)
  ierr=nf90_put_var(ncid,vid_psiou,psiary)
  ierr=nf90_put_var(ncid,vid_psiin,psivec)
  ierr=nf90_put_var(ncid,vid_Edcou,Edcary)
  ierr=nf90_put_var(ncid,vid_Edcin,Edcvec)
  ierr=nf90_put_var(ncid,vid_midin,midvec)
  ierr=nf90_put_var(ncid,vid_midou,midary)

  ierr=nf90_put_var(ncid,vid_js,js)
  ierr=nf90_put_var(ncid,vid_jsp,jsp)
  ierr=nf90_put_var(ncid,vid_pow,praytot)
  ierr=nf90_put_var(ncid,vid_pql,pqltot)
  ierr=nf90_put_var(ncid,vid_fe,fe)
  ierr=nf90_put_var(ncid,vid_vth,vtherm)
  ierr=nf90_put_var(ncid,vid_vel,vpar)
  ierr=nf90_put_var(ncid,vid_dc,dcoll)
!  if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr=nf90_put_var(ncid,vid_nuc,nucoll)
  ierr=nf90_put_var(ncid,vid_dql,Dql)
!  if (ierr /= nf90_noerr) call handle_err(ierr)
!  ierr=nf90_put_var(ncid,vid_powtr,powtsc)
!  ierr=nf90_put_var(ncid,vid_curtr,curtsc)
!  ierr=nf90_put_var(ncid,vid_x,x)
  ierr=nf90_put_var(ncid,vid_nzon,nz_ind)
!  ierr=nf90_put_var(ncid,vid_rbdry,rlcfs)
!  ierr=nf90_put_var(ncid,vid_zbdry,zlcfs)

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

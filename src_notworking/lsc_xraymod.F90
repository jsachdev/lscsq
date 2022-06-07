module lscsq_xraymod

  use iso_c_binding, only : fp => c_double
  implicit none

! parameters
  integer, parameter :: nmudim=21
  integer, parameter :: nendim=200
  integer, parameter :: maxpoints=1000
  integer, parameter :: spacedim=3
  integer, parameter :: npixdim=50  
  integer, parameter :: nchordim=6     
  integer, parameter :: nrdim=90       
  integer, parameter :: nzdim=90       
  integer, parameter :: nrzmxdim=90    
  integer, parameter :: nwkdix=1000
  integer, parameter :: npsidim=40  ! temporary
  integer, parameter :: nveldim=199 ! temporary
  

  integer n_pixel(SPACEDIM - 1)
  integer :: n_pixel_x, n_pixel_y
  integer :: nMUbins, nEbins, iAbsXray
  integer, dimension(nveldim) :: Eint
  integer, dimension(npixdim,npixdim) :: npoints
  integer :: nr_source, nz_source
  integer, dimension(nwkdix) :: iwkarx  

  real(fp), dimension(nwkdix,4) :: wkarx

  REAL(fp) :: ACCEPTABLE_ERROR
  REAL(fp) :: mu_min, mu_max, E_min, E_max, E_ph_min, dE_ph, dFoilTCM
  REAL(fp) :: dmu_inv
 
  REAL(fp) ::    mu_0, mu_width,                                           &
     &        Z_bound_min, Z_bound_max, R_bound_min, R_bound_max,       &
     &        Z_plasm_min, Z_plasm_max, R_plasm_min, R_plasm_max,       &
     &        r_source(NRDIM), z_source(NZDIM),                         &
     &        R_bound_min_sq, R_bound_max_sq
      REAL(fp) ::    PusherMajor, PusherMinor
      REAL(fp) ::    source_profile(NRDIM, NZDIM)

  real(fp), dimension(npixdim,npixdim) :: photon_count
  real(fp), dimension(3,maxpoints,npixdim,npixdim) :: chord
  REAL(fp) ::    pinhole_loc(SPACEDIM), camera_orientation(2),             &
              focal_length, screen_d,                                   &
              pinhole_size, pixel_size(SPACEDIM - 1),                   &
              R_tangent, z_tangent, RcntChrd,                           &
              mu_axis, phi_axis, Rpinhole,                              &
              Zpinhole, phi_pinhole, pix_size_x, pix_size_y,            &
              y_hat_camera(SPACEDIM, NPIXDIM, NPIXDIM),                 &
              chord_origin(SPACEDIM, NPIXDIM, NPIXDIM),                 &
              x_pix_vec(NPIXDIM), y_pix_vec(NPIXDIM),                   &
              r_tangent_actual(NPIXDIM, NPIXDIM),                       &
              z_tangent_actual(NPIXDIM, NPIXDIM)
  REAL(fp) :: sigtot(NMUDIM, NENDIM), E_incvec(NENDIM), XmnFac(NENDIM), &
              muvec(NMUDIM), Efrac(NVELDIM), Eofv(NVELDIM), dEofv(NVELDIM)
  REAL(fp) :: inten(NMUDIM, NPSIDIM) 

  CHARACTER(len=2) :: FoilCode

      EQUIVALENCE (n_pixel(1), n_pixel_x), (n_pixel(2), n_pixel_y)
      EQUIVALENCE (camera_orientation(1), mu_axis),                     &
     &     (camera_orientation(2), phi_axis)
      EQUIVALENCE (pinhole_loc(1), Rpinhole),                           &
     &     (pinhole_loc(2), Zpinhole),                                  &
     &     (pinhole_loc(3), phi_pinhole),                               &
     &     (pixel_size(1), pix_size_x), (pixel_size(2), pix_size_y)
!      COMMON /lscsq_camcom0/ n_pixel, npoints
!      COMMON /lscsq_camcom1/ pinhole_loc, camera_orientation,               &
!     &     focal_length, screen_d,                                      &
!     &     pinhole_size, pixel_size, R_tangent, z_tangent, RcntChrd,    &
!     &     y_hat_camera, x_pix_vec, y_pix_vec, chord_origin,            &
!     &     r_tangent_actual, z_tangent_actual
!      COMMON /lscsq_camcom2/ chord, photon_count


!     Xray.inc
!     E_photon        in MeV
!     E_incvec        in MeV
!     E_min E_max
!     COMMON /lscsq_Xray1a/ nMUbins, nEbins, iAbsXray
!     COMMON /lscsq_Xray2b/ Eint
!     COMMON /lscsq_Xray3c/ sigtot, E_incvec, XmnFac, muvec, inten, Efrac
!     COMMON /lscsq_Xray4d/ E_ph_min, dE_ph
!     COMMON /lscsq_Xray5e/ dFoilTCM
!     COMMON /lscsq_Xray6f/ Eofv, dEofv
!     COMMON /lscsq_Xray7g/ mu_min, mu_max, dmu_inv, E_min, E_max
!     COMMON /lscsq_Xray8h/ FoilCode
!                       ! dFoilTCM; or; dcu: thickness of foil (cm)
!                       ! FoilCode; or; FC: CU,AG,TA,MO, or 00 for no foil
!     Xray.inc
!     .         --------------------------------------------------------
!
!     COMMON /lscsq_emcom1/ nr_source, nz_source
!     COMMON /lscsq_emcom2/ mu_0, mu_width,                                 &
!    &        Z_bound_min, Z_bound_max, R_bound_min, R_bound_max,       &
!    &        Z_plasm_min, Z_plasm_max, R_plasm_min, R_plasm_max,       &
!    &     r_source, z_source, R_bound_min_sq,                          &
!    &     R_bound_max_sq, PusherMajor, PusherMinor
!     COMMON /lscsq_emcom3/ source_profile


      EQUIVALENCE (n_pixel(1), n_pixel_x), (n_pixel(2), n_pixel_y)
      EQUIVALENCE (camera_orientation(1), mu_axis),                     &
     &     (camera_orientation(2), phi_axis)
      EQUIVALENCE (pinhole_loc(1), Rpinhole),                           &
     &     (pinhole_loc(2), Zpinhole),                                  &
     &     (pinhole_loc(3), phi_pinhole),                               &
     &     (pixel_size(1), pix_size_x), (pixel_size(2), pix_size_y)
!      COMMON /lscsq_camcom0/ n_pixel, npoints
!      COMMON /lscsq_camcom1/ pinhole_loc, camera_orientation,               &
!     &     focal_length, screen_d,                                      &
!     &     pinhole_size, pixel_size, R_tangent, z_tangent, RcntChrd,    &
!     &     y_hat_camera, x_pix_vec, y_pix_vec, chord_origin,            &
!     &     r_tangent_actual, z_tangent_actual
!      COMMON /lscsq_camcom2/ chord, photon_count
!     camera.inc --------------------------------------------------------------


end module lscsq_xraymod

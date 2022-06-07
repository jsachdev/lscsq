SUBROUTINE lscsq_smooth(rmat, n1, n2, nsmoo, smvec)

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : nv
  implicit none

  integer, intent(in) :: nsmoo
  integer, intent(in) :: n1
  integer, intent(in) :: n2

  integer :: ipsi, nsm2 
  real(fp), dimension(n1,n2), intent(inout) :: rmat
  real(fp), dimension(n1), intent(in) :: smvec 
  real(fp), dimension(nv) :: tmp
!     smooth quasilinear diffusion coefficient
!     VecDmp in matr.F puts the Dql values for each v into wkv.
!
!     The calling logic is tricky in that it depends on storage order in f77.
!     convolve takes wkv, outputs Dql.  If qlsm is a delta fn, no change.
!     qlsm is integer vector gaussian of width nsmw, normalized to add to 1.

  nsm2 = (nsmoo-1)/2
  do ipsi=1, n2
     tmp(1:n1) = rmat(1:n1,ipsi)
     CALL lscsq_convolve(n1, nsm2, smvec, rmat(1:n1, ipsi), tmp)
  enddo
 
end subroutine lscsq_smooth
!     -----------------------------------------------------------------
subroutine lscsq_DqlGen
  ! deposit wave energy into dql matrix
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : power, ezsq, ivind, izind, rfudgdmp
  use lscsq_mod, only : fghz_now, vpar, fghz !, dvip 
  use lscsq_mod, only : dql, dqlnorm, qlsm, nv, npsi, nsmoo
  use lscsq_mod, only: nrays, npols, nzones, ind_ray, ok_ray
  implicit none

  integer :: iv, ip
  integer :: iry, izn !, itor
  real(fp) :: dv
!     deposit energy density on ray iray and zone izone
!     into appropriate
!     velocity bin.
!     See eq 22 and 28 of Valeo/Eder.
!     D_{ql} = \pi/2 (e/m)^2 Sum_j <E_z_j^2> \delta(\omega - k_par v_par)
!     where the j represents all rays and intersections with the flux surface
!     in question (index supressed here) and the <> mean flux surface avg.
!     The factor 2 is appropriate for E meaning an amplitude; it makes the
!     E^2 into an rms quantity.
!     The delta function connects the k-par dependence of E with the v
!     dependence of D.
!     \delta(w) = (2\pi h)^{-.5} exp(-(w/h)^2/2.) where h is a variable
!     width.  By writing in terms of velocity index and nsmoo
!
!     D(v_i) = \pi/2 (e/m)^2 Sum_j E_j^2 (v_par/(\omega \Delta v_par) *
!         (2\pi nsmoo)^{-.5} exp (- (i-j)^2/(2 nsmoo^2))
!
!     .                                 Clear the unsmoothed part of Dql

  dql(1:nv,1:npsi,1:2) = 0.0 
  ! Fill up the unsmoothed part of Dql
  do iry  = 1, nrays
!   if (ok_ray(iry).eq.1) then ! probably not needed
!      itor = ind_ray(1,iry)
      do izn   = 1, nzones-1
        iv = ivind(izn,iry)
        ! If  ivind=0, the ray was stopped before reaching this izone (izn).
        !    If Power=0, no contribution to Dql. No calculation is appropriate.
        if (iv.EQ.0 .or. Power(izn,iry).EQ.0.0) cycle    
        ip = izind(izn,iry)
        if (iv.eq.nv) dv = vpar(iv)-vpar(iv-1)
        if (iv.lt.nv) dv = vpar(iv+1)-vpar(iv)
        ! DMC: 1/frequency factor added here; see comments in DqlNorm subroutine
        Dql(iv,ip,1) = Dql(iv,ip,1)+Power(izn,iry)*rFudgDmp(izn,iry)   &
                     * Ezsq(izn,iry)*abs(vpar(iv))/dv/fghz(iry) 
      enddo        
!   endif
  enddo

  ! Multiply by normalization
  Dql = Dqlnorm*Dql  

  ! Copy unsmoothed into space for smoothed
  Dql(1:nv,1:npsi,2) = Dql(1:nv,1:npsi,1)

  ! Smooth the new copy
  call lscsq_Smooth(Dql(1:nv,1:npsi,2), nv, npsi, nsmoo, qlsm)

end subroutine lscsq_dqlgen    
!     -----------------------------------------------------------------
INTEGER FUNCTION lscsq_ivtabl(np)
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : vpar, ivzero, nv
  implicit none

  INTEGER :: iv
  REAL(fp) :: np, vv            
      
  vv = 1.0_fp / np

  iv = minloc(abs(vpar-vv),1)
  if (vpar(iv).lt.vv) iv=iv+1 

! v grid symmetry question: use vpar nearby of smaller energy
  if( iv .GT. ivZero) then
     lscsq_ivtabl = iv-1
  else
     lscsq_ivtabl = iv
  endif

end function lscsq_ivtabl
! -----------------------------------------------------------------
subroutine lscsq_DqlInit
  ! initialize d quasilinear constants
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : dqlnorm
  use lscsq_mod, only : pi, twopi, vc, qe_eV, me_Kg
  implicit none

!     Valeo original:
!     DqlNorm = (PI / 2.) * ((ECHARG / EMASS) ** 2) * EcgsbyEmks *
!    ^     EcgsbyEmks / (omega * omega) * (3.e10)
!     see Valeo-Eder, Eq. (22)
!     but why did he have omega twice
!     Here is an MKS version based on
!     Dql = \Sigma (\pi/2) q^2/m^2 E_z^2 /({\Delta k_\parallel} v_\parallel)
!         = \Sigma (\pi/2) q^2/m^2 E_z^2 / omega * v/delta-v
!         = \Sigma DqlNorm E_z^2 /(v/c *delta-n)
!
! DMC: replaced omega = 2*pi*1e9*fghz_now with 2*pi*1e9; fghz_now can vary
!      with ray now, so, this variation is now in the summations for the
!      individual dql array elements;
!         (pi/2) / (2*pi*1e9) -> (1/4e9)
!
  DqlNorm = 2.5e-10_fp*(qe_eV/me_Kg)**2 / vc**2
!      
end subroutine lscsq_DqlInit

!---------------------------------------------------
SUBROUTINE lscsq_convolve(n, nsm, smvec, vout, vin)

  use iso_c_binding, only : fp => c_double
  implicit none

  integer, intent(in) :: n, nsm

  real(fp), dimension(n), intent(in) :: vin
  real(fp), dimension(n), intent(out) :: vout
  real(fp), intent(in) :: smvec(-nsm:nsm)

  integer :: i, j, ilow, iup, jlow, jup
!
!     convolute source vin with smvec to yield vout
!
!     interior points
  ilow = nsm + 1
  iup = n - nsm
  jlow = - nsm
  jup = nsm
  do i = ilow, iup
     vout(i) = 0.0 
     do j = jlow, jup
        vout(i) = vout(i)+smvec(j)*vin(i-j)
     enddo        
  enddo        
!
!     lower boundary
  ilow = 1
  iup = nsm
  jlow = - nsm
  do i = ilow, iup
     vout(i) = 0.0 
     jup = i - 1
     do j = jlow, jup
        vout(i) = vout(i) + smvec(j) * vin(i - j)
     enddo          
  enddo       
!
!     upper boundary
  ilow = n - nsm + 1
  iup = n
  jup = nsm
  do i = ilow, iup
     vout(i) = 0.0 
    jlow = n - i
    do j = jlow, jup
       vout(i) = vout(i) + smvec(j) * vin(i - j)
    enddo         
  enddo       
!
end subroutine lscsq_convolve


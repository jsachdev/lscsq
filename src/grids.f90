SUBROUTINE lscsq_mkgrids

  use lscsq_mod

  CALL lscsq_psigrids ! generate psi grid
  CALL lscsq_vgrids   ! generate parallel velocity grid
  CALL lscsq_mkqlsm   ! generate smoothing function
END SUBROUTINE lscsq_mkgrids
!     ------------------------------------------------------------------
SUBROUTINE lscsq_mkqlsm

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : nsmoo, nsmw, qlsm, nsmsym
  implicit none

  INTEGER :: iv
  REAL(fp) :: earg, denomql

  ! initialize quasilinear smoothing function
  nsmsym = (nsmoo+1)/2
  denomql = 1.0_fp/REAL(nsmw,kind=fp)
  do iv = 1, nsmoo
     earg = REAL(iv-nsmsym,kind=fp)*denomql
     earg = 0.5_fp*earg*earg 
     qlsm(iv) = exp(-earg)
  enddo

  qlsm = qlsm/sum(qlsm) 

END SUBROUTINE lscsq_mkqlsm
!     ------------------------------------------------------------------
SUBROUTINE lscsq_psigrids
  ! generate psi grid
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : npsij, npsi, psiary, midary
  use lscsq_mod, only : delpsi !, psimin,psilim
  use lscsq_mod, only: lh_inp, lh_out
  implicit none

  REAL(fp) :: dpp
  INTEGER :: j
  real(fp), parameter :: ttiny=1.0e-5_fp

  dpp = (lh_inp%plflx(npsij)-lh_inp%plflx(1))*TTINY/REAL(npsi-1,kind=fp)
  ! generate psiary grid
  dpp = 0.0_fp
  CALL lscsq_ugrid(PsiAry, npsi, lh_inp%plflx(1)-dpp, lh_inp%plflx(npsij)+dpp) ! ENFORCES INTERPOLATION !!
  lh_out%psi = psiary

  delpsi = abs(PsiAry(2)-PsiAry(1))
  do j=1,Npsi-1
     MidAry(j) = 0.5_fp * (PsiAry(j+1)+PsiAry(j))
  enddo
  MidAry(Npsi) = PsiAry(Npsi)

  CALL lscsq_Volcalc


END SUBROUTINE lscsq_psigrids
!     ------------------------------------------------------------------
SUBROUTINE lscsq_VolCalc

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : npsi,npsij, midary, ivlary, dvol
  use lscsq_mod, only: lh_inp
  implicit none

  INTEGER :: ips
  REAL(fp) :: xlookup, yreturn

  iVlAry(npsi) = lh_inp%vol(npsij)

  ! fmp - I suspect there is an inconsistency between the grids that are input
  ! and the interpolated grids, in terms of zone center and zone boundary. TBC 
  ! compute integral of volume by interpolation
  do ips = 1, Npsi-1
     xlookup = MidAry(ips)
     CALL lscsq_linr1d(NpsiJ, lh_inp%plflx, lh_inp%vol, xlookup, yreturn)
     iVlAry(ips) = yreturn
  enddo

  ! compute dVol which is centered
  dVol(1) = iVlAry(1)
  do ips = 2, Npsi
     dVol(ips) = iVlAry(ips) - iVlAry(ips-1)
  enddo

END SUBROUTINE lscsq_VolCalc
!     ------------------------------------------------------------------
SUBROUTINE lscsq_vgrids

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : nv, vmin, vmax, ivzero, dvplus !, dvsym, dvip
  use lscsq_mod, only: lh_out
  implicit none      

  integer :: i
  real(fp), dimension(nv) :: dumr

  dvplus(1:nv)=0.0_fp

  ! set UniformGRID into Vpar
  CALL lscsq_ugrid(dumr, nv, Vmin, Vmax)
  lh_out%vpar = dumr
 
  ! compute forward difference on Vpar 
  CALL lscsq_mkdvp(nv, dvplus, lh_out%vpar)

  ! construct inverse, so dvisym is 1/dv
  IvZero = (nv+1)/2   ! requires nv to be odd !

end subroutine lscsq_vgrids

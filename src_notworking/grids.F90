SUBROUTINE lscsq_mkgrids

  use lscsq_mod
 
  CALL lscsq_psigrids ! generate psi grid
  CALL lscsq_vgrids   ! generate parallel velocity grid
  CALL lscsq_mkqlsm   ! generate smoothing function
  call lscsq_thgrid   ! generate array of poloidal angles
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
  use lscsq_mod, only : npsi, psiary, midary
  use lscsq_mod, only : delpsi, psimin,psilim
  implicit none

  REAL(fp) :: dpp
  INTEGER :: j
  real(fp), parameter :: ttiny=1.0e-5_fp

  dpp = (psilim-psimin)*TTINY/REAL(npsi-1,kind=fp)
  ! generate psiary grid
  CALL lscsq_ugrid(PsiAry, npsi, psimin-dpp, psilim+dpp) ! ENFORCES INTERPOLATION !!

  delpsi = abs(PsiAry(2)-PsiAry(1))
  do j=1,Npsi-1
     MidAry(j) = 0.5_fp * (PsiAry(j+1)+PsiAry(j))
  enddo
  MidAry(Npsi) = PsiAry(Npsi)

  CALL lscsq_Volcalc


END SUBROUTINE lscsq_psigrids
!     ------------------------------------------------------------------
SUBROUTINE lscsq_VolCalc
  use grap
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : npsi,npsij, midary, midvec, ivlvec, ivlary, dvol
  implicit none


!!!!!!     EXTERNAL ugrid
  INTEGER :: ips
  REAL(fp) :: xlookup, yreturn

  iVlAry(npsi) = iVlVec(NpsiJ)
 
  ! compute integral of volume by interpolation
  do ips = 1, Npsi-1
     xlookup = MidAry(ips)
     CALL lscsq_linr1d(NpsiJ, MidVec, iVlVec, xlookup, yreturn)
     iVlAry(ips) = yreturn
  enddo

  ! compute dVol which is centered
  dVol(1) = iVlAry(1)
  do ips = 2, Npsi
     dVol(ips) = iVlAry(ips) - iVlAry(ips-1)
  enddo

END SUBROUTINE lscsq_VolCalc
!     ------------------------------------------------------------------
subroutine lscsq_thgrid
   use Tracing
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : nth, thet0, dthet, thgrid
  implicit none      

  !integer :: i
  real(fp) :: thmin, thmax

  thmin = thet0-0.5_fp*dthet
  thmax = thet0+0.5_fp*dthet

  thgrid(1:nth)=thet0
  ! generate theta grid
  if ((nth.eq.1).or.(dthet.eq.0.0)) then
     thgrid(1) = thet0
  else
     if (2*(nth/2).EQ.nth) then
        CALL lscsq_LSCwarn (' require ODD nth ')
        nth = nth-1
     endif
     ! set UniformGRID into thgrid
     CALL lscsq_ugrid(thgrid, nth, thmin, thmax)  

  endif

end subroutine lscsq_thgrid
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
SUBROUTINE lscsq_vgrids
   use Tracing
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : nv, vmin, vmax, vpar, ivzero, dvplus !, dvsym, dvip
  implicit none      

 ! integer :: i

  dvplus(1:nv)=0.0_fp

  ! generate parallel velocity grid
  if (2*(nv/2).EQ.nv) then
     CALL lscsq_LSCwarn (' require ODD nv ')
     nv = nv-1
  endif
  ! set UniformGRID into Vpar
  CALL lscsq_ugrid(Vpar, nv, Vmin, Vmax)
 
  ! compute forward difference on Vpar 
  CALL lscsq_mkdvp(nv, dvplus, Vpar)

  ! construct inverse, so dvisym is 1/dv
  IvZero = (nv+1)/2   ! requires nv to be odd !

end subroutine lscsq_vgrids
!     ------------------------------------------------------------------
subroutine lscsq_MiscInit
  
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: mp_g, zcm3tom3, zcmtom, zev2kev
  use lscsq_mod, only : pi, twopi, zmtocm
  use lscsq_mod, only : qe_eV, vc, me_Kg, mp_Kg, me_g, mp_g
  use lscsq_mod, only : pe2fac,pe2fac14,pi2fac,aiofac
  use lscsq_mod, only : aelfac, omcfac, ceifac, cEparIK
  implicit none      

  pe2fac = 1.0e-5_fp*(qe_eV*vc)**2/(pi*me_Kg) 
  pe2Fac14= 1.0e-19_fp*(qe_eV*vc )**2/(PI*me_Kg)
  pi2fac  = pe2fac *(me_Kg/mp_Kg) 
  AioFac  = 3.0e-15_fp * pi2Fac * qe_eV/mp_Kg/TWOPI**2
  AelFac  = 0.75e-15_fp * pe2Fac * qe_eV/me_Kg/TWOPI**2 
  OmcFac  = 1.0e-9_fp*qe_eV/me_Kg/TWOPI 
  ceifac = 1.0e-18*(qe_eV/twopi)**2/me_Kg/mp_Kg
  cEparIK = 1.0e2_fp*twopi*twopi*qe_eV**2/me_g

 end subroutine lscsq_MiscInit
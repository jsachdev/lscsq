SUBROUTINE lscsq_RampPwr(arrys)
  use pl_types 
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : nrampup
  use lscsq_mod, only : iitr
  use lscsq_mod, only : pwrlevel
  use lscsq_mod, only : fe, dql, dfdv
  !use lscsq_mod, only : lh_out, fe
  implicit none
  type(storeAry) :: arrys
  integer :: iramp
  integer :: i1stcall=1
  integer :: iITRsave
  
  CALL lscsq_setpwrl(pwrlevel)

  iITR = 1
  if (i1stcall.eq.1) then
     i1stcall = 0
     CALL lscsq_FeAt0 ! this should be called only the first time LSC is called
  endif

  CALL lscsq_FePrime

  do iramp = 1, nrampup
   ! rescale power level
   call lscsq_rspwr(pwrlevel(iramp))
   call lscsq_RayDamp(arrys)
   call lscsq_DqlGen(arrys)
   iITR = mod(iITR,2) + 1
   call lscsq_FeCalc
   call lscsq_FeWeight
   call lscsq_FePrime
   call lscsq_FeConvrg(iramp)
enddo        
 
  iITRsave = iITR
  iITR = mod(iITR,2)+1
  call lscsq_FeAt0
  iITR = iITRsave

  !lh_out%fe = fe
  !lh_out%dql = dql
  !lh_out%dfdv = dfdv



end subroutine lscsq_rampPwr 
!-----------------------------------------------------------------
subroutine lscsq_FeConvrg(iramp)
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : fe, FeCvgAry, iitr, npsi, nv
  implicit none

  integer, intent(in) :: iramp
  integer :: ip,iv, inew, iold
  real(fp) :: dum

  inew = iITR
  iold = mod(inew,2) + 1

  dum = 0.0  
  do ip=1,npsi
     do iv=2,nv-1
        dum=dum+abs(fe(iv,ip,inew)-fe(iv,ip,iold)) /              &
                       (fe(iv,ip,inew)+fe(iv,ip,iold)+epsilon(1.0_fp))
     enddo
  enddo

  FeCvgAry(iramp) = dum/real(nv*npsi,kind=fp)

end subroutine lscsq_FeConvrg
!------------------------------------------------------------------
subroutine lscsq_rspwr(pwrlev)
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : power, totpwr, spec, npols,nth, ind_ray
  use lscsq_mod, only : nant, spec_ant, power_inp, ntors
  implicit none

  integer :: ir, itor, i1, i2, i
  real(fp), intent(in) :: pwrlev


  do i=1,nant
    Spec = Spec_ant(:,i)
    totpwr = power_inp(i)
    i1 = 1 + (i-1)*ntors*npols*nth
    i2 = i*ntors*npols*nth
    do ir=i1,i2  
!     if (ok_ray(ir).eq.1) then
        itor=ind_ray(1,ir) 
        power(1,ir) = TotPwr*pwrlev*Spec(itor)/real(nth,kind=fp)/real(npols,kind=fp)
!     endif
    enddo
  enddo

end subroutine lscsq_rspwr
!------------------------------------------------------------------
SUBROUTINE lscsq_setpwrl(pwrlevel)
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : nflat, nrampup
  implicit none

  REAL(fp), dimension(nrampup), intent(out) :: pwrlevel

  INTEGER :: i, nfl, n
  INTEGER :: newn 
  integer :: ialgorit = 3
  integer :: ioriginl = 1
  integer :: igeometr = 2
  integer :: ilinear = 3

  ! fmp - right now this is enforced. Consider introducing a new namelist
  ! variable to select the model
  n = nrampup

  if(nflat .LT. 1)then
     nfl = 1
  else
     nfl = nflat
  endif

  if(nflat .GE. n-1)then
     nfl = 1
  else
     nfl = nflat
  endif

  if (iALGORIT .EQ. iORIGINL) then  ! Original code.  Raise by 2 * each time
     pwrlevel(1) = 1.0_fp/(1.0_fp**(n-nfl))
     do i = 2, n-nfl+1
        pwrlevel(i) = 2.0_fp*pwrlevel(i-1)
     enddo
     do i = n-nfl+1, n
        pwrlevel(i) = pwrlevel(n-nfl+1)
     enddo
  else if (iALGORIT .EQ. iGEOMETR) then  ! New Geometric Ramp Up
     newn = n
     do i = newn, newn-nfl, -1
        pwrlevel(i) = 1.0_fp
     enddo

     do i = newn-nfl-1, 1, -1
        pwrlevel(i) = pwrlevel(i+1)/1.15_fp 
     enddo
  else if (iALGORIT .GE. iLINEAR) then  ! Linear Ramp Up
     do i = 1, n-nfl
        pwrlevel(i) = real(i,kind=fp)/real(n-nfl,kind=fp)  
     enddo
     pwrlevel(n-nfl+1:n) = 1.0_fp
  endif

end SUBROUTINE lscsq_setpwrl


subroutine lscsq_ugrid(vector,npts,vmin,vmax)
  use iso_c_binding, only : fp => c_double
  implicit none

  integer, intent(in) :: npts
  real(fp), intent(in) :: vmin, vmax
  real(fp), intent(out), dimension(npts) :: vector

  integer :: j
  real(fp) :: dv
   
  ! generate a uniform grid vector(j) from vmin to vmax
  if(npts.le.1)return
  dv=(vmax-vmin)/real(npts-1,kind=fp)
  do j=1,npts
     vector(j)=vmin+dv*real(j-1,kind=fp)
  enddo

end subroutine lscsq_ugrid

subroutine lscsq_ugridEXC(vector,npts,vmin,vmax,EXCLUDED)

  use iso_c_binding, only : fp => c_double
  implicit none

  integer, intent(in) :: npts
  real(fp), intent(in) :: vmin, vmax, excluded
  real(fp), intent(out), dimension(npts) :: vector

  integer :: ipts
  real(fp) :: dv 
  real(fp) :: vmaxPLUS, vminPLUS, vmaxMINU, vminMINU, v
 
  ! generate a uniform grid vector(j) from vmin to vmax
  ! excluding points with abs value .le. EXCLUDED
 
  if(npts .LE. 1 .or. vmax .LE. vmin ) return
!  rnpts1=real(npts-1,kind=fp)
  vmaxPLUS =  max ( abs(EXCLUDED) , vmax )
  vminPLUS =  max ( abs(EXCLUDED) , vmin )
  vmaxMINU =  min (-abs(EXCLUDED) , vmax )
  vminMINU =  min (-abs(EXCLUDED) , vmin )
  dv       =  max ( (vmaxPLUS-vminPLUS) , 0.0_fp) +                    &
              max ( (vmaxMINU-vminMINU) , 0.0_fp)
  dv = dv/real(npts-1,kind=fp)

  if (dv .LE. 0.0) return

! fmp - need to simplify the logic below and remove the gotos 
  v = vmin
  do ipts=1,npts
     vector(ipts) = v
     v = v + dv
     if (abs(v) .LT. abs(EXCLUDED)) go to 20
  enddo

  return
 
  ! return for normal exit when the excluded zone is not an issue
 
20   continue

  v = vmax
  do ipts=npts,1,-1
     vector(ipts) = v
     v = v - dv
     if (abs(v) .LT. abs(EXCLUDED)) return   
  enddo
 
  ! return for unusual exit when the positive used zone comes
  ! after a negative used zone

end subroutine lscsq_ugridEXC 
 
!------------------------------------------------------------
subroutine lscsq_mkdvp(n, delv, v)
  use iso_c_binding, only : fp => c_double
  implicit none

  integer, intent(in) :: n
  real(fp), intent(in), dimension(n) :: v
  real(fp), intent(out), dimension(n) :: delv

  integer :: i

  do i = 1, n-1
     delv(i) = v(i+1) - v(i)
  enddo

  delv(n) = delv(n-1)

end subroutine lscsq_mkdvp


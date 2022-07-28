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


subroutine spln1d(nr,xgrid,foo)

  use iso_c_binding, only: fp => c_double, c_int
  use lscsq_mod, only: lh_inp

  IMPLICIT NONE

  integer, intent(in) :: nr
  real(fp), dimension(nr), intent(in) :: xgrid
  real(fp), dimension(nr,3), intent(inout) :: foo

  INTEGER :: ispln
  INTEGER :: l,ld,lmax,lmax2,k,j,jm,jmm,jp,jpp,kl,kr
  INTEGER :: km,kmm,kp,kpp,ier,jd,lm,kmid,jsit,kn,nd,jl,jr
!============
  real(fp) :: fmax
  real(fp) :: dx,dum,dum1,ddf,err0
  real(fp) :: yc1,yc2,yc3,af,bf,cf

  real(fp), dimension(:), allocatable :: f1,f2,f3
      
  jl = 2
  jr = nr - 1
      
  ! need grid spacing ...
  dx = xgrid(2)-xgrid(1)

  !   Allocate arrays
  if(allocated(f1)) deallocate(f1,f2,f3)
  allocate(f1(nr),f2(nr),f3(nr))
      
  !   Fill in f1
  f1 = foo(1:nr,1)

  j = 1
  yc1 = f1(j+1)
  yc2 = f1(j)
  yc3 = f1(j+1)
  call get_coeffs(yc1,yc2,yc3,-dx,dx,af,bf,cf)
  f2(j) = bf
  f3(j) = cf
      
  do j = jl,jr
     yc1 = f1(j-1)
     yc2 = f1(j)
     yc3 = f1(j+1)
     call get_coeffs(yc1,yc2,yc3,-dx,dx,af,bf,cf)
     f2(j) = bf
     f3(j) = cf
  enddo
  f3(nr) = 0.0_fp
  f2(nr) = 2.0_fp*(f1(nr)-f1(nr-1))/dx
  ! f1,f2,f3-finished

  foo(:,2) = f2
  foo(:,3) = f3

end subroutine spln1d

!-----------------------
subroutine splnrz(foo)

  use iso_c_binding, only: fp => c_double, c_int
  use lscsq_mod, only: nx, nz, lh_inp

  IMPLICIT NONE

  INTEGER :: ispln
  INTEGER :: l,ld,lmax,lmax2,k,j,jm,jmm,jp,jpp,kl,kr
  INTEGER :: km,kmm,kp,kpp,ier,jd,lm,kmid,jsit,kn,nd,jl,jr
!============
  real(fp) :: fmax
  real(fp) :: dx,dz,dum,dum1,ddf,err0
  real(fp) :: yc1,yc2,yc3,af,bf,cf
  real(fp), dimension(nx,nz,9), intent(inout) :: foo

  real(fp), dimension(:,:), allocatable :: f1,f2,f3,f4,f5
  real(fp), dimension(:,:), allocatable :: f6,f7,f8,f9,wdum
      
! in domain j,k  f(x,z) = f1 + f2*dx + f3*dx**2 + f4*dz + f5*dz*dx
!  + f6*dz*dx**2 + f7*dz**2 + f8*dz**2*dx + f9*dz**2*dx**2
!  gelg solves bmat( l,m)x(m) = wk2(m), m=1,ndim
!   wk2 contains x after the call
!  bmat stored as vector bmat(lm),lm = l + (m-1)*ndim
!  before calling spln load function f into f1(j,k)=f(xj,zk)
!  after call,  f2(j,k),f3(j,k) etc will be determined
      
  jl = 2
  jr = nx - 1
  kl = 1
  kr = nz
      
  lmax= kr - kl + 1
  lmax2 = lmax*lmax
      
  dx = lh_inp%dx_grid
  dz = lh_inp%dz_grid

  !   Allocate arrays
  if(allocated(f1)) deallocate(f1,f2,f3,f4,f5,f6,f7,f8,f9,wdum)
  allocate(f1(nx,nz),f2(nx,nz),f3(nx,nz))
  allocate(f4(nx,nz),f5(nx,nz),f6(nx,nz))
  allocate(f7(nx,nz),f8(nx,nz),f9(nx,nz))
  allocate(wdum(nx,nz))
      
  !   Fill in f1
  f1 = foo(1:nx,1:nz,1)
      
  jp = 2
  do k = kl,kr
     f2(1,k) = 2.0_fp*(f1(2,k)-f1(1,k))/dx
  enddo
  f3 = 0.0_fp
      
  do k = 1,nz
     do j = jl,jr
        yc1 = f1(j-1,k)
        yc2 = f1(j,k)
        yc3 = f1(j+1,k)
        call get_coeffs(yc1,yc2,yc3,-dx,dx,af,bf,cf)
        f2(j,k) = bf
        f3(j,k) = cf
     enddo
  enddo    
  ! f1,f2,f3-finished

  !   find matrix for f4, f7
      
  kp = 2
  do j = 1,nx
     f4(j,1) = 20_fp*(f1(j,2)-f1(j,1))/dz
  enddo
  f7 = 0.0_fp
      
  do j = 1,nx
     do k = 2,nz-1
        yc1 = f1(j,k-1)
        yc2 = f1(j,k)
        yc3 = f1(j,k+1)
        call get_coeffs(yc1,yc2,yc3,-dz,dz,af,bf,cf)
        f4(j,k) = bf
        f7(j,k) = cf
     enddo
  enddo
  ! f4-f7-finished

  ! find matrix for f5,f8
  kp = 2
  do j = 1,nx
     f5(j,1) = 2.0_fp*(f2(j,2)-f2(j,1))/dz
  enddo
  f8 = 0.0_fp
      
  do j = 1,nx
     do k = 2,nz-1
        yc1 = f2(j,k-1)
        yc2 = f2(j,k)
        yc3 = f2(j,k+1)
        call get_coeffs(yc1,yc2,yc3,-dz,dz,af,bf,cf)
        f5(j,k) = bf
        f8(j,k) = cf
     enddo
  enddo
  ! f5-f8-finished

  ! find matrix for f6, f9
  kp = 2
  do j = 1,nx
     f6(j,1) = 2.0_fp*(f3(j,2)-f3(j,1))/dz
  enddo
  f9 = 0.0_fp
      
  do j = 1,nx
     do k = 2,nz-1
        yc1 = f3(j,k-1)
        yc2 = f3(j,k)
        yc3 = f3(j,k+1)
        call get_coeffs(yc1,yc2,yc3,-dz,dz,af,bf,cf)
        f6(j,k) = bf
        f9(j,k) = cf
     enddo
  enddo
  ! f6-f9-finished
      
  ! replace input arrays, then deallocate resources & return
  foo(1:nx,1:nz,1) = f1
  foo(1:nx,1:nz,2) = f2
  foo(1:nx,1:nz,3) = f3
  foo(1:nx,1:nz,4) = f4
  foo(1:nx,1:nz,5) = f5
  foo(1:nx,1:nz,6) = f6
  foo(1:nx,1:nz,7) = f7
  foo(1:nx,1:nz,8) = f8
  foo(1:nx,1:nz,9) = f9

  deallocate(f1,f2,f3,f4,f5,f6,f7,f8,f9,wdum)
      
end subroutine splnrz

!-----------------------------------------------------------
subroutine get_coeffs(f1,f2,f3,dx1,dx3,af,bf,cf)
      
  use iso_c_binding, only: fp => c_double, c_int
  implicit none

  real(fp) :: f1,f2,f3,dx1,dx3
  real(fp) :: v1,v2,v3
  real(fp) :: k,af,bf,cf
       
  ! assume dx1 = -dx3
  bf = 0.5_fp*(f3-f1)/dx3
  cf = 0.5_fp*(f3-2.0_fp*f2+f1)/(dx3*dx3)
  af = f2
        
end subroutine get_coeffs


!-------------------------------------------------------------
subroutine map1d(xdum,nr,xgrid,inp1d,out1d)

  use iso_c_binding, only: fp => c_double, c_int
  implicit none

  integer :: jd
  real(fp) :: dx,dx2
  integer, intent(in) :: nr
  real(fp), intent(out) :: out1d
  real(fp), intent(in) :: xdum
  real(fp), dimension(nr), intent(in) :: xgrid
  real(fp), dimension(nr,3), intent(in) :: inp1d

  jd=minloc(abs(xgrid-xdum),1)
      
  dx = xdum  - xgrid(jd)
  dx2 = dx*dx

  out1d = inp1d(jd,1) + inp1d(jd,2)*dx + inp1d(jd,3)*dx2 
     
        
end subroutine map1d

!-------------------------------------------------------------
subroutine psi_map(xdum,zdum,psi0)

  use iso_c_binding, only: fp => c_double, c_int
  use lscsq_mod, only: lh_inp, lh_coeff
  implicit none

  ! local b field 
  !============
  integer :: jd,idum,kd
  real(fp) :: dum,dx,dz,dx2,dz2
  real(fp), intent(in) :: xdum, zdum
  real(fp), intent(out) :: psi0

  jd=minloc(abs(lh_inp%rgrid-xdum),1)
  kd=minloc(abs(lh_inp%zgrid-zdum),1)
      
  dx = xdum  - lh_inp%rgrid(jd)
  dx2 = dx*dx
  dz = zdum  - lh_inp%zgrid(kd)
  dz2 = dz*dz

  psi0 = lh_coeff%psi(jd,kd,1) + lh_coeff%psi(jd,kd,2)*dx + lh_coeff%psi(jd,kd,3)*dx2 &
       + lh_coeff%psi(jd,kd,4)*dz + lh_coeff%psi(jd,kd,5)*dx*dz + lh_coeff%psi(jd,kd,6)*dz*dx2     &
       + lh_coeff%psi(jd,kd,7)*dz2 + lh_coeff%psi(jd,kd,8)*dz2*dx + lh_coeff%psi(jd,kd,9)*dz2*dx2
        
end subroutine psi_map 

!-------------------------------------------------------------
subroutine bfieldv_comp(xdum,zdum,bcompv)

  use iso_c_binding, only: fp => c_double, c_int
  use lscsq_mod, only: lh_inp, lh_coeff
  implicit none

  ! local b field 
  !============
  integer :: jd,idum,kd
  real(fp) :: dum,dx,dz,dx2,dz2
  real(fp), intent(in) :: xdum, zdum
  real(fp), dimension(3), intent(out) :: bcompv

  jd=minloc(abs(lh_inp%rgrid-xdum),1)
  kd=minloc(abs(lh_inp%zgrid-zdum),1)
      
  dx = xdum  - lh_inp%rgrid(jd)
  dx2 = dx*dx
  dz = zdum  - lh_inp%zgrid(kd)
  dz2 = dz*dz

  bcompv(1) = lh_coeff%Br(jd,kd,1) + lh_coeff%Br(jd,kd,2)*dx + lh_coeff%Br(jd,kd,3)*dx2 &
            + lh_coeff%Br(jd,kd,4)*dz + lh_coeff%Br(jd,kd,5)*dx*dz + lh_coeff%Br(jd,kd,6)*dz*dx2     &
            + lh_coeff%Br(jd,kd,7)*dz2 + lh_coeff%Br(jd,kd,8)*dz2*dx + lh_coeff%Br(jd,kd,9)*dz2*dx2
     
  bcompv(2) = lh_coeff%Bz(jd,kd,1) + lh_coeff%Bz(jd,kd,2)*dx + lh_coeff%Bz(jd,kd,3)*dx2 &
            + lh_coeff%Bz(jd,kd,4)*dz + lh_coeff%Bz(jd,kd,5)*dx*dz + lh_coeff%Bz(jd,kd,6)*dz*dx2     &
            + lh_coeff%Bz(jd,kd,7)*dz2 + lh_coeff%Bz(jd,kd,8)*dz2*dx + lh_coeff%Bz(jd,kd,9)*dz2*dx2
     
  bcompv(3) = lh_coeff%Bp(jd,kd,1) + lh_coeff%Bp(jd,kd,2)*dx + lh_coeff%Bp(jd,kd,3)*dx2 &
            + lh_coeff%Bp(jd,kd,4)*dz + lh_coeff%Bp(jd,kd,5)*dx*dz + lh_coeff%Bp(jd,kd,6)*dz*dx2     &
            + lh_coeff%Bp(jd,kd,7)*dz2 + lh_coeff%Bp(jd,kd,8)*dz2*dx + lh_coeff%Bp(jd,kd,9)*dz2*dx2
        
end subroutine bfieldv_comp


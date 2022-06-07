subroutine lscsq_PdepCalc
  implicit none

  CALL lscsq_RfDamp
  CALL lscsq_RfHeat

end subroutine lscsq_PdepCalc
 
!     -----------------------------------------------------------------
SUBROUTINE lscsq_RfDamp
  use lscsq_mod, only: power, pray, qlsm, ivind, izind, praytot, praysum, ppwrsum
  use lscsq_mod, only: npsi,nrays,nzones,npsi,nv,nsmoo, ok_ray
  implicit none

  integer :: izn, iry, iv, ips
!     uses decrement in ray energy along each ray to compute local
!     deposition and to deposit same into appropriate velocity
!     bin, differential volume element
!     assume power(izone, iray) is the power arriving at psi surface,
!     velocity zone izind(izone, iray), ivind(izone, iray) and use
!     differences to compute local deposition

  pray(1:nv,1:npsi) = 0.0 
  praytot(1:npsi) = 0.0 
  do iry = 1, nrays
    if (ok_ray(iry).eq.1) then
       do izn = 1, nzones-1
          iv = ivind(izn, iry)
          ! If  ivind=0, the ray was stopped before reaching this izone (izn).
          ! No calculation is appropriate.
          if (iv.EQ.0) cycle    
          ips = izind(izn, iry)
          PRay(iv,ips) = PRay(iv,ips)+(power(izn,iry)-power(izn+1,iry))
       enddo
    endif
  enddo

  CALL lscsq_Smooth(Pray, nv, npsi, nsmoo, qlsm)

  do ips = 1, npsi   
     praytot(ips) = sum(pray(1:nv,ips))
  enddo
 
  praysum = sum(praytot(1:npsi))

  ppwrsum = sum(power(1,1:nrays)-power(nzones,1:nrays))
 
end SUBROUTINE lscsq_RfDamp
!     -----------------------------------------------------------------
SUBROUTINE lscsq_RfHeat
  use lscsq_mod, only: ismo, vpar, dql, dfdv, dvol, pwrnorm, pqltot, pql, pqlsum, qlsm
  use lscsq_mod, only: npsi,nv,nsmoo, dvplus
  implicit none

  integer :: ips, iv, iSMOi

  ! compute power flow resulting from quasilinear diffusion
  pqltot(1:npsi) = 0.0 
  Pql(1:nv,1:npsi) = 0.0 
  iSMOi = mod(iSMO,2) + 1

  ! The quasilinear power deposited in each velocity bin in each psi bin is:
  !     Pql = m v Dql (-df/dv) dv dVol/dpsi dpsi
  do ips = 1, npsi
     do iv = 1, nv
        Pql(iv,ips) = -Vpar(iv)*Dql(iv,ips,iSMO)*dfdv(iv,ips,iSMOi) *  &
                       dvplus(iv) * dVol(ips) * PwrNorm
!                       dvsym(iv) * dVol(ips) * PwrNorm
     enddo
     Pqltot(ips) = sum(Pql(1:nv,ips))
  enddo

  CALL lscsq_Smooth(Pql, nv, npsi, nsmoo, qlsm)
  
  Pqlsum = sum(Pqltot(1:npsi))
 
end SUBROUTINE lscsq_RfHeat
!     -----------------------------------------------------------------
SUBROUTINE lscsq_RayDamp

  use iso_c_binding, only: fp => c_double, c_int
  use lscsq_mod, only : ivind, izind, nrays, nzones, ok_ray
  use lscsq_mod, only : rFudgDmp, dlnPdsX, power, dfdv, dlnPdsK, ivzero
  use lscsq_mod, only : pi, twopi
  implicit none

  character(len=70) :: ErrMsg
!  real(fp), dimension(nzones) :: xz, yx, yql 
  real(fp), dimension(:), allocatable :: xz, yx, yql 

  integer :: i, j, iry, ThisRay
  integer :: iv
  real(fp):: dum
  real(fp):: MxdPdZ=1.0e-5_fp
!  real(fp):: MxdPdZ=1.0e-10_fp

  allocate(xz(nzones))
  allocate(yx(nzones))
  allocate(yql(nzones))

  xz(1:nzones) = 0.0
  yx(1:nzones) = 0.0
  yql(1:nzones) = 0.0

  ! compute exponential decrement, deposit into yX (Maxwellian) and yQL (quasilinear)
  ! nzones is stored in, read from   ray.dat
 
  do iry = 1, nrays
    if (ok_ray(iry).eq.1) then
     xz(1)  = 1.0_fp
     yx(1)  = 1.0_fp
     yQL(1) = 1.0_fp
     ! First do Maxwellian damping dlnPdsX
     do i=1,nzones-1
        xz(i) = REAL(i,kind=fp)
        iv = ivind(i, iry)
        ! If ivind=0, the ray was stopped before reaching this izone
        ! No calculation is appropriate.
        if (iv.EQ.0) then
           yx(i+1) = yx(i)
           cycle
        endif

        dum = 1.0_fp + dlnPdsX(i,iry)
        ! Limit the damping per zone to the value MxdPdZ
        if (dum .LT. MxdPdz) then
           write(ErrMsg,'( '' Mxw damp! iry,izn,iv,dlnPdsX :'',i4,i4,i4,1pe11.3)')      &
                           iry, i, ivind(i,iry),dlnPdsX(i,iry)
           dum = MxdPdZ
        endif
        ! Guard against amplification of wave!
        if (dum .GT. 1.0_fp+MxdPdZ) then
           write(ErrMsg,'( '' Mxw ampl!iry,izn,iv,'',''dlnPdsX     :'',i4,i4,i4, 1pe11.3)')      &
                             iry, i, ivind(i,iry),dlnPdsX(i,iry)
           dum = 1.0_fp
        endif
        yx(i+1)=yx(i)*dum
        if(dum .EQ. MxdPdZ) then
          if (i.lt.nzones-1) yx(i+1:nzones) = 0.0 
          if (i.eq.nzones-1) yx(nzones) = 0.0  
          exit
        endif
     enddo
 
     ! Next do Quasilinear damping dlnPdsK with the Kernel
     do i=1,nzones-1
        rFudgDmp(i,iry) = 1.0_fp
        xz(i) = REAL(i,kind=fp)
        iv = ivind(i,iry)
        ! If  ivind=0, the ray was stopped before reaching this izone (izn).
        ! No calculation is appropriate.
        if (iv .EQ. 0) then
           yQL(i+1) = yQL(i)
           cycle
        endif
        dum = dlnPdsK(i,iry)*dfdv(ivind(i,iry),izind(i,iry),2)

        ! But if we are at negative velocity then dlnPdsK*dfdv needs a minus sign.
        if(iv.LT.IvZero) dum = -dum
        ! Make dum a positive number, hoped to be small compared to 1
        dum = -dum
        if (dum.LT.0.01_fp) then
           rFudgDmp(i,iry) = 1.0_fp-0.5_fp*dum+dum**2/6.0_fp
        else
           rFudgDmp(i,iry) = (1.0_fp-exp(-dum))/dum
        endif
        yQL(i+1) = yQL(i)*exp(-dum)

        if (yQL(i+1).LT.1.0e-20_fp) then
          if (i.lt.nzones-1) yql(i+1:nzones) = 0.0 
          if (i.eq.nzones-1) yql(nzones) = 0.0 
          exit
        endif
     enddo
     power(1:nzones,iry) = yql(1:nzones)*power(1,iry) !  compute power along ray
    else
      power(1:nzones,iry) = 0.0  
    endif
  enddo

  deallocate(xz, yx, yql)

end SUBROUTINE lscsq_RayDamp
!                                                                      |

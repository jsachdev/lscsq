subroutine lscsq_PdepCalc
  implicit none

  CALL lscsq_RfDamp
  CALL lscsq_RfHeat

end subroutine lscsq_PdepCalc
 
!     -----------------------------------------------------------------
SUBROUTINE lscsq_RfDamp
  use lscsq_mod, only: power, pray, qlsm, praytot, praysum, ppwrsum
  use lscsq_mod, only: npsi,nrays,nzones,npsi,nv,nsmoo
  use lscsq_mod, only: lh_out
  use iso_c_binding, only : fp => c_double
  implicit none

  integer :: izn, iry, iv, ips
!     uses decrement in ray energy along each ray to compute local
!     deposition and to deposit same into appropriate velocity
!     bin, differential volume element
!     assume power(izone, iray) is the power arriving at psi surface,
!     velocity zone izind(izone, iray), ivind(izone, iray) and use
!     differences to compute local deposition

  pray(1:nv,1:npsi) = 0.0_fp
  praytot(1:npsi) = 0.0_fp 
  do iry = 1, nrays
       do izn = 1, nzones-1
          iv = lh_out%ivind(izn, iry)
          ! If  ivind=0, the ray was stopped before reaching this izone (izn).
          ! No calculation is appropriate.
          if (iv.EQ.0) cycle    
          ips = lh_out%izind(izn, iry)
          PRay(iv,ips) = PRay(iv,ips)+(power(izn,iry)-power(izn+1,iry))
       enddo
  enddo

  CALL lscsq_Smooth(Pray, nv, npsi, nsmoo, qlsm)

  praytot = sum(pray,1)
 
  praysum = sum(praytot)

  ppwrsum = sum(power(1,1:nrays)-power(nzones,1:nrays))
 
end SUBROUTINE lscsq_RfDamp
!     -----------------------------------------------------------------
SUBROUTINE lscsq_RfHeat
  use lscsq_mod, only: ismo,  dvol, pqltot, pql, pqlsum, qlsm
  use lscsq_mod, only: npsi,nv,nsmoo, dvplus
  use lscsq_mod, only: lh_const, lh_out
  use iso_c_binding, only : fp => c_double
  implicit none

  integer :: ips, iv, iSMOi

  ! compute power flow resulting from quasilinear diffusion
  pqltot(1:npsi) = 0.0_fp 
  Pql(1:nv,1:npsi) = 0.0_fp 
  iSMOi = mod(iSMO,2) + 1

  ! The quasilinear power deposited in each velocity bin in each psi bin is:
  !     Pql = m v Dql (-df/dv) dv dVol/dpsi dpsi
  do ips = 1, npsi
     do iv = 1, nv
        Pql(iv,ips) = -lh_out%vpar(iv)*lh_out%Dql(iv,ips,iSMO)*lh_out%dfdv(iv,ips,iSMOi) *  &
                       dvplus(iv) * dVol(ips) * lh_const%PwrNorm
     enddo
     Pqltot(ips) = sum(Pql(1:nv,ips))
  enddo

  CALL lscsq_Smooth(Pql, nv, npsi, nsmoo, qlsm)
  
  Pqlsum = sum(Pqltot)
 
end SUBROUTINE lscsq_RfHeat
!     -----------------------------------------------------------------
SUBROUTINE lscsq_RayDamp

  use iso_c_binding, only: fp => c_double, c_int
  use lscsq_mod, only : nrays, nzones
  use lscsq_mod, only : rFudgDmp, dlnPdsX, power, dlnPdsK, ivzero
  use lscsq_mod, only : pi, twopi
  use lscsq_mod, only: lh_out
  implicit none

  character(len=70) :: ErrMsg
  real(fp), dimension(:), allocatable :: yx, yql 

  integer :: i, j, iry, ThisRay
  integer :: iv, ip
  real(fp):: dum
  real(fp):: MxdPdZ=1.0e-5_fp

  allocate(yx(nzones))
  allocate(yql(nzones))

  yx(1:nzones) = 0.0_fp
  yql(1:nzones) = 0.0_fp

  ! First do Maxwellian damping dlnPdsX
  ! compute exponential decrement, deposit into yX (Maxwellian) and yQL (quasilinear)
 
  do iry = 1, nrays
     if (lh_out%ok_ray(iry).eq.1) then
        yx(1)  = 1.0_fp
        yQL(1) = 1.0_fp
     endif
     ! First do Maxwellian damping dlnPdsX
     do i=1,nzones-1
        iv = lh_out%ivind(i, iry)
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
                           iry, i, iv,dlnPdsX(i,iry)
           dum = MxdPdZ
        endif
        ! Guard against amplification of wave!
        if (dum .GT. 1.0_fp+MxdPdZ) then
           write(ErrMsg,'( '' Mxw ampl!iry,izn,iv,'',''dlnPdsX     :'',i4,i4,i4, 1pe11.3)')      &
                             iry, i, iv,dlnPdsX(i,iry)
           dum = 1.0_fp
        endif
        yx(i+1)=yx(i)*dum
        if(dum .EQ. MxdPdZ) then
          if (i.lt.nzones-1) yx(i+1:nzones) = 0.0_fp 
          if (i.eq.nzones-1) yx(nzones) = 0.0_fp  
          exit
        endif
     enddo
 
     ! Next do Quasilinear damping dlnPdsK with the Kernel
     do i=1,nzones-1
        rFudgDmp(i,iry) = 1.0_fp
        iv = lh_out%ivind(i,iry)
        ip = lh_out%izind(i,iry)
        ! If  ivind=0, the ray was stopped before reaching this izone (izn).
        ! No calculation is appropriate.
        if (iv .EQ. 0) then
           yQL(i+1) = yQL(i)
           cycle
        endif
        dum = dlnPdsK(i,iry)*lh_out%dfdv(iv,ip,2)

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
          if (i.lt.nzones-1) yql(i+1:nzones) = 0.0_fp 
          if (i.eq.nzones-1) yql(nzones) = 0.0_fp 
          exit
        endif
     enddo
     power(1:nzones,iry) = yql(1:nzones)*power(1,iry) !  compute power along ray
  enddo

  deallocate(yx, yql)

end SUBROUTINE lscsq_RayDamp
!                                                                      |

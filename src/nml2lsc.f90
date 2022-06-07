subroutine nml2lsc
! read namelist and populates general parameters 
! will be replaced/integrated with IDS once completed 

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: nant,fghz, ntors, npols, nrays, nGrps,    &
              npeaks, centers_ant, widths_ant, powers_ant,     &
              fghz_ant, &
              nparmax,  nparmin,  npolmin,  npolmax,           &
              centers, couplers, widths, powers, phaseDeg,     &
              DoBram,  nslices, power_inp,                     &
              DiffuJrf, PrfSpred
  use lscsq_mod, only:  HstpLH,  nstep,  npsi,   nzones,         &
              nv,    nsmoo,     nsmw,                          &
              nRampUp,    nFlat, WeghtItr,                     &
              DoXcam,                                          &
              TailTeps, TailPeps, TailNeps,                    &
              ScatKdeg, TurnNegs,                              &
              thet0, dthet, nth
 
  use lscsq_mod, only: lscsq_par
  use lscsq_mod, only: lscsq_alloc, lscsq_allocrays
  use lscsq_mod, only: lh_inp
  implicit none

  character(len=200) :: zfile='input.lhh'
  logical :: exist
  integer :: i, i1, i2
  real(fp) :: tmpdum

! read first only the number of antennae and the number of peaks/groups
! allocate arrays, so that there are not zeros hanging around
! allocations are done in a separate subroutine
! arrays are deallocated at the end, before exiting

  namelist /inpsize/                                            &
              ntors, npols, ngrps, nslices,                     &
              npsi, nzones, nv, nth, nrampup, nant

  namelist /inpvalue/                                           &
              nparmax, nparmin, npolmin, npolmax,         &
              centers, couplers, widths, powers, phaseDeg,      &
              DoBram, DiffuJrf, PrfSpred, npeaks,               &
              fghz_ant,               &
              centers_ant, widths_ant, powers_ant,power_inp

  namelist /inpexprt/                                           &
              HstpLH, nstep, nsmoo, nsmw,                       &
              nFlat, WeghtItr,                                  &
              DoXcam,                                           &
              TailTeps, TailPeps, TailNeps,                     &
              ScatKdeg, TurnNegs,                               &
              thet0, dthet 
                                                                        
                                                                        
  inquire(file=zfile, exist=exist)
  if (exist) then
     open(72, file=zfile, status="old", action="read")
     read(72, nml=inpsize)
     call lscsq_checkdim
     call lscsq_alloc
     call lscsq_allocrays
     read(72, nml=inpvalue)
     read(72, nml=inpexprt)
     close(72)
  else
     write(*,*) 'template not found'
  endif

  fghz_ant = 1.0e-9*lh_inp%freqlh
  power_inp= lh_inp%powerlh
  do i=1,nant
!     ! check that normalizations add to one
     powers_ant(1:npeaks(i),i)=powers_ant(1:npeaks(i),i)/sum(powers_ant(:,i))
   enddo
  write(*,*) 'powers_ant:',powers_ant
  write(*,*) 'power_inp:',power_inp

  do i=1,nant
     i1 = 1+(i-1)*ntors*npols*nth
     i2 = i*ntors*npols*nth
     fghz(i1:i2) = fghz_ant(i)
  enddo
 
  if (DiffuJrf .LT. 0.0) DiffuJrf = 0.00
  if (PrfSpred.LT.0.0 .or. (DiffuJrf.EQ.0.0 .and. PrfSpred.GT.0.0)) PrfSpred = 0.00
  if (PrfSpred.GT.1.0)  PrfSpred = 1.0
  if (nparmin.GE.nparmax .and. DoBram .EQ. 0) then
     tmpdum = nparmax
     nparmax = nparmin+1.0e-3_fp
     nparmin = tmpdum -1.0e-3_fp
     CALL lscsq_LSCwarn(' nparmin/max reversed')
  endif

  if ( nsmoo.GT.nv/3 .or. nsmoo.GT.nzones/3 ) then
     CALL lscsq_LSCwarn (' nsmoo set to min(nv,nzones)/3 ')
     nsmoo = min(nv, nzones) / 3
  endif
 
  if(mod(nsmoo, 2) .EQ. 0) then
     CALL lscsq_LSCwarn (' nsmoo MUST BE ODD ')
     nsmoo = abs(nsmoo-1)
  endif
  if ( nsmw .LT. nsmoo/8) CALL lscsq_LSCwarn (' nsmoo-width seems too SMALL ')
  if ( nsmw .GE. nsmoo) CALL lscsq_LSCwarn (' nsmoo-width seems too LARGE ')
  if ( WeghtItr .GT. 1.0 .or. WeghtItr .LT. 0.0 ) then
     WeghtItr = 0.5
     CALL lscsq_LSCwarn (' WeghtItr set to 0.50 ')
  endif
  if ( TailTeps .GT. 0.00 .and. TailNeps .GT. 0.00 ) then
     TailPeps = TailNeps/TailTeps
  else if (TailTeps .GT. 0.00 .and. TailPeps .GT. 0.00 ) then
     TailNeps = TailPeps*TailTeps
  else if (TailNeps .GT. 0.00 .and. TailPeps .GT. 0.00 ) then
     TailTeps = TailNeps/TailPeps
  else
     TailTeps = 0.00
     TailNeps = 0.00
     TailPeps = 0.00
  endif
  if (TailTeps.GT.0.3 .or. TailPeps.GT.0.3 .or. TailNeps.GT.0.1 ) then
     CALL lscsq_LSCwarn(' Too much fast electron tail')
     TailTeps = 0.00
     TailNeps = 0.00
     TailPeps = 0.00
  endif


! fmp - need to define types and populate

contains

subroutine lscsq_checkdim

  use lscsq_mod, only: ind_ray
  implicit none

  integer :: itor, ipol, ith, ir, iant, i1, i2

  ! The if-else-if sequence below is to maintain compatibility w.
  ! old usage npols=1,ntors=nrays, ntors not explicitly used;
  ! IF ntors&nrays not given then use NTORDIM
  if (npols.LT.1) npols = 1
  if (ntors.LT.1) ntors = 1
  if (nth.LT.1) nth = 1
  
  nrays = ntors*npols*nth*nant
  ! define index for ntor and npol for each ray
  ! This is done to enable parallelization of the code rather than incrementing
  ! the rays indicator one by one in a loop.
  ! with multi-frequency code, this array needs to be extended to index all
  ! antennae
  if (.not.allocated(ind_ray)) allocate(ind_ray(4,nrays))
  if (.not.allocated(fghz)) allocate(fghz(nrays))

  ir = 1
  do iant=1,nant
    do ith=1,nth
      do itor = 1,ntors
        do ipol = 1,npols
           ind_ray(1,ir) = itor 
           ind_ray(2,ir) = ipol 
           ind_ray(3,ir) = ith 
           ind_ray(4,ir) = iant
           ir = ir+1
        enddo
      enddo
    enddo
  enddo

end subroutine lscsq_checkdim





end subroutine nml2lsc



subroutine nml2lsc
! read namelist and populates general parameters 
! will be replaced/integrated with IDS once completed 

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: nant,fghz, ntors, npols, nrays, nGrps,    &
              npeaks, centers_ant, widths_ant, powers_ant,     &
              fghz_ant, &
              nparmax,  nparmin,  npolmin,  npolmax,           &
              centers, widths, powers,     &
              power_inp,                     &
              DiffuJrf, PrfSpred
  use lscsq_mod, only:  HstpLH,  nstep,  npsi,   nzones,         &
              nv,    nsmoo,     nsmw,                          &
              nRampUp,    nFlat, WeghtItr,                     &
              TailTeps, TailPeps, TailNeps,                    &
              TurnNegs,                              &
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
              ntors, npols, ngrps,                      &
              npsi, nzones, nv, nth, nrampup, nant

  namelist /inpvalue/                                           &
              nparmax, nparmin, npolmin, npolmax,         &
              centers, widths, powers,       &
              DiffuJrf, PrfSpred, npeaks,               &
              fghz_ant,               &
              centers_ant, widths_ant, powers_ant,power_inp

  namelist /inpexprt/                                           &
              HstpLH, nstep, nsmoo, nsmw,                       &
              nFlat, WeghtItr,                                  &
              TailTeps, TailPeps, TailNeps,                     &
              TurnNegs,                               &
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
     ! check that normalizations add to one
     powers_ant(1:npeaks(i),i)=powers_ant(1:npeaks(i),i)/sum(powers_ant(:,i))
   enddo

  do i=1,nant
     i1 = 1+(i-1)*ntors*npols*nth
     i2 = i*ntors*npols*nth
     fghz(i1:i2) = fghz_ant(i)
  enddo
 
  if (DiffuJrf .LT. 0.0_fp) DiffuJrf = 0.0_fp
  if (PrfSpred.LT.0.0_fp .or. (DiffuJrf.EQ.0.0_fp .and. PrfSpred.GT.0.0_fp)) PrfSpred = 0.0_fp
  if (PrfSpred.GT.1.0_fp)  PrfSpred = 1.0_fp
  if (nparmin.GE.nparmax) then
     tmpdum = nparmax
     nparmax = nparmin+1.0e-3_fp
     nparmin = tmpdum -1.0e-3_fp
     CALL lscsq_LSCwarn(' nparmin/max reversed')
  endif

  ! check that nth is an odd number  
  if (2*(nth/2).EQ.nth) then
     CALL lscsq_LSCwarn (' require ODD nth ')
     nth = nth-1
  endif

  ! check that nv is an odd number  
  if (2*(nv/2).EQ.nv) then
     CALL lscsq_LSCwarn (' require ODD nv ')
     nv = nv-1
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
  if ( WeghtItr .GT. 1.0_fp .or. WeghtItr .LT. 0.0_fp ) then
     WeghtItr = 0.5_fp
     CALL lscsq_LSCwarn (' WeghtItr set to 0.50 ')
  endif
  if ( TailTeps .GT. 0.0_fp .and. TailNeps .GT. 0.0_fp ) then
     TailPeps = TailNeps/TailTeps
  else if (TailTeps .GT. 0.0_fp .and. TailPeps .GT. 0.0_fp ) then
     TailNeps = TailPeps*TailTeps
  else if (TailNeps .GT. 0.0_fp .and. TailPeps .GT. 0.0_fp ) then
     TailTeps = TailNeps/TailPeps
  else
     TailTeps = 0.0_fp
     TailNeps = 0.0_fp
     TailPeps = 0.0_fp
  endif
  if (TailTeps.GT.0.3_fp .or. TailPeps.GT.0.3_fp .or. TailNeps.GT.0.1_fp ) then
     CALL lscsq_LSCwarn(' Too much fast electron tail')
     TailTeps = 0.0_fp
     TailNeps = 0.0_fp
     TailPeps = 0.0_fp
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



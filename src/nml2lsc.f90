subroutine nml2lsc
! read namelist and populates general parameters 
! will be replaced/integrated with IDS once completed 

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: nant,fghz, ntors, npols, nrays, nGrps,    &
              npeaks, centers_ant, widths_ant, powers_ant,     &
              fghz_ant, &
              npolmin,  npolmax,           &
              centers, widths, powers,     &
              power_inp, Rant, Zant, Hant,                     &
              DiffuJrf, PrfSpred
  use lscsq_mod, only:  HstpLH,  nstep,  npsi,   nzones,         &
              nv,    nsmoo,     nsmw,                          &
              nRampUp,    nFlat, WeghtItr,                     &
              TailTeps, TailPeps, TailNeps,                    &
              TurnNegs,                              &
              nth
 
  use lscsq_mod, only: lscsq_par
  use lscsq_mod, only: lscsq_alloc, lscsq_allocrays
  use lscsq_mod, only: lh_inp
  implicit none

  character(len=200) :: zfile='input.lhh'
  logical :: exist
  integer :: i, j, i1, i2
  real(fp) :: tmpdum

! read first only the number of antennae and the number of peaks/groups
! allocate arrays, so that there are not zeros hanging around
! allocations are done in a separate subroutine
! arrays are deallocated at the end, before exiting

  namelist /inpsize/                                            &
              npols, npsi, nzones, nv, nrampup, nant

  namelist /inpvalue/                                           &
              ntors, nth, npolmin, npolmax,         &
              centers, widths, powers,       &
              DiffuJrf, PrfSpred, npeaks,               &
              fghz_ant,Rant,Zant,Hant,               &
              centers_ant, widths_ant, powers_ant,power_inp

  namelist /inpexprt/                                           &
              HstpLH, nstep, nsmoo, nsmw,                       &
              nFlat, WeghtItr,                                  &
              TailTeps, TailPeps, TailNeps,                     &
              TurnNegs                               
                                                                        
                                                                        
  inquire(file=zfile, exist=exist)
  if (exist) then
     open(72, file=zfile, status="old", action="read")
     read(72, nml=inpsize)
     call alloc_nant
     call lscsq_alloc
     read(72, nml=inpvalue)
     call lscsq_checkdim
!     call alloc_npeaks
     call lscsq_allocrays
     read(72, nml=inpexprt)
     close(72)
  else
     write(*,*) 'template not found'
  endif

!  fghz_ant = 1.0e-9*lh_inp%freqlh
  power_inp= lh_inp%powerlh
  do i=1,nant
     ! check that normalizations add to one
     powers_ant(1:npeaks(i),i)=powers_ant(1:npeaks(i),i)/sum(powers_ant(:,i))
  enddo

!  if (DiffuJrf .LT. 0.0_fp) DiffuJrf = 0.0_fp
!  if (PrfSpred.LT.0.0_fp .or. (DiffuJrf.EQ.0.0_fp .and. PrfSpred.GT.0.0_fp)) PrfSpred = 0.0_fp
!  if (PrfSpred.GT.1.0_fp)  PrfSpred = 1.0_fp

  if ( nsmoo.GT.nv/3 .or. nsmoo.GT.nzones/3 ) then
     CALL lscsq_LSCwarn (' nsmoo set to min(nv,nzones)/3 ')
     nsmoo = min(nv, nzones) / 3
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

  integer :: itor, ipol, ith, nt, ir, iant, ipk, npk, i1, i2

  ! The if-else-if sequence below is to maintain compatibility w.
  ! old usage npols=1,ntors=nrays, ntors not explicitly used;
  ! IF ntors&nrays not given then use NTORDIM
!  if (npols.LT.1) npols = 1
!  if (nth.LT.1) nth = 1
  
  nrays = 0
  do iant=1,nant
     nrays = nrays+npols*nth(iant)*sum(ntors(:,iant))
  enddo
  ! define index for ntor and npol for each ray
  ! This is done to enable parallelization of the code rather than incrementing
  ! the rays indicator one by one in a loop.
  if (.not.allocated(ind_ray)) allocate(ind_ray(5,nrays))
  if (.not.allocated(fghz)) allocate(fghz(nrays))

  fghz_ant = 1.0e-9*lh_inp%freqlh
  ir = 1
  do iant=1,nant
    npk = npeaks(iant)
    nt = nth(iant)
    do ith=1,nt
      do ipk=1,npk
        do ipol = 1,npols
          do itor = 1,ntors(ipk,iant)
             ind_ray(1,ir) = itor 
             ind_ray(2,ir) = ipol 
             ind_ray(3,ir) = ipk 
             ind_ray(4,ir) = ith 
             ind_ray(5,ir) = iant
             fghz(ir) = fghz_ant(iant)
             ir = ir+1
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine lscsq_checkdim

!-----------------------------------------------

subroutine alloc_nant

  use lscsq_mod, only: Rant, Zant, Hant, nant, ntors, npols, nth
  use lscsq_mod, only: fghz_ant, power_inp, npeaks
  use lscsq_mod, only: powers_ant, centers_ant, widths_ant
  use lscsq_mod, only: ntor_ant, spec_ant, npol_ant
  use lscsq_mod, only: ntor, npol, Spec

  implicit none

  integer :: np, nt


  if(.not.allocated(Rant)) allocate(Rant(nant))
  Rant = 0.0_fp
  if(.not.allocated(Zant)) allocate(Zant(nant))
  Zant = 0.0_fp
  if(.not.allocated(Hant)) allocate(Hant(nant))
  Hant = 0.0_fp
  if(.not.allocated(fghz_ant)) allocate(fghz_ant(nant))
  fghz_ant = 0.0_fp
  if(.not.allocated(power_inp)) allocate(power_inp(nant))
  power_inp = 0.0_fp
  if(.not.allocated(nth)) allocate(nth(nant))
  nth = 0
  if(.not.allocated(npeaks)) allocate(npeaks(nant))
  npeaks = 0

  np = 10 !  placeholder    
  if(.not.allocated(powers_ant)) allocate(powers_ant(np,nant))
  powers_ant = 0.0_fp
  if(.not.allocated(centers_ant)) allocate(centers_ant(np,nant))
  centers_ant = 0.0_fp
  if(.not.allocated(widths_ant)) allocate(widths_ant(np,nant))
  widths_ant = 0.0_fp
  if (.not.allocated(ntors)) allocate(ntors(np,nant)) 
  ntors = 0

  if(.not.allocated(npol)) allocate(npol(npols))
  npol(1:npols) = 0.0_fp

end subroutine alloc_nant

!------------------------------------------
subroutine alloc_npeaks

  use pl_types 
  implicit none

  type(ray_init), dimension(:,:), allocatable :: spec_ini  

  integer :: i, j, nt, np

  np = maxval(npeaks)
  if(.not.allocated(spec_ini)) allocate(spec_ini(np,nant))

  do i=1,nant ! test only 1 antenna for now
     do j=1,np
        spec_ini(j,i)%spec = 0.0_fp
        spec_ini(j,i)%ntor = 0.0_fp
        spec_ini(j,i)%ntors = 0
     enddo
!     do j=1,npeaks(i)
!        nt = ntors(j,i)
!        allocate(spec_ini(j,i)%spec(nt))
!        spec_ini(j,i)%spec = 0.0_fp
!        allocate(spec_ini(j,i)%ntor(nt))
!        spec_ini(j,i)%ntor = 0.0_fp
!        spec_ini(j,i)%ntors = ntors(j,i)
!     enddo
  enddo

end subroutine alloc_npeaks
!------------------------------------------






end subroutine nml2lsc



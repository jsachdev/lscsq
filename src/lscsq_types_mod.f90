module PL_TYPES
  use iso_c_binding, only : fp => c_double 

    type pl_data
      Real(fp) :: fghz
      Real(fp) :: Eper
      Real(fp) :: Epar
      Real(fp) :: Exy 
      Real(fp) :: omega
      Real(fp) :: woc
      Real(fp) :: woc2
      Real(fp) :: woc4
      Real(fp) :: epsq
      Real(fp) :: ecyc2
      Real(fp) :: ecyc
      Real(fp) :: wdddw
      Real(fp) :: dtdv
      Real(fp) :: pe2min
    end type pl_data
  
    type ray_init
       integer :: ntors
       integer :: npols
       real(fp), dimension(:), allocatable :: spec
       real(fp), dimension(:), allocatable :: ntor
       real(fp), dimension(:), allocatable :: npol
    end type ray_init

    type storeAry
      Real(fp):: rofray
      Real(fp):: zofray
      Real(fp):: Pofray
      Real(fp):: nparry
      Real(fp):: nperry
      Real(fp):: rtpsry
      Real(fp):: powrry
      Real(fp):: TimeRy
      Real(fp):: DistRy
      Real(fp):: DetrRy
      Real(fp):: NeofRy
      Real(fp):: BthRay   
      Real(fp):: dlnPdsK
      Real(fp):: dlnPdsX   
      Real(fp):: ezsq
      integer:: ivind
      integer:: izind
    end type storeAry

end module PL_TYPES

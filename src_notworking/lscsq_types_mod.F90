module PL_TYPES
  use iso_c_binding, only : fp => c_double 

  !$acc routine 
    type accum_data
    Real(fp) :: rayqt1
    Real(fp) :: accum1
    Real(fp) :: rayqt2
    Real(fp) :: accum2
    Real(fp) :: rayqt3
    Real(fp) :: accum3
    Real(fp) :: rayqt4
    Real(fp) :: accum4
    
      
      
    end type accum_data

    type old_data

    Real(fp) ::  oldRBphi
    

      Real(fp) :: RzindOld
      Real(fp) :: tOld 
      Real(fp) :: Rold
      Real(fp) :: Zold
      Real(fp) :: psiold       
      Real(fp) :: Brold
      Real(fp) :: Bzold
      Real(fp) :: omcold 
      Real(fp) :: RBphio
      Real(fp) :: pe2old
      Real(fp) :: Teeold
      Real(fp) :: pi2old
      Real(fp) :: Aioold
      Real(fp) :: Aelold
      Real(fp) :: tleave
      Real(fp) :: tenter
      integer :: izindold
      integer :: iread
      integer :: NinAc
      logical :: first_call
    end  type old_data

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


    type vecs
          Real(fp)  :: pe2vec
          Real(fp)  :: pi2vec
          Real(fp)  ::  teKev 
          Real(fp) ::  aiovec
          Real(fp)  :: aelvec 
          Real(fp)  ::  RBpvec 
          Real(fp)  ::  psivec
          ! Real(fp), allocatable, dimension(:,:) ::  psigrd
          ! Real(fp), dimension(:), allocatable :: vpar
          ! Real(fp),  :: dvol
    end type vecs

    ! type parameters
    !     ! Real(fp) :: qe_eV
    !     ! Real(fp) :: me_Kg
    !     Real(fp) :: pi
    !     !Real(fp) :: omcFac 
    !     Real(fp) :: rmin
    !     Real(fp) :: zmin
    !     Real(fp) :: eps0
    !     Real(fp) :: cEparIK
    !     Real(fp) :: delpsi
    !     Real(fp) :: psimin
    !     Real(fp) :: psilim
    !     Real(fp) :: pe2fac
    !     Real(fp) :: HstpLH
    !     ! Real(fp) :: deez
    !     ! Real(fp) :: deex
    !     integer :: nzones
    !    ! integer :: npsij
    !     ! integer :: nx 
    !     ! integer :: nv
    !     ! integer :: npsi
    !     ! integer :: nz 
    !     integer :: nrays
    !     !integer :: nCoeff
    !     integer :: ivZero

    !     real(fp) :: NparDead! = 9.8_fp 
    !     real(fp) :: NperDead != 147.0_fp
    !     real(fp) :: MaxwDead != 1.0e-12_fp
    !     real(fp) :: MxdPdz   != 0.90_fp
    !     real(fp) :: SMALL    != 1.0e-30_fp
    !     real(fp) :: large!=1.0e30_fp
    !     real(fp) :: Te2Ve
    !     real(fp) :: dr
        
    ! end type parameters

end module PL_TYPES
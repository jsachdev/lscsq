
subroutine lscsq_plasma2d(o_data, p_data, r,z,psi,Br,&
                            Bz,RBphi,omc,&
                            Tee,pe2,pi2,aio,ael)
    
    use pl_types 
    use plasma1d
    use grap
    use iso_c_binding, only : fp => c_double
    use lscsq_mod, only : omcfac
   
    implicit none

   ! type(vecs), intent(inout), dimension(npsij):: all_vecs
    Real(fp), intent(in) :: r
    Real(fp), intent(in) :: z
    Real(fp), intent(inout) :: Br
    Real(fp), intent(inout)  :: Bz
    Real(fp), intent(inout)  :: RBphi
    Real(fp), intent(inout)  :: Tee
    Real(fp), intent(inout)  :: pe2
    Real(fp), intent(inout)  :: pi2
    Real(fp), intent(inout)  :: aio
    Real(fp), intent(inout) :: ael
    Real(fp), intent(inout) :: psi
    Real(fp), intent(inout) :: omc

    Real(fp) :: large=1.0e30_fp
    Real(fp) :: Te2Ve=0.0445e-4_fp
    
    integer :: NpsiJ, jmin, i
    type(old_data), intent(inout) :: o_data
    type(pl_data), intent(inout) :: p_data
    !type(parameters), intent(inout) :: local_prmtrs

    Real(fp) :: psderiv(0:2, 0:2)

    !$acc routine     
    !$acc routine(lscsq_grapLSC)
    !$acc routine(lscsq_plasma1d)
    !$acc routine(lscsq_plasma1dACC)
    !$acc routine(lscsq_grnu1d)
    !$acc routine(lscsq_LSCwarn)


    if (o_data%iread == 1 ) go to 100
        o_data%iread = 1
        o_data%Rold  = 0.0_fp 
        o_data%Zold  = 0.0_fp 
        o_data%Psiold = 0.0_fp  
        o_data%Brold  = 0.0_fp 
        o_data%Bzold  = 0.0_fp 
        o_data%omcold = 0.0_fp 

    100  Continue

    if (r == o_data%Rold .and. z == o_data%Zold) then
        Psi    = o_data%psiold
        Br     = o_data%Brold
        Bz     = o_data%Bzold
        omc = o_data%omcold
    !  write(*,*) 'old'

    else
        o_data%Rold = r
        o_data%Zold = z

        call lscsq_grapLSC(r, z, psi, psderiv)
        Br = - psderiv(0, 1) / r
        Bz = psderiv(1, 0) / r

    
    Endif

    !call lscsq_plasma1dACC (psi,RBphi,Tee,pe2,pi2,aio,ael,o_data, p_data, all_vecs)
    call lscsq_plasma1d (psi, RBphi,Tee,pe2,pi2,aio,ael, o_data, p_data)
    
    !OmcFac  = 1.0e-9_fp*qe_eV/me_Kg/TWOPI 
    omc = OmcFac * sqrt( Br**2 + Bz**2 + (RBphi/r)**2 )
    Ael = Ael / omc**4
    !Te= Tee

    o_data%Brold  = Br
    o_data%Bzold  = Bz
    o_data%omcold = omc

end subroutine lscsq_plasma2d

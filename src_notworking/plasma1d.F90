module plasma1d
    contains

      subroutine lscsq_plasma1d (psi, RBphi,Tee,pe2,pi2,aio,ael, o_data, p_data )
        use pl_types
        use grap
        use Tracing
        use iso_c_binding, only : fp => c_double
        use lscsq_mod, only : pe2vec, pi2vec, teKev, aiovec, aelvec 
        !use lscsq_mod, only : psimin, psilim
        use lscsq_mod, only : RBpvec, psivec!, NpsiJ
        
        implicit none
        type(old_data), optional, intent(inout) :: o_data
        type(pl_data), optional, intent(inout) :: p_data
        Real(fp), intent(inout) :: psi 
        Real(fp), intent(out) :: RBphi
        Real(fp), intent(out) :: Tee  
        Real(fp), intent(out) :: pe2  
        Real(fp), intent(out) :: pi2  
        Real(fp), intent(out) :: aio  
        Real(fp), intent(out) :: ael  
        integer :: jmin, i, NpsiJ
        Real(fp) :: psl, psimin, psilim
        Real(fp) :: TeKevMin, min_val
        !Real(fp), dimension(NpsiJ) :: vec
        character(len=80) :: mes 

        !Real(fp), intent(inout), optional :: psimin1, psilim1
        !$acc routine seq
        !$acc routine(my_minloc)
        !$acc routine(lscsq_grnu1d)
        !acc routine(lscsq_LSCwarn)

        if(present(o_data)) then
          if(o_data%first_call) then
            o_data%first_call = .false.
            o_data%Psiold = 0.0
            o_data%RBphio = 0.0 
            o_data%pe2old = 0.0 
            o_data%Teeold = 0.0  
            o_data%pi2old = 0.0 
            o_data%Aioold = 0.0  
            o_data%Aelold = 0.0  
          endif
          if(psi == o_data%psiold) then
            RBphi  = o_data%RBphio
            pe2    = o_data%pe2old
            Tee    = o_data%Teeold
            pi2    = o_data%pi2old
            Aio    = o_data%Aioold
            Ael    = o_data%Aelold
         !   write(*,*) 'old1'
            return
          endif
        endif

        NpsiJ = size(psivec)

  
        !call my_minloc(abs(psivec-psi), npsij, jmin)
        min_val = 1e+30
        jmin = 1
        do i=1, npsij
            if (abs(psivec(i)-psi) < min_val) then
              jmin = I
              min_val = abs(psivec(i)-psi)
            endif
        enddo
  
        CALL lscsq_grnu1d(NpsiJ,PsiVec,RBpVec,jmin,psi,RBphi)
  
        ! enforce bounds on psi !!!!!!!!!!!!!!

        psilim = psivec(npsij)
        psimin = psivec(1)

        ! if (present(psimin1)) psimin1=psimin
        ! if (present(psilim1)) psilim1=psilim


        psi =  max (psimin, psi)
        psl =  min (psilim, psi)
      !     stop forcing bounds on psi on the maximum
      !     enforce bounds on psi !!!!!!!!!!!!!!
      !
        CALL lscsq_grnu1d(NpsiJ, PsiVec, pe2Vec, jmin, psl, pe2)
      !----------------------------------------------------
      ! DMC bugfix -- handle pe2 minimum value not at bdy
      !   (prevent LSC from executing premature termination of ray following)
  
        if(present(p_data)) then
          if ((psl.lt.PsiVec(NpsiJ-1)).AND.(pe2.le.p_data%pe2min)) &
          p_data%pe2min=pe2-1.0d-6*abs(pe2)
        endif
      
      !----------------------------------------------------
  
        CALL lscsq_grnu1d(NpsiJ, PsiVec, Tekev, jmin, psl, Tee)
        CALL lscsq_grnu1d(NpsiJ, PsiVec, pi2Vec, jmin, psl, pi2)
        CALL lscsq_grnu1d(NpsiJ, PsiVec, AioVec, jmin, psl, Aio)
        CALL lscsq_grnu1d(NpsiJ, PsiVec, AelVec, jmin, psl, Ael)
  
        if(present(o_data)) then
          o_data%Psiold = Psi
          o_data%RBphio = RBphi
          o_data%pe2old = pe2
          o_data%Teeold = Tee
          o_data%pi2old = pi2
          o_data%Aioold = Aio
          o_data%Aelold = Ael
          TeKevMin = 0.01_fp
        endif
        
        if(Tee.le.0.0) then

          write(*,*) ' negative temperature correction ', Tee
          !mes = ' negative temperature correction ' 
          !call lscsq_LSCwarn(mes)
          Tee = TeKevMin
        endif
      end subroutine lscsq_plasma1d

      ! subroutine set_consts(local_prmtrs)
      !   use pl_types
      !   use iso_c_binding, only : fp => c_double
      !   use lscsq_mod, only : qe_eV, me_Kg,pi,OmcFac,rmin,zmin,eps0,cEparIK,delpsi,psimin, psilim
      !   use lscsq_mod, only :  pe2fac,HstpLH,nzones, NpsiJ,nx,nv,npsi,nz,nrays,nCoeff,ivZero
      !   !use lscsq_mod, only :  pe2vec,pi2vec,teKev,aiovec,aelvec,RBpvec , psivec
        
        

      !   implicit none
      !   type(parameters), intent(inout) :: local_prmtrs
      !  ! type(vecs), dimension( NpsiJ), intent(inout):: all_vecs

 
      !   ! local_prmtrs%qe_eV = qe_eV
      !   ! local_prmtrs%me_Kg  = me_Kg
      !   local_prmtrs%pi   = pi
      !   local_prmtrs%OmcFac  = OmcFac
      !   local_prmtrs%rmin = rmin
      !   local_prmtrs%zmin = zmin
      !   local_prmtrs%eps0 = eps0
      !   local_prmtrs%cEparIK = cEparIK
      !   local_prmtrs%delpsi = delpsi
      !   local_prmtrs%psimin = psimin
      !   local_prmtrs%psilim = psilim
      !   local_prmtrs%pe2fac = pe2fac
      !   local_prmtrs%HstpLH = HstpLH
      !   ! local_prmtrs%deex =    5.0112905833333337E-003_fp  
      !   ! local_prmtrs%deez =  7.5430380120481933E-003_fp
      !   local_prmtrs%nzones = nzones
      !   !local_prmtrs%NpsiJ =  NpsiJ
      !   ! local_prmtrs%nx = nx
      !   ! local_prmtrs%nv = nv
      !   ! local_prmtrs%npsi = npsi
      !   ! local_prmtrs%nz = nz
      !   local_prmtrs%nrays = nrays
      !   !local_prmtrs%nCoeff = nCoeff
      !   local_prmtrs%ivZero = ivZero


        
      !   ! local_prmtrs%NparDead = 9.8_fp 
      !   ! local_prmtrs%NperDead = 147.0_fp
      !   ! local_prmtrs%MaxwDead = 1.0e-12_fp
      !   ! local_prmtrs%MxdPdz   = 0.90_fp
      !   ! local_prmtrs%SMALL    = 1.0e-30_fp
      !   ! local_prmtrs%large = 1.0e30_fp
      !   ! local_prmtrs%Te2Ve = 0.0445e-4_fp
      !   ! local_prmtrs%dr = 1.0e-4_fp


      !   end subroutine set_consts


  !     subroutine lscsq_plasma1dACC (psi, RBphi,Tee,pe2,pi2,aio,ael, o_data, p_data, all_vecs )
  !       use pl_types
  !       use grap
  !       use Tracing
  !       use iso_c_binding, only : fp => c_double
  !      ! use lscsq_mod, only : pe2vec, pi2vec, teKev, aiovec, aelvec 
  !       !use lscsq_mod, only : psimin, psilim
  !       !use lscsq_mod, only : RBpvec, psivec, NpsiJ
  !       use lscsq_mod, only : NpsiJ
        
  !       implicit none
  !       type(old_data), optional, intent(inout) :: o_data
  !       type(pl_data), optional, intent(inout) :: p_data
  !       Real(fp), intent(inout) :: psi 
  !       Real(fp), intent(out) :: RBphi
  !       Real(fp), intent(out) :: Tee  
  !       Real(fp), intent(out) :: pe2  
  !       Real(fp), intent(out) :: pi2  
  !       Real(fp), intent(out) :: aio  
  !       Real(fp), intent(out) :: ael  
  !       integer :: jmin, i!, NpsiJ
  !       Real(fp) :: psl, psimin, psilim
  !       Real(fp) :: TeKevMin, min_val
  !       Real(fp), dimension(NpsiJ) :: vec
  !       character(len=80) :: mes 

  !       type(vecs), dimension(:), intent(inout):: all_vecs

  !       !Real(fp), intent(inout), optional :: psimin1, psilim1

  !   !$acc routine seq
  !   !$acc routine(my_minloc)
  !   !$acc routine(lscsq_grnu1d)
  !   !$acc routine(lscsq_LSCwarn)

  !       ! do i=1,NpsiJ
  !       !   if(all_vecs(I)%pe2vec .ne. pe2vec(i)) write(*,*) 'peVec', i
  !       !   if(all_vecs(I)%pi2vec .ne. pi2vec(i)) write(*,*) 'piVec', i
  !       !   if(all_vecs(I)%teKev .ne. teKev(i)) write(*,*) 'teKev', i
  !       !   if(all_vecs(I)%aiovec .ne. aiovec(i)) write(*,*) 'aiovec', i
  !       !   if(all_vecs(I)%aelvec .ne. aelvec(i)) write(*,*) 'aelvec', i
  !       !   if(all_vecs(I)%RBpvec .ne. RBpvec(i)) write(*,*) 'RBpvec', i
  !       !   if(all_vecs(I)%psivec .ne. psivec(i)) write(*,*) 'psivec', i
  !       ! enddo
       

  !   !if(present(o_data)) then
  !     if(o_data%first_call) then
  !       o_data%first_call = .false.
  !       o_data%Psiold = 0.0
  !       o_data%RBphio = 0.0 
  !       o_data%pe2old = 0.0 
  !       o_data%Teeold = 0.0  
  !       o_data%pi2old = 0.0 
  !       o_data%Aioold = 0.0  
  !       o_data%Aelold = 0.0  
  !     endif
  !     if(psi == o_data%psiold) then
  !       RBphi  = o_data%RBphio
  !       pe2    = o_data%pe2old
  !       Tee    = o_data%Teeold
  !       pi2    = o_data%pi2old
  !       Aio    = o_data%Aioold
  !       Ael    = o_data%Aelold
  !    !   write(*,*) 'old1'
  !       return
  !     endif
  !  ! endif

  !   !NpsiJ = size(psivec)


  !   !call my_minloc(abs(psivec-psi), npsij, jmin)
  !   vec = abs(all_vecs(:)%psivec-psi)
  !   min_val = 1e+30
  !   jmin = 1
  !   do i=1, npsij
  !       if (vec(i) < min_val) then
  !         jmin = I
  !         min_val = vec(i)
  !       endif
  !   enddo

  !   CALL lscsq_grnu1d(NpsiJ,all_vecs(:)%PsiVec,all_vecs(:)%RBpVec,jmin,psi,RBphi)

  !   ! enforce bounds on psi !!!!!!!!!!!!!!

  !   psilim = all_vecs(npsij)%psivec
  !   psimin = all_vecs(1)%psivec

  !   ! if (present(psimin1)) psimin1=psimin
  !   ! if (present(psilim1)) psilim1=psilim


  !   psi =  max (psimin, psi)
  !   psl =  min (psilim, psi)
  ! !     stop forcing bounds on psi on the maximum
  ! !     enforce bounds on psi !!!!!!!!!!!!!!
  ! !
  !   CALL lscsq_grnu1d(NpsiJ, all_vecs(:)%PsiVec, all_vecs(:)%pe2Vec, jmin, psl, pe2)
  ! !----------------------------------------------------
  ! ! DMC bugfix -- handle pe2 minimum value not at bdy
  ! !   (prevent LSC from executing premature termination of ray following)

  !   !if(present(p_data)) then
  !     if ((psl.lt.all_vecs(NpsiJ-1)%PsiVec).AND.(pe2.le.p_data%pe2min)) &
  !     p_data%pe2min=pe2-1.0d-6*abs(pe2)
  !   !endif
  
  ! !----------------------------------------------------

  !   CALL lscsq_grnu1d(NpsiJ, all_vecs(:)%PsiVec, all_vecs(:)%Tekev, jmin, psl, Tee)
  !   CALL lscsq_grnu1d(NpsiJ, all_vecs(:)%PsiVec, all_vecs(:)%pi2Vec, jmin, psl, pi2)
  !   CALL lscsq_grnu1d(NpsiJ, all_vecs(:)%PsiVec, all_vecs(:)%AioVec, jmin, psl, Aio)
  !   CALL lscsq_grnu1d(NpsiJ, all_vecs(:)%PsiVec, all_vecs(:)%AelVec, jmin, psl, Ael)

  !   !if(present(o_data)) then
  !     o_data%Psiold = Psi
  !     o_data%RBphio = RBphi
  !     o_data%pe2old = pe2
  !     o_data%Teeold = Tee
  !     o_data%pi2old = pi2
  !     o_data%Aioold = Aio
  !     o_data%Aelold = Ael
  !     TeKevMin = 0.01_fp
  !  ! endif
    
  !   if(Tee.le.0.0) then

  !     !write(mes, '(a, e14.7)') ' negative temperature correction ', Tee
  !     mes = ' negative temperature correction ' 
  !    ! call lscsq_LSCwarn(mes)
  !     Tee = TeKevMin
  !   endif
  ! end subroutine lscsq_plasma1dACC

  !     subroutine set_vecs(all_vecs)
  !       use pl_types
  !       use iso_c_binding, only : fp => c_double
 
  !       use lscsq_mod, only :   NpsiJ
  !       use lscsq_mod, only :  pe2vec,pi2vec,teKev,aiovec,aelvec,RBpvec , psivec
        
        

  !       implicit none
  !       integer :: i
  !       !type(const), intent(inout) :: local_prmtrs
  !       type(vecs), dimension( NpsiJ), intent(inout) :: all_vecs

       
  !       do i=1, NpsiJ
  !           all_vecs(i)%pe2vec =pe2vec(i)
  !           all_vecs(i)%pi2vec =pi2vec(i)
  !           all_vecs(i)%teKev =teKev(i)
  !           all_vecs(i)%aiovec = aiovec(i)
  !           all_vecs(i)%aelvec  = aelvec(i)
  !           all_vecs(i)%RBpvec  = RBpvec(i)
  !           all_vecs(i)%psivec = psivec(i)
            
  !       enddo

  !       end subroutine set_vecs
 
end module plasma1d
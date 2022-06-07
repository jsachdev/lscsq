
module Integrator
    contains
SUBROUTINE RungeInit( yN1, yN2, yN3, yN4, yN5, yN6, yN7, yN8,&
                      yN1_1, yN2_1, yN3_1, yN4_1, yN5_1, yN6_1, yN7_1, yN8_1, &
                      fN1, fN2, fN3, fN4, fN5, fN6, fN7,&
                      fN1_1, fN2_1, fN3_1, fN4_1, fN5_1, fN6_1, fN7_1, &
                      arrys, lstop, izone, p_data, a_data, o_data)!, local_prmtrs)!, all_vecs)
                      ! , all_const, all_vecs, &
                      ! local_psigrd, local_dvol, local_vpar)

    use iso_c_binding, only : fp => c_double 
    use pl_types 
    use RayStore
    use rTrace
    use lscsq_mod, only : nzones, HstpLH
    
    implicit NONE
    type(pl_data), intent(inout) :: p_data
    type(accum_data), intent(inout) :: a_data
    type(old_data), intent(inout) :: o_data
    !type(parameters), intent(inout) :: local_prmtrs

    integer, intent(OUT) :: lstop, izone
    integer :: BoundsEr

    Real(fp), intent(inout) :: yN1, yN2, yN3, yN4, yN5, yN6, yN7, yN8
    Real(fp), intent(OUT) ::  yN1_1, yN2_1, yN3_1, yN4_1, yN5_1, yN6_1, yN7_1, yN8_1, &
                            fN1_1, fN2_1, fN3_1, fN4_1, fN5_1, fN6_1, fN7_1
    ! ,&
    ! fN1, fN2, fN3, fN4, fN5, fN6, fN7, yN1_1, yN2_1, yN3_1, yN4_1, yN5_1, yN6_1, yN7_1, yN8_1, fN1_1, fN2_1, fN3_1, fN4_1, fN5_1, fN6_1, fN7_1
    Real(fp) :: y1, y2, y3, y4, y5, y6, y7, y8, F1, F2, F3, F4, F5, F6, F7
    Real(fp), intent(inout) :: fN1, fN2, fN3, fN4, fN5, fN6, fN7
    Real(fp) :: cum_sum1, cum_sum2, cum_sum3, cum_sum4, cum_sum5, cum_sum6, cum_sum7
    !integer, intent(INOUT)  :: jstart
    integer :: ival, jval
    Real(fp) :: hstplh6, hstplh2, save_min

    ! type(const), intent(inout) :: all_const
    ! type(vecs), dimension(npsij), intent(inout) :: all_vecs
    ! Real(fp), dimension(nx,nz), intent(inout) :: local_psigrd
    ! Real(fp),  dimension(nv), intent(inout) :: local_vpar
    ! Real(fp),  dimension(npsi), intent(inout) :: local_dvol
    type(storeAry), dimension(nzones), intent(inout) :: arrys
    !Real(fp), dimension(5,5) :: small_psigrd


    !$acc routine seq
    !$acc routine(lscsq_E2byPr)
    !$acc routine(lscsq_ftion)

    hstplh6 = HstpLH / 6.0_fp
    hstplh2 = HstpLH / 2.0_fp
    BoundsEr = 0

    !call get_val(yN1, yN2, ival, jval, all_const)
    !local_psigrd = local_psigrd(ival-2:ival+2, jval-2:jval+2)
    !write(*,*) 'Inside5: ', ival, jval, small_psigrd(1,1)

    !write(*,*) 'Inside6: ', yN1, fN1
    
    y1 = yN1 + hstplh2 * fN1
    y2 = yN2 + hstplh2 * fN2
    !y3 = yN3 + hstplh2 * fN3
    y4 = yN4 + hstplh2 * fN4
    y5 = yN5 + hstplh2 * fN5
    y6 = yN6 + hstplh2 * fN6
    !y7 = yN7 + hstplh2 * fN7 
    cum_sum1 =  fN1
    cum_sum2 =  fN2
    cum_sum3 =  fN3
    cum_sum4 =  fN4
    cum_sum5 =  fN5
    cum_sum6 =  fN6
    cum_sum7 =  fN7

    !write(*,*) 'stop11', y1, YN1, hstplh2 * fN1, fN1,  p_data%pe2

    BoundsEr = 0
    CALL lscsq_ftion(BoundsEr, y1, y2, y4, y5, y6, f1, f2, f3, f4, f5, f6, f7,  o_data, p_data)!, local_prmtrs)!, all_vecs)
   
    !write(*,*) 'Inside5: ', y1, yN1, f1

    !if(BoundsEr .ne. 0) then
      !  lstop = 1
      !  write(*,*) 'stop12', y1, YN1, f1, fN1,  p_data%pe2
       
        !return
   ! endif

    y1 = yN1 + hstplh2 * f1
    y2 = yN2 + hstplh2 * f2
    !y3 = yN3 + hstplh2 * f3
    y4 = yN4 + hstplh2 * f4
    y5 = yN5 + hstplh2 * f5
    y6 = yN6 + hstplh2 * f6
    !y7 = yN7 + hstplh2 * f7 
    cum_sum1 =  cum_sum1 +  2.0_fp*f1
    cum_sum2 =  cum_sum2 +  2.0_fp*f2
    cum_sum3 =  cum_sum3 +  2.0_fp*f3
    cum_sum4 =  cum_sum4 +  2.0_fp*f4
    cum_sum5 =  cum_sum5 +  2.0_fp*f5
    cum_sum6 =  cum_sum6 +  2.0_fp*f6
    cum_sum7 =  cum_sum7 +  2.0_fp*f7  
    
    ! write(*,*) 'Inside1: ', yN1, y1

    !call get_val(y1, y2, ival, jval, all_const)
    !local_psigrd = local_psigrd(ival-2:ival+2, jval-2:jval+2)
    BoundsEr = 0
    CALL lscsq_ftion(BoundsEr, y1, y2, y4, y5, y6, f1, f2, f3, f4, f5, f6, f7, o_data, p_data)!, local_prmtrs)!, all_vecs)
  
    !write(*,*) 'Inside4: ', y1, yN1, f1

    !if(BoundsEr .ne. 0) then
     !   lstop = 1
       ! write(*,*) 'stop13', y1, yN1,  f1, fN1,  p_data%pe2
       ! return
    !endif

    y1 = yN1 + HstpLH * f1
    y2 = yN2 + HstpLH * f2
   ! y3 = yN3 + HstpLH * f3
    y4 = yN4 + HstpLH * f4
    y5 = yN5 + HstpLH * f5
    y6 = yN6 + HstpLH * f6
    !y7 = yN7 + HstpLH * f7 
    cum_sum1 =  cum_sum1 +  2.0_fp*f1
    cum_sum2 =  cum_sum2 +  2.0_fp*f2
    cum_sum3 =  cum_sum3 +  2.0_fp*f3
    cum_sum4 =  cum_sum4 +  2.0_fp*f4
    cum_sum5 =  cum_sum5 +  2.0_fp*f5
    cum_sum6 =  cum_sum6 +  2.0_fp*f6
    cum_sum7 =  cum_sum7 +  2.0_fp*f7   


    BoundsEr = 0
    CALL lscsq_ftion(BoundsEr, y1, y2, y4, y5, y6, f1, f2, f3, f4, f5, f6, f7, o_data, p_data)!, local_prmtrs)!, all_vecs)
    !,all_const, all_vecs, local_psigrd)

   ! write(*,*) 'Inside3: ', y1, yN1, f1

    !if(BoundsEr .ne. 0) then
      !  lstop = 1
       ! write(*,*) 'stop14', y1, yN1, f1, fN1,  p_data%pe2
       ! return
   ! endif

    yN1_1 = yN1 + hstplh6 * (f1 + cum_sum1)
    yN2_1 = yN2 + hstplh6 * (f2 + cum_sum2)
    yN3_1 = yN3 + hstplh6 * (f3 + cum_sum3)
    yN4_1 = yN4 + hstplh6 * (f4 + cum_sum4)
    yN5_1 = yN5 + hstplh6 * (f5 + cum_sum5)
    yN6_1 = yN6 + hstplh6 * (f6 + cum_sum6)
    yN7_1 = yN7 + hstplh6 * (f7 + cum_sum7) 
    yN8_1 = yN8 + HstpLH
    
    BoundsEr = 0

    CALL lscsq_ftion(BoundsEr, yN1_1, yN2_1, yN4_1, yN5_1, yN6_1, fN1_1, fN2_1, fN3_1, fN4_1, fN5_1, fN6_1, fN7_1,  &
    o_data, p_data)!, local_prmtrs)!, all_vecs)
    !,all_const, all_vecs, local_psigrd)

    if(BoundsEr > 0) then
        lstop = 1
        !write(*,*) 'stop after Runge', yN1_1, yN1, fN1_1, fN1
       ! read(*,*)
        return 
    endif

!     !write(*,*) 'Inside1: ', y1, y2, y3, y4, y5, y6, y7, yN8_1, izone, lstop, arrys, &
!     !p_data, a_data, o_data, all_const, all_vecs(1), local_psigrd(1,1), local_dvol(1), local_vpar(1)

    call lscsq_E2byPr(yN1_1, yN2_1, yN3_1, yN4_1, yN5_1, yN6_1, yN7_1, yN8_1, izone, lstop, arrys, &
                        p_data, a_data, o_data)!, local_prmtrs)!, all_vecs)
   
                        !, all_const, all_vecs, local_psigrd, local_dvol, local_vpar) 
   ! write(*,*)     arrys(izone)%RofRay
   ! write(*,*) 'finish Runge', yN1, yN1_1, fN1_1 , f1,fN1 

END SUBROUTINE RungeInit


! ------------------------------------------------------ 
! ------------------------------------------------------ 


! SUBROUTINE Runge4(BoundsEr, fN, yN,  o_data, p_data)

!     use iso_c_binding, only : fp => c_double  
!     use lscsq_mod, only : hstplh, NEQS
!     use pl_types 
!     use rTrace
!     implicit none
!     type(pl_data), intent(inout) :: p_data
!     type(old_data), intent(inout) :: o_data

!     !integer :: i
!     integer, intent(INOUT) :: BoundsEr
!     Real(fp), intent(INOUT):: fN, yN
!     Real(fp), dimension(NEQS) :: cum_sum
!     Real(fp), dimension(NEQS) :: f, y 
!     Real(fp) :: hstplh6, hstplh2

!     !$acc routine seq
!     !$acc routine(lscsq_ftion)
!     !$acc routine(shift)

!     hstplh6 = hstplh / 6.0_fp
!     hstplh2 = hstplh / 2.0_fp

!     BoundsEr = 0

!     CALL lscsq_ftion(BoundsEr, yN(1:NEQS,1), f,  &
!         o_data, p_data)  

!     if(BoundsEr .ne. 0) then
!         BoundsEr = 1
!         return
!     endif  

!     !do i=1, NEQS
!         y = yN(i,1) + hstplh2 * f(i)
!         cum_sum(i) =  f(i)
!     !enddo

!     CALL lscsq_ftion(BoundsEr, y, f,  o_data, p_data)

!     if(BoundsEr .ne. 0) then
!         BoundsEr = 2
!         return
!     endif

!     !do i=1, NEQS
!         y = yN(i,1) + hstplh2 * f(i)
!         cum_sum(i) = cum_sum(i) + 2.0_fp*f(i)
!     !enddo

!     CALL lscsq_ftion(BoundsEr, y, f,  o_data, p_data)

!     if(BoundsEr .ne. 0) then
!         BoundsEr = 3
!         return
!     endif

!     !do i=1, NEQS
!         y = yN(i,1) + hstplh * f(i)
!         cum_sum(i) = cum_sum(i) + 2.0_fp*f
!     !enddo

!     CALL lscsq_ftion(BoundsEr, y, f,  o_data, p_data)

!     if(BoundsEr .ne. 0) then
!         BoundsEr = 4
!         return
!     endif

!     !do i=1, NEQS
!         y = yN(i,1) + hstplh6*(f(i) + cum_sum(i))
!     !enddo

!     CALL lscsq_ftion(BoundsEr, y, f,  o_data, p_data)

!     if(BoundsEr .ne. 0) then
!         BoundsEr = 5
!         return 
!     endif

!     !call shift(yN, fN, y, f)

!     return

! END SUBROUTINE Runge4

! ------------------------------------------------------
! ------------------------------------------------------



! SUBROUTINE IntegratorPC(BoundsEr, yN1, yN2, yN3, yN4, yN5, yN6, yN7, yN8,&
!     fN1, fN2, fN3, fN4, fN5, fN6, fN7,  o_data, p_data)

!     use iso_c_binding, only : fp => c_double 
!     use lscsq_mod, only : NEQS,nCoeff , NEQSP1
!     use lscsq_mod, only :  corrector_tol, nCorrectIter, HstpLH 
!     use lscsq_mod, only :  bCorrect_1, bPredict_1, &
!                         bCorrect_2, bPredict_2, &
!                         bCorrect_3, bPredict_3, &
!                         bCorrect_4, bPredict_4, &
!                         bCorrect_5, bPredict_5, bCorrect_6
!     use pl_types 
!     use rTrace
!     implicit NONE
!     type(pl_data), intent(inout) :: p_data
!     type(old_data), intent(inout):: o_data
!     integer :: j, i, stopFlag
!     integer, intent(inout)  :: BoundsEr
!     Real(fp) :: c 
!     Real(fp), intent(INOUT), dimension(ncoeff) :: yN1, yN2, yN3, yN4, yN5, yN6, yN7, yN8,&
!     fN1, fN2, fN3, fN4, fN5, fN6, fN7
!     Real(fp) :: y1, y2, y3, y4, y5, y6, y7, y8, f1, f2, f3, f4, f5, f6, f7 
!     Real(fp) :: cConst

    
!     !$acc routine seq
!     !$acc routine(lscsq_ftion)
!     !$acc routine(shift)

!     BoundsEr = 0


!         do j=1,NEQS
!             y(j) = yN(J, 1)  + &
!             HstpLH  * (fN(J, 1) * bPredict_1 + &
!             fN(J, 2) * bPredict_2 + &
!             fN(J, 3) * bPredict_3 + &
!             fN(J, 4) * bPredict_4 + &
!             fN(J, 5) * bPredict_5)

!             cConst(j)= yN(j,1)  +  &
!             (bCorrect_2 * fN(j,1)+ &
!             bCorrect_3 * fN(j,2)+ &
!             bCorrect_4 * fN(j,3)+  &
!             bCorrect_5 * fN(j,4)+ &
!             bCorrect_6 * fN(j,5)) * HstpLH  
!         enddo


!         do i =1, nCorrectIter
!             CALL lscsq_ftion(BoundsEr, y, f, o_data, p_data)
!             if(BoundsEr .ne. 0) return

!             stopFlag = 0
!             do j=1,NEQS
!                 c = cConst(j) + bCorrect_1* f(j) * HstpLH  

!                 if (abs((c-y(j))/min(y(j),c)) < corrector_tol) then
!                     stopFlag = stopFlag + 1
!                 endif
!                 y(j) = c
!             enddo
!             if (stopFlag.eq.NEQS) exit
!         enddo

!         CALL lscsq_ftion(BoundsEr, y, f,  o_data, p_data)
!         if(BoundsEr .ne. 0) return


!         call shift(yN, fN, y, f)
! END SUBROUTINE IntegratorPC

! ------------------------------------------------------
! ------------------------------------------------------


! SUBROUTINE shift(yN, fN, y, f)

!     use iso_c_binding, only : fp => c_double 
!     use lscsq_mod, only : neqs, nCoeff, NEQSP1

!     implicit NONE
!     Real(fp), intent(INOUT), dimension(NEQS) :: f, y 
!     Real(fp), intent(INOUT), dimension(NEQSP1,nCoeff) :: fN, yN
!     integer :: i, J

!     !$acc routine seq

!     do j= nCoeff, 2, -1
!         do i = 1, NEQSP1
!             fN(i, j) = fN(i, j-1)
!             yN(i, j) = yN(i, j-1)
!         enddo
!     enddo

!     do i = 1, NEQS
!         fN(i, 1) = f(i)
!         yN(i, 1) = y(i)
!     enddo
! END SUBROUTINE shift
end module Integrator
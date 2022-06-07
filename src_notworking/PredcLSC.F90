
module PredcLSC
    contains
    subroutine lscsq_PredcLSC(arrys)
      

    use iso_c_binding, only : fp => c_double
    USE OMP_LIB
#ifdef _OPENACC
    use openacc
#endif
    use pl_types 
    use Integrator
    use rayIni
    use RayStore
    use Tracing
    use rayIni
    use rTrace
    use plasma1d
    
    use lscsq_mod, only : twopi, npsij, vc, ok_ray, iError, fghz, pe2vec, nzones, nrays, nCoeff, HstpLH, nstep

    implicit none
    integer, dimension(nrays) :: izone, lstop, next_j, step
 
    integer :: RayIniErr, iray
    type(pl_data), dimension(nrays) :: p_data
    type(accum_data), dimension(nrays) :: a_data
    type(old_data), dimension(nrays) :: o_data
    type(storeAry), dimension(nzones, nrays), intent(inout) :: arrys
    !type(parameters) :: local_prmtrs
    
    Real(fp),  dimension(nrays, 7):: fN1, yN1, fN2, yN2, fN3, yN3, fN4, yN4, fN5, yN5, fN6, yN6, fN7, yN7, yN8
    
    Real(fp),  dimension (nCoeff+1):: bCorrect, bPredict
    integer :: i, j, k
    integer, dimension(nrays) :: BoundsEr!, jstart
    !integer :: iBndsErr, iscatplt
    character(2) :: tmp

#ifdef _OPENMP      
      CALL OMP_SET_NUM_THREADS(omp_get_num_procs())   
      !CALL OMP_SET_NUM_THREADS(1)   
#endif

    !routine(lscsq_RayIni)
    !$acc routine(RungeInit)
    !$acc routine(lscsq_ftion)
    !$acc routine(lscsq_BounceIt)
    !$acc routine(lscsq_e2bypr)

    !nrays = 2

    !6th degree
    bCorrect(6) = 475.0_fp / 1440.0_fp
    bCorrect(5) = 1427.0_fp / 1440.0_fp
    bCorrect(4) = -798.0_fp / 1440.0_fp
    bCorrect(3) = 482.0_fp / 1440.0_fp
    bCorrect(2) = -173.0_fp / 1440.0_fp
    bCorrect(1) = 27.0_fp / 1440.0_fp
    !5h degree  
    ! bCorrect(6) = 251.0_fp / 720.0_fp
    ! bCorrect(5) = 646.0_fp / 720.0_fp
    ! bCorrect(4) = -264.0_fp / 720.0_fp
    ! bCorrect(3) = 106.0_fp / 720.0_fp
    ! bCorrect(2) = -19.0_fp / 720.0_fp
    ! bCorrect(1) = 0.0_fp
  
    !aPredict(1) = 1.0_fp
    bPredict(5) = 1901.0_fp / 720.0_fp
    bPredict(4) = -2774.0_fp  / 720.0_fp
    bPredict(3) = 2616.0_fp  / 720.0_fp
    bPredict(2) = -1274.0_fp  / 720.0_fp
    bPredict(1) = 251.0_fp  / 720.0_fp

 !OMP PARALLEL DO shared(p_data, o_data, a_data, iZone, next_j, step)
    do iray=1, nrays
          call lscsq_PlDataInit(iray,  p_data(iray), o_data(iray), a_data(iray))
          iZone(iray) = 1
          next_j(iray) = 6
          step(iray) = 1
    enddo
 !OMP END PARALLEL DO

    RayIniErr = 0

    !write(*,*) nrays, size(lstop)
    !OMP PARALLEL DO private(RayIniErr)  shared(yN1, yN2, yN3, yN4, yN5,  yN6, yN7, yN8,  arrys, p_data, o_data, BoundsEr, lstop, ok_ray)
    do iray=1, nrays
      call lscsq_RayIni(RayIniErr,   yN1(iray,step(iray)), yN2(iray, step(iray)), &
                                     yN3(iray, step(iray)), yN4(iray, step(iray)),&
                                     yN5(iray, step(iray)), yN6(iray, step(iray)),&
                                     yN7(iray, step(iray)), yN8(iray, step(iray)),&
                                     iray, arrys(:,iray), o_data(iray), p_data(iray))
      if (RayIniErr > 0) then
              ok_ray(iray) = 0
              RayIniErr = 0
              lstop(iray) = 1    
      else
              lstop(iray) = 0
              BoundsEr(iray) = 0
      endif
     ! write(*,*) lstop(iray)
    ENDDO
 !OMP END PARALLEL DO
    
   
    if (sum(ok_ray)==0) then
      write(*,*) 'None of ray were initialized'
      iError = 1
      return
    else
      iError = 0
    endif

 
    ! write(*,*) 'Staring 1st Runge', yN1(1,1), fN1(1,1)
    ! write(*,*) 'Not ok rays: ', sum(lstop)
    
    ! iscatplt = 1
    ! BoundsEr = 0
    ! iBndsErr = 0

    !$acc data copy(arrys) create(fN1, fN2, fN3, fN4, fN5, fN6, fN7, j) copyin(BoundsEr, yN1, yN2, yN3, yN4, yN5,  yN6, yN7, yN8, step, o_data, p_data, a_data, next_j, izone, lstop, bCorrect, bPredict) 


    !$acc parallel loop 
    !OMP PARALLEL DO shared(fN1, fN2, fN3, fN4, fN5, fN6, fN7, yN1, yN2, yN4, yN5,  yN6, p_data, o_data, BoundsEr, lstop)
    do iray=1, nrays

      if (lstop(iray) == 0) then 
        
        CALL lscsq_ftion(BoundsEr(iray), yN1(iray, step(iray)), yN2(iray, step(iray)), &
                         yN4(iray, step(iray)), yN5(iray, step(iray)),&
                         yN6(iray, step(iray)), fN1(iray, step(iray)),&
                         fN2(iray, step(iray)), fN3(iray, step(iray)),&
                         fN4(iray, step(iray)), fN5(iray, step(iray)),&
                         fN6(iray, step(iray)), fN7(iray, step(iray)),&
                                         o_data(iray), p_data(iray))
        if ((BoundsEr(iray) > 0)) then 
             ! write(*,*) 'Boundery error right 1st Runge for Ray# ', iray
              lstop(iray) = 1
        endif
      endif
    enddo
    !OMP PARALLEL DO

    do j=1, nCoeff-1
      !$acc parallel  loop 
      !OMP PARALLEL DO shared(fN1, fN2, fN3, fN4, fN5, fN6, fN7, yN1, yN2, yN3, yN4, yN5,  yN6, yN7, yN8, p_data, o_data, step, a_data, izone, lstop)
      do iray=1, nrays
        if (lstop(iray) == 0) then
              CALL RungeInit( yN1(iray, step(iray)),yN2(iray, step(iray)), &
                              yN3(iray, step(iray)), yN4(iray, step(iray)),&
                              yN5(iray, step(iray)), yN6(iray, step(iray)),&
                              yN7(iray, step(iray)), yN8(iray, step(iray)),&

                              yN1(iray, step(iray)+1), yN2(iray, step(iray)+1), &
                              yN3(iray, step(iray)+1), yN4(iray, step(iray)+1),&
                              yN5(iray, step(iray)+1), yN6(iray, step(iray)+1),&
                              yN7(iray, step(iray)+1), yN8(iray, step(iray)+1),&

                              fN1(iray, step(iray)), fN2(iray, step(iray)), &
                              fN3(iray, step(iray)), fN4(iray, step(iray)),&
                              fN5(iray, step(iray)), fN6(iray, step(iray)),&
                              fN7(iray, step(iray)),&

                              fN1(iray, step(iray)+1), fN2(iray, step(iray)+1), &
                              fN3(iray, step(iray)+1), fN4(iray, step(iray)+1),&
                              fN5(iray, step(iray)+1), fN6(iray, step(iray)+1),&
                              fN7(iray, step(iray)+1), arrys(:,iray), &
                              lstop(iray), izone(iray), p_data(iray), a_data(iray), o_data(iray))

              step(iray) = step(iray) + 1
          endif
        enddo  
        !OMP END PARALLEL DO
      enddo
      
      
     write(*,*) 'done with 1st Runge, starting the intergration'!, yN1(1,1), fN1(1,1)
     write(*,*) 'No ok rays:', sum(lstop)
    !read (*,*) tmp

    do j= 6,  nstep

      !if (j > 2000) then
          !$acc update host(lstop)
        !  write(*,*) 'step = ', j, sum(lstop)
          if (sum(lstop) == nrays) then
            write(*,*) 'All rays are dead'
            exit
          endif
       ! endif

      !$acc update device(j)

      !$acc parallel  loop async(1)
        !OMP PARALLEL DO shared(fN1, fN2, fN4, fN5, fN6, yN1, yN2, yN4, yN5,  yN6, bPredict, lstop, next_j)
        do iray=1,nrays 
            if ((next_j(iray) <= j).and.(lstop(iray) == 0)) then
                yN1(iray, 6) = yN1(iray, 5) + (fN1(iray, 1)   * bPredict(1)  + fN1(iray, 2) * bPredict(2) + &
                                              fN1(iray, 3) * bPredict(3) +  fN1(iray, 4) * bPredict(4)  + &
                                               fN1(iray, 5) * bPredict(5)) * HstpLH
                yN2(iray, 6) = yN2(iray, 5) + (fN2(iray, 1) * bPredict(1)  + fN2(iray, 2) * bPredict(2) + &
                                        fN2(iray, 3) * bPredict(3) + fN2(iray, 4) * bPredict(4)  + &
                                        fN2(iray, 5) * bPredict(5)) * HstpLH
                yN4(iray, 6) = yN4(iray, 5) + (fN4(iray, 1) * bPredict(1)  + fN4(iray, 2) * bPredict(2) + &
                                        fN4(iray, 3) * bPredict(3) + fN4(iray, 4) * bPredict(4)  + &
                                        fN4(iray, 5) * bPredict(5)) * HstpLH
                yN5(iray, 6) = yN5(iray, 5) + (fN5(iray, 1) * bPredict(1)  + fN5(iray, 2) * bPredict(2) + &
                                        fN5(iray, 3) * bPredict(3) + fN5(iray, 7) * bPredict(4)  + &
                                        fN5(iray, 5) * bPredict(5)) * HstpLH
                yN6(iray, 6) = yN6(iray, 5) + (fN6(iray, 1) * bPredict(1)  + fN6(iray, 2) * bPredict(2) + &
                                        fN6(iray, 3) * bPredict(3) + fN6(iray, 4) * bPredict(4)  + &
                                        fN6(iray, 5) * bPredict(5)) * HstpLH
            endif
        enddo  
        !OMP END PARALLEL DO

        !$acc  parallel loop async (2)
        !OMP PARALLEL DO shared(fN1, fN2, fN3, fN4, fN5, fN6, fN7, yN1, yN2, yN3, yN4, yN5,  yN6, yN7, bPredict, lstop, next_j)
        do iray=1,nrays 
          if ((next_j(iray) <= j).and.(lstop(iray) == 0)) then        
              yN3(iray, 6) = yN3(iray, 5) + (fN3(iray, 1) * bPredict(1)  + fN3(iray, 2) * bPredict(2) + &
                                        fN3(iray, 3) * bPredict(3) + fN3(iray, 4) * bPredict(4)  + &
                                        fN3(iray, 5) * bPredict(5)) * HstpLH
              yN3(iray, 7) = yN3(iray, 5) + (fN3(iray, 1)   * bCorrect(1)  + fN3(iray, 2) * bCorrect(2) + &
                                             fN3(iray, 3) * bCorrect(3) +  fN3(iray, 4) * bCorrect(4)  + &
                                             fN3(iray, 5) * bCorrect(5)) * HstpLH
              yN7(iray, 6) = yN7(iray, 5) + (fN7(iray, 1) * bPredict(1)  + fN7(iray, 2) * bPredict(2) + &
                                        fN7(iray, 3) * bPredict(3) + fN7(iray, 4) * bPredict(4)  + &
                                         fN7(iray, 5) * bPredict(5)) * HstpLH 
              yN7(iray, 7) = yN7(iray, 5) + (fN7(iray, 1)   * bCorrect(1)  + fN7(iray, 2) * bCorrect(2) + &
                                        fN7(iray, 3) * bCorrect(3) +  fN7(iray, 4) * bCorrect(4)  + &
                                        fN7(iray, 5) * bCorrect(5)) * HstpLH                                 
              yN1(iray, 7) = yN1(iray, 5) + (fN1(iray, 1)   * bCorrect(1)  + fN1(iray, 2) * bCorrect(2) + &
                                            fN1(iray, 3) * bCorrect(3) +  fN1(iray, 4) * bCorrect(4)  + &
                                            fN1(iray, 5) * bCorrect(5)) * HstpLH 
              yN2(iray, 7) = yN2(iray, 5) + (fN2(iray, 1)   * bCorrect(1)  + fN2(iray, 2) * bCorrect(2) + &
                                        fN2(iray, 3) * bCorrect(3) +  fN2(iray, 4) * bCorrect(4)  + &
                                        fN2(iray, 5) * bCorrect(5)) * HstpLH                        
              yN4(iray, 7) = yN4(iray, 5) + (fN4(iray, 1)   * bCorrect(1)  + fN4(iray, 2) * bCorrect(2) + &
                                        fN4(iray, 3) * bCorrect(3) +  fN4(iray, 4) * bCorrect(4)  + &
                                        fN4(iray, 5) * bCorrect(5)) * HstpLH 
              yN5(iray, 7) = yN5(iray, 5) + (fN5(iray, 1)   * bCorrect(1)  + fN5(iray, 2) * bCorrect(2) + &
                                        fN5(iray, 3) * bCorrect(3) +  fN5(iray, 4) * bCorrect(4)  + &
                                        fN5(iray, 5) * bCorrect(5)) * HstpLH
              yN6(iray, 7) = yN6(iray, 5) + (fN6(iray, 1)   * bCorrect(1)  + fN6(iray, 2) * bCorrect(2) + &
                                       fN6(iray, 3) * bCorrect(3) +  fN6(iray, 4) * bCorrect(4)  + &
                                       fN6(iray, 5) * bCorrect(5)) * HstpLH    
        endif
      enddo  
      !OMP END PARALLEL DO

      !$acc wait(1) async(20)
      !$acc parallel  loop async(20)
      !OMP PARALLEL DO  shared(fN1, fN2, fN3, fN4, fN5, fN6, fN7, yN1, yN2, yN4, yN5,  yN6, o_data, p_data, next_j, lstop, BoundsEr)
      do iray=1, nrays
        if ((next_j(iray) <= j).and.(lstop(iray) == 0)) then 
            CALL lscsq_ftion(BoundsEr(iray), yN1(iray, 6), yN2(iray, 6), yN4(iray, 6), yN5(iray, 6), yN6(iray, 6), &
                                              fN1(iray, 6), fN2(iray, 6), fN3(iray, 6),fN4(iray, 6), &
                                              fN5(iray, 6), fN6(iray, 6), fN7(iray, 6), &
                                              o_data(iray), p_data(iray))                     
          endif
      enddo
      !OMP END PARALLEL DO
      !$acc wait


      !$acc parallel  loop async(1)
      !OMP PARALLEL DO  shared(fN1, fN2, fN4, fN5, fN6, yN1, yN2, yN4, yN5,  yN6, bCorrect, next_j, lstop, BoundsEr)
      do iray=1, nrays
        if ((next_j(iray) <= j).and.(lstop(iray) == 0).and.(BoundsEr(iray) == 0)) then      
                yN1(iray, 6) = yN1(iray, 7) + bCorrect(6)* fN1(iray, 6) * HstpLH      
                yN2(iray, 6) = yN2(iray, 7) + bCorrect(6)* fN2(iray, 6) * HstpLH      
                yN4(iray, 6) = yN4(iray, 7) + bCorrect(6)* fN4(iray, 6) * HstpLH   
                yN5(iray, 6) = yN5(iray, 7) + bCorrect(6)* fN5(iray, 6) * HstpLH 
                yN6(iray, 6) = yN6(iray, 7) + bCorrect(6)* fN6(iray, 6) * HstpLH   
        endif
      enddo
      !OMP END PARALLEL DO


      !$acc  parallel loop async(2)
      !OMP PARALLEL DO  shared(fN3, fN7, yN3, yN7, bCorrect, next_j, lstop, BoundsEr)
      do iray=1,nrays 
        if ((next_j(iray) <= j).and.(lstop(iray) == 0).and.(BoundsEr(iray) == 0)) then     
                yN3(iray, 6) = yN3(iray, 7) + bCorrect(6)* fN3(iray, 6) * HstpLH      
                yN7(iray, 6) = yN7(iray, 7) + bCorrect(6)* fN7(iray, 6) * HstpLH  
        endif
      enddo

      !$acc wait(1) async(20)
      !$acc parallel loop async(20)
      !OMP PARALLEL DO  shared(fN1, fN2, fN3, fN4, fN5, fN6, fN7, yN1, yN2, yN4, yN5,  yN6, o_data, p_data, next_j, lstop, BoundsEr)
      do iray=1, nrays
        if ((next_j(iray) <= j).and.(lstop(iray) == 0).and.(BoundsEr(iray) == 0)) then 
          CALL lscsq_ftion(BoundsEr(iray), yN1(iray, 6), yN2(iray, 6), yN4(iray, 6), yN5(iray, 6), yN6(iray, 6), &
                                           fN1(iray, 6), fN2(iray, 6), fN3(iray, 6),fN4(iray, 6), &
                                           fN5(iray, 6), fN6(iray, 6), fN7(iray, 6), &
                                           o_data(iray), p_data(iray))!)!, all_vecs)!, all_const, all_vecs, local_psigrd)
        endif
      enddo
       !OMP END PARALLEL DO


      !$acc wait

      !$acc parallel  loop 
      !OMP PARALLEL DO shared(yN3, yN7, yN8, yN1, yN2, yN4, yN5,  yN6, o_data, p_data, next_j, lstop, BoundsEr,izone, arrys, a_data)
      do iray=1, nrays
        if ((next_j(iray) <= j).and.(lstop(iray) == 0).and.(BoundsEr(iray) == 0)) then 
          yN8(iray, 6) = yN8(iray, 5) + HstpLH
          call lscsq_E2byPr(yN1(iray, 6), yN2(iray, 6), yN3(iray, 6), yN4(iray, 6), yN5(iray, 6), &
                            yN6(iray, 6), yN7(iray, 6), yN8(iray, 6), &
                            izone(iray), lstop(iray), arrys(:, iray), &
                            p_data(iray), a_data(iray), o_data(iray))              
        endif
      enddo
      !OMP END PARALLEL DO

      
      
      do k=1, nCoeff
        !$acc parallel  loop 
        !OMP PARALLEL DO shared(fN1, fN2, fN3, fN4, fN5, fN6, fN7, yN1, yN2, yN3, yN4, yN5,  yN6, yN7, yN8, next_j, lstop, BoundsEr)
        do iray=1, nrays
          if ((next_j(iray) <= j).and.(lstop(iray) == 0) .and. (BoundsEr(iray) == 0)) then
            
            yN1(iray, k) = yN1(iray, k+1)
            yN2(iray, k) = yN2(iray, k+1)
            yN3(iray, k) = yN3(iray, k+1)
            yN4(iray, k) = yN4(iray, k+1)
            yN5(iray, k) = yN5(iray, k+1)
            yN6(iray, k) = yN6(iray, k+1)
            yN7(iray, k) = yN7(iray, k+1)
            yN8(iray, k) = yN8(iray, k+1)

            fN1(iray, k) = fN1(iray, k+1)
            fN2(iray, k) = fN2(iray, k+1)
            fN3(iray, k) = fN3(iray, k+1)
            fN4(iray, k) = fN4(iray, k+1)
            fN5(iray, k) = fN5(iray, k+1)
            fN6(iray, k) = fN6(iray, k+1)
            fN7(iray, k) = fN7(iray, k+1)
          endif
        enddo 
        !OMP END PARALLEL DO
      enddo

      !$acc parallel  loop 
      !OMP PARALLEL DO shared(yN1, yN2, yN4, yN5,  o_data, p_data, next_j, lstop, BoundsEr)
      do iray=1, nrays
        if ((next_j(iray) <= j).and.(lstop(iray) == 0) .and. (BoundsEr(iray) > 0)) then
            call lscsq_BounceIt(yN1(iray, 5), yN2(iray, 5), yN4(iray, 5), yN5(iray, 5), o_data(iray), p_data(iray))
        ENDIF
      enddo
       !OMP END PARALLEL DO


      !$acc parallel  loop 
      !OMP PARALLEL DO shared(yN3, yN7, yN8, yN1, yN2, yN4, yN5,  yN6, next_j, lstop, BoundsEr)
      do iray=1, nrays
        if ((next_j(iray) <= j).and.(lstop(iray) == 0) .and. (BoundsEr(iray) > 0)) then
            step(iray) = 1
            yN1(iray, step(iray)) = yN1(iray, 5)
            yN2(iray, step(iray)) = yN2(iray, 5)
            yN3(iray, step(iray)) = yN3(iray, 5)
            yN4(iray, step(iray)) = yN4(iray, 5)
            yN5(iray, step(iray)) = yN5(iray, 5)
            yN6(iray, step(iray)) = yN6(iray, 5)
            yN7(iray, step(iray)) = yN7(iray, 5)
            yN8(iray, step(iray)) = yN8(iray, 5)
        endif
      ENDDO
       !OMP END PARALLEL DO

      !$acc parallel loop 
      !OMP PARALLEL DO shared(fN1, fN2, fN3, fN4, fN5, fN6, fN7, yN1, yN2, yN4, yN5,  yN6, o_data, p_data, next_j, lstop, BoundsEr)
      do iray=1, nrays
        if ((next_j(iray) <= j).and.(lstop(iray) == 0).and.(BoundsEr(iray) > 0)) then 
            BoundsEr(iray) = 0
            CALL lscsq_ftion(BoundsEr(iray), yN1(iray, step(iray)), yN2(iray, step(iray)),&
                             yN4(iray, step(iray)), yN5(iray, step(iray)), yN6(iray, step(iray)), &
                             fN1(iray, step(iray)), fN2(iray, step(iray)), fN3(iray, step(iray)),&
                             fN4(iray, step(iray)), fN5(iray, step(iray)), fN6(iray, step(iray)),&
                             fN7(iray, step(iray)), o_data(iray), p_data(iray))
        !   if (iray == 24) then
           ! write(*,*) 'Got correction F:', yN1(iray, 6), fN1(iray, 6), BoundsEr(iray)
        !  endif 
            if ((BoundsEr(iray) > 0)) then 
                !  write(*,*) 'Boundery error right after Bounce!'
                  lstop(iray) = 1
            else
              BoundsEr(iray)  = 100
            endif
          endif
      enddo
       !OMP END PARALLEL DO
 
      do k=1, nCoeff-1 
          !$acc parallel  loop 
        !OMP PARALLEL DO shared(fN1, fN2, fN3, fN4, fN5, fN6, fN7, yN1, yN2, yN3, yN4, yN5,  yN6, yN7, yN8, p_data, o_data, step, a_data, izone, lstop, next_j, BoundsEr)
          do iray=1, nrays
            if ((next_j(iray) <= j).and.(lstop(iray) == 0) .and. (BoundsEr(iray)  == 100)) then

              if (k == 4) then
                next_j(iray) = j + 5
                BoundsEr(iray)= 0
              endif

              CALL RungeInit( yN1(iray, step(iray)), yN2(iray, step(iray)), &
                              yN3(iray, step(iray)), yN4(iray, step(iray)),&
                              yN5(iray, step(iray)), yN6(iray, step(iray)),&
                              yN7(iray, step(iray)), yN8(iray, step(iray)),&
  
                              yN1(iray, step(iray)+1), yN2(iray, step(iray)+1), &
                              yN3(iray, step(iray)+1), yN4(iray, step(iray)+1),&
                              yN5(iray, step(iray)+1), yN6(iray, step(iray)+1),&
                              yN7(iray, step(iray)+1), yN8(iray, step(iray)+1),&

                              fN1(iray, step(iray)), fN2(iray, step(iray)), &
                              fN3(iray, step(iray)), fN4(iray, step(iray)),&
                              fN5(iray, step(iray)), fN6(iray, step(iray)),&
                              fN7(iray, step(iray)),&
  
                              fN1(iray, step(iray)+1), fN2(iray, step(iray)+1), &
                              fN3(iray, step(iray)+1), fN4(iray, step(iray)+1),&
                              fN5(iray, step(iray)+1), fN6(iray, step(iray)+1),&
                              fN7(iray, step(iray)+1), arrys(:,iray), &
                              lstop(iray), izone(iray), p_data(iray), a_data(iray), o_data(iray))!, local_prmtrs)!, all_vecs)
                              ! ,all_const, &
                              ! all_vecs, local_psigrd, local_dvol, local_vpar)

              !if (lstop(iray) == 1) write(*,*) 'Out of Boundery after Runge for Ray# ', iray
              step(iray) = step(iray) + 1

            endif
          enddo
           !OMP END PARALLEL DO
      enddo  
      
      
    enddo 

    !$acc end data

    ! write(*,*) 'EXIT'
    ! do j=1, nrays 
    !   if (lstop(j)==0) write(*,*) j
    ! enddo

    ! do iray=1, nrays
    !   do iz=1, nzones
    !     if ((arrys(iz, iray)%ivind > nv)) write(*,*) '>nv', iray, iz, arrys(iz, iray)%ivind,arrys(iz, iray)%izind
    !     if ((arrys(iz, iray)%izind > npsi)) write(*,*) '<npsi', iray, iz, arrys(iz, iray)%ivind,arrys(iz, iray)%izind
    !     if ( arrys(iz,iray)%dlnPdsK .ne. arrys(iz,iray)%dlnPdsK )  write(*,*)  'Nan1',iz,iray
    !     if ( arrys(iz,iray)%dlnPdsx .ne. arrys(iz,iray)%dlnPdsx )  write(*,*)  'Nan2',iz,iray
    !   !  if (iz < 5)  write(*,*) arrys(iz, iray)%ivind,arrys(iz, iray)%izind,  arrys(iz,iray)%dlnPdsK , arrys(iz,iray)%dlnPdsK 
    !   ENDDO
    ! enddo
    end subroutine lscsq_PredcLSC



    !
    !     -----------------------------------------------------------------
    !     -----------------------------------------------------------------
    !

    ! subroutine Bounce(y1, y2, y4, y5, o_data, p_data)
    !   !, all_const, all_vecs, local_psigrd)
    !   !subroutine Bounce(y, yok, iBndsErr, o_data, p_data)
    !   use pl_types 
    !   use iso_c_binding, only : fp => c_double
    !   ! use lscsq_mod, only : neqsp1
    !   ! use lscsq_mod, only : incithet, scatthet!, iscatplt 
    !   ! use lscsq_mod, only : deg2rad
    !   ! use lscsq_mod, only : scatkdeg

    !   implicit NONE
    !   integer :: nBndsErr = 25
    !   !integer, intent(INOUT) :: iBndsErr, iscatplt
    !   Real(fp), intent(INOUT):: y1, y2, y4, y5
    !   type(old_data), intent(inout) :: o_data
    !   type(pl_data), intent(inout) :: p_data

    !   !Real(fp), intent(INOUT) :: yok(NEQSP1) 
    !   !Real(fp), intent(INOUT) :: thi, tho

    !   ! type(const), intent(inout) :: all_const
    !   ! type(vecs), dimension(npsij), intent(inout) :: all_vecs
    !   ! !Real(fp), dimension(5,5), intent(inout) :: small_psigrd
    !   ! Real(fp), dimension(nx,nz), intent(in) :: local_psigrd

    !   !$acc routine seq
    !   !$acc routine(lscsq_BounShft)
    !   !$acc routine(lscsq_BounceIt)
     
    !   !if (ScatKdeg.gt.1.0e-3_fp) then
    !    ! CALL lscsq_BounShft(y, NEQSP1, thi, tho, o_data, p_data)
    !    ! inciThet(iscatplt) = thi * deg2rad
    !    ! scatThet(iscatplt) = tho * deg2rad
    !    ! iscatplt=iscatplt+1
    !   !else

    !     CALL lscsq_BounceIt(y1, y2, y4, y5, o_data, p_data )!, all_const, all_vecs, local_psigrd)
    !  ! endif

    !   ! iBndsErr = iBndsErr + 1
    !   ! if (iBndsErr .GE. nBndsErr) then
    !   !   write(*, *) ' Rays out of bounds; Times of recovery: ', iBndsErr
    !   !   iBndsErr = 0
    !   ! endif
   
    ! END subroutine Bounce
    !
    !     -----------------------------------------------------------------
    !     -----------------------------------------------------------------
    !
    SUBROUTINE lscsq_BounceIt (y1, y2, y4, y5, o_data, p_data)!, local_prmtrs)!, all_vecs)!, all_const, all_vecs, local_psigrd)
      use pl_types 
      use iso_c_binding, only : fp => c_double
      implicit none
     
      type(old_data), intent(inout) :: o_data
      type(pl_data), intent(inout) :: p_data
      !type(parameters), intent(inout) :: local_prmtrs

      Real(fp), intent(in):: y1, y2
      Real(fp), intent(inout):: y4, y5
      Real(fp) ::    r,z,   psi, Br, Bz, RBphi, omc, Tee, pe2, pi2, aio, ael
      Real(fp) ::    B2, crr, crz, czr, czz, KrNew, KzNew
      !     Go back to the last ok solution point, and re-arrange kr and kz
      !     such that k \dot B is unchanged but k \cross B changes sign.  Then
      !     put this new information into the starting condition for y.
      !     y(1)  y(2)  y(3)  y(4)  y(5)  y(6)
      !      r     z    phi   k_r   k_z    n

      ! type(const), intent(inout) :: all_const
      ! type(vecs), dimension(npsij), intent(inout) :: all_vecs
      ! !Real(fp), dimension(5,5), intent(inout) :: small_psigrd
      ! Real(fp), dimension(nx,nz), intent(in) :: local_psigrd

      !$acc routine seq
      !$acc routine(lscsq_plasma2d)



      !r  = yok(1)
      !z  = yok(2)
      r  = y1
      z  = y2
!all_const, all_vecs, local_psigrd, 
      call  lscsq_plasma2d   (o_data, p_data, r,z,psi,Br,&
      Bz,RBphi,omc,&
     Tee,pe2,pi2,aio,ael)!, all_vecs)

      

      B2  =  Br*Br + Bz*Bz
      crr = (Br*Br - Bz*Bz)/B2
      czz = - crr
      crz = 2.0_fp*Br*Bz/B2
      czr = crz

      KrNew = crr * y4 + crz * y5
      KzNew = czr * y4 + czz * y5

      y4 = KrNew
      y5 = KzNew

      !do i=1, n
       ! y(i) = yok(i)
      !enddo

    end subroutine lscsq_bounceit
    !
    !     -----------------------------------------------------------------
    !     -----------------------------------------------------------------
    !
    ! subroutine lscsq_BounShft (y, n, thi, tho, o_data, p_data)
    !   use pl_types
    !   use iso_c_binding, only : fp => c_double
    !   use lscsq_mod, only : deg2rad
    !   use lscsq_mod, only : scatkdeg
    !   implicit none
    !   type(old_data), intent(inout) :: o_data
    !   type(pl_data), intent(inout) :: p_data

    !   integer, intent(in) :: n
    !   !Real(fp), dimension(n), intent(inout) :: yok
    !   Real(fp), dimension(n), intent(inout) :: y

    !   integer :: i
    !   Real(fp):: tauFWP
    !   Real(fp):: r,z, psi, Br, Bz, RBphi, omc, Tee, pe2, pi2, aio, ael
    !   Real(fp):: Bph, kro, kzo, kpho, thi, tho

    !   !$acc routine seq
    !   !$acc routine(lscsq_GetKscat)
    !   !$acc routine(lscsq_plasma2d)
      
    !   r   = y(1)
    !   z   = y(2)
    !   call lscsq_plasma2d (o_data, p_data, r,z,psi,Br,&
    !   Bz,RBphi,omc,&
    !     Tee,pe2,pi2,aio,ael)

    !   Bph = RBphi/r
    !   tauFWP =  (ScatKdeg*deg2rad)**2
    !   CALL lscsq_GetKscat(y(4), y(5), y(6)/r, Br, Bz, Bph, tauFWP, kro, kzo, kpho, thi, tho)
    !   y(4) = kro
    !   y(5) = Kzo
    !   y(6) = kpho*r

    !   !do i=1, n
    !   !y(i) = yok(i)
    !   !enddo
    
    ! end subroutine lscsq_bounshft
      !
      !     -----------------------------------------------------------------
      !     -----------------------------------------------------------------
      !
    subroutine lscsq_GetKscat( kri, kzi, kphi, Br, Bz, Bph, tau, kro, kzo, kpho, thi, tho )
      !     Get K scattered by fluctuations.
      !     Copyright 1993 by F. W. Perkins, H. Takahashi, D. W. Ignat, & E. J. Valeo
      !     Given the initial k in r,z,phi space and the B at that point
      !     return the out k
      !     by transforming to a space in grad psi, transverse, and parallel
      !     referred to as                  rr    ,   tt      ,     ll
      !     rotating the k randomly around the ll (field) direction thus preserving
      !     k_parallel (kll)  and k_perp (kperp) but not k-phi
      !     The parameter tau controls the randomness.
      !     If tau is large then the scattered k peaks perpendicular to the wall
      !     If tau is small, then the reflection is specular, without randomization
      !
      !     The angle is measured with respect to the tt direction

      !      tt ( perp to B and to grad psi )
      !      ^
      !      |  theta(i)
      !      |          .
      !      |     .                           (parallel to B directionout of paper)
      !      | .
      !      +------- >  rr (grad poloidal flux
      !
      !     Following Perkins the Prob(ability) of scattering from thetai to theta is
      !
      !          Prob ~ sin(theta) exp{ - (theta-thetai)^2/tau }
      !     and one forms the integral
      !          ProbIntl(theta,thetai,tau) which is zero at theta=0 and 1 at theta=PI
      !     so that if a random number is chosen between 0 and 1, a unique angle
      !     is determined.  If tau is big, PI/2 is most likely; if tau is small
      !     then thetai is most likely
      !
      !
      use iso_c_binding, only : fp => c_double
      use lscsq_mod, only : rad2deg
      implicit none

      INTEGER i,j,  isimp, ifirst, NumThets
      INTEGER ithetai, itheta
      INTEGER idum
      Real(fp) ::    kri, kzi, kphi
      Real(fp) ::     kro, kzo, kpho, thi, tho
      Real(fp) ::     krr, ktt, kll, kperp, signth
      Real(fp) ::     ktoti, ktoto, ktott
      Real(fp) ::     Br, Bz, Bph, Bp, B, tau
      Real(fp) ::     TAUMAX, TAUMIN, tauold, thetai, dtheta
      Real(fp) ::     cosalph, sinalph, cosbeta, sinbeta, simp(2)
      Real(fp) ::     PI, PIby2, SMALL, ERR
      Real(fp) ::     random
!!!      Real(fp) ::     lscsq_ran3
      PARAMETER (NumThets=101)
      !      PARAMETER (TAUMAX=100.0_fp, TAUMIN=1.0e-05)
      DATA       TAUMAX ,       TAUMIN /                                &
      !&     1.0d2, 1.0d-5/
      &     1.0_fp, 1.0_fp/
      !     if  exp{ - (th_out_deg-th_in_deg)^2/del_th_deg^2 }
      !     and exp{ - (th_out_rad-th_in_rad)^2/tauFWPerkins }
      !     then    (del_th_deg,tauFWP) = (0.1deg,2.e-6) (1.0deg,3e-4) (360deg,40.0_fp)
      Real(fp) ::     ProbIntl(NumThets, NumThets),theta(NumThets)
      !Real(fp) ::      tauDEG
      DATA ifirst,       SMALL   ,       ERR     /                      &
      !&        1  , 1.0d-25, 1.0d-5 /
      & 1, 1.0_fp, 1.0_fp /
      DATA       simp(1),       simp(2) /                               &
      ! &     4.0d0 , 2.0d0  /
      & 4.0_fp, 2.0_fp/
      Real(fp) :: re41

      !---------------------------------------------
      SAVE :: PI,PIby2,dtheta,theta,probintl,tauold
      !---------------------------------------------
      !     If first call, set up PI and theta array
      !     Also if first CALL lscsq_set the ProbIntl to the result for large tau

      !$acc routine seq
      !$acc routine(lscsq_ran3)

      if (ifirst .EQ. 1) then
        ifirst = 0
        PIby2  = asin(1.0_fp)
        PI     = PIby2 + PIby2
        dtheta = PI / real(NumThets - 1,kind=fp)
        do i = 1,NumThets
          theta(i) = dtheta * real(i-1,kind=fp)
          ProbIntl(i,1) = 0.5_fp * ( 1.0_fp - cos (theta(i)) )


          ProbIntl(i,2:numthets) = ProbIntl(i,1)
        enddo
        tauold = TAUMAX
      endif
      !
      !     Fill the array of ProbIntl unless the values from the last call
      !     or from the first CALL lscsq_set up are ok
      !
      if ( tau .ne. tauold .and. (tau .LE. tauold .or. tauold .ne. TAUMAX) ) then

        ProbIntl(1,1:numthets) = 0.0_fp

        do j=1,NumThets
          thetai = theta(j)
          do i=2,NumThets
            isimp = mod(i,2) + 1
            ProbIntl(i,j) = ProbIntl(i-1,j)+simp(isimp)*sin(theta(i))*exp(-(theta(i)-thetai)**2/tau)
          enddo
        enddo

        !       do 30 j=1,NumThets
        !       do 30 i=2,NumThets
        do j=1,NumThets
          do i=2,NumThets
            ProbIntl(i,j) =  ProbIntl(i,j)/ProbIntl(NumThets,j)
          enddo
        enddo
        ! 30     continue

        !cccc  Contour plot begins
        !     ------------------------------------------------------------------
        !     kclev1: >0 --> clevelin contains levels to be used;
        !                    dots for index less than kclev2;
        !                    solid for index greater or equal kclev2
        !     kclev1: =0 --> first contour at clevelin(1) with
        !                    next one up by clevelin(2), and so on and so on
        !     kclev1: <0 --> rcontr to choose -kclev1 equally spaced values between
        !                    clevelin(1) and clevelin(2)
        !     clevelin:      array of contour levels; this is output if kclev1<0
        !     kclev2:        separates dots from solid lines
        !     CALL lscsq_EZrcon(ix1(1),ix2(1), jy1(1), jy2(1),
        !    ^             xa,    xb,     ya,     yb   ,
        !    ^  kclev1,clevelin,kclev2,
        !    ^  PsiContr,
        !    ^  i1stdim,
        !    ^  xAry, ixmin, ixmax, ixstep,
        !    ^  yAry, jymin, jymax, jystep )
        !      do i=1,9
        !         clevelin(i)=0.1*REAL(i,kind=fp)
        !      enddo
        !      CALL lscsq_EZrcon(150, 500, 250, 600, ZERO, 3.2_fp, ZERO, 3.2_fp,   &
        !     &  9,clevelin,5,                                                   &
        !     &  ProbIntl,                                                       &
        !     &  NumThets,                                                       &
        !     &  theta, 1, NumThets-1, 1,                                        &
        !     &  theta, 1, NumThets-1, 1)
        !      CALL lscsq_EZwrit( 150, 150,                                        &
        !     &            'ProbIntl; theta_in ordinate$',0,0)
        !      CALL lscsq_EZwrit( 150, 125,                                        &
        !     &            'theta_out abcsissa(radians)$',0,0)
        !      tauDEG = sqrt(tau)*rad2deg     
        !      write(MyString,'(''tauDEG, tauFWP: '',                            &
        !     &                   f5.1,1x,1pe9.2,''$'')')tauDEG,tau
        !      CALL lscsq_EZwrit( 150, 100,MyString,0,0)
        !      CALL lscsq_EZfini(0,0)
        !      CALL lscsq_MkGrfLst   (' ProbIntl plot ')
        !cccc  Contour plot ends
      endif
      !     .
      !     .                                 Here is the normal starting point !!
      !     .
      Bp  = sqrt(Br*Br + Bz*Bz) + SMALL
      B   = sqrt(Br*Br + Bz*Bz + Bph*Bph) + SMALL
      cosalph = abs(Bph)/B
      sinalph =     Bp  /B
      cosbeta =     Br  /Bp
      sinbeta =     Bz  /Bp

      krr =         sinbeta*kri -         cosbeta*kzi !+    0.0_fp*kphi
      ktt = cosalph*cosbeta*kri + cosalph*sinbeta*kzi - sinalph*kphi
      kll = sinalph*cosbeta*kri + sinalph*sinbeta*kzi + cosalph*kphi

      kperp = sqrt(krr*krr + ktt*ktt)
      ktott = sqrt(krr*krr + ktt*ktt + kll*kll)
      if ( krr  .LT. 0.0_fp ) then
        signth = - 1.0_fp
        krr    = abs(krr)
        thetai = asin (krr/kperp)
        if (ktt .LT. 0.0_fp ) thetai = PI - thetai
      else
        signth = + 1.0_fp
        thetai = asin (krr/kperp)
        if (ktt .LT. 0.0_fp ) thetai = PI - thetai
      endif
      RE41 = thetai/dtheta
      ithetai = int(RE41) + 1
      !     ithetai = ifix(thetai/dtheta) + 1
      thi = theta(ithetai)
      !
      if (tau .LT. TAUMIN) then
        krr = -signth*krr
        ktt = +       ktt
        tho = thi
      else
      !
      idum=1
      !call random_seed()
      !call random_number(random)
      random = lscsq_ran3(idum)
      !
      !       do 50 i=2,NumThets
        do i=2,NumThets
          itheta = i
          !         if (ProbIntl(i,ithetai) .GT. random ) go to 51
          if (ProbIntl(i,ithetai) .GT. random ) exit     
        enddo

      ! 50     continue
      ! 51     continue


        krr = -signth*kperp* sin(theta(itheta))
        ktt =         kperp* cos(theta(itheta))
        tho = theta(itheta)

      endif

      kro =  sinbeta*krr + cosalph*cosbeta*ktt + sinalph*cosbeta*kll
      kzo = -cosbeta*krr + cosalph*sinbeta*ktt + sinalph*sinbeta*kll
      kpho=     0.0_fp*krr - sinalph*        ktt + cosalph*        kll

      ktoti = sqrt(kri*kri + kzi*kzi + kphi*kphi)
      ktoto = sqrt(kro*kro + kzo*kzo + kpho*kpho)
      !
      !     Take this out June 2000, below
      !      if ( abs(ktoti-ktoto) .gt. ERR*ktoti .or.
      !     ^     abs(ktoti-ktott) .gt. ERR*ktoti .or.
      !     ^     abs(ktott-ktoto) .gt. ERR*ktoto  ) then
      !cc       CALL lscsq_LSCpause
      !        write(6,'(''GetKscat error; ktoti, ktoto, ktott:'',
      !     ^                           3(1x,1pe10.3))') ktoti, ktoto, ktott
      !cc       CALL lscsq_LSCpause
      !      endif
      !     Take this out June 2000, above
      !
      tauold = tau
      return
    END subroutine lscsq_GetKscat
    !
    !-----------------------------------------------------------------
    !-----------------------------------------------------------------

    !Real(fp) 
    FUNCTION lscsq_ran3(idum)
      !     Transcribed 1993 by D. W. Ignat
      !     References:
      !     W. H. Press et al, Numerical Recipes in Fortran, p199
      !     D. Knuth, Seminumerical Algorithms, Vol 2 of The Art of Computer
      !          Programming
      !     Returns a uniform random deviate between 0.0_fp and 1.0.
      !     Set idum to any negative value to initialize or reinitialize the
      !     sequence.
      !     Substitute the CCCommented lines for the ones following to
      !     render the routine entirely floating point.
      use iso_c_binding, only : fp => c_double
      implicit none
      INTEGER idum
      integer :: iff, mj,mk
      INTEGER :: i,ii,k, inext, inextp

      real(fp) :: lscsq_ran3
      ! 55 dimension is special; no changes!
      integer :: ma(55)
      integer, parameter :: mbig=1.0e9
      integer, parameter :: mseed=161803398
      integer, parameter :: mz=0
      Real(fp), parameter :: fac=1.0_fp/mbig
      DATA iff /0/
      !  Initialize

      !$acc routine seq

      if(idum .LT. 0 .or. iff .EQ. 0) then
      iff=1
      ! Initialize ma(55) using the seed idum
      ! and the large number MSEED
      mj=MSEED-abs(idum)
      mj=mod(mj,MBIG)
      ma(55)=mj
      mk=1
      ! Now initialize the rest of the table
      ! in a slightly random order
      ! with numbers that are not especially
      ! random
      do i=1,54
      ii=mod(21*i,55)
      ma(ii)=mk
      mk=mj-mk
      if(mk .LT. MZ) mk=mk+MBIG
      mj=ma(ii)
      enddo

      ! We randomize them by
      ! 'warming up the generator'
      do k=1,4
        do i=1,55
          ma(i)=ma(i)-ma(1+mod(i+30,55))
          if(ma(i) .LT. MZ) ma(i)=ma(i)+MBIG
        enddo
      enddo
      ! Prepare indices for our first
      ! generated number.   The constant 31
      ! is special...see Knuth
      inext=0
      inextp=31
      idum=1
      endif
      !
      !     .                                 Here is where we start usually
      !     .                                 Increment inext, wrap around 56 to 1
      inext=inext+1
      if(inext .EQ. 56)inext=1
      inextp=inext+1
      !     .                                 Ditto for inextp
      if(inextp .EQ. 56)inextp=1
      !     .                                 Now generate a new random number
      mj=ma(inext)-ma(inextp)
      !     .                                 Be sure it is in the range
      if(mj .LT. MZ)mj=mj+MBIG
      !     .                                 Store it, and output ran3
      ma(inext)=mj
      lscsq_ran3=mj*FAC
      return
    END FUNCTION lscsq_ran3
    !
end module PredcLSC

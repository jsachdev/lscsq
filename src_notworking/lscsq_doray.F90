subroutine lscsq_DoRay(arrys)
    use pl_types 
    !USE OMP_LIB
    
    use PredcLSC
    use Tracing
    use rayIni

    use iso_c_binding, only : fp => c_double

    use lscsq_mod, only : iError, rmin, psimin, psilim
    use lscsq_mod, only:  cEparIK, delpsi, &
                          vpar, omcfac, pe2fac, nzones, nrays, &
                          zmin, psigrd,  &
                          hstplh, nstep, &
                          lfast, begin, ntors,&
                          pe2vec, pi2vec, teKev, aiovec, aelvec, RBpvec, psivec, &
                          dvol, ivzero!, NpsiJ, twopi, vc, deg2rad,ntor, npol, thgrid, ind_ray, lfast,&
                        ! begin, rmax, zmax, rmaj, SEARCHINCR, fghz, ok_ray, npar
    
    implicit none
    type(storeAry), dimension(nzones, nrays), intent(inout) :: arrys
    !integer :: iray, iErrCount


   !type(const),  intent(inout) :: all_const
  !  type(pl_data), optional, intent(inout) :: p_data
  !  type(old_data), dimension(nrays), intent(inout) :: o_data
   !type(vecs), dimension(npsij), intent(inout) :: all_vecs

    !omcfac, psimin, psilim,npsij, nx, nz, nv, 
  !global arrays
  !$acc update device(pe2vec, pi2vec, teKev, aiovec, aelvec, RBpvec, psivec)
  !$acc update device(dvol, vpar, psigrd)

  !global scalars
  !$acc update device(psimin, psilim, hstplh, nzones, nrays)
  !$acc update device(OmcFac, zmin, rmin, ivzero, pe2Fac , cEparIK, delpsi)

  !! $acc update device(twopi, vc, deg2rad,ntor, npol, thgrid, ind_ray, lfast)
  !! $acc update device(begin, rmax, zmax, rmaj, SEARCHINCR, fghz,ok_ray, npar)

    
  ! !$acc update device(nzones)
  ! !$acc update device(lfast, begin, ntors, hstplh)
  ! !$acc update device(zmin, rmin)
  ! !$acc update device(ivzero)
 

  !variables for BounceShft - not using right now
  !acc update device(incithet, scatthet, scatkdeg)
  
  !variable used by lscsq_LSCendr - in the master this subroutine was called if Runge failed. 
  !If !=0 then main routine will through a message, 
  !we can remove this global variable and through message in subroutine lscsq_PredcLSC(arrys)
  !acc update device(iEndRy)

  !No need to keep global variable for formatting printout in
  !acc update device( nTSCscrn)

  !Remove this variable for now. During integaration it used in lscsq_plasma1d
  !for message ' negative temperature correction ' - never had this message, so we can revisit it later
  !acc update device(iRayTrsi)
  
    

  call lscsq_RayDataInit(arrys)

  call lscsq_PredcLSC(arrys)

  write(*,*) 'Error after lscsq_PredcLSC - ', iError

 
  end subroutine lscsq_doray
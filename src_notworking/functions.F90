module functions
    contains


              

  
  ! ------------------------------------------------------
  !     OMPN  function for integration       
  
  
  
  
  subroutine lscsq_ftion(BoundsEr, y, f, o_data, p_data, vecs)
    use grap
    use grapgrd
    use plasma1d
 
    use pl_types 
    use iso_c_binding, only : fp => c_double
    use lscsq_mod, only : neqs, npsij
    implicit none
    type(pl_data), intent(inout) :: p_data
    type(old_data), intent(inout) :: o_data
    type(vec_data), dimension(npsij), intent(inout)  :: vecs
    Real(fp) :: D11er, D11w0, D33w0, D12w0
    !, D33er, D12er, D11ar, D33ar, D12ar,&
                
    Real(fp) :: denom, dkpar, dkper 
    Real(fp), intent(inout), dimension(neqs) :: f, y 
  
  !!!!!!     EXTERNAL plasma2d, eps, epsdr
  !     KdB   k \cdot B = ( k_r B_r + k_z B_z + (n/r) B_{\phi} )
  !     KdK   k \cdot k = ( k_r^2   + k_z^2   + (n/r)^2 )
  !     \p = \partial , the Greek for partial derivative
  !     \R = {\bf r}  , the vector r
  !     \K = {\bf k}  , the vector k
  !     \abs{#1} = \mid #1 \mid , the absolute value
  !     s  =  s       , the abs value of \R; ds = path length
  !     D  = D(\R,\K,\omega) ; D = 0 is the dispersion relation
  !     D  = Kper4 * Eper
  !        + Kper2 * [ Qpar*(Epar+Eper) + woc2*Exy**2 ]  ! Kper2 * bb
  !        +  Epar * [(Qpar )**2 - (woc2*Exy)**2 ]       ! Epar  * cc
  !     bkb(3) = 0.5 \p k_{\parallel}^2 / \p (k_r,k_z,n)
  !     per(3) = 0.5 \p k_{\perp}^2     / \p (k_r,k_z,n)
  !     dRdwt  = -\frac{\p D/\p\K}{\omega \p D/\p\omega}
  !     dKdwt  = +\frac{\p D/\p\R}{\omega \p D/\p\omega}
  !     dsdwt  = \abs{dRdwt}
  !     dRds   = dRdwt/dsdwt      == f(1,2,3)
  !     dKds   = dKdwt/dsdwt      == f(4,5,6)
  !     dwtds  = 1./dsdwt
  !     wdDdw  = \omega \p D/\p \omega
  !     DKpar ==  \p D/\p k_{\parallel}^2
  !     DKper ==  \p D/\p k_{\perp}^2
  !     Qpar = Kpar2 - woc2*Eper
  !     QparE= \p Qpar/\p Kper2
  !     QparA= \p Qpar/\p Kpar2
  !     QparW= \p Qpar/\p \omega^2  \cdot \omega^2
  !
    integer, intent(inout) :: BoundsEr
    Real(fp) :: aion, ipsq,  max_d!, lscsq_DispRela
    Real(fp) :: KdB,KdK,  Kper2,Kper4,Kpar2
    Real(fp) :: Btot2, Qpar, QparE,  QparW, bb, cc!, QparA
    Real(fp) :: bkb(3), per(3), dRdwt(3), dKdwt(3),  dsdwt
    Real(fp) :: r, z, RBphi, omc, Tee, pi2, aio, ael, &
    psi, pe2
    Real(fp) :: Br, Bz
    Real(fp) :: dr = 1.0e-4_fp      
    Real(fp) :: dz = 1.0e-4_fp  
    
      !acc routine vector
      !acc routine(lscsq_eps)
      !acc routine(lscsq_plasma2d)  
      !acc routine(lscsq_DispRela)
  
  
    r = Y(1)
    z = Y(2)
  
    call lscsq_plasma2d (o_data, p_data, vecs, r,z,psi,Br,&
                          Bz,RBphi,omc,&
                          Tee,pe2,pi2,aio,ael)
    
    If (pe2.le.p_data%pe2min) then
       BoundsEr = 1
       return
    endif
    !BoundsEr = 0
  
    Btot2 =  Br*Br + Bz*Bz + (RBphi/r)**2
    KdK = Y(4)**2 + Y(5)**2 + (Y(6)/r)**2
    KdB = Y(4)*Br + Y(5)*Bz + Y(6)/r**2 *RBphi
    Kpar2 = KdB**2/Btot2
    Kper2 = KdK - Kpar2
  
    CALL lscsq_eps(kper2, psi, omc, pe2, pi2, tee, &
    aio, aion, ael, ipsq, p_data)
  
    d11er = -(aion + ael)
    d11w0 = ipsq+2.0_fp*aion*kper2
    d12w0 = -0.5_fp*p_data%exy
    d33w0 = p_data%epsq
  
    bkb(1) = Br*KdB/Btot2
    bkb(2) = Bz*KdB/Btot2
    bkb(3) = rBphi*KdB/Btot2 /r**2 
  
    per(1) = Y(4) - bkb(1)
    per(2) = Y(5) - bkb(2)
    per(3) = Y(6)/r**2 - bkb(3)
   
    Kper4 = Kper2**2
    Qpar = Kpar2 - p_data%woc2*p_data%Eper
    QparE= -p_data%woc2*D11er
    QparW= -p_data%woc2*(p_data%Eper+D11w0)
  
    bb = Qpar*(p_data%Epar+p_data%Eper) + p_data%woc2*p_data%Exy**2
    cc = Qpar*Qpar - p_data%woc4*p_data%Exy**2
   
    denom = Kper4*D11w0                                               &
             + Kper2*(QparW*(p_data%Epar+p_data%Eper) + Qpar*(D33w0+D11w0)&
             + p_data%woc2*p_data%Exy*p_data%Exy + 2.0_fp*p_data%woc2*p_data%EXY*D12w0 )&
             + cc*D33w0&
             + p_data%Epar*(2.0_fp*Qpar*QparW -2.0_fp*p_data%woc4*p_data%Exy *p_data%Exy&
             -2.0_fp*p_data%woc4*p_data%Exy*D12w0)  
                                               
    denom = 2.0_fp*denom
    p_data%wdDdw = denom
  
    DKper = 2.0_fp*p_data%Eper*Kper2 + bb &
              + D11er *Kper4&
              + Kper2*(QparE*(p_data%Epar+p_data%Eper) + Qpar*D11er)&
              + p_data%Epar*Qpar*QparE*2.0_fp
  
    DKpar = Kper2*(p_data%Epar+p_data%Eper) + 2.0_fp*p_data%Epar*Qpar
  
    dRdwt(1) = -2.0_fp* ( DKper*per(1) + DKpar*bkb(1) )!/denom
    dRdwt(2) = -2.0_fp* ( DKper*per(2) + DKpar*bkb(2) )!/denom
    dRdwt(3) = -2.0_fp* ( DKper*per(3) + DKpar*bkb(3) )!/denom
  
    dsdwt    = sqrt (dRdwt(1)**2 + dRdwt(2)**2 + r*r*dRdwt(3)**2)
  
    dKdwt(1) =(lscsq_DispRela(Y(1)+DR,Y(2), Y(4),Y(5),Y(6),  o_data, p_data, vecs, max_d ) -      &
               lscsq_DispRela(Y(1)-DR,Y(2), Y(4),Y(5),Y(6),  o_data, p_data, vecs, max_d  ) ) /(2.0_fp*DR)!/denom
    dKdwt(2) =(lscsq_DispRela(Y(1),Y(2)+DZ, Y(4),Y(5),Y(6),  o_data, p_data, vecs, max_d  ) -      &
               lscsq_DispRela(Y(1),Y(2)-DZ, Y(4),Y(5),Y(6),  o_data, p_data, vecs, max_d  ) ) /(2.0_fp*DZ)!/denom
    !dKdwt(3) = 0.0
  
    f(1)     = dRdwt(1)/dsdwt
    f(2)     = dRdwt(2)/dsdwt
    f(3)     = dRdwt(3)/dsdwt
    f(4)     = dKdwt(1)/dsdwt
    f(5)     = dKdwt(2)/dsdwt
    f(6)     = 0.0
    f(7)     = denom/dsdwt
  
  end subroutine lscsq_ftion

  

  end module functions

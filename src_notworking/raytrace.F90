module rTrace
  contains
  
   ! eps -- dielectric tensor              
    SUBROUTINE lscsq_eps(kper2, psi, omc, pe2, pi2, tee, &
      aio, aion, ael, ipsq, p_data)

      use iso_c_binding, only : fp => c_double
      use pl_types
      implicit none
      type(pl_data), intent(inout) :: p_data
      Real(fp) :: kper2, veow2,elam,elamc, ecyc 
      Real(fp), intent(INOUT) :: aio,ael, ipsq, Aion
      Real(fp), intent(IN) :: omc,Tee, pe2, pi2, psi
      Real(fp) :: large=1.0e30_fp
      Real(fp) :: Te2Ve=0.0445e-4_fp

      !$acc routine
    
      !call lscsq_plasma2d (o_data, p_data, r,z,psi,Br,Bz,RBphi,omc,&
      !Tee,pe2,pi2,aio,ael)
      If (psi>=large) return
          ecyc  = omc/p_data%fghz
          p_data%ecyc2 = ecyc*ecyc
          p_data%epsq  = pe2/p_data%fghz**2
          ipsq  = pi2/p_data%fghz**2
          veow2 = Te2Ve*tee/p_data%fghz**2
          elamc = veow2/p_data%ecyc2
          elam  = elamc*kper2
          ! This is the electron Lambda - - (Kper RhoE)^2
          ! Now form the quantities exp(-elam) I0,1 (elam)
          ! emli0= lscsq_bsi0(elam)
          ! emli1= lscsq_bsi1(elam)
          p_data%Epar = 1.0_fp - p_data%epsq
          aion = aio/p_data%fghz**4
          p_data%Eper = 1.0_fp+p_data%epsq/p_data%ecyc2-ipsq
        
          p_data%Eper = p_data%Eper-kper2*(aion+ael)
          p_data%Exy  = p_data%epsq/ecyc !Exy0

    end subroutine lscsq_eps

    ! ------------------------------------------------------
    !     OMPN  function for integration     

    subroutine lscsq_ftion(BoundsEr, y1, y2, y4, y5, y6, f1, f2, f3, f4, f5, f6, f7,  o_data, p_data)!, local_prmtrs)!, all_vecs)

        use pl_types 
        use iso_c_binding, only : fp => c_double
        !use lscsq_mod, only : npsij
        implicit none
        type(pl_data), intent(inout) :: p_data
        type(old_data), intent(inout) :: o_data
        !type(parameters), intent(inout) :: local_prmtrs

        type(pl_data) :: tmp_p_data
        type(old_data) :: tmp_o_data
        
        Real(fp) :: D11er, D11w0, D33w0, D12w0
        !, D33er, D12er, D11ar, D33ar, D12ar,&
                    
        Real(fp) :: denom, dkpar, dkper, d1, d2, d3, d4
        Real(fp), intent(inout) :: y1, y2, y4, y5, y6, f1, f2, f3, f4, f5, f6, f7

        integer, intent(inout) :: BoundsEr
        Real(fp) ::     dis, aion, ipsq
        Real(fp) ::     KdB,KdK,  kper2,Kper4,Kpar2
        Real(fp) ::     Btot2, Qpar, QparE,  QparW, bb, cc!, QparA
        Real(fp) ::     bkb(3), per(3), dRdwt(3), dKdwt(3),  dsdwt
        Real(fp) ::     r,z, RBphi,omc,Tee,pi2,aio,ael, &
        psi, pe2
        Real(fp) :: Br,Bz
        Real(fp) :: dr = 1.0e-4_fp      
        Real(fp) :: dz = 1.0e-4_fp     
        
        !$acc routine
        !$acc routine(lscsq_plasma2d)
        !$acc routine(lscsq_eps)
        !$acc routine(lscsq_DispRela)

        r = y1
        z = y2

        call lscsq_plasma2d(o_data, p_data, r,z,psi,Br,&
                            Bz,RBphi,omc,&
                            Tee,pe2,pi2,aio,ael)!,all_vecs)

        If (pe2 <= p_data%pe2min) then
          BoundsEr = 1
          return
        endif
        BoundsEr = 0

        Btot2 =  Br*Br + Bz*Bz + (RBphi/r)**2
        KdK = y4**2 + y5**2 + (y6/r)**2
        KdB = y4*Br + y5*Bz + y6/r**2 *RBphi
        Kpar2 = KdB**2/Btot2
        kper2 = KdK - Kpar2

        CALL lscsq_eps(kper2, psi, omc, pe2, pi2, tee, &
        aio, aion, ael, ipsq, p_data)

        d11er = -(aion + ael)
        d11w0 = ipsq+2.0_fp*aion*kper2
        d12w0 = -0.5_fp*p_data%exy
        d33w0 = p_data%epsq

        bkb(1) = Br*KdB/Btot2
        bkb(2) = Bz*KdB/Btot2
        bkb(3) = rBphi*KdB/Btot2 /r**2 

        per(1) = y4 - bkb(1)
        per(2) = y5 - bkb(2)
        per(3) = y6/r**2 - bkb(3)

        Kper4 = kper2**2
        Qpar = Kpar2 - p_data%woc2*p_data%Eper
        QparE= -p_data%woc2*D11er
        QparW= -p_data%woc2*(p_data%Eper+D11w0)

        bb = Qpar*(p_data%Epar+p_data%Eper) + p_data%woc2*p_data%Exy**2
        cc = Qpar*Qpar - p_data%woc4*p_data%Exy**2
      
        denom = Kper4*D11w0                                               &
                + kper2*(QparW*(p_data%Epar+p_data%Eper) + Qpar*(D33w0+D11w0)&
                + p_data%woc2*p_data%Exy*p_data%Exy + 2.0_fp*p_data%woc2*p_data%EXY*D12w0 )&
                + cc*D33w0&
                + p_data%Epar*(2.0_fp*Qpar*QparW -2.0_fp*p_data%woc4*p_data%Exy *p_data%Exy&
                -2.0_fp*p_data%woc4*p_data%Exy*D12w0)  
                                                  
        denom = 2.0_fp*denom
        p_data%wdDdw = denom

        DKper = 2.0_fp*p_data%Eper*kper2 + bb &
                  + D11er *Kper4&
                  + kper2*(QparE*(p_data%Epar+p_data%Eper) + Qpar*D11er)&
                  + p_data%Epar*Qpar*QparE*2.0_fp

        DKpar = kper2*(p_data%Epar+p_data%Eper) + 2.0_fp*p_data%Epar*Qpar

        dRdwt(1) = -2.0_fp* ( DKper*per(1) + DKpar*bkb(1) )!/denom
        dRdwt(2) = -2.0_fp* ( DKper*per(2) + DKpar*bkb(2) )!/denom
        dRdwt(3) = -2.0_fp* ( DKper*per(3) + DKpar*bkb(3) )!/denom
        dsdwt    = sqrt (dRdwt(1)**2 + dRdwt(2)**2 + r*r*dRdwt(3)**2)

        tmp_p_data = p_data
        tmp_o_data = o_data

        dKdwt(1) =(lscsq_DispRela((y1+DR),y2, y4,y5,y6, o_data, p_data) -      &
                  lscsq_DispRela((y1-DR),y2, y4,y5,y6, o_data, p_data) ) /(2.0_fp*DR)!/denom
        dKdwt(2) =(lscsq_DispRela(y1,(y2+DR), y4,y5,y6, o_data, p_data) -      &
                  lscsq_DispRela(y1,(y2-DR), y4,y5,y6, o_data, p_data) ) /(2.0_fp*DR)!/denom

        p_data = tmp_p_data
        o_data = tmp_o_data      

        f1     = dRdwt(1)/dsdwt
        f2     = dRdwt(2)/dsdwt
        f3     = dRdwt(3)/dsdwt
        f4     = dKdwt(1)/dsdwt
        f5     = dKdwt(2)/dsdwt
        f6     = 0.0
        f7     = denom/dsdwt

    end subroutine lscsq_ftion


      Function lscsq_disprela(r, z, Kr, Kz, Kphi,  &
        o_data, p_data, max_d)
        use iso_c_binding, only : fp => c_double
        use pl_types 

        implicit none
        type(pl_data), intent(inout) :: p_data
        type(old_data), intent(inout) :: o_data
        
        Real(fp) :: d1,d2,d4 
        Real(fp), optional, intent(inout) :: max_d
        Real(fp) :: lscsq_disprela

        Real(fp) ::   psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio, &
        aion, ael, ipsq, Kpar2, kper2, Kdk
        Real(fp),  intent(in) ::  r,z,Kr, Kz, Kphi
        Real(fp) ::  Btot2, Qpar

        !$acc routine
        !$acc routine(lscsq_plasma2d)
        !$acc routine(lscsq_eps)

        KdK  = Kr**2 + (Kz)**2 + (KPhi/r)**2
        call lscsq_plasma2d(o_data, p_data, r,z,psi,Br,&
                            Bz,RBphi,omc,&
                            Tee,pe2,pi2,aio,ael)

        Btot2 =  Br**2 + Bz**2 + (RBphi/r)**2
        Kpar2 = ( Kr*Br + Kz*Bz + Kphi*rBphi/r/r )**2/Btot2
        kper2 = KdK - Kpar2

        CALL lscsq_eps (kper2, psi, omc, pe2, pi2, &
        tee, aio, aion, ael, ipsq, p_data)

        Qpar = Kpar2 - p_data%woc2*p_data%Eper
        D1 = kper2**2*p_data%Eper
        D2 = kper2*( (p_data%Epar+p_data%Eper)*Qpar + p_data%woc2*p_data%Exy**2 )
        D4 = p_data%Epar*( Qpar**2 - p_data%woc4*p_data%Exy**2 )

        if (present(max_d)) max_d = max(abs(d1),abs(d2),abs(d4))  
        lscsq_disprela  = D1 + D2 + D4

    END function lscsq_disprela
    !                             
    !                                                                      
    !     DispRela ends                         

    ! ------------------------------------------------------      

end module rTrace

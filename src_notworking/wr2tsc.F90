subroutine lscsq_PJrfIgrl
  ! Computes the integrals of Prf, Pql, Jrf, as function of psi
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : pi, twopi
  use lscsq_mod, only: npsi
  use lscsq_mod, only: neary
  use lscsq_mod, only: jrundot, nrundot, js, jsp,IrIntgrl
  use lscsq_mod, only: IpIntgrl,PqIntgrl,PqlTot
  use lscsq_mod, only: PrIntgrl, praytot, dvol
  implicit none

  integer :: ips
  real(fp), dimension(npsi) :: nRdotNorm, jRdotNorm
 
  ! Integrate Pray and Pql in LSC flux surfaces
  PrIntgrl(1) = PRayTot(1)
  PqIntgrl(1) = PqlTot (1)
  do ips = 2, npsi
     PrIntgrl(ips) = PrIntgrl(ips-1) + PRayTot(ips)
     PqIntgrl(ips) = PqIntgrl(ips-1) + PqlTot(ips)
  enddo

  ! Integrate Jrf and Jrf-plus in LSC flux surfaces
  ! Divide by 2 pi R to make an area element !
  ! the integral is multiplied back to 2 pi R at teh exit from this subroutine
  ! and these variables are used only to calculate the current profile by taking
  ! the difference between adjacent zones. We can probably remove the divide by
  ! here and the multiply by in lscsq_output.
  IrIntgrl(1) = js(1)*dVol(1)
  IpIntgrl(1) = jsp(1)*dVol(1)
  do ips = 2, npsi
     IrIntgrl(ips) = IrIntgrl(ips-1) + js(ips)*dVol(ips)
     IpIntgrl(ips) = IpIntgrl(ips-1) +jsp(ips)*dVol(ips)
  enddo
!  IrIntgrl = IrIntgrl/(twopi*xmag)
!  IpIntgrl = IpIntgrl/(twopi*xmag)

  do ips=1,Npsi
     if (PrIntgrl(ips) .LE. 0.05_fp*PrIntgrl(Npsi) .or.              &
        PrIntgrl(ips) .GE. 0.95_fp*PrIntgrl(Npsi)) then
        nRunDot(ips) = 0.00
        jRunDot(ips) = 0.00
     endif
     nRdotNorm(ips) = nRunDot(ips)/(NeAry(ips) + 1.e10_fp)
     jRdotNorm(ips) = jRunDot(ips)/(   js(ips) + 1.e01_fp)
  enddo

end subroutine lscsq_PJrfIgrl
!------------------------------------------------------------ 
subroutine lscsq_output(Pelfnd,Jrffnd, Pqlfnd)
  use grap
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : zero, pi, twopi, zm3tocm3, zm2tocm2
  use lscsq_mod, only : zcm3tom3, zcm2tom2
  use lscsq_mod, only: edcvec
  use lscsq_mod, only: IpIntgrl,IrIntgrl,PqIntgrl,PrIntgrl
  use lscsq_mod, only: midvec,midary,npsij
  use lscsq_mod, only: dEdcAmnt,dVlVec
  use lscsq_mod, only: npsi, powtsc, curtsc, dJdE,dlJdlE
  implicit none

  integer ::  l
  real(fp)::  PelFnd, JrfFnd, PqlFnd!, PioFnd
  real(fp)::  xlookup, yreturn  
 
  real(fp), dimension(:), allocatable :: curtscp,powdql, powtsci, curtsci, curtscIp, powdqli

!------------------------------
!  dmc -- interpolation repairs (axial singularity)
  integer :: inside
  real(fp) :: zvin
!------------------------------
  !integer :: idbg = 0

  if(.not.allocated(curtscp)) allocate(curtscp(npsij))
  if(.not.allocated(curtsci)) allocate(curtsci(npsij))
  if(.not.allocated(powtsci)) allocate(powtsci(npsij))
  if(.not.allocated(curtscip)) allocate(curtscip(npsij))
  if(.not.allocated(powdql)) allocate(powdql(npsij))
  if(.not.allocated(powdqli)) allocate(powdqli(npsij))

  powtsc(1:npsij) = zero
  powtsci(1:npsij) = zero
  curtsc(1:npsij) = zero
  curtscp(1:npsij) = zero
  curtsci(1:npsij) = zero
  curtscip(1:npsij) = zero
  powdql(1:npsij) = zero
  powdqli(1:npsij) = zero
  dlJdlE(1:npsij) = zero
  dJdE(1:npsij) = zero

  CALL lscsq_PJrfIgrl

  !dmc fixup block
  inside=0
  zvin=0
  !end dmc fixup block

  do l = 1, NpsiJ-1
     xlookup =  MidVec(l)
     !dmc fixup block
     if(xlookup.lt.MidAry(1)) then
        inside=l
        zvin=zvin+dVlVec(l)
        cycle
     endif
     !end dmc fixup block
 
     ! interpolate for power integral
     CALL lscsq_linr1d(Npsi, MidAry, PrIntgrl, xlookup, yreturn)
     powtscI(l) =  abs(yreturn) ! power in watts inside the l flux surface boundary
 
     ! interpolate for current integral
     CALL lscsq_linr1d(Npsi, MidAry, IrIntgrl, xlookup, yreturn)
     curtscI (l)= yreturn  ! current in amps inside the l flux surface boundary

     ! interpolate for current integral at enhanced E [ E plus dE]
     CALL lscsq_linr1d(Npsi, MidAry, IpIntgrl, xlookup, yreturn)
     curtscIp(l)= yreturn
 
     ! interpolate for the ql power deposited
     CALL lscsq_linr1d(Npsi, MidAry, PqIntgrl, xlookup, yreturn)
     powDqlI(l) = yreturn
  enddo

  powtscI (NpsiJ) = PrIntgrl(Npsi)
  curtscI (NpsiJ) = IrIntgrl(Npsi)
  curtscIp(NpsiJ) = IpIntgrl(Npsi)
  powDqlI (NpsiJ) = PqIntgrl(Npsi)

!  powtscI (NpsiJ-1) = PrIntgrl(Npsi)
!  curtscI (NpsiJ-1) = IrIntgrl(Npsi)
!  curtscIp(NpsiJ-1) = IpIntgrl(Npsi)
!  powDqlI (NpsiJ-1) = PqIntgrl(Npsi)

!dmc fixup block -- flatten current and power inside innermost grid
  if(inside.gt.0) then
     zvin=zvin+dVlVec(inside+1)
     powtscI(1)=powtscI(inside+1)*dvlvec(1)/zvin
     powDqlI(1)=powDqlI(inside+1)*dvlvec(1)/zvin
     curtscI(1)=curtscI(inside+1)*dvlvec(1)/zvin
     curtscIp(1)=curtscIp(inside+1)*dvlvec(1)/zvin
     do l=2,inside
        powtscI(l)=powtscI(l-1)+powtscI(inside+1)*dvlvec(l)/zvin
        powDqlI(l)=powDqlI(l-1)+powDqlI(inside+1)*dvlvec(l)/zvin
        curtscI(l)=curtscI(l-1)+curtscI(inside+1)*dvlvec(l)/zvin
        curtscIp(l)=curtscIp(l-1)+curtscIp(inside+1)*dvlvec(l)/zvin
     enddo
  endif
!end dmc fixup block

  !     divide by volume for power, area for current.
  !     Note backward volume convention

!  powtsc(1) = powtscI(1)/dVlVec(1)
!  powDql(1) = powDqlI(1)/dVlVec(1)
!  curtsc(1) = curtscI(1)/dVlVec(1)
!  curtscp(1)= curtscIp(1)/dVlVec(1)

  do l = 2, NpsiJ
     powtsc(l) = (powtscI(l)-powtscI(l-1))/dVlVec(l)
     powDql(l) = (powDqlI(l)-powDqlI(l-1))/dVlVec(l)
     curtsc(l) = (curtscI(l)-curtscI(l-1))/dVlVec(l)
     curtscp(l)= (curtscIp(l)-curtscIp(l-1))/dVlVec(l)
  enddo

  powtsc(1)=powtsc(2)
  powDql(1)=powDql(2)
  curtsc(1)=curtsc(2)
  curtscp(1)=curtscp(2)

  powtsc(1:npsij) = zcm3tom3*powtsc(1:npsij)
  powDql(1:npsij) = zcm3tom3*powDql(1:npsij)
  ! I have to mult by 2 pi R to make an area! (comment by DMC)
  ! fmp - I do not like this multiplication, since we are using a volume
  ! integration by flux surface and multiplying by a constant major radius
!  curtsc(1:npsij) = zcm2tom2*curtsc(1:npsij)*(twopi*xmag)
!  curtscp(1:npsij) = zcm2tom2*curtscp(1:npsij)*(twopi*xmag)
  curtsc(1:npsij) = zcm2tom2*curtsc(1:npsij)
  curtscp(1:npsij) = zcm2tom2*curtscp(1:npsij)

  do l = 1, NpsiJ
     if(curtsc(l).eq.zero .or. dEdcAmnt.eq.zero .or. EdcVec(l).eq.zero) then
        dlJdlE(l) = zero
     else
        dlJdlE(l) = (curtscp(l)-curtsc(l))/curtsc(l)*EdcVec(l)/dEdcAmnt
     endif
     ! TRANSP REQUEST:  
     dJdE(l) = (curtscp(l)-curtsc(l))/dEdcAmnt
  enddo
  
  ! powtsc is in watts/cc
  ! curtsc is in amps/cm**2
  ! powtscI is in watts
  ! curtscI is in amps

  PelFnd = powtscI(NpsiJ)
  JrfFnd = curtscI(NpsiJ)
  PqlFnd = powDqlI(NpsiJ)

  call lscsq_FastFrac  ! Call FastFrac every time, but write it out only if the proper flag is set.
                     ! New FastFrac logic at request of TRANSP

end subroutine lscsq_output

!     ------------------------------------------------------------------
subroutine lscsq_FastFrac
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : zero
  use lscsq_mod, only: vpar, vtherm,fe,ivzero,iitr, dvplus
  use lscsq_mod, only: fenorm,FstFracN,FstFracE
  use lscsq_mod, only: npsi, nv
  implicit none

      INTEGER i,j, jwt
      REAL(fp)    ExpMax, RsltMin
      REAL(fp)    duNorm, duVth2, duFracN, duFracE, duMaxwN, duMaxwE
      REAL(fp)    dumFe, dumMx, dumV2, duDelV, exp1
      REAL(fp)    wt(2)
      DATA         wt(1)          ,       wt(2)             /           &
     &     0.666666666666d0, 1.333333333333d0 /
      DATA    ExpMax / 100.0d0 /

!     FastFrac: Fast Fraction of N FstFracN and of Energy FstFracE
!     are computed and writtend in Subroutine FastFrac(npsi)
!     du        prefix meaning dummy
!     duNorm    Normalization constant for this flux surface
!               such that f_e(v) = duNorm * exp(-v^2/(2 v_T^2) for
!               Maxwellian
!     duVth2    v_T^2 normalizing the exponent for this flux surface
!     duDelV    The width of the velocity grid at the center
!     duFracN   The fraction of particles in the heated portion of the tail
!     duFracE   The fraction of energy in the heated portion of the tail
!     duMaxwN   "particles" in Maxwellian; (local units not relevant)
!     duMaxwE   "energy" in Maxwellian; (local units not relevant)
!     FeNorm    array dimensioned NPSIDIM containing NeAry/(2 PI)/Vtherm(ip)
!
!     Simpsons Rule is (f1 + 4f2 + f3)h/3 or
!                     (2f1 + 4f2 + 2f3 + 4f4 + 2f5 )h/3 less
!                   - ( f1                      f5 )h/3 which can be ignored
!     because the end points are exponentially small.  The weight (wt) is
!     from the expression  [ 2 + 2( j / 2 ) - j ]  = 1, {j=1,3,5...}
!                                                  = 2, {j=2,4,6...}
!
      duDelV  = abs ( Vpar(IvZero+1) - Vpar(IvZero) )
      RsltMin = exp ( -ExpMax )
      do 20 i = 1, npsi
        duFracN = ZERO
        duFracE = ZERO
        duMaxwN = ZERO
        duMaxwE = ZERO
        duNorm  = FeNorm(i)
        duVth2  = Vtherm(i)*Vtherm(i)
        if ( vtherm(i) .GE. duDelV ) then
          do 10 j = 1, nv
            jwt =  ( 2 + 2 * ( j / 2 ) - j )
!            dumFe  =              fe(j,i,iITR) * wt(jwt) * dvsym(j)
!            dumV2  = vpar(j)**2 * fe(j,i,iITR) * wt(jwt) * dvsym(j)
            dumFe  =              fe(j,i,iITR) * wt(jwt) * dvplus(j)
            dumV2  = vpar(j)**2 * fe(j,i,iITR) * wt(jwt) * dvplus(j)
            duFracN = duFracN + dumFe
            duFracE = duFracE + dumV2

            exp1 =  vpar(j)*vpar(j)/(2.*duVth2)
            if(exp1 .LT. ExpMax) then
              exp1 = exp( -exp1 )
            else
              exp1 = RsltMin
            endif
            exp1 = exp1 * duNorm
!           exp1 = exp( -vpar(j)*vpar(j)/(2.*duVth2) )* duNorm
!             dumMx  =              exp1 * wt(jwt) * dvsym(j)
!             dumV2  = vpar(j)**2 * exp1 * wt(jwt) * dvsym(j)
              dumMx  =              exp1 * wt(jwt) * dvplus(j)
              dumV2  = vpar(j)**2 * exp1 * wt(jwt) * dvplus(j)
              duMaxwN = duMaxwN + dumMx
              duMaxwE = duMaxwE + dumV2
 10       continue

          FstFracN(i) = ( duFracN - duMaxwN ) / duMaxwN
            FstFracE(i) = ( duFracE - duMaxwE ) / duMaxwE
        else
          FstFracN(i) = ZERO
          FstFracE(i) = ZERO
        endif
 20   continue
!      if (PrFlg(2) .GE. 1.0)  then
!!       CALL lscsq_LSCpause
!        write(nLSCcomm,1000)
! 1000       format (                                                    &
!     &'  Psi Index   Fraction Tail Particles  Fraction Tail Energy',/)
!         write(nLSCcomm,1001)                                           &
!     &                     (  i, FstFracN(i), FstFracE(i)  , i=1,npsi)
! 1001       format(t8, i4, t28, 1pe10.3, t50, 1pe10.3)
!!       CALL lscsq_LSCpause
!      endif
      return
END
!
!     ------------------------------------------------------------------

module Tracing
    contains
  SUBROUTINE lscsq_LSCstop ( ErrMsg )
    use iso_c_binding, only : fp => c_double
    use lscsq_mod, only: nTSCscrn
    use lscsq_mod, only: ierror
    implicit none
        CHARACTER*(*) ErrMsg
        !acc routine seq
       ! write(nTSCscrn,'('' LSCstop: '',a40 )')                           &
      !&                     ErrMsg

        write(nTSCscrn,*) ErrMsg
        iError = iError + 1
        return
      END  SUBROUTINE lscsq_LSCstop     
  !     ------------------------------------------------------------------
  SUBROUTINE lscsq_LSCwarn ( ErrMsg )
    use iso_c_binding, only : fp => c_double
    use lscsq_mod, only:  iRayTrsi!, iendry
    implicit none

    character(len=*) :: ErrMsg
    integer, save :: nwarning = 0
    integer, save :: maxWarn  = 1000

    !acc routine seq
    !acc routine(lscsq_LSCstop)
    
    !integer :: i, ilen
    !ilen = len(trim(ErrMsg))
    !write(0,*) ilen, ErrMsg

    nWarning = nWarning + 1
    if (nWarning .GT. maxWarn ) then
      ! write(nTSCscrn,'('' LSCquit: '',70a1 )') (ErrMsg(i:i), i=1,ilen)
      write(0,*) 'LSCquit:', ErrMsg
      CALL lscsq_LSCstop ( ' too many warnings issued! ' )
      nWarning = 0
    else if (iRayTrsi .ne. 0) then
      !write(nTSCscrn,'('' LSCwarn: '',70a1 )') (ErrMsg(i:i), i=1,ilen)
      write(0,*) 'LSCquit:', ErrMsg
    else
      ! Do not issue warnings if only re-doing the J of E calculation:  iRayTrsi=0
      return
    endif
  end subroutine lscsq_LSCwarn
  !     ------------------------------------------------------------------
  SUBROUTINE lscsq_LSCendr ( ErrMsg, lstop)
  !     Trys to end ray tracing asap because an error is encountered
    use iso_c_binding, only : fp => c_double
    use lscsq_mod, only: iendry
    !use lscsq_mod, only: lstop
    implicit none
    integer, intent(INOUT):: lstop
    CHARACTER*(*) ErrMsg
    !acc routine seq
    !acc routine(lscsq_LSCwarn)
        
        Lstop = 1
        iEndRy = iEndRy + 1
        CALL lscsq_LSCwarn ( ErrMsg )
        return
  !     ------------------------------------------------------------------
  !      ENTRY lscsq_LSCbigK ( ErrMsg )
  !      Lstop = 1
  !      CALL lscsq_LSCwarn ( ErrMsg )
  !      return

  END SUBROUTINE lscsq_LSCendr 
  !     ------------------------------------------------------------------
  SUBROUTINE lscsq_LSCtrace ( ErrMsg )
    use iso_c_binding, only : fp => c_double
    use lscsq_mod, only: nTSCscrn
    implicit none
        CHARACTER*(*) ErrMsg
        write(nTSCscrn,'('' LSCtrace called on exiting: '',a40 )')        &
      &                     ErrMsg
        return
  END SUBROUTINE lscsq_LSCtrace

  

  
  end module Tracing
!
!     ------------------------------------------------------------------

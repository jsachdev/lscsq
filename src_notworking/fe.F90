SUBROUTINE lscsq_FeCalc
! compute electron distribution function, given quasilinear
! diffusion coefficient and plamsa profiles from psi index ipsiL to ipsiU

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: npsi,nv,ivzero,vpar,dql, iitr
  use lscsq_mod, only: nucoll, dcoll, fenorm, fe
  implicit none

  integer :: ip
  
  ! set to zero to start
  fe(1:nv,1:npsi,iitr) = 0.0

  ! copy normalization into fe(v = vpar(ivZero), psi)
  ! fmp - I think there is a problem here. The values of Fe at v=0 are huge in
  ! the output. 
  fe(ivzero,1:npsi,iitr) = Fenorm(1:npsi)

   do ip = 1, npsi
     ! solve for positive velocity
     CALL lscsq_FePlus(fe(1:nv,ip,iitr), nuColl(1:nv, ip), Dcoll(1:nv,ip), Dql(1:nv,ip,2), Vpar, ivZero)
     ! solve for negative velocity
     CALL lscsq_FeMinus(fe(1:nv,ip,iitr), nuColl(1:nv,ip), Dcoll(1:nv,ip), Dql(1:nv,ip,2), Vpar, ivZero)
   enddo

end subroutine lscsq_FeCalc
!     -----------------------------------------------------------------
subroutine lscsq_FeMkNorm
!     determine normalized value of fe at ivZero vs psi
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: npsi,tailneps,tailteps
  use lscsq_mod, only: vtherm, fe0, fenorm, neary
  implicit none

  integer :: ip
  real(fp)  :: TailFact

  TailFact = (1.0_fp+TailNeps*sqrt(TailTeps))/(1.0_fp+TailNeps)
  do ip = 1, npsi
     FeNorm(ip) = fe0*NeAry(ip)/Vtherm(ip)*TailFact
  enddo

end subroutine lscsq_FeMkNorm
!     -----------------------------------------------------------------
subroutine lscsq_FePlus(fe, nuCollx, Dcollx, Dqlx, Vpar, ivZero)
  ! solve for fe for positive velocity, using the convention that
  ! the new value is the old value plus the mean integrand
  ! times the delta-v to the new value
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: nv
  implicit none

  integer :: ivZero, iv
  real(fp):: RsltMin
  real(fp), dimension(nv), intent(inout) :: fe
  real(fp), dimension(nv), intent(in) :: nucollx
  real(fp), dimension(nv), intent(in) :: Dcollx
  real(fp), dimension(nv), intent(in) :: Dqlx
  real(fp), dimension(nv), intent(in) :: vpar
  real(fp) :: expmax=100.0_fp
  real(fp), dimension(nv)  :: lgfe


  lgFe(ivZero) = 0.0_fp  
  RsltMin = exp(-ExpMax)
  do iv = ivZero, nv-1
     lgFe(iv+1)=lgFe(iv)+                                       &
                0.5_fp*(nuCollx(iv+1)/(Dqlx(iv+1)+Dcollx(iv+1))+  &
                      nuCollx(iv)/(Dqlx(iv)+Dcollx(iv)))*       &
                     (Vpar(iv+1)-Vpar(iv))
  enddo
  
  do iv = ivZero+1, nv-2
     if(lgFe(iv) .LT. ExpMax) then
        fe(iv) = fe(ivZero)*exp(-lgFe(iv))
     else
        fe(iv) = fe(ivZero)*RsltMin
     endif
   enddo
   fe(nv-1) = 0.0_fp  
   fe(nv)   = 0.0_fp   

end subroutine lscsq_feplus

!     -----------------------------------------------------------------

subroutine lscsq_FeMinus(fe, nuCollx, Dcollx, Dqlx, Vpar, ivZero)
!     solve for fe for negative velocity, using the convention that
!     the new value is the old value plus the the mean integrand
!     times the delta-v to the new value
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: nv
  implicit none

  real(fp), dimension(nv), intent(inout) :: fe
  real(fp), dimension(nv), intent(in) :: nucollx
  real(fp), dimension(nv), intent(in) :: Dcollx
  real(fp), dimension(nv), intent(in) :: Dqlx
  real(fp), dimension(nv), intent(in) :: vpar
  integer :: ivZero, iv
  real(fp):: RsltMin
  real(fp), dimension(nv) :: lgFe
!                                       solve for fe for negative velocity
  real(fp) :: expmax=100.0_fp

  ! note:  nu < 0 for iv < ivZero
  lgFe(ivZero) = 0.0_fp   
  RsltMin = exp (-ExpMax)
  do iv = ivZero-1, 1, -1
     lgFe(iv) = lgFe(iv+1) +                                      &
                0.5_fp*(nuCollx(iv+1) / (Dqlx(iv+1)+Dcollx(iv+1)) + &
                      nuCollx(iv)/(Dqlx(iv)+Dcollx(iv)))*  &
                (Vpar(iv)-Vpar(iv+1))
  enddo

  do iv = ivZero-1, 3, -1
     if( lgFe(iv) .LT. ExpMax) then
        fe(iv) = fe(ivZero)*exp(-lgFe(iv))
     else
        fe(iv) = fe(ivZero)*RsltMin
     endif
  enddo

  fe(1) = 0.0_fp  
  fe(2) = 0.0_fp  

end subroutine lscsq_FeMinus

!     -----------------------------------------------------------------

subroutine lscsq_FeWeight
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: npsi, nv, weghtitr, fe, iitr
  implicit none

  INTEGER iv,ip,new,old
  real(fp) :: CtrWeght

  CtrWeght = 1.0_fp - WeghtItr
  new = iITR
  old = mod(iITR,2)+1
  do iv = 1, nv
     do ip = 1, npsi
        fe(iv,ip,new) = fe(iv,ip,new)*WeghtItr + fe(iv,ip,old)*CtrWeght
     enddo
  enddo

end SUBROUTINE lscsq_FeWeight

!     -----------------------------------------------------------------

subroutine lscsq_FeInit
  ! initialize fe solver, fe
  use iso_c_binding, only : fp => c_double
  implicit none

  integer :: ifirstcall = 1

  if (ifirstcall.EQ.1) then
     ifirstcall = 0
     CALL lscsq_FeConst
  endif

  CALL lscsq_FeArrays

end subroutine lscsq_FeInit
 
!     -----------------------------------------------------------------

subroutine lscsq_FeArrays

  CALL lscsq_mkvth     ! generate thermal velocity array (vthermal AND Vperpsq vs psi)
  CALL lscsq_FeMkNorm  ! array of normalization for Fe (Fe(v = 0, psi))
  CALL lscsq_Fecvecs   ! initialize collisional vectors

end subroutine lscsq_FeArrays
!     -----------------------------------------------------------------
SUBROUTINE lscsq_FeCvecs
  !    initialize collisional diffusion and drag terms Dcoll, nuColl with velocity bins and psi bins.
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: tailteps, tailneps, tailpeps, tailvtrn
  use lscsq_mod, only: nu0psi, nucoll, dcoll, vtherm
  use lscsq_mod, only: neary,lnlary,betzary, npsi, nv, vpar
  use lscsq_mod, only: dcollnorm, nunorm
  implicit none

  INTEGER :: ip, iv
  real(fp) :: vthsq, vth3, vth5, harg, hvpar, vpnorm, v12
  real(fp) :: TransVel, LogThing, TailT32, TailT12

  TailT12 = sqrt (TailTeps)
  TailT32 = sqrt (TailTeps) * TailTeps

  if (TailNeps * TailTeps .GT. 0.0_fp ) then
     LogThing = (-log(TailNeps) - 0.5_fp*log(TailTeps)) /(1.0_fp-TailTeps)
  else
     LogThing = 1.0e06_fp
  endif
  TailVtrn = sqrt(2.0_fp*LogThing)
  do ip = 1, npsi
     vthsq = vtherm(ip)**2
     TransVel = vtherm(ip)*TailVtrn
     vth3 = vthsq*vtherm(ip)
     vth5 = vth3 * vthsq
     nu0psi(ip) = NeAry(ip)/vth3
     do iv = 1, nv
        v12 = vpar(iv)
        vpnorm = (v12 / vtherm(ip))
        harg = 1.0_fp + vpnorm**2
        hvpar = 1.0_fp/(harg*sqrt(harg))

        Dcoll(iv, ip) = DcollNorm * nu0psi(ip) * vthsq * hvpar * LnlAry(ip) * BetZAry(ip)
        nuColl(iv, ip) = nuNorm * nu0psi(ip) * hvpar * v12 * LnlAry(ip) * BetZAry(ip)
 
        if ( TailPeps .GT. 0.0 .and. abs(vpar(iv)) .GT. TransVel) then
           Dcoll (iv,ip) = Dcoll(iv,ip) * TailT12
           nuColl(iv,ip) = nuColl(iv,ip)* TailT32
        endif
     enddo
  enddo
      
end subroutine lscsq_FeCvecs 

!     -----------------------------------------------------------------

subroutine lscsq_mkVth
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: npsi, teary,vtherm,vperpsq,vthnorm
  implicit none

  integer :: ip

  do ip = 1, npsi
     Vtherm(ip) = VthNorm*sqrt(TeAry(ip))
     VperpSq(ip) = Vtherm(ip)**2
  enddo

end subroutine lscsq_mkVth

!     -----------------------------------------------------------------

subroutine lscsq_FeAt0
  ! value of fe at t = 0.
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: npsi
  implicit none

  integer :: i

  ! initial conditions for fe(v, psi)
  do i = 1, npsi
     CALL lscsq_FeMaxw(i)   ! set fe at iipsi vs velocity
  enddo

end subroutine lscsq_FeAt0
!     -----------------------------------------------------------------
subroutine lscsq_FeMaxw(iipsi)
  ! contruct Maxwell distribution at Psi(iipsi)
  ! set fe at iipsi vs velocity
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : pi, twopi
  use lscsq_mod, only: tailteps, tailneps, nv
  use lscsq_mod, only: vtherm, vpar, fenorm, fe, iitr, ivzero
  implicit none

  integer, intent(in) :: iipsi
  integer ::  iv
  real(fp) :: VthsInv, argument
  real(fp) :: Factor, FactorM1, TransVel, RsltMin !, REXPMIN
  real(fp) :: LoVfac, HiVfac, LoVterm, HiVterm

  real(fp) :: toobig=200.0_fp
  real(fp) :: rexpmin=100.0_fp
!
!  fmp - Oct 2021 -- maybe it is time to replace this? the old fix was to
!  replace 100 with 85 (???)
!  dmc 13 Jan 1995 -- replaced "100." with REXPMIN in a couple of
!  places; this avoids compiler problem under osf/1 v3.0, f77 v3.6:
!  evaluation of exp( -100. ) gave an arithmetic error.  In the TRANSP
!  environment, however, it should be noted that exp( -REXPMIN )
!  evaluates by silent underflow to zero.
!
!     produce Maxwellian fe(ipsi) normalized to unity (MKS units)
!     TailPeps  Regarding a ficticious fast electron tail in both
!     TailNeps  directions, the fraction eps(ilon) for Pressure, Ne and Te
!     TailTeps  such that
!     .                  TailNeps = TailPeps * TailTeps
!     .                  TailTeps = T_thermal/T_fast
!     .                  TailPeps = (n_fast T_fast)/(n_thermal T_thermal)
!     .         and
!     .              f(v) = (2 pi v_t^2)^-.5 1. n_e exp[- v^2/(2 v_t^2)]; all v
!     .         f_fast(v) = (2 pi v_f^2)^-.5 1. n_f exp[- v^2/(2 v_f^2)]; all v
!     .
!     .         where n_f/n_e = TailNeps
!     .               v_t/v_f = TailTeps^.5
!     .               n_f v_f^2 / [ n_t v_t^2 ] = TailPeps
!     .
!     TailVtrn  is the transition velocity relative to v_t at which the fast
!               electron tail becomes more important
!                        TailVtrn = v/v_t | transition
!     .                           = sqrt[2 ln(1/Neps 1/Teps^.50)/(1 - Teps^2)]
!

  VthsInv = 0.5_fp/Vtherm(iipsi)**2
  Factor  = 1.0_fp + sqrt(TailTeps)*TailNeps
  FactorM1= sqrt(TailTeps)*TailNeps
 
  if (FactorM1 .EQ. 0.0) then
     do iv = 3, nv-2
        argument = VthsInv*vpar(iv)**2
        if (argument .LT. TOOBIG) then
          fe(iv,iipsi,iITR) = FeNorm(iipsi)*exp(-argument)
        else
          fe(iv,iipsi,iITR) = 0.0_fp 
        endif
     enddo
 
     fe(1   ,iipsi,iITR) = 0.0_fp 
     fe(2   ,iipsi,iITR) = 0.0_fp 
     fe(nv-1,iipsi,iITR) = 0.0_fp 
     fe(nv  ,iipsi,iITR) = 0.0_fp 
  endif
 
  if (FactorM1 .ne. 0.0 ) then
     RsltMin = exp ( - REXPMIN )
     TransVel = vtherm(iipsi)*sqrt(2.0_fp*(-log(TailNeps)-0.5_fp*log(TailTeps))/(1.0_fp-TailTeps))
      
     fe(ivZero,iipsi,iITR) = 0.0_fp
     LoVfac = VthSInv
     HiVfac = LoVfac*TailTeps
     do  iv = ivZero - 1, 3, - 1
        if(abs(vpar(iv)) .LT. TransVel) then
          HiVterm = vpar(iv)*LoVfac
        else
          HiVterm = vpar(iv)*HiVfac
        endif
        if(abs(vpar(iv+1)) .LT. TransVel) then
          LoVterm = vpar(iv+1)*LoVfac
        else
          LoVterm = vpar(iv+1)*HiVfac
        endif
        fe(iv,iipsi,iITR) = fe(iv+1,iipsi,iITR)+(HiVterm + LoVterm)*(Vpar(iv)-Vpar(iv+1))
     enddo
 
     do  iv = ivZero + 1, nv-2,  1
        if(abs(vpar(iv)) .LT. TransVel) then
          HiVterm = vpar(iv)*LoVfac
        else
          HiVterm = vpar(iv)*HiVfac
        endif
        if(abs(vpar(iv-1)) .LT. TransVel) then
          LoVterm = vpar(iv-1)*LoVfac
        else
          LoVterm = vpar(iv-1)*HiVfac
        endif
        fe(iv,iipsi,iITR) = fe(iv-1,iipsi,iITR) + (HiVterm+LoVterm)*(Vpar(iv)-Vpar(iv-1))
     enddo
     do iv = 3, nv-2
        if( fe(iv,iipsi,iITR) .LT. REXPMIN ) then
          fe(iv,iipsi,iITR) = FeNorm(iipsi)*exp(-fe(iv,iipsi,iITR))
        else
          fe(iv,iipsi,iITR) = FeNorm(iipsi) * RsltMin
        endif
      enddo
      fe(1,iipsi,iITR) = 0.0_fp 
      fe(2,iipsi,iITR) = 0.0_fp 
      fe(nv-1,iipsi,iITR) = 0.0_fp 
      fe(nv  ,iipsi,iITR) = 0.0_fp 
    endif

end subroutine lscsq_FeMaxw
!     -----------------------------------------------------------------
subroutine lscsq_FeConst
!                                       constants needed in solution of fe
!     In the Valeo-Eder paper
!     \nu_0 = \beta_z lnLambda 4 \pi e^4 n /( m^2 v_t^3)
!
!     where \beta_z is of order 1/2 and is about (1+Z)/5
!     where lnLambda is the Coulomb logarithm
!     so
!     \nu_0 = \beta_z lnLambda \cdot
!             TWOPI/5  1.6^4  3.0^1 / 9.11^2   n_{13}/(v/c)^3

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : pi, twopi, vc, me_g, me_Kg, zel
  use lscsq_mod, only : zkeV2eV, zmtocm, zcm3tom3
  use lscsq_mod, only: vthnorm,fe0,nu0,pwrnorm,nunorm,dcollnorm
  implicit none

  VthNorm = sqrt(1.6e-12_fp*zkev2eV/me_g)/(vc*zmtocm)       

  fe0 = 1.0_fp/sqrt(twopi)
  nu0 = 4.0_fp*pi*zel**4/me_g**2

  PwrNorm = 1.0e6_fp*me_Kg*vc**2

!     normalization for computation of QL power deposition
!     (see files [FW]power.F, subroutine RFheat[Ele])
!     The heating rate is:
!     3/2 n dT/dt = \int S_w \p \eps/\p v dv^3
!     where S_w is the wave induced flux in velocity space, and
!     where \eps is the energy per particle  == 1/2 mv^2
!     S_w = - D_{QL} \cdot \p f/\p v
!     so heating  = \int D_{QL} mv df/dv dv^1 (integrating over v-perps)
!                 = \int mc^2 (D_{QL} v d(cf)/dv dv (normalizing D and v to c)
!     The expression for PwrNorm gives mc^2 and then a conversion to m^-3
!     because f has density in it at cm-3 units.
 
  nuNorm =    zcm3tom3*nu0/vc**3
  DcollNorm = nuNorm

end subroutine lscsq_Feconst
! -----------------------------------------------------------------
subroutine lscsq_FePrime
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: npsi, nv, nsmoo, dfdv, vpar, qlsm, ivzero
  implicit none

  integer :: ips, iv
 
  call lscsq_FePru
  ! compute unsmoothed derivative

  do ips = 1, npsi
     do iv = 1, nv
        dfdv(iv,ips,2) = vpar(iv)*dfdv(iv,ips,1)
     enddo
  enddo

  call lscsq_smooth(dfdv(1:nv,1:npsi,2), nv, npsi, nsmoo, qlsm)

!     convolve with smoothing function
!
!     NOTE: For energy conservation, the smoothing done here must
!           be identical to that done in constructing Dql.
!           Conceptually equivalent to inclusion of resonance broadening
!           in wave-particle interaction. Loops added at 5,14,15,16
!           make smoothing over v df/dv, rather than over df/dv alone.
 
  do ips = 1, npsi
     do iv = 1, ivZero-1
        dfdv(iv,ips,2) = dfdv(iv,ips,2)/vpar(iv)
     enddo        
     dfdv(ivZero,ips,2) = 0.0_fp
     do iv = ivZero+1 , nv
        dfdv(iv,ips,2) = dfdv(iv,ips,2)/vpar(iv)
     enddo        
  enddo        
end subroutine lscsq_feprime
 
subroutine lscsq_FePrU
  ! Fe Prime Unsmoothed
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only: npsi, nv, vpar
  use lscsq_mod, only: ivzero, fe, dfdv, iitr ! dvip
  implicit none

  integer :: ips, iv
 
  do ips = 1, npsi
     do iv = ivZero+1, nv-1
        dfdv(iv,ips,1)=(fe(iv+1,ips,iITR)-fe(iv,ips,iITR))/(vpar(iv+1)-vpar(iv))
      enddo
     ! v grid symmetry make dfdv by looking out to high abs(v) 
     do iv = 2, ivZero-1
        dfdv(iv,ips,1)=-(fe(iv-1,ips,iITR)-fe(iv,ips,iITR))/(vpar(iv+1)-vpar(iv))
     enddo
     ! compute unsmoothed derivative
     dfdv(ivZero, ips,1) = 0.0_fp
     dfdv(nv,ips,1) = dfdv(nv-1,ips,1)
     dfdv(1,ips,1) = dfdv(2,ips,1)
  enddo         

end subroutine lscsq_fePru

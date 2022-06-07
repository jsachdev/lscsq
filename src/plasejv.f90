SUBROUTINE lscsq_ProfInit

  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : npsi, psiary, neary, teary, zbrary, lnlary, betzary
  implicit none

  integer :: ip

  do ip = 1, npsi
     CALL lscsq_GePlPar(PsiAry(ip), NeAry(ip), TeAry(ip), ZbrAry(ip), LnlAry(ip), BetZAry(ip)) ! get plasma parameters at given psi
  enddo

end subroutine lscsq_ProfInit


SUBROUTINE lscsq_GePlPar(psi, Negl, Tegl, Zbrgl, Lnlgl, BetZgl)
  ! get plasma parameters
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : pe2fac, pi2fac
  implicit none

  real(fp), intent(in) :: psi
  real(fp), intent(out):: Negl      ! electron density in cm^-3 
  real(fp), intent(out):: Tegl      ! electron temperature in Kev
  real(fp), intent(out):: Zbrgl     ! Z bar for e-i collison rate
  real(fp), intent(out):: Lnlgl     ! Coulomb logarithm
  real(fp), intent(out):: BetZgl    ! Beta value of Eder & Valeo

  real(fp) :: RBphi, pe2, pi2, aio, ael
  real(fp) :: psiold = -100.0_fp

  call lscsq_plasma1d(psi, psiold, RBphi, Tegl, pe2, pi2, aio, ael)
  
  pe2 = max(pe2,1.0e-3_fp*pe2Fac)  ! ensure ne > 1.0e11
  Negl = 1.0e14_fp* pe2/pe2Fac

  Zbrgl = pe2/pi2*(pi2Fac/pe2Fac)
  Zbrgl = max(0.0_fp, min(1.0e2_fp, Zbrgl))  ! ensure 0 < Zbar < 100

  if (Tegl .LE. 0.01_fp) then
     Lnlgl = 23.0_fp-log(sqrt(Negl)*(1.0e-3_fp)**(-1.5_fp)/(Tegl*sqrt(Tegl)))
  else
     Lnlgl = 24.0_fp-log(sqrt(Negl)*1.0e-3_fp/Tegl)
  endif
  BetZgl = 0.2_fp*(1.0_fp+Zbrgl)

end subroutine lscsq_gePlPar

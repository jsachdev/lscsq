subroutine lscsq_grapLSC(xx, zz, psval, psderiv, isw, deex, deez, isave, jsave, fmat)
!***************************************************************
!     *
!     GRid APproximation in the Lower hybrid Simulation Code
!     local bivariate interpolation*
!     *
!***************************************************************
!     isw=1 forces recalculation of coefficient matrix
!     =0 may not recalculate coefficient matrix if same
!     i,j as last call
!
  use iso_c_binding, only : fp => c_double
  use lscsq_mod, only : nx, nz, rmin, zmin, psigrd
 
  implicit none

      INTEGER isw, jval, ival, m, k, l, ii, i, j,        &
           iplus, iminus, jplus, jminus, isave, jsave, frstcall
      INTEGER LLC, LRC, ULC, URC, LEFT,                                 &
           RIGHT, TOP, BOTTOM, INTERIOR 
      PARAMETER(                                                        &
           INTERIOR = 1,                          &
           LLC = 5,  LRC = 6,    ULC = 7,   URC = 8,                    &
           LEFT = 9, RIGHT = 10, TOP = 11,  BOTTOM = 12)
      REAL(fp)                                                              &
           zz, xx, gradsq, dpsidx, dpsidz, psval, psixz,                &
           psixx, psizz, deez, deex,                         &
           rival, rjval, cx, cy, psderiv

  real(fp), dimension(4,4) :: vmat

      REAL(fp)                                                              &
           fmat, poly, xpa, zpa, rhs1, rhs2, rhs3, rhs4, sum1,    &
           sum2, sum3, sum4, sum5, sum6, xp, zp, ff, f1
      DATA frstcall / 1 /
      DIMENSION fmat(0:3,0:*),poly(0:3,0:3)
      DIMENSION xpa(0:3),zpa(0:3)
      DIMENSION psderiv(0:2, 0:*)
      DATA xpa(0),zpa(0)/ 1.d0, 1.d0/

      rjval = (zz-zmin)/deez+1.0_fp
      rival = (xx-rmin )/deex+1.0_fp
      jval = rjval
      ival = rival
 
      if(ival .LT. 2)then
         ival = 2
      else if(ival.GT. nx - 2)then
         ival = nx - 2
      endif
      if(jval .LT. 2)then
         jval = 2
      else if(jval.GT. nz - 2)then
         jval = nz - 2
      endif
  !    if((frstcall.EQ.0).and.(ival.EQ.isave .and. jval.EQ.jsave .and. isw.EQ.0)) go to 200
      frstcall = 0
      isave = ival
      jsave = jval
 
!     begin INTERIOR point
  do ii=1,4
         if (ii.eq.1) then  
!       go to(101,102,103,104),ii
! 101     i      = ival
            i      = ival
            j      = jval
            iplus = i + 1
            jplus = j + 1
            iminus = i - 1
            jminus = j - 1
            cx = 0.5_fp
            cy = 0.5_fp
         elseif (ii.eq.2) then 
!         go to 105
! 102     i      = ival
         i      = ival
         j      = jval+1
         iplus = i + 1
         jplus = j + 1
         iminus = i - 1
         jminus = j - 1
         cx = 0.5_fp
         cy = 0.5_fp
!         go to 105
         elseif (ii.eq.3) then 
! 103     i      = ival+1
         i      = ival+1
         j      = jval
         iplus = i + 1
         jplus = j + 1
         iminus = i - 1
         jminus = j - 1
         cx = 0.5_fp
         cy = 0.5_fp
!         go to 105
! 104     i      = ival+1
      elseif (ii.eq.4) then 
         i      = ival+1
         j      = jval+1
         iplus = i + 1
         jplus = j + 1
         iminus = i - 1
         jminus = j - 1
         cx = 0.5_fp
         cy = 0.5_fp
      endif
! 105     continue
!
!...psi:
!
         vmat(1, ii) =      psigrd(i, j)     - psigrd(ival, jval)
         vmat(2, ii) = cx * (psigrd(iplus, j  ) - psigrd(iminus, j) )
         vmat(3, ii) = cy * (psigrd(i, jplus) - psigrd(i, jminus)  )
         vmat(4, ii) = cx * cy * (psigrd(iplus, jplus) -                   &
     &        psigrd(iminus, jplus) - psigrd(iplus, jminus) +                 &
     &        psigrd(iminus, jminus) )
!
  enddo
! 10   continue
!     end INTERIOR point
!
!
!..........................................................
!
!     calculate function values and function derivatives
!     at four corners of the reference cell
!
!..........................................................
!..........................................................c
!     determine coefficients by evaluating polynomial
!     at four corner points c
!..........................................................c
!...point(i,j):c
      fmat(0,0) = 0.0  
      fmat(1,0) = vmat(2,1)
      fmat(0,1) = vmat(3,1)
      fmat(1,1) = vmat(4,1)
!
!...point(i,j+1):
!
      fmat(0,2) = 3.*vmat(1,2) - 2.*fmat(0,1)                           &
     &     -    vmat(3,2)
      fmat(0,3) =    vmat(3,2) +    fmat(0,1)                           &
     &     - 2.*vmat(1,2)
      fmat(1,2) = 3.*vmat(2,2) - 3.*fmat(1,0)                           &
     &     -    vmat(4,2) - 2.*fmat(1,1)
      fmat(1,3) =    vmat(4,2) + 2.*fmat(1,0)                           &
     &     - 2.*vmat(2,2) +    fmat(1,1)
!
!...point(i+1,j):
!
      fmat(2,0) = 3.*vmat(1,3) - 2.*fmat(1,0)                           &
     &     -    vmat(2,3)
      fmat(3,0) =    vmat(2,3) +    fmat(1,0)                           &
     &     - 2.*vmat(1,3)
      fmat(2,1) = 3.*vmat(3,3) - 3.*fmat(0,1)                           &
     &     -    vmat(4,3) - 2.*fmat(1,1)
      fmat(3,1) =    vmat(4,3) + 2.*fmat(0,1)                           &
     &     - 2.*vmat(3,3) +    fmat(1,1)
!
!...point(i+1,j+1):
!
      rhs1=   vmat(1,4)-vmat(1,3)-vmat(3,3)                             &
     &     -  (fmat(0,2)+fmat(1,2)                                      &
     &     +   fmat(0,3)+fmat(1,3) )
      rhs2=   vmat(2,4)-vmat(2,3)-vmat(4,3)                             &
     &     -  (fmat(1,2)+fmat(1,3) )
      rhs3=   vmat(3,4)-vmat(3,3)                                       &
     &     -2*(fmat(0,2)+fmat(1,2) )                                    &
     &     -3*(fmat(0,3)+fmat(1,3) )
      rhs4=   vmat(4,4)-vmat(4,3)                                       &
     &     -2*(fmat(1,2) )                                              &
     &     -3*(fmat(1,3) )
!
      fmat(2,2) =  9*rhs1-3*rhs2-3*rhs3+rhs4
      fmat(3,2) = -6*rhs1+3*rhs2+2*rhs3-rhs4
      fmat(2,3) = -6*rhs1+2*rhs2+3*rhs3-rhs4
      fmat(3,3) =  4*rhs1-2*rhs2-2*rhs3+rhs4
 200  continue
!...........................................................
!
!     evaluate function and derivatives at (xx,zz)
!
!...........................................................
      xp     = rival - ival
      zp     = rjval - jval
 
!     begin evaluate function
      do m=1,3
         xpa(m) = xp*xpa(m-1)
         zpa(m) = zp*zpa(m-1)
      enddo
      do l=0,3
         do k=0,3
            poly(k,l) = xpa(k)*zpa(l)
         enddo
      enddo
 
      sum1   = 0.0  
      do k = 0,3
         do l = 0,3
            f1     = fmat(k,l)
            sum1   = sum1 + f1*poly(k,l)
         enddo
      enddo
      psval  = sum1 + psigrd(ival,jval)
!     end evaluate function
 
!     begin evaluate derivatives
      sum2   = 0.0 
      do l = 0,3
         do k = 0,2
            ff   = (k+1)*fmat(k+1,l)
            sum2 = sum2 + ff*poly(k,l)
         enddo
      enddo
      dpsidx = sum2/deex
 
      sum3   = 0.0 
      do k = 0,3
         do l = 0,2
            ff   = (l+1)*fmat(k,l+1)
            sum3 = sum3 + ff*poly(k,l)
         enddo
      enddo
      dpsidz = sum3/deez
      gradsq = dpsidx**2 + dpsidz**2
!
      sum5   = 0.0 
      do l = 0,3
         do k = 0,1
            ff   = (k+1)*(k+2)*fmat(k+2,l)
            sum5 = sum5 + ff*poly(k,l)
         enddo
      enddo
      psixx  = sum5/deex/deex
!
      sum6   = 0.0 
      do k = 0,3
         do l = 0,1
            ff   = (l+1)*(l+2)*fmat(k,l+2)
            sum6 = sum6 + ff*poly(k,l)
         enddo
      enddo
      psizz  = sum6/deez/deez
 
      sum4   = 0.0  
      do k = 0,2
         do l = 0,2
            ff   = (k+1)*(l+1)*fmat(k+1,l+1)
            sum4 = sum4 + ff*poly(k,l)
         enddo
      enddo
      psixz  = sum4/deex/deez
!     end evaluate derivatives
!
      psderiv(0, 0) = gradsq
      psderiv(1, 0) = dpsidx
      psderiv(2, 0) = psixx
      psderiv(0, 1) = dpsidz
      psderiv(0, 2) = psizz
      psderiv(1, 1) = psixz

end subroutine lscsq_grapLSC
!
!     ------------------------------------------------------------------
!
SUBROUTINE lscsq_grnu1d(nx, xdata, ydata, jold, coef, x, y, yp )

  use iso_c_binding, only : fp => c_double
  implicit none

  integer, intent(in) :: nx
  integer, intent(inout) :: jold
  real(fp), dimension(nx), intent(in) :: xdata
  real(fp), dimension(nx), intent(in) :: ydata
  real(fp), intent(in) :: x
  real(fp), intent(out) :: y 
  real(fp), dimension(4), intent(out) :: coef
  real(fp), intent(out) :: yp 

  integer :: j

  real(fp) :: x21,      x32
  real(fp) :: d21, d31, d32, d42, dx

!     GRid, NonUniform, interpolation in 1D
!     cubic interpolation yint(x) for points 1, 2, 3, 4 set by conditions
!     yint(2) = ydata(2)
!     yint(3) = ydata(3)
!     d yint(2) / dx = (ydata(3) - ydata(1)) / (xdata(3) - xdata(1) )
!     d yint(3) / dx = (ydata(4) - ydata(2)) / (xdata(4) - xdata(2) )
!
!     here: yint = sum(1 - 4) * coef(j) * ((x - x(2)) ** (j - 1))

  j = jold
  if (xdata(j).gt.x) j=j-1

  coef(1:4) = 0.0 

  if (j .EQ. 0 ) then
     y  = ydata(1)
     yp = 0.0 
     return
  endif
  if ( j.EQ. nx) then
     y  = ydata(nx)
     yp = 0.0 
     return
  endif

  !     compute new coefficients
  jold = j
  if((j.GT.1) .and. (j.LT.(nx-1))) then
     x32     = (xdata(j+1) - xdata(j))
     d31     = (ydata(j+1) - ydata(j-1)) / (xdata(j+1) - xdata(j-1))
     d32     = (ydata(j+1) - ydata(j)) / (xdata(j+1) - xdata(j))
     d42     = (ydata(j+2) - ydata(j)) / (xdata(j+2) - xdata(j))

     coef(1) = ydata(j)
     coef(2) = d31
     coef(3) = (3.0_fp*(d32-d31) -    (d42-d31)) / x32
     coef(4) = (   (d42-d31) - 2.0_fp*(d32-d31)) / (x32*x32)
  else if(j.EQ.1)then
     x21     = (xdata(2) - xdata(1))
     d21     = (ydata(2) - ydata(1)) /(xdata(2) - xdata(1))
     d31     = (ydata(3) - ydata(1)) /(xdata(3) - xdata(1))
     coef(1) = ydata(1)
     coef(2) = (2.0_fp*d21 - d31)
     coef(3) = (d31 - d21) / x21
     coef(4) = 0.0 
  else if(j.EQ.(nx-1))then
     x32     = (xdata(nx) - xdata(nx-1))
     d31     = (ydata(nx) - ydata(nx-2)) / (xdata(nx) - xdata(nx-2))
     d32     = (ydata(nx) - ydata(nx-1)) / (xdata(nx) - xdata(nx-1))
     coef(1) = ydata(nx-1)
     coef(2) =  d31
     coef(3) = (d32-d31) / x32
     coef(4) = 0.0 
  endif
  !     compute y, yp
  dx = x - xdata(j)
  y  = ((coef(4)*dx + coef(3))*dx + coef(2))*dx + coef(1)
  yp = ((3.0_fp*coef(4)*dx + 2.0_fp*coef(3))*dx + coef(2))

end subroutine lscsq_grnu1d 
 
!     ------------------------------------------------------------------
subroutine lscsq_linr1d(nx, xdata, ydata, x, y )
  
  use iso_c_binding, only : fp => c_double
  implicit none

  integer, intent(in) :: nx

  integer :: j, jguess

  real(fp), dimension(nx) :: xdata, ydata
  real(fp), intent(out) :: x, y
  real(fp) :: coef(2)

  real(fp) :: d21, dx
!     LINeaR interpolation in 1D
!
!     Try to look much like grnu1d, which see above.
!     The derivative is guaranteed not continuous.
!     Given to make consistent with grnu1d.
!     Linear interpolation of interest to avoid fake oscillations
!     in the interpolant.
!
!     here: yint = sum(1 - 2) * coef(j) * ((x - x(1)) ** (j - 1))
!
!                                       find which interval we are in by
!                                       Numerical Recipes HUNT routine
  j = minloc(abs(xdata-x),1)
  if (xdata(j).gt.x) j=j-1

  if (j .EQ. 0 ) then
     y  = ydata(1)
     return
  endif
  if ( j.EQ. nx) then
     y  = ydata(nx)
     return
  endif
 
  ! This is the computation of the slope, and delta-x
  d21 = (ydata(j+1)-ydata(j)) / (xdata(j+1)-xdata(j))

  coef(1) = ydata(j)
  coef(2) = d21
  dx = x - xdata(j)
 
  y  = coef(2)*dx + coef(1)

end subroutine lscsq_linr1d 
 
!     ------------------------------------------------------------------
!

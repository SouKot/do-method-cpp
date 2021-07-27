!* ============================================================================ *
!*
!* Multi-Layer Shallow Water Model
!*
!* Version 1.2 - 20/02/2008 - ATvS
!*
!* ============================================================================ *
! NOTES BY KURT LUST:
! - I'm still using fixed-format sources, but I've been carefull to make sure
!   that conversion to free-format should be easy.
! - These routines interface with fitpack and are more or less equivalent
!   to the NAG routines used in older versions of the code.
!
      subroutine spfit( psi, npsi, theta, ntheta, field, knotpsi, knottheta, c, ifail )
!
! subroutine spfit
!
! This routine computes an interpolating spline function through
! the data. We assume that the data is on a sphere and that no values
! are given at the poles.
!
! The spline returned by this routine is defined for psi in [0, 2 pi]
! and theta in [-pi/2, pi/2]. The knots are periodic in the 
! psi-direction. The computed spline s(psi,theta) satisfies the following 
! boundary conditions:
!
!  (1) s(psi,+/-pi/2) = s(0,+/-pi/2)
!
!      d s(psi,+/-pi/2)
!  (2) ---------------- 
!          d theta
!                    d s(0,+/-pi/2)              d s(pi/2,+/-pi/2)
!       = cos(psi) * -------------- + sin(theta) -----------------
!                       d theta                      d theta
!
! for all psi in [0, 2 pi].
! 
! Input variables:
! - psi(npsi): The psi-coordinates of the data points.
! - theta(ntheta): The theta-coordinates of the data points.
! - field(npsi,ntheta): The field to interpolate.
!
! Output variables:
! - knotpsi(npsi+7): Vector with the knots in the psi direction. 
! - knottheta(ntheta+8): Vector with the knots in the theta-direction.
! - c(npsi+3,ntheta+4): Matrix with the B-spline coefficients computed 
!   with fitpack.
! - ifail: Error code. 0 when successful.
!
      implicit none
      include 'constants.com'
!
      integer, intent(in):: npsi, ntheta
      real, intent(in):: psi(npsi), theta(ntheta), field(npsi,ntheta)
      real, intent(out):: knotpsi(npsi+7), knottheta(ntheta+8)
      real, intent(out):: c(npsi+3,ntheta+4)
      integer, intent(out):: ifail
!
      integer:: i, j
      integer:: iopt(3)
      integer:: ider(4)
      real:: thetashift(ntheta)
      integer:: nknpsi, nkntheta
      integer:: nknpsiest, nknthetaest
      real r0, r1, s, fp
! Size of work space for fitpack:
!   nknotpsi = npsi + 7
!   nknottheta = ntheta + 8
!   lwrk = 12 + nknottheta * (npsi + nknotpsi + 3)
!          + 24 * nknotpsi + 4 * ntheta + 8 * npsi
!          + npsi + nknotpsi + nknottheta
!        = 50 * npsi + 15 * ntheta + 2 * npsi * ntheta + 275
!   kwrk = 5 + npsi + ntheta + nknotpsi + nknottheta
!        = 2*npsi + 2*ntheta + 20
      integer:: lwrk, kwrk
      integer:: iwrk(2*npsi + 2*ntheta + 20)
      real:: wrk(50*npsi + 15*ntheta + 2*npsi*ntheta + 275)
!
! Since this piece of code is by no means time-critical, it makes
! perfect sense to test whether a number of assumptions that have been
! made when writing the code, are met. This prevents nasty bugs.
!
! - Are all psi-values strictly increasing and in [0, 2 pi)?
      if( ( psi(1) < 0. ) .or. ( psi(npsi) >= tpi ) ) then
        write(stderr,*) 'Some psi-values are out-of-range'
        ifail = 1
        return
      end if
      do i = 1, npsi - 1
        if ( psi(i+1) <= psi(i) ) then
          write(stderr,*)&
      'The strict ordering of the psi-values is not satisfied.'
          ifail = 1
          return
        end if
      end do
! - Are all theta-values strictly ordered and in (-pi/2, pi/2)?
      if( ( theta(1) <= -pio2 ) .or. ( theta(ntheta) >= pio2 ) ) then
        write(stderr,*) 'Some theta-values are out-of-range'
        ifail = 1
        return
      end if
      do j = 1, ntheta - 1
        if ( theta(j+1) <= theta(j) ) then
          write(stderr,*)& 
           'The strict ordering of the theta-values is not satisfied.'
          ifail = 1
          return
        end if
      end do
!
! Now we'll shift the theta values to bring to adapt them to the coordinate
! system used by spgrid.
!
      thetashift = theta + pio2
!
! Next we construct the knot vector in the psi-direction. We use the 
! data values for this purpose.
!
      knotpsi(1:3) = psi(npsi-2:npsi) - tpi
      knotpsi(4:npsi+3) = psi
      knotpsi(npsi+4:npsi+7) = psi(1:4) + tpi
!
! Next we construct the knot vector in the theta direction. Here we
! have to take into account that fitpack uses a different convention 
! than the QG code.
!
      knottheta(1:4) = 0.
      knottheta(5:ntheta+4) = thetashift
      knottheta(ntheta+5:ntheta+8) = 0.    
!
! Other preparations for fitpack: Set fitpack input variables.
!
      iopt = (/ -1, 1, 1 /)
      ider = (/ -1, 0, -1, 0 /)
      nknpsi = npsi + 7
      nknpsiest = nknpsi
      nkntheta = ntheta + 8
      nknthetaest = nkntheta
      r0 = 0.
      r1 = 0.
      s = 0.
      lwrk = 50*npsi + 15*ntheta + 2*npsi*ntheta + 275
      kwrk = 2*npsi + 2*ntheta + 20
!
! Call the fitpack spgrid routine.
! Note that we first put iwrk to 0 to work around a bug in fitpack.
!
      iwrk = 0
      call spgrid( iopt, ider, ntheta, thetashift, npsi, psi, field,&
                  r0, r1, s, nknthetaest, nknpsiest,&
                  nkntheta, knottheta, nknpsi, knotpsi, c, fp,&
                  wrk, lwrk, iwrk, kwrk, ifail )
!
! Returned from fitpack. 
! - Check for errors.
      if ( ifail > 0 ) then
        write(stderr,*) 'Fitpack failed to compute the spline.'
        write(stderr,*) 'Error return code: ', ifail
        return
      end if
! - Convert the theta-knots back to the coordinate system favored by 
!   the ocean model.
      knottheta = knottheta - pi / 2.
!
      end
! end spfit
!
!

      subroutine speval( ncoefpsi, ncoeftheta, knotpsi, knottheta, c,&
                        gridpsi, ngridpsi, gridtheta, ngridtheta,&
                        field, ifail )
!
! subroutine speval
!
! This routine evaluates a bicubic spline function on a grid of points.
!
! Input arguments:
! - ncoefpsi, ncoeftheta: The number of spline coefficients in the psi-
!   and theta-directions.
! - knotpsi(ncoefpsi+4): Knots in the psi-direction.
! - knottheta(ncoeftheta+4): Knots in the theta-direction.
! - c(ncoefpsi,ncoeftheta): The coefficients of the bicubic spline in the
!   B-spline basis.
! - gridpsi(ngridpsi): psi-values at which the spline should be evaluated.
! - gridtheta(ngridtheta): theta-values at which the spline should be 
!   evaluated.
! 
! Output arguments:
! - field(ngridpsi,nthetapsi): Computed approximation.
! - ifail: Error code. The routine was successfull if ifail=0.
!
! Assumptions made:
! - The values in gridtheta should be strictly increasing and in [-pi/2, pi/2].
! - This implementation makes no assumptions at all about the psi variables.
!   Since we want to allow values outside of [0, 2 pi] to allow for a 
!   computational domain crossing the 0-meridian, we cannot exploit the 
!   ordering of the psi variables much.
! - The spline is defined for psi in [0, 2 pi] and theta in [-pi/2, pi/2].
!
      implicit none
      include 'constants.com'
!
      integer, intent(in):: ncoefpsi, ncoeftheta
      real, intent(in):: knotpsi(ncoefpsi+4)
      real, intent(in):: knottheta(ncoeftheta+4)
      real, intent(in):: c(ncoefpsi,ncoeftheta)
      integer, intent(in):: ngridpsi, ngridtheta
      real, intent(in):: gridpsi(ngridpsi), gridtheta(ngridtheta)
      real, intent(out):: field(ngridpsi,ngridtheta)
      integer, intent(out):: ifail
!
      integer, parameter:: degree = 3
!
      integer i, j
      real psib
      real dfield(ngridtheta)
      integer lwrk, kwrk
      real wrk((degree+1)*(ngridtheta+1))
      integer iwrk(ngridtheta+1)
!
! First we'll check the assumption about the theta-coordinates of the
! grid points: are they in [-pi/2, pi/2] and strictly ordered?
!
      if( ( gridtheta(1) < -pio2 ) .or.& 
         ( gridtheta(ngridtheta) > pio2 ) ) then
        write(stderr,*)& 
         'Some theta-values of the grid points are out-of-range'
        ifail = 1
        return
      end if
      do j = 1, ngridtheta - 1
        if ( gridtheta(j+1) <= gridtheta(j) ) then
          write(stderr,*)& 
           'The strict ordering of the theta-values (grid) is not OK.'
          ifail = 1
          return
        end if
      end do
!
! Set the other input variables for bispev.
!
      lwrk = (degree+1) * (ngridtheta + 1)
      kwrk = ngridtheta + 1
!
! Loop over all psi-values of the grid:
!
      do i = 1, ngridpsi
! * Bring into [0, 2 pi].
        psib = gridpsi(i)
        do while ( psib < 0. )
          psib = psib + tpi
        end do
        do while ( psib > tpi )
          psib = psib - tpi
        end do
! * Evaluate the spline for this value of psi and all values of gridtheta.
        call bispev( knottheta, ncoeftheta+4, knotpsi, ncoefpsi+4,&
                    c, degree, degree,& 
                    gridtheta, ngridtheta, psib, 1, dfield,&
                    wrk, lwrk, iwrk, kwrk, ifail )
! * Test the return code.
        if ( ifail /= 0 ) then
          write(stderr,*)& 
           'Evaluation of the spline failed.',&
           'There might be a problem with the arguments of bispev.'
          return
        end if
! * Copy to the matrix receiving the result.
        field(i,:) = dfield
      end do
!
      end
! end speval
!
!

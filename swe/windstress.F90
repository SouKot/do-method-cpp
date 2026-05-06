!* ============================================================================ *
!*
!* Multi-Layer Shallow Water Model
!*
!* Version 1.2-20/02/2008-ATvS
!*
!* ============================================================================ *
! This file contains all routines that have to do with the initialisation
! of the wind stress field.
!
! NOTES BY KURT LUST:
! - I'm still using fixed-format sources, but I've been carefull to make sure
!   that conversion to free-format should be easy.
! - I've parameterised windfit since it needs so little from the common
!   blocks that it is preferable to pass the variables as arguments.
!

!*******************************************************************************
!
! windfun
!
!*******************************************************************************

      real function windfun( xx, yy, v1, curl )
!
!  function windfun
!
!  This function gives the wind stress field.
!
!     ! implicit none
!      use m_par
!      use m_usr
!      implicit none
!! Arguments
!      real, intent(in)    :: xx, yy
!      integer, intent(in):: v1
!      real, intent(in)    :: curl
!! Local variables
!      real:: pi
!      real:: alp, bet, falp, rl, ys, ym
!
!      pi = 4.0*atan(1.0)
!!     alp = 9.0
!!     bet = pi/4.*(1.0-curl)
!!     falp = 1./(alp*alp-1)
!!     rl = r0dim/1.e+06
!!     ys = -39.5*pi/180.
!
!!     windfun = -cos(2*pi*(yy-ymin)/(ymax-ymin))
!      !ym = tghmax*pi/180.0
!      ys = -(40.0+par(16))*pi/180.
!      if (yy .ge. ys) then
!         windfun = cos(pi*(yy-ys)/(ygmax-ys))-1.
!      !  windfun = 0.03
!      else
!         windfun = 0.0
!      endif


        use m_par
        use m_usr
        implicit none
! Arguments
        real, intent(in)    :: xx, yy
        integer, intent(in):: v1
        real, intent(in)    :: curl
! Local variables
        real:: pi, k1, k2, earth_radius, phi_b, dd, bx, aphi
        real:: alp, bet, falp, rl, ys, ym, yyd

        pi = 4.0*atan(1.0)

        k1 = 2.0*pi/38.0
        k2 = 2.0*pi/70.0
        phi_b = 37.0-pi/k1
        dd = r0dim*cos(yy)*(xx-xgmin)
        bx = sin(pi*dd/(2*10700000) + pi/4)
        yyd = yy*180/pi
        if (yyd .lt. 5.0) then
          aphi = -cos(k2*(5.0-phi_b))
        elseif (yyd .ge. 5.0 .and. yyd .lt. phi_b ) then
          aphi = -cos(k2*(yyd-phi_b))
        elseif (yyd .ge. phi_b .and. yyd .lt. 37.0) then
          aphi = cos(k1*(yyd-37.0))
        elseif (yyd .ge. 37.0) then
          aphi = cos(k2*(yyd-37.0))
        endif
        windfun = 0.046*aphi*bx-0.004
        !windfun = 0.0

        end


!*******************************************************************************
!
! windforcing
!
!*******************************************************************************

        subroutine windforcing( windfromdata, gridpsi, ngridpsi, gridtheta, ngridtheta, taux, tauy, tmax )
!
! subroutine windforcing
!
! Initialises the shape of the wind-forcing
!
          implicit none
! Arguments.
          logical, intent(in):: windfromdata
          integer, intent(in):: ngridpsi, ngridtheta
          real, intent(in):: gridpsi(ngridpsi), gridtheta(ngridtheta)
          real, intent(out):: taux(ngridpsi, ngridtheta)
          real, intent(out):: tauy(ngridpsi, ngridtheta)
          real, intent(out):: tmax
! Local variables.
          real:: curl = 0.0

          if ( windfromdata ) then
            call windfit( gridpsi, ngridpsi, gridtheta, ngridtheta, taux, tauy, tmax )    
          else
! TODO: We only call windfit to get the same initialisation of tmax
! as in the original code, but this initialisation is probably wrong!
!       call windfit( gridpsi, ngridpsi, gridtheta, ngridtheta, 
!    &                taux, tauy, tmax )    
            curl = 0.0
            call windfromfun( gridpsi, ngridpsi, gridtheta, ngridtheta, curl, taux, tauy, tmax )    
          endif

          end

!*******************************************************************************
!
! windfromfun
!
!*******************************************************************************

          subroutine windfromfun( gridpsi, ngridpsi, gridtheta, ngridtheta, curl, taux, tauy, tmax )
!
! subroutine windfromfun.
!
! This routine computes the wind stress field from the given user
! function windfun (declared at the top of this file).
!
! Input parameters:
! - gridpsi(ngridpsi) and gridtheta(ngridtheta): The psi-and theta
!   components of the grid points. 
! - curl: A parameter which is passed on to the function windfun.
!
! Output parameters:
! - taux(ngridpsi, ngridtheta) and tauy(ngridpsi, ngridtheta): 
!   The two components of the windstress.
! - tmax: Scaling factor used at the end to rescale taux and tauy.
!
! Note that I have chosen to add a parameter list to windfromfun so that
! I could get rid of the common blocks. This reduces the chance for
! errors and makes this routine more easily portable to other codes.
!
! Kurt Lust, summer 2006.
!
            implicit none
! Arguments.
            integer, intent(in):: ngridpsi, ngridtheta
            real, intent(in):: gridpsi(ngridpsi), gridtheta(ngridtheta)
            real, intent(in):: curl
            real, intent(out):: taux(ngridpsi, ngridtheta)
            real, intent(out):: tauy(ngridpsi, ngridtheta)
            real, intent(out):: tmax
! Local variables
            integer:: i, j
! Functions used.
            real    :: windfun

            do j = 1, ngridtheta
            do i = 1, ngridpsi
            taux(i, j) = windfun( gridpsi(i), gridtheta(j), 1, curl )
            tauy(i, j) = 0.
            enddo
            enddo

! Note that the following line is only valid for the current windfun.
! With a more complicated function, the lines from the end of windfit
! should be used.

! TODO: The next line seems to be a more approriate initilisation of tmax!
            tmax = 1. ! ATvS-01/09

            end

!*******************************************************************************
!
! windfit
!
!*******************************************************************************

            subroutine windfit( gridpsi, ngridpsi, gridtheta, ngridtheta, taux, tauy, tmax )
!
! subroutine windfit.
!
! This routine reads the windstress data and interpolates the data to
! the regular computational grid.
!
! Input parameters:
! - gridpsi(ngridpsi) and gridtheta(ngridtheta): The psi-and theta
!   components of the grid points. 
!
! Output parameters:
! - taux(ngridpsi, ngridtheta) and tauy(ngridpsi, ngridtheta): 
!   The two components of the windstress.
! - tmax: Scaling factor used at the end to rescale taux and tauy.
!
! Note that I have chosen to add a parameter list to windfit so that I
! could get rid of the common blocks. This reduces the chance for errors
! and makes this routine more easily portable to other codes.
!
! Kurt Lust, summer 2006.
!
              implicit none
              include 'constants.com'
! Arguments.
              integer, intent(in):: ngridpsi, ngridtheta
              real, intent(in):: gridpsi(ngridpsi), gridtheta(ngridtheta)
              real, intent(out):: taux(ngridpsi, ngridtheta)
              real, intent(out):: tauy(ngridpsi, ngridtheta)
              real, intent(out):: tmax
! Constants
              integer, parameter:: uin = 10
! Local variables
              integer npsi, ntheta
              integer i, j
              integer ifail
              real, allocatable:: psi(:), theta(:)
              real, allocatable:: fieldu(:,:), fieldv(:,:)
              real, allocatable:: knotpsi(:), knottheta(:)
              real, allocatable:: c(:,:)
!
! Read trtau.dat.
!
              open(uin, file='trtau.dat')
! - Read the size of the data grid.
!   This is not yet robust  ! I wonder why I have to omit the "=" after m.
10    format( ' trenberth: tx, ty; n=',I0, ',m',I0)
              read(uin, 10) npsi, ntheta
              print *, 'npsi = ', npsi, ', nntheta = ', ntheta
              npsi = npsi-1  ! We don't want to read the excess psi value.
! - Allocate memory.
              allocate( psi(npsi), theta(ntheta), fieldu(npsi, ntheta), fieldv(npsi, ntheta), stat = ifail )
              psi = 0; theta = 0; fieldu = 0; fieldv = 0
              if ( ifail /= 0 ) then
                write(stderr, *) "windfit: Failed to allocate memory."
                stop
              end if
! - Read the psi-and theta-values.
              do i = 1, npsi
              read(uin, *) psi(i)
              end do
              read(uin, *)  ! Skip the last psi-value.
              do j = 1, ntheta
              read(uin, *) theta(j)
              end do
! - Read the windstresses. We'll skip the last line.
              do i = 1, npsi
              do j = 1, ntheta
              read(uin, *) fieldu(i, j), fieldv(i, j)
              end do
              end do
              close(uin)
!
! Convert psi and theta to radians.
!
              psi   = psi*pi/180.
              theta = theta*pi/180.
!
! Allocate memory for the spline knots and coefficients.
!
              allocate( knotpsi(npsi+7), knottheta(ntheta+8), c(npsi+3, ntheta+4), stat = ifail )
              knotpsi = 0; knottheta = 0; c = 0
              if ( ifail /= 0 ) then
                write(stderr, *) "windfit: Failed to allocate memory."
                stop
              end if
!
! Now, for each field, compute the spline approximation and interpolate
! on the grid.
!
              call spfit( psi, npsi, theta, ntheta, fieldu, knotpsi, knottheta, c, ifail )
              if ( ifail /= 0 ) then
                write(stderr, *)& 
                  "windfit: Failed to compute the spline approximation."
                stop
              end if
              call speval( npsi+3, ntheta+4, knotpsi, knottheta, c, gridpsi, ngridpsi, gridtheta, ngridtheta, taux, ifail )
              if ( ifail /= 0 ) then
                write(stderr, *)& 
                  "windfit: Failed to evaluate the spline approximation."
                stop
              end if
              call spfit( psi, npsi, theta, ntheta, fieldv, knotpsi, knottheta, c, ifail )
              if ( ifail /= 0 ) then
                write(stderr, *)& 
                  "windfit: Failed to compute the spline approximation."
                stop
              end if
              call speval( npsi+3, ntheta+4, knotpsi, knottheta, c, gridpsi, ngridpsi, gridtheta, ngridtheta, tauy, ifail )
              if ( ifail /= 0 ) then
                write(stderr, *)& 
                  "windfit: Failed to evaluate the spline approximation."
                stop
              end if
!
! Rescale the stresses.
!
              tmax = max( maxval(taux), maxval(tauy) )
              taux = taux/tmax
              tauy = tauy/tmax
!
! Clean-up the workspace.
!
              deallocate( psi, theta, fieldu, fieldv )
              deallocate( knotpsi, knottheta, c )
!     
              write(uout, *) 'fit of windstress done, tmax  = ', tmax
!
              end


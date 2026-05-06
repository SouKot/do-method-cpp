! module implementing the discretization for the 3D lid driven cavity problem
module m_LidDrivenCavity

#define _SUB2IND_(nn, mm, ll, dd, ii, jj, kk, vv) ((((kk)-1)*(mm) + (jj)-1)*(nn) + (ii)-1)*(dd) + (vv)

use, intrinsic:: iso_c_binding
use m_types, m_types_stat=>istat
use m_tools
!use m_constants
!use m_spf

implicit none

TYPE(Model_t), save:: Model
!TYPE(Matrix_t), save:: jac

! we define some quantities that essentially point
! to data structures in the model to simplify the
! programming
integer, save:: nx, ny, nz, ndim, ndof, npBiL, np2, npar2  ! np and npar is removed!!!!!!!!!
integer, save:: stxb, stxe, styb, stye, stzb, stze

REAL(kind = dp), DIMENSION(:,:,:,:,:,:), POINTER, save:: Al
REAL(kind = dp), DIMENSION(:,:,:,:,:,:), POINTER, save:: BiL
REAL(kind = dp), DIMENSION(:,:,:,:,:,:), POINTER, save:: AnlF
REAL(kind = dp), DIMENSION(:,:,:,:,:,:), POINTER, save:: AnlJ
REAL(kind = dp), DIMENSION(:), POINTER, save:: mass
REAL(kind = dp), DIMENSION(:,:,:,:), POINTER, save:: Frc
REAL(kind = dp), DIMENSION(:), POINTER, save:: x, y, z
REAL(kind = dp), DIMENSION(:,:,:), POINTER, save:: u, v, w, p
REAL(kind = dp), save:: nlin, oldparam

! define ordering of continuation parameters
integer, parameter ::  TIME = 1

contains
!==========================================================================================
subroutine initialize(nx_in, ny_in, nz_in, phi_in, theta_in, z_in, ngx, ngy, ngz,  xmini, xmaxi, ymini, &
      ymaxi, xgmini, xgmaxi, ygmini, ygmaxi, material, nrows, sol, ierr)&
      bind(C, name='model_init')

 use m_par
 use m_mat
 use m_usr
integer(kind = c_int), intent(in):: nx_in, ny_in, nz_in, ngx, ngy, ngz
real(kind = c_double):: xmini, xmaxi, xgmini, xgmaxi, dphi, dtheta, dz
real(kind = c_double):: ymini, ymaxi, ygmini, ygmaxi
real(kind = c_double), dimension(nx_in+1):: phi_in
real(kind = c_double), dimension(ny_in+1):: theta_in
real(kind = c_double), dimension(nz_in+1):: z_in
integer(kind = c_int), dimension((nx+2)*(ny+2)*(nz+2)), intent(in):: material
integer(kind = c_int), intent(out):: ierr
integer(kind = c_int), intent(in):: nrows
real(kind = c_double), dimension(nrows), intent(out):: sol
! local variables
integer:: i, j, k, pos
real(kind = c_double):: zero, one, pi
!ndim = 3
ndof = ndim
npar2 = 1  ! just Reynolds right now
np2 = 0  ! number of possible neighbours in stencil
npBiL = 0  !width of 1D stencils for Bilinear form
nx = nx_in
ny = ny_in
nz = nz_in
pmin = xmini
pmax = xmaxi
thmin = ymini
thmax = ymaxi
pgmin = xgmini
pgmax = xgmaxi
tghmin = ygmini
tghmax = ygmaxi

!print *,' theta_in ', theta_in
!allocate_grid(xmin, xmax, ymin, ymax)
!m = ny
!nl = nz
call allocate_usr(nx, ny, nz)
one = 1.0
pi = 4*atan(one)
!x = phi_in*pi/180
!y = theta_in*pi/180
!do j = 1, ny-1
! yv(j)=(y(j)+y(j+1))/2
!enddo
!print *, 'ngx, ngy ',ngx, ngy
!dy=(ygmaxi-ygmini)*pi/(ngy*180)
!dx=(xgmaxi-xgmini)*pi/(ngx*180)
!Print *, 'value of pmin pmax thmin thmax ',pmin, pmax, thmin, thmax
!print *,'value of n m and nl are : ', n, m, nl
call allocate_mat()
call AllocateModel(Model, ndim, ndof, npar2, nx, ny, nz, np2, npBiL)
!
if (m_types_stat /= 0) then
  ierr = -1
  return
end if
!-----------------------------------------------------------
!The stencils for the nonlinear part immediately use x
!so these are not used currently
!call AllocateField(Model%Field(UU), 0, n, 1, m, 0, l+1)
!call AllocateField(Model%Field(VV), 1, n, 0, m, 0, l+1)
!call AllocateField(Model%Field(WW), 0, n+1, 0, m+1, 0, nz)
!call AllocateField(Model%Field(PP), 0, n+1, 0, m+1, 0, l+1)
!

if (m_types_stat /= 0) then
  ierr = -1
  return
end if


!Al => Model%Stencil%Al
!BiL => Model%Stencil%BiL
!AnlF => Model%Stencil%AnlF
!AnlJ => Model%Stencil%AnlJ
!Frc => Model%Forcing%Frc
!x => Model%Grid%x
!y => Model%Grid%y
!z => Model%Grid%z
!stzb = Model%Stencil%stzb
!styb = Model%Stencil%styb
!stxb = Model%Stencil%stxb
!stze = Model%Stencil%stze
!stye = Model%Stencil%stye
!stxe = Model%Stencil%stxe

!u => Model%Field(UU)%Values
!v => Model%Field(VV)%Values
!w => Model%Field(WW)%Values
!p => Model%Field(PP)%Values

!Model % Header % Label = 'Lid-Driven Cavity'

!dond do
!o j = 0, ny
! Model%Grid%y(j) = y_in(j+1)
!nd do
!o k = 0, nz
! Model%Grid%z(k) = z_in(k+1)
!nd do
!
! we need some extra coordinates at the boundaries
!Model%Grid%x(-1) = 2.d0*Model%Grid%x(0) - Model%Grid%x(1)
!Model%Grid%x(nx+1) = 2.d0*Model%Grid%x(nx) - Model%Grid%x(nx-1)
!
!Model%Grid%y(-1) = 2.d0*Model%Grid%y(0) - Model%Grid%y(1)
!Model%Grid%y(ny+1) = 2.d0*Model%Grid%y(ny) - Model%Grid%y(ny-1)
!
!Model%Grid%z(-1) = 2.d0*Model%Grid%z(0) - Model%Grid%z(1)
!Model%Grid%z(nz+1) = 2.d0*Model%Grid%z(nz) - Model%Grid%z(nz-1)
!
!do k = 0, nz+1
!  do j = 0, ny+1
!    do i = 0, nx+1
!      pos = _SUB2IND_(nx+2, ny+2, nz+2, 1, i+1, j+1, k+1, 1)
!      Model % Grid % Material(i, j, k) = material(pos)
!     end do
!  end do
!end do
!
!!     Adaptation for LidDrivenCavity
!Model % Grid % Material(1:nx, 1:ny, nz+1)=MLID

  !call AddModelParameter(Model, REYN, 'Reynolds Number',1.0D0)
 zero = 0.00
  call AddModelParameter(Model, TIME, 'Time',zero)
!
!! build linear part of matrix coefficients
!call lin
!Frc = 0.D0  ! don't forget this line: forces are added
!call boundaries(Model%Stencil, Model%Forcing,  Model%Grid%Material)
!
!!set oldparam which tracks whether coefficient has changed.
!oldparam = Model%Param(REYN)%value
!PRINT *, oldparam
!!build coefficients from bilinear form
!!with nlin we can indicate whether non linear terms will be added or not
!! nlin .ne. 0 will add these terms
!nlin = 1.0D0
!call bilin
!call BndBiL(Model%Stencil, Model%Grid%Material)
!AnlF = 0.0D0
!AnlJ = 0.0D0
!#ifdef DEBUGGING
!! for debugging purposes...
!open(unit = 42, file='forcing.txt')
!call DumpForcing(Model%Forcing, 42)
!close(42)
!#endif

 call stpnt(sol)
 call massB

 mass => coB

! print *,(sol(i), i = 1, nx*ny*nz*ndof)
 call init


end subroutine initialize
!
subroutine free() bind(C, name='model_free')

implicit none

  call DeleteModel(Model)

end subroutine free
!=================================================================================
subroutine get_num_params(npar) bind(C, name='model_get_num_params')

implicit none
integer(kind = c_int), intent(out):: npar

npar = Model%Header%npar

end subroutine get_num_params

! get parameter name i (i 1-based), null-terminated so it can be used in C)
!============================================================================

subroutine get_param_name(idx, maxlen, name, default_value, ierr) &
 bind(C, name='model_get_param_name')

integer(kind = c_int), intent(in):: idx, maxlen
character(kind = c_char), dimension(maxlen), intent(out):: name
real(kind = c_double), intent(out):: default_value
integer(kind = c_int), intent(out):: ierr

integer:: i, length

default_value = 0.0
do i = 1, maxlen
  name(i)=''
end do
if (idx > Model%Header%npar) then
  ierr = -1
  return
end if
length = min(maxlen-1, LEN_TRIM(Model%Param(idx)%name))
if (length > maxlen-1) then
  ierr = 1
end if
do i = 1, length
  name(i) = Model%Param(idx)%name(i:i)
end do
name(length+1)=C_NULL_CHAR

default_value = Model%Param(idx)%default_value

end subroutine get_param_name


!======================================================================================

subroutine set_params(npar, pvalues, ierr) bind(C, name='model_setparams')

integer(kind = c_int), intent(in):: npar
real(kind = c_double), dimension(npar), intent(in):: pvalues
integer(kind = c_int):: ierr

integer:: i

ierr = 0

if (npar /= Model%Header%npar) then
  ierr = -1
  return
end if

do i = 1, npar
  if (Model % Param(i) % Idx /= i) then
    ierr = -2
    return
  end if
  Model % Param(i) % Value = pvalues(i)
end do
return
end subroutine set_params
!*****************************************************************************
     subroutine boundaries_bilin
!*     insert conditions at the 'real' boundaries of the domain
!*     changes in boundary conditions here must be implemented
!*     in 'solu' as well!
     !implicit none
     use m_par
     use m_usr
     use m_mat
     implicit none
!*     local
     integer i, j, v, ms, l, k
!*
!*     2  5 (8)  
!*     1  4  7
!*        3  6
!*
!*     pos 2 wordt alleen gebruikt in vu
!*     pos 6 wordt alleen gebruikt in uv
!*
     if (pper) then
     ms = mper
     Alvu(1, mper, 2, :,:) = 0.0
   else
     ms = 1
   endif	
     do k = 1, nl
     do j = ms, m
     ! x = 0
     ! -> v(0, i)+v(1, i)=0
     aluv(1, j, 4, k, k)=aluv(1, j, 4, k, k)-aluv(1, j, 1, k, k)
     aluv(1, j, 1, k, :)=0.0
     aluv(1, j, 5, k, k)=aluv(1, j, 5, k, k)-aluv(1, j, 2, k, k)
     aluv(1, j, 2, k, :)=0.0
     alvv(1, j, 4, k, k)=alvv(1, j, 4, k, k)-alvv(1, j, 1, k, k)
     alvv(1, j, 1, k, :)=0.0
     alvv(1, j, 5, k, k)=alvv(1, j, 5, k, k)-alvv(1, j, 2, k, k)
     alvv(1, j, 2, k, :)=0.0
     alhv(1, j, 4, k, k)=alhv(1, j, 4, k, k)-alhv(1, j, 1, k, k)
     alhv(1, j, 1, k, :)=0.0
     alhv(1, j, 5, k, k)=alhv(1, j, 5, k, k)-alhv(1, j, 2, k, k)
     alhv(1, j, 2, k, :)=0.0
     ! -> u(1, i)=0
     aluu(1, j, 1, k, :)=0.0
     aluu(1, j, 2, k, :)=0.0
     alvu(1, j, 1, k, :)=0.0
     alvu(1, j, 2, k, :)=0.0
     alhu(1, j, 1, k, :)=0.0
     alhu(1, j, 2, k, :)=0.0
     ! -> h(1, i)=0
     aluh(1, j, 1, k, :)=0.0
     aluh(1, j, 2, k, :)=0.0
     alvh(1, j, 1, k, :)=0.0
     alvh(1, j, 2, k, :)=0.0
     alhh(1, j, 1, k, :)=0.0
     alhh(1, j, 2, k, :)=0.0

     cpnluh(1, j, 1, k, :)=0.0
     cpnluh(1, j, 2, k, :)=0.0
     cpnlvh(1, j, 1, k, :)=0.0
     cpnlvh(1, j, 2, k, :)=0.0

     ! x = 1:
     ! -> v(n, i)+v(n+1, i)=0
     alvv(n, j, 4, k, k)=alvv(n, j, 4, k, k)-alvv(n, j, 7, k, k)
     alvv(n, j, 7, k, :)=0.0
     alvv(n, j, 3, k, k)=alvv(n, j, 3, k, k)-alvv(n, j, 6, k, k)
     alvv(n, j, 6, k, :)=0.0
     alhv(n, j, 4, k, k)=alhv(n, j, 4, k, k)-alhv(n, j, 7, k, k)
     alhv(n, j, 7, k, :)=0.0
     alhv(n, j, 3, k, k)=alhv(n, j, 3, k, k)-alhv(n, j, 6, k, k)
     alhv(n, j, 6, k, :)=0.0
     ! -> u(n, i)=0
     aluu(n, j, :,k, :)=0.0
     aluu(n, j, 4, k, k)=0.0
     aluv(n, j, :,k, :)=0.0
     aluh(n, j, :,k, :)=0.0
     ! -> h(n, i)=0
     aluh(n, j, 6, k, :)=0.0
     aluh(n, j, 7, k, :)=0.0
     alvh(n, j, 6, k, :)=0.0
     alvh(n, j, 7, k, :)=0.0
     alhh(n, j, 6, k, :)=0.0
     alhh(n, j, 7, k, :)=0.0

     cpnluh(n, j, 6, k, :)=0.0
     cpnluh(n, j, 7, k, :)=0.0
     cpnlvh(n, j, 6, k, :)=0.0
     cpnlvh(n, j, 7, k, :)=0.0

     enddo

     ! y = 0:
     ! -> u(i, 0)+u(i, 1)=0
     aluu(:,1, 4, k, k)=aluu(:,1, 4, k, k)-aluu(:,1, 3, k, k)
     aluu(:,1, 3, k, :)=0.0
     aluu(:,1, 7, k, k)=aluu(:,1, 7, k, k)-aluu(:,1, 6, k, k)
     aluu(:,1, 6, k, :)=0.0
     alvu(:,1, 4, k, k)=alvu(:,1, 4, k, k)-alvu(:,1, 3, k, k)
     alvu(:,1, 3, k, :)=0.0
     alvu(:,1, 7, k, k)=alvu(:,1, 7, k, k)-alvu(:,1, 6, k, k)
     alvu(:,1, 6, k, :)=0.0
     alhu(:,1, 4, k, k)=alhu(:,1, 4, k, k)-alhu(:,1, 3, k, k)
     alhu(:,1, 3, k, :)=0.0
     alhu(:,1, 7, k, k)=alhu(:,1, 7, k, k)-alhu(:,1, 6, k, k)
     alhu(:,1, 6, k, :)=0.0
     ! -> v(i, 0)=0
     aluv(:,1, 3, k, :)=0.0
     aluv(:,1, 6, k, :)=0.0
     alvv(:,1, 3, k, :)=0.0
     alvv(:,1, 6, k, :)=0.0
     alhv(:,1, 3, k, :)=0.0
     alhv(:,1, 6, k, :)=0.0
     ! -> h(i, 0)=0
     aluh(:,1, 3, k, :)=0.0
     aluh(:,1, 6, k, :)=0.0
     alvh(:,1, 3, k, :)=0.0
     alvh(:,1, 6, k, :)=0.0
     alhh(:,1, 3, k, :)=0.0
     alhh(:,1, 6, k, :)=0.0

     cpnluh(:,1, 3, k, :)=0.0
     cpnluh(:,1, 6, k, :)=0.0
     cpnlvh(:,1, 3, k, :)=0.0
     cpnlvh(:,1, 6, k, :)=0.0

     ! y = 1:
     ! -> u(i, m)+u(i, m+1)=0	
     aluu(:,m, 4, k, k)=aluu(:,m, 4, k, k)-aluu(:,m, 5, k, k)
     aluu(:,m, 5, k, :)=0.0
     aluu(:,m, 1, k, k)=aluu(:,m, 1, k, k)-aluu(:,m, 2, k, k)
     aluu(:,m, 2, k, :)=0.0
     alhu(:,m, 4, k, k)=alhu(:,m, 4, k, k)-alhu(:,m, 5, k, k)
     alhu(:,m, 5, k, :)=0.0
     alhu(:,m, 1, k, k)=alhu(:,m, 1, k, k)-alhu(:,m, 2, k, k)
     alhu(:,m, 2, k, :)=0.0
     ! -> v(i, m)=0
     alvu(:,m, :,k, :)=0.0
     alvv(:,m, :,k, :)=0.0
     alvv(:,m, 4, k, k)=0.0
     alvh(:,m, :,k, :)=0.0
     ! -> h(i, m)=0
     aluh(:,m, 2, k, :)=0.0
     aluh(:,m, 5, k, :)=0.0
     alvh(:,m, 2, k, :)=0.0
     alvh(:,m, 5, k, :)=0.0
     alhh(:,m, 2, k, :)=0.0
     alhh(:,m, 5, k, :)=0.0

     cpnluh(:,m, 2, k, :)=0.0
     cpnluh(:,m, 5, k, :)=0.0
     cpnlvh(:,m, 2, k, :)=0.0
     cpnlvh(:,m, 5, k, :)=0.0

     enddo



     if (pper) then
     ms = mper
     llvu(1, mper, 2, :,:) = 0.0
   else
     ms = 1
   endif	
     do k = 1, nl
     do j = ms, m
     lluv(1, j, 4, k, k)=lluv(1, j, 4, k, k)-lluv(1, j, 1, k, k)
     lluv(1, j, 1, k, :)=0.0
     lluv(1, j, 5, k, k)=lluv(1, j, 5, k, k)-lluv(1, j, 2, k, k)
     lluv(1, j, 2, k, :)=0.0
     llvv(1, j, 4, k, k)=llvv(1, j, 4, k, k)-llvv(1, j, 1, k, k)
     llvv(1, j, 1, k, :)=0.0
     llvv(1, j, 5, k, k)=llvv(1, j, 5, k, k)-llvv(1, j, 2, k, k)
     llvv(1, j, 2, k, :)=0.0
     llhv(1, j, 4, k, k)=llhv(1, j, 4, k, k)-llhv(1, j, 1, k, k)
     llhv(1, j, 1, k, :)=0.0
     llhv(1, j, 5, k, k)=llhv(1, j, 5, k, k)-llhv(1, j, 2, k, k)
     llhv(1, j, 2, k, :)=0.0
     ! -> u(1, i)=0
     lluu(1, j, 1, k, :)=0.0
     lluu(1, j, 2, k, :)=0.0
     llvu(1, j, 1, k, :)=0.0
     llvu(1, j, 2, k, :)=0.0
     llhu(1, j, 1, k, :)=0.0
     llhu(1, j, 2, k, :)=0.0
     ! -> h(1, i)=0
     lluh(1, j, 1, k, :)=0.0
     lluh(1, j, 2, k, :)=0.0
     llvh(1, j, 1, k, :)=0.0
     llvh(1, j, 2, k, :)=0.0
     llhh(1, j, 1, k, :)=0.0
     llhh(1, j, 2, k, :)=0.0

     ! x = 1:
     ! -> v(n, i)+v(n+1, i)=0
     llvv(n, j, 4, k, k)=llvv(n, j, 4, k, k)-llvv(n, j, 7, k, k)
     llvv(n, j, 7, k, :)=0.0
     llvv(n, j, 3, k, k)=llvv(n, j, 3, k, k)-llvv(n, j, 6, k, k)
     llvv(n, j, 6, k, :)=0.0
     llhv(n, j, 4, k, k)=llhv(n, j, 4, k, k)-llhv(n, j, 7, k, k)
     llhv(n, j, 7, k, :)=0.0
     llhv(n, j, 3, k, k)=llhv(n, j, 3, k, k)-llhv(n, j, 6, k, k)
     llhv(n, j, 6, k, :)=0.0
     ! -> u(n, i)=0
     lluu(n, j, :,k, :)=0.0
     lluu(n, j, 4, k, k)=0.0
     lluv(n, j, :,k, :)=0.0
     lluh(n, j, :,k, :)=0.0
     ! -> h(n, i)=0
     lluh(n, j, 6, k, :)=0.0
     lluh(n, j, 7, k, :)=0.0
     llvh(n, j, 6, k, :)=0.0
     llvh(n, j, 7, k, :)=0.0
     llhh(n, j, 6, k, :)=0.0
     llhh(n, j, 7, k, :)=0.0
     enddo

     ! y = 0:
     ! -> u(i, 0)+u(i, 1)=0
     lluu(:,1, 4, k, k)=lluu(:,1, 4, k, k)-lluu(:,1, 3, k, k)
     lluu(:,1, 3, k, :)=0.0
     lluu(:,1, 7, k, k)=lluu(:,1, 7, k, k)-lluu(:,1, 6, k, k)
     lluu(:,1, 6, k, :)=0.0
     llvu(:,1, 4, k, k)=llvu(:,1, 4, k, k)-llvu(:,1, 3, k, k)
     llvu(:,1, 3, k, :)=0.0
     llvu(:,1, 7, k, k)=llvu(:,1, 7, k, k)-llvu(:,1, 6, k, k)
     llvu(:,1, 6, k, :)=0.0
     llhu(:,1, 4, k, k)=llhu(:,1, 4, k, k)-llhu(:,1, 3, k, k)
     llhu(:,1, 3, k, :)=0.0
     llhu(:,1, 7, k, k)=llhu(:,1, 7, k, k)-llhu(:,1, 6, k, k)
     llhu(:,1, 6, k, :)=0.0
     ! -> v(i, 0)=0
     lluv(:,1, 3, k, :)=0.0
     lluv(:,1, 6, k, :)=0.0
     llvv(:,1, 3, k, :)=0.0
     llvv(:,1, 6, k, :)=0.0
     llhv(:,1, 3, k, :)=0.0
     llhv(:,1, 6, k, :)=0.0
     ! -> h(i, 0)=0
     lluh(:,1, 3, k, :)=0.0
     lluh(:,1, 6, k, :)=0.0
     llvh(:,1, 3, k, :)=0.0
     llvh(:,1, 6, k, :)=0.0
     llhh(:,1, 3, k, :)=0.0
     llhh(:,1, 6, k, :)=0.0

     ! y = 1:
     ! -> u(i, m)+u(i, m+1)=0	
     lluu(:,m, 4, k, k)=lluu(:,m, 4, k, k)-lluu(:,m, 5, k, k)
     lluu(:,m, 5, k, :)=0.0
     lluu(:,m, 1, k, k)=lluu(:,m, 1, k, k)-lluu(:,m, 2, k, k)
     lluu(:,m, 2, k, :)=0.0
     llhu(:,m, 4, k, k)=llhu(:,m, 4, k, k)-llhu(:,m, 5, k, k)
     llhu(:,m, 5, k, :)=0.0
     llhu(:,m, 1, k, k)=llhu(:,m, 1, k, k)-llhu(:,m, 2, k, k)
     llhu(:,m, 2, k, :)=0.0
     ! -> v(i, m)=0
     llvu(:,m, :,k, :)=0.0
     llvv(:,m, :,k, :)=0.0
     llvv(:,m, 4, k, k)=0.0
     llvh(:,m, :,k, :)=0.0
     ! -> h(i, m)=0
     lluh(:,m, 2, k, :)=0.0
     lluh(:,m, 5, k, :)=0.0
     llvh(:,m, 2, k, :)=0.0
     llvh(:,m, 5, k, :)=0.0
     llhh(:,m, 2, k, :)=0.0
     llhh(:,m, 5, k, :)=0.0
     enddo

     end
     
!====================================================================================

subroutine land2mat_bilin(F)
  !implicit none
  use m_par
  use m_usr
  use m_mat
  implicit none

  !*     IMPORT/EXPORT
  real    F(ndim)
  !*     LOCAL
  integer idx, v, row, k, mod3
  !*
  !*     par(10) = 1  land is 'down'
  !*     par(10) = 0  land is 'up'
  do idx = 1, ndir
    row = dir(idx)
  do v = bega(row), bega(row+1)-1
    coa(v) = par(10)*coa(v)
    if (jcoa(v).eq.row) coa(v) = coa(v)
  enddo
  if (mod(row, 3).ne.0) then
    F(row) = par(10)*F(row)
  else
    do k = 1, nl-1
      if ((row .gt. 3*(k-1)*n*m).and.(row .le. 3*k*n*m)) then
        F(row) = par(10)*F(row)
      endif
    enddo  
    if (row .gt. 3*(nl-1)*n*m) then
      F(row) = par(10)*F(row)
    endif 
  endif
  enddo
  END

!===================================================================================
      subroutine nluh_bilin(hdet, h)
!     produce local matrices for nonlinear operators for calc of rhs

     use m_par
     use m_usr
     implicit none
!     import/export
      real    h(0:n+1, 0:m+1, 1:nl), hdet(0:n+1, 0:m+1, 1:nl)
!     local
      real    txh(n, m, np)
      real    sig, eps2
      integer i, j, k

      sig  = par(1)*par(11)
      eps2 = par(3)
! TODO KURT: Was it correct to disable the next call?
!      call forcing

      nluh = 0.0

      do k = 1, nl
        txh = 0.0

        if (k .eq. 1) then
         DO j = 1, m
          DO i = 1, n-1
           !txh(i, j, 4) = (h(i, j)+h(i+1, j))*(tx(i, j)+tx(i+1, j))/(hdet(i, j)+hdet(i+1, j))**3
           !txh(i, j, 7) = (h(i, j)+h(i+1, j))*(tx(i, j)+tx(i+1, j))/(hdet(i, j)+hdet(i+1, j))**3
           txh(i, j, 4) = h(i, j, k)*(tx(i, j)+tx(i+1, j))/(hdet(i, j, k)+hdet(i+1, j, k))**3
           txh(i, j, 7) = h(i, j, k)*(tx(i, j)+tx(i+1, j))/(hdet(i, j, k)+hdet(i+1, j, k))**3
          ENDDO
          txh(n, j, 4) = h(n, j, k)*(tx(n, j)+tx(1, j))/(hdet(n, j, k)+hdet(1, j, k))**3
          txh(n, j, 7) = h(n, j, k)*(tx(n, j)+tx(1, j))/(hdet(n, j, k)+hdet(1, j, k))**3
         ENDDO
          nluh(:,:,:,k, k) = nluh(:,:,:,k, k)-sig*txh
        endif
        enddo 
        !print *, "norm of tx = ", Norm2(tx)
        !print *, "norm of txh = ", Norm2(txh)
        !print *, "norm of h = ", Norm2(h)
        !print *, "norm of hdet = ", Norm2(hdet), sig
      end


!=======================================================================================
      subroutine bilin(udet, un, vn, bil) bind(C, name='model_bil')

        use m_par
        use m_mat
        use m_usr
    
        implicit none
    
!*     import/export
        real(kind = c_double), dimension(ndim), intent(in) ::    un, vn, udet
        real(kind = c_double), dimension(ndim), intent(out) ::    bil
        integer(kind = c_int) ::    i, j, k, row
        real(kind = c_double)::    u(0:n, 0:m+1, 1:nl), v(0:n+1, 0:m, 1:nl), h(0:n+1, 0:m+1, 1:nl)
        real(kind = c_double)::    u1(0:n, 0:m+1, 1:nl), v1(0:n+1, 0:m, 1:nl), h1(0:n+1, 0:m+1, 1:nl)
        !real(kind = c_double) ::    hbx(n, m), hby(n, m), eps2, mix(ndim)
        !           eps2 = par(3)
        !   print *, "in SWE_interface.F90_bilin n, m, nl =",n, m, nl
        !
        !c      call depth_cor(un)
        !print *, "norm of udet =", Norm2(udet)
        !print *, "norm of un =", Norm2(un)
        !print *, "norm of vn =", Norm2(vn)
        call usol(udet, u1, v1, h1)      ! conversion state -> u, v, h
        call solu(udet, u1, v1, h1)      ! conversion state -> u, v, h
        call usol(vn, u, v, h)      ! conversion state -> u, v, h
        call solu(vn, u, v, h)
        call usol(un, u, v, h)      ! conversion state -> u, v, h
        call solu(un, u, v, h)
    
        call nlin_rhs(u, v, h)      ! local matrix, non-linear
        call nluh_bilin(h1, h)
        !call nlin_rhs(u1, v1, h1)      ! local matrix, non-linear
        !if (Norm2(un) > Norm2(vn)) then
        !call usol(udet, u1, v1, h1)      ! conversion state -> u, v, h
        !call solu(udet, u1, v1, h1)      ! conversion state -> u, v, h
          !call usol(un, u1, v1, h1)      ! conversion state -> u, v, h
          !call solu(un, u1, v1, h1)      ! conversion state -> u, v, h
          !call usol(vn, u, v, h)      ! conversion state -> u, v, h
          !call nluh_bilin(h1, h)
          !aluh = nluh
        !else
          !call usol(vn, u1, v1, h1)    ! conversion state -> u, v, h 
          !call solu(vn, u1, v1, h1)    ! conversion state -> u, v, h
          !call usol(un, u, v, h)
          !aluh = 0
        !endif
        !print *, "norm of nluh =", Norm2(nluh)

        !Print *,'lluu=', lluu
        aluu =  nluu        ! total local matrix
        aluv =  nluv
        aluh =  nluh
        alvu =  nlvu
        alvv =  nlvv
        alvh =  0
        !alvh =  nlvh
        alhu =  nlhu
        alhv =  nlhv
        alhh =  nlhh
        !print *, "swe_interface: norm of nlhh ", Norm2(nlhh)
        !print *, nlhu(1, 2, :,:,:)+nlhv(1, 2, :,:,:)+nlhh(1, 2, :,:,:)
        !print *, nlhh(1, 2, :,:,:)
        call boundaries_bilin           ! implement boundary conditions
        call assemble             ! assemble global matrix from local
        !print *, "assemble: norm of coA", Norm2(coA)
    
        if (mix_flag) then
          call mix_jac(un)
          stop  "WARNING: mix_flag is on!!!!"
        endif
      !*************************************************************************
      ! The following line is a hack and would require the change if 
      ! 'cont' is false (which makes par(10)=1). My current guess is 
      ! we might have to silence all lines which have the followng 
      ! F(row) = par(10)*F(row)
      ! It is to pass the following test for Jacobian and bilnear term relation 
      ! J(0)*v + <x, v> + <v, x> + J(x)*v
      !*************************************************************************
        call land2mat_bilin(frc)        ! remove land
        !print *, "land2mat: norm of coA", Norm2(coA)
        if (pmode .eq. 3) call intcond(Frc)
        !print *, "intcond: norm of coA", Norm2(coA)
        call packa                ! sparse row storage and removal of small elements
        !print *, "packa: norm of coA", Norm2(coA)
        !stop
    
        !    print *, cpnluh
        ! print *, "sig = ",sig
        !   pause
        !print *, "bil = ", bil
        !print *, "inside bilin func"
        !call writemat
        !pause
        call matavec(vn, bil)
        !print *, "bil = ", bil
        !
        !Print *, 'norm of aluu=',Norm2(aluu)
        !Print *, 'norm of aluv=',Norm2(aluv)
        !Print *, 'norm of aluh=',Norm2(aluh)
        !Print *, 'norm of alvu=',Norm2(alvu)
        !Print *, 'norm of alvv=',Norm2(alvv)
        !Print *, 'norm of alvh=',Norm2(nlvh)
        !Print *, 'norm of alhu=',Norm2(alhu)
        !Print *, 'norm of alhv=',Norm2(alhv)
        !Print *, 'norm of alhh=',Norm2(alhh)
        bil = -bil
        !do i = 1, 30
        !print *, bil(i)
        !end do
        100   format(3i6, 2es12.4)
        end

!=================================================================
    subroutine matrix(nrows, nzmax, val, rows, cols, x_in, ierr) &
        bind(C, name='model_jac')
      use m_Tools, only: toString
      !implicit none
      use m_par
      use m_mat
      use m_usr
      implicit none
      integer(kind = c_int), intent(in):: nrows, nzmax
      integer(kind = c_int), dimension(nrows+1), intent(out):: rows
      integer(kind = c_int), dimension(nzmax), intent(out):: cols
      real(kind = c_double), dimension(nzmax), intent(out):: val
      real(kind = c_double), dimension(nrows), intent(in):: x_in
      integer(kind = c_int), intent(out):: ierr
      !real sig
      !integer type
      !#ifdef DEBUGGING
      !character(len = 4):: label
      !#endif
      integer:: i, j, k, nnz


      call matrix_swe(x_in, 1, 0.0)

      !nzmax = begA(ndim+1)-1
      !PRINT *, nrows, ndim, nzmax, begA(ndim+1)-1
      !PRINT *, cols(10713), val(10713)
      !PRINT *, (jcoA(begA(i)), i = 5398, 5405)
      !print *, (begA(i), i = 5398, 5405)
      !return
      nnz = begA(ndim+1)-1

      rows = begA
      val(1:nnz)=-coA(1:nnz)
      cols(1:nnz)=jcoA(1:nnz)

      ierr = 0
    end subroutine matrix
    !=======================================================================================
    subroutine matrixlin(nrows, nzmax, ierr, linval, linrows, lincols) &
        bind(C, name='model_jaclin')
      !We assume that this call is after that of "matrix"
      use m_Tools, only: toString
      !implicit none
      use m_par
      use m_mat
      use m_usr
      implicit none
      integer(kind = c_int), intent(in):: nrows, nzmax
      integer(kind = c_int), dimension(nrows+1), intent(out):: linrows
      integer(kind = c_int), dimension(nzmax), intent(out) ::  lincols
      real(kind = c_double), dimension(nzmax), intent(out):: linval
      integer(kind = c_int), intent(out):: ierr
      !real sig
      !integer type
      !#ifdef DEBUGGING
      !character(len = 4):: label
      !#endif
      integer:: i, j, k, nnz

      ! call matrix_swe(x_in, 1, 0.0)
      nnz = ccbegA(ndim+1)-1
      linrows = ccbegA
      linval(1:nnz)=-ccoA(1:nnz)
      lincols(1:nnz)=ccjcoA(1:nnz)
    end subroutine matrixlin
    !==============================================================================
    subroutine massmat(nrows, x_in, massdiag) &
        bind(C, name='model_massmat')
      implicit none
      integer(kind = c_int), intent(in):: nrows
      real(kind = c_double), dimension(nrows), intent(in):: x_in
      real(kind = c_double), dimension(nrows), intent(out):: massdiag
      massdiag(1:nrows) = mass(1:nrows)
    end subroutine massmat
    !==============================================================================
    subroutine mkrhs(nrows, rhs, x_in, ierr) bind(C, name='model_rhs')
      !This routine has much resemblance to the assemble routine.
      !An explicit conversion to U, V, W, P is not necessary ones everything is defined in terms of stencils.

      implicit none

      integer(kind = c_int), intent(in):: nrows
      real(kind = c_double), dimension(nrows), intent(out):: rhs
      real(kind = c_double), dimension(nrows), intent(in):: x_in
      integer(kind = c_int), intent(out):: ierr
      ! local
      integer row
      integer:: line, plane
      ! cout<<i <<"\n"<< rows[i] << "\t" << rows[i+1] <<"\t"<< nrows<<"\t"<<cols[j] i
      !print *, 'x_in is : ',x_in
      ierr = 1
      call mrhs(x_in, rhs)
      !print *, 'norm of fortran rhs ',Norm2(rhs)
      !rhs = rhs
      ierr = 0
    end subroutine mkrhs
    !==============================================================================
    subroutine mkrhs_s(nrows, nzmax, stress_term, srows, scols, ierr) bind(C, name='model_rhs_s')
      !Assumed to be called immedeately after mkrhs, if called.
      !This routine has much resemblance to the assemble routine.
      !An explicit conversion to U, V, W, P is not necessary ones everything is defined in terms of stencils.
      use m_par
      use m_mat
      use m_usr

      implicit none

      integer(kind = c_int), intent(in):: nrows, nzmax
      integer(kind = c_int), intent(out):: ierr
      integer(kind = c_int), dimension(nrows+1), intent(out):: srows
      integer(kind = c_int), dimension(nzmax), intent(out):: scols
      real(kind = c_double), dimension(nzmax), intent(out):: stress_term

      ! local
      integer row, nnz
      ierr = 1

      !call mrhs(x_in, rhs, stress_term, cpjcoA, cpbegA)

      nnz = cpbegA(ndim+1)-1
      srows = cpbegA
      stress_term(1:nnz)=-cpcoA(1:nnz)
      scols(1:nnz)=cpjcoA(1:nnz)
      ierr = 0
    end subroutine mkrhs_s

    subroutine memory_estimate(ints, doubles) bind(C, name='model_memory_estimate')
      use, intrinsic:: iso_c_binding
      implicit none

      REAL(kind = c_double), intent(out):: ints, doubles

      call mem_required(Model, ints, doubles)

    end subroutine memory_estimate
    !=================================================================================

    subroutine StochFrc(nrows, ret, col, ierr) bind(C, name='model_stoch_frc')
      use m_usr
      implicit none

      integer(kind = c_int), intent(in):: nrows, col
      real(kind = c_double), dimension(nrows), intent(out):: ret
      integer(kind = c_int), intent(out):: ierr
      ! local
      real(kind = c_double):: etax, etay, l, C, xc, yc
      integer(kind = c_int):: i, j, row
      l = 0.125; 
      C = 0.1
      ierr = 1
      ! We expect to column numbers: 0 and 1. These are used to define the tow vectors
      !etax = C*col; 
      !etay = C*(1-col); 
      xc=(xgmin+xgmax)/2
      yc=(ygmin+ygmax)/2
      row=(nz-1)*nx*ny*3+1 
      ret = 0
      !if (nz .ne. 2) STOP 'number of layers not two'
      !open(23, file = 'data1.dat', status = 'new')
      do j = 1, ny
      do i = 1, nx
      !ret(row)= C*exp(-2. * ((x(i)-xc) * (x(i)-xc) + (y(j)-yc) * (y(j)-yc)) / (4. * l*l))*tx(i, j)
      ret(row)= C*exp(-2. * ((x(i)-xc) * (x(i)-xc) + (y(j)-yc) * (y(j)-yc)) / (4. * l*l))
      ret(row+1)=0
      ret(row+2)=0
      !  write(23, *) ret(row)
      row = row+3; 
      enddo
      enddo
      !close(23)
      ierr = 0
    end subroutine StochFrc

    !=================================================================================
    subroutine get_num_stochvectors(nVec) bind(C, name='model_get_num_stochvectors')

      implicit none
      integer(kind = c_int), intent(out):: nVec
      nVec = 1

    end subroutine get_num_stochvectors


  end module m_LidDrivenCavity

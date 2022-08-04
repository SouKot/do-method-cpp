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
!=======================================================================================
subroutine bilin(un, vn, bil) bind(C, name='model_bil')

use m_par
use m_mat
use m_usr

implicit none

!*     import/export
  real(kind = c_double), dimension(ndim), intent(in) ::    un, vn
  real(kind = c_double), dimension(ndim), intent(out) ::    bil
  integer(kind = c_int) ::    i, j, k, row
  real(kind = c_double)::    u(0:n, 0:m+1, 1:nl), v(0:n+1, 0:m, 1:nl), h(0:n+1, 0:m+1, 1:nl)
  !real(kind = c_double) ::    hbx(n, m), hby(n, m), eps2, mix(ndim)
!           eps2 = par(3)
  !   print *, "in SWE_interface.F90_bilin n, m, nl =",n, m, nl
  !   print *, "HELOO!!!"
  !   print *, " un = ", un
  !
!c      call depth_cor(un)
  call usol(un, u, v, h)      ! conversion state -> u, v, h
     !print *, "un = ",un
     !print *, "h = ",h
     !pause
  !print *,'HHH =',h
  !call lin
  call nlin_rhs(u, v, h)      ! local matrix, non-linear
  !Print *,'lluu=', lluu
  aluu =  nluu        ! total local matrix
  aluv =  nluv
  aluh =  0
  alvu =  nlvu
  alvv =  nlvv
  alvh =  0
  alhu =  nlhu
  alhv =  nlhv
  alhh =  nlhh
   ! Print *, 'norm of aluu=',Norm2(aluu)
   ! !
   ! Print *, 'norm of aluv=',Norm2(aluv)
   !
   !Print *, 'norm of aluh=',Norm2(aluh)
   !Print *, 'norm of alvu=',Norm2(alvu)
   !Print *, 'norm of alvv=',Norm2(alvv)
   !Print *, 'norm of alvh=',Norm2(nlvh)
   !Print *, 'norm of alhu=',Norm2(alhu)
   !Print *, 'norm of alhv=',Norm2(alhv)
   !Print *, 'norm of alhh=',Norm2(alhh)

  call boundaries           ! implement boundary conditions
  call assemble             ! assemble global matrix from local

  if (mix_flag) then
  call mix_jac(un)
  stop  "WARNING: mix_flag is on!!!!"
endif

  call land2mat(frc)        ! remove land
  if (pmode .eq. 3) call intcond(Frc)
  call packa                ! sparse row storage and removal of small elements
  !stop

  !    print *, cpnluh
  ! print *, "sig = ",sig
  !   pause
!print *, "bil = ", bil
  call matavec(vn, bil)
!print *, "bil = ", bil
!
 bil = -bil
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
#ifdef DEBUGGING
character(len = 4):: label
#endif
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
#ifdef DEBUGGING
character(len = 4):: label
#endif
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
do j = 1, ny
  do i = 1, nx
    !ret(row)= C*exp(-2. * ((x(i)-xc) * (x(i)-xc) + (y(j)-yc) * (y(j)-yc)) / (4. * l*l))*tx(i, j)
    ret(row)= C*exp(-2. * ((x(i)-xc) * (x(i)-xc) + (y(j)-yc) * (y(j)-yc)) / (4. * l*l))
    ret(row+1)=0
    ret(row+2)=0
    row = row+3; 
  enddo
enddo
ierr = 0
end subroutine StochFrc

!=================================================================================
subroutine get_num_stochvectors(nVec) bind(C, name='model_get_num_stochvectors')

implicit none
integer(kind = c_int), intent(out):: nVec
nVec = 1

end subroutine get_num_stochvectors


end module m_LidDrivenCavity

MODULE m_Types

use m_Tools

implicit none

INTEGER, PARAMETER :: MAX_NAME_LEN = 1024
INTEGER, PARAMETER :: DP = 8

integer :: istat

! basic description of the model (independent of the actual grid used)
TYPE ModelDescription_t

CHARACTER(LEN=MAX_NAME_LEN) :: Label
INTEGER :: ndim ! number of spatial dimensions
INTEGER :: ndof ! number of degrees of freedom per grid cell
INTEGER :: npar=0 ! number of model parameters

INTEGER :: pid=0 ! global processor rank (for parallel computations,
                 ! can be useful when writing files in parallel etc.)

END TYPE ModelDescription_t

! basic type for model parameters such as the Reynolds number
TYPE Param_t

CHARACTER(len=MAX_NAME_LEN) :: name ! name of the parameter

INTEGER :: idx  ! position in parameter vector
REAL(KIND=dp) :: value=0.0
REAL(KIND=dp) :: default_value=0.0

END TYPE Param_t

! basic type for defining a stencil in each grid point
TYPE Stencil_t

integer :: nx=0, ny=0, nz=0 ! number of grid cells
integer :: np          ! max number of nonzero neighbours
integer :: ndof       ! number of dofs per grid cell

integer :: stxb,stxe,styb,stye,stzb,stze ! stencil width, defaults to [-1, +1]

! stencil for the linear part of the problem
! The dimensions are Al(1:nx,1:ny,1:nz,1:np,1:ndof,1:ndof)
REAL(KIND=dp), pointer, dimension(:,:,:,:,:,:) :: Al

! stencil for the bilinear part
! The dimensions are BiL(1:nx,1:ny,1:nz,1:np,0:npBiL-1,1:2*(ndof-1),1:ndof-1)
REAL(KIND=dp), pointer, dimension(:,:,:,:,:,:) ::  BiL

! stencil for the additional terms in F
! The dimensions are AnlF(1:nx,1:ny,1:nz,1:np,1:ndof-1,1:ndof-1)
REAL(KIND=dp), pointer, dimension(:,:,:,:,:,:) :: AnlF

! stencil for the additional terms in the jacobian
! The dimensions are Al(1:nx,1:ny,1:nz,1:np,1:ndof-1,1:ndof-1)
REAL(KIND=dp), pointer, dimension(:,:,:,:,:,:) :: AnlJ

! vector to store diagonal mass matrix
REAL(KIND=dp), pointer, dimension(:) :: Mass

END TYPE Stencil_t

TYPE Forcing_t
integer :: nx=0, ny=0, nz=0 ! number of grid cells
integer :: ndof       ! number of dofs per grid cell
! Forcing array  
! The dimensions are Frc(0:nx,0:ny,0:nz,1:ndof)
REAL(KIND=dp), pointer, dimension(:,:,:,:) :: Frc

END TYPE Forcing_t

TYPE Matrix_t

! standard sparse matrix in CRS format with 
! memory allocated for nzmax nonzeros and nnz
! nonzeros actually filled in.
integer :: nrows, nnz, nzmax
integer, pointer, dimension(:) :: rows
integer, pointer, dimension(:) :: cols
real(kind=dp), pointer, dimension(:) :: values


END TYPE Matrix_t

TYPE Grid_t

! number of grid cells
integer nx,ny,nz

! vertex coordinates
real(kind=dp), pointer, dimension(:) :: x,y,z

! indicates what kind of properties the cell has (solid, fluid...)
integer, pointer, dimension(:,:,:) :: material

END TYPE Grid_t

TYPE ScalarField_t

! index bounds
integer :: i0,i1,j0,j1,k0,k1

! scalar field
real(kind=dp), pointer, dimension(:,:,:) :: values

END TYPE ScalarField_t

! a type to contain the complete model description
TYPE Model_t

TYPE(ModelDescription_t) :: header

TYPE(Param_t), POINTER, DIMENSION(:) :: param

TYPE(Grid_t) :: grid

TYPE(Stencil_t) :: Stencil

TYPE(Forcing_t) :: Forcing

TYPE(ScalarField_t), POINTER, DIMENSION(:) :: Field

logical :: initialized = .false.

END TYPE Model_t

CONTAINS

! allocate memory and set parameters in all the data structures of the model.
! input: model, string to describe it, dimension, dof/cell,number of parameters, 
! grid size in 3D, max num neighbours in stencil
SUBROUTINE AllocateModel(M,ndim,ndof,npar,nx,ny,nz,np,npBiL)

#ifdef HAVE_MPI
use MPI
#endif
implicit none
TYPE(Model_t) :: M
integer, intent(in) :: ndim,ndof,npar,nx,ny,nz,np,npBiL

integer :: i

M%Header%label = '<model description not set>'
M%Header%ndim = ndim
M%Header%ndof = ndof
M%Header%npar = npar

#ifdef HAVE_MPI
call MPI_Comm_rank(MPI_COMM_WORLD,M%Header%pid, istat)
#else
M%Header%pid = 0
#endif
allocate(M%param(npar))
allocate(M%Field(ndof))
do i=1,ndof
  nullify(M%Field(i)%values)
end do
call AllocateGrid(M%Grid,nx,ny,nz)
call AllocateStencil(M%Stencil,nx,ny,nz,ndof,np,npBiL)
call AllocateForcing(M%Forcing,nx,ny,nz,ndof)

if (istat/=0) then
  call Error('Memory allocation failed!',__FILE__,__LINE__)
else
  M%initialized=.true.
end if
END SUBROUTINE AllocateModel


SUBROUTINE AllocateGrid(G,nx,ny,nz)

TYPE(Grid_t) :: G
integer :: nx,ny,nz

G%nx = nx
G%ny = ny
G%nz = nz

ALLOCATE(G%x(-1:nx+1),stat=istat)
ALLOCATE(G%y(-1:ny+1),stat=istat)
ALLOCATE(G%z(-1:nz+1),stat=istat)

ALLOCATE(G%Material(0:nx+1,0:ny+1,0:nz+1),stat=istat)

END SUBROUTINE AllocateGrid

SUBROUTINE AllocateStencil(S,nx,ny,nz,ndof,np,npBiL)

TYPE(Stencil_t) :: S
integer :: nx,ny,nz,ndof, np, npBiL

S%nx=nx
S%ny=ny
S%nz=nz

S%np=np
S%ndof=ndof

S%stxb=-1
S%stxe=+1

S%styb=-1
S%stye=+1

S%stzb=-1
S%stze=+1

ALLOCATE(S%Al(nx,ny,nz,1:np,ndof,ndof),stat=istat)
ALLOCATE(S%AnlF(nx,ny,nz,1:np,ndof-1,ndof-1),stat=istat)
ALLOCATE(S%AnlJ(nx,ny,nz,1:np,ndof-1,ndof-1),stat=istat)
ALLOCATE(S%BiL(nx,ny,nz,0:npBiL-1,2*(ndof-1),ndof-1),stat=istat)
ALLOCATE(S%Mass(ndof*nx*ny*nz),stat=istat)


END SUBROUTINE AllocateStencil

subroutine AllocateForcing(F,nx,ny,nz,ndof)

TYPE(Forcing_t) :: F
integer :: nx,ny,nz,ndof

F%nx=nx
F%ny=ny
F%nz=nz

F%ndof=ndof
ALLOCATE(F%Frc(nx,ny,nz,ndof),stat=istat)

F%Frc(:,:,:,:) = 0.D0

END SUBROUTINE AllocateForcing

SUBROUTINE AllocateField(F,i0,i1,j0,j1,k0,k1)

IMPLICIT NONE

TYPE(ScalarField_t) :: F

INTEGER, intent(in) :: i0,i1,j0,j1,k0,k1

F%i0=i0
F%i1=i1
F%j0=j0
F%j1=j1
F%k0=k0
F%k1=k1

ALLOCATE(F%values(i0:i1,j0:j1,k0:k1),stat=istat)

END SUBROUTINE AllocateField

subroutine AllocateMatrix(A,nrows,nzmax)

TYPE(Matrix_t) :: A
integer :: nrows,nzmax

A%nrows = nrows
A%nzmax = nzmax
A%nnz = 0
ALLOCATE(A%rows(nrows+1), &
         A%cols(nzmax), &
         A%values(nzmax), stat=istat)

end subroutine AllocateMatrix

! create a 'view' of the matrix, e.g. set the internal
! pointers in Matrix_t A to those passed in. This is useful
! e.g. for filling C++ arrays using functions that accept an
! Matrix_t opbject.
subroutine CreateView(A,nrows,nzmax,beg,co,jco)

TYPE(Matrix_t) :: A
integer :: nrows,nzmax
integer, dimension(nrows+1), target :: beg
integer, dimension(nzmax), target :: jco
real(kind=dp), dimension(nzmax), target :: co

A%nrows = nrows
A%nzmax = nzmax
A%nnz = beg(nrows+1)
A%rows => beg
A%cols => jco
A%values => co

end subroutine CreateView

subroutine AddModelParameter(M,idx,name,default)

TYPE(Model_t) :: M
integer, intent(in) :: idx
character(len=*),intent(in) :: name
real(kind=DP),intent(in) :: default

if (.not. M%initialized) then
  call Error('model not initialized',__FILE__,__LINE__)
  return
end if

if (idx>M%Header%npar) return

M%Param(idx)%idx=idx
M%Param(idx)%name=name
M%Param(idx)%value=default
M%Param(idx)%default_value=default

end subroutine AddModelParameter

subroutine DumpMatrix(A,fid)

implicit none

TYPE(Matrix_t), intent(in) :: A
integer,intent(in) :: fid
integer :: i,j

write(fid,*) 'n=',A%nrows
write(fid,*) 'nnz=',A%nnz
write(fid,*) 'nzmax=',A%nzmax

write(fid,*) 'DAT=[...'
do i=1,A%nrows
  do j=A%rows(i),A%rows(i+1)-1
    write(fid,*) i,A%cols(j),A%values(j)
  end do
end do
write(fid,*) '];'
write(fid,*) 'A=sparse(DAT(:,1),DAT(:,2),DAT(:,3));'

end subroutine DumpMatrix

subroutine DumpForcing(Frc,fid)

implicit none

TYPE(Forcing_t), intent(in) :: Frc
integer, intent(in) :: fid

integer i,j,k,xx

write(fid,*) '% forcing field'

do xx=1,Frc%ndof
  write(42,*) '% variable ',xx
  do k=1,Frc%nz
    do j=1,Frc%ny
      write(fid,'(100E16.8)') (Frc%Frc(i,j,k,xx),i=1,Frc%nx)
    end do
  end do
end do

end subroutine DumpForcing

subroutine mem_required(M,ints,doubles)

implicit none

TYPE(Model_t) :: M

real(kind=dp) :: ints, doubles
integer :: i

ints=0.0
doubles=0.0

! we use this somewhat awkward implementation to avoid integer overflows
  if (associated(M%Stencil%Al)) then
    doubles = doubles + product(dble(shape(M%Stencil%Al)))
  end if
  if (associated(M%Stencil%BiL)) then
    doubles = doubles + product(dble(shape(M%Stencil%BiL)))
  end if
  if (associated(M%Stencil%AnlF)) then
    doubles = doubles + product(dble(shape(M%Stencil%AnlF)))
  end if
  if (associated(M%Stencil%AnlJ)) then
    doubles = doubles + product(dble(shape(M%Stencil%AnlJ)))
  end if
  if (associated(M%Stencil%Mass)) then
    doubles = doubles + product(dble(shape(M%Stencil%Mass)))
  end if

  if (associated(M%Grid%x)) then
    doubles = doubles + product(dble(shape(M%Grid%x)))
  end if
  if (associated(M%Grid%y)) then
    doubles = doubles + product(dble(shape(M%Grid%y)))
  end if
  if (associated(M%Grid%z)) then
    doubles = doubles + product(dble(shape(M%Grid%z)))
  end if
  if (associated(M%Grid%material)) then
    ints = ints + product(dble(shape(M%Grid%material)))
  end if

  if (associated(M%Forcing%Frc)) then
    doubles = doubles + product(dble(shape(M%Forcing%Frc)))
  end if

if (associated(M%Field)) then
  do i=1,size(M%Field)
    if (associated(M%Field(i)%values)) then
      doubles = doubles + product(dble(shape(M%Field(i)%values)))
    end if
  end do
end if
end subroutine mem_required

subroutine DeleteModel(M)

implicit none

TYPE(Model_t) :: M

integer :: i

!call DeleteStencil(M%Stencil)
!call DeleteGrid(M%Grid)
 call DeleteForcing(M%Forcing)
if (associated(M%Field)) then
do i=1,size(M%Field)
    call DeleteField(M%Field(i))
  end do
end if
end subroutine DeleteModel

SUBROUTINE DeleteGrid(G)

TYPE(Grid_t) :: G

if (G%nx+G%ny+G%nz==0) return;

G%nx = 0
G%ny = 0
G%nz = 0

DEALLOCATE(G%x,G%y,G%z)

DEALLOCATE(G%Material)

END SUBROUTINE DeleteGrid

SUBROUTINE DeleteStencil(S)

TYPE(Stencil_t) :: S

if (S%nx+S%ny+S%nz==0) return;

S%nx=0
S%ny=0
S%nz=0

S%np=0
S%ndof=0

S%stxb=0
S%stxe=0

S%styb=0
S%stye=0

S%stzb=0
S%stze=0

DEALLOCATE(S%Al,S%AnlF,S%AnlJ,S%BiL,S%Mass)

END SUBROUTINE DeleteStencil

subroutine DeleteForcing(F)

TYPE(Forcing_t) :: F

if (F%nx+F%ny+F%nz==0) return;

F%nx=0
F%ny=0
F%nz=0

F%ndof=0
DEALLOCATE(F%Frc)

END SUBROUTINE DeleteForcing

SUBROUTINE DeleteField(F)

IMPLICIT NONE

TYPE(ScalarField_t) :: F

F%i0=0
F%i1=0
F%j0=0
F%j1=0
F%k0=0
F%k1=0

if (associated(F%values)) DEALLOCATE(F%values)

END SUBROUTINE DeleteField

END MODULE m_Types

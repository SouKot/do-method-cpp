!* ============================================================================ *
!*
!* Multi-Layer Shallow Water Model
!*
!* Version 1.2 - 20/02/2008 - ATvS
!*
!* ============================================================================ *
!*     CONSTANT
      module m_mat
      !use m_usr
      integer   ::   adim=0
      !* parameter(adim=ndim*(3*8*nl+1))
      real(kind=8),  dimension(:),allocatable  :: rmx
      real(kind=8),   dimension(:),allocatable  :: coA,cpcoA,ccoA
      real(kind=8),   dimension(:),allocatable, target  :: coB
      integer,dimension(:),allocatable  :: jcoA,cpjcoA,ccjcoA
      integer,dimension(:),allocatable  :: jcoB
      integer,dimension(:),allocatable  :: begA,begB,cpbegA,ccbegA

      contains
  
      subroutine allocate_mat()
      use m_usr
      adim=ndim*(3*8*nl+1)
      allocate(rmx(ndim),coA(adim),coB(ndim),jcoA(adim),cpjcoA(adim),cpcoA(adim),ccjcoA(adim),ccoA(adim))
      allocate(jcoB(ndim),begA(ndim+1),begB(ndim+1),cpbegA(ndim+1),ccbegA(ndim+1))
      rmx=0;coA=0;coB=0;jcoA=0;cpjcoA=0;cpcoA=0;ccjcoA=0;ccoA=0
      jcoB=0;begA=0;begB=0;cpbegA=0;ccbegA=0
      end subroutine allocate_mat
      !common /MAT/ coA,rmx,jcoA,begA
      !common /MATB/coB,jcoB,begB
      end module m_mat

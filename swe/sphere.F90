!* ============================================================================ *
!*
!* Multi-Layer Shallow Water Model
!*
!* Version 1.2-20/02/2008-ATvS
!*
!* ============================================================================ *
!
!* dx = coef(i, j, pos, unkn)
!*
!*position: 4 central(i, j)
!*          1 west   (i-1)           2  5   
!*          7 east   (i+1)           1  4  7
!*          3 south  (j-1)              3  6
!*          5 north  (j+1)
!********************************************************
      SUBROUTINE unlin(type, atom, u, v, h)
      use m_par
      !implicit none
      use m_usr
       
!*     nonlinear terms for the u-equation                      u
!*       uu = u_ij*u_ij                                      v   v   j
!*       uv = u_ij*(v_ij+v_ij-1+v_i+1j-1+v_i+1j)*0.25      u h u h  j
!*       uux = u_ij*(u_i+1j-u_i-1j)*tdxi                    v   v   j-1
!*       vuy= (u_ij+1-u_ij-1)*(v_ij+v_ij-1+v_i+1j-1+v_i+1j)    u    j-1
!*                                                           i  i+1
!*     1:  urU                                                 i 
!*     2:  urUx /(r cos y)
!*     3:  Uurx /(r cos y)
!*     4:  vrUy/r
!*     5:  Vury/r
!*     6:  urV  (tan y)/r
!*     7:  Uvr  (tan y)/r
!*     IMPORT/EXPORT
      integer type
      real    atom(n, m, np)
      real    u(0:n, 0:m+1), v(0:n+1, 0:m), h(0:n+1, 0:m+1)
      integer i, j
      real    costdxi(m), tdyi, tanr(m)
!*
      atom = 0.0
!*
      SELECT CASE(type)
         CASE(1)            ! urU
           DO j = 1, m
            DO i = 1, n
             atom(i, j, 4) = u(i, j)
            ENDDO
           ENDDO
         CASE(2)            ! urUx
           costdxi = 1.0/(cos(y)*2*dx)
           DO j = 1, m
            DO i = 1, n-1
             atom(i, j, 4) = (u(i+1, j)-u(i-1, j))*costdxi(j)
            ENDDO
            atom(n, j, 4) = (u(1, j)-u(n-1, j))*costdxi(j)
           ENDDO
         CASE(3)            ! Uurx
           costdxi = 1.0/(cos(y)*2*dx)
           DO j = 1, m
            DO i = 1, n
             atom(i, j, 1) =-u(i, j)*costdxi(j)
             atom(i, j, 7) = u(i, j)*costdxi(j)
            ENDDO
           ENDDO
         CASE(4)            ! vrUy
           tdyi = 1.0/(8*dy)
           DO j = 1, m
            DO i = 1, n
             atom(i, j, 3) = tdyi*(u(i, j+1)-u(i, j-1))
             atom(i, j, 4) = tdyi*(u(i, j+1)-u(i, j-1))
             atom(i, j, 6) = tdyi*(u(i, j+1)-u(i, j-1))
             atom(i, j, 7) = tdyi*(u(i, j+1)-u(i, j-1))
            ENDDO
           ENDDO
         CASE(5)            ! Vury
          tdyi = 1.0/(8*dy)
          DO j = 1, m
           DO i = 1, n
            atom(i, j, 3) =-tdyi*(v(i, j-1)+v(i, j)+v(i+1, j-1)+v(i+1, j) )
            atom(i, j, 5) = tdyi*(v(i, j-1)+v(i, j)+v(i+1, j-1)+v(i+1, j) )
           ENDDO
          ENDDO
         CASE(6)            ! urV
          tanr = tan(y)/4.0
          DO j = 1, m
           DO i = 1, n
            atom(i, j, 4) =(v(i, j)+v(i, j-1)+v(i+1, j-1)+v(i+1, j))*tanr(j)
           ENDDO
          ENDDO
         CASE(7)            ! Uvr
           tanr = tan(y)/4.0
           DO j = 1, m
            DO i = 1, n
             atom(i, j, 3) = u(i, j)*tanr(j)
             atom(i, j, 4) = u(i, j)*tanr(j)
             atom(i, j, 6) = u(i, j)*tanr(j)
             atom(i, j, 7) = u(i, j)*tanr(j)
            ENDDO
           ENDDO
         END SELECT
!*
      END
!********************************************************
      SUBROUTINE vnlin(type, atom, u, v, h)
           use m_par
      use m_usr
      implicit none

!*     nonlinear terms for the v-equation                  u   u    j+1
!*       vv = v_ij*v_ij                                      v   v   j
!*       uu = ((u_i-1j+u_i-1j+1+u_ij+u_ij+1)*0.25)^2       u h u h  j
!*       uvx= (u_i-1j+u_i-1j+1+u_ij+u_ij+1)*(v_ij+1-v_ij-1)  v   v   j-1
!*       vvy = v_ij*(vij+1-v_ij-1)                              u    j-1
!*                                                           i  i+1 
!*     1:  vrV                                            i-1  i 
!*     2:  urVx /(r cos y)
!*     3:  Uvrx /(r cos y)   2 5            
!*     4:  vrVy/r           1 4 7
!*     5:  Vvry/r             3 6
!*     6:  urU  (tan y)/r
      integer type
      real    atom(n, m, np)
      real    u(0:n, 0:m+1), v(0:n+1, 0:m), h(0:n+1, 0:m+1)
      integer i, j
      real    costdxi(m), tdyi, tanr(m)
!*
      atom = 0.0
!*
      SELECT CASE(type)
       CASE(1)            ! vrV
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 4) = v(i, j)
         ENDDO
        ENDDO
       CASE(2)            ! urVx/(cos y)
        costdxi = 1.0/(cos(yv)*8*dx)
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 1) = (v(i+1, j)-v(i-1, j))*costdxi(j)
          atom(i, j, 2) = (v(i+1, j)-v(i-1, j))*costdxi(j)
          atom(i, j, 4) = (v(i+1, j)-v(i-1, j))*costdxi(j)
          atom(i, j, 5) = (v(i+1, j)-v(i-1, j))*costdxi(j)
         ENDDO
        ENDDO
       CASE(3)            ! Uvrx/(cos y)
        costdxi = 1.0/(cos(yv)*8*dx)
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 1) =-(u(i-1, j)+u(i-1, j+1)+u(i, j)+u(i, j+1))*costdxi(j)
          atom(i, j, 7) = (u(i-1, j)+u(i-1, j+1)+u(i, j)+u(i, j+1))*costdxi(j)
         ENDDO
        ENDDO
       CASE(4)            ! vrVy
        tdyi = 1.0/(2*dy)
        DO j = 1, m-1
         DO i = 1, n
          atom(i, j, 4) = tdyi*(v(i, j+1)-v(i, j-1))
         ENDDO
        ENDDO
       CASE(5)            ! Vvry
        tdyi = 1.0/(2*dy)
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 3) =-v(i, j)*tdyi
          atom(i, j, 5) = v(i, j)*tdyi
         ENDDO
        ENDDO
       CASE(6)            ! urU (tan y)
        DO j = 1, m
         DO i = 1, n
          tanr = tan(yv)/16.0
          atom(i, j, 1) =(u(i-1, j)+u(i-1, j+1)+u(i, j)+u(i, j+1))*tanr(j)
          atom(i, j, 2) =(u(i-1, j)+u(i-1, j+1)+u(i, j)+u(i, j+1))*tanr(j)
          atom(i, j, 4) =(u(i-1, j)+u(i-1, j+1)+u(i, j)+u(i, j+1))*tanr(j)
          atom(i, j, 5) =(u(i-1, j)+u(i-1, j+1)+u(i, j)+u(i, j+1))*tanr(j)
         ENDDO
        ENDDO
      END SELECT
!*
      END
!*****************************************************
      SUBROUTINE uderiv(type, atom)
      !implicit none
      use m_par
      use m_usr
       implicit none

!*     1:  u
!*     2:  u   /cos^2 y
!*     3:  uy  tan y
!*     4:  uxx/cos^2 y
!*     5:  uyy 
!*     6:  vy  2sin y/cos^2 y
!*     IMPORT/EXPORT
      integer type, i, j
      real    atom(n, m, np), fac(n, m)
      real    tdxi, tandyi(m), cosdx2i(m), rdy2i
      real    cos2i(m), sincos2i(m)
!*
      !print *, dx, dy
      atom = 0.0
      SELECT CASE(type)
       CASE(1)
        atom(:,:,4) = 1.0
       CASE(2)
        cos2i = 1.0/cos(y)**2
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 4) = cos2i(j)
         ENDDO
        ENDDO
       CASE(3)
        tandyi = tan(y)/(2*dy)
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 3) =-tandyi(j)
          atom(i, j, 5) = tandyi(j)
         ENDDO
        ENDDO
       CASE(4)
         
        cosdx2i = (1.0/(cos(y)*dx))**2
        !print *, "hey!!!!!"
        !print *, cos(y)*dx
        !print *, dx 
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 1) = cosdx2i(j)
          atom(i, j, 4) =-2*cosdx2i(j)
          atom(i, j, 7) = cosdx2i(j)
         ENDDO
        ENDDO
       CASE(5)
        rdy2i = (1.0/dy)**2
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 3) = rdy2i
          atom(i, j, 4) =-2*rdy2i
          atom(i, j, 5) = rdy2i
         ENDDO
        ENDDO
       CASE(6)
        sincos2i = sin(y)/(dx*cos(y)**2)
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 3) =-sincos2i(j)
          atom(i, j, 4) =-sincos2i(j)
          atom(i, j, 6) = sincos2i(j)
          atom(i, j, 7) = sincos2i(j)
         ENDDO
        ENDDO
      END SELECT
!*
      END
!*****************************************************
      SUBROUTINE vderiv(type, atom)
      !implicit none
      use m_par
      use m_usr
      implicit none
!*     1:  v
!*     2:  v   /cos^2 y
!*     3:  vy  tan y
!*     4:  vxx/cos^2 y
!*     5:  vyy
!*     6:  uy  2sin y/cos^2 y
!*     IMPORT/EXPORT
      integer type, i, j
      real    atom(n, m, np), fac(n, m)
      real    tdxi, tandyi(m), cosdx2i(m), dy2i
      real    cos2i(m), sincos2i(m)
!*
      atom = 0.0
      SELECT CASE(type)
       CASE(1)
        atom(:,:,4) = 1.0
       CASE(2)
        cos2i = 1.0/cos(yv)**2
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 4) = cos2i(j)
         ENDDO
        ENDDO
       CASE(3)
        tandyi = tan(yv)/(2*dy)
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 3) =-tandyi(j)
          atom(i, j, 5) = tandyi(j)
         ENDDO
        ENDDO
       CASE(4)
        cosdx2i = (1.0/(cos(yv)*dx))**2
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 1) = cosdx2i(j)
          atom(i, j, 4) =-2*cosdx2i(j)
          atom(i, j, 7) = cosdx2i(j)
         ENDDO
        ENDDO
       CASE(5)
        dy2i = (1.0/dy)**2
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 3) = dy2i
          atom(i, j, 4) =-2*dy2i
          atom(i, j, 5) = dy2i
         ENDDO
        ENDDO
       CASE(6)
        sincos2i = sin(y)/(dx*cos(y)**2)
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 1) =-sincos2i(j)
          atom(i, j, 2) =-sincos2i(j)
          atom(i, j, 4) = sincos2i(j)
          atom(i, j, 5) = sincos2i(j)
         ENDDO
        ENDDO
      END SELECT
!*
      END
!**********************************************************
      SUBROUTINE coriolis(type, atom)
     ! implicit none
      use m_par
      use m_usr
       implicit none

!*     1: fv in the u momentum equation     *   u   u   j+1
!*     2: fu in the v momentum equation     *     v   v    j
!*     1: v_ij+v_i+1j+v_ij-1 + v_i+1j-1 *   u h u h  j  
!*     2: u_ij+u_ij+1 + u_i-1j+u_i-1j+1 *     v   v    j-1
!*                                          *     i  i+1
!*     IMPORT/EXPORT                        *  i-1  i   
      integer type
      real    atom(n, m, np)
!*     LOCAL
      integer i, j
      real    cor(m)
!*
!*     
      IF (type .EQ. 1) THEN
       cor = sin(y)
       DO j = 1, m
        DO i = 1, n
         atom(i, j, 3) = cor(j)*0.25
         atom(i, j, 4) = cor(j)*0.25
         atom(i, j, 6) = cor(j)*0.25
         atom(i, j, 7) = cor(j)*0.25
        ENDDO
       ENDDO
      ELSE IF (type .EQ. 2) THEN
       cor = sin(yv)
       DO j = 1, m-1        !!!! potential problems !!!
        DO i = 1, n
         atom(i, j, 1) = cor(j)*0.25
         atom(i, j, 2) = cor(j)*0.25
         atom(i, j, 4) = cor(j)*0.25
         atom(i, j, 5) = cor(j)*0.25
        ENDDO
       ENDDO
      END IF
!*
      END
!**********************************************************
      SUBROUTINE gradh(type, atom)
      !implicit none
      use m_par
      use m_usr
       implicit none

!*     1: hx in the u momentum equation     *     h     j+1
!*     2: hy in the v momentum equation     *     v        j
!*     1: h_i+1j-h_ij                     *     h u h  j  
!*     2: h_ij+1 - h_ij                     *         
!*                                          *     i  i+1
!*     IMPORT/EXPORT                        *       i   
      integer type
      real    atom(n, m, np)
!*     LOCAL
      integer i, j
      real    cosdxi(m), dyi, grav
!*
      atom = 0.0
      IF (type .EQ. 1) THEN
       cosdxi = 1./(cos(y)*dx)
       DO j = 1, m
        DO i = 1, n
         atom(i, j, 4) =-cosdxi(j)
         atom(i, j, 7) = cosdxi(j)
        ENDDO
       ENDDO
      ELSE IF (type .EQ. 2) THEN
       dyi = 1./dy
       atom(:,:,4) =-dyi
       atom(:,:,5) = dyi
      END IF
!*
      END
!********************************************************
      SUBROUTINE hnlin(type, atom, u, v, h)
      !implicit none
      use m_par
      use m_usr
       implicit none

!*     nonlinear terms for the h-equation
!*                                      
!*     (              u_{i, j}~(h_{i+1, j}+h_{i, j}) - 
!*                               u_{i-1, j}~(h_{i, j}+h_{i-1, j}) )*costdxi_j
!*     (cos(yv_{i, j})~v_{i, j}~(h_{i, j+1}+h_{i, j}) - 
!*               cos(yv_{i, j-1})~v_{i, j-1}~(h_{i, j}+h_{i, j-1}) )*costdyi_j
!*
!*     1:  (urH)_x/cos y
!*     2:  (Uhr)_x/cos y
!*     3:  (vrH cos y)_y/cos y   2 5            
!*     4:  (Vhr cos y)_y/cos y   1 4 7
!*                                  3 6
      integer type
      real    atom(n, m, np)
      real    u(0:n, 0:m+1), v(0:n+1, 0:m), h(0:n+1, 0:m+1)
      integer i, j
      real    costdxi(m), costdyi(m), cov(0:m)
!*
      atom = 0.0
      cov(0) = cos(ygmin)
      !print *, 'YMIN', ygmin
      !print *, '***********************'
      DO j = 1, m
       cov(j) = cos(yv(j))
      ENDDO
!*
      SELECT CASE(type)
       CASE(1)            ! (urH)_x
        costdxi = 1.0/(2*cos(y)*dx)
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 1) =-(h(i-1, j)+h(i, j))*costdxi(j)
          atom(i, j, 4) = (h(i+1, j)+h(i, j))*costdxi(j)
         ENDDO
        ENDDO
       CASE(2)            ! (Uhr)_x
        costdxi = 1.0/(2*cos(y)*dx)
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 1) = (-u(i-1, j)         )*costdxi(j)
          atom(i, j, 4) = ( u(i, j)  -u(i-1, j))*costdxi(j)
          atom(i, j, 7) = ( u(i, j)           )*costdxi(j)
         ENDDO
        ENDDO
       CASE(3)            ! (vrH)_y
        costdyi = 1.0/(2*cos(y)*dy)
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 3) =-(h(i, j-1)+h(i, j))*costdyi(j)*cov(j-1)
          atom(i, j, 4) = (h(i, j+1)+h(i, j))*costdyi(j)*cov(j)
         ENDDO
        ENDDO
       CASE(4)            ! (Vhr)_y
        costdyi = 1.0/(2*cos(y)*dy)
        DO j = 1, m
         DO i = 1, n
          atom(i, j, 3)=(             -v(i, j-1)*cov(j-1))*costdyi(j)
          atom(i, j, 4)=(v(i, j)*cov(j)-v(i, j-1)*cov(j-1))*costdyi(j)
          atom(i, j, 5)=(v(i, j)*cov(j)                  )*costdyi(j)
         ENDDO
        ENDDO
      END SELECT
!*
      END
!**************************************************************************
      SUBROUTINE wind(type, atom, u, v, h)
      !implicit none
      use m_par
      use m_usr
       implicit none

!*     wind-forcing terms in u-and v-equations
!*       txh = taux/(h_ij+hi+1j)
!*       tyh = taux/(h_ij+hij+1)
!*
!*     1:  txh
!*     2:  tyh
!*     IMPORT/EXPORT
      integer type
      real    atom(n, m, np)
      real    u(0:n, 0:m+1), v(0:n+1, 0:m), h(0:n+1, 0:m+1)
      integer i, j, k
!*

      atom = 0.0
      SELECT CASE(type)
       CASE(1)            ! txh
         DO j = 1, m
          DO i = 1, n-1
           atom(i, j, 4) =(tx(i, j)+tx(i+1, j))/(h(i, j)+h(i+1, j))**2
           atom(i, j, 7) =(tx(i, j)+tx(i+1, j))/(h(i, j)+h(i+1, j))**2
      !    print *, 'i, j, tx=',i, j, tx(i, j), tx(i+1, j), (tx(n, j)+tx(1, j))
        !print *, 'sphere.f  h =',h(i, j), h(i+1, j), (h(i, j)+h(i+1, j))
       ! print *, 'sphere.f l_hth=',l_hth
  
          ENDDO
!x        atom(n, j, 1) = (tx(i, j)+tx(i-1, j))/(h(i, j)+h(i-1, j))**2
!x        atom(n, j, 4) = (tx(i, j)+tx(i-1, j))/(h(i, j)+h(i-1, j))**2
          atom(n, j, 4) =(tx(n, j)+tx(1, j))/(h(n, j)+h(1, j))**2
          atom(n, j, 7) =(tx(n, j)+tx(1, j))/(h(n, j)+h(1, j))**2
       ! print *, 'n, j, tx=', n, j, tx(n, j), tx(1, j), (tx(n, j)+tx(1, j))
       ! print *, 'h =',h(n, j), h(1, j), (h(n, j)+h(1, j))
        !print *, 'x-stress=',(tx(n, j)+tx(1, j))/(h(n, j)+h(1, j))
       
          ENDDO
               CASE(2)            ! tyh
         DO j = 1, m-1     !!!! potential problems !!!!
          DO i = 1, n
           atom(i, j, 4) =(ty(i, j)+ty(i, j+1))/(h(i, j)+h(i, j+1))**2
           atom(i, j, 5) =(ty(i, j)+ty(i, j+1))/(h(i, j)+h(i, j+1))**2
         !  print *, 'i, j, ty=',i, j, ty(i, j), ty(i+1, j), (ty(n, j)+ty(1, j))
       ! print *, 'h=',h(i, j), h(i+1, j), (h(i, j)+h(i+1, j))
        !print *, 'y-stress=',(ty(i, j)+ty(i+1, j))/(h(i, j)+h(i+1, j))
         ENDDO
         ENDDO
        
      END SELECT
     !print *, "tx = "
      !    print *,tx
          !print *, "cpnluh = "
          !print *, cpnluh
       !   print *, "h = ",h
        !  pause

100   format(2i8, 4es16.8)
      END
!**************************************************************************
      SUBROUTINE gradhb(hbx, hby)
      !implicit none
      use m_par
      use m_usr
       implicit none

!*     1: hbx in the u momentum equation: hb_i+1j-hb_ij
!*     2: hby in the v momentum equation: hb_ij+1 - hb_ij
!*     IMPORT/EXPORT
      real    hbx(n, m), hby(n, m)
!*     LOCAL
      integer i, j
      real    cosdxi(m), dyi
!*
      cosdxi = 1./(cos(y)*dx)
      dyi = 1./dy
      DO j = 1, m
       DO i = 1, n-1
        hbx(i, j) = cosdxi(j)*(hb(i+1, j)-hb(i, j))
       ENDDO
      ENDDO
      IF (periodic) THEN
       hbx(n, :) = cosdxi(:)*(hb(1, :)-hb(n, :))
      ELSE
       hbx(n, :) = 0.0
      ENDIF
      DO j = 1, m-1
       DO i = 1, n
        hby(i, j) = dyi*(hb(i, j+1)-hb(i, j))
       ENDDO
      ENDDO
      hby(:,m) = 0.0
!*
      END
!**************************************************************************
      SUBROUTINE layer(type, atom, u, v, h)
     ! implicit none
      use m_par
      use m_usr
      implicit none
 
!*     wind-forcing terms in u-and v-equations
!*       txh = taux/(h_ij+hi+1j)
!*       tyh = taux/(h_ij+hij+1)
!*
!*     1:  txh
!*     2:  tyh
!*     IMPORT/EXPORT
      integer type
      real    atom(n, m, np)
      real    u(0:n, 0:m+1), v(0:n+1, 0:m), h(0:n+1, 0:m+1)
      integer i, j, k, exp
      
      exp = int(par(15))
!*
      atom = 0.0
      SELECT CASE(type)
       CASE(1)            ! txh
         DO j = 1, m
          DO i = 1, n-1
           atom(i, j, 4) =(par(14)**exp)/(h(i, j)+h(i+1, j))**exp
           atom(i, j, 7) =(par(14)**exp)/(h(i, j)+h(i+1, j))**exp
          ENDDO
!x        atom(n, j, 1) = (tx(i, j)+tx(i-1, j))/(h(i, j)+h(i-1, j))**2
!x        atom(n, j, 4) = (tx(i, j)+tx(i-1, j))/(h(i, j)+h(i-1, j))**2
          atom(n, j, 4) =(par(14)**exp)/(h(n, j)+h(1, j))**exp
          atom(n, j, 7) =(par(14)**exp)/(h(n, j)+h(1, j))**exp
         ENDDO
       CASE(2)            ! tyh
         DO j = 1, m-1     !!!! potential problems !!!!
          DO i = 1, n
           atom(i, j, 4) =(par(14)**exp)/(h(i, j)+h(i, j+1))**exp
           atom(i, j, 5) =(par(14)**exp)/(h(i, j)+h(i, j+1))**exp
          ENDDO
         ENDDO
      END SELECT

100   format(2i8, 4es16.8)
      END
!**************************************************************************

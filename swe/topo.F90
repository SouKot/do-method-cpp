!* ============================================================================ *
!*
!* Multi-Layer Shallow Water Model
!*
!* Version 1.2-20/02/2008-ATvS
!*
!* ============================================================================ *
      SUBROUTINE removeland
        !implicit none
        use m_par
        use m_usr
        implicit none

!*     IMPORT/EXPORT
!*     LOCAL
        integer i, j, row(4+nl), rw, idx, is, js, ks, p, k, k2
!*     FUNCTION
        logical new
!*
!*     numbering u, v, h (from i, j)
!*     at every land point (central position of grid box)
!*     remove u_ij, u_i-1j, v_ij, v_ij-1, h_ij
        idx = 0
        is = 3
        js = 3*N
        ks = 3*n*m
        do k2 = 1, nl
        do j = 1, m
        do i = 1, n

        if (lnd(i, j).eq.1) then
!             Compute the position of the data for cell (i, j) in the
!             1D-array.
          rw  = 3*(n*m*(k2-1)+n*(j-1) + i-1)
          row(1) = rw-is   +1    !  u_(i-1, j) 
          row(2) = rw      +1    !  u_(i, j)
          row(3) = rw   -js+2    !  v_(i, j-1)
          row(4) = rw      +2    !  v_(i, j)
          row(5) = rw      +3

          do k = 1, 5
          if (new(dir, idx, row(k))) then
            if ((row(k).ge.3*n*m*(k2-1)+1).and. (row(k).le.3*n*m*k2)) then
              idx = idx+1
              dir(idx) = row(k)
            endif
          endif
          enddo

        endif

        enddo
        enddo
        enddo
        ndir = idx
        write(99, *) 'number of dirichl pnts: ',ndir
!*
        END
!********************************************************
        FUNCTION new(arr, nel, val)
!*     checks if val is a new value in de array of length nel
          implicit none
!*     IMPORT/EXPORT
          integer nel, arr(nel), val
!*     LOCAL
          integer i
!*     FUNCTION
          logical new
!*
          new = .true.
          DO i = 1, nel
          IF (val .EQ. arr(i)) new = .false.
          ENDDO
!*
          END
!********************************************************
          subroutine land2mat(F)
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
            if (jcoa(v).eq.row) coa(v) = coa(v)+(1.0-par(10))
            enddo
            if (mod(row, 3).ne.0) then
              F(row) = par(10)*F(row)
            else
              do k = 1, nl-1
              if ((row .gt. 3*(k-1)*n*m).and.(row .le. 3*k*n*m)) then
                F(row) = par(10)*F(row)+(1.0-par(10))*(l_hth(k)-l_hth(k+1))
              endif
              enddo  
              if (row .gt. 3*(nl-1)*n*m) then
                F(row) = par(10)*F(row)+(1.0-par(10))* l_hth(nl)
              endif 
            endif
            enddo
            END
!********************************************************
            SUBROUTINE topography
!     implicit none
              use m_par
              use m_usr
              implicit none

!*     IMPORT/EXPORT
!*     LOCAL
              integer i, j, ifile, ice, icw, jcs, jcn
              real    pi, dum, x0, y0, x1, y1, alp
              integer*2 land(n, m)
! SP = specific
              real    yp(5), xp(5), f1, f2, ap, bp
              integer ns

              logical curve
!*
              curve= .true.
              lnd = 0
              IF (curve) THEN
                pi = 4.0*atan(1.0)
! SP = specific
                f1 = pi/180.
                xp(1) = f1*122.
                xp(2) = f1*134.
                xp(2) = f1*134.
                xp(3) = f1*141.
                xp(4) = f1*143.
                xp(5) = f1*171.
                yp(1) = f1*22.
                yp(2) = f1*33.
                yp(3) = f1*35.
                yp(4) = f1*42.
                yp(5) = f1*55.
                DO j = 1, m
                DO ns = 1, 4
                IF ((yp(ns).le.y(j)).and.(y(j).le.yp(ns+1))) THEN
                  ap  = (yp(ns+1)-yp(ns))/(xp(ns+1) - xp(ns))
                  bp  = yp(ns) - ap*xp(ns)
                  f2 = (y(j)- bp)/ap
                  DO i = 1, n
                  IF (x(i).lt.f2) THEN
                    lnd(i, j) = 1
                  ENDIF
                  ENDDO
                ENDIF
                ENDDO
                ENDDO
                hb  = 0.0
                hb1 = 0.0
                davg = 0.0
!*       lnd = 0
              ELSE
                call topofit
              ENDIF
              call removeland
!     do j = 1, m
!        write(90, 101) (lnd(i, j), i = 1, n)
!     enddo
!101  format(1x, 160(i1, 1x))
!     stop
!*
              END
!****************************************************************
              SUBROUTINE topofit
     ! implicit none
                use m_par
                use m_usr
                implicit none

!*     CONSTANT
                !integer, parameter:: nd = 6*360, md = 6*180
     ! integer  lwrk, liwrk 
                !*     LOCAL
                integer(2):: d10(nd, md)
                real    :: rd(md, nd), xd(nd), yd(md), dxd, dyd, wrk2(lwrk)
                integer:: ifail, px, py, iwrk(liwrk)
                real    :: lambda(nd+4), mu(md+4), cspl(nd*md), wrk((nd+6)*(md+6))
                real    :: dumd(n*m)
                real    :: xb(n)

                real    :: depth(n, m), pi, dav, hbs(n, m), dmin, dmax
                real    :: davnd, davd, tpi
                integer:: i, j, io, sfend, sfstart, l, tel
                logical:: isopen
!*
!      open(20, file='topo10.bin', form='unformatted',
!     &     access='sequential', status='old', recordtype='stream')
! TODO: Check this change. RECORDTYPE is non-standard, thus probably
! not needed!
                !     lwrk = 4*n+nd+4
                !     liwrk = n+nd

                open(20, file='topo10.bin', form='unformatted', access='sequential', status='old')
                inquire( unit = 20, OPENED = isopen )
                if ( .not. isopen ) then
                  write(*,*) 'Failed to open topo10.bin'
                  stop 
                endif      
                read(20) d10
                close(20)

                DO j = 1, md
                DO i = 1, nd
                rd(md+1-j, i) = real(d10(i, j))
                ENDDO
                ENDDO
!*
                pi  = 4.0*atan(1.0)
                tpi = 2*pi
                dxd = 2.0*pi/nd
                dyd = pi/md
                DO i = 1, nd
                xd(i) = (real(i)-0.5)*dxd
                ENDDO
                DO j = 1, md
                yd(j) = (real(j)-0.5)*dyd-0.5*pi
                ENDDO
!*
!TODO: Remove NAG routines from topofit.
                ifail = 0
                write (*,*) 'Routine topofit: NAG nog verwijderen!'
                stop
!      call e01daf(nd, md, xd, yd, rd, px, py, lambda, mu, cspl, wrk, ifail)
                DO i = 1, n
                IF (x(i).gt.tpi) THEN
                  xb(i) = x(i) - tpi
                ELSE
                  xb(i) = x(i)
                END IF 
                ENDDO      
!      call e02dff(n, m, px, py, xb, y, lambda, mu, cspl, dumd, wrk2, lwrk, 
!     &            iwrk, liwrk, ifail)
                DO i = 1, n
                DO j = 1, m
                depth(i, j) = dumd(m*(i-1)+j)
                ENDDO
                ENDDO
!*
                call depth2land(depth)
                tel = 0
                DO j = 1, m
                DO i = 1, n
                IF (lnd(i, j).GT.0) then
                  lnd(i, j) = 1
                  tel = tel+1
                ENDIF 
                ENDDO
                ENDDO
                write(99, *) 'no land points: ',tel
                nopo = n*m - tel 
!*
!* Determine the average depth of unsmoothed topography 
!* Use only ocean points in bathymetry
!* 
                dav = 0.0
                DO i = 1, n
                DO j = 1, m
                dav = dav+depth(i, j)*(1. - lnd(i, j))
                ENDDO
                ENDDO
                dav = dav/nopo      
                dmax = maxval(depth)
                dmin = minval(depth)
                write(99, *) 'min, max, average depth:', dmin, dmax, dav
!* H = hdim is the characteristic vertical scale       
!* real topography 
!      hb =  (hdim+depth)/hdim
!*     hb =  (hdim+dav)/hdim
!* test
!*     DO i = 1, n
!*       DO j = 1, m
!*          hb(i, j) = -1.0 + (i-1)*0.5/(n-1)
!*          hb(i, j) = -1.0 + (i-1)*(j-1)*0.5/((n-1)*(m-1))
!*       ENDDO
!*     ENDDO
!* If ism > 0, smooth the topography 
                if (ism .ge. 1) then
                  do l = 1, ism
                  hbs = hb
                  do i = 2, n-1
                  do j = 2, m-1
                  hbs(i, j) = (0.5*(hb(i+1, j) + hb(i-1, j) +  hb(i, j+1) + hb(i, j-1)) + 2*hb(i, j))/4
                  enddo
                  enddo
                  hb = hbs
                  enddo
                  do j = 1, m
                  do i = 1, n
!*    minimum depth at 50 m
! TODO: hdim = 4000, so 200m is the limit?
                  if (hbs(i, j).gt.0.95) hbs(i, j) = 0.95
                  if (lnd(i, j).gt.0)    hbs(i, j) = 0.0
                  enddo
                  enddo
!* determine the average depth of the smoothed topography
                  dav = 0.0
                  do i = 1, n
                  do j = 1, m
                  dav = dav+hbs(i, j)
                  enddo
                  enddo
                  davnd = dav/nopo
                  davd =  hdim*davnd-hdim         
                  dmax = hdim*maxval(hbs) - hdim
                  dmin = hdim*minval(hbs) - hdim
                  write(99, *) 'min, max, average depth [smoothed]:', dmin, dmax, davd, davnd
                  do i = 1, n
                  do j = 1, m
                  davg(i, j) = davnd * (1-lnd(i, j))
                  enddo
                  enddo
                  hb1 = hbs
                  hb = hbs
                else
                  do i = 1, n
                  do j = 1, m
!*    minimum depth at 50 m
                  if (hb(i, j).gt.0.95) hb(i, j) = 0.95
                  if (lnd(i, j).gt.0)   hb(i, j) = 0.0
                  enddo
                  enddo
! TODO: strange order of assignements!
                  hb1 = hb     ! unsmoothed topography 
                  hb = hb1
                endif 

                write(99, *) maxval(hb1), minval(hb1)
!*
                END
!******************************************************************
                SUBROUTINE depth2land(depth)
                  !implicit none
                  use m_par
                  use m_usr
                  implicit none

!*     IMPORT/EXPORT
                  real    depth(n, m)
!*     LOCAL
                  integer i, j, iter, itmax
                  real    pi, x0, x1, y0, y1
!*
                  pi = 4.0*atan(1.0)
                  lnd = 1
                  DO j = 1, m
                  DO i = 1, n
                  IF (depth(i, j).LT.0) THEN
                    lnd(i, j) = 2
!*       'missing' middle america
                    x0 = 275*pi/180.
                    x1 = 282*pi/180.
                    y0 =   8*pi/180.
                    y1 =  10*pi/180.
                    IF ( (x(i).ge.x0).AND.(x(i).le.x1) .AND. (y(j).ge.y0).AND.(y(j).le.y1) ) lnd(i, j) = 1
!*       'close the mediterran
                    x0 =  -5*pi/180.
                    x1 =  10*pi/180.
                    y0 =  35*pi/180.
                    y1 =  45*pi/180.
                    IF ( (x(i).ge.x0).AND.(x(i).le.x1) .AND. (y(j).ge.y0).AND.(y(j).le.y1) ) lnd(i, j) = 1
                  ENDIF
!*       remove south orkney islands
                  x0 = 295*pi/180.
                  x1 = 305*pi/180.
                  y0 = -65*pi/180.
                  y1 = -60*pi/180.
                  IF ( (x(i).ge.x0).AND.(x(i).le.x1) .AND. (y(j).ge.y0).AND.(y(j).le.y1) ) lnd(i, j) = 2
!*       remove kerguelen islands
                  x0 =  68*pi/180.
                  x1 =  72*pi/180.
                  y0 = -50*pi/180.
                  y1 = -45*pi/180.
                  IF ( (x(i).ge.x0).AND.(x(i).le.x1) .AND. (y(j).ge.y0).AND.(y(j).le.y1) ) lnd(i, j) = 2
!*       remove falkland islands
                  x0 = 297*pi/180.
                  x1 = 306*pi/180.
                  y0 = -55*pi/180.
                  y1 = -50*pi/180.
                  IF ( (x(i).ge.x0).AND.(x(i).le.x1) .AND. (y(j).ge.y0).AND.(y(j).le.y1) ) lnd(i, j) = 2
                  ENDDO
                  ENDDO
!*
                  lnd(n/2, m/2) = 3
                  itmax = 400
                  DO iter = 1, itmax
                  DO j = 1, m
                  DO i = 1, n
                  IF (lnd(i, j).EQ.3) THEN
                    IF (i .LT. n) THEN
                      IF (lnd(i+1, j).EQ.2) lnd(i+1, j) = 3
                    ELSE IF (periodic) THEN
                      IF (lnd(1, j).EQ.2)   lnd(1, j) = 3
                    ENDIF
                    IF (i .GT. 1) THEN
                      IF (lnd(i-1, j).EQ.2) lnd(i-1, j) = 3
                    ELSE IF (periodic) THEN
                      IF (lnd(n, j).EQ.2)   lnd(n, j) = 3
                    ENDIF
                    IF (j .LT. m) THEN
                      IF (lnd(i, j+1).EQ.2) lnd(i, j+1) = 3
                    ENDIF
                    IF (j .GT. 1) THEN
                      IF (lnd(i, j-1).EQ.2) lnd(i, j-1) = 3
                    ENDIF
                    lnd(i, j) = 0
                  ENDIF
                  ENDDO
                  ENDDO
                  ENDDO
!*
                  END

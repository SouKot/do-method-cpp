!* ============================================================================ *
!*
!* Multi-Layer Shallow Water Model
!*
!* Version 1.2-20/02/2008-ATvS
!*
!* ============================================================================ *
      subroutine init

        use m_par
        use m_usr
        implicit none

        real    pi
        integer i, j

        integer info, iwa(6*ndim), liwa  ! mixing
        integer minrow, maxrow
        real dnsm

        pi = 4.0*atan(1.0)
        !print *, 'pmin, pmax, thmin, thmax :',pgmin, pgmax, tghmin, tghmax
        !pause
        ! conversion deg -> rad
        xmin = pmin*pi/180.0
        xmax = pmax*pi/180.0
        ymin = thmin*pi/180.0
        ymax = thmax*pi/180.0
        xgmin = pgmin*pi/180.0
        xgmax = pgmax*pi/180.0
        ygmin = tghmin*pi/180.0
        ygmax = tghmax*pi/180.0

!TODO: Why not test on abs( pmax-pmin-360. )? Less roundoff problems!
        ! test for global grid
        if (abs(xgmax-xgmin-360.*pi/180).lt.1.0e-2) periodic = .true.

        write (99, *) periodic
        ! test if xmax < xmin
        if (xmax .lt. xmin) then
          xmax = 2*pi+xmax
        endif 

        ! define grid
        dx = (xmax-xmin)/n
        dy = (ymax-ymin)/m
        do i = 1, n
        x(i) = (real(i)-0.5)*dx+xmin
        enddo
        do j = 1, m
        y(j) = (real(j)-0.5)*dy+ymin
        yv(j)= (real(j)    )*dy+ymin
        enddo
! read winds from file and interpolate on grid
        call windforcing( windfromdata, x, n, y, m, tx, ty, tmax )

!OLD:       call windfit( x, n, y, m, taux, tauy, tmax )
! read bottom topography from file and interpolate on grid
        call topography

        call renumber

        ! mixing
        if (mix_flag) then
          call mix_el(irow, icol, idim)

          liwa = 6*ndim
          call dsm(ndim, ndim, idim, irow, icol, ngrp, maxgrp, mingrp, info, ipntr, jpntr, iwa, liwa)
          if (info .le. 0) then
            write (99, *) 'Error in DSM'
            stop
          endif  
          maxrow = 0
          minrow = ndim
          do i = 1, ndim
          maxrow = max(maxrow, ipntr(i+1)-ipntr(i))
          minrow = min(minrow, ipntr(i+1)-ipntr(i))
          enddo
          dnsm = real(idim)/(real(ndim)**2)
          write (99, '(3(i10, 3x), es16.8, 4x, 5(i10, 5x))')&
            ndim, idim, info, dnsm, minrow, maxrow, mingrp, maxgrp
        endif
        ! mixing

        end

!*******************************************************************************
!
! Test the constraints on variables and constants.
!
!*******************************************************************************

        subroutine testconstconstraints

          use m_par
          use m_usr
          implicit none

          ! If we have a n.5-layer model, then bfric = ifric and
          ! tot = .false.
          if ( mmode .eq. 1 ) then 
            if ( bfric .ne. ifric ) then
              print *, 'In a n.5-layer model, we require bfric = ifric.'
              stop
            endif   
            if ( topo ) then      
              print *, 'In a n.5-layer model, topo should be .false.'
              stop    
            endif      
          endif

          end


!*******************************************************************************
!
! Variable initialisations.
!
!*******************************************************************************

!
! subroutine iniparconst
!
! Initialise the vector par from the constants defined in param.com
!
          subroutine iniparconst


            use m_par
            use m_usr    ! Needed: par, tmax.
implicit none
            par = 0.0

            ! parameters: 
            ! r : mean radius of earth (6.37 10^6 m)
            ! om: rotation rate of earth (7.292 10^-5 s-1)
            ! u : velocity scale (0.1 m/s)
            ! g : acceleration due to  gravity (9.8  m/s)
            ! a : laplace friction (control)
            ! r : bottom friction (0 s-1)
            ! h : layer depth (1000 m)
            ! t0: wind-stress 0.1 n/m (pa) * tmax uit fit 
            ! t = r/u t, x = r x, u = u u
            ! 1: wind stress sigma (t0/(2 om h u rho) = 1.4 10^-2
            ! 2: bottom friction (r/2 om) = 0.0
            ! 3: eps (u/2 om r) = 1.02 10 ^-4
            ! 4: ek  (a/2 om r^2) = 1.5 10^-10*a 
            ! 6: f (gh/u^2) = 1e6
            ! 7: h (equilibrium value) 
            ! 9: sine/data windstress (0/1)
            ! 10: land/no land (0/1)
            ! 12: topography/no topography (1/0)  
            ! 13: smoothed/non-smoothed topo (0/1)
            ! tmax is the windstress scaling factor determined in windfit.
            ! All other constants from which the entries of par are computed 
            ! are defined in param.com.
            tmax = 1.
            par(1) = taudim*tmax/(2*omegadim*hdim*udim*rhodim)    ! windstress coef.
            !par(1) = taudim/rhodim    ! windstress coef.
            par(2) = bfric                                        ! bottom friction
            par(3) = udim/(2*omegadim*r0dim)                      ! rossby number
            par(4) = Ahdim/(2.0*omegadim*r0dim**2)                ! ekman number
            par(5) = r0dim/udim                 ! time scale
            par(6) = hdim*gdim/(udim*udim)  ! xoxgdim*hdim/(udim*udim)  ! froud number
            par(7) = 1.0e+00                                      ! eq. value h
            par(8) = ifric                                        ! interface friction
            if ( windfromdata ) then  ! sine/data windstress [0/1]
              par(9) = 1.0e+00 
            else
              par(9) = 0.0e+00
            endif
            if (cont) then            ! land/no land [0/1]
              par(10) = 0.0e+00
            else
              par(10) = 1.0e+00
            endif
            par(11)= 1.0e+00          ! unused
            if (topo) then            ! topography/no topography [1/0]
              par(12) = 1.0e+00
            else
              par(12) = 0.0e+00
            endif
            if (smoothed) then        ! smoothed/non smoothed topo. [0/1]
              par(13) = 0.0e+00
            else
              par(13) = 1.0e+00
            endif
            if (potential) then
              par(14)=h0
              par(15)=hn
            else
              par(14)=0.0
              par(15)=0.0
            endif  

            end

!
! subroutine initrivstate( un )
!
! Initialise the state with a `trivial` initialisation: zero velocity 
! and layer heights from a constant.
!
            subroutine initrivstate( un )


              use m_par
              use m_usr
              implicit none
! Arguments
              real, dimension(ndim), intent(out):: un
! Local variables
              real, dimension(0:n, 0:m+1, 1:nl)   :: u
              real, dimension(0:n+1, 0:m, 1:nl)   :: v
              real, dimension(0:n+1, 0:m+1, 1:nl):: hh
              integer:: i

              l_hth(1)= 1.0
              !l_hth(2)=0.825 

              u = 0.0
              v = 0.0
              if ( nl .gt. 1 ) then
                do i = 1, nl-1
                hh(:,:,i) = l_hth(i) - l_hth(i+1)

                enddo
              endif
              hh(:,:,nl) = l_hth(nl)
              hh(1:n, 1:m, nl) = l_hth(nl) - par(12)*hb
              !print *,'h = ',hh
              !print *, 'l_hth = ',l_hth
              call solu(un, u, v, hh)

              end

!
! subroutine inicoef
!
! Initialises the common block variable l_coef(user.com), the layer
! coefficients.
!
! Uses: par. Note that in some cases, values in par are changed!
!
              subroutine inilcoef


                use m_par
                use m_usr
                implicit none
! Local variables
                integer:: i, j


                select case(mmode)

                case(0)  ! full model with topography
                  l_coef(1, 1:nl) = par(3)*par(6)
                  if (topo) l_coef(1, nl+1) = par(3)*par(6)
                  if ( nl .gt. 1 ) then
                    do i = 2, nl
                    l_coef(i, :) = l_coef(i-1, :)
                    do j = i, nl
                    l_coef(i, j) = l_coef(i, j)&              
                      + par(3)*par(6)*(l_den(i)-l_den(i-1))/rhodim
                    enddo
                    if (topo) l_coef(i, nl+1) = l_coef(i, nl+1)&
                      +par(3)*par(6)*(l_den(i)-l_den(i-1))/rhodim
                    enddo
                  endif

                case(1)  ! n.5 layer model (n+1-th layer in rest)
                  ! note: par(12) = 0.0 (no bottom topography)
                  ! note: par(2)  = par(8)
                  if ( par(12) .ne. 0.0 ) then
                    print *, '!!!!!'
                    print *, 'WARNING: par(12) (topo) was not zero, while ',&
                      'this is assumed to be the case for a n.5-layer model.',&
                      'Its value will be changed to zero.'
                    print *, '!!!!!'
                  endif
                  par(12) = 0.0
                  if ( par(2) .ne. par(8) ) then
                   ! print *, '!!!!!'
                   ! print *, 'For a n.5-layer model, par(2) (bfric) and ',&
                   !   'par(8) (ifric) should be the same. par(2) will be ',&
                   !   'changed!'
                   ! print *, '!!!!!'
                  endif
                  par(2) = par(8)
                  l_coef(nl, 1:nl)=par(3)*par(6)*(l_den(nl+1)-l_den(nl))/rhodim
                  if ( nl .gt. 1 ) then
                    do i = nl-1, 1, -1
                    l_coef(i, 1:nl) = l_coef(i+1, 1:nl)
                    do j = 1, i
                    l_coef(i, j) = l_coef(i, j)&
                      + par(3)*par(6)*(l_den(i+1)-l_den(i))/rhodim
                    enddo
                    enddo
                  endif

                end select    

                end  ! subroutine inilcoef

!
! subroutine stpnt( un )
!
! Initialise the state and parameter vector from constants.
! Initialise l_coef.
!
                subroutine stpnt(un)

                  use m_par
                  use m_usr
implicit none

! Arguments
                  real, dimension(ndim), intent(out):: un
! Local variables
                  integer:: i, j       

                  ! Initialise the par vector using the constants in param.com
                  call iniparconst

                  ! layer coefficients for discretization
                  call inilcoef
!l_hth(1)=1.0; 
 !    l_hth(2)=0.825; 
 ! Initial solution.
 ! This is done after the call to inilcoef since that routine 
 ! might change par(12) which is used in initrivstate!
 call initrivstate( un )

!print *, "Hello there!!!!!"
 ! Write some output to the screen.
 do i = 1, nl
 write (99, 100) i
 write (99, 102) l_den(i)
 !PRINT *, 'help density', l_den
 write (99, 103) l_hth(i)
 !PRINT *, 'help height', l_hth
 !STOP

 if ( mmode .eq. 1 ) then
   if ( i .eq. nl)  then
     write (99, 104) gdim*(rhodim-l_den(i))/rhodim
   else
     write (99, 104) gdim*(l_den(i+1)-l_den(i))/rhodim
   endif  
   endif  
     write (99, *) ' - Coefficients pressure term:'
     do j = 1, nl
     write (99, 101) j, l_coef(i, j)
     enddo
     write (99, 105) l_coef(i, nl+1)
     write (99, *)
     enddo
!print *,'density',l_den
!print *, 'l_hth',l_hth        
100     format('Layer:', i2)
101     format('   ',i2, ':',es12.4)
102     format('  - Density:', es12.4)
103     format('  - Interface eq. heigth:', es12.4)
104     format('  - Reduced gravity:', es12.4)
105     format('    b:',es12.4)
!stop
     end  ! subroutine stpnt



!*******************************************************************************
!
! Data conversion functions: Conversion of the state vector.
!
!*******************************************************************************

     subroutine solu(un, u, v, h)
!*     go from u, v, h to un

     use m_par
     use m_usr
!*     import/export
     implicit none
     real    un(ndim)
     real    u(0:n, 0:m+1, 1:nl), v(0:n+1, 0:m, 1:nl), h(0:n+1, 0:m+1, 1:nl)
!*     local
     integer i, j, k, row
!*
     do k = 1, nl
     do j = 1, m
     do i = 1, n
     row =  3*(n*m*(k-1)+n*(j-1) + i-1)
     un(row+1) = u(i, j, k)
     un(row+2) = v(i, j, k)
     un(row+3) = h(i, j, k)
     enddo
     enddo
     enddo

     end

!*******************************************************************************

     subroutine usol(un, u, v, h)


     use m_par
     use m_usr
 implicit none
     real, dimension(ndim), intent(in):: un
     real, dimension(0:n, 0:m+1, 1:nl), intent(out):: u
     real, dimension(0:n+1, 0:m,  1:nl), intent(out):: v
     real, dimension(0:n+1, 0:m+1, 1:nl), intent(out):: h
!      real    un(ndim)
!      real    u(0:n, 0:m+1, 1:nl), v(0:n+1, 0:m, 1:nl), h(0:n+1, 0:m+1, 1:nl)

     integer i, j, k, row

     u = 0.0
     v = 0.0
     !h = 0.0
     do k = 1, nl
     do j = 1, m
     do i = 1, n
     row =  3*(n*m*(k-1)+n*(j-1)+i-1)
     u(i, j, k) = un(row+1)
     v(i, j, k) = un(row+2)
     h(i, j, k) = un(row+3)
     !print *, i, j, k, h(i, j, k)
     enddo
     enddo
     enddo
!Print *, 'h in solu =',h
! Boundary conditions psi-direction.
     do j = 1, m
     if (periodic) then
     u(0, j, :)  = u(n, j, :)
     v(0, j, :)  = v(1, j, :)
     v(n+1, j, :)= v(n, j, :)
     h(0, j, :)  = h(n, j, :)
     h(n+1, j, :)= h(1, j, :)
   else
     u(0, j, :)  = 0.0
     v(0, j, :)  =-v(1, j, :)
     v(n+1, j, :)=-v(n, j, :)
     h(0, j, :)  = 0.0
     h(n+1, j, :)= 0.0
   endif
     enddo
     if (pper) then
     do j = 1, mper
     u(0, j, :)  = u(n, j, :)
     v(0, j, :)  = v(1, j, :)
     v(n+1, j, :)= v(n, j, :)
     h(0, j, :)  = h(n, j, :)
     h(n+1, j, :)= h(1, j, :)
     enddo
     do j = mper+1, m
     u(0, j, :)  = 0.0
     v(0, j, :)  =-v(1, j, :)
     v(n+1, j, :)=-v(n, j, :)
     h(0, j, :)  = 0.0
     h(n+1, j, :)= 0.0
     enddo
   endif

! Boundary conditions theta-direction.
     do i = 1, n
     u(i, 0, :)  =-u(i, 1, :)
     u(i, m+1, :)=-u(i, m, :)
     v(i, 0, :)  = 0.0
     h(i, 0, :)  = 0.0
     h(i, m+1, :)= 0.0
     enddo


     end

!*******************************************************************************

     subroutine shuflvec(rl)


     use m_par
     use m_usr
implicit none
     real rl(ndim), rdum(ndim)
     integer i

     do i = 1, ndim
     rdum(i)=rl(i)
     enddo
     do i = 1, ndim
     rl(i)=rdum(num(i))
     enddo

     end

!*******************************************************************************

     subroutine reshuflvec(rl)

     use m_par
     use m_usr
implicit none

     real rl(ndim), rdum(ndim)
     integer i

     do i = 1, ndim
     rdum(i)=rl(i)
     enddo
     do i = 1, ndim
     rl(i)=rdum(inv(i))
     enddo

     end

!*******************************************************************************
!
! Data conversion: Jacobian matrix.
!
!*******************************************************************************

     subroutine shuflmat

     use m_par
     use m_usr
     use m_mat
     implicit none

     integer i, tel, v, vv, jdum, k, j
     real dum
     real, allocatable, dimension(:):: storco
     integer, allocatable, dimension(:):: storjco
     integer storbeg(ndim+1), nnz

     nnz = begA(ndim+1)-1
     allocate(storco(nnz))
     allocate(storjco(nnz))
     storco = 0; storjco = 0

     storco = coA(1:nnz)
     storjco = jcoA(1:nnz)
     storbeg = begA

     ! interchange rows
     tel = 1
     do i = 1, ndim
     begA(i)=tel
     do v = storbeg(num(i)), storbeg(num(i)+1)-1
     coA(tel)=storco(v)
     jcoA(tel)=inv(storjco(v))
     tel = tel+1
     enddo     
     enddo
     begA(ndim+1)=tel

     ! Sorting
     do i = 1, ndim
     do v = begA(i), begA(i+1)-2
     k = v
     do vv = k+1, begA(i+1)-1
     if (jcoA(vv).lt.jcoA(k)) then
     k = vv
   endif
     dum     = coA(k)
     jdum    = jcoA(k)
     coA(k)  = coA(v)
     jcoA(k) = jcoA(v)
     coA(v)  = dum
     jcoA(v) = jdum
     enddo
     enddo
     enddo



     end



!*******************************************************************************
!
! Data conversion: Auxiliary routine.
!
!*******************************************************************************


     subroutine renumber

     use m_par
     use m_usr
     implicit none

     integer:: k, i, j, l
     integer:: row
     !print *, 'value of m n and nl is : ', m, n, nl

     l = 0
     do j = 1, m
     do i = 1, n
     do k = 1, nl
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     num(l+1) = row+1
     num(l+2) = row+2
     num(l+3) = row+3
     inv(row+1) = l+1
     inv(row+2) = l+2
     inv(row+3) = l+3
     l = l+3
     enddo  
     enddo
     enddo
     !stop
     end



!*******************************************************************************
!
!
!
!*******************************************************************************


     subroutine assemble
!*     assemble the global matrix a from the local matrices
     use m_par
     use m_usr
     use m_mat
!*     locl 
     implicit none

     integer i, j, k, l, v, row
!*
!*     2 5  
!*     1 4 7
!*       3 6
!*
     cpcoA = 0.0
     ccoA = 0.0
     coA = 0.0
     call fillcola
     do k = 1, nl
     do j = 1, m
     do i = 1, n
     row = 3*8*3*nl*(n*m*(k-1)+n*(j-1)+i-1)
     do l = 1, nl
     v = row+24*(l-1)
     coA(v+1) = aluu(i, j, 1, k, l)
     coA(v+2) = aluv(i, j, 1, k, l)
     coA(v+3) = aluh(i, j, 1, k, l)
     cpcoA(v+3)=cpnluh(i, j, 1, k, l)
     coA(v+4) = aluu(i, j, 2, k, l)
     coA(v+5) = aluv(i, j, 2, k, l)
     coA(v+6) = aluh(i, j, 2, k, l)
     cpcoA(v+6)=cpnluh(i, j, 2, k, l)
     coA(v+7) = aluu(i, j, 3, k, l)
     coA(v+8) = aluv(i, j, 3, k, l)
     coA(v+9) = aluh(i, j, 3, k, l)
     cpcoA(v+9)=cpnluh(i, j, 3, k, l)
     coA(v+10)= aluu(i, j, 4, k, l)
     coA(v+11)= aluv(i, j, 4, k, l)
     coA(v+12)= aluh(i, j, 4, k, l)
     cpcoA(v+12)=cpnluh(i, j, 4, k, l)
     coA(v+13)= aluu(i, j, 5, k, l)
     coA(v+14)= aluv(i, j, 5, k, l)
     coA(v+15)= aluh(i, j, 5, k, l)
     cpcoA(v+15)=cpnluh(i, j, 5, k, l)
     coA(v+16)= aluu(i, j, 6, k, l)
     coA(v+17)= aluv(i, j, 6, k, l)
     coA(v+18)= aluh(i, j, 6, k, l)
     cpcoA(v+18)=cpnluh(i, j, 6, k, l)
     coA(v+19)= aluu(i, j, 7, k, l)
     coA(v+20)= aluv(i, j, 7, k, l)
     coA(v+21)= aluh(i, j, 7, k, l)
     cpcoA(v+21)=cpnluh(i, j, 7, k, l)


     ccoA(v+1) = lluu(i, j, 1, k, l)
     ccoA(v+2) = lluv(i, j, 1, k, l)
     ccoA(v+3) = lluh(i, j, 1, k, l)

     ccoA(v+4) = lluu(i, j, 2, k, l)
     ccoA(v+5) = lluv(i, j, 2, k, l)
     ccoA(v+6) = lluh(i, j, 2, k, l)

     ccoA(v+7) = lluu(i, j, 3, k, l)
     ccoA(v+8) = lluv(i, j, 3, k, l)
     ccoA(v+9) = lluh(i, j, 3, k, l)

     ccoA(v+10)= lluu(i, j, 4, k, l)
     ccoA(v+11)= lluv(i, j, 4, k, l)
     ccoA(v+12)= lluh(i, j, 4, k, l)

     ccoA(v+13)= lluu(i, j, 5, k, l)
     ccoA(v+14)= lluv(i, j, 5, k, l)
     ccoA(v+15)= lluh(i, j, 5, k, l)

     ccoA(v+16)= lluu(i, j, 6, k, l)
     ccoA(v+17)= lluv(i, j, 6, k, l)
     ccoA(v+18)= lluh(i, j, 6, k, l)

     ccoA(v+19)= lluu(i, j, 7, k, l)
     ccoA(v+20)= lluv(i, j, 7, k, l)
     ccoA(v+21)= lluh(i, j, 7, k, l)

     enddo

     do l = 1, nl
     v = row+24*nl+24*(l-1)
     coA(v+1 )= alvu(i, j, 1, k, l)
     coA(v+2 )= alvv(i, j, 1, k, l)
     coA(v+3 )= alvh(i, j, 1, k, l)
     coA(v+4 )= alvu(i, j, 2, k, l)
     coA(v+5 )= alvv(i, j, 2, k, l)
     coA(v+6 )= alvh(i, j, 2, k, l)
     coA(v+7 )= alvu(i, j, 3, k, l)
     coA(v+8 )= alvv(i, j, 3, k, l)
     coA(v+9 )= alvh(i, j, 3, k, l)
     coA(v+10)= alvu(i, j, 4, k, l)
     coA(v+11)= alvv(i, j, 4, k, l)
     coA(v+12)= alvh(i, j, 4, k, l)
     coA(v+13)= alvu(i, j, 5, k, l)
     coA(v+14)= alvv(i, j, 5, k, l)
     coA(v+15)= alvh(i, j, 5, k, l)
     coA(v+16)= alvu(i, j, 6, k, l)
     coA(v+17)= alvv(i, j, 6, k, l)
     coA(v+18)= alvh(i, j, 6, k, l)
     coA(v+19)= alvu(i, j, 7, k, l)
     coA(v+20)= alvv(i, j, 7, k, l)
     coA(v+21)= alvh(i, j, 7, k, l)
     cpcoA(v+3)=cpnlvh(i, j, 1, k, l)
     cpcoA(v+6)=cpnlvh(i, j, 2, k, l)
     cpcoA(v+9)=cpnlvh(i, j, 3, k, l)
     cpcoA(v+12)=cpnlvh(i, j, 4, k, l)
     cpcoA(v+15)=cpnlvh(i, j, 5, k, l)
     cpcoA(v+18)=cpnlvh(i, j, 6, k, l)
     cpcoA(v+21)=cpnlvh(i, j, 7, k, l)


     ccoA(v+1 )= llvu(i, j, 1, k, l)
     ccoA(v+2 )= llvv(i, j, 1, k, l)
     ccoA(v+3 )= llvh(i, j, 1, k, l)
     ccoA(v+4 )= llvu(i, j, 2, k, l)
     ccoA(v+5 )= llvv(i, j, 2, k, l)
     ccoA(v+6 )= llvh(i, j, 2, k, l)
     ccoA(v+7 )= llvu(i, j, 3, k, l)
     ccoA(v+8 )= llvv(i, j, 3, k, l)
     ccoA(v+9 )= llvh(i, j, 3, k, l)
     ccoA(v+10)= llvu(i, j, 4, k, l)
     ccoA(v+11)= llvv(i, j, 4, k, l)
     ccoA(v+12)= llvh(i, j, 4, k, l)
     ccoA(v+13)= llvu(i, j, 5, k, l)
     ccoA(v+14)= llvv(i, j, 5, k, l)
     ccoA(v+15)= llvh(i, j, 5, k, l)
     ccoA(v+16)= llvu(i, j, 6, k, l)
     ccoA(v+17)= llvv(i, j, 6, k, l)
     ccoA(v+18)= llvh(i, j, 6, k, l)
     ccoA(v+19)= llvu(i, j, 7, k, l)
     ccoA(v+20)= llvv(i, j, 7, k, l)
     ccoA(v+21)= llvh(i, j, 7, k, l)
     enddo
     do l = 1, nl
     v = row+24*nl*2+24*(l-1)
     coA(v+1 )= alhu(i, j, 1, k, l)
     coA(v+2 )= alhv(i, j, 1, k, l)
     coA(v+3 )= alhh(i, j, 1, k, l)
     coA(v+4 )= alhu(i, j, 2, k, l)
     coA(v+5 )= alhv(i, j, 2, k, l)
     coA(v+6 )= alhh(i, j, 2, k, l)
     coA(v+7 )= alhu(i, j, 3, k, l)
     coA(v+8 )= alhv(i, j, 3, k, l)
     coA(v+9 )= alhh(i, j, 3, k, l)
     coA(v+10)= alhu(i, j, 4, k, l)
     coA(v+11)= alhv(i, j, 4, k, l)
     coA(v+12)= alhh(i, j, 4, k, l)
     coA(v+13)= alhu(i, j, 5, k, l)
     coA(v+14)= alhv(i, j, 5, k, l)
     coA(v+15)= alhh(i, j, 5, k, l)
     coA(v+16)= alhu(i, j, 6, k, l)
     coA(v+17)= alhv(i, j, 6, k, l)
     coA(v+18)= alhh(i, j, 6, k, l)
     coA(v+19)= alhu(i, j, 7, k, l)
     coA(v+20)= alhv(i, j, 7, k, l)
     coA(v+21)= alhh(i, j, 7, k, l)


     ccoA(v+1 )= llhu(i, j, 1, k, l)
     ccoA(v+2 )= llhv(i, j, 1, k, l)
     ccoA(v+3 )= llhh(i, j, 1, k, l)
     ccoA(v+4 )= llhu(i, j, 2, k, l)
     ccoA(v+5 )= llhv(i, j, 2, k, l)
     ccoA(v+6 )= llhh(i, j, 2, k, l)
     ccoA(v+7 )= llhu(i, j, 3, k, l)
     ccoA(v+8 )= llhv(i, j, 3, k, l)
     ccoA(v+9 )= llhh(i, j, 3, k, l)
     ccoA(v+10)= llhu(i, j, 4, k, l)
     ccoA(v+11)= llhv(i, j, 4, k, l)
     ccoA(v+12)= llhh(i, j, 4, k, l)
     ccoA(v+13)= llhu(i, j, 5, k, l)
     ccoA(v+14)= llhv(i, j, 5, k, l)
     ccoA(v+15)= llhh(i, j, 5, k, l)
     ccoA(v+16)= llhu(i, j, 6, k, l)
     ccoA(v+17)= llhv(i, j, 6, k, l)
     ccoA(v+18)= llhh(i, j, 6, k, l)
     ccoA(v+19)= llhu(i, j, 7, k, l)
     ccoA(v+20)= llhv(i, j, 7, k, l)
     ccoA(v+21)= llhh(i, j, 7, k, l)

     enddo
     enddo   
     enddo
     enddo

     end
!*****************************************************************************
     subroutine jacob(un, sig)

     use m_par
     use m_usr
     implicit none
!*    import/export
     real    un(ndim)
     real    sig
!*    local
     real    u(0:n, 0:m+1, 1:nl), v(0:n+1, 0:m, 1:nl), h(0:n+1, 0:m+1, 1:nl)
     real    eps2
     integer i, j, k

     eps2 = par(3)
     cpnluh = 0
     cpnlvh = 0

     Frc = 0.0
     call usol(un, u, v, h)

     call lin
     call nlin_rhs(u, v, h)
     call nlin_jac(u, v, h)


     Aluu = Lluu+Nluu
     Aluv = Lluv+Nluv
     Aluh = Lluh+Nluh
     Alvu = Llvu+Nlvu
     Alvv = Llvv+Nlvv
     Alvh = Llvh+Nlvh
     Alhu = Llhu+Nlhu
     Alhv = Llhv+Nlhv
     Alhh = Llhh+Nlhh

     do k = 1, nl
     Aluu(:,:,4, k, k) = Aluu(:,:,4, k, k) + sig*eps2
     Alvv(:,:,4, k, k) = Alvv(:,:,4, k, k) + sig*eps2
     Alhh(:,:,4, k, k) = Alhh(:,:,4, k, k) + sig
     enddo

     do k = 1, nl
     lluu(:,:,4, k, k) = lluu(:,:,4, k, k) + sig*eps2
     llvv(:,:,4, k, k) = llvv(:,:,4, k, k) + sig*eps2
     llhh(:,:,4, k, k) = llhh(:,:,4, k, k) + sig
     enddo

     call boundaries
     call assemble    ! Convert to a global sparse matrix.
     !call writemat
     if (mix_flag) call mix_jac(un)  ! mix
     call land2mat(Frc)
     if (pmode .eq. 3) call intcond(Frc)
     call packa       ! Remove very small elements.
100   format(3i4, 2es12.4)
     end
!*****************************************************************************
     subroutine mrhs(un, b)
     use m_par
     use m_usr
     use m_mat
implicit none

!*     import/export
     real    un(ndim), b(ndim)
!*     local
     integer i, j, k, row
     real    au(ndim), cpau(ndim)
     real    u(0:n, 0:m+1, 1:nl), v(0:n+1, 0:m, 1:nl), h(0:n+1, 0:m+1, 1:nl)
     real    hbx(n, m), hby(n, m), eps2, mix(ndim)
     eps2 = par(3)
     !print *, "in user.f_mrhs l_hth =",l_hth
!      call depth_cor(un)
     call usol(un, u, v, h)       ! conversion state -> u, v, h

     call lin
     call nlin_rhs(u, v, h)      ! local matrix, non-linear
     if (mmode .eq. 0) then
     call gradhb(hbx, hby)
     do k = 1, nl
     do j = 1, m
     do i = 1, n 
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     frc(row+1)=-l_coef(k, nl+1)*par(12)*hbx(i, j)
     frc(row+2)=-l_coef(k, nl+1)*par(12)*hby(i, j)
     frc(row+3)= 0.0
     enddo
     enddo
     enddo
   else
     Frc = 0.0
   endif

     aluu = lluu+nluu        ! total local matrix
     aluv = lluv+nluv
     aluh = lluh+nluh
     alvu = llvu+nlvu
     alvv = llvv+nlvv
     alvh = llvh+nlvh
     alhu = llhu+nlhu
     alhv = llhv+nlhv
     alhh = llhh+nlhh

     call boundaries           ! implement boundary conditions

     call assemble             ! assemble global matrix from local
     if (mix_flag) call mix_jac(un)
     call land2mat(frc)        ! remove land
     if (pmode .eq. 3) call intcond(Frc)
     call packa                ! sparse row storage and removal of small elements

     call matavec(un, au)

     b = -au+frc

100   format(3i6, 2es12.4)
     end
!*****************************************************************************
     real function l2nrm(f, ndim)       ! 2-norm
     implicit none
     integer i, ndim
     real    f(ndim), nrm
!*
     nrm = 0.0
     do i = 1, ndim
     nrm = f(i)*f(i) + nrm
     enddo
     l2nrm = sqrt(nrm)
!*
     end
!*****************************************************************************
     real function linrm(f, ndim)       ! max-norm
     implicit none
     integer i, ndim
     real    f(ndim), nrm
!*
     nrm = 0.0
     do i = 1, ndim
     if ( abs(f(i)).gt.nrm ) nrm = abs(f(i))
     enddo
     linrm = nrm
!*
     end
!*****************************************************************************
     subroutine lin
!*     shallow water equations
!*     produce local element matrices for linear operators
     use m_par
     use m_usr
     implicit none

!*     local
     integer k, k2
     real    u(n, m, np), uxx(n, m, np), uyy(n, m, np), fv(n, m, np), hx(n, m, np)
     real    uy(n, m, np), vy(n, m, np)
     real    ucsi(n, m, np), vcsi(n, m, np), uxsncsi(n, m, np), vxsncsi(n, m, np)
     real    v(n, m, np), vxx(n, m, np), vyy(n, m, np), fu(n, m, np), hy(n, m, np)
     real    ek, ri, rb, dum
!*
     rb  = par(2)  ! bottom friction
     ri  = par(8)  ! interface friction
     ek  = par(4)  ! Ekman number

     lluu = 0.0
     lluv = 0.0
     lluh = 0.0
     llvu = 0.0
     llvv = 0.0
     llvh = 0.0
     llhu = 0.0
     llhv = 0.0
     llhh = 0.0

     u = 0.0
     ucsi = 0.0
     uy = 0.0
     uxx = 0.0
     uyy = 0.0
     vxsncsi = 0.0
     fv = 0.0
     hx = 0.0

     v = 0.0
     vcsi = 0.0
     vy = 0.0
     vxx = 0.0
     vyy = 0.0
     vxsncsi = 0.0
     fu = 0.0
     hy = 0.0

     if (mix_flag) dum = 0.0

! Equation for u
     call uderiv(1, u )
     if(.not.mix_flag) then
     call uderiv(2, ucsi)
     call uderiv(6, vxsncsi)
     call uderiv(3, uy)
     call uderiv(4, uxx)
     call uderiv(5, uyy)
   endif
     call coriolis(1, fv)
     call gradh(1, hx)
     do k = 1, nl
     lluu(:,:,:,k, k) = -ek*(-uy+uyy+uxx-ucsi)
     lluv(:,:,:,k, k) = -fv-ek*(- vxsncsi)
     ! pressure term
     do k2 = 1, nl
     lluh(:,:,:,k, k2) = l_coef(k, k2)*hx
     enddo
     enddo
! TODO: Restructure the following part of code. The nl .eq. 1 case can be 
! treated separately, avoing the other if-statements and improving 
! readablility.
     ! friction terms
     ! - top layer
     if (nl .gt. 1) then 
     lluu(:,:,:,1, 1) = lluu(:,:,:,1, 1) + ri*u
     lluu(:,:,:,1, 2) = -ri*u
   endif
     ! - intermediate layers (only for more than 2 layers)
! TODO: The if here is not needed since the loop will not execute if 
! nl <= 2, unless Fortran 66 behaviour would be selected.
     if (nl .gt. 2) then
     do k = 2, nl-1
     lluu(:,:,:,k, k-1) = -ri*u
     lluu(:,:,:,k, k)   =  lluu(:,:,:,k, k) + 2.0*ri*u
     lluu(:,:,:,k, k+1) = -ri*u
     enddo
   endif
     ! - bottom layer
     if (nl .gt. 1) then
     lluu(:,:,:,nl, nl-1) = -ri*u
     lluu(:,:,:,nl, nl) = lluu(:,:,:,nl, nl) + (ri+rb)*u
   else
     ! Single layer, nl = 1.
     lluu(:,:,:,nl, nl) = lluu(:,:,:,nl, nl) + rb*u
   endif

! Equation for v
     call vderiv(1, v  )
     if (.not.mix_flag) then
     call vderiv(2, vcsi)
     call vderiv(6, uxsncsi)
     call vderiv(3, vy )
     call vderiv(4, vxx)
     call vderiv(5, vyy)
   endif
     call coriolis(2, fu)
     call gradh(2, hy)
     do k = 1, nl
     llvu(:,:,:,k, k) = fu  - ek*( uxsncsi)
     llvv(:,:,:,k, k) = -ek*(-vy+vyy+vxx-vcsi)
     ! pressure term
     do k2 = 1, nl
     llvh(:,:,:,k, k2) = l_coef(k, k2)*hy
     enddo
     enddo
! TODO: Restructure the following part of code. The nl .eq. 1 case can be 
! treated separately, avoing the other if-statements and improving 
! readablility.
     ! friction terms
     ! - top layer
     if (nl .gt. 1) then
     llvv(:,:,:,1, 1) = llvv(:,:,:,1, 1) + ri*v
     llvv(:,:,:,1, 2) = -ri*v
   endif
     ! - intermediate layers (only for more than 2 layers)
     if (nl .gt. 2) then
     do k = 2, nl-1
     llvv(:,:,:,k, k-1) = -ri*u
     llvv(:,:,:,k, k)   =  llvv(:,:,:,k, k) + 2.0*ri*v
     llvv(:,:,:,k, k+1) = -ri*v
     enddo
   endif
     ! - bottom layer
     if (nl .gt. 1) then
     llvv(:,:,:,nl, nl-1) = -ri*v
     llvv(:,:,:,nl, nl) = llvv(:,:,:,nl, nl) + (ri+rb)*v 
   else
     llvv(:,:,:,nl, nl) = llvv(:,:,:,nl, nl) + rb*v
   endif

! no linear terms in equation for h

     end
!*****************************************************************************
     subroutine nlin_rhs(u, v, h)
!*     produce local matrices for nonlinear operators for calc of rhs

     use m_par
     use m_usr
     implicit none
!*     import/export
     real    u(0:n, 0:m+1, 1:nl), v(0:n+1, 0:m, 1:nl), h(0:n+1, 0:m+1, 1:nl)
!*     local
     real    uux(n, m, np), vuy(n, m, np)
     real    uvx(n, m, np), vvy(n, m, np)
     real     uv(n, m, np), uu(n, m, np)
     real    uhx(n, m, np), vhy(n, m, np)
     real    txh(n, m, np), tyh(n, m, np)
     real    lxh(n, m, np), lyh(n, m, np)
     real    sig, eps2
     integer i, j, k

     sig  = par(1)*par(11)
     eps2 = par(3)
! TODO KURT: Was it correct to disable the next call?
!      call forcing

     nluu = 0.0
     nluv = 0.0
     nluh = 0.0
     nlvu = 0.0
     nlvv = 0.0
     nlvh = 0.0
     nlhu = 0.0
     nlhv = 0.0
     nlhh = 0.0

     do k = 1, nl
     uux = 0.0
     vuy = 0.0
     uv = 0.0
     txh = 0.0
     uu = 0.0
     uvx = 0.0
     vvy = 0.0
     tyh = 0.0
     uhx = 0.0
     vhy = 0.0
     call unlin(3, uux, u(:,:,k), v(:,:,k), h(:,:,k))
     call unlin(5, vuy, u(:,:,k), v(:,:,k), h(:,:,k))
     call unlin(6, uv, u(:,:,k), v(:,:,k), h(:,:,k))
     nluu(:,:,:,k, k) =eps2*(-uv+uux+vuy)
     nluv(:,:,:,k, k) = 0.0

     if (potential) then
     call layer( 1, lxh, u(:,:,k), v(:,:,k), h(:,:,k))
     nluh(:,:,:,k, k) = l_coef(k, k)*lxh
   endif
     if (k .eq. 1) then
     call wind( 1, txh, u(:,:,k), v(:,:,k), h(:,:,k))
     nluh(:,:,:,k, k) = nluh(:,:,:,k, k)-sig*txh
     cpnluh(:,:,:,k, k)=-sig*txh
   endif
     call vnlin(6, uu, u(:,:,k), v(:,:,k), h(:,:,k))
     call vnlin(3, uvx, u(:,:,k), v(:,:,k), h(:,:,k))
     call vnlin(5, vvy, u(:,:,k), v(:,:,k), h(:,:,k))
     nlvu(:,:,:,k, k) = eps2*(uu)
     nlvv(:,:,:,k, k) = eps2*(uvx+vvy)
     if (potential) then
     call layer( 2, lyh, u(:,:,k), v(:,:,k), h(:,:,k))
     nlvh(:,:,:,k, k) = l_coef(k, k)*lyh
   endif
     if (k .eq. 1) then
     call wind( 2, tyh, u(:,:,k), v(:,:,k), h(:,:,k))
     nlvh(:,:,:,k, k) = nlvh(:,:,:,k, k)-sig*tyh
     cpnlvh(:,:,:,k, k)=-sig*tyh
   endif
!*
     call hnlin(2, uhx, u(:,:,k), v(:,:,k), h(:,:,k))
     call hnlin(4, vhy, u(:,:,k), v(:,:,k), h(:,:,k))
     nlhh(:,:,:,k, k) = uhx+vhy
     enddo
     end
!*****************************************************************************
     subroutine nlin_jac(u, v, h)
!*     produce local matrices for nonlinear operators for calc of jacobian
     !implicit none
     use m_par
     use m_usr
     implicit none
!*     import/export
     real    u(0:n, 0:m+1, 1:nl), v(0:n+1, 0:m, 1:nl), h(0:n+1, 0:m+1, 1:nl)
!*     local
     real    urux(n, m, np), uurx(n, m, np), vruy(n, m, np), vury(n, m, np)
     real    urvx(n, m, np), uvrx(n, m, np), vrvy(n, m, np), vvry(n, m, np)
     real    urhx(n, m, np), uhrx(n, m, np), vrhy(n, m, np), vhry(n, m, np)
     real    urv(n, m, np), uvr(n, m, np), uru(n, m, np)
     real    hrux(n, m, np), hurx(n, m, np), hrvy(n, m, np), hvry(n, m, np)
     real    txh(n, m, np), tyh(n, m, np)
     real    lxh(n, m, np), lyh(n, m, np)
     real    sig, eps2
     integer i, j, k, k2, l
     real dum
!*
     sig = par(1)*par(11)
     eps2 = par(3)
! TODO Kurt: Was it correct to disable the next call?
!      call forcing

     nluu = 0.0
     nluv = 0.0
     nluh = 0.0
     nlvu = 0.0
     nlvv = 0.0
     nlvh = 0.0
     nlhu = 0.0
     nlhv = 0.0
     nlhh = 0.0

     do k = 1, nl
     urux = 0.0
     uurx = 0.0
     vruy = 0.0
     vury = 0.0
     urv = 0.0
     uvr = 0.0
     txh = 0.0
     urvx = 0.0
     uvrx = 0.0
     vrvy = 0.0
     vvry = 0.0
     uru = 0.0
     tyh = 0.0
     urhx = 0.0
     uhrx = 0.0
     vrhy = 0.0
     vhry = 0.0

     call unlin(2, urux, u(:,:,k), v(:,:,k), h(:,:,k))
     call unlin(3, uurx, u(:,:,k), v(:,:,k), h(:,:,k))
     call unlin(4, vruy, u(:,:,k), v(:,:,k), h(:,:,k))
     call unlin(5, vury, u(:,:,k), v(:,:,k), h(:,:,k))
     call unlin(6, urv, u(:,:,k), v(:,:,k), h(:,:,k))
     call unlin(7, uvr, u(:,:,k), v(:,:,k), h(:,:,k))
     nluu(:,:,:,k, k) =eps2*(-urv+uurx+urux+vury)
     nluv(:,:,:,k, k) =eps2*(-uvr+vruy)
     if (potential) then
     call layer( 1, lxh, u(:,:,k), v(:,:,k), h(:,:,k))
     nluh(:,:,:,k, k) = -l_coef(k, k)*lxh
   endif
     if (k .eq. 1) then 
     call wind( 1, txh, u(:,:,k), v(:,:,k), h(:,:,k))
     nluh(:,:,:,k, k) = nluh(:,:,:,k, k) + sig*txh
   endif

     call vnlin(2, urvx, u(:,:,k), v(:,:,k), h(:,:,k))
     call vnlin(3, uvrx, u(:,:,k), v(:,:,k), h(:,:,k))
     call vnlin(4, vrvy, u(:,:,k), v(:,:,k), h(:,:,k))
     call vnlin(5, vvry, u(:,:,k), v(:,:,k), h(:,:,k))
     call vnlin(6, uru, u(:,:,k), v(:,:,k), h(:,:,k))
     nlvu(:,:,:,k, k) = eps2*(2*uru+urvx)
     nlvv(:,:,:,k, k) = eps2*(uvrx+vvry+vrvy)
     if (potential) then
     call layer( 2, lyh, u(:,:,k), v(:,:,k), h(:,:,k))
     nlvh(:,:,:,k, k) = -l_coef(k, k)*lyh
   endif
     if (k .eq. 1) then
     call wind( 2, tyh, u(:,:,k), v(:,:,k), h(:,:,k))
     nlvh(:,:,:,k, k) =nlvh(:,:,:,k, k)+sig*tyh
   endif

     call hnlin(1, urhx, u(:,:,k), v(:,:,k), h(:,:,k))
     call hnlin(2, uhrx, u(:,:,k), v(:,:,k), h(:,:,k))
     call hnlin(3, vrhy, u(:,:,k), v(:,:,k), h(:,:,k))
     call hnlin(4, vhry, u(:,:,k), v(:,:,k), h(:,:,k))
     nlhu(:,:,:,k, k) = urhx
     nlhv(:,:,:,k, k) = vrhy
     nlhh(:,:,:,k, k) = uhrx+vhry
     enddo
     end
!*****************************************************************************
     subroutine mixing(un, mix)
!*     produce local matrices for nonlinear operators for calc of jacobian
     !implicit none
     use m_par
     use m_usr
     implicit none
!*     import/export
     real    un(ndim), mix(ndim)
     real    u(0:n, 0:m+1, 1:nl), v(0:n+1, 0:m, 1:nl), h(0:n+1, 0:m+1, 1:nl)
     integer i, j, k, row
     real AH_h

     call usol(un, u, v, h)       ! conversion state -> u, v, h

     mix = 0.0
     print *, nl,l_hth(1),l_hth(2)
      
     do i = 1, nl-1
     h(0  ,:  ,i)=l_hth(i)-l_hth(i+1)
     h(n+1, :  ,i)=l_hth(i)-l_hth(i+1)
     h(:  ,0  ,i)=l_hth(i)-l_hth(i+1)
     h(:  ,m+1, i)=l_hth(i)-l_hth(i+1)
     enddo
     h(0  ,:  ,nl)=l_hth(nl)
     h(n+1, :  ,nl)=l_hth(nl)
     h(:  ,0  ,nl)=l_hth(nl)
     h(:  ,m+1, nl)=l_hth(nl)

     do k = 1, nl
     do j = 1, m
     do i = 1, n
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     ! *** U-momentum: row+1
     ! AH_x*u_x
     mix(row+1)=mix(row+1)+&
     (AH_h(h(i+1, j, k))-Ah_h(h(i, j, k)))*&  ! E
     (u(i+1, j, k)-u(i-1, j, k))/&
     (2.0*dx**2*cos(y(j))**2)
     ! AH   * u_xx
     mix(row+1)=mix(row+1)+&
     (AH_h(h(i+1, j, k))+Ah_h(h(i, j, k)))*&  ! E
     (u(i+1, j, k)-2.0*u(i, j, k)+u(i-1, j, k))/&
     (2.0*dx**2*cos(y(j))**2)
     ! AH_y*u_y
     mix(row+1)=mix(row+1)+&
     (AH_h(h(i+1, j+1, k))-Ah_h(h(i+1, j-1, k))+&     ! E, N, S
     AH_h(h(i, j+1, k))-Ah_h(h(i, j-1, k)))*&
     (u(i, j+1, k)-u(i, j-1, k))/&
     (8.0*dy**2)
     ! AH   * u_yy
     mix(row+1)=mix(row+1)+&
     (AH_h(h(i+1, j, k))+Ah_h(h(i, j, k)))*&  ! E
     (u(i, j+1, k)-2.0*u(i, j, k)+u(i, j-1, k))/&
     (2.0*dy**2)
     ! AH   * u_y  * tan y
     mix(row+1)=mix(row+1)-&
     (AH_h(h(i+1, j, k))+AH_h(h(i, j, k)))*&  ! E
     (u(i, j+1, k)-u(i, j-1, k))*&
     tan(y(j))/(4.0*dy)

     ! *** V-momentum: row+1
     ! AH_x*v_x
     mix(row+2)=mix(row+2)+&
     (AH_h(h(i+1, j+1, k))-Ah_h(h(i-1, j+1, k))+&     !E, N, W
     AH_h(h(i+1, j, k))-Ah_h(h(i-1, j, k)))*&
     (v(i+1, j, k)-v(i-1, j, k))/&
     (8.0*dx**2*cos(yv(j))**2)
     ! AH   * v_xx
     mix(row+2)=mix(row+2)+&
     (AH_h(h(i, j+1, k))+Ah_h(h(i, j, k)))*&  ! N
     (v(i+1, j, k)-2.0*v(i, j, k)+v(i-1, j, k))/&
     (2.0*dx**2*cos(yv(j))**2)
     ! AH_y*v_y
     mix(row+2)=mix(row+2)+&
     (AH_h(h(i, j+1, k))-Ah_h(h(i, j, k)))*&  ! N
     (v(i, j+1, k)-v(i, j-1, k))/&
     (2.0*dy**2)
     ! AH   * v_yy
     mix(row+2)=mix(row+2)+&              ! N
     (AH_h(h(i, j+1, k))+Ah_h(h(i, j, k)))*&
     (v(i, j+1, k)-2.0*v(i, j, k)+v(i, j-1, k))/&
     (2.0*dy**2)
     ! AH   * v_y  * tan y
     mix(row+2)=mix(row+2)-&
     (AH_h(h(i, j+1, k))+AH_h(h(i, j, k)))*&       ! N
     (v(i, j+1, k)-v(i, j-1, k))*&
     tan(yv(j))/(4.0*dy)
     enddo
     enddo
     enddo

     end
!*****************************************************************************
     subroutine mix_el(imix, jmix, nmix)
!*     produce local matrices for nonlinear operators for calc of jacobian
     !implicit none
     use m_par
     use m_usr
     implicit none
!*     import/export
     integer  imix(ndim*25), jmix(ndim*25), nmix
     integer i, j, k, row, l

     imix = 0
     jmix = 0
     nmix = 0
     l    = 0

     do k = 1, nl
     ! interior
     do j = 2, m-1
     do i = 2, n-1
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     imix(l+1:l+11)=row+1
     jmix(l+1)=row-3*n     + 1
     jmix(l+2)=row-3*n     + 3
     jmix(l+3)=row-3*n+3 + 3
     jmix(l+4)=row       - 3+1
     jmix(l+5)=row           + 1
     jmix(l+6)=row           + 3
     jmix(l+7)=row       + 3+1
     jmix(l+8)=row       + 3+3
     jmix(l+9)=row+3*n     + 1
     jmix(l+10)=row+3*n     + 3
     jmix(l+11)=row+3*n+3 + 3

     imix(l+12:l+22)=row+2
     jmix(l+12)=row-3*n     + 2
     jmix(l+13)=row       - 3+2
     jmix(l+14)=row       - 3+3
     jmix(l+15)=row           + 2
     jmix(l+16)=row           + 3
     jmix(l+17)=row       + 3+2
     jmix(l+18)=row       + 3+3
     jmix(l+19)=row+3*n-3 + 3
     jmix(l+20)=row+3*n     + 2
     jmix(l+21)=row+3*n     + 3
     jmix(l+22)=row+3*n+3 + 3

     l = l+22
     enddo
     enddo
     ! Western Boundary
     do j = 2, m-1
     i = 1
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     imix(l+1:l+10)=row+1
     jmix(l+1)=row-3*n     + 1
     jmix(l+2)=row-3*n     + 3
     jmix(l+3)=row-3*n+3 + 3
     jmix(l+4)=row           + 1
     jmix(l+5)=row           + 3
     jmix(l+6)=row       + 3+1
     jmix(l+7)=row       + 3+3
     jmix(l+8)=row+3*n     + 1
     jmix(l+9)=row+3*n     + 3
     jmix(l+10)=row+3*n+3 + 3

     imix(l+11:l+18)=row+2
     jmix(l+11)=row-3*n     + 2
     jmix(l+12)=row           + 2
     jmix(l+13)=row           + 3
     jmix(l+14)=row       + 3+2
     jmix(l+15)=row       + 3+3
     jmix(l+16)=row+3*n     + 2
     jmix(l+17)=row+3*n     + 3
     jmix(l+18)=row+3*n+3 + 3

     l = l+18
     ! Eastern Boundary
     i = n
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)            
     imix(l+1:l+8)=row+2
     jmix(l+1)=row-3*n     + 2
     jmix(l+2)=row       - 3+2
     jmix(l+3)=row       - 3+3
     jmix(l+4)=row           + 2
     jmix(l+5)=row           + 3
     jmix(l+6)=row+3*n-3 + 3
     jmix(l+7)=row+3*n     + 2
     jmix(l+8)=row+3*n     + 3

     l = l+8
     enddo
     do i = 2, n-1
     ! Southern Boundary
     j = 1
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     imix(l+1:l+8)=row+1
     jmix(l+1)=row       - 3+1
     jmix(l+2)=row         + 1
     jmix(l+3)=row         + 3
     jmix(l+4)=row       + 3+1
     jmix(l+5)=row       + 3+3
     jmix(l+6)=row+3*n     + 1
     jmix(l+7)=row+3*n     + 3
     jmix(l+8)=row+3*n+3 + 3

     imix(l+9:l+18)=row+2
     jmix(l+9)=row       - 3+2
     jmix(l+10)=row       - 3+3
     jmix(l+11)=row         + 2
     jmix(l+12)=row         + 3
     jmix(l+13)=row       + 3+2
     jmix(l+14)=row       + 3+3
     jmix(l+15)=row+3*n-3 + 3
     jmix(l+16)=row+3*n     + 2
     jmix(l+17)=row+3*n     + 3
     jmix(l+18)=row+3*n+3 + 3

     l = l+18

     ! Northern Boundary
     j = m
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     imix(l+1:l+8)=row+1
     jmix(l+1)=row-3*n     + 1
     jmix(l+2)=row-3*n     + 3
     jmix(l+3)=row-3*n+3 + 3
     jmix(l+4)=row       - 3+1
     jmix(l+5)=row         + 1
     jmix(l+6)=row         + 3
     jmix(l+7)=row       + 3+1
     jmix(l+8)=row       + 3+3

     l = l+8
     enddo
     ! South-West Corner
     i = 1
     j = 1
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     imix(l+1:l+7)=row+1
     jmix(l+1)=row         + 1
     jmix(l+2)=row         + 3
     jmix(l+3)=row       + 3+1
     jmix(l+4)=row       + 3+3
     jmix(l+5)=row+3*n     + 1
     jmix(l+6)=row+3*n     + 3
     jmix(l+7)=row+3*n+3 + 3

     imix(l+8:l+14)=row+2
     jmix(l+8)=row         + 2
     jmix(l+9)=row         + 3
     jmix(l+10)=row       + 3+2
     jmix(l+11)=row       + 3+3
     jmix(l+12)=row+3*n     + 2
     jmix(l+13)=row+3*n     + 3
     jmix(l+14)=row+3*n+3 + 3

     l = l+14
     ! South-East Corner
     i = n
     j = 1
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     imix(l+1:l+7)=row+2
     jmix(l+1)=row       - 3+2
     jmix(l+2)=row       - 3+3
     jmix(l+3)=row         + 2
     jmix(l+4)=row         + 3
     jmix(l+5)=row+3*n-3 + 3
     jmix(l+6)=row+3*n     + 2
     jmix(l+7)=row+3*n     + 3

     l = l+7
     ! North-West Corner
     i = 1
     j = m
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     imix(l+1:l+7)=row+1
     jmix(l+1)=row-3*n     + 1
     jmix(l+2)=row-3*n     + 3
     jmix(l+3)=row-3*n+3 + 3
     jmix(l+4)=row           + 1
     jmix(l+5)=row           + 3
     jmix(l+6)=row       + 3+1
     jmix(l+7)=row       + 3+3

     l = l+7
     enddo

     nmix = l

     end
!*****************************************************************************
     subroutine mix_jac(un)
     !implicit none

     use m_par
     use m_usr
     use m_mat
     implicit none      
     integer i, j, k, numgrp
     real un(ndim), mix(ndim), d(ndim), und(ndim), mixd(ndim)
     real fjac(idim)
     logical col
     real time1, time2, time3

     call mixing(un, mix)

     do numgrp = 1, maxgrp
     do j = 1, ndim
     d(j)=0.0
     if (ngrp(j) .eq. numgrp) d(j) = 1.0e-4
     und(j)=un(j)+d(j)
     enddo
     call mixing(und, mixd)
     do i = 1, ndim
     mixd(i)=mixd(i)-mix(i)
     enddo
     col=.false.
     if (col) then
     call fdjs(ndim, ndim, col, irow, jpntr, ngrp, numgrp, d, mixd, fjac)
   else
     call fdjs(ndim, ndim, col, icol, ipntr, ngrp, numgrp, d, mixd, fjac)
   endif
     enddo

     ! note: row = i, columns are icol(ipntr(i):ipntr(i+1)-1)  

     numgrp = 0
     do i = 1, ndim  ! loop over rows
     do k = ipntr(i), ipntr(i+1)-1  ! loop over columns (in approx.)
     do j = begA(i), begA(i+1)-1
     if (icol(k).eq.jcoA(j)) then  ! Add+Exit
     coA(j)=coA(j)-fjac(k)*par(4)
     exit
   endif  
     enddo
     enddo
     enddo

     end
!*****************************************************************************
     real function AH_h(val)
     ! implicit none
     use m_par
     use m_usr
     implicit none
     real  val

     AH_h = 1.1/(0.1+val)
     end

!*****************************************************************************
     subroutine dirset(f, row, fd)
!*     inforce dirichlet condition
     ! implicit none
     use m_par
     use m_usr
     use m_mat
     implicit none
!*     import/export
     integer row
     real    f(ndim), fd
!*     local
     integer v
!*
     f(row)  =  fd
     do v = begA(row), begA(row+1)-1
     coA(v)= 0.0
     if (jcoA(v).eq.row) coA(v) = 1.0
     enddo
!*
     end
!*****************************************************************************
     subroutine matrix_swe(un, type, sig)
!     construct the jacobian a or the 'mass' matrix b
!     type = 1: put b in coB and a in coA 
!     type = 2: put b in coB and a-sig*b in coA
!     type = 3: Only compute b, put in coB.
     !implicit none
     use m_par
     use m_usr
     implicit none
!     Arguments
     integer type
     real    un(ndim)
     real    sig
     call massB
     if (type .eq. 2) then
!       put a-sig*b in coA
     call jacob(un, sig)

   elseif (type .eq. 1) then
!       put a in coA

     call jacob(un, 0.0)

   endif

     end
!*****************************************************************************
     subroutine massB
!     Construct the mass matrix B.
     !implicit none
     use m_par
     use m_usr
     use m_mat
 implicit none
!     Local variabeles
     integer i, j, k, row, idx
     real    eps2, dirb(ndim)
!*
     eps2 = par(3)  ! Rossby number.
!*     put b in coB
!*     b is a diagonal matrix. (once suffices)
     dirb = 1.0
     do idx = 1, ndir
     dirb(dir(idx)) = par(10)
     enddo
!*
     do k = 1, nl
     do j = 1, m
     do i = 1, n
     row =  3*(n*m*(k-1)+n*(j-1) + i-1)
     coB(row+1) = -eps2*dirb(row+1)
     coB(row+2) = -eps2*dirb(row+2)
     coB(row+3) =  -1.0*dirb(row+3)
     enddo
     enddo
     enddo
     do k = 1, nl
     i = n
     do j = 1, m
     if (.not.(periodic)) then
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     coB(row+1)=0.0
   endif
     enddo
     j = m
     do i = 1, n
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     coB(row+2)=0.0
     enddo
     enddo

     if (pmode .eq. 3) then
     do k = 1, nl
     row = 3*(n*m*(k-1)+n*(m-1)+n-1)+3      
     coB(row) = 0.0            ! zero in b for integral condition
     enddo
   endif
     end

!****************************************************************************
     subroutine intcond(f)
!*     impose integral condition
     !implicit none
     use m_par
     use m_usr
     use m_mat
     implicit none
!*     import/export
     real    f(ndim)
!*     local
     integer i, j, v, k, row, k2
     real    atot, btot
     integer row2
!*
!*     replace h equation at ndim with an 'integral' condition for h
!*     write(99, *) 'int: ',row, ndim

     atot = 0.0
     btot = 0.0
     do j = 1, m
     do i = 1, n
     atot = atot+cos(y(j))*(1.0-real(lnd(i, j))*(1.0-par(10)))
     btot = btot+cos(y(j))*par(12)*hb(i, j)*(1.0-real(lnd(i, j))*(1.0-par(10)))
     enddo
     enddo

     do k = 1, nl
     row = 3*(n*m*(k-1)+n*(m-1)+n-1)+3
     if (k .eq. nl) then
     f(row)=l_hth(k)*atot-btot
   else  
     f(row)=(l_hth(k)-l_hth(k+1))*atot
   endif
     ! Make room in storage for constraint
     j = n*m-(bega(row+1)-bega(row))  ! shift

     if (k .ne. nl) then
     do i = bega(ndim+1)-1, bega(row+1), -1
     jcoa(i+j)=jcoa(i)
     coa(i+j)=coa(i)
     enddo
     bega(row+1:ndim+1)=bega(row+1:ndim+1)+j
   else
     bega(ndim+1)=bega(ndim+1)+j
   endif

     ! Set equal to zero
     do i = bega(row), bega(row+1)-1
     coa(i)=0.0
     jcoa(i)=0
     enddo

     v = bega(row)
     do j = 1, m
     do i = 1, n
     row2 =  3*(n*m*(k-1)+n*(j-1)+i-1)
     jcoa(v)= row2+3
     coa(v) = cos(y(j))*(1.0-real(lnd(i, j))*(1.0-par(10)))
     v = v+1
     enddo
     enddo
     enddo

     end
!***************************************************************************
     subroutine writemat
!*     write the matrix a to a file
     !implicit none
     use m_par
     use m_usr
     use m_mat
     implicit none
!*     local
     integer i

     write(99, *) 'write matrix'
!*     write(10, *) ndim
     open(11, file='c.beg')
     open(12, file='c.jco')
     open(13, file='c.co')
     do i = 1, ndim+1
     write(11, *) begA(i)
     enddo
     do i = 1, begA(ndim+1)-1
     write(12, *) jcoA(i)
     enddo
     do i = 1, begA(ndim+1)-1
     write(13, *) coA(i)
     enddo



     end
!*****************************************************************************
     subroutine boundaries
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
     aluu(n, j, 4, k, k)=1.0
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
     alvv(:,m, 4, k, k)=1.0
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
     lluu(n, j, 4, k, k)=1.0
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
     llvv(:,m, 4, k, k)=1.0
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
!*****************************************************************************
     subroutine fillcola
!*     fill the collomns of a
     !implicit none
     use m_par
     use m_usr
     use m_mat
     implicit none
!*     local
     integer i, j, k, l, v, is, js, row, bis, ks, f, mfin
!*
!*     2 5  
!*     1 4 7
!*       3 6
!      call writemat()
     !print *,'fillcola',n, m, nl
     is = 3      
     js = 3*n
     ks = 3*n*m
     bis = 3*(n-1)
     do k = 1, nl    ! loop over layers
     do j = 1, m   ! loop in y direction
     do i = 1, n  ! loop in x direction
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     begA(row+1)= 3*3*8*nl*(n*m*(k-1)+n*(j-1)+i-1)+1
     begA(row+2)= 3*3*8*nl*(n*m*(k-1)+n*(j-1)+i-1)+24*nl+1
     begA(row+3)= 3*3*8*nl*(n*m*(k-1)+n*(j-1)+i-1)+24*nl*2+1
     do f = 1, 3  ! loop over variables
     do l = 1, nl
     v = 3*3*8*nl*(n*m*(k-1)+n*(j-1)+i-1)+(f-1)*24*nl+24*(l-1)
     jcoA(v+1) = row-is   +ks*(l-k)+1
     jcoA(v+2) = row-is   +ks*(l-k)+2
     jcoA(v+3) = row-is   +ks*(l-k)+3
     jcoA(v+4) = row-is+js+ks*(l-k)+1
     jcoA(v+5) = row-is+js+ks*(l-k)+2
     jcoA(v+6) = row-is+js+ks*(l-k)+3
     jcoA(v+7) = row   -js+ks*(l-k)+1
     jcoA(v+8) = row   -js+ks*(l-k)+2
     jcoA(v+9) = row   -js+ks*(l-k)+3
     jcoA(v+10)= row      +ks*(l-k)+1
     jcoA(v+11)= row      +ks*(l-k)+2
     jcoA(v+12)= row      +ks*(l-k)+3
     jcoA(v+13)= row   +js+ks*(l-k)+1
     jcoA(v+14)= row   +js+ks*(l-k)+2
     jcoA(v+15)= row   +js+ks*(l-k)+3
     jcoA(v+16)= row+is-js+ks*(l-k)+1
     jcoA(v+17)= row+is-js+ks*(l-k)+2
     jcoA(v+18)= row+is-js+ks*(l-k)+3
     jcoA(v+19)= row+is   +ks*(l-k)+1
     jcoA(v+20)= row+is   +ks*(l-k)+2
     jcoA(v+21)= row+is   +ks*(l-k)+3
     jcoA(v+22)= row+is+js+ks*(l-k)+1
     jcoA(v+23)= row+is+js+ks*(l-k)+2
     jcoA(v+24)= row+is+js+ks*(l-k)+3
     enddo
     enddo
     enddo
     enddo
     enddo
     begA(ndim+1)=24*3*nl*n*m*nl+1
     if ((periodic).or.(pper)) then  ! periodic boundary conditions
     if (pper) then
     mfin = mper
   else
     mfin = m
   endif
     do k = 1, nl    ! loop over layers
     do j = 1, mfin   ! loop in y direction
     do i = 1, n, n-1  ! loop in x direction
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)
     do f = 1, 3  ! loop over variables
     do l = 1, nl
     v = 3*3*8*nl*(n*m*(k-1)+n*(j-1)+i-1)+(f-1)*24*nl+24*(l-1)
     if (i .eq. 1) then
     jcoa(v+1) = row+bis   +ks*(l-k)+1
     jcoa(v+2) = row+bis   +ks*(l-k)+2
     jcoa(v+3) = row+bis   +ks*(l-k)+3
     jcoa(v+4) = row+bis+js+ks*(l-k)+1
     jcoa(v+5) = row+bis+js+ks*(l-k)+2
   endif  
     if (i .eq. n) then
     jcoa(v+16)= row-bis-js+ks*(l-k)+1
     jcoa(v+17)= row-bis-js+ks*(l-k)+2
     jcoa(v+18)= row-bis-js+ks*(l-k)+3
     jcoa(v+19)= row-bis   +ks*(l-k)+1
     jcoa(v+20)= row-bis   +ks*(l-k)+2
     jcoa(v+21)= row-bis   +ks*(l-k)+3
     jcoa(v+22)= row-bis+js+ks*(l-k)+1
     jcoa(v+23)= row-bis+js+ks*(l-k)+2
     jcoa(v+24)= row-bis+js+ks*(l-k)+3
   endif
     enddo
     enddo
     enddo
     enddo
     enddo
   endif
!print *, (jcoa(i), i = 1700, 1799)
     end
!*****************************************************************************
     subroutine packa
!*     remove zero entries in a
     ! implicit none
     use m_par
     use m_usr
     use m_mat
     implicit none
!*     local
     integer vv, v, begin, i, j, k, l, row, oldsize, cpoldsize, cpvv, cpbegin, ccvv, ccoldsize, ccbegin
     integer bn, bs
     logical q
!*
!*     indikken van a door alleen elementen > 1e-12 mee
!*     te nemen
!*
 !Print *, '**********coA is*************= '



     ! remove small elements
     vv = 1
     cpvv = 1
     ccvv = 1
     cpbegA = begA
     ccbegA = begA
     cpjcoA = jcoA
     ccjcoA = jcoA
     oldsize = begA(ndim+1)-1
     cpoldsize = cpbegA(ndim+1)-1
     ccoldsize = ccbegA(ndim+1)-1 
     do k = 1, nl
     do j = 1, m
     do i = 1, n
     do l = 1, 3
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)+l
     begin = vv
     do v = begA(row), begA(row+1)-1
     if (abs(coA(v)).gt.1.0e-12) then
     coA(vv) = coA(v)
     jcoA(vv)=jcoA(v)
     vv = vv+1
   endif
     enddo
     begA(row)= begin
     enddo
     enddo
     enddo
     enddo
     begA(ndim+1)=vv
     do k = 1, nl
     do j = 1, m
     do i = 1, n
     do l = 1, 3
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)+l
     cpbegin = cpvv
     do v = cpbegA(row), cpbegA(row+1)-1
     if (abs(cpcoA(v)).gt.1.0e-12) then
     cpcoA(cpvv) = cpcoA(v)
     cpjcoA(cpvv)=cpjcoA(v)
     cpvv = cpvv+1
   endif
     enddo
     cpbegA(row)= cpbegin
     enddo
     enddo
     enddo
     enddo



     do k = 1, nl
     do j = 1, m
     do i = 1, n
     do l = 1, 3
     row = 3*(n*m*(k-1)+n*(j-1)+i-1)+l
     ccbegin = ccvv
     do v = ccbegA(row), ccbegA(row+1)-1
     if (abs(ccoA(v)).gt.1.0e-12) then
     ccoA(ccvv) = ccoA(v)
     ccjcoA(ccvv)=ccjcoA(v)
     ccvv = ccvv+1
   endif
     enddo
     ccbegA(row)= ccbegin
     enddo
     enddo
     enddo
     enddo
     cpbegA(ndim+1)=cpvv 
     ccbegA(ndim+1)=ccvv



     end
!*******************************************************************************
     SUBROUTINE matAvec(v1, v2)
!*     This multiplies sparse matrix A and vector v1 to vector v2
     !implicit none
     use m_par
     use m_usr
     use m_mat
     implicit none
     real     v1(ndim), v2(ndim) !,v3(ndim), v4(ndim), v5(ndim)
!*     LOCAL
     integer  i, v
!*
     v2 = 0.0
!           v3 = 0.0
!           v4 = 0.0
!           v5 = 0.0
     DO i = 1, ndim
     DO v = begA(i), begA(i+1)-1
     v2(i) = coA(v)*v1(jcoA(v)) + v2(i)
     ENDDO
     ! DO v = cpbegA(i), cpbegA(i+1)-1
     !!!!!!IS  THIS CORRECT?????? TO GET THE tau./h^2 = txh
     ! v4(i) = cpcoA(v)*1.0+v4(i)
     !v3(i) = cpcoA(v)*v1(cpjcoA(v)) + v3(i)
     ! v5(i) = cpcoA(v)*2.0+v5(i)
!print *,cpcoA(v), i
     ! ENDDO

     ENDDO
     ! do i = 1, ndim        
     ! print *,v1(i), v2(i)
     ! enddo
     ! !print *, '**************************'
     ! pause


     !  DO i = 1, ndim
     ! print *, v3(i), v4(i), v5(i), v1(i)
     !ENDDO
     !pause 
100   format(2i7, 2es12.4)
     END
!*******************************************************************************
     SUBROUTINE matBvec(v1, v2)
!*     This multiplies diagonal matrix B and vector v1 to vector v2
!*     B is a diagonal matrix
     !implicit none
     use m_par
     use m_usr
     use m_mat
     implicit none
     real     v1(ndim), v2(ndim)
!*     LOCAL
     integer  i
!*
     DO i = 1, ndim
     v2(i) = coB(i)*v1(i)
     ENDDO
!*
     END
!* ============================================================================ *
     subroutine check_int_cond(un, type)
     !implicit none

     use m_par
     use m_usr
     use m_mat
     implicit none
     real un(ndim)
     real dum
     integer i, j, k, l, row, type

     if (type .eq. 1) call matrix_swe (un, 1, 0.0)

     do l = 1, nl
     dum = 0.0
     row = 3*(n*m*(l-1)+n*(m-1)+n-1)+3
     do k = begA(row), begA(row+1)-1
     dum = dum+coA(k)*un(jcoA(k))
     enddo   
     write (99, 100) l, dum-frc(row), frc(row)
     enddo

100   format('- int.cond. layer',i2, ':',es12.4, ' diff. wrt.',es12.4)
     end      
!* ============================================================================ *
     subroutine check_int_cond2(un, out)
     !implicit none

     use m_par
     use m_usr
     use m_mat
     implicit none
     real un(ndim)
     real out(nl)
     integer i, j, k, l, row, type

     call matrix_swe (un, 1, 0.0)

     do l = 1, nl
     out(l)=0.0
     row = 3*(n*m*(l-1)+n*(m-1)+n-1)+3
     do k = begA(row), begA(row+1)-1
     out(l)=out(l)+coA(k)*un(jcoA(k))
     enddo   
     enddo

100   format('- int.cond. layer',i2, ':',es12.4, ' diff. wrt.',es12.4)
     end      
!* ============================================================================ *
!********************************************************************
     SUBROUTINE outreg(un, ifile)
!*     write solution
     !implicit none
     use m_par
     use m_usr
     implicit none
!*     IMPORT/EXPORT
     integer ifile
     real    un(ndim)
!*     LOCAL
     integer i, j, k, ifile2
     real    u(0:n, 0:m+1, 1:nl), v(0:n+1, 0:m, 1:nl), h(0:n+1, 0:m+1, 1:nl)
!*
     call usol(un, u, v, h)

     do k = 1, nl
     ifile2 = ifile+k-1
     write(ifile2, 999) n, m
     write(ifile2, 997) xmin, xmax, ymin, ymax
     write(ifile2, 997) l_hth
     DO i = 1, n
     DO j = 1, m
     write(ifile2, 998)&
     0.5*(u(i, j, k)+u(i-1, j, k)), &
     0.5*(v(i, j, k)+v(i, j-1, k)), &
     h(i, j, k), &
     tx(i, j), ty(i, j), &
     lnd(i, j), hb(i, j)
     ENDDO
     ENDDO
     enddo
 999  format(2i4)
 998  format(5es18.8E3, i5, es18.8E3)
 997  format(8es18.8E3)
     END
!********************************************************************
     SUBROUTINE outmodes(unr, uni, iMode, sigma, para)
!*     write solution
     !implicit none
     use m_par
     use m_usr
     implicit none
!*     IMPORT/EXPORT
     integer iMode
     real    unr(ndim), uni(ndim), para(npar)
!*     LOCAL
     integer i, j, k, ifile
     real    ur(0:n, 0:m+1, 1:nl), vr(0:n+1, 0:m, 1:nl), hr(0:n+1, 0:m+1, 1:nl)
     real    ui(0:n, 0:m+1, 1:nl), vi(0:n+1, 0:m, 1:nl), hi(0:n+1, 0:m+1, 1:nl)
     real    sigma(2)
     character*50 filename
!*
     call usol(unr, ur, vr, hr)
     call usol(uni, ui, vi, hi)

     ifile = 50
     filename(1:5) = 'Mode_'
     write(filename(6:7), '(i2.2)') iMode
     filename(8:11) = '_nl_'

     do k = 1, nl
     write(filename(12:13), '(i2.2)') k
     open(unit = ifile, file = filename)
     write(ifile, 999) n, m
     write(ifile, 997) xmin, xmax, ymin, ymax
     write(ifile, 997) l_hth
     write(ifile, 998) sigma
     DO i = 1, n
     DO j = 1, m
     write(ifile, 998)&
     0.5*(ur(i, j, k)+ur(i-1, j, k)), &
     0.5*(vr(i, j, k)+vr(i, j-1, k)), &
     hr(i, j, k), &
     0.5*(ui(i, j, k)+ui(i-1, j, k)), &
     0.5*(vi(i, j, k)+vi(i, j-1, k)), &
     hi(i, j, k), &
     tx(i, j), ty(i, j), &
     lnd(i, j), par(12)*hb(i, j)
     ENDDO
     ENDDO
     close(ifile)
     enddo
 999  format(2i4)
 998  format(8es18.8E3, i5, es18.8E3)
 997  format(8es18.8E3)
     END
!********************************************************************
     SUBROUTINE depth_cor(un)
!*     write solution
     !implicit none
     use m_par
     use m_usr
     implicit none
!*     IMPORT/EXPORT
     integer ifile
     real    un(ndim)
!*     LOCAL
     integer i, j, k, ifile2
     real    u(0:n, 0:m+1, 1:nl), v(0:n+1, 0:m, 1:nl), h(0:n+1, 0:m+1, 1:nl)
!*
     call usol(un, u, v, h)

     do k = 1, nl
     DO i = 1, n
     DO j = 1, m
     if (h(i, j, k).lt.par(14)) h(i, j, k)=par(14)
     ENDDO
     ENDDO
     enddo
     END

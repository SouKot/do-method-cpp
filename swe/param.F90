!* ============================================================================ *
!*
!* Multi-Layer Shallow Water Model
!*
!* Version 1.4-20/05/2008-ATvS
!*
!* ============================================================================ *
! Model:
      module m_par
        integer, parameter:: mmode = 1			! 0: full model; 1: n.5 layer model
        integer   np, npar
        parameter (np = 7, npar = 20)

! Grid: 
        logical, parameter:: pper=.false.		! partial boundary conditions
        integer, parameter:: mper = 5			! number of point in partial periodic boundary
        real:: pmin = 0.0 			! psi min [in degrees] -> conversion 2 rads in init
        real:: pmax = 0.0 			! psi max [in degrees] -> conversion 2 rads in init
        real:: thmin = 0.0			! theta min [in degrees] -> conversion 2 rads in init
        real:: thmax = 0.0		 	! theta max [in degrees] -> conversion 2 rads in init
        real:: pgmin = 0.0 			! psi min [in degrees] -> conversion 2 rads in init
        real:: pgmax = 0.0 			! psi max [in degrees] -> conversion 2 rads in init
        real:: tghmin = 0.0			! theta min [in degrees] -> conversion 2 rads in init
        real:: tghmax = 0.0		 	! theta max [in degrees] -> conversion 2 rads in init


! Conversion quantities:
        real,    parameter:: hdim = 5.0e+02		! Layer heigth
        real,    parameter:: omegadim = 7.292e-05		! Angular velocity of the earth
        real,    parameter:: r0dim = 6.37e+06		! Radius of the earth
        real,    parameter:: udim = 1.0e-01		! Velocity scale
        real,    parameter:: gdim = 9.8		! Gravity
        real,    parameter:: rhodim = 1.0235e+03		! Densiity
        real,    parameter:: taudim = 1.0		! Amplitude of the windstress
        real,    parameter:: Ahdim = 2.20e+4
        real,    parameter:: bfric = 0.0e-0		! Bottom friction
        real,    parameter:: ifric = 5.0e-4		! Interfacial friction
        real,    parameter:: h0 = 10.0/hdim		! Minumum layer thickness when using potential for layer thicknesses
        real,    parameter:: hn = 12.0			! Exponent in potential formulaation for layer thicknesses

! Flags and switches:
        logical, parameter:: windfromdata=.false.	! windstress from data yes/no
        logical, parameter:: potential=.false.		! potential for layer thinknesses yes/no
        logical, parameter:: chk_intcond=.false.		! checking of integral conditions yes/no
        logical, parameter:: mix_flag=.false.		! Obsolete

! Newton iteration:
        integer, parameter:: NewtonMethod = 2		! Newton-Method: 1) standard, 2) adaptive shamanskii, 3) Newton-Chord
        integer, parameter:: OptNewtonIter = 6		! Optimal number of Newton iterations
        integer, parameter:: predictor = 1			! predictor: 1) Euler, 2) secant

! Topography & continents: 
        logical, parameter:: cont=.true.		! continents yes/no
        logical, parameter:: topo=.false.		! topography yes/no
        logical, parameter:: smoothed=.false.		! apply smoothing yes/no
        integer, parameter:: ism = 4            		! number of smoothings of topography data

! Interface equilibrium heigth and densities
        ! equilibrium heigth of interfaces
        !   note: bottom or top n+1-th layer at z = 0
        !   note: interface with atmosphere at z = 1
        !   note: l_hth(i) is interface at top of each layer
        !   note: for full model with topography, make sure topography remains within the bottom layer

     !double precision, dimension(2):: l_hth =(/ 1.0, 0.875/)
     double precision, dimension(1):: l_hth =(/ 1.0/)


        ! layer densities
        !   note: for n.5 layer model, density of the n+1-th layer is equal to rhodim
        !   note: for n.5 layer model, make sure rhodim > density of n-th layer
        double precision, dimension(2):: l_den = (/ 1.0235e+3, 1.02810575e+03/)
!contains

        !subroutine allocate_grid(Xmin, Xmax, Ymin, Ymax)
     ! real:: Xmin, Xmax, Ymin, Ymax
        !pmin = Xmin
        !pmax = Xmax
        !thmin = Ymin
        !thmax = Ymax
        !end subroutine


      end module m_par



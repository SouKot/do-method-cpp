cmake_minimum_required(VERSION 3.8)
project(DO C CXX Fortran)
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
#set(CMAKE_CXX_COMPILER_ID "Intel")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  message(STATUS "Using GNU compiler flags")

  if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	  #set(CMAKE_C_COMPILER "scorep mpicc")
	  #set(CMAKE_CXX_COMPILER "scorep mpic++")
	  #set(CMAKE_Fortran_COMPILER "scorep gfortran")
    set(CMAKE_C_FLAGS        "${CMAKE_C_FLAGS} -O2  -pg  -g -lpthread -lrt  -Wno-deprecated-declarations")
    set(CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS} -O2  -pg  -g -DGNU -std=c++17 -Wno-literal-suffix -lpthread -lrt -Wno-deprecated-declarations")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -std=legacy -O2  -pg  -g -ffixed-line-length-132 -lpthread -lrt -fdefault-real-8 -Wno-deprecated-declarations")
##    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ffree-form -O3 -g -lpthread -lrt -fdefault-real-8 -DTESTING -Wno-deprecated-declarations")
  endif()

  set(CMAKE_C_FLAGS_RELEASE        "-O3 -ffast-math -lpthread -fopenmp -lrt -march=native -Wno-deprecated-declarations")
  set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -ffast-math -DGNU  -std=c++17 -Wno-literal-suffix -lpthread -fopenmp -lrt -march=native -Wno-deprecated-declarations")
  set(CMAKE_Fortran_FLAGS_RELEASE  "-std=legacy -O3 -ffast-math -ffixed-line-length-132 -fopenmp -lpthread -lrt -fdefault-real-8 -march=native -Wno-deprecated-declarations")

  set (CMAKE_C_FLAGS_DEBUG       "-O0 -ggdb3 -lpthread -fopenmp -lrt -Wno-deprecated-declarations")
  set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -ggdb3  -DGNU  -std=c++17 -Wno-literal-suffix -lpthread -fopenmp -lrt  -Wno-deprecated-declarations")
  set (CMAKE_Fortran_FLAGS_DEBUG "-std=legacy -O0 -ggdb3 -lpthread -fopenmp -lrt -ffixed-line-length-132 -fdefault-real-8  -Wno-deprecated-declarations")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO       "${CMAKE_C_FLAGS_RELEASE}       -pg  -g")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELEASE}     -pg  -g")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -pg  -g")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  message(STATUS "Using Intel compiler flags")

  if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_C_FLAGS        "${CMAKE_C_FLAGS} -O2 -g -qmkl")
    set(CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS} -O2 -g -DINTEL -std=c++17 -qmkl")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -std=legacy -O2 -g -ffixed-line-length-132 -lpthread -lrt -fdefault-real-8 -Wno-deprecated-declarations")
    #set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -O2 -g -mkl -DTESTING -fbounds-check")
  endif()

  set(CMAKE_C_FLAGS_RELEASE        "-O3 -qmkl -march=native")
  set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -qmkl -DINTEL -std=c++17 -march=native")
  set(CMAKE_Fortran_FLAGS_RELEASE  "-O3 -qmkl -march=native -fbounds-check")

  set (CMAKE_C_FLAGS_DEBUG       "-O0 -g -qmkl -DTESTING -DDEBUGGING")
  set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -qmkl -DINTEL -DTESTING -std=c++17 -DDEBUGGING")
  set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -qmkl -DTESTING -DDEBUGGING")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO       "${CMAKE_C_FLAGS_RELEASE}       -g")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELEASE}     -g")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")

  add_definitions(-DUSE_MKL)
endif()

#set_property(SOURCE FVM_LocaInterface.C APPEND PROPERTY OBJECT_DEPENDS "/home/p266319/Documents/fredwubs/fvm/stochSWE/SWE_test/FVM_model_interface.h")

set(GEN_SOURCES
  general/Types.F90 
  general/Tools.F90
  general/FVM_Domain.C
  general/NOX_Epetra_LinearSystem_Belos.C
  general/GaleriExt_Utils.C
### general/OceanOutputXdmf.C needs Trios_Domain.H
  general/Filestreams.C)

# Sources shared by all DO targets (stochastic solver, basis solver)
# Stochastic solver core: basis solver, coefficient solver, linear solver back-end
set(DO_COMMON_SOURCES
  DO/StochSys.cpp
  DO/BasisSolver.cpp
  DO/LinearSolverWrapper.cpp
)

# Burgers and QG targets: direct Newton solver, no NOX/LOCA
set(DO_NONSWE_SOURCES
  DO/timedependent_do.cpp
  ${DO_COMMON_SOURCES}
)

# SWE target: NOX/LOCA implicit time stepper, adaptive step-size control
set(DO_SWE_SOURCES
  DO/timedependent_swe.cpp
  ${DO_COMMON_SOURCES}
)

set(BRGR_SOURCES
  Burger/Mean.cpp
  Burger/PDEAssembler.cpp)

set(QG_SOURCES
  QG/QG_Data.cpp
  QG/QG.cpp
  QG/Mean.cpp
  QG/PDEAssembler.cpp
)

set(SWE_SOURCES 
##  swe/Types.F90 
##  swe/Tools.F90
  swe/SWE_Interface.F90 
  swe/param.F90 
  swe/m_usr.F90 
  swe/m_mat.F90 
  swe/user.F90 
  swe/compleet.F90 
  swe/sphere.F90 
  swe/topo.F90 
  swe/windstress.F90 
  swe/splinefit.F90 
  swe/fitpack-dierckx.F90 
  swe/FVM_LocaInterface.C
  swe/ImplicitTimeStepper.cpp
  swe/ThetaStepperEvaluator.C
  swe/NOX_Epetra_LinearSystem_Amesos.cpp
)

include_directories(includes)
#find_package(MPI)
#set(PYTHON_INCLUDE_DIRS "/opt/intel/oneapi/intelpython/python3.9/include/python3.9")
#set(PYTHON_LIBRARIES "/usr/lib/x86_64-linux-gnu/libpython3.10.so")
#find_package(Python3 COMPONENTS Interpreter NumPy Development)
#PYHTON_INCLUDE_DIRS+="/usr/lib/python3/dist-packages/numpy/core/include/numpy/"
#include_directories(${PYTHON_INCLUDE_DIRS})
#link_directories(${PYTHON_LIBRARIES})
#message("Python_FOUND:${PYTHON_FOUND}")
#message("Python_VERSION:${PYTHON_VERSION}")
#message("Python_Development_FOUND:${PYTHON_DEVELOPMENT_FOUND}")
#message("Python_LIBRARIES:${PYTHON_LIBRARIES}")
find_package(MPI REQUIRED)
set(istrng "0")
if (DEFINED ENV{TRNG_HOME})
  include_directories($ENV{TRNG_HOME}/include)
  link_directories($ENV{TRNG_HOME}/lib)
  set(trng_lib "trng4")
  set(istrng "1")
endif()

if (DEFINED ENV{TRILINOS_HOME})
 set(Trilinos_DIR "$ENV{TRILINOS_HOME}/lib/cmake/Trilinos")
endif()
find_package(Trilinos REQUIRED)

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Trilinos_TPL_INCLUDE_DIRS})

link_directories(${Trilinos_LIBRARY_DIRS})
link_directories(${Trilinos_TPL_LIBRARY_DIRS})

message(STATUS "istrng=${istrng}")
message(STATUS "trilinos_LIBRARIES=${Trilinos_LIBRARIES}")
#string(REPLACE "${Trilinos_LIBRARIES}" "${Trilinos_LIBRARIES};trilinosss" tril_libs "${Trilinos_LIBRARIES}")
#string(REPLACE "pytrilinos;" "" tril_libs "${Trilinos_LIBRARIES}")
#set(Trilinos_LIBRARIES ${tril_libs})
#####message(STATUS "trilinos_TPL_LIBRARIES=${Trilinos_TPL_LIBRARIES}")

#message(STATUS "trilinos_LIBRARIES : ${Trilinos_LIBRARIES}")

link_directories("$ENV{TRILINOS_HOME}/lib")
#####include_directories("$ENV{HYMLS_HOME}/include")
#####link_directories("$ENV{HYMLS_HOME}/lib")
if (DEFINED ENV{HYMLS_HOME})
  set(HYMLS_DIR "$ENV{HYMLS_HOME}/lib/cmake/HYMLS")
endif()
find_package(HYMLS REQUIRED)
######target_link_libraries(SWEtimedep "libhymls.a" -lgfortran )
include_directories(${HYMLS_INCLUDE_DIRS})

add_library(gen SHARED ${GEN_SOURCES})
target_link_libraries(gen ${HYMLS_LIBRARIES})
target_link_libraries(gen ${Trilinos_TPL_LIBRARIES})
target_link_libraries(gen ${Trilinos_LIBRARIES})

add_library(BRGR SHARED ${BRGR_SOURCES} )
target_link_libraries(BRGR gen)

list(FIND Trilinos_TPL_LIST SCALAPACK SCALAPACK_found)

if(SCALAPACK_found GREATER -1)
  add_executable(BRGRDO ${DO_NONSWE_SOURCES})
  target_compile_definitions(BRGRDO PRIVATE brgr=1 quasi_geo=0 need_locaInterface=0 use_trng=${istrng})
  target_link_libraries(BRGRDO BRGR)
  target_link_libraries(BRGRDO ${HYMLS_LIBRARIES})
  target_link_libraries(BRGRDO ${Trilinos_TPL_LIBRARIES})
  target_link_libraries(BRGRDO ${Trilinos_LIBRARIES})
  target_link_libraries(BRGRDO ${PYTHON_LIBRARIES})
  target_link_libraries(BRGRDO ${trng_lib})
  target_link_libraries(BRGRDO gen)
endif()

add_library(QG SHARED ${QG_SOURCES} )
target_link_libraries(QG gen)

if(SCALAPACK_found GREATER -1)
  add_executable(QGDO ${DO_NONSWE_SOURCES})
  target_compile_definitions(QGDO PRIVATE brgr=0 quasi_geo=1 need_locaInterface=0 use_trng=${istrng})
  target_link_libraries(QGDO QG)
  target_link_libraries(QGDO ${HYMLS_LIBRARIES})
  target_link_libraries(QGDO ${Trilinos_TPL_LIBRARIES})
  target_link_libraries(QGDO ${Trilinos_LIBRARIES})
  target_link_libraries(QGDO ${PYTHON_LIBRARIES})
  target_link_libraries(QGDO ${trng_lib})
  target_link_libraries(QGDO gen)
endif()

add_library(SWE SHARED ${SWE_SOURCES} )
set_target_properties (SWE PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/swe
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/general
  )
target_link_libraries(SWE gen)

if(SCALAPACK_found GREATER -1)
  add_executable(SWEDO ${DO_SWE_SOURCES})
  target_compile_definitions(SWEDO PRIVATE brgr=0 quasi_geo=0 need_locaInterface=1 use_trng=${istrng})
  target_link_libraries(SWEDO SWE)
  target_link_libraries(SWEDO ${HYMLS_LIBRARIES})
  target_link_libraries(SWEDO ${Trilinos_TPL_LIBRARIES})
  target_link_libraries(SWEDO ${Trilinos_LIBRARIES})
  target_link_libraries(SWEDO "libtrilinosss.so" -lgfortran)
  target_link_libraries(SWEDO ${PYTHON_LIBRARIES})
  #target_link_libraries(SWEDO ${trng_lib})
  target_link_libraries(SWEDO gen)
endif()

find_package(Doxygen)
#set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${PROJECT_SOURCE_DIR}/README.md")
if(DOXYGEN_FOUND)
    set(ORIGINAL_README "${PROJECT_SOURCE_DIR}/README.md")
    set(DOXYGEN_README "${CMAKE_CURRENT_BINARY_DIR}/README_doxygen.md")

    if(EXISTS "${ORIGINAL_README}")
        file(READ "${ORIGINAL_README}" README_CONTENTS)
        
        # Pass 1: Translate GitHub block math ($$...$$) to Doxygen block math (\f[...\f])
        string(REGEX REPLACE "\\$\\$([^$]+)\\$\\$" "\\\\f[\\1\\\\f]" README_CONTENTS "${README_CONTENTS}")
        
        # Pass 2: Secure standard currency syntax (e.g., $10) so they don't break
        # This prevents isolated numeric dollar instances from turning into math blocks
        string(REGEX REPLACE " \\$([0-9])" " \\\\$\\1" README_CONTENTS "${README_CONTENTS}")
        
        # Pass 3: Translate GitHub inline math ($...$) to Doxygen inline math (\f$... \f$)
        # The [^$]+ rule ensures it matches only true enclosed inline segments
        string(REGEX REPLACE "\\$([^$ \n][^$\n]*[^$ \n])\\$" "\\\\f$\\1\\\\f$" README_CONTENTS "${README_CONTENTS}")
        
        file(WRITE "${DOXYGEN_README}" "${README_CONTENTS}")
    else()
        set(DOXYGEN_README "${ORIGINAL_README}") 
    endif()

    # Pass the converted file to your Doxygen processor
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${DOXYGEN_README}")

    # Set Doxygen configuration variables
    set(DOXYGEN_PROJECT_NAME "DO Project")
    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_RECURSIVE YES)
    set(DOXYGEN_USE_MATHJAX YES)
    set(DOXYGEN_MATHJAX_FORMAT "HTML-CSS") 
    # Crucial for your Fortran source files (swe/*.F90)
    set(DOXYGEN_OPTIMIZE_FOR_FORTRAN YES)
    set(DOXYGEN_EXTRACT_ALL YES)
    set(DOXYGEN_EXTRACT_PRIVATE YES)
    set(DOXYGEN_EXTRACT_STATIC YES)

    # Tell Doxygen where to look for source code. 
    # This covers your general, DO, Burger, QG, and swe directories.
    set(DOXYGEN_SOURCE_DIR 
        ${PROJECT_SOURCE_DIR}/general
        ${PROJECT_SOURCE_DIR}/DO
        ${PROJECT_SOURCE_DIR}/Burger
        ${PROJECT_SOURCE_DIR}/QG
        ${PROJECT_SOURCE_DIR}/swe
        ${PROJECT_SOURCE_DIR}/includes
    )

    # Create the 'docs' target
    doxygen_add_docs(docs 
        ${DOXYGEN_SOURCE_DIR}
        "${DOXYGEN_README}"
        #"${PROJECT_SOURCE_DIR}/README.md"
        COMMENT "Generating API documentation with Doxygen"
    )
    
    message(STATUS "Doxygen target 'docs' has been created.")
else()
    message(WARNING "Doxygen not found - documentation target will not be created.")
endif()

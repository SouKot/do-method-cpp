cmake_minimum_required(VERSION 3.8)
project(DO C CXX Fortran)
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  message(STATUS "Using GNU compiler flags")

  if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	  #set(CMAKE_C_COMPILER "scorep mpicc")
	  #set(CMAKE_CXX_COMPILER "scorep mpic++")
	  #set(CMAKE_Fortran_COMPILER "scorep gfortran")
    set(CMAKE_C_FLAGS        "${CMAKE_C_FLAGS} -O2  -pg  -g -lpthread -lrt  -Wno-deprecated-declarations")
    set(CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS} -O2  -pg  -g -DGNU  -std=c++17 -Wno-literal-suffix -lpthread -lrt -Wno-deprecated-declarations")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -std=legacy -O2  -pg  -g -ffixed-line-length-132 -lpthread -lrt -fdefault-real-8 -Wno-deprecated-declarations")
##    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ffree-form -O3 -g -lpthread -lrt -fdefault-real-8 -DTESTING -Wno-deprecated-declarations")
  endif()

  set(CMAKE_C_FLAGS_RELEASE        "-O3 -ffast-math -lpthread -lrt -march=native -Wno-deprecated-declarations")
  set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -DGNU -ffast-math -std=c++17 -Wno-literal-suffix -lpthread -lrt -march=native -Wno-deprecated-declarations")
  set(CMAKE_Fortran_FLAGS_RELEASE  "-std=legacy -O3 -ffixed-line-length-132 -ffast-math -lpthread -lrt -fdefault-real-8 -march=native -Wno-deprecated-declarations")

  set (CMAKE_C_FLAGS_DEBUG       "-O0 -ggdb3 -lpthread -lrt -Wno-deprecated-declarations")
  set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -ggdb3  -DGNU  -std=c++17 -Wno-literal-suffix -lpthread -lrt  -Wno-deprecated-declarations")
  set (CMAKE_Fortran_FLAGS_DEBUG "-std=legacy -O0 -ggdb3 -lpthread -lrt -ffixed-line-length-132 -fdefault-real-8  -Wno-deprecated-declarations")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO       "${CMAKE_C_FLAGS_RELEASE}       -pg  -g")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELEASE}     -pg  -g")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -pg  -g")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  message(STATUS "Using Intel compiler flags")

  if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_C_FLAGS        "${CMAKE_C_FLAGS} -O2 -g -mkl -DTESTING")
    set(CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS} -O2 -g -std=c++11 -mkl -DTESTING")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -O2 -g -mkl -DTESTING")
  endif()

  set(CMAKE_C_FLAGS_RELEASE        "-O3 -mkl -ffast-math -march=native")
  set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -mkl -ffast-math -std=c++11 -march=native")
  set(CMAKE_Fortran_FLAGS_RELEASE  "-O3 -mkl -ffast-math -march=native")

  set (CMAKE_C_FLAGS_DEBUG       "-O0 -g -mkl -DTESTING -DDEBUGGING")
  set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -mkl -DTESTING -std=c++11 -DDEBUGGING")
  set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -mkl -DTESTING -DDEBUGGING")

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

set(DO_SOURCES 
  DO/timedependent.cpp 
  DO/StochSys.cpp
  DO/Problem_Interface.cpp
 )

set(BRGR_SOURCES
  Burger/Mean.cpp)

set(QG_SOURCES
  QG/QG_Data.cpp
  QG/QG.cpp
  QG/Mean.cpp
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

#find_package(PythonLibs 2.7)
#include_directories(${PYTHON_INCLUDE_DIRS})
#link_directories(${PYTHON_LIBRARIES})

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
find_package(Trilinos)

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Trilinos_TPL_INCLUDE_DIRS})

link_directories(${Trilinos_LIBRARY_DIRS})
link_directories(${Trilinos_TPL_LIBRARY_DIRS})

message(STATUS "istrng=${istrng}")
message(STATUS "trilinos_LIBRARIES=${Trilinos_LIBRARIES}")
string(REPLACE "teko;" "" tril_libs "${Trilinos_LIBRARIES}")
set(Trilinos_LIBRARIES ${tril_libs})
message(STATUS "trilinos_LIBRARIES : ${Trilinos_LIBRARIES}")
#####include_directories("$ENV{HYMLS_HOME}/include")
#####link_directories("$ENV{HYMLS_HOME}/lib")
if (DEFINED ENV{HYMLS_HOME})
  set(HYMLS_DIR "$ENV{HYMLS_HOME}/lib/cmake/HYMLS")
endif()
find_package(HYMLS)

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
  add_executable(BRGRDO ${DO_SOURCES})
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
  add_executable(QGDO ${DO_SOURCES})
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
  add_executable(SWEDO ${DO_SOURCES})
  target_compile_definitions(SWEDO PRIVATE brgr=0 quasi_geo=0 need_locaInterface=1 use_trng=${istrng})
  target_link_libraries(SWEDO SWE)
  target_link_libraries(SWEDO ${HYMLS_LIBRARIES})
  target_link_libraries(SWEDO ${Trilinos_TPL_LIBRARIES})
  target_link_libraries(SWEDO ${Trilinos_LIBRARIES})
  target_link_libraries(SWEDO ${PYTHON_LIBRARIES})
  target_link_libraries(SWEDO ${trng_lib})
  target_link_libraries(SWEDO gen)
endif()


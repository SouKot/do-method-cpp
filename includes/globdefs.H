/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef GLOBDEFS_H
#define GLOBDEFS_H

#ifdef TRIOS_MACROS_H
#error                                                                         \
  "you should include trilinos_thcm's globdefs.H before trios' TRIOS_Macros.H"
#endif

#include "Teuchos_RCP.hpp"
#include "Teuchos_Time.hpp"
#include <cmath>
#include <iostream>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "my_f2c.H"

// some global definitions for the C++ part of this project
// They should be consistent with usr.F90 (formerly usr.com)

// pi
#define PI_ 3.14159265358979323846

//! these streams are linked to globally unique filenames
//! so that output from different processes is placed into
//! different files:

#ifdef DEBUGGING
#ifndef TESTING
#define TESTING
#endif
extern Teuchos::RCP<std::ostream> debug;
#ifndef DEBUG
#define DEBUG(s) (*debug) << s << std::endl;
#define DEBVAR(s) (*debug) << #s << " = " << s << std::endl;
#endif
#ifndef INFO
#define INFO(s)                                                                \
  {                                                                            \
    (*info) << s << std::endl;                                                 \
    (*debug) << s << std::endl;                                                \
  }
#endif
#else
#ifndef INFO
#define INFO(s) (*info) << s << std::endl;
#endif
#ifndef DEBUG
#define DEBUG(s)
#define DEBVAR(s)
#endif
#endif

extern Teuchos::RCP<std::ostream> info;
extern Teuchos::RCP<Teuchos::Time> coeffTime;
extern Teuchos::RCP<Teuchos::Time> basisTime;
extern Teuchos::RCP<Teuchos::Time> meanTime;
extern Teuchos::RCP<Teuchos::Time> meanSolverTime;
extern Teuchos::RCP<Teuchos::Time> meanMassMatTime;
extern Teuchos::RCP<Teuchos::Time> meanRhsTime;
extern Teuchos::RCP<Teuchos::Time> meanJacTime;
// macros for checking Epetra calls
#ifndef CHECK_ZERO
#define CHECK_ZERO(funcall)                                                    \
  {                                                                            \
    int ierr = funcall;                                                        \
    if (ierr) {                                                                \
      INFO("Trilinos Error " << ierr << " returned from call " << #funcall);   \
    }                                                                          \
  }
#define CHECK_TRUE(funcall)                                                    \
  {                                                                            \
    bool berr = funcall;                                                       \
    if (!berr) {                                                               \
      INFO("Trilinos call " << #funcall << " returned false");                 \
    }                                                                          \
  }
#define CHECK_NONNEG(funcall)                                                  \
  {                                                                            \
    int ierr = funcall;                                                        \
    if (ierr < 0) {                                                            \
      INFO("Trilinos Error " << ierr << " returned from call " << #funcall);   \
    }                                                                          \
  }
#endif

//! this function prints an error message and terminates the program.
//! It is implemented in Filestreams.C
void
Error(std::string msg, std::string file, int line);

#ifndef MOD
//! our own 'modulo' function, which behaves like mod in matlab.
//! the C++ built-in '%' operator returns -1%n=-1 and is therefore
//! not very useful for periodic boundaries...
#define MOD(x, y)                                                              \
  (((double)(y) == 0.0)                                                        \
     ? (double)(x)                                                             \
     : ((double)(x)-floor((double)(x) / ((double)(y))) * ((double)(y))))
#endif

#endif

/*
 * =====================================================================================
 *
 *       Filename:  Problem_Interface.hpp
 *
 *    Description:  Backward-compatibility shim for the DO basis solver.
 *
 *                  The class was renamed to BasisSolver in Phase 2d of the
 *                  refactoring.  This header provides an alias so that all
 *                  existing call sites continue to compile without modification.
 *
 *                  New code should include BasisSolver.hpp directly.
 *
 *        Version:  2.0  (shim; formerly the full Problem_Interface declaration)
 *        Created:  2026-05-13
 *       Revision:  none
 *       Compiler:  gcc / clang (C++17)
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:
 *
 * =====================================================================================
 */
#ifndef Problem_Interface_H
#define Problem_Interface_H

#include "BasisSolver.hpp"

/// Alias kept for backward compatibility — prefer BasisSolver in new code.
using Problem_Interface = BasisSolver;

#endif // Problem_Interface_H

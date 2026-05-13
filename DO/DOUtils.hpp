/*
 * =====================================================================================
 *
 *       Filename:  DOUtils.hpp
 *
 *    Description:  Lightweight shared utilities for all DO targets (BRGRDO, QGDO,
 *                  SWEDO).  Consolidates the identical printnormMV helpers that
 *                  previously existed separately in Problem_Interface, Y_Stoch,
 *                  and QG::Mean into a single free function.
 *
 *        Version:  1.0
 *        Created:  2026-05-06
 *       Revision:  none
 *       Compiler:  gcc / clang (C++17)
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:
 *
 * =====================================================================================
 */
#ifndef DO_UTILS_HPP
#define DO_UTILS_HPP

#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include <iostream>
#include <string>
#include <vector>

/**
 * @brief Print the 1- or 2-norm of each column of an Epetra_MultiVector.
 *
 * Output is written to stdout only from MPI rank 0 to avoid duplicate lines
 * in parallel runs.  Replaces the four identical per-class implementations
 * that previously existed in Problem_Interface, Y_Stoch, and QG::Mean.
 *
 * @param[in,out] mv       The multi-vector whose column norms are printed.
 * @param[in]     normType 1 for the 1-norm, 2 for the 2-norm.
 * @param[in]     str      Label string printed before the norm values.
 */
inline void printnormMV(Epetra_MultiVector& mv, int normType, const std::string& str)
{
    std::vector<double> nrm(mv.NumVectors());
    if (normType == 1)
        mv.Norm1(nrm.data());
    else if (normType == 2)
        mv.Norm2(nrm.data());
    if (mv.Comm().MyPID() == 0) {
        std::cout << "\n" << str << std::endl;
        for (int i = 0; i < mv.NumVectors(); i++)
            std::cout << "  " << nrm[i];
    }
}

#endif // DO_UTILS_HPP

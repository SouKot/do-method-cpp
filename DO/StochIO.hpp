/*
 * =====================================================================================
 *
 *       Filename:  StochIO.hpp
 *
 *    Description:  Centralised MatrixMarket I/O helpers for DO stochastic solvers.
 *                  Wraps EpetraExt::MultiVectorToMatrixMarketFile,
 *                  EpetraExt::RowMatrixToMatrixMarketFile, and
 *                  EpetraExt::MatrixMarketFileToMultiVector behind a thin
 *                  namespace so that callers never see raw EpetraExt returns or
 *                  raw-pointer ownership rules.
 *
 *                  All write helpers throw std::runtime_error on I/O failure.
 *                  readMV returns an RCP<Epetra_MultiVector> that owns the object
 *                  so callers need no manual delete.
 *
 *        Version:  1.0
 *        Created:  2026-05-13
 *       Revision:  none
 *       Compiler:  gcc / clang (C++17)
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:
 *
 * =====================================================================================
 */
#ifndef DO_STOCH_IO_HPP
#define DO_STOCH_IO_HPP

#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include <string>

// =====================================================================================
/// @name StochIO — MatrixMarket I/O for DO stochastic data
// =====================================================================================
/// @{

/**
 * @brief Centralised MatrixMarket I/O namespace for the DO method.
 *
 * All functions are inline; the namespace adds no object files.  Include this
 * header wherever DO multi-vector or matrix data is read or written.
 */
namespace StochIO {

/**
 * @brief Write an Epetra_MultiVector to a MatrixMarket file.
 *
 * @param fname Destination file path (created or overwritten).
 * @param mv    Multi-vector to serialise.
 * @throws std::runtime_error if EpetraExt reports an error.
 */
inline void writeMV(const std::string& fname, const Epetra_MultiVector& mv)
{
    int err = EpetraExt::MultiVectorToMatrixMarketFile(fname.c_str(), mv);
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error,
        "StochIO::writeMV: failed to write \"" + fname + "\"");
}

/**
 * @brief Write an Epetra_RowMatrix to a MatrixMarket file.
 *
 * @param fname Destination file path (created or overwritten).
 * @param mat   Row matrix to serialise.
 * @throws std::runtime_error if EpetraExt reports an error.
 */
inline void writeMatrix(const std::string& fname, const Epetra_RowMatrix& mat)
{
    int err = EpetraExt::RowMatrixToMatrixMarketFile(fname.c_str(), mat);
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error,
        "StochIO::writeMatrix: failed to write \"" + fname + "\"");
}

/**
 * @brief Read a MatrixMarket file into a new, RCP-owned Epetra_MultiVector.
 *
 * Replaces the raw-pointer EpetraExt pattern:
 * @code
 *   Epetra_MultiVector* raw = nullptr;
 *   EpetraExt::MatrixMarketFileToMultiVector(fname.c_str(), map, raw);
 *   // ... use raw ...
 *   delete raw;
 * @endcode
 * with the safer:
 * @code
 *   auto mv = StochIO::readMV(fname, someObj.Map());
 *   // ... use mv as RCP<Epetra_MultiVector> ...
 * @endcode
 *
 * @param fname    Source file path.
 * @param blockmap Epetra map describing the distribution of the returned vector.
 *                 Accepts Epetra_BlockMap so callers can pass the result of
 *                 Epetra_MultiVector::Map() or Epetra_Vector::Map() directly.
 *                 The underlying object must be an Epetra_Map (all Epetra
 *                 distributed objects that hold row data use Epetra_Map).
 * @return RCP owning the newly allocated Epetra_MultiVector.
 * @throws std::runtime_error if the file cannot be read.
 */
inline Teuchos::RCP<Epetra_MultiVector>
readMV(const std::string& fname, const Epetra_BlockMap& blockmap)
{
    // EpetraExt takes Epetra_Map; the actual object behind any Epetra
    // distributed vector map is always an Epetra_Map, so the cast is safe.
    const Epetra_Map& map = static_cast<const Epetra_Map&>(blockmap);
    Epetra_MultiVector* raw = nullptr;
    int err = EpetraExt::MatrixMarketFileToMultiVector(fname.c_str(), map, raw);
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0 || raw == nullptr, std::runtime_error,
        "StochIO::readMV: failed to read \"" + fname + "\"");
    return Teuchos::rcp(raw);
}

} // namespace StochIO

/// @}

#endif // DO_STOCH_IO_HPP

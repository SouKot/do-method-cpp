/*
 * =====================================================================================
 *
 *       Filename:  CoeffSolver.hpp
 *
 *    Description:  class for solving stochastic system
 *
 *        Version:  1.0
 *        Created:  01/15/2016 04:12:43 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */

#ifndef DO_COEFF_SOLVER_HPP
#define DO_COEFF_SOLVER_HPP
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Vector.h"
#include "Ifpack_Preconditioner.h"
#include <BelosEpetraAdapter.hpp>
#include <BelosIMGSOrthoManager.hpp>
#include <BelosMultiVec.hpp>
#include <Epetra_SerialDenseSVD.h>
#include <Epetra_SerialDenseSolver.h>
#include <Teuchos_SerialDenseMatrix.hpp>
#include "NoiseGenerator.hpp"
#include "StochasticState.hpp"
// #include "Sacado.hpp"
// #include "Sacado_Fad_BLAS.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_RCP.hpp"
#include "globdefs.H"
#include "DOUtils.hpp"
#if need_locaInterface == 1
#include "FVM_Domain.H"
using DomainPtr = Teuchos::RCP<FVM::Domain>;
#else /* -----  not NEED_LOCAINTERFACE  ----- */
#include "Interface.hpp"
/// @brief Placeholder type used when FVM::Domain is not available (non-LOCA builds).
struct DomainPlaceholder {};
using DomainPtr = Teuchos::RCP<DomainPlaceholder>;
#endif /* -----  not NEED_LOCAINTERFACE  ----- */
// Include header to define eigenproblem Ax = \lambda*x
#include "AnasaziBasicEigenproblem.hpp"
// Include header to provide Anasazi with Epetra adapters.  If you
// plan to use Tpetra objects instead of Epetra objects, include
// AnasaziTpetraAdapter.hpp instead; do analogously if you plan to use
// Thyra objects instead of Epetra objects.
#include "AnasaziEpetraAdapter.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
/**
 * @brief Stochastic coefficient solver for the Dynamically Orthogonal (DO) method.
 *
 * Solves the evolution equation for the stochastic coefficient matrix
 * @f$ Y(t) \in \mathbb{R}^{m \times S} @f$, where @e m is the number of
 * DO basis modes and @e S is the number of stochastic realisations.
 *
 * The DO decomposition represents a stochastic field as
 * @f$ u(x,t,\omega) = \bar{u}(x,t) + V(x,t)\,Y(t,\omega) @f$,
 * where @f$ \bar{u} @f$ is the mean (solved by the Mean class) and @f$ V @f$
 * is the orthonormal basis (solved by BasisSolver).  This class handles:
 *   - Assembly of the coefficient RHS and Jacobian from bilinear/linear terms.
 *   - A Newton iteration (with optional backtracking) to advance Y one time step.
 *   - Computation of second-moment statistics @f$ E[yy^T] @f$ and
 *     @f$ E[V y V^T y y^T] @f$ needed by the basis and mean solvers.
 *   - Generation of Wiener increments via NoiseGenerator.
 *
 * @see BasisSolver  The companion class solving for the spatial basis V.
 * @see StochasticState  The shared data contract coupling CoeffSolver and BasisSolver.
 *
 * @class CoeffSolver CoeffSolver.hpp "CoeffSolver.hpp"
 */
class CoeffSolver
{
public:
  /**
   * @brief Construct the stochastic coefficient solver.
   *
   * @param NumStochIter  Number of stochastic realisations.
   * @param num_Subtime_Step  Sub-time-steps per outer step (subdt = dt / num_Subtime_Step).
   * @param m  Number of DO basis modes.
   * @param dt  Pointer to the current time-step size.
   * @param A  Jacobian of the deterministic RHS.
   * @param uav  Mean-field solution vector.
   * @param domain  FVM domain (LOCA builds) or Teuchos::null (non-LOCA builds).
   * @param Vn  Stochastic basis multi-vector (N × m).
   * @param Wb  Stochastic forcing multi-vector.
   * @param comm  Epetra communicator.
   * @param CoefParams  Parameter list for this class.
   * @param sharedState  Shared state object visible to BasisSolver.
   * @param maxNumIter  Maximum Newton iterations.
   * @param useBacktracking  Enable backtracking line-search.
   * @param numBackTrackingSteps  Number of backtracking steps.
   * @param toleranceRHS  Newton convergence tolerance.
   * @param NormRHS  Initial RHS norm.
   */
  CoeffSolver(int NumStochIter,
          int num_Subtime_Step,
          int m,
          double* dt,
          const Teuchos::RCP<Epetra_CrsMatrix>& A,
          const Teuchos::RCP<Epetra_Vector>& uav,
          const DomainPtr& domain,
          const Teuchos::RCP<Epetra_MultiVector>& Vn,
          const Teuchos::RCP<Epetra_MultiVector>& Wb,
          const Teuchos::RCP<Epetra_Comm>& comm,
          const Teuchos::RCP<Teuchos::ParameterList>& CoefParams,
          const Teuchos::RCP<StochasticState>& sharedState,
          int maxNumIter,
          bool useBacktracking,
          double numBackTrackingSteps,
          double toleranceRHS,
          double NormRHS);
  //** \name Overridden from EpetraExt::ModelEvaluator . */
  Teuchos::RCP<Epetra_Map> get_x_map();
  Teuchos::RCP<Epetra_Map> get_f_map();
  Teuchos::RCP<Epetra_MultiVector> get_x_init();
  void set_x(const Teuchos::RCP<Epetra_MultiVector>& x0_temp);
  /*!
   \brief set dW (wiener process)  for the stochastic iteartions.
   \fn setDwiener
   \param z0_temp : pointer to stochastic random gaussion matrix of size
           NumStochForcingVecotr*NumStochiter.

    * NOTE: We set dW=z0_temp and then mutiply dW with \f$\sqrt{dt}$\f.
    * It implies that the address of variable(in our case, z_) which is
    * paased on to zo_temp also gets multiplied by \f$\sqrt{dt}\f.
  */
  void setDwiener(const Teuchos::RCP<Epetra_MultiVector>& z0_temp);
  Teuchos::RCP<Epetra_Vector> getEVyVy();
  /**
   * @brief Evaluate the bilinear form B(ud, u_i) for each column of u.
   *
   * For the non-LOCA path this delegates to the StochasticState::bilinearTerm
   * callback set by the driver; for the LOCA path it calls the Fortran model_bil.
   *
   * @param ud  Mean-field vector.
   * @param u   Multi-vector of basis or coefficient columns.
   * @param v   Second argument multi-vector.
   * @param uv  Result multi-vector, overwritten with B(u_i, v_i).
   */
  void Bilinear(const Teuchos::RCP<Epetra_Vector>& ud,
                const Teuchos::RCP<Epetra_MultiVector>& u,
                const Teuchos::RCP<Epetra_MultiVector>& v,
                const Teuchos::RCP<Epetra_MultiVector>& uv);
  void HBilinV();
  // void BilinTerm();
  void jacBilinTerm();
  void NonLinRHS();
  void LinCoeff();
  void computeRepEVyVy();
  void y_rhs();
  void y_jac();
  void Solve();
  void StochasticIterations();
  void computeExpectations();
  int BlockDavidsonMethod();
  void PostProcess(Epetra_Vector& second_mmnt);
  void Newton();
  void RunBackTracking();
  void CreateLocMultiVec(const std::string& WhichOne);
  void CreateDistTransMultivec();
  bool Converged() { return isConverged_; }
  int Iterations() { return iter_; }
  Teuchos::RCP<Epetra_SerialDenseMatrix> getJacobian() { return Yjac; }
  Teuchos::RCP<Epetra_MultiVector> getEyy()         { return Eyy; }
  Teuchos::RCP<Epetra_MultiVector> getEDEyy()   { return EDEyy; }
  Teuchos::RCP<Epetra_MultiVector> getY()           { return y_; }
  Epetra_MultiVector&              getYTrans()       { return *yTrans_; }

  int numSubTimeStep, MyLDA;
  int m_, rszyy;
  int stochiter;
  int N_, MyPID, itrtr = 0;
  double subdt_;
  double* dt_;
  Teuchos::RCP<Epetra_CrsMatrix> A_, EYY;
  // Teuchos::RCP<Epetra_CrsMatrix> mass_;
  Teuchos::RCP<Epetra_MultiVector> Vnew, fortVn, fortviv,
    viv; //*********************************************FILL IT With the
         // value of bases
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*****************
  Teuchos::RCP<Epetra_Vector> udet,fortudet;
  Teuchos::RCP<Epetra_MultiVector> AV;
  Teuchos::RCP<Epetra_MultiVector> VAV;
  Teuchos::RCP<Epetra_MultiVector> udet_mtimes;
  Teuchos::RCP<Epetra_Vector> ones;
  Teuchos::RCP<Epetra_MultiVector> Vudet, udetV;
  Teuchos::RCP<Epetra_MultiVector> VVudet;
  Teuchos::RCP<Epetra_MultiVector> lin_coeff, JacNonLin;
  Teuchos::RCP<Epetra_MultiVector> repEVyVy;
  Teuchos::RCP<Epetra_MultiVector> Vtemp;
  Teuchos::RCP<Epetra_Vector>
    EVyVy; // map of udet from deterministic part( n*1 length vector )
  Teuchos::RCP<Epetra_MultiVector> const_coeff;
  Teuchos::RCP<Epetra_MultiVector> rhs;
  Teuchos::RCP<Epetra_MultiVector> B;
  Teuchos::RCP<Epetra_MultiVector> VB;
  Teuchos::RCP<Epetra_MultiVector> dW;
  Teuchos::RCP<Epetra_MultiVector> VBdW;
  Teuchos::RCP<Epetra_Vector> xndiff;
  Teuchos::RCP<Epetra_SerialDenseMatrix> Yjac;
  Teuchos::RCP<Epetra_MultiVector> jac;
  Epetra_SerialDenseMatrix DenseDx_, DenseRHS;
  // Teuchos::RCP<Epetra_SerialDenseVector> DenseRHS;
  Teuchos::RCP<Epetra_MultiVector> z_, y_, rhsNonLin, yTrans_, zTrans_,
    EVyVyyDEyy, EDEyy, eye;
  Epetra_SerialDenseSVD y_prob;
  Teuchos::RCP<Epetra_MultiVector> H, Hn, VHn, Rv, YY;
  Teuchos::RCP<Epetra_Map> Tmap, map_mm;
  Teuchos::RCP<Epetra_MultiVector> Ezy, EzyPrev, EzyDEyy, EVyVyy;
  Teuchos::RCP<Epetra_MultiVector> Eyy, LocEyyy, GlobEyyy,
    EyyTyT, Utmp;
  Epetra_SerialDenseMatrix ru;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>> Ru;

  // Epetra_MultiVector* viv;
  // Teuchos::RCP<Epetra_Vector>Teuchos::RCP<Epetra_CrsMatrix> W_jac;
  // Teuchos::RCP<Epetra_MultiVector> y_jac_ ;
  // void model_bil(double* u, double* v, double* uv){};
  // typedef Sacado::Fad::DFad<double> FadType;
  // Teuchos::RCP< Epetra_Vector> x_old ;
  Teuchos::RCP<Epetra_MultiVector> f_out;
  Teuchos::RCP<Epetra_MultiVector> x0_, x_, dx_;
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP>> problem;
  std::vector<Anasazi::Value<double>> evals;
  Teuchos::RCP<MV> evecs;
  Anasazi::Eigensolution<double, MV> EigSol;

  void printTransNormMV(Epetra_MultiVector& mv, int normType, const std::string& str);
  void computeEVyVy();
  void computeCrossVariance();
  void computeEDEyy();
  void computeEVyVyy();
  void computeEyyTyT();
  void SymMatPseudoInverse(Epetra_MultiVector& mat, Epetra_MultiVector& matInv);

  void test_jac_billin();
  
private:
  // /////////////////////////////////////
  // Private member data
  NoiseGenerator noiseGen_;
  int iter_;
  int maxNumIterations_;
  int backTrack_;
  int numBackTrackingSteps_;
  double toleranceRHS_;
  double NormRHS_;
  double NormRHStest_;
  int stochiter_;
  int NumStochIter_;
  //!
  bool test_, debug_, useNwtn_, shouldEquil;
  bool isConverged_;
  bool backTracking_; // perhaps call this enableBacktracking_
  Teuchos::RCP<Epetra_Comm> Comm_;
  Teuchos::RCP<Epetra_LocalMap> map_z_, map_x_;
  Teuchos::RCP<Epetra_MultiVector> EyyOld_;
  Teuchos::RCP<StochasticState> sharedState_;
  DomainPtr domain_;
};
#endif

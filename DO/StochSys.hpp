/*
 * =====================================================================================
 *
 *       Filename:  StochSys.hpp
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

#ifndef do_method_h
#define do_method_h
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
#else /* -----  not NEED_LOCAINTERFACE  ----- */
#include "Interface.hpp"
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
/*!
 \brief Class for calculating stochastic coefficients.

 \class Y_Stoch StochSys.hpp "StochSys.hpp"
*/
class Y_Stoch
{
public:
#if need_locaInterface == 1
  Y_Stoch(int NumStochIter,
          int num_Subtime_Step,
          int m,
          double* dt,
          const Teuchos::RCP<Epetra_CrsMatrix>& A,
          const Teuchos::RCP<Epetra_Vector>& uav,
          const Teuchos::RCP<FVM::Domain>& domain,
          const Teuchos::RCP<Epetra_MultiVector>& Vn,
          const Teuchos::RCP<Epetra_MultiVector>& Wb,
          const Teuchos::RCP<Epetra_Comm>& comm,
          const Teuchos::RCP<Teuchos::ParameterList>& CoefParams,
          int maxNumIter,
          bool useBacktracking,
          double numBackTrackingSteps,
          double toleranceRHS,
          double NormRHS);

#else
  /*!
   \brief constructor with no links to LocaInterface

   \fn Y_Stoch
   \param NumStochIter : number of stochastic iteration
   \param num_Subtime_Step : subdt= dt/num_Subtime_step
   \param m : is the number of stochastic basis
   \param dt :
   \param A : is jacobian of RHS of deterministic system's (i.e., system without
   stochastics) \param uav \param Vn : stocahstic basis (of size N*m) \param Wb
   : stochastic forcing \param comm : \param CoefParams : a parameterlist
   containing parameters for this class \param model : Mean class object \param
   maxNumIter : maximum iterations newton solver \param useBacktracking :
   whether to use backtracking or not \param numBackTrackingSteps \param
   toleranceRHS : tolerance for newton solver \param NormRHS
  */
  Y_Stoch(int NumStochIter,
          int num_Subtime_Step,
          int m,
          double* dt,
          const Teuchos::RCP<Epetra_CrsMatrix>& A,
          const Teuchos::RCP<Epetra_Vector>& uav,
          const Teuchos::RCP<Epetra_MultiVector>& Vn,
          const Teuchos::RCP<Epetra_MultiVector>& Wb,
          const Teuchos::RCP<Epetra_Comm>& comm,
          const Teuchos::RCP<Teuchos::ParameterList>& CoefParams,
          const Teuchos::RCP<Mean>& model,
          int maxNumIter,
          bool useBacktracking,
          double numBackTrackingSteps,
          double toleranceRHS,
          double NormRHS);
#endif /* -----  not NEED_LOCAINTERFACE  ----- */
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
   * ===  FUNCTION
   *======================================================================
   *         Name:  Bilinear
   *  Description:  Bilinear form should be defined here!!! This function should
   *  be called for computing a bilinear term!!!
   *		  (TODO:: how to call it in y_jac for non-linear term????!!!)
   * =====================================================================================
   **/
  void Bilinear(const Teuchos::RCP<Epetra_Vector>& ud,
                const Teuchos::RCP<Epetra_MultiVector>& u,
                const Teuchos::RCP<Epetra_MultiVector>& v,
                const Teuchos::RCP<Epetra_MultiVector>& uv);
  void HBilinV();
  // void BilinTerm();
  void jacBilinTerm();
  void NonLinRHS();
  void LinCoeff();
  void computeExpVVyVy();
  void y_rhs();
  void y_jac();
  void Solve();
  void StochasticIterations();
  void computeExpVal();
  int BlockDavidsonMethod();
  void PostProcess(Epetra_Vector& second_mmnt);
  void Newton();
  void RunBackTracking();
  void CreateLocMultiVec(const std::string& WhichOne);
  void CreateDistTransMultivec();
  bool Converged() { return isConverged_; }
  int Iterations() { return iter_; }
  Teuchos::RCP<Epetra_SerialDenseMatrix> getJacobian() { return Yjac; }
  Teuchos::RCP<Epetra_MultiVector> getEyy()         { return Exp_yy_; }
  Teuchos::RCP<Epetra_MultiVector> getExpDExpyy()   { return ExpDExpyy; }
  Teuchos::RCP<Epetra_MultiVector> getY()           { return y_; }
  Epetra_MultiVector&              getYTrans()       { return *yTrans_; }

  int numSubTimeStep, MyLDA;
  int m_, ttt, rszyy;
  int stochiter;
  int N_, MyPID, itrtr = 0;
  double subdt_;
  double* dt_;
  Teuchos::RCP<Epetra_CrsMatrix> A_, ExpYY;
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
  Teuchos::RCP<Epetra_MultiVector> ExpVar, Vexpv4, rep_exp_vyvy, row_one;
  Teuchos::RCP<Epetra_MultiVector> Vtemp;
  Teuchos::RCP<Epetra_Vector>
    expv4; // map of udet from deterministic part( n*1 length vector )
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
    ExpVyVyyDExpyy, ExpDExpyy, eye;
  Epetra_SerialDenseSVD y_prob;
  Teuchos::RCP<Epetra_MultiVector> H, Hn, VHn, Rv, YY;
  Teuchos::RCP<Epetra_Map> map_expv4, Tmap, map_mm;
  Teuchos::RCP<Epetra_MultiVector> Exp_zy, Exp_zy_, ExpzyDExpyy, Exp_VyVyy;
  Teuchos::RCP<Epetra_MultiVector> Exp_yy, Exp_yy_, LocExpyyy, GlobExpyyy,
    EyyTyT, Utmp;
  Teuchos::RCP<Epetra_SerialDenseMatrix> sol;
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
  void computeExpDExpyy();
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
  int numPrecRecomputes_;
  double sig_;
  int backTrack_;
  int numBackTrackingSteps_;
  double toleranceRHS_;
  double NormRHS_;
  double NormRHStest_;
  int stochiter_;
  int NumStochIter_;
  //!
  bool isInitialized_, test_, debug_, useNwtn_, shouldEquil;
  bool isConverged_;
  bool backTracking_; // perhaps call this enableBacktracking_
  bool showGetInvalidArg_;
  Teuchos::RCP<Epetra_Comm> Comm_;
  Teuchos::RCP<Epetra_LocalMap> map_z_, map_x_;
  Teuchos::RCP<Epetra_Vector> p_;
  Teuchos::RCP<Epetra_CrsGraph> W_graph_;
  Teuchos::RCP<Epetra_CrsMatrix> stress_;
  Teuchos::RCP<Epetra_MultiVector> exp_yy_;
#if need_locaInterface == 1
  Teuchos::RCP<FVM::Domain> domain_;
#else  /* -----  not NEED_LOCAINTERFACE  ----- */
  Teuchos::RCP<Mean> model_;
#endif /* -----  not NEED_LOCAINTERFACE  ----- */
};
#endif

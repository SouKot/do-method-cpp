/*
 * =====================================================================================
 *
 *
 *       Filename:  StochCoeff.cpp
 *
 *    Description:  Class to solve the stochastic system
 *
 *        Version:  1.0
 *        Created:  01/15/2016 04:09:52 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sourabh Kotnala & Eric Mulder(Newton method and backtracking
 *                  functions in this class are largely a slight modification of
 *                  his template class)
 *
 *   Organization:
 *
 * =====================================================================================
 */
#include "CoeffSolver.hpp"
#include "StochIO.hpp"
#include "Epetra_LinearProblem.h"
#include "Epetra_Operator.h"
#include "Epetra_Time.h"
#include "HYMLS_MatrixUtils.hpp"
#include "Teuchos_LAPACK.hpp"
//#include "matplotlibcpp.h"
#include <Teuchos_RCPDecl.hpp>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>
//#include <Epetra_SerialDenseSVD.h>
#if need_locaInterface == 1
#include "FVM_model_interface.h"
#endif /* -----  not NEED_LOCAINTERFACE  ----- */

CoeffSolver::CoeffSolver(int numSamplesIn,
                 int numSubSteps,
                 int m,
                 double* dt,
                 const Teuchos::RCP<Epetra_CrsMatrix>& A,
                 const Teuchos::RCP<Epetra_Vector>& uMeanIn,
                 const DomainPtr& domain,
                 const Teuchos::RCP<Epetra_MultiVector>& Vn,
                 const Teuchos::RCP<Epetra_MultiVector>& W_init,
                 const Teuchos::RCP<Epetra_Comm>& comm,
                 const Teuchos::RCP<Teuchos::ParameterList>& CoefParams,
                 const Teuchos::RCP<StochasticState>& sharedState,
                 int maxNumIter,
                 bool useBacktracking,
                 double numBackTrackingSteps,
                 double toleranceRHS,
                 double initialRHSNorm)
  : Comm_(comm)
  , NumStochIter_(numSamplesIn)
  , numSubTimeStep(numSubSteps)
  , m_(m)
  , dt_(dt)
  , A_(A)
  , uMean(uMeanIn)
  , domain_(domain)
  , sharedState_(sharedState)
  , V(Vn)
  , W(W_init)
  , nDOF(uMeanIn->GlobalLength())
  , y_prob()
  , isConverged_(false)
  , backTracking_(useBacktracking)
  , iter_(0)
  , maxNumIterations_(maxNumIter)
  , toleranceRHS_(toleranceRHS)
  , NormRHS_(initialRHSNorm)
  , numBackTrackingSteps_(numBackTrackingSteps)
  , noiseGen_(A->Comm().NumProc(), A->Comm().MyPID())
{
  MyPID = A->Comm().MyPID();
  using Teuchos::rcp;
  localMapY = rcp(new Epetra_LocalMap(m_, 0, *Comm_));
  localMapZ = rcp(new Epetra_LocalMap(W->NumVectors(), 0, A_->Comm()));
  localMapYY = rcp(new Epetra_LocalMap(m_ * m_, 0, *Comm_));
  stochMap = rcp(new Epetra_Map(NumStochIter_, 0, *Comm_));

  yTrans_ = Teuchos::rcp(new Epetra_MultiVector(*stochMap, m_));
  zTrans_ = Teuchos::rcp(new Epetra_MultiVector(*stochMap, W->NumVectors()));
  numSamples = yTrans_->MyLength();
  y_ = Teuchos::rcp(new Epetra_MultiVector(*localMapY, numSamples));
  z_ = Teuchos::rcp(new Epetra_MultiVector(*localMapZ, numSamples));
  yyOuter = Teuchos::rcp(new Epetra_MultiVector(*localMapYY, numSamples));
  sizeYY = m_ * m_;
  yOld = rcp(new Epetra_MultiVector(*localMapY, numSamples));
  yCurr = rcp(new Epetra_MultiVector(*localMapY, numSamples));
  yDelta = rcp(new Epetra_MultiVector(*localMapY, numSamples));
  residual = Teuchos::rcp(new Epetra_Vector(*localMapY, numSamples));
  int n;
  n = A->NumMyRows();
  bilinWork = rcp(new Epetra_MultiVector(*V));
  rhs = rcp(new Epetra_MultiVector(*localMapY, numSamples));
#if need_locaInterface == 1
  // Allocate Fortran-map work vectors for the LOCA/SWE bilinear interface.
  Teuchos::RCP<Epetra_Map> fortMap = domain_->GetAssemblyMap();
  fortudet = Teuchos::rcp(new Epetra_Vector(*fortMap));
  fortVn = Teuchos::rcp(new Epetra_MultiVector(*fortMap, m_));
  fortviv = Teuchos::rcp(new Epetra_MultiVector(*fortMap, m_));
#endif
  AV = rcp(new Epetra_MultiVector(*V));
  VAV = rcp(new Epetra_MultiVector(*localMapY, m_));
  udet_mtimes = rcp(new Epetra_MultiVector(*V));
  ones = rcp(new Epetra_Vector(*localMapY));
  Vudet = rcp(new Epetra_MultiVector(*V));
  udetV = rcp(new Epetra_MultiVector(*V));
  VVudet = rcp(new Epetra_MultiVector(*localMapY, m_));
  lin_coeff = rcp(new Epetra_MultiVector(*localMapY, m_));
  JacNonLin = rcp(new Epetra_MultiVector(*localMapY, m_));
  rhsNonLin = rcp(new Epetra_MultiVector(*localMapY, numSamples));
  Ezy = rcp(new Epetra_MultiVector(*localMapZ, m_));
  EzyPrev = rcp(new Epetra_MultiVector(*localMapZ, m_));
  Eyy = rcp(new Epetra_MultiVector(*localMapY, m_));
  EyyTyT = rcp(new Epetra_MultiVector(*localMapYY, m_));
  EzyDEyy = rcp(new Epetra_MultiVector(*localMapZ, m_));
  EVyVyyDEyy = rcp(new Epetra_MultiVector(*V));
  EDEyy = rcp(new Epetra_MultiVector(*V));
  LocEyyy = rcp(new Epetra_MultiVector(*localMapY, V->NumVectors()));

  EVyVyy = rcp(new Epetra_MultiVector(*V));

  workMxM = Teuchos::rcp(new Epetra_MultiVector(*localMapY, m_));
  Epetra_Map Glob_x_map(m_, 0, *Comm_);
  EYY = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, Glob_x_map, m_));
  int GlobInd;
  double val;
  for (int i = 0; i < m_; i++) {
    for (int j = i; j < m_; j++) {
      val = 1.0;
      GlobInd = EYY->InsertGlobalValues(i, 1, &val, &j);
    }
  }

  EYY->FillComplete();
  identity = Teuchos::rcp(new Epetra_MultiVector(Epetra_Map(m_, 0, A_->Comm()), m_));
  Vwork = rcp(new Epetra_MultiVector(*V));
  const_coeff = rcp(new Epetra_MultiVector(
    *localMapY, 1)); // map of uMean from deterministic part( n*1 length vector )
  repEVyVy = rcp(new Epetra_MultiVector(*localMapY, numSamples));
  for (int i = 0; i < numSamples; i++)
    (*repEVyVy)[i] = (*const_coeff)[0];

  EVyVy = rcp(new Epetra_Vector(
    *uMean)); // map of uMean from deterministic part( n*1 length vector )
  EVyVy->PutScalar(0.0);
  xndiff = Teuchos::rcp(new Epetra_Vector(*localMapY, false));
  VB = rcp(new Epetra_MultiVector(*localMapY, W->NumVectors()));
  dW = rcp(new Epetra_MultiVector(*z_));
  VBdW = rcp(new Epetra_MultiVector(*yCurr));
  jacView = Teuchos::rcp(new Epetra_MultiVector(*lin_coeff));
  {
    double* jac_val;
    jacView->ExtractView(&jac_val, &m_);
    jacDense = rcp(new Epetra_SerialDenseMatrix(Epetra_DataAccess::View, jac_val, m_, m_, m_));
  }
  BilinTensor = rcp(new Epetra_MultiVector(*localMapY, m_ * m_));
  BilinTensorN = rcp(new Epetra_MultiVector(V->Map(), m_ * m_));
  VtBilinTensorN = rcp(new Epetra_MultiVector(*localMapY, m_ * m_));
  EyyInvView = Teuchos::rcp(new Epetra_MultiVector(*Eyy));
  {
    double* r_val;
    EyyInvView->ExtractView(&r_val, &m_);
    EyyInvTeuchos = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double>(
      Teuchos::View, r_val, m_, m_, m_));
    EyyInvEpetra = Epetra_SerialDenseMatrix(Epetra_DataAccess::View, r_val, m_, m_, m_);
  }

  // WARNING: At present the 'TypeCoeffFile' can be either 'None' or
  // 'Variance'. Using 'CoeffMatrix' will result in an error !!!!
  std::string CoefFile;
  std::string TypeCoeffFile = CoefParams->get("Type of Coeff File", "None");
  test_ = CoefParams->get("Class Testing", false);
  debug_ = CoefParams->get("Class Debugging", false);
  useNwtn_ = CoefParams->get("Use Newton", true);
  if (TypeCoeffFile == "None") {
    if (MyPID == 0)
      std::cout << "Initializing Stoch. Coeff to zero \n";
    y_->PutScalar(0.0);
  }
  if (TypeCoeffFile != "None")
    CoefFile = CoefParams->get("StochCoefFile", "y.mm");
  if (TypeCoeffFile == "CoeffMatrix")
  {
    if (MyPID == 0)
      std::cout << " Initializing Coeff.'s using " << CoefFile << ".\n";
    auto y_coeff = StochIO::readMV(CoefFile, *stochMap);
    *yTrans_ = *y_coeff;
    CreateLocalMultiVec("y");
  }
  if (TypeCoeffFile == "Variance") {
    if (MyPID == 0)
      std::cout << " Using provided variance vector for initializing Stoch. "
                << "Coeff.'s \n";
    auto y_coeff = StochIO::readMV(CoefFile, *localMapY);
    for (int i = 0; i < yTrans_->NumVectors(); i++)
      for (int j = 0; j < yTrans_->MyLength(); j++)
        (*(*yTrans_)(i))[j] = sqrt((*(*y_coeff)(0))[i]) * noiseGen_.sample();
    CreateLocalMultiVec("y");
  }
  double temval = 1.0, one = 1.0;
  if (debug_) {
    std::cout << "\nnum vectors in locexpyy = " << LocEyyy->NumVectors()
              << "\n";
  }
  for (int i = 0; i < identity->MyLength(); i++) {
    int gr = identity->Map().GID(i);
    identity->ReplaceGlobalValue(gr, gr, one);
  }
  if (debug_)
    printnormMV(*V, 2, "nrm of each vector in Vn:");

  // Publish owned data into the shared state so BasisSolver sees it.
  sharedState_->y         = y_;
  sharedState_->EDEyy = EDEyy;

} /* end of constructor */

#if 0
  // Test on consistency rhs and bilinear form after pitchfork bifurcation
  // If xp, xn are the solutions of the two stable branches then the average
  // x=(xp+xn)/2 and
  // half the difference v=(xp-xn)/2 satisfy f(x)+<v,v>=0.
  // We create xp and xn by solving for the steady states using appropriate
  // initial conditions
  // Next we compute in matlab the average and difference. Here we read them and
  // perform the test.
  Teuchos::RCP<Epetra_Map> xmap = rcp(new Epetra_Map(n, 0, *Comm_));
  Teuchos::RCP<Epetra_Vector> vt, diff, x, v;
  // Epetra_MultiVector *xptr, *vptr;
  vt = rcp(new Epetra_Vector(*xmap));
  diff = rcp(new Epetra_Vector(*xmap));
  x = rcp(new Epetra_Vector(*xmap));
  v = rcp(new Epetra_Vector(*xmap));
  double nrm;
  Epetra_CrsMatrix testJac(*A);
  testJac.FillComplete();
  //    testJac.ExtractCrsDataPointers(rows,cols,values);
  x->PutScalar(0.0);
  v->Random();
  model->computeJacobian(*x, testJac);
  Epetra_Vector J0v(*x);
  Epetra_Vector Jxv(*x);
  testJac.Multiply(false, *v, J0v); // compute J(0)*v
  x->Random();
  model->computeJacobian(*x, testJac);
  testJac.Multiply(false, *v, Jxv); // compute J(x)*v
  // compute rhs = J(0)*v - Bilin(x,v) - Bilin(v,x) - J(x)*v
  Bilinear(v, x, vt);
  diff->Update(-1.0, *vt, 0.0);
  Bilinear(x, v, vt);
  diff->Update(-1.0, *vt, 1.0);
  diff->Update(1.0, J0v, -1.0, Jxv, 1.0);
  diff->Norm2(&nrm);
  std::cout << "\n===========================================================";
  std::cout << "\n TEST FOR BILINEAR AND JACOBIAN";
  std::cout << "\n Norm of J(0)*v - Bilin(x,v) - Bilin(v,x) - J(x)*v = " << nrm;
  std::cout << "\n===========================================================";
  namespace plt = matplotlibcpp;
  double arr[diff->GlobalLength()];
  diff->ExtractCopy(arr);
  std::vector<double> vec(arr, arr + sizeof(arr) / sizeof(double));
  plt::plot(vec);
  plt::show();
  getchar();
#endif

/* ******************************************************** */

void
CoeffSolver::CreateLocalMultiVec(const std::string& WhichOne)
{
  if (WhichOne == "y") {
    for (int i = 0; i < yTrans_->MyLength(); i++) {
      for (int j = 0; j < yTrans_->NumVectors(); j++) {
        (*y_)[i][j] = (*yTrans_)[j][i];
      }
    }
  }
  if (WhichOne == "z") {
    for (int i = 0; i < zTrans_->MyLength(); i++) {
      for (int j = 0; j < zTrans_->NumVectors(); j++) {
        (*z_)[i][j] = (*zTrans_)[j][i];
      }
    }
  }
}

/* ******************************************************** */

void
CoeffSolver::localToDistributed()
{
  for (int i = 0; i < yTrans_->MyLength(); i++) {
    for (int j = 0; j < yTrans_->NumVectors(); j++) {
      (*yTrans_)[j][i] = (*y_)[i][j];
    }
  }
  for (int i = 0; i < zTrans_->MyLength(); i++) {
    for (int j = 0; j < zTrans_->NumVectors(); j++) {
      (*zTrans_)[j][i] = (*z_)[i][j];
    }
  }
}

/* ******************************************************** */

void
CoeffSolver::Bilinear(const Teuchos::RCP<Epetra_Vector>& ud,
                  const Teuchos::RCP<Epetra_MultiVector>& u,
                  const Teuchos::RCP<Epetra_MultiVector>& v,
                  const Teuchos::RCP<Epetra_MultiVector>& uv)
{
  int mu = u->NumVectors();
  int mv = v->NumVectors();
  int len;
  //  u->Random();
  //  v->Random();
  len = v->GlobalLength();
  double *udet_ptr, *u_ptr, *v_ptr, *uv_ptr;
  if (mu >= mv) {
    for (int i = 0; i < mu; i++) {
      for (int j = 0; j < mv; j++) {
#if need_locaInterface == 1
	ud->ExtractView(&udet_ptr);
        (*u)(i)->ExtractView(&u_ptr);
        (*v)(j)->ExtractView(&v_ptr);
        (*uv)(i)->ExtractView(&uv_ptr);
        model_bil(udet_ptr, u_ptr, v_ptr, uv_ptr);
#else
        /* -----  not NEED_LOCAINTERFACE  ----- */
        sharedState_->bilinearTerm(Teuchos::rcpFromRef(*((*u)(i))),
                             Teuchos::rcpFromRef(*((*v)(j))),
                             Teuchos::rcpFromRef(*((*uv)(i))));
#endif /* -----  not NEED_LOCAINTERFACE  ----- */
      }
    }
  } else {
    for (int i = 0; i < mu; i++) {
      for (int j = 0; j < mv; j++) {
#if need_locaInterface == 1
	ud->ExtractView(&udet_ptr);
        (*u)(i)->ExtractView(&u_ptr);
        (*v)(j)->ExtractView(&v_ptr);
        (*uv)(j)->ExtractView(&uv_ptr);
        model_bil(udet_ptr, u_ptr, v_ptr, uv_ptr);
#else  /* -----  not NEED_LOCAINTERFACE  ----- */
        sharedState_->bilinearTerm(Teuchos::rcpFromRef(*((*u)(i))),
                             Teuchos::rcpFromRef(*((*v)(j))),
                             Teuchos::rcpFromRef(*((*uv)(j))));
#endif /* -----  not NEED_LOCAINTERFACE  ----- */
      }
    }
  }
}

/* ******************************************************** */

/* BilinTensor - V'E<Vy,Vy>
 * BilinTensorN - E<Vy,Vy>
 */
void
CoeffSolver::computeBilinTensor()
{
  //  Epetra_MultiVector temp(*Eyy);
  int k;
  for (int i = 0; i < m_; i++) {

#if need_locaInterface == 1
    CHECK_ZERO(domain_->Solve2Assembly(*V, *fortVn));
    CHECK_ZERO(domain_->Solve2Assembly(*uMean, *fortudet));
    Bilinear(fortudet, Teuchos::rcp((*fortVn)(i), false), fortVn, fortviv);
    CHECK_ZERO(domain_->Assembly2Solve(*fortviv, *bilinWork));

#else  /* -----  not NEED_LOCAINTERFACE  ----- */
    Bilinear(uMean,Teuchos::rcp((*V)(i), false), V, bilinWork);
#endif /* -----  not NEED_LOCAINTERFACE  ----- */
    //   temp.Multiply('T', 'N', 1.0, *V, *bilinWork, 0.0);
    k = 0;
    for (int j = i * m_; j < (i + 1) * m_; j++) {
      //   (*BilinTensor)(j)->Update(1.0, *(temp(k)), 0.0);
      (*BilinTensorN)(j)->Update(1.0, *((*bilinWork)(k)), 0.0);
      k = k + 1;
    }
  }
  VtBilinTensorN->Multiply('T', 'N', 1.0, *V, *BilinTensorN, 0.0);
  if (debug_)
  {
    printnormMV(*V, 2, "nrm of each vector in V:");
    printnormMV(*BilinTensorN, 2, "nrm of each vector in BilinTensorN:");
  }
}

/* ******************************************************** */

void
CoeffSolver::computeRepEVyVy()
{
  /********* EVyVy = E[<Vy,Vy>] ********/
  // EVyVy->Scale((-subdt_));
  /********* repEVyVy = repmat(V' * (EVyVy), stochiter)*****/
  int err = const_coeff->Multiply('T', 'N', 1.0, *V, *EVyVy, 0.0);

  if (debug_) {
    std::cout << "\nError in calculating V*expvyvy : " << err << std::endl;
  }
  //  const_coeff->Scale(-1.0);
}

/* ******************************************************** */

void
CoeffSolver::computeNonlinearRHS()
{
  if (!useNwtn_) {
    rhsNonLin->Multiply('N', 'N', 1.0, *VtBilinTensorN, *yyOuter, 0.0);
  } else {
    rhsNonLin->Multiply('N', 'N', 1.0, *VtBilinTensorN, *yyOuter, 0.0);
  }
  if (debug_) {
    (*rhsNonLin)(0)->Print(std::cout << "\nrhsnonlin(0) values are \n");
  }
}

/* ******************************************************** */

void
CoeffSolver::computeBilinearJacobian()
{
  JacNonLin->Scale(0.0);
  for (int i = 0; i < m_; i++) {
    for (int j = 0; j < m_; j++) {
      (*JacNonLin)(i)->Update((*(*yOld)(i))[j], *((*BilinTensor)(i * m_ + j)), 1.0);
    }
  }
  for (int j = 0; j < m_; j++) {
    for (int i = 0; i < m_; i++) {
      (*JacNonLin)(j)->Update((*(*yOld)(i))[i], *((*BilinTensor)(i * m_ + j)), 1.0);
    }
  }
}

/* ******************************************************** */

void
CoeffSolver::computeLinearCoeff()
{
  A_->Multiply(false, *V, *AV);
  VAV->Multiply('T', 'N', 1.0, *V, *AV, 0.0);
  // lin_coeff = I - dt * (VAV)
  lin_coeff->Update(
    1.0, *VAV, 0.0); /* bilinear terms are not added in lin_coeff!!!!!! */
  lin_coeff->Scale(-1.0 * (subdt_));
  for (int i = 0; i < m_; i++) {
    lin_coeff->SumIntoMyValue(i, i, 1.00);
  }
  if (debug_) {
    std::cout << "\n frobenius norm of A : " << A_->NormFrobenius() << "\n";
    printnormMV(*lin_coeff, 2, "norm of y-jacView cols:");
    printnormMV(*VAV, 2, "norm of VAV cols:");
    printnormMV(*V, 2, "norm of V cols:");
  }
  if (debug_) {
    std::cout << "\n In CoeffSolver::computeLinearCoeff" << std::endl;
    double nrm1[V->NumVectors()];
    lin_coeff->Norm1(&nrm1[0]);
    std::cout << "One norm of I-dt*VAV" << std::endl;
    for (int i = 0; i < V->NumVectors(); ++i) {
      std::cout << nrm1[i] << "  ";
    }
    std::cout << std::endl;
    lin_coeff->MaxValue(&nrm1[0]);
    std::cout << "max value of I-dt*VAV" << std::endl;
    for (int i = 0; i < V->NumVectors(); ++i) {
      std::cout << nrm1[i] << "  ";
    }
    std::cout << std::endl;
    std::cout << "Norm of A = " << A_->NormFrobenius() << std::endl;
  }
  if (debug_) {

    getchar();
  }
}

/* ******************************************************** */

/*
 * ===  FUNCTION
 * ======================================================================
 *         Name:  assembleJacobian
 *  Description:  computes the jacobian. It uses SACADO package computing for
 * nonlinear
 *                part of jacobian.
 * =====================================================================================
 */
void
CoeffSolver::assembleJacobian()
{
  // computeBilinearJacobian();
  *jacView = *lin_coeff;
  // jacView->Update(subdt_,*JacNonLin,1.0);
}
/*
 * ===  FUNCTION
 * ======================================================================
 *         Name:  assembleRHS
 *  Description:  calculates rhs of the system!!!
 *  yCurr - the value obtained at previous time step
 *  yOld- the value to be obtained for current time step
 * =====================================================================================
 */
void
CoeffSolver::assembleRHS()
{
  if (!useNwtn_) {
    computeNonlinearRHS();
    rhs->PutScalar(0.0);
    rhs->Update(1.0, *yCurr, 1.0);
    // rhs=-x+lin_coeff*x0
    std::cout.precision(16);
    // rhs= x + dt*(V'<Vx,Vx>-V'E[<Vx,Vx>])
    rhs->Update(subdt_, *repEVyVy, 1.0);
    rhs->Update(-subdt_, *rhsNonLin, 1.0);
    rhs->Update(1.0, *VBdW, 1.0);
    if (debug_ && stochiter_ == 0) {
      double nrm;
      double nrm1[V->NumVectors()], nrm2[y_->NumVectors()];
      printnormMV(*V, 2, "norm of V: ");
      printnormMV(*rhsNonLin, 2, "norm of rhsnonlin: ");
      printnormMV(*zTrans_, 2, "norm of dW: ");
      printnormMV(*EyyOld_, 2, "norm of EyyOld_: ");
    }

    if (debug_ && stochiter_ == 0) {
    }
  } else {
    computeNonlinearRHS();
    rhs->PutScalar(0.0);
    rhs->Multiply('N', 'N', 1.0, *lin_coeff, *yOld, 1.0); // rhs=-x+lin_coeff*x0
    rhs->Update(-1.0, *yCurr, 1.0);
    rhs->Update(-subdt_, *repEVyVy, 1.0); // rhs=-x+lin_coeff*x0+const_coeff
    // rhs=-x + lin_coeff*x0 + const_coeff + rhs_nonlin
    rhs->Update(subdt_, *rhsNonLin, 1.0);
    rhs->Update(-1.0, *VBdW, 1.0);
  }
}

Teuchos::RCP<Epetra_Map>
CoeffSolver::get_x_map()
{
  return localMapY;
}

void
CoeffSolver::set_x(const Teuchos::RCP<Epetra_MultiVector>& x0_temp)
{
  yOld = x0_temp;
  *yCurr = *x0_temp;
}

void
CoeffSolver::setDwiener(const Teuchos::RCP<Epetra_MultiVector>& z0_temp)
{
  dW = z0_temp;

  if ((0)) { /* TODO: Put it in debug_ or testing_*/
    debugWienerCount += 1;
    std::string flnm = "wnr" + Teuchos::toString(debugWienerCount) + ".mm";
    std::cout << "\n filename is: " << flnm << "\n";
    auto tmpmv = StochIO::readMV(flnm, *localMapZ);
    *dW = *tmpmv;
  } else {
    dW->Scale(sqrt(subdt_));
  }
}

Teuchos::RCP<Epetra_Map>
CoeffSolver::get_f_map()
{
  return localMapY;
}

Teuchos::RCP<Epetra_MultiVector>
CoeffSolver::get_x_init()
{
  return yOld;
}

Teuchos::RCP<Epetra_Vector>
CoeffSolver::getEVyVy()
{
  return EVyVy;
}
void
CoeffSolver::Solve()
{
  double* Viewx;
  double* ViewRHS;
  yDelta->ExtractView(&Viewx, &m_);
  rhs->ExtractView(&ViewRHS, &m_);
  DenseDx_ =
    Epetra_SerialDenseMatrix(Epetra_DataAccess::Copy, Viewx, m_, m_, numSamples);
  DenseRHS =
    Epetra_SerialDenseMatrix(Epetra_DataAccess::Copy, ViewRHS, m_, m_, numSamples);

  y_prob.SetVectors(DenseDx_, DenseRHS);
  y_prob.Solve();
  if (debug_) {
    y_prob.Matrix()->Print(std::cout << "inverted matrix \n");
    jacDense->Print(std::cout << "jacDense \n");
    std::cout << "\n norm one jacDense: " << jacDense->NormOne() << std::endl;
    Epetra_SerialDenseMatrix res;
    res = DenseRHS;
    res.Multiply('N', 'N', 1.0, *jacDense, DenseDx_, -1.0);
    //res.Print(std::cout << "res \n");
    std::cout << "\n If Solved : " << y_prob.Solved() << std::endl;
    std::cout << "\n Debug_ value : " << debug_<< std::endl;
    getchar();
  }

  for (int i = 0; i < numSamples; i++)
    (*yDelta)[i] = DenseDx_[i];
}

void
CoeffSolver::StochasticIterations()
{
  subdt_ = *dt_ / numSubTimeStep;
  double time, LocTime;
  Epetra_Time timer1(V->Comm());

  computeEVyVy();

  if (test_) {
    LocTime = timer1.ElapsedTime();
    Comm_->SumAll(&LocTime, &time, 1);
    if (Comm_->MyPID() == 0) {
      std::cout << "\n*******************************************\n";
      std::cout << "average time per proc. for computeBilinTensor: "
                << (time) / (A_->Comm().NumProc()) << " sec\n";
      std::cout << "*******************************************\n\n";
    }
  }
  Epetra_Time timer2(timer1);

  computeLinearCoeff();

  if (test_) {
    LocTime = timer2.ElapsedTime();
    Comm_->SumAll(&LocTime, &time, 1);
    if (Comm_->MyPID() == 0) {
      std::cout << "\n*******************************************\n";
      std::cout << "average time per proc. for computeLinearCoeff: "
                << (time) / (A_->Comm().NumProc()) << " sec\n";
      std::cout << "*******************************************\n\n";
    }
  }
  Epetra_Time timer3(timer1);
  assembleJacobian();

  if (debug_) {
    std::cout.precision(18);
    std::cout << "\n frob-norm of det-jacView: " << A_->NormFrobenius();
    printnormMV(*jacView, 2, "2-norm of y-jacView");
  }
  Epetra_MultiVector JAC(*jacView);
  JAC = *jacView;

  if (test_) {
    LocTime = timer3.ElapsedTime();
    Comm_->SumAll(&LocTime, &time, 1);
    if (Comm_->MyPID() == 0) {
      std::cout << "\n*******************************************\n";
      std::cout << "average time per proc. for y-jacView() "
                << (time) / (A_->Comm().NumProc()) << " sec\n";
      std::cout << "*******************************************\n\n";
    }
  }
  Epetra_Time timer4(timer1);
  y_prob.SetMatrix(*jacDense);
  if (debug_) {
    printnormMV(*V, 1, "1-norm of V in CoeffSolver::StochasticIterations():");
    std::cout << "One Norm of jacDense = " << jacDense->NormOne() << std::endl;
  }
  // shouldEquil=y_prob.ShouldEquilibrate();
  //  if(shouldEquil)
  //  {
  //    y_prob.FactorWithEquilibration(true);
  //  }
  //    y_prob.SolveToRefinedSolution(true);
  //    y_prob.EstimateSolutionErrors(true);
  //
  y_prob.Invert(10e-10);
  double* valptr;
  valptr = y_prob.AI();
  Epetra_MultiVector inv_jac(
    Epetra_DataAccess::Copy, jacView->Map(), valptr, m_, m_);
  if (test_) {
    LocTime = timer4.ElapsedTime();
    Comm_->SumAll(&LocTime, &time, 1);
    if (Comm_->MyPID() == 0) {
      std::cout << "\n*******************************************\n";
      std::cout << "average time per proc. for factorization "
                << (time) / (A_->Comm().NumProc()) << " sec\n";
      std::cout << "*******************************************\n\n";
    }
  }

  VB->Multiply('T', 'N', 1.0, *V, *W, 0.0);
  Epetra_Time timer5(timer1);
  for (int substep = 1; substep <= numSubTimeStep; substep++) {
    double avrg[zTrans_->NumVectors()];
    Epetra_MultiVector avgMV(*zTrans_);
    for (int i = 0; i < W->NumVectors(); i++) {
      avrg[i] = 0;
      for (stochiter_ = 0; stochiter_ < yTrans_->MyLength(); stochiter_++) {
        (*(*zTrans_)(i))[stochiter_] = noiseGen_.sample();
        avrg[i] += (*(*zTrans_)(i))[stochiter_];
      }
      avrg[i] = avrg[i] / NumStochIter_;
      avgMV(i)->PutScalar(avrg[i]);
    }
    // Make expectation of 'z' exactly zero.
    zTrans_->Update(-1.0, avgMV, 1.0);

    std::string type = "z";
    CreateLocalMultiVec(type);
    set_x(y_);
    setDwiener(z_);
    /********* stoch_coeff = VB*dW ****************/
    VBdW->Multiply('N', 'N', 1.0, *VB, *dW, 0.0);
    if (!useNwtn_) {
      yDelta->PutScalar(0.0);
      assembleRHS();
      yDelta->Multiply('N', 'N', 1.0, inv_jac, *rhs, 0.0);
      yOld->Update(1.0, *yDelta, 0.0);
    } else
      Newton();
    if (debug_ && stochiter_ == 0) {
      double nrm1[V->NumVectors()];
      inv_jac.Norm2(&nrm1[0]);
      std::cout << "\n norm of inv(jacDense)\n";
      for (int i = 0; i < V->NumVectors(); ++i) {
        std::cout << nrm1[i] << "  ";
      }
      std::cout << std::endl;

      std::cout << "\n frob. norm of inv(assembleJacobian)  = " << jacDense->NormOne()
                << std::endl;
      getchar();
    }
  }
  if (test_) {
    LocTime = timer5.ElapsedTime();
    Comm_->SumAll(&LocTime, &time, 1);
    if (Comm_->MyPID() == 0) {
      std::cout << "\n*******************************************\n";
      std::cout << "average time per proc. for Iterations "
                << (time) / (A_->Comm().NumProc()) << " sec\n";
      std::cout << "*******************************************\n" << std::endl;
    }
  }
  computeExpectations();
}

void
CoeffSolver::computeExpectations()
{
  computeCrossVariance();
  computeEyyTyT();
  computeEVyVyy();
  computeEDEyy();
}

void
CoeffSolver::computeEVyVy()
{
  EVyVy->Scale(0.0);
  for (int i = 0; i < m_; i++) {
    for (int j = 0; j < m_; j++) {
      EVyVy->Update((*Eyy)[j][i], (*((*BilinTensorN)(i * m_ + j))), 1.0);
    }
  }
}

void
CoeffSolver::computeCrossVariance()
{
  localToDistributed(); /* Create zTrans and yTrans which have distributed
                                map */
  int tmp = Eyy->Multiply(
    'T', 'N', 1.0 / NumStochIter_, *yTrans_, *yTrans_, 0.0); /* Calculate yy' */
}

void
CoeffSolver::computeEyyTyT()
{
  double* ViewYY;
  double* ViewLocExpyyy;
  yyOuter->ExtractView(&ViewYY, &sizeYY);
  LocEyyy->ExtractView(&ViewLocExpyyy, &m_);
  Epetra_Vector tmpy(*(*y_)(0));
  int err1;
  for (int i = 0; i < numSamples; i++) {
    tmpy = (*(*y_)(i));
    err1 = LocEyyy->Multiply('N', 'T', 1.0, tmpy, tmpy, 0.0);
    for (int j = 0; j < yyOuter->MyLength(); j++) {
      ViewYY[i * yyOuter->MyLength() + j] = ViewLocExpyyy[j];
    }
  }
  double cpyEyyy[EyyTyT->NumVectors()*EyyTyT->MyLength()];
  const int sz=EyyTyT->NumVectors()*EyyTyT->MyLength();
  double globSumEyyy[EyyTyT->NumVectors()*EyyTyT->MyLength()];
  EyyTyT->Multiply('N', 'T', 1.0 / NumStochIter_, *yyOuter, *y_, 0.0);

  EyyTyT->ExtractCopy(cpyEyyy, EyyTyT->MyLength());
  int count = EyyTyT->NumVectors() * EyyTyT->MyLength();
  Comm_->SumAll(cpyEyyy, globSumEyyy, count);
  for (int i=0; i< EyyTyT->NumVectors(); i++)
  {
      for (int j=0; j< EyyTyT->MyLength(); j++)
      {
        (*(*EyyTyT)(i))[j]=globSumEyyy[i*EyyTyT->MyLength()+j];
      }
  }
}

void
CoeffSolver::computeEVyVyy()
{
  EVyVyy->Multiply('N', 'N', 1.0, *BilinTensorN, *EyyTyT, 0.0);
}
// ************************************************************************
void
CoeffSolver::computeEDEyy()
// ************************************************************************
{
  Epetra_MultiVector Inv_Exp_yy(Eyy->Map(), Eyy->NumVectors());
  SymMatPseudoInverse(*Eyy, Inv_Exp_yy); /* compute pseudo-inverse */
  Epetra_MultiVector zyT(*EzyDEyy), WzyT(*V);
  /* Calculate zy' */
  zyT.Multiply('T', 'N', 1.0 / NumStochIter_, *zTrans_, *yTrans_, 0.0);
  /* Calculate W*E[zy'] */
  WzyT.Multiply('N', 'N', numSubTimeStep, *W, zyT, 0.0);
  Epetra_MultiVector Wzy_P_EVyVyy(*V);
  /* Calculate W*E[zy']+E[<Vy,Vy>] */
  Wzy_P_EVyVyy.Update(*dt_, *EVyVyy, 1.0, WzyT, 0.0);
  /* Calculate (W*E[zy']+E[<Vy,Vy>])/E[yy'] */
  EDEyy->Multiply('N', 'N', 1.0, Wzy_P_EVyVyy, Inv_Exp_yy, 0.0);

  if (debug_) {
    Eyy->Print(std::cout << "Eyy \n");
    Inv_Exp_yy.Print(std::cout << "Inv_Exp_yy \n");
    Epetra_MultiVector identity(*Eyy);
    identity.Multiply('N', 'N', 1.0, Inv_Exp_yy, *Eyy, 0.0);
    identity.Print(std::cout << "Eyy * Inv(Eyy) \n");
    printnormMV(*V, 2, "norm of V : ");
    printnormMV(*yTrans_, 2, "norm of Y^T : ");
    printnormMV(*W, 2, "norm of W : ");
    printnormMV(*VB, 2, "norm of VW : ");
    printTransposedNorms(*rhs, 2, "norm of rhs^T: ");
    printTransposedNorms(*dW, 2, "norm of dw^T: ");
    printnormMV(*EyyTyT, 2, "norm of Eyyy: ");
    printnormMV(Wzy_P_EVyVyy, 2, "norm of WE[zyT]+E[VyVy]: ");
    printnormMV(*EDEyy, 2, "norm of EDEyy: ");
  }
}

void
CoeffSolver::SymMatPseudoInverse(Epetra_MultiVector& mat,
                             Epetra_MultiVector& matInv)
{
  int myLDA = mat.GlobalLength();
  double matView[myLDA * myLDA];
  mat.ExtractCopy(matView, myLDA);
  Epetra_SerialDenseMatrix matSDM(
    Epetra_DataAccess::Copy, matView, myLDA, myLDA, myLDA);
  Epetra_SerialDenseSVD invsol;
  invsol.SetMatrix(matSDM);
  invsol.Factor();
  double *u, *vt, *s;
  u = invsol.U_;
  vt = invsol.Vt_;
  s = invsol.S_;
  Epetra_MultiVector U(Epetra_DataAccess::View, mat.Map(), u, myLDA, myLDA);
  Epetra_MultiVector Vt(Epetra_DataAccess::View, mat.Map(), vt, myLDA, myLDA);
  Epetra_Vector S(Epetra_DataAccess::View, mat.Map(), s);
  double tol = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < S.MyLength(); ++i) {
    if (abs(S[i]) <= tol)
      S[i] = 0.0;
    else
      S[i] = 1.0 / S[i];
    U(i)->Scale(S[i]);
  }
  matInv.Multiply('N', 'N', 1.0, U, Vt, 0.0);
}
// ************************************************************************
void
CoeffSolver::PostProcess(Epetra_Vector& second_mmnt)

// ************************************************************************
{
  using namespace std;
  localToDistributed();
  double A[m_][m_];
  int LWORK = 1 + 5 * m_ * m_, intdum = 1;
  double WORK[LWORK], singValues[m_], U[m_][m_], Vdum;
  char COMPZ = 'V', UPLO = 'U';
  int INFO;
  char JOBU = 'A', JOBVT = 'O';
  Epetra_LAPACK lapack;
  Epetra_MultiVector yy(*localMapY, m_);

  yy.Multiply('T', 'N', 1.0, *yTrans_, *yTrans_, 0.0);
  yy.ExtractCopy(&A[0][0], m_);
  lapack.GESVD(JOBU,
               JOBVT,
               m_,
               m_,
               &A[0][0],
               m_,
               &singValues[0],
               &U[0][0],
               m_,
               &Vdum,
               intdum,
               &WORK[0],
               &LWORK,
               &INFO);
  if (A_->Comm().MyPID() == 0) {
    std::ofstream DeviationFile("DevFile.mm", ios::out);
    std::cout << "\nDeviations:\n";
    DeviationFile << "%%MatrixMarket matrix array real general\n"
                  << m_ << " " << 1 << "\n";
    DeviationFile << scientific;
    DeviationFile.precision(16);
    for (int i = 0; i < m_; i++) {
      std::cout << "  " << sqrt(singValues[i] / float(NumStochIter_)) << "\n";
      DeviationFile << sqrt(singValues[i] / float(NumStochIter_)) << "\n";
    }
    std::cout << "\nVariances:\n";
    std::ofstream VarianceFile("VarFile.mm", ios::out);
    VarianceFile << "%%MatrixMarket matrix array real general\n"
                 << m_ << " " << 1 << "\n";
    VarianceFile << scientific;
    VarianceFile.precision(16);
    for (int i = 0; i < m_; i++) {
      std::cout << "  " << singValues[i] / NumStochIter_ << "\n";
      VarianceFile << singValues[i] / NumStochIter_ << "\n";
    }
  }
  int err;
  second_mmnt.PutScalar(0.0);
  for (int i = 0; i < V->NumVectors(); i++) {
    for (int j = 0; j < V->NumVectors(); j++) {
      second_mmnt.Multiply(
        yy[i][j] / NumStochIter_, *(*V)(i), *(*V)(j), 1.0);
    }
  }
}

void
CoeffSolver::Newton()
{
  isConverged_ = false;
  yDelta->PutScalar(0.0);
  // assembleJacobian();
  assembleRHS();
  rhs->Scale(-1.0);
  std::vector<double> normsRHS(rhs->NumVectors());
  rhs->Norm2(normsRHS.data());
  NormRHS_ = *std::max_element(normsRHS.begin(), normsRHS.end());
  //std::cout << "\n norm of rhs = " << NormRHS_ << std::endl;
  //assembleJacobian();
  //y_prob.SetMatrix(*jacDense); // 2. Give the updated matrix to the SVD solver
  //y_prob.Factor();
  INFO("Newton:      norm: " << NormRHStest_ );
  for (iter_ = 0; iter_ != maxNumIterations_; ++iter_) {
    Solve();
    if (debug_)
      yOld->Update(1.0, *yDelta, 0.0);
    else {
      yOld->Update(1.0, *yDelta, 1.0);
      assembleRHS();
      rhs->Scale(-1.0);
    }
    std::vector<double> normsTest(rhs->NumVectors());
    rhs->Norm2(normsTest.data());
    NormRHStest_ = *std::max_element(normsTest.begin(), normsTest.end());
    if (debug_) {
      std::cout << "\n RHS NORM = " << NormRHStest_ << std::endl;
      Epetra_Vector Err(*(*y_)(0));
      Err.Multiply('N', 'N', 1.0, *jacView, *yOld, 0.0);
      Err.Update(1.0, *rhs, 1.0);
      double nrm;
      Err.Norm2(&nrm);
      std::cout << "\n norm of err = jacView*x0-rhs = " << nrm << std::endl;
    } else {
      if (NormRHStest_ < toleranceRHS_) {
        DEBUG("Success...");
        break;
      }
      if (backTracking_ and (NormRHS_ < NormRHStest_))
        runBackTracking();
      NormRHS_ = NormRHStest_;
    }
  }
  if (iter_ == maxNumIterations_) {
    std::cout << "Newton: ---> TROUBLE" << __FILE__ << __LINE__;
    std::cout << "\nNewton not converged after Max number of iter!!!!!!\n";
    std::cout << " \n Stoch. iter. no. = " << stochiter_ << std::flush;
    std::cout << "\nrhs norm = " << NormRHStest_;
    getchar();
  } else {
    isConverged_ = true;
  }
}
//======================================================================
void
CoeffSolver::runBackTracking()
{
  // Initialize reduction with -1/2
  double reduction = -1.0 / 2;
  for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_) {
    if (NormRHStest_ < NormRHS_) {
      DEBUG("Success...");
      break;
    }
    // Apply reduction to the state vector
    yOld->Update(reduction, *yDelta, 1.0);
    assembleRHS();
    std::vector<double> normsBT(rhs->NumVectors());
    rhs->Norm2(normsBT.data());
    NormRHStest_ = *std::max_element(normsBT.begin(), normsBT.end());
    if (debug_) {

      std::cout << "Newton: --> backtracking:\n "
                << " step: " << backTrack_ << "\n reduction: " << reduction
                << "\n norm: " << NormRHStest_;
    }
    // Update reduction
    reduction /= 2.0;
  }
  if (backTrack_ == numBackTrackingSteps_)
    std::cout << "\nNewton: --> BACKTRACKING FAILED" << __FILE__ << __LINE__;
}

void
CoeffSolver::printTransposedNorms(Epetra_MultiVector& mv, int normType, const std::string& str)
{
  Epetra_Map mp(mv.NumVectors(), 0, mv.Comm());
  Epetra_MultiVector mvT(mp, mv.MyLength());
  for (int i = 0; i < mvT.MyLength(); i++) {
    for (int j = 0; j < mvT.NumVectors(); j++) {
      mvT[j][i] = mv[i][j];
    }
  }
  double nrm1[mvT.NumVectors()];
  if (normType == 1)
    mvT.Norm1(&nrm1[0]);
  else if (normType == 2)
    mvT.Norm2(&nrm1[0]);
  if (MyPID==0)
  {
    std::cout << "\n" << str << std::endl;
    for (int i = 0; i < mvT.NumVectors(); i++)
      std::cout << "  " << nrm1[i];
  }
}

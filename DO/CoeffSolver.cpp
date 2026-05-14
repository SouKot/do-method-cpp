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

CoeffSolver::CoeffSolver(int NumStochIter,
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
                 double NormRHS)
  : Comm_(comm)
  , NumStochIter_(NumStochIter)
  , numSubTimeStep(num_Subtime_Step)
  , m_(m)
  , dt_(dt)
  , A_(A)
  , udet(uav)
  , domain_(domain)
  , sharedState_(sharedState)
  , Vnew(Vn)
  , B(Wb)
  , N_(uav->GlobalLength())
  , y_prob()
  , isConverged_(false)
  , backTracking_(useBacktracking)
  , iter_(0)
  , maxNumIterations_(maxNumIter)
  , toleranceRHS_(toleranceRHS)
  , NormRHS_(NormRHS)
  , numBackTrackingSteps_(numBackTrackingSteps)
  , noiseGen_(A->Comm().NumProc(), A->Comm().MyPID())
{
  MyPID = A->Comm().MyPID();
  using Teuchos::rcp;
  map_x_ = rcp(new Epetra_LocalMap(m_, 0, *Comm_));
  map_z_ = rcp(new Epetra_LocalMap(B->NumVectors(), 0, A_->Comm()));
  map_mm = rcp(new Epetra_LocalMap(m_ * m_, 0, *Comm_));
  Tmap = rcp(new Epetra_Map(NumStochIter_, 0, *Comm_));

  yTrans_ = Teuchos::rcp(new Epetra_MultiVector(*Tmap, m_));
  zTrans_ = Teuchos::rcp(new Epetra_MultiVector(*Tmap, B->NumVectors()));
  MyLDA = yTrans_->MyLength();
  y_ = Teuchos::rcp(new Epetra_MultiVector(*map_x_, MyLDA));
  z_ = Teuchos::rcp(new Epetra_MultiVector(*map_z_, MyLDA));
  YY = Teuchos::rcp(new Epetra_MultiVector(*map_mm, MyLDA));
  rszyy = m_ * m_;
  x0_ = rcp(new Epetra_MultiVector(*map_x_, MyLDA));
  x_ = rcp(new Epetra_MultiVector(*map_x_, MyLDA));
  dx_ = rcp(new Epetra_MultiVector(*map_x_, MyLDA));
  f_out = Teuchos::rcp(new Epetra_Vector(*map_x_, MyLDA));
  int n;
  n = A->NumMyRows();
  viv = rcp(new Epetra_MultiVector(*Vnew));
  rhs = rcp(new Epetra_MultiVector(*map_x_, MyLDA));
#if need_locaInterface == 1
  // Allocate Fortran-map work vectors for the LOCA/SWE bilinear interface.
  Teuchos::RCP<Epetra_Map> fortMap = domain_->GetAssemblyMap();
  fortudet = Teuchos::rcp(new Epetra_Vector(*fortMap));
  fortVn = Teuchos::rcp(new Epetra_MultiVector(*fortMap, m_));
  fortviv = Teuchos::rcp(new Epetra_MultiVector(*fortMap, m_));
#endif
  AV = rcp(new Epetra_MultiVector(*Vnew));
  VAV = rcp(new Epetra_MultiVector(*map_x_, m_));
  udet_mtimes = rcp(new Epetra_MultiVector(*Vnew));
  ones = rcp(new Epetra_Vector(*map_x_));
  Vudet = rcp(new Epetra_MultiVector(*Vnew));
  udetV = rcp(new Epetra_MultiVector(*Vnew));
  VVudet = rcp(new Epetra_MultiVector(*map_x_, m_));
  lin_coeff = rcp(new Epetra_MultiVector(*map_x_, m_));
  JacNonLin = rcp(new Epetra_MultiVector(*map_x_, m_));
  rhsNonLin = rcp(new Epetra_MultiVector(*map_x_, MyLDA));
  Ezy = rcp(new Epetra_MultiVector(*map_z_, m_));
  EzyPrev = rcp(new Epetra_MultiVector(*map_z_, m_));
  Eyy = rcp(new Epetra_MultiVector(*map_x_, m_));
  EyyTyT = rcp(new Epetra_MultiVector(*map_mm, m_));
  EzyDEyy = rcp(new Epetra_MultiVector(*map_z_, m_));
  EVyVyyDEyy = rcp(new Epetra_MultiVector(*Vnew));
  EDEyy = rcp(new Epetra_MultiVector(*Vnew));
  LocEyyy = rcp(new Epetra_MultiVector(*map_x_, Vnew->NumVectors()));

  EVyVyy = rcp(new Epetra_MultiVector(*Vnew));

  Utmp = Teuchos::rcp(new Epetra_MultiVector(*map_x_, m_));
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
  eye = Teuchos::rcp(new Epetra_MultiVector(Epetra_Map(m_, 0, A_->Comm()), m_));
  Vtemp = rcp(new Epetra_MultiVector(*Vnew));
  const_coeff = rcp(new Epetra_MultiVector(
    *map_x_, 1)); // map of udet from deterministic part( n*1 length vector )
  repEVyVy = rcp(new Epetra_MultiVector(*map_x_, MyLDA));
  for (int i = 0; i < MyLDA; i++)
    (*repEVyVy)[i] = (*const_coeff)[0];

  EVyVy = rcp(new Epetra_Vector(
    *udet)); // map of udet from deterministic part( n*1 length vector )
  EVyVy->PutScalar(0.0);
  xndiff = Teuchos::rcp(new Epetra_Vector(*map_x_, false));
  VB = rcp(new Epetra_MultiVector(*map_x_, B->NumVectors()));
  dW = rcp(new Epetra_MultiVector(*z_));
  VBdW = rcp(new Epetra_MultiVector(*x_));
  jac = Teuchos::rcp(new Epetra_MultiVector(*lin_coeff));
  {
    double* jac_val;
    jac->ExtractView(&jac_val, &m_);
    Yjac = rcp(new Epetra_SerialDenseMatrix(Epetra_DataAccess::View, jac_val, m_, m_, m_));
  }
  H = rcp(new Epetra_MultiVector(*map_x_, m_ * m_));
  Hn = rcp(new Epetra_MultiVector(Vnew->Map(), m_ * m_));
  VHn = rcp(new Epetra_MultiVector(*map_x_, m_ * m_));
  Rv = Teuchos::rcp(new Epetra_MultiVector(*Eyy));
  {
    double* r_val;
    Rv->ExtractView(&r_val, &m_);
    Ru = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double>(
      Teuchos::View, r_val, m_, m_, m_));
    ru = Epetra_SerialDenseMatrix(Epetra_DataAccess::View, r_val, m_, m_, m_);
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
    auto y_coeff = StochIO::readMV(CoefFile, *Tmap);
    *yTrans_ = *y_coeff;
    CreateLocMultiVec("y");
  }
  if (TypeCoeffFile == "Variance") {
    if (MyPID == 0)
      std::cout << " Using provided variance vector for initializing Stoch. "
                << "Coeff.'s \n";
    auto y_coeff = StochIO::readMV(CoefFile, *map_x_);
    for (int i = 0; i < yTrans_->NumVectors(); i++)
      for (int j = 0; j < yTrans_->MyLength(); j++)
        (*(*yTrans_)(i))[j] = sqrt((*(*y_coeff)(0))[i]) * noiseGen_.sample();
    CreateLocMultiVec("y");
  }
  double temval = 1.0, one = 1.0;
  if (debug_) {
    std::cout << "\nnum vectors in locexpyy = " << LocEyyy->NumVectors()
              << "\n";
  }
  for (int i = 0; i < eye->MyLength(); i++) {
    int gr = eye->Map().GID(i);
    eye->ReplaceGlobalValue(gr, gr, one);
  }
  if (debug_)
    printnormMV(*Vnew, 2, "nrm of each vector in Vn:");

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
CoeffSolver::CreateLocMultiVec(const std::string& WhichOne)
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
CoeffSolver::CreateDistTransMultivec()
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

/* H - V'E<Vy,Vy>
 * Hn - E<Vy,Vy>
 */
void
CoeffSolver::HBilinV()
{
  //  Epetra_MultiVector temp(*Eyy);
  int k;
  for (int i = 0; i < m_; i++) {

#if need_locaInterface == 1
    CHECK_ZERO(domain_->Solve2Assembly(*Vnew, *fortVn));
    CHECK_ZERO(domain_->Solve2Assembly(*udet, *fortudet));
    Bilinear(fortudet, Teuchos::rcp((*fortVn)(i), false), fortVn, fortviv);
    CHECK_ZERO(domain_->Assembly2Solve(*fortviv, *viv));

#else  /* -----  not NEED_LOCAINTERFACE  ----- */
    Bilinear(udet,Teuchos::rcp((*Vnew)(i), false), Vnew, viv);
#endif /* -----  not NEED_LOCAINTERFACE  ----- */
    //   temp.Multiply('T', 'N', 1.0, *Vnew, *viv, 0.0);
    k = 0;
    for (int j = i * m_; j < (i + 1) * m_; j++) {
      //   (*H)(j)->Update(1.0, *(temp(k)), 0.0);
      (*Hn)(j)->Update(1.0, *((*viv)(k)), 0.0);
      k = k + 1;
    }
  }
  VHn->Multiply('T', 'N', 1.0, *Vnew, *Hn, 0.0);
  if (debug_)
  {
    printnormMV(*Vnew, 2, "nrm of each vector in V:");
    printnormMV(*Hn, 2, "nrm of each vector in Hn:");
  }
}

/* ******************************************************** */

void
CoeffSolver::computeRepEVyVy()
{
  /********* EVyVy = E[<Vy,Vy>] ********/
  // EVyVy->Scale((-subdt_));
  /********* repEVyVy = repmat(V' * (EVyVy), stochiter)*****/
  int err = const_coeff->Multiply('T', 'N', 1.0, *Vnew, *EVyVy, 0.0);

  if (debug_) {
    std::cout << "\nError in calculating V*expvyvy : " << err << std::endl;
  }
  //  const_coeff->Scale(-1.0);
}

/* ******************************************************** */

void
CoeffSolver::NonLinRHS()
{
  if (!useNwtn_) {
    rhsNonLin->Multiply('N', 'N', 1.0, *VHn, *YY, 0.0);
  } else {
    rhsNonLin->Multiply('N', 'N', 1.0, *VHn, *YY, 0.0);
  }
  if (debug_) {
    (*rhsNonLin)(0)->Print(std::cout << "\nrhsnonlin(0) values are \n");
  }
}

/* ******************************************************** */

void
CoeffSolver::jacBilinTerm()
{
  JacNonLin->Scale(0.0);
  for (int i = 0; i < m_; i++) {
    for (int j = 0; j < m_; j++) {
      (*JacNonLin)(i)->Update((*(*x0_)(i))[j], *((*H)(i * m_ + j)), 1.0);
    }
  }
  for (int j = 0; j < m_; j++) {
    for (int i = 0; i < m_; i++) {
      (*JacNonLin)(j)->Update((*(*x0_)(i))[i], *((*H)(i * m_ + j)), 1.0);
    }
  }
}

/* ******************************************************** */

void
CoeffSolver::LinCoeff()
{
  A_->Multiply(false, *Vnew, *AV);
  VAV->Multiply('T', 'N', 1.0, *Vnew, *AV, 0.0);
  // lin_coeff = I - dt * (VAV)
  lin_coeff->Update(
    1.0, *VAV, 0.0); /* bilinear terms are not added in lin_coeff!!!!!! */
  lin_coeff->Scale(-1.0 * (subdt_));
  for (int i = 0; i < m_; i++) {
    lin_coeff->SumIntoMyValue(i, i, 1.00);
  }
  if (debug_) {
    std::cout << "\n frobenius norm of A : " << A_->NormFrobenius() << "\n";
    printnormMV(*lin_coeff, 2, "norm of y-jac cols:");
    printnormMV(*VAV, 2, "norm of VAV cols:");
    printnormMV(*Vnew, 2, "norm of V cols:");
  }
  if (debug_) {
    std::cout << "\n In CoeffSolver::LinCoeff" << std::endl;
    double nrm1[Vnew->NumVectors()];
    lin_coeff->Norm1(&nrm1[0]);
    std::cout << "One norm of I-dt*VAV" << std::endl;
    for (int i = 0; i < Vnew->NumVectors(); ++i) {
      std::cout << nrm1[i] << "  ";
    }
    std::cout << std::endl;
    lin_coeff->MaxValue(&nrm1[0]);
    std::cout << "max value of I-dt*VAV" << std::endl;
    for (int i = 0; i < Vnew->NumVectors(); ++i) {
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
 *         Name:  y_jac
 *  Description:  computes the jacobian. It uses SACADO package computing for
 * nonlinear
 *                part of jacobian.
 * =====================================================================================
 */
void
CoeffSolver::y_jac()
{
  // jacBilinTerm();
  *jac = *lin_coeff;
  // jac->Update(subdt_,*JacNonLin,1.0);
}
/*
 * ===  FUNCTION
 * ======================================================================
 *         Name:  y_rhs
 *  Description:  calculates rhs of the system!!!
 *  x_ - the value obtained at previous time step
 *  x0_- the value to be obtained for current time step
 * =====================================================================================
 */
void
CoeffSolver::y_rhs()
{
  if (!useNwtn_) {
    NonLinRHS();
    rhs->PutScalar(0.0);
    rhs->Update(1.0, *x_, 1.0);
    // rhs=-x+lin_coeff*x0
    std::cout.precision(16);
    // rhs= x + dt*(V'<Vx,Vx>-V'E[<Vx,Vx>])
    rhs->Update(subdt_, *repEVyVy, 1.0);
    rhs->Update(-subdt_, *rhsNonLin, 1.0);
    rhs->Update(1.0, *VBdW, 1.0);
    if (debug_ && stochiter_ == 0) {
      double nrm;
      double nrm1[Vnew->NumVectors()], nrm2[y_->NumVectors()];
      printnormMV(*Vnew, 2, "norm of V: ");
      printnormMV(*rhsNonLin, 2, "norm of rhsnonlin: ");
      printnormMV(*zTrans_, 2, "norm of dW: ");
      printnormMV(*EyyOld_, 2, "norm of EyyOld_: ");
    }

    if (debug_ && stochiter_ == 0) {
    }
  } else {
    NonLinRHS();
    rhs->PutScalar(0.0);
    rhs->Multiply('N', 'N', 1.0, *lin_coeff, *x0_, 1.0); // rhs=-x+lin_coeff*x0
    rhs->Update(-1.0, *x_, 1.0);
    rhs->Update(-subdt_, *repEVyVy, 1.0); // rhs=-x+lin_coeff*x0+const_coeff
    // rhs=-x + lin_coeff*x0 + const_coeff + rhs_nonlin
    rhs->Update(subdt_, *rhsNonLin, 1.0);
    rhs->Update(-1.0, *VBdW, 1.0);
  }
}

Teuchos::RCP<Epetra_Map>
CoeffSolver::get_x_map()
{
  return map_x_;
}

void
CoeffSolver::set_x(const Teuchos::RCP<Epetra_MultiVector>& x0_temp)
{
  x0_ = x0_temp;
  *x_ = *x0_temp;
}

void
CoeffSolver::setDwiener(const Teuchos::RCP<Epetra_MultiVector>& z0_temp)
{
  dW = z0_temp;

  if ((0)) { /* TODO: Put it in debug_ or testing_*/
    itrtr += 1;
    std::string flnm = "wnr" + Teuchos::toString(itrtr) + ".mm";
    std::cout << "\n filename is: " << flnm << "\n";
    auto tmpmv = StochIO::readMV(flnm, *map_z_);
    *dW = *tmpmv;
  } else {
    dW->Scale(sqrt(subdt_));
  }
}

Teuchos::RCP<Epetra_Map>
CoeffSolver::get_f_map()
{
  return map_x_;
}

Teuchos::RCP<Epetra_MultiVector>
CoeffSolver::get_x_init()
{
  return x0_;
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
  dx_->ExtractView(&Viewx, &m_);
  rhs->ExtractView(&ViewRHS, &m_);
  DenseDx_ =
    Epetra_SerialDenseMatrix(Epetra_DataAccess::Copy, Viewx, m_, m_, MyLDA);
  DenseRHS =
    Epetra_SerialDenseMatrix(Epetra_DataAccess::Copy, ViewRHS, m_, m_, MyLDA);

  y_prob.SetVectors(DenseDx_, DenseRHS);
  y_prob.Solve();
  if (debug_) {
    y_prob.Matrix()->Print(std::cout << "inverted matrix \n");
    Yjac->Print(std::cout << "Yjac \n");
    std::cout << "\n norm one Yjac: " << Yjac->NormOne() << std::endl;
    Epetra_SerialDenseMatrix res;
    res = DenseRHS;
    res.Multiply('N', 'N', 1.0, *Yjac, DenseDx_, -1.0);
    //res.Print(std::cout << "res \n");
    std::cout << "\n If Solved : " << y_prob.Solved() << std::endl;
    std::cout << "\n Debug_ value : " << debug_<< std::endl;
    getchar();
  }

  for (int i = 0; i < MyLDA; i++)
    (*dx_)[i] = DenseDx_[i];
}

void
CoeffSolver::StochasticIterations()
{
  subdt_ = *dt_ / numSubTimeStep;
  double time, LocTime;
  Epetra_Time timer1(Vnew->Comm());

  computeEVyVy();

  if (test_) {
    LocTime = timer1.ElapsedTime();
    Comm_->SumAll(&LocTime, &time, 1);
    if (Comm_->MyPID() == 0) {
      std::cout << "\n*******************************************\n";
      std::cout << "average time per proc. for HBilinV: "
                << (time) / (A_->Comm().NumProc()) << " sec\n";
      std::cout << "*******************************************\n\n";
    }
  }
  Epetra_Time timer2(timer1);

  LinCoeff();

  if (test_) {
    LocTime = timer2.ElapsedTime();
    Comm_->SumAll(&LocTime, &time, 1);
    if (Comm_->MyPID() == 0) {
      std::cout << "\n*******************************************\n";
      std::cout << "average time per proc. for LinCoeff: "
                << (time) / (A_->Comm().NumProc()) << " sec\n";
      std::cout << "*******************************************\n\n";
    }
  }
  Epetra_Time timer3(timer1);
  y_jac();

  if (debug_) {
    std::cout.precision(18);
    std::cout << "\n frob-norm of det-jac: " << A_->NormFrobenius();
    printnormMV(*jac, 2, "2-norm of y-jac");
  }
  Epetra_MultiVector JAC(*jac);
  JAC = *jac;

  if (test_) {
    LocTime = timer3.ElapsedTime();
    Comm_->SumAll(&LocTime, &time, 1);
    if (Comm_->MyPID() == 0) {
      std::cout << "\n*******************************************\n";
      std::cout << "average time per proc. for y-jac() "
                << (time) / (A_->Comm().NumProc()) << " sec\n";
      std::cout << "*******************************************\n\n";
    }
  }
  Epetra_Time timer4(timer1);
  y_prob.SetMatrix(*Yjac);
  if (debug_) {
    printnormMV(*Vnew, 1, "1-norm of V in CoeffSolver::StochasticIterations():");
    std::cout << "One Norm of Yjac = " << Yjac->NormOne() << std::endl;
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
    Epetra_DataAccess::Copy, jac->Map(), valptr, m_, m_);
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

  VB->Multiply('T', 'N', 1.0, *Vnew, *B, 0.0);
  Epetra_Time timer5(timer1);
  for (int substep = 1; substep <= numSubTimeStep; substep++) {
    double avrg[zTrans_->NumVectors()];
    Epetra_MultiVector avgMV(*zTrans_);
    for (int i = 0; i < B->NumVectors(); i++) {
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
    CreateLocMultiVec(type);
    set_x(y_);
    setDwiener(z_);
    /********* stoch_coeff = VB*dW ****************/
    VBdW->Multiply('N', 'N', 1.0, *VB, *dW, 0.0);
    if (!useNwtn_) {
      dx_->PutScalar(0.0);
      y_rhs();
      dx_->Multiply('N', 'N', 1.0, inv_jac, *rhs, 0.0);
      x0_->Update(1.0, *dx_, 0.0);
    } else
      Newton();
    if (debug_ && stochiter_ == 0) {
      double nrm1[Vnew->NumVectors()];
      inv_jac.Norm2(&nrm1[0]);
      std::cout << "\n norm of inv(Yjac)\n";
      for (int i = 0; i < Vnew->NumVectors(); ++i) {
        std::cout << nrm1[i] << "  ";
      }
      std::cout << std::endl;

      std::cout << "\n frob. norm of inv(y_jac)  = " << Yjac->NormOne()
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
      EVyVy->Update((*Eyy)[j][i], (*((*Hn)(i * m_ + j))), 1.0);
    }
  }
}

void
CoeffSolver::computeCrossVariance()
{
  CreateDistTransMultivec(); /* Create zTrans and yTrans which have distributed
                                map */
  int tmp = Eyy->Multiply(
    'T', 'N', 1.0 / NumStochIter_, *yTrans_, *yTrans_, 0.0); /* Calculate yy' */
}

void
CoeffSolver::computeEyyTyT()
{
  double* ViewYY;
  double* ViewLocExpyyy;
  YY->ExtractView(&ViewYY, &rszyy);
  LocEyyy->ExtractView(&ViewLocExpyyy, &m_);
  Epetra_Vector tmpy(*(*y_)(0));
  int err1;
  for (int i = 0; i < MyLDA; i++) {
    tmpy = (*(*y_)(i));
    err1 = LocEyyy->Multiply('N', 'T', 1.0, tmpy, tmpy, 0.0);
    for (int j = 0; j < YY->MyLength(); j++) {
      ViewYY[i * YY->MyLength() + j] = ViewLocExpyyy[j];
    }
  }
  double cpyEyyy[EyyTyT->NumVectors()*EyyTyT->MyLength()];
  const int sz=EyyTyT->NumVectors()*EyyTyT->MyLength();
  double globSumEyyy[EyyTyT->NumVectors()*EyyTyT->MyLength()];
  EyyTyT->Multiply('N', 'T', 1.0 / NumStochIter_, *YY, *y_, 0.0);

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
  EVyVyy->Multiply('N', 'N', 1.0, *Hn, *EyyTyT, 0.0);
}
// ************************************************************************
void
CoeffSolver::computeEDEyy()
// ************************************************************************
{
  Epetra_MultiVector Inv_Exp_yy(Eyy->Map(), Eyy->NumVectors());
  SymMatPseudoInverse(*Eyy, Inv_Exp_yy); /* compute pseudo-inverse */
  Epetra_MultiVector zyT(*EzyDEyy), WzyT(*Vnew);
  /* Calculate zy' */
  zyT.Multiply('T', 'N', 1.0 / NumStochIter_, *zTrans_, *yTrans_, 0.0);
  /* Calculate W*E[zy'] */
  WzyT.Multiply('N', 'N', numSubTimeStep, *B, zyT, 0.0);
  Epetra_MultiVector Wzy_P_EVyVyy(*Vnew);
  /* Calculate W*E[zy']+E[<Vy,Vy>] */
  Wzy_P_EVyVyy.Update(*dt_, *EVyVyy, 1.0, WzyT, 0.0);
  /* Calculate (W*E[zy']+E[<Vy,Vy>])/E[yy'] */
  EDEyy->Multiply('N', 'N', 1.0, Wzy_P_EVyVyy, Inv_Exp_yy, 0.0);

  if (debug_) {
    Eyy->Print(std::cout << "Eyy \n");
    Inv_Exp_yy.Print(std::cout << "Inv_Exp_yy \n");
    Epetra_MultiVector eye(*Eyy);
    eye.Multiply('N', 'N', 1.0, Inv_Exp_yy, *Eyy, 0.0);
    eye.Print(std::cout << "Eyy * Inv(Eyy) \n");
    printnormMV(*Vnew, 2, "norm of V : ");
    printnormMV(*yTrans_, 2, "norm of Y^T : ");
    printnormMV(*B, 2, "norm of W : ");
    printnormMV(*VB, 2, "norm of VW : ");
    printTransNormMV(*rhs, 2, "norm of rhs^T: ");
    printTransNormMV(*dW, 2, "norm of dw^T: ");
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
  CreateDistTransMultivec();
  double A[m_][m_];
  int LWORK = 1 + 5 * m_ * m_, intdum = 1;
  double WORK[LWORK], W[m_], U[m_][m_], Vdum;
  char COMPZ = 'V', UPLO = 'U';
  int INFO;
  char JOBU = 'A', JOBVT = 'O';
  Epetra_LAPACK lapack;
  Epetra_MultiVector yy(*map_x_, m_);

  yy.Multiply('T', 'N', 1.0, *yTrans_, *yTrans_, 0.0);
  yy.ExtractCopy(&A[0][0], m_);
  lapack.GESVD(JOBU,
               JOBVT,
               m_,
               m_,
               &A[0][0],
               m_,
               &W[0],
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
      std::cout << "  " << sqrt(W[i] / float(NumStochIter_)) << "\n";
      DeviationFile << sqrt(W[i] / float(NumStochIter_)) << "\n";
    }
    std::cout << "\nVariances:\n";
    std::ofstream VarianceFile("VarFile.mm", ios::out);
    VarianceFile << "%%MatrixMarket matrix array real general\n"
                 << m_ << " " << 1 << "\n";
    VarianceFile << scientific;
    VarianceFile.precision(16);
    for (int i = 0; i < m_; i++) {
      std::cout << "  " << W[i] / NumStochIter_ << "\n";
      VarianceFile << W[i] / NumStochIter_ << "\n";
    }
  }
  int err;
  second_mmnt.PutScalar(0.0);
  for (int i = 0; i < Vnew->NumVectors(); i++) {
    for (int j = 0; j < Vnew->NumVectors(); j++) {
      second_mmnt.Multiply(
        yy[i][j] / NumStochIter_, *(*Vnew)(i), *(*Vnew)(j), 1.0);
    }
  }
}

void
CoeffSolver::Newton()
{
  isConverged_ = false;
  dx_->PutScalar(0.0);
  // y_jac();
  y_rhs();
  rhs->Scale(-1.0);
  std::vector<double> normsRHS(rhs->NumVectors());
  rhs->Norm2(normsRHS.data());
  NormRHS_ = *std::max_element(normsRHS.begin(), normsRHS.end());
  //std::cout << "\n norm of rhs = " << NormRHS_ << std::endl;
  //y_jac();
  //y_prob.SetMatrix(*Yjac); // 2. Give the updated matrix to the SVD solver
  //y_prob.Factor();
  INFO("Newton:      norm: " << NormRHStest_ );
  for (iter_ = 0; iter_ != maxNumIterations_; ++iter_) {
    Solve();
    if (debug_)
      x0_->Update(1.0, *dx_, 0.0);
    else {
      x0_->Update(1.0, *dx_, 1.0);
      y_rhs();
      rhs->Scale(-1.0);
    }
    std::vector<double> normsTest(rhs->NumVectors());
    rhs->Norm2(normsTest.data());
    NormRHStest_ = *std::max_element(normsTest.begin(), normsTest.end());
    if (debug_) {
      std::cout << "\n RHS NORM = " << NormRHStest_ << std::endl;
      Epetra_Vector Err(*(*y_)(0));
      Err.Multiply('N', 'N', 1.0, *jac, *x0_, 0.0);
      Err.Update(1.0, *rhs, 1.0);
      double nrm;
      Err.Norm2(&nrm);
      std::cout << "\n norm of err = jac*x0-rhs = " << nrm << std::endl;
    } else {
      if (NormRHStest_ < toleranceRHS_) {
        DEBUG("Success...");
        break;
      }
      if (backTracking_ and (NormRHS_ < NormRHStest_))
        RunBackTracking();
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
CoeffSolver::RunBackTracking()
{
  // Initialize reduction with -1/2
  double reduction = -1.0 / 2;
  for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_) {
    if (NormRHStest_ < NormRHS_) {
      DEBUG("Success...");
      break;
    }
    // Apply reduction to the state vector
    x0_->Update(reduction, *dx_, 1.0);
    y_rhs();
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
CoeffSolver::printTransNormMV(Epetra_MultiVector& mv, int normType, const std::string& str)
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

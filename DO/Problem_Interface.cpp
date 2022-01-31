// ----------   Includes   ---
// ------- #include <iostream>
#include "Problem_Interface.hpp"

#include <EpetraExt_MatrixMatrix.h>

#include "Amesos_Klu.h"
#include "Amesos_Lapack.h"
#include "Amesos_Scalapack.h"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AztecOO.h"
#include "BelosSolverFactory.hpp"
#include "BelosSolverManager.hpp"
#include "BelosTypes.hpp"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_LAPACK.h"
#include "HYMLS_MatrixUtils.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#if need_locaInterface == 1
#include "FVM_model_interface.h"
#endif
using namespace std;
// using Teuchos::RCP;
// using Teuchos::rcp; // Save some typing
//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(Teuchos::RCP<Epetra_CrsMatrix> A,
                                     Teuchos::RCP<Teuchos::ParameterList> SolverParams,
                                     Teuchos::RCP<Epetra_Vector> udet, int &m, double *t, double &dt,
                                     Teuchos::RCP<Epetra_MultiVector> W, Teuchos::RCP<Epetra_CrsMatrix> mass,
                                     int iter)
    : // **************************************************************************************
      A_(A), SolverParams_(SolverParams), udet_(udet), m_(m), t_(t), dt_(dt), W_(W), stochiter(iter), mass_(mass),
      v_prob(), v_solve()
{
  MyPID = A_->Comm().MyPID();
  InitFile_ = SolverParams_->get("StochBasisFile", "v2.mm");
  solver_type = SolverParams_->get("Solver Package", "any");
  test_ = SolverParams_->get("Class Testing", false);
  debug_ = SolverParams_->get("Class Debugging", false);
  DbgLvl_ = SolverParams_->get("Debug Level", 0);
  Teuchos::ParameterList &SolverSublist = SolverParams_->sublist(solver_type);
  std::string solver_name = SolverSublist.get("Solver name", "any");
  const Teuchos::RCP<Teuchos::ParameterList> SolParams = Teuchos::getParametersFromXmlFile(solver_name + ".xml");
  if (MyPID == 0)
    cout << solver_type << " solver(" << solver_name << ") has been chosen as stoch basis solver";
  // Set up the LU factorization for M-dt*A
  std::flush(cout << "\n no. of global rows in jacobian " << A_->NumMyRows());
  Epetra_Map rowmap(A_->RowMap());
  V = Teuchos::rcp(new Epetra_MultiVector(rowmap, m_));
  map_ = Teuchos::rcp(new Epetra_Map(m_, 0, A_->Comm()));
  locmap_ = Teuchos::rcp(new Epetra_LocalMap(m_, 0, A_->Comm()));
  y_map = Teuchos::rcp(new Epetra_Map(W_->NumVectors(), 0, A_->Comm()));
  iteration = 1;
  expv3 = Teuchos::rcp(new Epetra_MultiVector(*V));
  eye = Teuchos::rcp(new Epetra_MultiVector(*map_, m_));
  Exp_zy = Teuchos::rcp(new Epetra_MultiVector(*y_map, m_));
  Exp_yy = Teuchos::rcp(new Epetra_MultiVector(*map_, m_));
  Exp_yyy = Teuchos::rcp(new Epetra_MultiVector(*map_, m_));
  exp_yy_inv = Teuchos::rcp(new Epetra_MultiVector(*map_, m_));
  RHS_block_1 = Teuchos::rcp(new Epetra_MultiVector(A_->RowMap(), m_));
  ExpDExpyy = rcp(new Epetra_MultiVector(A->RowMap(), m_));
  LHS_block_1_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, A_->RowMap(), 3));
  Rvec = Teuchos::rcp(new Epetra_MultiVector(*locmap_, m_));
  Rvec->ExtractView(&r_val, &m_);
  R = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double>(Teuchos::DataAccess::View, r_val, m_, m_, m_));
  /****** LHS_block_1 = M-dt*A *********/
  LHS_block_1 = LHS_block_1_.get();
  // EpetraExt::MatrixMatrix::Add(*mass_, false, 1.0, *LHS_block_1_, 1.0);
  // EpetraExt::MatrixMatrix::Add(*A_, false, -dt_, *LHS_block_1_, 1.0);
  // LHS_block_1_->FillComplete();
  EpetraExt::MatrixMatrix::Add(*mass_, false, 1.0, *A_, false, -(dt_), LHS_block_1);
  LHS_block_1->FillComplete();
  v_prob = Teuchos::rcp(new Epetra_LinearProblem);
  if (solver_type == "Amesos")
  {
    v_prob->SetOperator(LHS_block_1_.get());
    Amesos Factory;
    v_solve = Teuchos::rcp(Factory.Create(solver_name.c_str(), *v_prob), false);
    // v_solve = Teuchos::rcp(new Amesos_Klu(*v_prob));
    v_solve->SetUseTranspose(false);
    v_solve->SetParameters(*SolParams);
    Epetra_Time timer(A_->Comm());
    v_solve->SymbolicFactorization();
    double time;
    if (test_)
    {
      double LocTime = timer.ElapsedTime();
      A_->Comm().SumAll(&LocTime, &time, 1);
      if (A_->Comm().MyPID() == 0)
      {
        cout << "\n*******************************************\n";
        cout << "Average time per proc. for factorization"
             << " = " << (time) / (A_->Comm().NumProc()) << " sec\n";
        cout << "*******************************************\n" << std::endl;
      }
    }
  }
  else if (solver_type == "Amesos2")
  {
    amesos2_solver = Amesos2::create<MAT,MV>(solver_name,LHS_block_1_);
    amesos2_solver->setParameters(SolParams); 
    Epetra_Time timer(A_->Comm());
    amesos2_solver->symbolicFactorization();
    double time;
    if (test_)
    {
      double LocTime = timer.ElapsedTime();
      A_->Comm().SumAll(&LocTime, &time, 1);
      if (A_->Comm().MyPID() == 0)
      {
        cout << "\n*******************************************\n";
        cout << "Average time per proc. for factorization"
             << " = " << (time) / (A_->Comm().NumProc()) << " sec\n";
        cout << "*******************************************\n" << std::endl;
      }
    }
  }
  else
  {
    // allocates an IFPACK factory. No data is associated
    // to this object (only method Create()).
    Ifpack Factory;

    // create the preconditioner. For valid precType values,
    // please check the documentation
    std::string precType = SolParams->get("Ifpack Preconditioner Name", "any");
    int OverlapLevel = SolParams->get("Overlap Level", 0);
    // it is ignored.

    prec = Teuchos::rcp(Factory.Create(precType, &*LHS_block_1_, OverlapLevel));
    assert(prec != Teuchos::null);
    
    Teuchos::ParameterList ifpackList=SolParams->sublist(precType);

    (prec->SetParameters(ifpackList));

    // initialize the preconditioner. At this point the matrix must
    // have been FillComplete()'d, but actual values are ignored.
    (prec->Initialize());

    // Create the Belos preconditioned operator from the Ifpack preconditioner.
    belosPrec = rcp(new Belos::EpetraPrecOp(prec));
    bool success = true;
    bool leftprec = true; // left preconditioning or right.
    int numrhs= m_;
    std::string solverType = SolParams->get("Solver","any");
    ParameterList belosList = SolParams->sublist(solverType);
    MT tol = belosList.get("Convergence Tolerance",1.0);
    belosList.set("Convergence Tolerance", tol); // Relative convergence tolerance requested
    if (numrhs > 1)
    {
      belosList.set("Show Maximum Residual Norm Only",
                    true); // Show only the maximum residual norm
    }
    bool verbose = belosList.get("Verbose",false) ;
    if (verbose)
    {
      belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
      int frequency = belosList.get("output frequency",0);   // frequency of status test output.
      if (frequency > 0)
        belosList.set("Output Frequency", frequency);
    }
    else
      belosList.set("Verbosity", Belos::Errors + Belos::Warnings);

    problem = rcp(new Belos::LinearProblem<double, MV, OP>());
    problem->setOperator(LHS_block_1_); // set the operator for belos
    
    if (leftprec)
    {
      problem->setLeftPrec(belosPrec);
    }
    else
    {
      problem->setRightPrec(belosPrec);
    }
    Teuchos::writeParameterListToXmlFile(*SolParams,"storedBelosList.xml"); 
    // Create an iterative solver manager.
    
    Belos::SolverFactory<double,MV,OP> belosFactory;
    v_solve_iter = belosFactory.create(solverType, rcp(&belosList,false));
    //v_solve_iter = rcp(new Belos::BlockGmresSolMgr<double, MV, OP>(problem, rcp(&belosList, false)));
  }
  double one = 1.0;
  for (int i = 0; i < m_; i++)
    eye->ReplaceGlobalValue(i, i, one);
  // Currently the forcing is linear, so we compute it only once!!!!
  int n;
  n = A->NumMyRows();

#if need_locaInterface == 1
//  double* W_ptr;
//  int ierr;
//  for (int col = 0; col < W_->NumVectors(); col++) {
//    (*W_)(col)->ExtractView(&W_ptr);
//    model_stoch_frc(&n, W_ptr, &col, &ierr);
//  }
#endif
}
// *************************************************************************
Problem_Interface::~Problem_Interface()
// ************************************************************************
{
}
// ************************************************************************
int Problem_Interface::SolveV(Epetra_MultiVector& LHS, Epetra_MultiVector& RHS, std::string label)
// ************************************************************************
{
  double time, LocTime;
  Epetra_Time timer1(A_->Comm());
  if (solver_type == "Amesos")
  {
    v_prob->SetLHS(&LHS);
    v_prob->SetRHS(&RHS);
    int err = (v_solve->Solve());
    if (debug_)
    {
      v_solve->PrintStatus();
    }
  }
  else if (solver_type == "Amesos2")
  {
    amesos2_solver->solve(&LHS,&RHS);
  }
  else
  {
    // Builds the preconditioners, by looking for the values of
    // the matrix.
    //v_solve_iter->reset(Belos::Problem);
    problem->setOperator(LHS_block_1_);
   // problem->setProblem(Teuchos::rcp(LHS), Teuchos::rcp(RHS));
    problem->setLHS(Teuchos::rcpFromRef(LHS));
    problem->setRHS(Teuchos::rcpFromRef(RHS));
    problem->setProblem();
    v_solve_iter->setProblem(problem);
    v_solve_iter->solve(); // TODO: Parameteres Should come from xml file!!!!

    if (v_solve_iter->isLOADetected())
      std::flush(std::cout << "\n Loss of accuracy detected!!! \n" 
	  <<"achieved tol = " << v_solve_iter->achievedTol());

    if (debug_ && MyPID == 0)
      std::cout << "\nIterative solution\n"
                << "Solver performed " << v_solve_iter->getNumIters() << " iterations." << std::endl
                << "Achieved tolerance = " << v_solve_iter->achievedTol() << std::endl;
  }
  if (debug_ && DbgLvl_ % 2 == 0)
  {
    Epetra_MultiVector rh(LHS);
    LHS_block_1_->Multiply(false, LHS, rh);
    rh.Update(-1.0, RHS, 1.0);
    double nrm1[rh.NumVectors()];
    rh.Norm2(&nrm1[0]);
    std::cout << "\nNorm of residual (" << label << ")" << std::endl;
    for (int i = 0; i < rh.NumVectors(); i++)
    {
      cout << "  " << nrm1[i];
    }
  }
  if (test_)
  {
    LocTime = timer1.ElapsedTime();
    A_->Comm().SumAll(&LocTime, &time, 1);
    if (A_->Comm().MyPID() == 0)
    {
      cout << "\n*******************************************\n";
      cout << label << " = " << (time) / (A_->Comm().NumProc()) << " sec\n";
      cout << "*******************************************\n" << std::endl;
    }
  }
  return 0;
}
// ************************************************************************
bool Problem_Interface::computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac)
// ************************************************************************
{
  Jac = Teuchos::dyn_cast<Epetra_Operator>(getJacobian());
  return true;
}
//**************************************************************************************
bool Problem_Interface::computeShiftedMatrix(double alpha, double beta, const Epetra_Vector &x, Epetra_Operator &A)
//*************************************************************************
{
  return false;
}
// **************************************************************************************
void Problem_Interface::init_v(double t)
// **************************************************************************************
{
  // Initialize V such that if intial time=0 and W={w_i}
  // then  V=[w_i,randomvec_i+1,...,randomvec_n]
  // else read from the file.
  if (InitFile_ == "None")
  {
    V->Random();
//    int myln = V->GlobalLength();
//    V->Scale(1.0/myln);
/*     for (int i = 0; i < V->NumVectors(); i++)
 *     {
 *       for (int j = 0; j < myln; j+=2)
 *       {
 * 	(*((*V)(i)))[j]=cos(j*M_PI* myln/2);
 * 	(*((*V)(i)))[j+1]=1.0;
 *       }
 * 
 *     }
 */
/*     if (abs(stchFrcStren_)>=10e-12)
 *     {
 *       for (int i = 0; i < std::min(W_->NumVectors(),V->NumVectors()) ; i++)
 *         *((*V)(i))=*((*W_)(i));
 *     }
 */

#if need_locaInterface == 1
    Epetra_Vector massDiag(V->Map());
    mass_->ExtractDiagonalCopy(massDiag);
    for (int i = 0; i < massDiag.MyLength(); i++)
    {
      if (abs(massDiag[i]) < 10e-10)
      {
        for (int j = 0; j < V->NumVectors(); j++)
        {
          (*V)[j][i] = 0.0;
        }
      }
    }
#endif
    //for (int i = 0; i < W_->NumVectors(); i++)
      //mass_->Multiply(false, *((*W_)(i)), *((*V)(i)));
    // if (solver_type == "direct") {
    //  v_prob->SetLHS(((*V)(0)));
    //  v_prob->SetRHS(((*V)(0)));
    //  v_solve->Solve();
    //} else {
    //  v_solve_iter->SetLHS(((*V)(0)));
    //  v_solve_iter->SetRHS(((*V)(0)));
    //  v_solve_iter->Iterate(
    //    MaxIter_, Tol_); // TODO: Parameteres Should come from xml file!!!!
    //}
    // for (int i = 1; i < V->NumVectors(); i++) {
    //  mass_->Multiply(false, *((*V)(i - 1)), *((*V)(i)));
    //  if (solver_type == "direct") {
    //    v_prob->SetLHS(((*V)(i)));
    //    v_prob->SetRHS(((*V)(i)));
    //    v_solve->Solve();
    //  } else {
    //    v_solve_iter->SetLHS(((*V)(i)));
    //    v_solve_iter->SetRHS(((*V)(i)));
    //    v_solve_iter->Iterate(
    //      MaxIter_, Tol_); // TODO: Parameteres Should come from xml file!!!!
    //  }
    //}
    MOrth(); // V is global variable for this class. It need not to be
    // passed.
    if (debug_ && DbgLvl_ % 3 == 0)
    {
      double nrm1[V->NumVectors()];
      V->Norm2(&nrm1[0]);
      std::cout << "\nNorm of V" << std::endl;
      for (int i = 0; i < V->NumVectors(); i++)
      {
        cout << "  " << nrm1[i];
      }
    }
  }
  else
  {
    Epetra_MultiVector *vt;
    EpetraExt::MatrixMarketFileToMultiVector(InitFile_.c_str(), (V->Map()), vt);
    V = Teuchos::rcp(vt);
    if (debug_)
    {
      printnormMV(*V, 2, "initial norm of V without m-orthogonalization:");
    }
    MOrth(); // V is global variable for this class. It need not to be passed.
    // It is handy not to assume a certain orthogonalization. So we should scale
    // y's after they are read too: Ry -> y
  }
}
// *************************************************************************
void Problem_Interface::MOrth()
// *************************************************************************
{
  Epetra_MultiVector Mmv(*V);
  double time, LocTime;
  Epetra_Time timer2(A_->Comm());
  R->putScalar(0.0);
  typedef Epetra_MultiVector mv;
  typedef Epetra_Operator OP;
  // mmv = Teuchos::rcp(new Epetra_MultiVector(*V));
  mass_->Multiply(false, *V, Mmv);
  A_->Comm().Barrier();
  Belos::IMGSOrthoManager<double, mv, OP> orthman;
  orthman.setOp(mass_);
  int rank = orthman.normalize(*V, Teuchos::rcpFromRef(Mmv), R);
  double orthoerr = orthman.orthonormError(*V);
  //  Epetra_MultiVector MV(*((*V)(0)));
  //  double nrm;
  //  mass_->Multiply(false, *((*V)(0)), MV);
  //  MV.Dot(*((*V)(0)), &nrm);
  //  (*R)(0, 0) = sqrt(nrm);
  //  (*V)(0)->Scale(1 / (*R)(0, 0));

  //  if (debug_)
  //    cout << "\nNORM\n";

  //  for (int i = 1; i < m_; i++) {
  //    for (int j = 0; j < i; j++) {
  //      mass_->Multiply(false, *((*V)(i)), MV);
  //      MV.Dot(*((*V)(j)), &((*R)(j, i)));
  //      (*V)(i)->Update(-(*R)(j, i), *((*V)(j)), 1.0);
  //    }
  //    mass_->Multiply(false, *((*V)(i)), MV);
  //    MV.Dot(*((*V)(i)), &nrm);
  //    (*R)(i, i) = sqrt(nrm);
  //    (*V)(i)->Scale(1 / (*R)(i, i));
  //  }

  if (test_)
  {
    LocTime = timer2.ElapsedTime();
    A_->Comm().SumAll(&LocTime, &time, 1);
    if (A_->Comm().MyPID() == 0)
    {
      cout << "\n*******************************************\n";
      cout << " Avergage time per proc. for Orthogonalization = " << (time) / (A_->Comm().NumProc()) << " sec\n";
      cout << "*******************************************\n" << std::endl;
    }
  }
  if (debug_ && DbgLvl_ % 5 == 0)
  {
    typedef Anasazi::MultiVecTraits<double, mv> MVT;
    Teuchos::RCP<mv> mmv, tem, temv, temv2;
    temv = Teuchos::rcp(new mv(*V));
    double nrm1[V->NumVectors()];
    // test M-orthogonality
    mass_->Multiply(false, *V, *temv);
    temv->Dot(*V, &nrm1[0]);
    // mest constraint satisfaction
    cout << "\n Norm of V'MV\n";
    for (int i = 0; i < V->NumVectors(); i++)
      std::cout << "  " << nrm1[i];
  }
}
// **************************************************************************************
int Problem_Interface::v_stoch_init(Epetra_MultiVector *Vold)
// **************************************************************************************
{
  Epetra_MultiVector X1(*RHS_block_1);
  Epetra_MultiVector X2(*V);
  std::string label;
  /*****SOlve X1=J^-1*RHS_block_1*************/
  if (debug_)
    label = "J*X1-RHS_block_1";
  else
    label = "average time per proc. for soving X1=J^-1*RHS_block_1";
  if (solver_type == "Amesos")
  {
    // std::flush(std::cout<<"\n HELLO1!! \n");
    // std::flush(std::cout<<"\n HELLO1!!
    // "<<v_solve->NumericFactorization()<<"\n");
    AMESOS_CHK_ERR(v_solve->NumericFactorization());
  }
  else if (solver_type == "Amesos2")
  {
    // std::flush(std::cout<<"\n HELLO1!! \n");
    // std::flush(std::cout<<"\n HELLO1!!
    // "<<v_solve->NumericFactorization()<<"\n");
    (amesos2_solver->numericFactorization());
  }
  else
  {
    IFPACK_CHK_ERR(prec->Compute());
  }
  SolveV(X1, *RHS_block_1, label);
  /*******************************
    SOlve X2=J^-1*(MV)
   ********************************/
  Epetra_MultiVector MV(*V);
  if (debug_)
    label = "J*X2-MV";
  else
    label = "average time per proc. for soving X2=J^-1*MV";
  mass_->Multiply(false, *V, MV);
  SolveV(X2, MV, label);
  /*******************************
    SOlve Y=(MV'*X2)^-1*(MV'*X1)
   ********************************/
  /***** TODO: THE MATRIX BELOW SHOULD BE MADE DENSE **/
  int numvecV = V->NumVectors();
  double one = 1.0;
  Epetra_LocalMap MVtX_map(numvecV, 0, V->Comm());
  Epetra_CrsMatrix MVtX2(Epetra_DataAccess::Copy, MVtX_map, numvecV);
  for (int i = 0; i < numvecV; i++)
  {
    for (int j = 0; j < numvecV; j++)
    {
      MVtX2.InsertGlobalValues(i, 1, &one, &j);
    }
  }
  MVtX2.FillComplete();
  Epetra_MultiVector MVtX2_temp(MVtX_map, numvecV);
  Epetra_MultiVector MVtX1(MVtX2_temp), Y_sol(MVtX2_temp);
  int err1 = MVtX2_temp.Multiply('T', 'N', 1.0, MV, X2, 0.0);
  int err2 = MVtX1.Multiply('T', 'N', 1.0, MV, X1, 0.0);

  double jacval;
  for (int i = 0; i < numvecV; i++)
  {
    for (int j = 0; j < numvecV; j++)
    {
      jacval = MVtX2_temp.operator()(j)->operator[](i);
      MVtX2.ReplaceGlobalValues(i, 1, &jacval, &j);
    }
  }
  Epetra_LinearProblem v_prob2;
  Amesos_Lapack v_solve2(v_prob2);
  v_prob2.SetOperator(&MVtX2);
  v_prob2.SetLHS(&Y_sol);
  v_prob2.SetRHS(&MVtX1);
  v_solve2.SetUseTranspose(false);
  (v_solve2.SymbolicFactorization());
  (v_solve2.NumericFactorization());
  (v_solve2.Solve());
  if (debug_ && DbgLvl_ % 2 == 0)
  {
    Epetra_MultiVector rh(MVtX1);
    HYMLS::MatrixUtils::Dump(MVtX1, "MVtX1");
    HYMLS::MatrixUtils::Dump(MVtX2, "MVtX2");
    HYMLS::MatrixUtils::Dump(Y_sol, "Y_sol");
    MVtX2.Multiply(false, Y_sol, rh);
    rh.Update(-1.0, MVtX1, 1.0);
    double nrm1[rh.NumVectors()];
    rh.Norm2(&nrm1[0]);
    std::cout << "\nNorm of residual (MVtX2*Y_sol-MVtX1)" << std::endl;
    for (int i = 0; i < rh.NumVectors(); i++)
    {
      cout << "  " << nrm1[i];
    }
  }

  // SOlve X=X1-(X2*Y)
  Epetra_MultiVector X_sol(X1);
  X_sol = X1;
  int errcode = X_sol.Multiply('N', 'N', -1.0, X2, Y_sol, 1.0);
  V->Update(1.0, *Vold, 0.0);
  V->Update(1.0, X_sol, 1.0);
  // Orthonormalize V
  MOrth();
  if (debug_ && DbgLvl_ % 3 == 0)
  {
    printnormMV(*V, 2, "norm of V after m-orthogonalization:");
  }
  return (EXIT_SUCCESS);
}
// **************************************************************************************
void Problem_Interface::computeBlocks(double dt)
// **************************************************************************************
{
  dt_ = dt;
  /****** LHS_block_1 = M-dt*J *********/
  // EpetraExt::MatrixMatrix::Add(*mass_, false, 1.0, *LHS_block_1_, 1.0);
  // EpetraExt::MatrixMatrix::Add(*A_, false, -dt_, *LHS_block_1_, 1.0);
  EpetraExt::MatrixMatrix::Add(*mass_, false, 1.0, *A_, false, -(dt_), LHS_block_1);
  if (debug_ && DbgLvl_ % 7 == 0)
  {
    std::cout << "\nNorm of deterministic Jacobian A = " << A_->NormFrobenius() << std::endl;
  }
  Epetra_Map zeromap_(m_, 0, A_->Comm());
  RCP<Epetra_MultiVector> zeros;
  zeros = Teuchos::rcp(new Epetra_MultiVector(zeromap_, m_));
  zeros->PutScalar(0.0);
  int N = A_->NumGlobalRows();
  /***********************************************************************************/
  // A_ should be the copy of deterministic linear jacobian ....not pointer to
  // it !!!!!!!!!!!!!!!!!!!!!!!
  // LHS
  // |(M-dt*J) MV |
  // |            | =
  // | V'M      0 |
  // RHS_block_1 = Fv = dt*JV + (dt*E[<VY,VY>Y^T] + W*E[dw*Y^T])/E[YY^T])
  Teuchos::RCP<Epetra_CrsMatrix> Acopy = Teuchos::rcp(new Epetra_CrsMatrix(*A_));
  *Acopy = *A_;
  /*********** AV= A*V*****************/
  Epetra_MultiVector AV(*V);
  Acopy->Multiply(false, *V, AV);
  Epetra_MultiVector Fv(AV);

  /* ExpDExpyy= (numsubtimesteps*W*Exp[zy]+dt*Exp[<Vy,Vy>y])/Exp[yy] */
  Fv.Update(dt_, AV, 1.0, *ExpDExpyy, 0.0);
  if (debug_)
  {
    HYMLS::MatrixUtils::Dump(AV, "AV");
    HYMLS::MatrixUtils::Dump(Fv, "FV");
    HYMLS::MatrixUtils::Dump(*ExpDExpyy, "Expect");
    double *FVarray;
    int Ldim;
    Fv.ExtractView(&FVarray, &Ldim);
    double nrm = 0;
    for (int i = 1; i < Fv.MyLength(); i = i + 2)
    {
      nrm = nrm + abs(FVarray[i]) / dt_;
    }
    std::ofstream NormConstrFile("Constraint.txt", std::ios::out | std::ios::app);
    NormConstrFile << std::scientific;
    NormConstrFile.precision(16);
    NormConstrFile << nrm << "\n";
  }
  if (debug_ && DbgLvl_ % 11 == 0)
  {
    double nrm1[m_];
    AV.Norm2(&nrm1[0]);
    std::cout << "\n nrm of AV" << std::endl;
    for (int i = 0; i < V->NumVectors(); i++)
    {
      cout << "  " << nrm1[i];
    }
    ExpDExpyy->Norm2(&nrm1[0]);
    std::cout << "\n nrm of (dt*E[<Vy,Vy>y']+numsubtimestep*E[zy'])/E[yy']" << std::endl;
    for (int i = 0; i < V->NumVectors(); i++)
    {
      cout << "  " << nrm1[i];
    }
  }

  RHS_block_1->Update(1.0, Fv, 0.0);

  if (debug_ && DbgLvl_ % 11 == 0)
  {
    double nrm1[m_];
    RHS_block_1->Norm2(&nrm1[0]);
    std::cout << "\n  nrm of RHS_block_1" << std::endl;
    for (int i = 0; i < V->NumVectors(); i++)
    {
      cout << "  " << nrm1[i];
    }
  }
}
// ************************************************************************
void Problem_Interface::TransferNorm()
// ************************************************************************
{
  // y(:,:)=R*y(:,:)
  Epetra_MultiVector ycpy(*y);
  ycpy = *y;
  y->Multiply('N', 'N', 1.0, *Rvec, ycpy, 0.0);
}
void Problem_Interface::printnormMV(Epetra_MultiVector &mv, int normType, string str)
{
  double nrm1[mv.NumVectors()];
  if (normType == 1)
    mv.Norm1(&nrm1[0]);
  else if (normType == 2)
    mv.Norm2(&nrm1[0]);
  std::flush(std::cout << "\n" << str << std::endl);
  for (int i = 0; i < mv.NumVectors(); i++)
    std::flush(std::cout << "  " << nrm1[i]);
}

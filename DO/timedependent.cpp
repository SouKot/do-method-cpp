#include <iostream>
using std::cerr;
using std::endl;
#include <fstream>
#include <sstream>
using std::ofstream;
#include "AztecOO.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "HYMLS_HyperCube.hpp"
#include "HYMLS_Tools.hpp"
#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX.H"
#include "NOX_Epetra_LinearSystem_Belos.H"
#include "NOX_Epetra_LinearSystem_Hymls.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "globdefs.H"
#include <EpetraExt_MatrixMatrix.h>
#include <iomanip>
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_OperatorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorIn.h"
#include "Epetra_LinearProblem.h"
#include "FVM_model_interface.h"
#include "HYMLS_MatrixUtils.hpp"
#include "Interface.hpp"
#include "Problem_Interface.hpp"
#include "StochSys.hpp"
//#include "matplotlibcpp.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <sstream>
#ifndef INFO
#define INFO(s) std::cout << s << std::endl;
#endif
// Global Timers
Teuchos::RCP<Teuchos::Time> coeffTime =
  Teuchos::TimeMonitor::getNewTimer("my total time spent on Stoch. Coeff. part");
Teuchos::RCP<Teuchos::Time> coeffTime2 =
  Teuchos::TimeMonitor::getNewTimer("my Bilinear form computation time");
Teuchos::RCP<Teuchos::Time> basisTime =
  Teuchos::TimeMonitor::getNewTimer("my total time spent on Stoch. Basis part");
Teuchos::RCP<Teuchos::Time> meanTime =
  Teuchos::TimeMonitor::getNewTimer("my total time spent on Mean part");
  Teuchos::RCP<Teuchos::Time> meanSolverTime = Teuchos::TimeMonitor::getNewTimer("my time spent on mean solver");
Teuchos::RCP<Teuchos::Time> meanMassMatTime = Teuchos::TimeMonitor::getNewTimer("my time spent on mass matrix computation");
Teuchos::RCP<Teuchos::Time> meanRhsTime =
  Teuchos::TimeMonitor::getNewTimer("my total time spent on RHS computation");
Teuchos::RCP<Teuchos::Time> meanJacTime =
  Teuchos::TimeMonitor::getNewTimer("my total time spent on Jacobian computation");
using namespace std;
#if need_locaInterface == 0
// namespace plt = matplotlibcpp;
// void
// plot(const Epetra_Vector& v, int mmnt)
//{
//  std::cout << "\n HEELLOO" << std::endl;
//  double arr[v.GlobalLength()];
//  double localcpy[v.MyLength()];
//  v.ExtractCopy(localcpy);
//  v.Map().Comm().GatherAll(localcpy, arr, v.MyLength());
//  if (mmnt == 2) {
//    vector<double> vec(arr, arr + sizeof(arr) / sizeof(double));
//    plt::plot(vec);
//    plt::show();
//  }
//}
#endif
void printnormMV(Epetra_MultiVector &mv, int normType, string str) {
  double nrm1[mv.NumVectors()];
  if (normType == 1)
    mv.Norm1(&nrm1[0]);
  else if (normType == 2)
    mv.Norm2(&nrm1[0]);
  std::cout << "\n" << str << std::endl;
  for (int i = 0; i < mv.NumVectors(); i++)
    std::cout << "  " << nrm1[i];
}
const double r0dim = 6370.0e3;
const double udim = 0.1;
const double comp_time_scaling() {
#if need_locaInterface == 1
  return r0dim / udim;
#else
  return 1.0;
#endif
} // r0dim/udim;}  // time scale [s]
const double timesc = comp_time_scaling();

// a function for nice output of time data
string time_out(double t) {
  string unit;
  double val;

  // try to figure out which unit gives a result closest to 1
  double sec = t * timesc;
  double diff = std::abs(1.0 - sec);
  double diffmin = diff;
  val = sec;
  unit = "s";

  std::stringstream s;
  s << val << " [" << unit << "]";
  return s.str();
}

int main(int argc = 0, char *argv[] = NULL) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  bool status = true;

  ////////////////////////////////////////
  // Create the Parallel  Communicator  //
  ////////////////////////////////////////

#ifdef HAVE_MPI
  Teuchos::RCP<Epetra_MpiComm> Comm =
    Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  // HYMLS::HyperCube Topology;
  // Teuchos::RCP<Epetra_MpiComm> Comm = Topology.GetComm();
#else
  Teuchos::RCP<Epetra_SerialComm> Comm = rcp(new Epetra_SerialComm());
#endif

  Epetra_Time timer(*Comm);
  // Get process ID and total number of processes
  int MyPID = Comm->MyPID();
  int NumProc = Comm->NumProc();

  if (MyPID == 0) {
    std::ostringstream infofile;
    infofile << "info_" << MyPID << ".txt";
    std::string myinfo = infofile.str();
    const char* Myinfo = myinfo.c_str();
    cout << "info is written to " << Myinfo << endl;
    info = Teuchos::rcp(new std::ofstream(Myinfo));

  } else {
    info = Teuchos::rcp(new Teuchos::oblackholestream());
  }

  (*info) << std::setw(15) << std::setprecision(15);

  Teuchos::RCP<std::ostream> outstream;

  if (MyPID == 0) {
    outstream = Teuchos::rcp(new std::ofstream("timestepping.out"));
  } else {
    outstream = Teuchos::rcp(new Teuchos::oblackholestream());
  }
  // this allows us to use the DEBUG macro, the Tools::Out function etc.
  // HYMLS::Tools::InitializeIO_std(Comm, outstream);
  (*outstream) << std::setw(15) << std::setprecision(15);

  // HYMLS::Tools::out() << "FVM revision: "<<FVM_REVISION<<std::endl;
  // HYMLS::Tools::out() << "HYMLS revision:
  // "<<HYMLS::Tools::Revision()<<std::endl;

  //#ifdef HAVE_MPI
  // HYMLS::Tools::out() << Topology << std::endl;
  //#endif

  DEBUG("*********************************************")
  DEBUG("* Debugging output for process " << MyPID)
  DEBUG("* To prevent this file from being written,  ")
  DEBUG("* omit the -DDEBUGGING flag when compiling. ")
  DEBUG("*********************************************")

  INFO(*Comm);

  bool stat = true;

  try {

    ////////////////////////////////////////////////////////
    // Setup Parameter Lists                              //
    ////////////////////////////////////////////////////////

    // read parameters from files
    // we create one big parameter list and further down we will set some things
    // that are not straight-forwars in XML (verbosity etc.)
    Teuchos::Ptr<Teuchos::ParameterList> paramListptr =
      Teuchos::ptr(new Teuchos::ParameterList);

    // override default settings with parameters from user input file
    Teuchos::updateParametersFromXmlFile("params.xml", paramListptr);

    Teuchos::RCP<Teuchos::ParameterList> paramList = rcpFromPtr(paramListptr);

    // extract the final sublists:
#if need_locaInterface == 1
    // Get the model sublist
    Teuchos::ParameterList& modelList = paramList->sublist("Model");

    // Get the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

    Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
    printParams.set("MyPID", MyPID);
#ifdef OUTPUT_TO_FILE
    printParams.set("Output Stream", outstream);
    printParams.set("Error Stream", outstream);
#endif
    printParams.set("Output Process", 0);
    printParams.set("Output Information",0
#ifdef DEBUGGING
                   + NOX::Utils::Details + NOX::Utils::OuterIteration +
                      NOX::Utils::InnerIteration +
                      NOX::Utils::OuterIterationStatusTest +
                      NOX::Utils::LinearSolverDetails +
                      NOX::Utils::Debug + NOX::Utils::StepperParameters +
                      NOX::Utils::Warning + NOX::Utils::StepperDetails +
                      NOX::Utils::StepperIteration
#endif
		    );
    // Create the "Direction" sublist for the "Line Search Based" solver
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");

    // Create the "Line Search" sublist for the "Line Search Based" solver
    Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");

    // Create the "Direction" sublist for the "Line Search Based" solver
    Teuchos::ParameterList& newtParams = dirParams.sublist("Newton");

    // Create the "Linear Solver" sublist for the "Direction' sublist
    Teuchos::ParameterList& lsParams = newtParams.sublist("Linear Solver");
#endif

    // Get the "Time Stepping" sublist
    Teuchos::ParameterList& transParams = paramList->sublist("Time Stepping");
    double t_start = transParams.get("Start Time", 0.0) / timesc;
    double t_end = transParams.get("End Time", 1.0) / timesc;
    double dt = transParams.get("Step Size", 0.01) / timesc;
    double t = t_start;

#if need_locaInterface == 1
    string stepperScheme = transParams.get("Scheme", "Theta");
    bool gradual_startup = transParams.get("Gradual Start-Up", true);
    string startup_param = "None";
    double startup_rate = 0.0;
    double startup_value, startup_max;

    if (gradual_startup) {
      startup_param = transParams.get("Start-Up Parameter", "Wind Forcing");
      startup_rate = transParams.get("Start-Up Rate", 0.1);
      startup_max =
        modelList.sublist("Starting Parameters").get(startup_param, 0.0);
    }
    transParams.set("Output Stream", outstream);
    int nslow = 0;
    double dt_min =
      transParams.get("Minimum Step Size", (double)(dt / 10)) / timesc;
    double dt_max = transParams.get("Maximum Step Size", dt) / timesc;
    int maxsteps =
      transParams.get("Max Num Steps", (int)((t_end - t_start) / dt));
    double lteTol = transParams.get("Error Tolerance", 1.0e-3);
    string steadyMon = transParams.get("Steady-State Monitor", "Norm of F");
    double steadyTol = transParams.get("Steady-State Tolerance", 0.0);

    string stepChoice = transParams.get("Step Size Control", "Constant");

    // step size change factors for "Constant" mode
    double red_fac = transParams.get("Failed Step Reduction Factor", 0.5);
    double inc_fac = transParams.get("Successful Step Increase Factor", 1.5);

    int step = 1;
    double maxEffort = 0.9; // prevent increasing step size if more than
    // 90% of the maximum Newton steps were needed.
    double PsiMaxOld = -1.0;
#else
    double dtIncrmnt = transParams.get("time step increment", 1.0e-4);
    Teuchos::ParameterList MeanSolveList;
    MeanSolveList = (paramList->sublist("Model Parameters"));
#endif
    bool timeProf = transParams.get("Time Profiling", false);

    /////////////////////////////////////////////////////////////////////
    // courant-number for CFL check
    // double CN = transParams.get("Courant Number",0.9);

    // if the flow stays under cflmin for maxslow steps we increase the
    // time-step
    // int maxslow = transParams.get("Slow Flow Delay",10);
    // double slw = transParams.get("Slow Flow Condition",0.4);
    // double cflmin = slw*std::abs(CN);
#if need_locaInterface == 1
    // these LOCA data structures aren't really used by the OceanModel,
    // but are required for the ModelEvaluatorInterface

    // Create Epetra factory
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
      Teuchos::rcp(new LOCA::Epetra::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, epetraFactory);

    // the model serves as Preconditioner factory if "User Defined" is selected
    std::string PrecType = lsParams.get("Preconditioner", "Ifpack");

    Teuchos::RCP<Teuchos::ParameterList> myPrecList = Teuchos::null;
    if (PrecType == "User Defined") {
      myPrecList = Teuchos::rcp(&lsParams, false);
    }

    // for some purposes it's good to know which one is the continuation
    // parameter (i.e. backup in regular intervals)
    modelList.set("Parameter Name", "Time");

    // this is the LOCA interface (LOCA::Epetra::Interface::TimeDependent) to
    // our EpetraExt ModelEvaluator class 'OceanModel'.
    Teuchos::RCP<FVM::LocaInterface> model = Teuchos::rcp(
      new FVM::LocaInterface(modelList, Comm, globalData, myPrecList));
    Teuchos::RCP<FVM::Domain> domain = model->getDomain();
    Teuchos::RCP<LOCA::ParameterVector> pVector = model->getParameterVector();
    pVector->setValue("Time", t_start);
#else
    Teuchos::RCP<Mean> model =
      Teuchos::rcp(new Mean(rcpFromRef(MeanSolveList), &t, &dt, Comm));
#endif

    // Get the vector from the problem
    Teuchos::RCP<Epetra_Vector> soln = model->getSolution();
    Epetra_Vector soln_old(soln->Map(), Epetra_DataAccess::Copy);

#if need_locaInterface == 1
    // check for starting solution
    string StartConfigFile = modelList.get("Starting Solution File", "None");
    if (StartConfigFile != "None") {
      // HYMLS::Tools::Error("not implemented",__FILE__,__LINE__);
      INFO("Read Start File...");
      soln = model->ReadConfiguration(StartConfigFile, *pVector);
      try {
        t_start = pVector->getValue("Time");
      } catch (...) {
        Error("no time info found in start file", __FILE__, __LINE__);
      }
      transParams.set("Start Time", t_start);
      model->setParameters(*pVector);
    }
    //  INFO("Construct the Implicit Time Stepper...");
    Teuchos::RCP<AbstractTimeStepper> stepper =
      Teuchos::rcp(new ImplicitTimeStepper(model, transParams, nlParams));
#if 1
    Epetra_Vector F(*soln);
    model->computeF(*soln, F, FVM::LocaInterface::Residual);
    double nrm;
    soln->Norm2(&nrm);
    if (MyPID == 0)
      std::cout << "||u|| of starting solution: " << nrm << std::endl;
    F.Norm2(&nrm);
    if (MyPID == 0)
      std::cout << "||f|| of starting solution: " << nrm << std::endl;
#endif

#endif

    int glen = soln->GlobalLength();

  Teuchos::RCP<std::ostream> outstream2;
  if (MyPID == 0) {
    if (fabs(t_start-0.0)<.000001)
      outstream2 = Teuchos::rcp(new std::ofstream("time_dt_norm_Evyvy.txt"));
    else{
      cout<<"\n***********************************************************************\n";
      cout<<"NOTE : Will APPEND \"time_dt_norm_Evyvy.txt\" file if already present !!! \n";
      cout<<"***********************************************************************\n";
      outstream2 = Teuchos::rcp(new std::ofstream("time_dt_norm_Evyvy.txt", ios::app));
    }
  } else {
    outstream2 = Teuchos::rcp(new Teuchos::oblackholestream());
  }

#if quasi_geo == 1
  // NOTE: following settings are for very specific QG problem and when no
  // initial solution is provided! 
  double Re=MeanSolveList.get("Reynolds Number", 0.0);
  int stableSol=MeanSolveList.get("Stable Solution", 0);
  if(Re==40.0 && stableSol==1)
  {
    std::cout<<"\n setting the initial solution for Re = "<<Re;
    std::cout<<" and looking for stable solution.\n";
    for (int i=0; i<soln->MyLength();i++)
      (*soln)[i]= exp(cos(M_PI*double(i)/double(glen)));
  }
  else if(Re==40.0 && stableSol==0)
  {
    std::cout<<"\n setting the initial solution for Re = "<<Re;
    std::cout<<" and looking for unstable solution.\n";
    soln->PutScalar(0.0);
  }
  else
    soln->PutScalar(0.0);
#endif
// when initial solution is provided.
#if need_locaInterface==0
  std::string initSolnFile=MeanSolveList.get("Initial Solution File", "None");
  if(initSolnFile.compare("None")!=0)
  {
    Epetra_MultiVector *initsol;
    EpetraExt::MatrixMarketFileToMultiVector(initSolnFile.c_str(), soln->Map(), initsol);
    *soln=*(*initsol)(0);

  }
#endif



    bool increase_dt = false;

    /////////////////////////////////////////////////////////////////////////////
    Teuchos::Ptr<Teuchos::ParameterList> stochParamListptr =
      Teuchos::ptr(new Teuchos::ParameterList);
    if (MyPID == 0)
      cout << "file name to be read is StochasticParams.xml" << std::endl;
    // override default settings with parameters from user input file
    Teuchos::updateParametersFromXmlFile("StochasticParams.xml",
                                         stochParamListptr);

    Teuchos::RCP<Teuchos::ParameterList> stochParamList =
      rcpFromPtr(stochParamListptr);
    // extract the final sublists:

    std::string ProbName = stochParamList->get("Problem Type", "MOC");
    // Get the model sublist
    Teuchos::ParameterList& stochModelList =
      stochParamList->sublist("StochModel");
    double prntintvl =
      stochModelList.get("Save Solution Interval", 1.0) / timesc;

    // std::string MassFile =     stochModelList.get("MassMtxFile","mass.mm");
    std::string FrcFile = stochModelList.get("StochFrcFile", "frc.mm");
    double StochFrcStren = stochModelList.get("StochFrc Strength", 1.0);
    bool Scaling = stochModelList.get("Scaling", false);
    // std::string JacFile  =     stochModelList.get("JacDetFile","jac.mm");
    // Get the "Time Stepping" sublist
    Teuchos::ParameterList& BasisParams =
      stochModelList.sublist("Stochastic Basis");

    /* No.of bases .....Should come from xml file!!!!!!!!!!*/
    int numvecV = BasisParams.get("Max. Stoch subspace Dimension", 10);

    Teuchos::ParameterList& CoefParams =
      stochModelList.sublist("Stochastic System");
    std::string CoefFile = CoefParams.get("StochCoefFile", "y.mm");
    int Stochit = CoefParams.get("Stochastic Iterations", 100);
    int numSubTimeStep = CoefParams.get("Number of Sub Time Steps", 1);
    int maxnumiter = CoefParams.get("Max Num of Iter", 10);
    bool usebacktrack = CoefParams.get("Use BackTracking", true);
    double backtrackstep = CoefParams.get("Max Num of BackTracking Steps", 10);
    double tolRHS = CoefParams.get("Tolerance of RHS", 10e-6);
    double normRHS = CoefParams.get("Norm of RHS", 1.0);

    // Initialize class for solving basis
    Epetra_Vector V1(*soln);
    Epetra_MultiVector* frc;
    // No.of bases .....Should come from xml file!!
    int numvecW = model->get_dim_W();
    double pi = 3.14159265359;
    // Hard Coded forcing on h equation!!
    // Only work for 10*20*2*3 grid!!
    /*if(FrcFile != "None")
      {
      EPETRA_CHK_ERR(EpetraExt::MatrixMarketFileToMultiVector(FrcFile.c_str(),soln->Map(),frc));

      numvecW=frc->NumVectors();
      cout<< "frc read \n";
      }
      else
      {
      if(MyPID==0)
      {
      cout<<"\nWARNING!!!-->> No File for Stochastic forcing is provided."
      <<"Assumed to be hard coded with one column(after this message).\n";
      }
      */

    frc = new Epetra_MultiVector(V1.Map(), numvecW);
    frc->PutScalar(0.0);
    /*******************mass matrix is negative************** */
    Teuchos::RCP<Epetra_CrsMatrix> massmat;
#if need_locaInterface == 1
    massmat = Teuchos::rcp(new Epetra_CrsMatrix(*(model->getMassMatrix())));
    *massmat = *(model->getMassMatrix());
    massmat->Scale(-1.0);
#else
    massmat = (model->getMassMatrix());
#endif
    if (!massmat->Filled())
      massmat->FillComplete();

    const char* massfile = "mass.mm";
    EpetraExt::RowMatrixToMatrixMarketFile(massfile, *massmat);
    Epetra_Vector Fu(*soln);
#if need_locaInterface == 1
    model->computeF(*soln, Fu, FVM::LocaInterface::Residual);
#else
    model->computeF(*soln, Fu);
#endif
    Teuchos::RCP<Epetra_CrsMatrix> A = model->getJacobian();
    if (!A->Filled()) {
      std::cout << "\n A is not filled\n";
      A->FillComplete();
    }

    //    soln->PutScalar(0.0);
    //    for (int i = 2; i < soln->MyLength(); i += 3)
    //      (*soln)[i] = 0.1 * sin(i * pi / glen);
    // const char* initSolFile = "initSol.mm";
    // EpetraExt::MultiVectorToMatrixMarketFile(initSolFile, *soln);
    Epetra_CrsMatrix detA(*A);
    model->computeJacobian(
      *soln, detA); // FIXME :detA will be overwritten, so what is this???
    detA.FillComplete();
    const char* jac = "jac.mm";
    EpetraExt::RowMatrixToMatrixMarketFile(jac, detA);
    EpetraExt::MatrixMatrix::Add(*A, false, 1.0, detA, 0.0); // detA=A
    Epetra_Vector ScR(*soln);
    Epetra_Vector ScL(*soln);
    // Epetra_CrsMatrix MassMat(Copy,massmat->RowMap(),1,false);
    // Epetra_CrsMatrix Jac(Copy,massmat->RowMap(),5,false);
    if (Scaling && ProbName == "SWE") {
      Epetra_MultiVector* Frc;
      double val;
      for (int i = 0; i < massmat->NumMyRows(); i = i + 3) {
        ScR[i] = 1;
        ScR[i + 1] = 1;
        ScR[i + 2] = 1e-4 / 8;
        ScL[i] = 1;
        ScL[i + 1] = 1;
        ScL[i + 2] = 1e-2;
      }
      Frc = new Epetra_MultiVector(*frc);
      massmat->RightScale(ScR);
      massmat->LeftScale(ScL);
      detA.RightScale(ScR);
      detA.LeftScale(ScL);
      Frc->Multiply(1.0, ScL, *frc, 1.0);
      *frc = *Frc;
      delete Frc;
    } else if (Scaling && ProbName != "SWE" && MyPID == 0) {
      cout << "\nWARNING : Scaling is ON but it has NOT been implemented"
           << "\n for the current Problem Type. Apply it manually!!!!\n";
    }
    Teuchos::RCP<Epetra_MultiVector> Wbase;
    Wbase = Teuchos::rcp(new Epetra_MultiVector(soln->Map(), numvecW));
#if need_locaInterface == 1
    Teuchos::RCP<Epetra_Map> fortMap = domain->GetAssemblyMap();
    Epetra_MultiVector fortStochFrc(*fortMap, 1);
    double* W_ptr;
    int ierr;
    std::cout << "\nnumber of rows in W_ = " << fortStochFrc.GlobalLength()
              << "\n";
    int n = fortStochFrc.MyLength();
    for (int col = 0; col < fortStochFrc.NumVectors(); col++) {
      fortStochFrc(col)->ExtractView(&W_ptr);
      model_stoch_frc(&n, W_ptr, &col, &ierr);
    }
    CHECK_ZERO(domain->Assembly2Solve(fortStochFrc, *Wbase));
    HYMLS::MatrixUtils::Dump(*Wbase, "Forcing.txt");
    //*Wbase = *frc;
    // bool success3 = model->computeShiftedMatrix(1.0, 0.0, *soln, detA);
#else
    Wbase = rcp(new Epetra_MultiVector(soln->Map(), numvecW));
    Wbase = model->get_W();
    numvecW = Wbase->NumVectors();
#endif
    delete frc;
    Wbase->Scale(StochFrcStren);
    std::flush(std::cout << "\n frob. norm of detA = " << detA.NormFrobenius()
                         << std::endl);
    Teuchos::RCP<Problem_Interface> Vstoch =
      Teuchos::rcp(new Problem_Interface(Teuchos::rcpFromRef(detA),
                                         Teuchos::rcpFromRef(BasisParams),
                                         model->getSolution(),
                                         numvecV,
                                         &t,
                                         dt,
                                         Wbase,
                                         massmat,
                                         Stochit));

    // std::string frcfile = "frc.mm";
    //    EpetraExt::MultiVectorToMatrixMarketFile(frcfile.c_str(), Wbase);
    Vstoch->set_frcStrength(StochFrcStren);
    Vstoch->init_v(t);

    /*****************************************************************************************************/
    // initialization for Y_stoch class
    /*****************************************************************************************************/
    int NumGlobalElements =
      numvecV; // should be equal to the no of bases given to Vstoch class

    Teuchos::RCP<Epetra_MultiVector> Vn = Vstoch->V;

    if (MyPID == 0)
      cout << "\nInitializing Stoc. Coeff. class..... " << std::endl;
#if need_locaInterface == 1
    Teuchos::RCP<Y_Stoch> y_interface =
      Teuchos::rcp(new Y_Stoch(Stochit,
                               numSubTimeStep,
                               NumGlobalElements,
                               &dt,
                               Teuchos::rcpFromRef(detA),
                               soln,
                               domain,
                               Vn,
                               Wbase,
                               Comm,
                               Teuchos::rcpFromRef(CoefParams),
                               maxnumiter,
                               usebacktrack,
                               backtrackstep,
                               tolRHS,
                               normRHS));
#else
    Teuchos::RCP<Y_Stoch> y_interface =
      Teuchos::rcp(new Y_Stoch(Stochit,
                               numSubTimeStep,
                               NumGlobalElements,
                               &dt,
                               Teuchos::rcpFromRef(detA),
                               soln,
                               Vn,
                               Wbase,
                               Comm,
                               Teuchos::rcpFromRef(CoefParams),
                               model,
                               maxnumiter,
                               usebacktrack,
                               backtrackstep,
                               tolRHS,
                               normRHS));
#endif
    if (MyPID == 0)
      cout << "\nDONE" << std::endl;
    Vstoch->ExpDExpyy = y_interface->ExpDExpyy;
    Vstoch->y = y_interface->y_;
    double nrm1[Vstoch->V->NumVectors()];
#if need_locaInterface == 1
    model->printSolution(*soln, t);
#else
    string fname = "MeanSol_" + std::to_string(t) + ".text";
    model->WriteSolution(fname, t, *soln);
#endif
    y_interface->HBilinV();
    y_interface->computeEyyTyT();
    model->setExpVyVy(y_interface->getEVyVy());
    double count = t;
    int yycount = 0;
    if (MyPID == 0) {
      cout << "\n end time = " << t_end;
      cout << "\n start time = " << t;
      cout << "\n starting time iterations:"
           << "\n";
    }
    printnormMV(*soln, 2, "inital norm of soln");
#if quasi_geo == 1
    model->setSolution(*soln);
#endif
    Teuchos::RCP<Epetra_MultiVector> expyy = y_interface->getEyy();
    double saveEyyintrvl=0.2;
    int nv = round(saveEyyintrvl/dt);
    std::cout<< "\n nv = "<<nv<<"\n";
    int ttlelmnt=expyy->NumVectors()*expyy->GlobalLength();
    Epetra_MultiVector Eyy(Epetra_Map(ttlelmnt,0,*Comm),nv);
    INFO("Start Time integration");
    /*****************************************************************************************************/
    // START time integration loop
    /*****************************************************************************************************/
    while (t < t_end) {
#if need_locaInterface == 1
      if (step > maxsteps) {
        (*outstream) << "Maximum number of steps exceeded, stopping...\n";
        cout << "Maximum number of steps exceeded, stopping...\n";
        break;
      }
      (*outstream) << std::endl;
      (*outstream) << "##############################################\n";
      (*outstream) << "start time step " << step << ": t=" << time_out(t)
                   << std::endl;
#endif
      // make sure that we hit t_end exactly
      // If the time-step would become too
      // small, we allow the stepper to
      // go a little too far (may want to fix
      // this)
      soln_old = *soln;
      
#if need_locaInterface == 1
      if (t + dt > t_end) {
        dt = std::max(dt_min, t_end - t);
      }

      (*outstream) << "step size dt = " << time_out(dt) << std::endl;
      (*outstream) << "##############################################\n";
      (*outstream) << std::endl;

      // bool success2 = model->computeShiftedMatrix(1.0, 0.0, soln_old, detA);
#endif
      /***************************************************************
        Solve For stochastic coefficients
       *****************************************************************/
      // if (MyPID == 0)
      // std::flush(std::cout << "\n solving for Y");
      if (timeProf)
        coeffTime->start();

      y_interface->StochasticIterations();

      if (timeProf)
        coeffTime->stop();
      coeffTime->incrementNumCalls();
      /************************************************************
        Solve For Stochastic Bases
       *************************************************************/
      // Crude iteration scheme for V-equation
      // copy Vn in a temporary variable
      // if (MyPID == 0)
      // std::flush(std::cout << "\n solving for V");
      // Note:Temp Disable

      Epetra_MultiVector Vtemp(*Vn);
      //// Do a loop
      if (timeProf)
        basisTime->start();

      for (int i = 0; i < 1; i = i + 1) {
        Vstoch->computeBlocks(dt);
        Vstoch->v_stoch_init(&Vtemp);
      }
      Vstoch->TransferNorm();


      if (timeProf)
      {
        basisTime->stop();
	basisTime->incrementNumCalls();
      }

      if (timeProf)
        coeffTime2->start();

      y_interface->HBilinV();
      y_interface->computeEyyTyT();
      y_interface->computeEVyVy();
      

      if (timeProf)
      {
        coeffTime2->stop();
	coeffTime2->incrementNumCalls();
      }
        /******************************************************
          SOLVE FOR Udet!!!!!!
         *****************************************************/
#if need_locaInterface == 1
      // if (MyPID == 0)
      // std::flush(std::cout << "\n solving for mean");
      if (timeProf)
        meanTime->start();

      bool success = stepper->Step(soln_old, t, *soln, dt);
      
      if (timeProf)
      {
        meanTime->stop();
	meanTime->incrementNumCalls();
      }
        // std::flush(cout << "\nsuccess = " << success);
#else
      if (timeProf)
        meanTime->start();
      bool success = model->NewtonSolver();
      //bool success = model->newtonLineSearchSolve(*soln);
      if (timeProf)
      {
      	meanTime->stop();
	meanTime->incrementNumCalls();
      }
#endif
     EpetraExt::MatrixMatrix::Add(*A, false, 1.0, detA, 0.0);
      if (Scaling && ProbName == "SWE") {
        detA.RightScale(ScR);
        detA.LeftScale(ScL);
      }
      increase_dt = success;
#if need_locaInterface == 1
      if (!success) {
        if (dt == dt_min) {
          (*outstream)
            << "Time stepping failed at minimum step-size, stopping run...\n";
          break;
        } else {
          (*outstream) << "Step failed!\n";
        }
      }
#else
      if (!success) {
        (*outstream) << "Time stepping failed, stopping run...\n";
      }
#endif
      (*outstream) << std::endl;
      (*outstream) << "##############################################\n";
#if need_locaInterface == 1
      double lte = -1.0;

      double effortEst = stepper->getEffortEstimate();
      (*outstream) << "Step effort: " << effortEst * 100 << "%\n";

      // success only indicates wether the nonlinear solver converged,
      // if we have more criteria we must check them now:
      if (success && (stepper->hasErrorEstimate())) {
        lte = stepper->getErrorEstimate();
        (*outstream) << "Local truncation error tau: " << lte << std::endl;

        if (success && stepChoice == "Adaptive") {
          if (lte > lteTol) {
            if (dt > dt_min) {
              success = false;
              (*outstream) << "Step rejected: estimate of LTE (" << lte
                           << ")\n";
              (*outstream) << "               does not satisfy tolerance ("
                           << lteTol << ")\n";
            } else {
              (*outstream)
                << "WARNING: step failed to achieve desired accuracy, \n";
              (*outstream)
                << "         but as the step-size is already at its   \n";
              (*outstream)
                << "         minimum, we accept it anyway.            \n";
            }
          }
        }
      }

      // advance one step:
      if (success) {
        if (MyPID == 0) {
          std::flush(cout << "+");
//          std::flush(cout << yycount<<" ");
        }
        t += dt;
        pVector->setValue("Time", t);
        model->setParameters(*pVector);
        std::string Filename;
	double cpyeyy[ttlelmnt];
	expyy->ExtractCopy(cpyeyy,expyy->MyLength());
	for(int i=0; i<expyy->NumVectors();i++)
	{
	  for(int j=0; j<expyy->MyLength();j++)
	  {
	    Eyy.ReplaceGlobalValue(i*expyy->MyLength()+j,yycount,(*expyy)[i][j]);
	  }
	}
	yycount=yycount+1;
	if (yycount==nv ) {
	  Filename = "tsEyy_" + Teuchos::toString(float(t)) + ".mm";
	  HYMLS::MatrixUtils::mmwrite(Filename, Eyy);
	  yycount=0;
	}
        if ((t - count - prntintvl) >= 0) {
	  Filename = "v_" + Teuchos::toString(float(t)) + ".mm";
	  HYMLS::MatrixUtils::mmwrite(Filename, *Vn);
          // EpetraExt::MultiVectorToMatrixMarketFile(Filename.c_str(), *Vn);
          Filename = "yT_" + Teuchos::toString(float(t)) + ".mm";
	  HYMLS::MatrixUtils::mmwrite(Filename, *(y_interface->yTrans_));
          // EpetraExt::MultiVectorToMatrixMarketFile(Filename.c_str(),
          //                                         *(y_interface->yTrans_));
	  Filename = "mean_" + Teuchos::toString(float(t)) + ".mm";
	  HYMLS::MatrixUtils::mmwrite(Filename, *soln);
	  //model->WriteConfiguration(Filename, *pVector, *soln);
          //EpetraExt::MultiVectorToMatrixMarketFile(Filename.c_str(), *soln);
          count = count + prntintvl; // getchar();
        }
        step++;
        (*outstream) << "Step was successful!\n";
        const Epetra_Vector& fn = stepper->getF();
        double nrmF, nrmexp;
        CHECK_ZERO(fn.Norm2(&nrmF));
        (y_interface->getEVyVy())->Norm2(&nrmexp);
        (*outstream2) << t << "\t" << dt << "\t" << nrmF << "\t" << nrmexp
                      << std::endl;
      }

      if (effortEst < maxEffort) {
        increase_dt = success;
        nslow = 0;
      } else {
        (*outstream)
          << "step size not increased because system to hard to solve\n";
        increase_dt = false;
      }

      if (stepChoice == "Constant") // step-size control only based on nonlinear
                                    // solver convergence
      {
        if (increase_dt) {
          if (dt < dt_max) {
            (*outstream) << "increasing step size for next step\n";
          }
          dt = min(dt_max, inc_fac * dt);
        } else if (!success) // reduce step size and try again
        {
          (*outstream) << "reducing step size\n";
          dt = max(dt_min, red_fac * dt);
          nslow = 0;
        }
      } else if (stepChoice ==
                 "Adaptive") // adaptivity based on error estimate,
      {
        if (lte != -1.0) {
          if ((lte > lteTol) || (lte < lteTol / 4)) {
            // compute new step size
            double newdt = sqrt(lteTol / (2 * lte)) * dt;
            newdt = max(min(newdt, 2 * dt), dt / 2);
            newdt = min(newdt, dt_max);
            newdt = max(newdt, dt_min);
            DEBVAR(newdt);
            dt = newdt;
          }
        } else // same as "Constant". hasErrorEstimate==false may
        {      // simply mean that this was the first step or a restart
          if (!success) {
            (*outstream) << "reducing step size\n";
            dt = max(dt_min, dt / 2);
          }
        }
      }

      // check for convergence to a steady state (when dU/dt->0)
      if (success && steadyMon == "Norm of F") {
        const Epetra_Vector& fn = stepper->getF();
        double nrmF;
        CHECK_ZERO(fn.Norm2(&nrmF));
        (*outstream) << "Norm of F(u): " << nrmF << std::endl;
        if (nrmF < steadyTol) {
          (*outstream)
            << "Convergence to a steady-state has been detected, stopping here!"
            << std::endl;
          break;
        }
      }
#else
      if (success) {
        if (MyPID == 0) {
          std::flush(cout << "+");
        }
        t += dt;
        dt = min(dtIncrmnt + dt, t_end - t);

        std::string Filename;
	double cpyeyy[ttlelmnt];
	expyy->ExtractCopy(cpyeyy,expyy->MyLength());
	for(int i=0; i<expyy->NumVectors();i++)
	{
	  for(int j=0; j<expyy->MyLength();j++)
	  {
	    Eyy.ReplaceGlobalValue(i*expyy->MyLength()+j,yycount,(*expyy)[i][j]);
	  }
	}
	yycount=yycount+1;
	if (yycount==nv ) {
	  Filename = "tsEyy_" + Teuchos::toString(float(t)) + ".mm";
	  HYMLS::MatrixUtils::mmwrite(Filename, Eyy);
	  yycount=0;
	}
        if ((t - count - prntintvl) >= 0) {
	  
          Filename = "v_" + Teuchos::toString(float(t)) + ".mm";
          EpetraExt::MultiVectorToMatrixMarketFile(Filename.c_str(), *Vn);
          Filename = "yT_" + Teuchos::toString(float(t)) + ".mm";
          EpetraExt::MultiVectorToMatrixMarketFile(Filename.c_str(),
                                                   *(y_interface->yTrans_));
          Filename = "mean_" + Teuchos::toString(float(t)) + ".mm";
          EpetraExt::MultiVectorToMatrixMarketFile(Filename.c_str(), *soln);
          count = count + prntintvl; // getchar();
        }
        const Epetra_Vector& fn = model->getF();
        double nrmF;
        CHECK_ZERO(fn.Norm2(&nrmF));
        (*outstream) << "Step was successful!\n";
        (*outstream2) << t << "\t" << dt << "\t" << nrmF << std::endl;
      }
#endif

    } // while loop
    if (Scaling && ProbName == "SWE") {
      Epetra_MultiVector Vtemp(*Vn);
      Vtemp.Multiply(1.0, ScL, *Vn, 1.0);
      *Vn = Vtemp;

      ScL.Reciprocal(ScL);
      ScR.Reciprocal(ScR);
      massmat->RightScale(ScR);
      massmat->LeftScale(ScL);
      Vstoch->MOrth();
      Vstoch->TransferNorm();
    }

    Epetra_Vector scnd_mmnt(*soln);
    y_interface->PostProcess(scnd_mmnt);
    double time = timer.ElapsedTime();
    (*info) << "\n*******************************************\n";
    (*info) << "Elapsed time during entire run: " << (int)(time / 60)
            << " min\n";
    (*info) << "*******************************************\n\n";

    (*outstream) << std::endl;
    (*outstream)
      << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    (*outstream)
      << "Time-stepping run finished, store final solution and finish...\n";
    (*outstream) << "Elapsed time during entire run: " << (int)(time / 60)
                 << " min\n";
    (*outstream)
      << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
    (*outstream) << std::endl;

    // TODO: the openDX file is finished automatically by the OceanModel
    // destructor,
    //       but the final state is generally not stored!

#if need_locaInterface == 1
    model->WriteConfiguration("MeanSolnFinal.txt", *pVector, *soln);
#else

    HYMLS::MatrixUtils::Dump(*soln, "MeanSolnFinal.mm");
#endif
    INFO(" Time-stepping run finished, store solution...");
    INFO("done!");
    std::string flnm = "V_BASE.mm";
    HYMLS::MatrixUtils::mmwrite(flnm, *Vn);
    flnm = "YTrans_COEFF.mm";
    HYMLS::MatrixUtils::mmwrite(flnm, *(y_interface->yTrans_));

    (*outstream) << "Final Parameters: " << std::endl;
    (*outstream) << *paramList << std::endl;
    if (timeProf) {
      std::ofstream timeFile("timeProf.yml", ios::out);
      Teuchos::RCP<Teuchos::ParameterList> reportTimeParams =
        parameterList(*(Teuchos::TimeMonitor::getValidReportParameters()));
      // pars->set("Report format", "YAML");
      // pars->set("YAML style", "compact");
      // Get a summary from the time monitor.
      Teuchos::TimeMonitor::report(timeFile, "my", reportTimeParams);
    }
    // std::cout << timeFile.str() << std::endl;
    // The HDF5 file is closed when OceanModel is deleted => OceanOutput is
    // deleted, but this happens only at the end of main because of the stepper
    // hook. we need to make sure everyone closes the file before calling
    // MPI_Finalize()
    // TODO: This is a bug, we probably violated one or more of the 10
    // Teuchos::RCP-commandments
    // model->finishOutput();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, status);
  //}catch(...){std::cout << "caught something!\n";}
  ///////////////////////////////////////////////////////////////

  // end main

  Comm->Barrier();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return status ? EXIT_SUCCESS : EXIT_FAILURE;
}

/*
 * =====================================================================================
 *
 *       Filename:  timedependent_swe.cpp
 *
 *    Description:  DO time-integration driver for the Shallow Water Equations
 *                  (SWEDO target).  Uses NOX/LOCA for implicit time stepping
 *                  with LTE-based adaptive step-size control.  Compiled
 *                  exclusively with need_locaInterface=1; contains no
 *                  brgr / quasi_geo conditionals.
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

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using std::ofstream;
using namespace std;

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "AztecOO.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_OperatorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#include "StochIO.hpp"
#include "EpetraExt_VectorIn.h"
#include "HYMLS_HyperCube.hpp"
#include "HYMLS_MatrixUtils.hpp"
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
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "FVM_model_interface.h"
#include "globdefs.H"
#include "Interface.hpp"        // provides FVM::LocaInterface + ImplicitTimeStepper
#include "Problem_Interface.hpp"
#include "StochSys.hpp"
#include "DOUtils.hpp"
#include "DOTimeLoop.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#ifndef INFO
#define INFO(s) std::cout << s << std::endl;
#endif

static const double r0dim = 6370.0e3;
static const double udim = 0.1;
static const double timesc = r0dim / udim; // [s]

// Non-static: libSWE.so calls this symbol directly.
const double comp_time_scaling() { return timesc; }

static string time_out(double t) {
    std::stringstream s;
    s << t * timesc << " [s]";
    return s.str();
}

// Global timers: referenced as extern symbols inside libSWE.so
// (ImplicitTimeStepper, ThetaStepperEvaluator).
Teuchos::RCP<Teuchos::Time> meanSolverTime = Teuchos::TimeMonitor::getNewTimer("my time spent on mean solver");
Teuchos::RCP<Teuchos::Time> meanMassMatTime = Teuchos::TimeMonitor::getNewTimer(
    "my time spent on mass matrix computation");
Teuchos::RCP<Teuchos::Time> meanRhsTime = Teuchos::TimeMonitor::getNewTimer("my total time spent on RHS computation");
Teuchos::RCP<Teuchos::Time> meanJacTime = Teuchos::TimeMonitor::getNewTimer(
    "my total time spent on Jacobian computation");

//====================================================================================

#ifdef DEBUGGING
static void test_jac_billin(Teuchos::RCP<FVM::LocaInterface> model,
                            Teuchos::RCP<Epetra_CrsMatrix> A,
                            Teuchos::RCP<Y_Stoch> Ystoch) {
    std::cout << "\n inside test_jac_bilin()\n";
    int n = A->NumMyRows();
    std::cout << "number of elements, n= " << n << "\n";
    Teuchos::RCP<Epetra_Map> xmap = Teuchos::rcp(new Epetra_Map(n, 0, A->Comm()));
    Teuchos::RCP<Epetra_Vector> vx, xv, diff, x, v, xdet;
    vx = Teuchos::rcp(new Epetra_Vector(*xmap));
    xv = Teuchos::rcp(new Epetra_Vector(*xmap));
    diff = Teuchos::rcp(new Epetra_Vector(*xmap));
    x = Teuchos::rcp(new Epetra_Vector(*xmap));
    xdet = Teuchos::rcp(new Epetra_Vector(*xmap));
    v = Teuchos::rcp(new Epetra_Vector(*xmap));
    double nrm;
    Epetra_CrsMatrix *linJac;
    EpetraExt::MatrixMarketFileToCrsMatrix("linjac.mm", A->Comm(), linJac, 0, 1);
    linJac->Scale(-1.0);
    xdet->PutScalar(10.0);
    x->PutScalar(10.0);
    v->PutScalar(0.1);
    Epetra_CrsMatrix testJac(*A);
    testJac.FillComplete();
    Epetra_Vector J0v(*x);
    linJac->Multiply(false, *v, J0v);
    model->computeJacobian(*x, testJac);
    Epetra_Vector Jxv(*x);
    testJac.Multiply(false, *v, Jxv);
    Ystoch->Bilinear(xdet, v, x, vx);
    Ystoch->Bilinear(xdet, x, v, xv);
    diff->Update(1.0, *xv, 1.0, *vx, 0.0);
    diff->Update(1.0, J0v, -1.0, Jxv, 1.0);
    diff->MaxValue(&nrm);
    StochIO::writeMatrix("jacobian.mm", testJac);
    StochIO::writeMV("diff.mm", *diff);
    diff->Norm2(&nrm);
    std::cout << "\n Norm of J(0)*v - Bilin(x,v) - Bilin(v,x) - J(x)*v = " << nrm << "\n";
    getchar();
}
#endif

//====================================================================================

int main(int argc = 0, char *argv[] = NULL) {
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif
    bool status = true;

    // ===== PHASE 1: Communicator and output streams =====

#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Teuchos::RCP<Epetra_SerialComm> Comm = rcp(new Epetra_SerialComm());
#endif
    Epetra_Time timer(*Comm);
    int MyPID = Comm->MyPID();
    int NumProc = Comm->NumProc();

    if (MyPID == 0) {
        std::ostringstream fname;
        fname << "info_" << MyPID << ".txt";
        info = Teuchos::rcp(new std::ofstream(fname.str().c_str()));
    } else {
        info = Teuchos::rcp(new Teuchos::oblackholestream());
    }
    (*info) << std::setw(15) << std::setprecision(15);

    Teuchos::RCP<std::ostream> outstream;
    if (MyPID == 0)
        outstream = Teuchos::rcp(new std::ofstream("timestepping.out"));
    else
        outstream = Teuchos::rcp(new Teuchos::oblackholestream());
    (*outstream) << std::setw(15) << std::setprecision(15);

    DEBUG("*********************************************")
    DEBUG("* Debugging output for process " << MyPID)
    DEBUG("*********************************************")

    INFO(*Comm);

    Teuchos::RCP<Teuchos::Time> coeffTime = Teuchos::TimeMonitor::getNewTimer("Stoch. Coeff. time");
    Teuchos::RCP<Teuchos::Time> coeffTime2 = Teuchos::TimeMonitor::getNewTimer("Bilinear form time");
    Teuchos::RCP<Teuchos::Time> basisTime = Teuchos::TimeMonitor::getNewTimer("Stoch. Basis time");
    Teuchos::RCP<Teuchos::Time> meanTime = Teuchos::TimeMonitor::getNewTimer("Mean solve time");
    StochTimers stochTimers = {coeffTime, coeffTime2, basisTime};

    try {
        // ===== PHASE 2: Read parameter lists =====

        Teuchos::Ptr<Teuchos::ParameterList> paramListptr = Teuchos::ptr(new Teuchos::ParameterList);
        Teuchos::updateParametersFromXmlFile("params.xml", paramListptr);
        Teuchos::RCP<Teuchos::ParameterList> paramList = rcpFromPtr(paramListptr);

        Teuchos::ParameterList &modelList = paramList->sublist("Model");
        Teuchos::ParameterList &nlParams = paramList->sublist("NOX");
        Teuchos::ParameterList &printParams = nlParams.sublist("Printing");
        printParams.set("MyPID", MyPID);
#ifdef OUTPUT_TO_FILE
        printParams.set("Output Stream", outstream);
        printParams.set("Error Stream", outstream);
#endif
        printParams.set("Output Process", 0);
        printParams.set("Output Information", 0
#ifdef DEBUGGING
                        +NOX::Utils::Details + NOX::Utils::OuterIteration +
                                NOX::Utils::InnerIteration + NOX::Utils::OuterIterationStatusTest +
                                NOX::Utils::LinearSolverDetails + NOX::Utils::Debug +
                                NOX::Utils::StepperParameters + NOX::Utils::Warning +
                                NOX::Utils::StepperDetails + NOX::Utils::StepperIteration
#endif
        );
        Teuchos::ParameterList &dirParams = nlParams.sublist("Direction");
        Teuchos::ParameterList &searchParams = nlParams.sublist("Line Search");
        Teuchos::ParameterList &newtParams = dirParams.sublist("Newton");
        Teuchos::ParameterList &lsParams = newtParams.sublist("Linear Solver");

        Teuchos::ParameterList &transParams = paramList->sublist("Time Stepping");
        double t_start = transParams.get("Start Time", 0.0) / timesc;
        double t_end = transParams.get("End Time", 1.0) / timesc;
        double dt = transParams.get("Step Size", 0.01) / timesc;
        double t = t_start;
        bool timeProf = transParams.get("Time Profiling", false);
        string stepperScheme = transParams.get("Scheme", "Theta");

        bool gradual_startup = transParams.get("Gradual Start-Up", true);
        string startup_param = "None";
        double startup_rate = 0.0;
        double startup_max = 0.0;
        if (gradual_startup) {
            startup_param = transParams.get("Start-Up Parameter", "Wind Forcing");
            startup_rate = transParams.get("Start-Up Rate", 0.1);
            startup_max = modelList.sublist("Starting Parameters").get(startup_param, 0.0);
        }

        transParams.set("Output Stream", outstream);
        double dt_min = transParams.get("Minimum Step Size", dt / 10.0) / timesc;
        double dt_max = transParams.get("Maximum Step Size", dt) / timesc;
        int maxsteps = transParams.get("Max Num Steps", (int) ((t_end - t_start) / dt));
        double lteTol = transParams.get("Error Tolerance", 1.0e-3);
        string steadyMon = transParams.get("Steady-State Monitor", "Norm of F");
        double steadyTol = transParams.get("Steady-State Tolerance", 0.0);
        string stepChoice = transParams.get("Step Size Control", "Constant");
        double red_fac = transParams.get("Failed Step Reduction Factor", 0.5);
        double inc_fac = transParams.get("Successful Step Increase Factor", 1.5);
        double maxEffort = 0.8;
        int step = 1;

        // ===== PHASE 3: Initialize FVM model and implicit time stepper =====

        Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
                Teuchos::rcp(new LOCA::Epetra::Factory);
        Teuchos::RCP<LOCA::GlobalData> globalData =
                LOCA::createGlobalData(paramList, epetraFactory);

        string PrecType = lsParams.get("Preconditioner", "Ifpack");
        Teuchos::RCP<Teuchos::ParameterList> myPrecList = Teuchos::null;
        if (PrecType == "User Defined")
            myPrecList = Teuchos::rcp(&lsParams, false);

        modelList.set("Parameter Name", "Time");
        Teuchos::RCP<FVM::LocaInterface> model =
                Teuchos::rcp(new FVM::LocaInterface(modelList, Comm, globalData, myPrecList));
        Teuchos::RCP<FVM::Domain> domain = model->getDomain();
        Teuchos::RCP<LOCA::ParameterVector> pVector = model->getParameterVector();
        pVector->setValue("Time", t_start);

        Teuchos::RCP<Epetra_Vector> soln = model->getSolution();
        Epetra_Vector soln_old(soln->Map(), Epetra_DataAccess::Copy);

        {
            string StartConfigFile = modelList.get("Starting Solution File", "None");
            if (StartConfigFile != "None") {
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
        }

        Teuchos::RCP<AbstractTimeStepper> stepper =
                Teuchos::rcp(new ImplicitTimeStepper(model, transParams, nlParams));

        {
            Epetra_Vector F(*soln);
            model->computeF(*soln, F, FVM::LocaInterface::Residual);
            double nrm;
            soln->Norm2(&nrm);
            if (MyPID == 0) cout << "||u|| of starting solution: " << nrm << "\n";
            F.Norm2(&nrm);
            if (MyPID == 0) cout << "||f|| of starting solution: " << nrm << "\n";
        }

        // ===== PHASE 4: Initialize basis and forcing =====

        Teuchos::Ptr<Teuchos::ParameterList> stochParamListptr =
                Teuchos::ptr(new Teuchos::ParameterList);
        if (MyPID == 0) cout << "Reading StochasticParams.xml\n";
        Teuchos::updateParametersFromXmlFile("StochasticParams.xml", stochParamListptr);
        Teuchos::RCP<Teuchos::ParameterList> stochParamList = rcpFromPtr(stochParamListptr);

        std::string ProbName = stochParamList->get("Problem Type", "MOC");
        Teuchos::ParameterList &stochModelList = stochParamList->sublist("StochModel");
        double prntintvl = stochModelList.get("Save Solution Interval", 1.0) / timesc;
        double saveEyyIntrvl = stochModelList.get("Save Eyy Interval", 0.4) / timesc;
        double StochFrcStren = stochModelList.get("StochFrc Strength", 1.0);
        bool Scaling = stochModelList.get("Scaling", false);
        bool isStochOn = stochModelList.get("Use Stochastic", true);

        Teuchos::ParameterList &BasisParams = stochModelList.sublist("Stochastic Basis");
        int numvecV = BasisParams.get("Max. Stoch subspace Dimension", 10);

        Teuchos::ParameterList &CoefParams = stochModelList.sublist("Stochastic System");
        int Stochit = CoefParams.get("Stochastic Iterations", 100);
        int numSubTimeStep = CoefParams.get("Number of Sub Time Steps", 1);
        int maxnumiter = CoefParams.get("Max Num of Iter", 10);
        bool usebacktrack = CoefParams.get("Use BackTracking", true);
        double backtrackstep = CoefParams.get("Max Num of BackTracking Steps", 10);
        double tolRHS = CoefParams.get("Tolerance of RHS", 10e-6);
        double normRHS = CoefParams.get("Norm of RHS", 1.0);

        // Mass matrix: negated to match FVM sign convention
        Teuchos::RCP<Epetra_CrsMatrix> massmat =
                Teuchos::rcp(new Epetra_CrsMatrix(*(model->getMassMatrix())));
        *massmat = *(model->getMassMatrix());
        massmat->Scale(-1.0);
        if (!massmat->Filled()) massmat->FillComplete();

        {
            Epetra_Vector Fu(*soln);
            model->computeF(*soln, Fu, FVM::LocaInterface::Residual);
        }
        Teuchos::RCP<Epetra_CrsMatrix> A = model->getJacobian();
        if (!A->Filled()) A->FillComplete();

        Epetra_CrsMatrix detA(*A);
        model->computeJacobian(*soln, detA);
        detA.FillComplete();
        dumpInitialJacobian(A, detA, *massmat);

        // Scaling (SWE-specific)
        Epetra_Vector ScR(*soln), ScL(*soln);
        if (Scaling) {
            for (int i = 0; i < massmat->NumMyRows(); i += 3) {
                ScR[i] = 1;
                ScR[i + 1] = 1;
                ScR[i + 2] = 1e-4 / 8.0;
                ScL[i] = 1;
                ScL[i + 1] = 1;
                ScL[i + 2] = 1e-2;
            }
            massmat->RightScale(ScR);
            massmat->LeftScale(ScL);
            detA.RightScale(ScR);
            detA.LeftScale(ScL);
        }

        // Stochastic forcing via Fortran interface
        int numvecW = model->get_dim_W();
        Teuchos::RCP<Epetra_MultiVector> Wbase =
                Teuchos::rcp(new Epetra_MultiVector(soln->Map(), numvecW));
        {
            Teuchos::RCP<Epetra_Map> fortMap = domain->GetAssemblyMap();
            Epetra_MultiVector fortStochFrc(*fortMap, 1);
            double *W_ptr;
            int ierr, n = fortStochFrc.MyLength();
            for (int col = 0; col < fortStochFrc.NumVectors(); col++) {
                fortStochFrc(col)->ExtractView(&W_ptr);
                model_stoch_frc(&n, W_ptr, &col, &ierr);
            }
            CHECK_ZERO(domain->Assembly2Solve(fortStochFrc, *Wbase));
            HYMLS::MatrixUtils::mmwrite("Forcing.txt", *Wbase);
        }
        Wbase->Scale(StochFrcStren);

        // ===== Shared stochastic state =====
        auto sharedState = Teuchos::rcp(new StochasticState());
        sharedState->udet = soln;
        sharedState->A    = Teuchos::rcpFromRef(detA);
        sharedState->W    = Wbase;

        Teuchos::RCP<Problem_Interface> Vstoch =
                Teuchos::rcp(new Problem_Interface(Teuchos::rcpFromRef(detA),
                                                   Teuchos::rcpFromRef(BasisParams),
                                                   model->getSolution(),
                                                   numvecV, &t, dt,
                                                   Wbase, massmat, Stochit,
                                                   sharedState));
        Vstoch->set_frcStrength(StochFrcStren);
        Vstoch->init_v(t);

        // ===== PHASE 5: Initialize stochastic coefficient solver =====

        Teuchos::RCP<Epetra_MultiVector> Vn = Vstoch->getBasisV();
        printnormMV(*Vn, 2, "initial norm of V");
        if (MyPID == 0) cout << "\nInitializing Stoch. Coeff. class...\n";

        Teuchos::RCP<Y_Stoch> y_interface =
                Teuchos::rcp(new Y_Stoch(Stochit, numSubTimeStep, numvecV,
                                         &dt, Teuchos::rcpFromRef(detA), soln, domain,
                                         Vn, Wbase, Comm,
                                         Teuchos::rcpFromRef(CoefParams),
                                         sharedState,
                                         maxnumiter, usebacktrack, backtrackstep,
                                         tolRHS, normRHS));

#ifdef DEBUGGING
        test_jac_billin(model, Teuchos::rcpFromRef(detA), y_interface);
#endif

        if (MyPID == 0) cout << "DONE\n";
        // No sync needed — both solvers share data via sharedState.

        model->printSolution(*soln, t);
        y_interface->HBilinV();
        y_interface->computeEyyTyT();
        model->setExpVyVy(y_interface->getEVyVy());

        if (MyPID == 0 && !isStochOn)
            cout << "\n WARNING: Stochastic classes initialized but Use Stochastic is FALSE.\n";
        printnormMV(*soln, 2, "initial norm of soln");

        Teuchos::RCP<Epetra_MultiVector> expyy = y_interface->getEyy();
        int nv = std::max(1, (int) round(saveEyyIntrvl / dt));
        if (MyPID == 0) cout << "\nEyy snapshots per file: " << nv << "\n";
        int ttlelmnt = expyy->NumVectors() * expyy->GlobalLength();
        Epetra_MultiVector Eyy(Epetra_Map(ttlelmnt, 0, *Comm), nv);
        int eyyBufferIdx = 0;
        double count = t;

        Teuchos::RCP<std::ostream> outstream2;
        if (MyPID == 0) {
            if (std::fabs(t_start) < 1e-6)
                outstream2 = Teuchos::rcp(new std::ofstream("time_dt_norm_Evyvy.txt"));
            else {
                cout << "\nNOTE: Appending to time_dt_norm_Evyvy.txt\n";
                outstream2 = Teuchos::rcp(new std::ofstream("time_dt_norm_Evyvy.txt",
                                                            ios::app));
            }
        } else {
            outstream2 = Teuchos::rcp(new Teuchos::oblackholestream());
        }

        // ===== PHASE 6: Time integration (adaptive, NOX-driven) =====

        INFO("Start Time integration");
        if (MyPID == 0)
            cout << "\n start time = " << time_out(t)
                    << "\n end time   = " << time_out(t_end) << "\n";

        while (t < t_end) {
            if (step > maxsteps) {
                (*outstream) << "Maximum number of steps exceeded, stopping...\n";
                cout << "Maximum number of steps exceeded, stopping...\n";
                break;
            }
            (*outstream) << "\n##############################################\n";
            (*outstream) << "start time step " << step << ": t=" << time_out(t) << "\n";

            if (t + dt > t_end) dt = std::max(dt_min, t_end - t);
            (*outstream) << "step size dt = " << time_out(dt) << "\n";
            (*outstream) << "##############################################\n\n";

            soln_old = *soln;

            //model->createW(t);
            //Wbase = model->get_W();

            runStochStep(isStochOn, timeProf, dt, y_interface, Vstoch, Vn, stochTimers);

            if (timeProf) meanTime->start();
            bool success = stepper->Step(soln_old, t, *soln, dt);
            if (timeProf) {
                meanTime->stop();
                meanTime->incrementNumCalls();
            }

            EpetraExt::MatrixMatrix::Add(*A, false, 1.0, detA, 0.0);
            if (Scaling) {
                detA.RightScale(ScR);
                detA.LeftScale(ScL);
            }

            (*outstream) << "\n##############################################\n";

            double lte = -1.0;
            double effortEst = stepper->getEffortEstimate();
            (*outstream) << "Step effort: " << effortEst * 100 << "%\n";

            if (success && stepper->hasErrorEstimate()) {
                lte = stepper->getErrorEstimate();
                (*outstream) << "Local truncation error tau: " << lte << "\n";
                if (stepChoice == "Adaptive" && lte > lteTol) {
                    if (dt > dt_min) {
                        success = false;
                        (*outstream) << "Step rejected: LTE=" << lte
                                << " > tol=" << lteTol << "\n";
                    } else {
                        (*outstream) << "WARNING: step failed accuracy target but "
                                "step-size is already at minimum.\n";
                    }
                }
            }

            if (success) {
                if (MyPID == 0) std::flush(cout << "+");
                t += dt;
                pVector->setValue("Time", t);
                model->setParameters(*pVector);

                saveTimestepOutputs(t, eyyBufferIdx, nv, Eyy, *expyy,
                                    isStochOn, *Vn, y_interface->getYTrans(),
                                    *soln, prntintvl, count);

                step++;
                (*outstream) << "Step was successful!\n";
                double nrmF, nrmexp;
                CHECK_ZERO(stepper->getF().Norm2(&nrmF));
                (y_interface->getEVyVy())->Norm2(&nrmexp);
                (*outstream2) << t << "\t" << dt << "\t" << nrmF << "\t" << nrmexp << "\n";
            } else {
                if (dt == dt_min)
                    (*outstream) << "Time stepping failed at minimum step-size, stopping.\n";
                else
                    (*outstream) << "Step failed!\n";
                if (dt == dt_min) break;
            }

            // Step-size control
            bool increase_dt = (effortEst < maxEffort) ? success : false;
            if (effortEst >= maxEffort)
                (*outstream) << "Step size not increased: system too hard to solve.\n";

            if (stepChoice == "Constant") {
                if (increase_dt)
                    dt = std::min(dt_max, inc_fac * dt);
                else if (!success)
                    dt = std::max(dt_min, red_fac * dt);
            } else if (stepChoice == "Adaptive") {
                if (lte != -1.0) {
                    if (lte > lteTol || lte < lteTol / 4.0) {
                        double newdt = std::sqrt(lteTol / (2.0 * lte)) * dt;
                        newdt = std::max(std::min(newdt, 2.0 * dt), dt / 2.0);
                        newdt = std::min(newdt, dt_max);
                        newdt = std::max(newdt, dt_min);
                        dt = newdt;
                    }
                } else if (!success) {
                    (*outstream) << "Reducing step size.\n";
                    dt = std::max(dt_min, dt / 2.0);
                }
            }

            // Steady-state convergence check
            if (success && steadyMon == "Norm of F") {
                double nrmF;
                CHECK_ZERO(stepper->getF().Norm2(&nrmF));
                (*outstream) << "Norm of F(u): " << nrmF << "\n";
                if (nrmF < steadyTol) {
                    (*outstream) << "Convergence to steady-state detected, stopping.\n";
                    break;
                }
            }
        }

        // ===== PHASE 7: Post-processing and final output =====

        if (Scaling) {
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

        double elapsed = timer.ElapsedTime();
        (*info) << "\nElapsed: " << (int) (elapsed / 60) << " min\n";
        (*outstream) << "\nTime-stepping run finished.\nElapsed: "
                << (int) (elapsed / 60) << " min\n";

        model->WriteConfiguration("MeanSolnFinal.txt", *pVector, *soln);
        INFO("done!");
        HYMLS::MatrixUtils::mmwrite("V_BASE.mm", *Vn);
        HYMLS::MatrixUtils::mmwrite("YTrans_COEFF.mm", y_interface->getYTrans());

        (*outstream) << "Final Parameters:\n" << *paramList << "\n";

        if (timeProf) {
            std::ofstream timeFile("timeProf.yml");
            Teuchos::TimeMonitor::report(timeFile, "my",
                                         parameterList(*(Teuchos::TimeMonitor::getValidReportParameters())));
        }
    } // end try
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, status);

    Comm->Barrier();
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return status ? EXIT_SUCCESS : EXIT_FAILURE;
}

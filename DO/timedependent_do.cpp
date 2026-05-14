/*
 * =====================================================================================
 *
 *       Filename:  timedependent_do.cpp
 *
 *    Description:  DO time-integration driver for Burgers (brgr=1) and
 *                  Quasi-Geostrophic (quasi_geo=1) targets.  Implements a
 *                  direct Newton time-stepping loop without NOX/LOCA.
 *                  The Mean class is resolved at compile time by Interface.hpp
 *                  according to the brgr / quasi_geo preprocessor definitions.
 *                  need_locaInterface is always 0 for this translation unit.
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
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "EpetraExt_MatrixMatrix.h"
#include "StochIO.hpp"
#include "HYMLS_MatrixUtils.hpp"
#include "HYMLS_Tools.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "globdefs.H"
#include "Interface.hpp"        // resolves Mean → Burger::Mean or QG::Mean
#include "Problem_Interface.hpp"
#include "CoeffSolver.hpp"
#include "DOUtils.hpp"
#include "DOTimeLoop.hpp"

#ifndef INFO
#define INFO(s) std::cout << s << std::endl;
#endif

// Burgers/QG operate in non-dimensionalized simulation time; no scaling needed.
static const double timesc = 1.0;

// =====================================================================================
// MeanSolver concept (C++17: duck-typed via template; upgrade to C++20 Concept when
// Trilinos requires C++20).
//
// Any type passed as MeanType to runSimulation<> must provide:
//   MeanType(RCP<ParameterList>, double* t, double* dt, RCP<Epetra_Comm>)
//   void             refreshForcing(double t)
//   bool             NewtonSolver()
//   void             computeF(Epetra_Vector& u, Epetra_Vector& F)
//   void             computeJacobian(Epetra_Vector& x, Epetra_CrsMatrix& A)
//   Epetra_Vector&   getF()
//   RCP<Epetra_Vector>       getSolution()
//   RCP<Epetra_CrsMatrix>    getJacobian()
//   RCP<Epetra_CrsMatrix>    getMassMatrix()
//   RCP<Epetra_MultiVector>  get_W()
//   void             setEVyVy(RCP<Epetra_Vector>)
//   void             WriteSolution(string, double, const Epetra_Vector&)
// =====================================================================================

// =====================================================================================
/// @brief Run the full DO simulation for a given Mean model type.
///
/// Encapsulates the full simulation pipeline (parameter reading, model construction, stochastic
/// solver initialisation, time loop, post-processing).  The model type is
/// resolved at compile time via the MeanType template parameter; no runtime
/// branching on model type occurs inside this function.
///
/// @tparam MeanType  A type satisfying the MeanSolver concept documented above.
///                   Interface.hpp provides the correct alias (Burger::Mean or
///                   QG::Mean) via its own brgr / quasi_geo compile-time dispatch.
///
/// @param[in] Comm   Epetra communicator (serial or MPI).
/// @param[in] timer  Wall-clock timer started in main() before this call.
// =====================================================================================
template<typename MeanType>
void runSimulation(Teuchos::RCP<Epetra_Comm> Comm, Epetra_Time& timer)
{
    int MyPID = Comm->MyPID();

    // ----- Output stream for time-stepping log -----
    Teuchos::RCP<std::ostream> outstream;
    if (MyPID == 0)
        outstream = Teuchos::rcp(new std::ofstream("timestepping.out"));
    else
        outstream = Teuchos::rcp(new Teuchos::oblackholestream());
    (*outstream) << std::setw(15) << std::setprecision(15);

    // ----- Performance timers -----
    Teuchos::RCP<Teuchos::Time> coeffTime  = Teuchos::TimeMonitor::getNewTimer("Stoch. Coeff. time");
    Teuchos::RCP<Teuchos::Time> coeffTime2 = Teuchos::TimeMonitor::getNewTimer("Bilinear form time");
    Teuchos::RCP<Teuchos::Time> basisTime  = Teuchos::TimeMonitor::getNewTimer("Stoch. Basis time");
    Teuchos::RCP<Teuchos::Time> meanTime   = Teuchos::TimeMonitor::getNewTimer("Mean solve time");
    StochTimers stochTimers = { coeffTime, coeffTime2, basisTime };

    // ===== Load simulation parameters from XML =====

    Teuchos::Ptr<Teuchos::ParameterList> paramListptr = Teuchos::ptr(new Teuchos::ParameterList);
    Teuchos::updateParametersFromXmlFile("params.xml", paramListptr);
    Teuchos::RCP<Teuchos::ParameterList> paramList = rcpFromPtr(paramListptr);

    Teuchos::ParameterList& transParams   = paramList->sublist("Time Stepping");
    Teuchos::ParameterList  MeanSolveList = paramList->sublist("Model Parameters");

    double t_start   = transParams.get("Start Time",           0.0);
    double t_end     = transParams.get("End Time",             1.0);
    double dt        = transParams.get("Step Size",           0.01);
    double dtIncrmnt = transParams.get("time step increment", 1.0e-4);
    bool   timeProf  = transParams.get("Time Profiling",      false);
    double t         = t_start;

    // ===== Construct and initialise the mean-field model =====

    Teuchos::RCP<MeanType> model =
        Teuchos::rcp(new MeanType(rcpFromRef(MeanSolveList), &t, &dt, Comm));
    Teuchos::RCP<Epetra_Vector> soln = model->getSolution();
    Epetra_Vector soln_old(soln->Map(), Epetra_DataAccess::Copy);
    int glen = soln->GlobalLength();

    // QG requires a model-specific initial solution; Burgers starts from zero.
#if quasi_geo == 1
    {
        double Re     = MeanSolveList.get("Reynolds Number", 0.0);
        int stableSol = MeanSolveList.get("Stable Solution", 0);
        if (Re == 40.0 && stableSol == 1) {
            if (MyPID == 0)
                cout << "\n Setting initial solution for Re=40, stable branch.\n";
            for (int i = 0; i < soln->MyLength(); i++)
                (*soln)[i] = exp(cos(M_PI * double(i) / double(glen)));
        } else {
            if (MyPID == 0)
                cout << "\n Setting initial solution for Re=" << Re << " (zero).\n";
            soln->PutScalar(0.0);
        }
    }
#endif

    {
        std::string initSolnFile = MeanSolveList.get("Initial Solution File", "None");
        if (initSolnFile != "None") {
            auto initsol = StochIO::readMV(initSolnFile, soln->Map());
            *soln = *(*initsol)(0);
        }
    }

    // ===== Initialise stochastic basis (V) and forcing (W) =====

    Teuchos::Ptr<Teuchos::ParameterList> stochParamListptr =
        Teuchos::ptr(new Teuchos::ParameterList);
    if (MyPID == 0) cout << "Reading StochasticParams.xml\n";
    Teuchos::updateParametersFromXmlFile("StochasticParams.xml", stochParamListptr);
    Teuchos::RCP<Teuchos::ParameterList> stochParamList = rcpFromPtr(stochParamListptr);

    Teuchos::ParameterList& stochModelList = stochParamList->sublist("StochModel");
    double prntintvl     = stochModelList.get("Save Solution Interval", 1.0);
    double saveEyyIntrvl = stochModelList.get("Save Eyy Interval",      0.4);
    double StochFrcStren = stochModelList.get("StochFrc Strength",      1.0);
    bool   isStochOn     = stochModelList.get("Use Stochastic",         true);

    Teuchos::ParameterList& BasisParams = stochModelList.sublist("Stochastic Basis");
    int numvecV = BasisParams.get("Max. Stoch subspace Dimension", 10);

    Teuchos::ParameterList& CoefParams = stochModelList.sublist("Stochastic System");
    int    Stochit        = CoefParams.get("Stochastic Iterations",         100);
    int    numSubTimeStep = CoefParams.get("Number of Sub Time Steps",        1);
    int    maxnumiter     = CoefParams.get("Max Num of Iter",                10);
    bool   usebacktrack   = CoefParams.get("Use BackTracking",             true);
    double backtrackstep  = CoefParams.get("Max Num of BackTracking Steps",  10);
    double tolRHS         = CoefParams.get("Tolerance of RHS",           10e-6);
    double normRHS        = CoefParams.get("Norm of RHS",                  1.0);

    { Epetra_Vector Fu(*soln); model->computeF(*soln, Fu); }

    Teuchos::RCP<Epetra_CrsMatrix> A       = model->getJacobian();
    Teuchos::RCP<Epetra_CrsMatrix> massmat = model->getMassMatrix();
    if (!A->Filled())       A->FillComplete();
    if (!massmat->Filled()) massmat->FillComplete();

    Epetra_CrsMatrix detA(*A);
    model->computeJacobian(*soln, detA);
    detA.FillComplete();
    dumpInitialJacobian(A, detA, *massmat);

    Teuchos::RCP<Epetra_MultiVector> Wbase = model->get_W();
    int numvecW = Wbase->NumVectors();
    Wbase->Scale(StochFrcStren);

    // ===== Shared stochastic state =====
    auto sharedState = Teuchos::rcp(new StochasticState());
    sharedState->uMean = soln;
    sharedState->A    = Teuchos::rcpFromRef(detA);
    sharedState->W    = Wbase;
    sharedState->bilinearTerm = [&model](const Teuchos::RCP<Epetra_Vector>& u,
                                         const Teuchos::RCP<Epetra_Vector>& v,
                                         const Teuchos::RCP<Epetra_Vector>& uv) {
        model->BilinearTerm(u, v, uv);
    };

    Teuchos::RCP<Problem_Interface> Vstoch =
        Teuchos::rcp(new Problem_Interface(Teuchos::rcpFromRef(detA),
                                           Teuchos::rcpFromRef(BasisParams),
                                           model->getSolution(),
                                           numvecV, &t, dt,
                                           Wbase, massmat, Stochit,
                                           sharedState));
    Vstoch->set_frcStrength(StochFrcStren);
    Vstoch->init_v(t);

    // ===== Initialise stochastic coefficient solver (Y) =====

    Teuchos::RCP<Epetra_MultiVector> Vn = Vstoch->getBasisV();
    printnormMV(*Vn, 2, "initial norm of V");
    if (MyPID == 0) cout << "\nInitializing Stoch. Coeff. class...\n";

    Teuchos::RCP<CoeffSolver> y_interface =
        Teuchos::rcp(new CoeffSolver(Stochit, numSubTimeStep, numvecV,
                                 &dt, Teuchos::rcpFromRef(detA), soln,
                                 Teuchos::null, Vn, Wbase, Comm,
                                 Teuchos::rcpFromRef(CoefParams),
                                 sharedState,
                                 maxnumiter, usebacktrack, backtrackstep,
                                 tolRHS, normRHS));

    if (MyPID == 0) cout << "DONE\n";
    // No sync needed — both solvers share data via sharedState.

    model->WriteSolution("MeanSol_" + std::to_string(t) + ".text", t, *soln);
    y_interface->HBilinV();
    y_interface->computeEyyTyT();
    model->setEVyVy(y_interface->getEVyVy());

#if quasi_geo == 1
    model->setSolution(*soln);
#endif

    if (MyPID == 0 && !isStochOn)
        cout << "\n WARNING: Stochastic classes initialized but Use Stochastic is FALSE.\n";
    printnormMV(*soln, 2, "initial norm of soln");

    Teuchos::RCP<Epetra_MultiVector> expyy = y_interface->getEyy();
    int nv = std::max(1, (int)round(saveEyyIntrvl / dt));
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

    // ===== Time integration loop =====

    INFO("Start Time integration");
    if (MyPID == 0)
        cout << "\n start time = " << t << "\n end time = " << t_end << "\n";

    while (t < t_end) {
        soln_old = *soln;

        model->refreshForcing(t);

        runStochStep(isStochOn, timeProf, dt, y_interface, Vstoch, Vn, stochTimers);

        if (timeProf) meanTime->start();
        bool success = model->NewtonSolver();
        if (timeProf) { meanTime->stop(); meanTime->incrementNumCalls(); }

        EpetraExt::MatrixMatrix::Add(*A, false, 1.0, detA, 0.0);

        (*outstream) << "\n##############################################\n";

        if (!success) {
            (*outstream) << "Time stepping failed, stopping run...\n";
            break;
        }

        if (MyPID == 0) std::flush(cout << "+");
        t  += dt;
        dt  = std::min(dtIncrmnt + dt, t_end - t);

        saveTimestepOutputs(t, eyyBufferIdx, nv, Eyy, *expyy,
                            isStochOn, *Vn, y_interface->getYTrans(),
                            *soln, prntintvl, count);

        double nrmF;
        CHECK_ZERO(model->getF().Norm2(&nrmF));
        (*outstream) << "Step was successful!\n";
        (*outstream2) << t << "\t" << dt << "\t" << nrmF << "\n";
    }

    // ===== Post-processing and final output =====

    Epetra_Vector scnd_mmnt(*soln);
    y_interface->PostProcess(scnd_mmnt);

    double elapsed = timer.ElapsedTime();
    (*info) << "\nElapsed: " << (int)(elapsed / 60) << " min\n";
    (*outstream) << "\nTime-stepping run finished.\nElapsed: "
                 << (int)(elapsed / 60) << " min\n";

    HYMLS::MatrixUtils::Dump(*soln, "MeanSolnFinal.mm");
    INFO("done!");
    HYMLS::MatrixUtils::mmwrite("V_BASE.mm",       *Vn);
    HYMLS::MatrixUtils::mmwrite("YTrans_COEFF.mm", y_interface->getYTrans());

    (*outstream) << "Final Parameters:\n" << *paramList << "\n";

    if (timeProf) {
        std::ofstream timeFile("timeProf.yml");
        Teuchos::TimeMonitor::report(timeFile, "my",
            parameterList(*(Teuchos::TimeMonitor::getValidReportParameters())));
    }
}

// =====================================================================================
// main: MPI setup, global streams, and model dispatch.
// All simulation logic lives in runSimulation<MeanType>().
// Interface.hpp resolves Mean → Burger::Mean (brgr=1) or QG::Mean (quasi_geo=1)
// at compile time, so no #if guard is needed here.
// =====================================================================================

int main(int argc = 0, char* argv[] = NULL)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif
    bool status = true;

    // ===== MPI communicator and global output stream =====

#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Teuchos::RCP<Epetra_SerialComm> Comm = rcp(new Epetra_SerialComm());
#endif
    Epetra_Time timer(*Comm);
    int MyPID = Comm->MyPID();

    if (MyPID == 0) {
        std::ostringstream fname;
        fname << "info_" << MyPID << ".txt";
        info = Teuchos::rcp(new std::ofstream(fname.str().c_str()));
    } else {
        info = Teuchos::rcp(new Teuchos::oblackholestream());
    }
    (*info) << std::setw(15) << std::setprecision(15);

    INFO(*Comm);

    // Mean is the type alias provided by Interface.hpp for this build target.
    // No #if guard needed: the alias already encodes the brgr / quasi_geo choice.
    try {
        runSimulation<MeanSolver>(Comm, timer);
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, status);

    Comm->Barrier();
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return status ? EXIT_SUCCESS : EXIT_FAILURE;
}

/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Uppsala 20012.                  *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.uu.se                                          *
 **********************************************************************/
#ifndef FVM_LOCA_INTERFACE_H
#define FVM_LOCA_INTERFACE_H

#include "EpetraExt_ModelEvaluator.h"
#include "LOCA_Epetra_ModelEvaluatorInterface.H"
#include "LOCA_SaveEigenData_AbstractStrategy.H"
#include "NOX.H"
#include "NOX_Epetra_Interface_Preconditioner.H"

#include "Epetra_Vector.h" 
#include "LOCA_Epetra.H"
#include <fstream>
#include <iostream>
#include <sstream>

#include "Teuchos_Array.hpp"
#include "Ifpack_Preconditioner.h"

#include "HYMLS_Tools.hpp"

//#include "HYMLS_MatrixUtils.H"
//#include "HYMLS_Preconditioner.H"
class Epetra_IntVector;
class OceanOutputXdmf;
namespace LOCA{
  class ParameterVector;
}
const double comp_time_scaling();

namespace FVM
{

class Domain;

//! class to interface between Trilinos and OpenFOAM

/*! This class implements the EpetraExt ModelEvaluator class 
    (which is also compatible with the Thyra package).    
    */
class ModelEvaluator :
 public EpetraExt::ModelEvaluator,
 public Epetra_Object
{

public:

  //! constructor. The model details are taken from the command line. The
  //! contents of the parameter list are only used to set up the solver, so
  //! the "Linear Solver"->"HYMLS" parameter list should be consistent with
  //! the openFOAM model setup passed on the command line.
  ModelEvaluator(Teuchos::RCP<Epetra_Comm> comm, Teuchos::ParameterList& plist);

  //! destructor
  ~ModelEvaluator();

//! \name EpetraExt/Thyra interface 

//@{
  Teuchos::RCP<Domain> getDomain();
  //! get the map of the variable vector [u,v,w,p,T,S]'
  Teuchos::RCP<const Epetra_Map> get_x_map() const;

  //! get the map of the result vector from a function evaluation f(x) (same as the x_map)
  Teuchos::RCP<const Epetra_Map> get_f_map() const;

  //! starting guess
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  //! get the Jacobian
  Teuchos::RCP<Epetra_Operator> create_W() const;

  //! get the names of parameters
  Teuchos::RCP<const Teuchos::Array<std::string> >
  get_p_names(int l) const {return p_names_;}
  
  //! get initial values for parameters
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const
    {return p_values_;}

  //! get map for parameters
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const
    {return p_map_;}

  //!
  InArgs createInArgs() const;
  //! 
  OutArgs createOutArgs() const;


  //! compute F, A, alpha*B+beta*A, ...
  //! See Thyra_ModelEvaluator for info on how to use this
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  Teuchos::RCP<Epetra_Vector> ReadConfiguration(std::string filename, LOCA::ParameterVector& params);

   void WriteConfiguration(std::string filename, const LOCA::ParameterVector& pvec,
            const Epetra_Vector& soln);
  //! get matrix
  Teuchos::RCP<Epetra_CrsMatrix> getJacobian() {return jac_;}
  
  //! get solution vector
  Teuchos::RCP<Epetra_Vector> getSolution() {return sol_;}

  //! get stress terms alphs*tau/h
  // Teuchos::RCP<Epetra_CrsMatrix> getStress() {return stress;}

  //! get parameters
  Teuchos::RCP<LOCA::ParameterVector> getParameterVector() {return pVector_;}
  
  //! get scaling object
  Teuchos::RCP<NOX::Epetra::Scaling> getScaling() {return scaling_;}
 
  void setEVyVy(Teuchos::RCP<Epetra_Vector>EVyVy){
  expv4_=EVyVy;}

  //! get dimension of stoch forcing (constant)
  int get_dim_W() {return nWvec_;} 
  


protected:  

  //! communicator
  Teuchos::RCP<Epetra_Comm> comm_;
  
  //! domain object
  Teuchos::RCP<Domain> domain_;
  
  //! vector of material info (LAND/OCEAN etc in THCM)
  Teuchos::RCP<Epetra_IntVector> material_;
  
  //! parameter names (this is redundant, but needed for the ModelEvaluator)
  Teuchos::RCP<Teuchos::Array<std::string> > p_names_;

  //! parameter values (this is redundant, but needed for the ModelEvaluator)
  Teuchos::RCP<Epetra_Vector> p_values_;
  
  //! simple replicated map for the p_values
  Teuchos::RCP<Epetra_Map> p_map_;
  
  //! global shared parameter list
  Teuchos::ParameterList& paramList_;
  
  //! backup interval
  double backup_interval_;

  //! when was the last backup written?
  double last_backup_;
  
  //! time series output interval
  double output_interval_;

  //! when was the last series output written?
  double last_output_;

  //! when was the last backup written?
  double last_backup;
  



  //! computes some scalar quantities based on the solution
  //! and stores a file that can be used for plotting bifurcation
  //! diagrams.
  void Monitor(Epetra_Vector& x,double conParam,int dump=0);



  //! transform parameter list (read from XML) to parameter vector (continuation params)
  int ReadParameters(Teuchos::ParameterList& List, LOCA::ParameterVector& pVector);

  //! \name model definition and state

  //!@{
  
  //! problem dimension (i.e. 2 or 3D)
  int dim_;

  //! deg. of freedom per grid cell
  int dof_;
  
  //! global grid size
  int nx_, ny_, nz_, ncells_;
    
  //! physical bounds of the global domain
  double xmin_,xmax_,ymin_,ymax_,zmin_,zmax_;
    
  //! loca parameter vector indicating the current model state
  Teuchos::RCP<LOCA::ParameterVector> pVector_;
  
  //!@}

  //! \name linear algebra objects
  //@{
  
  //! map defining the variable ordering (parallel)
  Teuchos::RCP<Epetra_Map> map_;
  
  //! distributed solution and rhs
  Teuchos::RCP<Epetra_Vector> rhs_, sol_;

  //! distributed solution and rhs with overlap between subdomains
  Teuchos::RCP<Epetra_Vector> localRhs_, localSol_;
  Teuchos::RCP<Epetra_Vector> expv4_, localexpv4_;
  
  //! distributed sparse matrix (Jacobian)
  Teuchos::RCP<Epetra_CrsMatrix> jac_;
  
  //! this is our 'mass matrix', a diagonal matrix with ones for u/v and 0 for p
  Teuchos::RCP<Epetra_Vector> mass_;

  //! null space
  Teuchos::RCP<Epetra_MultiVector> nullSpace_;
  
  //! left and right scaling vectors for the linear solver
  mutable Teuchos::RCP<Epetra_Vector> sca_left_, sca_right_;
  
  //! scaling object
  mutable Teuchos::RCP<NOX::Epetra::Scaling> scaling_;
  
  //!@}

  //! which parameter are we 'continuating' in?
  std::string cont_param_;

  //! this is for creating a bifurcation diagram
  int bif_counter_;
  

  //! number of parameters (constant)
  int npar_; 
  
  //! number of parameters (constant)
  int nWvec_; 
  

  //! function to construct the scaling vectors and object
  //! The scaling is done so that the p-rows/cols contain all ones.
  int RecomputeScaling() const;
  
  //! generate the default geometry, a closed box with solid cells around
  //! each subdomain (the interior boundaries are removed by the assembly
  //! process). This function allocates the material_ vector and fills it,
  //! if geomfile!="None" a material mask is read in from a file.
  int CreateGeometry(std::string geomfile);  

  //! construct vertex coordinates in one dimension 
  //! including overlap.
  Teuchos::RCP<Epetra_Vector> Vertices1D(std::string dim);

   //! write 'light-weight' data (parameters and binary filename)
  int write_state(std::string asciifilename, std::string binfilename);

  //! read 'light-weight' data (parameters and binary filename)
  int read_state(std::string asciifilename, std::string& binfilename);

  //! store a vector (currently in ASCAII MatrixMarket)
  int write_vector(const Epetra_MultiVector& vec, std::string filename);

  //! read vector stored by write_vector()
  int read_vector(Epetra_Vector& vec, std::string filename);

  //! label
  std::string label_;
  
  private:

  //! read parameter entry from ascii stream
  void read_parameter_entry(Teuchos::RCP<std::istream>, std::string&, double&);
    

  };

//! this class is an adaptor for LOCA: it inherits everything      
//! from ModelEvaluator and ModelEvaluatorInterface to form a 
//! complete Interface::TimeDependent.                             
//! This interface uses the following formulation of the model:    
//! B dx/dt = f(x,t), W = alpha*df/dx + beta*B                         

class LocaInterface : 
  public ModelEvaluator,
  public LOCA::Epetra::ModelEvaluatorInterface,
  public NOX::Abstract::PrePostOperator,
  public LOCA::SaveEigenData::AbstractStrategy
  {
  public:
  
  //! constructor:
  //! if preclist!=null, a preconditioner will be          
  //! available through computePreconditioner() and getPreconditioner().
  //! plist is the "THCM" sublist, preclist is the                      
  //! "Linear Solver"->"Block Preconditioner" sublist.                  
  LocaInterface(Teuchos::ParameterList& plist, Teuchos::RCP<Epetra_Comm> comm,
        const Teuchos::RCP<LOCA::GlobalData>& globalData,
        Teuchos::RCP<Teuchos::ParameterList> preclist=Teuchos::null);

  // destructor
  ~LocaInterface();


   //! give the model access to the linear system object
  void SetSharedData(Teuchos::Ptr<NOX::Epetra::LinearSystem> linSys,
                     Teuchos::Ptr<Teuchos::ParameterList> globalParams)
    {
    linSys_=linSys;
    sharedParams_=globalParams;
    return;
    }

//! \name NOX/LOCA interface 

//!@{
  
/*!
          \brief Provides data to application for output files.

          This routine is called from Interface::xyzt::printSolution() just
          before the call to Interface::Required::printSolution(x,param),
          and gives the application some indices that can be used for
          creating a unique name/index for the output files.
        */
        void dataForPrintSolution(const int conStep, const int timeStep,
                                  const int totalTimeSteps);
  
  
  //! Call user's own print routine for vector-parameter pair
  void printSolution(const Epetra_Vector& x_,
                           double conParam);
                           
  //! save eigenvalue/-vector data
  NOX::Abstract::Group::ReturnType
      save(Teuchos::RCP< std::vector<double> >& evals_r,
           Teuchos::RCP< std::vector<double> >& evals_i,
           Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
           Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i);

                                         

  //!@}

   //! \name Pre/Post operations
  //@{

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::iterate().
  void runPreIterate(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::iterate().
  void runPostIterate(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::solve().
  void runPreSolve(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::solve().
  void runPostSolve(const NOX::Solver::Generic& solver);

  // LOCA preprocessing step (before each continuation step)
  void preProcessContinuationStep(
                             LOCA::Abstract::Iterator::StepStatus stepStatus,
                             LOCA::Epetra::Group& group);

  // LOCA postprocessing step (after each continuation step)
 // void postProcessContinuationStep(
   //                          LOCA::Abstract::Iterator::StepStatus stepStatus,
     //                        LOCA::Epetra::Group& group);

  //@}
                           
  //! close output files (happens automatically during destruction but
  //! it may have to be done by hand before calling MPI_Finalize)
  //! This means that the solution is no longer written in printSolution().
 // void finishOutput() {xdmfWriter =Teuchos::null; output_interval=-1;}

//!\name NOX Preconditioner interface

//!@{

//! Compute preconditioner \f$M\f$.

//! note: this may cause an error if you didn't give a preclist 
//!       to the constructor.                                   
virtual bool computePreconditioner(const Epetra_Vector& x,
                                         Epetra_Operator& Prec,
                                         Teuchos::ParameterList* p = NULL);                                                                                              
                                                                                              
//!@}

  //! return the preconditioner operator. Will only be non-null    
  //! if you passed a preclist to the constructor. Before using    
  //! the preconditioner, computePreconditioner() should be called.
  Teuchos::RCP<Epetra_Operator> getPreconditioner()
    {
//    return Teuchos::rcp_dynamic_cast<Epetra_Operator>(precPtr_);
    return precPtr_;
    }

  //! get the mass matrix
  Teuchos::RCP<Epetra_CrsMatrix> getMassMatrix() const;

   //! get the null space
  Teuchos::RCP<Epetra_MultiVector> getNullSpace() const {return nullSpace_;}

  //! always backup, regardless of parameter value (mainly used for saving finl config)
  void setForceBackup(bool value){force_backup_=value;}
  
  //! print performance to screen
  void printPerformance() const;

 int write_vector(const Epetra_MultiVector& vec, std::string filename);
protected:

  //! filename of backup file (by default it's "IntermediateConfig.txt", this only
  //! changes in XYZT mode)
 std::string backup_filename_;
  
  //! Object for writing Xdmf data files at regular intervals
  Teuchos::RCP<OceanOutputXdmf> xdmfWriter;
  
  //! this is usually false, but if it's true a backup file is always written
  bool force_backup_;

   //! pointer to the NOX object for linear systems, can be set by SetSharedLinSys()
  Teuchos::Ptr<NOX::Epetra::LinearSystem> linSys_;

  //! pointer to the parameter list shared with NOX/LOCA
  Teuchos::Ptr<Teuchos::ParameterList> sharedParams_;
  
  //! custom preconditioner, constructed when you first call computePreconditioner
  Teuchos::RCP<Ifpack_Preconditioner> precPtr_;
 
  //! time series output interval
  double output_interval;
  
  //! dump all vectors and matrices in the Newton process
  bool dumpLinearSystems_;
 //! use our own eigenvalue strategy after a Newton solve
  bool ownEigStrat_;
  
  //! counters
  int numNewtonSolves_;
  int numNewtonIters_;  
  
private:

  //! call fortran routine to store a vector in binary format
 

  //! call fortran routine to read a vector in binary format
  int read_vector(Epetra_Vector& vec, std::string filename);
  
  //! write 'light-weight' data (parameters and binary filename)
  int write_state(std::string asciifilename, std::string binfilename);

  //! read 'light-weight' data (parameters and binary filename)
  int read_state(std::string asciifilename, std::string& binfilename);


  };

} //namespace Foam2Loca
#endif

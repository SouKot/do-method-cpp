/********************************************************************** * Copyright by Jonas Thies, Univ. of Uppsala 2012.                *  *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.uu.se                                          *
 **********************************************************************/
#include <iostream>
using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <iomanip>
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Epetra_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_BlockMapOut.h"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Epetra_MultiVector.H"

#include "LOCA_Parameter_SublistParser.H"
#include "GaleriExt_Utils.H"
#include "HYMLS_Tools.hpp"
#include "HYMLS_MatrixUtils.hpp"
#include "HYMLS_Preconditioner.hpp"
#include "HYMLS_Solver.hpp"
#include "FVM_Enums.H"
#include "FVM_LocaInterface.H"
#include "FVM_Domain.H"
#include "FVM_model_interface.h"
#include "OceanOutputXdmf.H"
#include "NOX_Epetra_LinearSystem_Hymls.hpp"
//#include "FVM_Eigs.H"

namespace FVM {

// for testing purposes, mostly:
static int step_counter;

ModelEvaluator::ModelEvaluator(Teuchos::RCP<Epetra_Comm> comm, 
                Teuchos::ParameterList& plist):
  comm_(comm),
  paramList_(plist)
  {
  step_counter=0;
  bif_counter_=0;
  this->SetLabel("FVM::ModelEvaluator");
  label_=string(Label());
  
  HYMLS_PROF3(label_,"Constructor");
  
  // construct the domain
  nx_ = paramList_.get("nx",16);
  ny_ = paramList_.get("ny",16);
  nz_ = paramList_.get("nz",16);
  dim_= paramList_.get("dim",3);
  dof_ = paramList_.get("dof",1);
  xmin_=paramList_.get("xmin",0.0);
  xmax_=paramList_.get("xmax",1.0);
  ymin_=paramList_.get("ymin",0.0);
  ymax_=paramList_.get("ymax",1.0);
  zmin_=paramList_.get("zmin",0.0);
  zmax_=paramList_.get("zmax",1.0); 
  // this does most of the parallel setup like domain decomposition,
  // map generation etc.  
  domain_ = Teuchos::rcp(new Domain(nx_,ny_,nz_,dof_,
        xmin_, xmax_, ymin_, ymax_, zmin_, zmax_, 
        comm_));

// do the domain decomposition
bool decomp2D = paramList_.get("2D Domain Decomposition",false);
if (decomp2D)
  {
  CHECK_ZERO(domain_->Decomp2D());
  }
else
  {
  CHECK_ZERO(domain_->Decomp3D());
  }
  // create a parallel map.
  map_=domain_->GetSolveMap();
  Teuchos::RCP<Epetra_Map> localMap=domain_->GetAssemblyMap();    
  //std::flush(std::cout<<"\n!!!!!!*****"<<__FILE__<<" "<<__LINE__);
  //DEBVAR(*map_);
  //std::flush(std::cout<<"\n!!!!!!*****"<<__FILE__<<" "<<__LINE__);
  //DEBVAR(*localMap);

  sol_=Teuchos::rcp(new Epetra_Vector(*map_));
  rhs_=Teuchos::rcp(new Epetra_Vector(*map_));
  expv4_=Teuchos::rcp(new Epetra_Vector(*map_));
  localSol_=Teuchos::rcp(new Epetra_Vector(*localMap));
  localRhs_=Teuchos::rcp(new Epetra_Vector(*localMap));
  // call model's initialize function to allocate memory etc. on each subdomain
  int nloc = domain_->LocalN();
  int mloc = domain_->LocalM();
  int lloc = domain_->LocalL();
  double xminloc = domain_->XminLoc();
  double xmaxloc = domain_->XmaxLoc();
  double yminloc = domain_->YminLoc();
  double ymaxloc = domain_->YmaxLoc();
  double zminloc = domain_->ZminLoc();
  double zmaxloc = domain_->ZmaxLoc();

  std::string geomfile = "None";
  geomfile = paramList_.get("Material Mask File",geomfile);
  CHECK_ZERO(this->CreateGeometry(geomfile));
  
  int* material_ptr;
  CHECK_ZERO(material_->ExtractView(&material_ptr));
  
  // generate coordinates of cell corners and cell centers for a cartesian grid.
  // There are two ways to use these coordinates
  // - simply set "xmin"=0,"xmax=1","cell ratio (x)"=1 etc and define the grid
  // - in the model implementation in Fortran based on the resulting regular grid
  // - use the parameters to define a simple stretched grid in a cartesian domain
  Teuchos::RCP<Epetra_Vector> x_coords=this->Vertices1D("x");
  Teuchos::RCP<Epetra_Vector> y_coords=this->Vertices1D("y");
  Teuchos::RCP<Epetra_Vector> z_coords=this->Vertices1D("z");

 if (y_coords->Map().Comm().MyPID()==0 && z_coords->Map().Comm().MyPID()==0)
    {
    HYMLS::MatrixUtils::Dump(*x_coords,"x_coords.txt",false,HYMLS::MatrixUtils::GATHER);
    }
  if (x_coords->Map().Comm().MyPID()==0 && z_coords->Map().Comm().MyPID()==0)
    {
    HYMLS::MatrixUtils::Dump(*y_coords,"y_coords.txt",false,HYMLS::MatrixUtils::GATHER);
    }
  if (x_coords->Map().Comm().MyPID()==0 && y_coords->Map().Comm().MyPID()==0)
    {
    HYMLS::MatrixUtils::Dump(*z_coords,"z_coords.txt",false,HYMLS::MatrixUtils::GATHER);
    }

  double *x_ptr, *y_ptr, *z_ptr;
  CHECK_ZERO(x_coords->ExtractView(&x_ptr));
  CHECK_ZERO(y_coords->ExtractView(&y_ptr));
  CHECK_ZERO(z_coords->ExtractView(&z_ptr));
    double *sol;
    CHECK_ZERO(localSol_->ExtractView(&sol));
  int ierr=0;
   int nrows = localSol_->MyLength();
  model_init(&nloc,&mloc,&lloc,x_ptr,y_ptr,z_ptr,&nx_,&ny_,&nz_,
        &xminloc,&xmaxloc,&yminloc,&ymaxloc,&xmin_,&xmax_,&ymin_,&ymax_,
        material_ptr, &nrows, sol, &ierr);
  domain_->Assembly2Solve(*localSol_, *sol_);
  Teuchos::RCP<Epetra_Vector> SOL=Teuchos::rcp(new Epetra_Vector(*map_));
  domain_->Assembly2Solve(*localSol_, *sol_);//Teuchos::RCP<Epetra_Vector> SOL(sol_);
  if (ierr) HYMLS::Tools::Error("fortran error initializing model",__FILE__,__LINE__);

  // how many continuation parameters are there in the model?
  model_get_num_params(&npar_);

  // continuation parameter (may be "Time" in transient mode)
  // This is typically set by the main program
  cont_param_ = paramList_.get("Parameter Name","Undefined");
  double start_val = paramList_.get("Parameter Initial Value",0.0);    
 
  pVector_=Teuchos::rcp(new LOCA::ParameterVector());

  Teuchos::ParameterList& startList = paramList_.sublist("Starting Parameters");
  CHECK_ZERO(this->ReadParameters(startList,*pVector_));
  // update parameter names and initial values for ModelEvaluator
  Epetra_SerialComm scomm;
  p_map_ = Teuchos::rcp(new Epetra_Map(npar_,0,scomm));
  p_values_ = Teuchos::rcp(new Epetra_Vector(*p_map_));
  p_names_=rcp(new Teuchos::Array<string>(npar_));
  
  bool status=true;    
  try {
  pVector_->setValue(cont_param_,start_val);
  } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,status);
  if (!status) 
    {
    HYMLS::Tools::Error("'Continuation Parameter' not in 'Starting Parameters' sublist",
        __FILE__,__LINE__);
    }
  for (int i=0;i<npar_;i++)
    {
    (*p_names_)[i] = pVector_->getLabel(i);
    (*p_values_)[i] = pVector_->getValue((*p_names_)[i]);
    }

  backup_interval_ = paramList_.get("Backup Interval",-1.0);
  last_backup_ = start_val-1e-12;// we subtract eps because of the way 
                                  // backuping is treated in XYZT mode. 
                                  // otherwise the initial solution would
                                  // be backuped in that case.
  const double timesc=comp_time_scaling();
  output_interval_=paramList_.get("Output Frequency",-1.0)/timesc;
  last_output_=last_backup_-1.1*output_interval_;
  
  double* par_ptr = pVector_->getDoubleArrayPointer();

  model_setparams(&npar_, par_ptr, &ierr);
  if (ierr) HYMLS::Tools::Error("fortran error setting parameters",__FILE__,__LINE__);

  int nnzPerRow=128; // estimate
  Teuchos::DataAccess c1 = Teuchos::Copy;
  Epetra_DataAccess c2 = static_cast<Epetra_DataAccess>(c1);
  jac_=Teuchos::rcp(new Epetra_CrsMatrix( c2,*map_,nnzPerRow));

  mass_=Teuchos::rcp(new Epetra_Vector(*map_));
  
  sca_left_=Teuchos::rcp(new Epetra_Vector(*map_));
  sca_right_=Teuchos::rcp(new Epetra_Vector(*map_));
  
  CHECK_ZERO(sca_left_->PutScalar(1.0));
  CHECK_ZERO(sca_right_->PutScalar(1.0));

  scaling_=Teuchos::rcp(new NOX::Epetra::Scaling());

  scaling_->addUserScaling(NOX::Epetra::Scaling::Right,sca_right_);
  scaling_->addUserScaling(NOX::Epetra::Scaling::Left,sca_left_);
          
  // we do one evaluation of the Jacobian with a random input vector
  // and random parameters to get the pattern right:
  EpetraExt::ModelEvaluator::OutArgs outargs = this->createOutArgs();
  EpetraExt::ModelEvaluator::InArgs inargs = this->createInArgs();

  model_get_num_stochvectors(&nWvec_);
  LOCA::ParameterVector pCopy = *pVector_;
  
  HYMLS::MatrixUtils::Random(*SOL);
  SOL->PutScalar(1.0);
  for (int i=0;i<npar_;i++)
    {
    (*p_values_)[i]=Teuchos::ScalarTraits<double>::random();
    }

  inargs.set_x(SOL);
  inargs.set_p(0,p_values_);
  inargs.set_beta(1.0);
  inargs.set_alpha(1.0);
  outargs.set_W(jac_);    

  DEBUG(" compute initial matrix");
  status=true;

  try {
    this->evalModel(inargs, outargs);
      } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,status);
 
  if (!status) HYMLS::Tools::Fatal("caught an exception when trying to evaluate model",__FILE__,__LINE__);
  *pVector_ = pCopy;

  for (int i=0;i<npar_;i++)
    {
    (*p_values_)[i] = pVector_->getValue((*p_names_)[i]);
    }
  std::string restartFile = paramList_.get("Restart File","None");
  if (restartFile!="None")
    {
    std::string vecFile;
    CHECK_ZERO(this->read_state(restartFile,vecFile));
    if (vecFile=="undefined")
      {
      HYMLS::Tools::Error("restart failed - xml file has no 'State Vector' entry",
        __FILE__,__LINE__);
      }
    if(comm_->MyPID()==0)
      HYMLS::Tools::out() <<" RESTART FROM: "<<vecFile<<std::endl;
    
    CHECK_ZERO(this->read_vector(*sol_,vecFile));
    }  
  
  double ints, doubles;
  model_memory_estimate(&ints,&doubles);

  Teuchos::RCP<Epetra_Map> MAP;
  MAP=HYMLS::MatrixUtils::Gather(*map_,0);

  return;
  }

  ModelEvaluator::~ModelEvaluator()
  {
  HYMLS_PROF3(label_,"Destructor");
  model_free();
  }

  void ModelEvaluator::read_parameter_entry(RCP<std::istream> in, string& key, double& value)
  {
  int j;
  string tmp;
  (*in) >> j;
  (*in) >> key;
  while (1)
    {
    (*in) >> tmp;
    if (tmp=="=")
      {
      (*in) >> value;
      break;
      }
    else
      {
      key = key + " "+tmp;
      }
    }
  DEBVAR(j);
  DEBVAR(key);
  DEBVAR(value);
  }

  Teuchos::RCP<Epetra_Vector> ModelEvaluator::ReadConfiguration(std::string filename ,LOCA::ParameterVector& pVec)
  {

  Teuchos::RCP<std::istream> in;
  in = Teuchos::rcp(new std::ifstream(filename.c_str()) );
  std::string s1,s2,s3;
  (*in) >> s1;
  DEBVAR(s1);
  if (s1!="LOCA::ParameterVector")
    {
    Error("Error reading start config",__FILE__,__LINE__);
    }

  // read THCM Parameter vector
  int npar;
  (*in) >> s1 >> s2 >> npar >>s3;
  DEBVAR(npar)

  int j;
  string key;
  double value;
  for (int i=0;i<npar;i++)
    {
    read_parameter_entry(in,key,value);
    if (pVec.isParameter(key))
      {
      pVec.setValue(key, value);
      } 
    else
      {
      pVec.addParameter(key,value);
      }
    }

  // read current solution
  Teuchos::RCP<Epetra_Vector> dsoln = getSolution();
  Teuchos::RCP<Epetra_Map> dmap = domain_->GetSolveMap();

  Teuchos::RCP<Epetra_MultiVector> gsoln = HYMLS::MatrixUtils::Gather(*dsoln,0);

  if (comm_->MyPID()==0)
    {
    (*in) >> s1;
    if (s1!="Epetra::Vector")
      {
      INFO("Bad Vector label: should be Epetra::Vector, found "<<s1<<std::endl);
      Error("Error reading start config",__FILE__,__LINE__);
      }
    (*in) >> s1 >> s2 >> s3;
    if (s1+s2+s3!="MyPIDGIDValue")
      {
      Error("Error reading start config",__FILE__,__LINE__);
      }
    int pid,gid;
    double val;
    for (int i=0;i<gsoln->GlobalLength();i++)
      {
      (*in) >> pid >> gid >> val;
      (*gsoln)[0][gid]=val;
      }
    }
    
  Teuchos::RCP<Epetra_MultiVector> tmp = HYMLS::MatrixUtils::Scatter(*gsoln,*dmap);
  *dsoln = *((*gsoln)(0));
  
  try {
    last_backup=pVec.getValue(cont_param_);
    } catch (...) {
    Error("Missing continuation parameter in starting file!",__FILE__,__LINE__);
    }
  return dsoln;
  }

  void ModelEvaluator::WriteConfiguration(std::string filename , const LOCA::ParameterVector& 
                   pVector, const Epetra_Vector& soln)
  {
  Teuchos::RCP<std::ostream> out;
  if (comm_->MyPID()==0)
    {
    out = Teuchos::rcp(new std::ofstream(filename.c_str()) );
    }
  else
    { // dummy stream
    out = Teuchos::rcp(new Teuchos::oblackholestream());
    }
  (*out) << std::setw(15) << std::setprecision(15);
  out->setf(ios::scientific);
  (*out) << pVector;
  (*out) << *(HYMLS::MatrixUtils::Gather(soln,0));
  }
  // compute and store streamfunction in 'fort.7'
  void ModelEvaluator::Monitor(Epetra_Vector& x, double conParam,int dump)
  {
  HYMLS_PROF3(label_,"Monitor");
  bif_counter_++;

  //TODO: output something more meaningful than then vector norms
  double norm1,norm2, normInf;
  CHECK_ZERO(x.Norm1(&norm1));  
  CHECK_ZERO(x.Norm2(&norm2));  
  CHECK_ZERO(x.NormInf(&normInf));  
  
  if (comm_->MyPID()==0)  
    {
    std::ofstream os("bif_data.m",ios::app);

    os << cont_param_ <<"("<<bif_counter_<< ") = "<<conParam  << "; ";
    os << "norm1_sol("<<bif_counter_<<") = " << " " << norm1 << "; ";
    os << "norm2_sol("<<bif_counter_<<") = " << " " << norm2 << "; ";
    os << "normInf_sol("<<bif_counter_<<") = " << " " << normInf << "; ";
    os.close();
    }
  return;
  }

  Teuchos::RCP<Domain> ModelEvaluator::getDomain()
  {
    return domain_;
  }
//////////// EpetraExt::ModelEvaluator interface //////////////////////

  // get the map of our variable vector [u,v,w,p,T,S]'
  Teuchos::RCP<const Epetra_Map> ModelEvaluator::get_x_map() const
    {
    return map_;
    }

///////////////////////////////////////////////////////////////////////

  // get the map of our 'model response' F(u)
  Teuchos::RCP<const Epetra_Map> ModelEvaluator::get_f_map() const
    {
    return map_;
    }
  
///////////////////////////////////////////////////////////////////////

  // get initial guess (all zeros)
  Teuchos::RCP<const Epetra_Vector> ModelEvaluator::get_x_init() const
    {
    return sol_;
    }
    
///////////////////////////////////////////////////////////////////////

  // create the Jacobian
  Teuchos::RCP<Epetra_Operator> ModelEvaluator::create_W() const
    {
    return jac_;
    }
    
///////////////////////////////////////////////////////////////////////

  EpetraExt::ModelEvaluator::InArgs ModelEvaluator::createInArgs() const
    {
  HYMLS_PROF3(label_,"createInArgs");
    EpetraExt::ModelEvaluator::InArgsSetup inArgs;
    inArgs.setModelEvalDescription("2D Driven Cavity");
    inArgs.setSupports(IN_ARG_x,true);
    inArgs.setSupports(IN_ARG_x_dot,true); 
    inArgs.setSupports(IN_ARG_alpha,true); 
    inArgs.setSupports(IN_ARG_beta,true);
    inArgs.setSupports(IN_ARG_t,true);
    inArgs.set_Np(1); // note: there are actually npar_ parameters,
                      // but we store them in a single Epetra_Vector
    return inArgs;
    }

///////////////////////////////////////////////////////////////////////
  
  EpetraExt::ModelEvaluator::OutArgs ModelEvaluator::createOutArgs() const
    {
  HYMLS_PROF3(label_,"createOutArgs");
    EpetraExt::ModelEvaluator::  OutArgsSetup outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(OUT_ARG_f,true);
    outArgs.setSupports(OUT_ARG_W,true);
    // TODO: is this correc? I think I just copied it and left it there...
    outArgs.set_W_properties(
    DerivativeProperties(DERIV_LINEARITY_NONCONST
                        ,DERIV_RANK_FULL,true // supportsAdjoint
                    ));
    return outArgs;
    }
    
///////////////////////////////////////////////////////////////////////
// this function is from the EpetraExt::ModelEvaluator interface, it is
// called by LOCA's computeJacobian, computeF and computeShiftedMatrix 
// functions. Note that if LOCA wants alpha A + beta B (e.g. compute-  
// ShiftedMatrix(alpha,beta) is called, it passes alpha^=-beta and     
// beta^=alpha to this function. The ModelEvaluator interface is some- 
// what outdated and not well-documented, but it was chosen here a     
// long time ago so we stick to it. So, this function computes
// beta^*A-alpha^*B
void ModelEvaluator::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
  {
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  //
  // Get the input arguments
  //
  Teuchos::RCP<const Epetra_Vector> x = inArgs.get_x();
  Teuchos::RCP<const Epetra_Vector> xdot = inArgs.get_x_dot();
  double Alpha = inArgs.get_alpha();
  double beta = inArgs.get_beta();
  const Epetra_Vector &x_dot = *inArgs.get_x_dot();
  double t = inArgs.get_t();
  DEBVAR(t)
  Teuchos::RCP<const Epetra_Vector> p_values=inArgs.get_p(0);
 
  if (t==0.0) t =(*p_values)[0];
  
  if (p_values.get()!=p_values_.get())
    {
    (*p_values_) = (*p_values);
    }

  for (int i=0;i<p_values->MyLength();i++)
    {
    if (pVector_->isParameter((*p_names_)[i]))
      {
      pVector_->setValue((*p_names_)[i],(*p_values_)[i]);
      }
    else
      {
      pVector_->addParameter((*p_names_)[i],(*p_values_)[i]);
      }
    }
 
  // set continuation parameters in the model
  int npar=npar_;
  if (pVector_->length()!=npar_)
    {
    HYMLS::Tools::Warning("number of parameters changed???",__FILE__,__LINE__);
    npar=pVector_->length();
    }
  double *param_ptr = pVector_->getDoubleArrayPointer();
  int ierr=0;
  DEBUG("call model_setparams");
  model_setparams(&npar,param_ptr,&ierr);
  if (ierr) HYMLS::Tools::Error("fortran error setting parameters",
        __FILE__,__LINE__);  

  // Get the output arguments
  //
  Teuchos::RCP<Epetra_Vector> f_out = outArgs.get_f();
  
  Teuchos::RCP<Epetra_Operator>     W_out = outArgs.get_W();
 

  if (Alpha!=0.0)
  {
  // we retrieve the mass matrix every time so that in principle
  // the model can have a parameter-dependence in M as well

  //
  int nrows = localSol_->MyLength();
  Epetra_Vector localMass(localSol_->Map());
  //double Mass[nrows];
  double *mass, *xptr ;
  //Epetra_Comm* rcomm=comm_.getRawPtr();
  //Epetra_Map map1(nrows,0,*rcomm);
  CHECK_ZERO(localMass.ExtractView(&mass));
  CHECK_ZERO(localSol_->ExtractView(&xptr));
  meanMassMatTime->start();
  model_massmat(&nrows, xptr, mass);
  meanMassMatTime->stop();
  meanMassMatTime->incrementNumCalls();
    // discard overlap
    CHECK_ZERO(domain_->Assembly2Solve(localMass,*mass_));
  }

  DEBVAR(Alpha)
  DEBVAR(beta)
  
  bool want_A = ((W_out!=Teuchos::null)&&(beta!=0.0));
  bool want_F = (f_out!=Teuchos::null);
  if (W_out!=Teuchos::null)
    {
    // if any of beta*A+alpha*B is requested and the matrix
    // pattern is defined, zero out the matrix.
    if (jac_->Filled())
      {
      jac_->PutScalar(0.0);
      
      }
    else if (want_A==false)
      {
      // mass matrix requested, but matrix not filled -> Error
      HYMLS::Tools::Error("trying to compute mass matrix before Jacobian, not allowed here",
        __FILE__, __LINE__);
      }
    }
  if (want_F || want_A)
    {
    // import overlap
    CHECK_ZERO(domain_->Solve2Assembly(*x,*localSol_));
    }
 

 if (want_F)
    {
 
    double *rhs_ptr,*x_ptr;
    CHECK_ZERO(localRhs_->ExtractView(&rhs_ptr));
    CHECK_ZERO(localSol_->ExtractView(&x_ptr));
    int nloc = localRhs_->MyLength();
    int ierr=0;

    DEBUG("call model_rhs");
    
    meanRhsTime->start();
    model_rhs(&nloc,rhs_ptr,x_ptr,&ierr);
   step_counter ++ ; 

    if (ierr) HYMLS::Tools::Error("fortran error computing rhs",
        __FILE__,__LINE__);
  
    CHECK_ZERO(domain_->Assembly2Solve(*localRhs_,*rhs_));
    meanRhsTime->stop();
    meanRhsTime->incrementNumCalls();
    //**************ADDING THE EXPV4 VALUE***********************/  
    double nrm; 
    expv4_->Norm2(&nrm);
   /*     int MyPID=jac_->Comm().MyPID();
       	if(MyPID==0)
	  std::cout<<"norm of expv4: "<<nrm;
        rhs_->Norm2(&nrm);
	if(MyPID==0)
	  std::cout<<",  of rhs: "<<nrm;

	rhs_->Update(1.0,*expv4_,1.0);

	rhs_->Norm2(&nrm);
	if(MyPID==0) 
	  std::cout<<",  of rhs+expv4: "<<nrm<<"\n";*/ 
	
    rhs_->Update(1.0,*expv4_,1.0);   
   }
  if (want_A)
    {
    double *x_ptr;
    CHECK_ZERO(localSol_->ExtractView(&x_ptr));
    int ierr=0;
    int nrows = localSol_->MyLength();
    int nzmax = nrows*dof_*27;//let's hope that's enough, otherwise we should
                          // get an ierr/=0 from the fortran code.
    int *rows = new int[nrows+1];
    int *cols = new int[nzmax];
    double *values = new double[nzmax];

    DEBUG("call model_jac");
    meanJacTime->start();
    model_jac(&nrows, &nzmax, values, rows, cols,x_ptr, &ierr);  
    if (ierr!=0)
      {
      HYMLS::Tools::Error("non-zero error code "+Teuchos::toString(ierr)+
        " returned from call model_jac(...)",__FILE__,__LINE__);
      }
    DEBUG("fill matrix in Trilinos");
    DEBVAR(nrows);
    // copy the 'real' matrix rows into the global matrix
    for (int i = 0; i<nrows; i++)
      {
      if (!domain_->IsGhost(i))
        {
        int row = localSol_->Map().GID(i);
        int numentries =  rows[i+1] - rows[i];
#ifdef TESTING
        if (map_->MyGID(row)==false) HYMLS::Tools::Error("invalid row index",__FILE__,__LINE__);
#endif        
        // the rows and cols arrays come from Fortran, so the indexing is 1-based
         for (int j = rows[i]-1; j <  rows[i+1]-1; j++)
          { 
          cols[j] = localSol_->Map().GID(cols[j]-1);
#ifdef TESTING
          if (cols[j] == -1) { 
               HYMLS::Tools::Error("invalid column index",
                __FILE__,__LINE__);}
#endif          
#ifdef DEBUGGING_
          HYMLS::Tools::deb() << "("<<cols[j]<<","<<values[j]<<") ";
#endif          
          }
        //return;
        if (jac_->Filled())
          {
          int ierr = jac_->ReplaceGlobalValues(row, numentries,
                            &(values[rows[i]-1]),&(cols[rows[i]-1]));
          //ierr==3 means not all entries replaced, but we zeroed out jac_ first
          // so it doesn't matter.
          if ((ierr!=0) && (ierr!=3))
            {
            HYMLS::Tools::Error("non-zero error code "+Teuchos::toString(ierr)+
                " returned from call jac_->ReplaceGlobalValues",__FILE__,__LINE__);
            }
          }
        else
          {
        CHECK_ZERO(jac_->InsertGlobalValues(row, numentries,
                            &(values[rows[i]-1]),&(cols[rows[i]-1])));
        }
        }
      }
    
    CHECK_ZERO(jac_->FillComplete());
    meanJacTime->stop();
    meanJacTime->incrementNumCalls();
    DEBVAR(nrows);

    delete [] rows;
    delete [] cols;
    delete [] values;
    }


  if (want_F)
    {
    *f_out = *rhs_;
#ifdef STORE_MATRICES
    HYMLS::MatrixUtils::Dump(*rhs_,"rhs.txt");
#endif
    if (beta!=1.0) CHECK_ZERO(f_out->Scale(beta));

#ifdef TESTING
    double divloc,divglob;
    divloc=0.0;
    
    for (int i=dim_;i<rhs_->MyLength();i+=dof_)
      {
      divloc+= (*rhs_)[i] * (*rhs_)[i];
      }
    CHECK_ZERO(comm_->SumAll(&divloc,&divglob,1));
#endif   
    
    }
    

  if (W_out!=Teuchos::null) 
    {
    // pass it on to the OutArgs:
    Teuchos::RCP<Epetra_CrsMatrix> W = 
    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out,true);
    if (beta!=1.0) CHECK_ZERO(W->Scale(beta)); 
    
     if (Alpha!=0.0)
      {
  	if (comm_->MyPID()==0)
	  cout <<" Adapting diagonal " << std::endl; 
      Epetra_Vector diag(jac_->RowMap());
      CHECK_ZERO(jac_->ExtractDiagonalCopy(diag));
      CHECK_ZERO(diag.Update(-Alpha,*mass_,1.0));
      // for Navier-Stokes this returns 1 for the P-rows
      // and ignores the entry in the mass matrix (which
      // should be zero anyway)
      CHECK_NONNEG(jac_->ReplaceDiagonalValues(diag));
      }
     
      if (W.get()!=jac_.get())
        {
        *W=*jac_;
        }
     
    // compute new scaling
    CHECK_ZERO(this->RecomputeScaling());
    }
  return ;
  }

//////////////////////////////////////////////////////////////////////

// after the Jacobian is available, construct left- and right scaling
// so that it has all ones in the p-rows/cols.
int ModelEvaluator::RecomputeScaling() const
  {
  //HYMLS_PROF2(label_,"RecomputeScaling");

  EPETRA_CHK_ERR(sca_left_->PutScalar(1.0));
  EPETRA_CHK_ERR(sca_right_->PutScalar(1.0));
  
  // the matrix hasx been computed and 
  // distributed, and is now in jac_.    
  int *indices;
  double *values;
  int len;
  int lrid,var;
  
  if (jac_->Filled()==false)
    {
    HYMLS::Tools::Error("matrix not filled!",__FILE__,__LINE__);    
    }
#ifdef TTFW
  
  for (int i=0;i<jac_->NumMyRows();i+=dof_)
    {    
    // loop over velocities and look at p-columns
    for (var=0;var<dof_-1;var++)
      {
      lrid = i+var; // local row ID
      EPETRA_CHK_ERR(jac_->ExtractMyRowView(lrid, len, values, indices));
      // loop over matrix row
      for (int j=0;j<len;j++)
        {
        int gcid = jac_->GCID(indices[j]); // global column ID
        // if this column index is a pressure node
        if (MOD(gcid+1,dof_)==0)
          {
          double val = abs(values[j]);
          if (val>1.0e-6)
            {
            (*sca_left_)[lrid] = val; // NOX computes the reciprocal of this value
            break;
            }
          }
        }
      }

    // loop over p-row and pick scaling values for the divergence operator
    var=dof_-1;
    lrid = i+var;
    EPETRA_CHK_ERR(jac_->ExtractMyRowView(lrid, len, values, indices));
    // loop over matrix row
    for (int j=0;j<len;j++)
      {
      int gcid = jac_->GCID(indices[j]); // global column ID
      // if this column index is a velocity node
      if (MOD(gcid+1,dof_)!=0)
        {
        if (map_->MyGID(gcid))
          {
          double val = abs(values[j]);
          if (val>1.0e-6)
            {
            (*sca_right_)[map_->LID(gcid)] = val; // NOX computes the reciprocal of this value
            }
          }
        }
      }
    }
#endif
  
#ifdef TESTING
    HYMLS::MatrixUtils::Dump(*sca_left_,"left_scaling.txt");
    HYMLS::MatrixUtils::Dump(*sca_right_,"right_scaling.txt");
#endif  
  
  return 0;
  }

int ModelEvaluator::ReadParameters(Teuchos::ParameterList& List, LOCA::ParameterVector& pVector)
  {
  HYMLS_PROF3(label_,"ReadParameters");
  if (pVector.length()!=0)
    {
    HYMLS::Tools::Error("initial pVector should be empty",
        __FILE__,__LINE__);
    }
  DEBVAR(npar_);
  // first get the names of the parameters that the model accepts
  int ierr=0;
  int maxlen=1024;
  char buf[maxlen];
  double default_value;
  for (int i=1;i<=npar_;i++)
    {
    model_get_param_name(&i,&maxlen,buf,&default_value,&ierr);
    if (ierr) HYMLS::Tools::Error("fortran call get_param_name returned error",
        __FILE__,__LINE__);
    std::string pname(buf);
    DEBVAR(pname);
    pVector.addParameter(pname,default_value);
    }

  DEBVAR(npar_);
  DEBVAR(pVector);
  
  // now override the default values using the input list
  for (Teuchos::ParameterList::ConstIterator i=List.begin(); i!=List.end(); i++)
    {
    std::string label = i->first;
    int pos;
    try {
    pos = pVector.getIndex(label);
    } catch(...){pos=-1;}
    
    if (pos<0) 
      {
      //HYMLS::Tools::Warning("Your parameter '"+label+"' is not among the model's"
      //  " valid continuation parameters and will be ignored.",__FILE__,__LINE__);
      }
    else
      {
      double value = List.get(label,0.0);
      pVector.setValue((unsigned int)pos,value);
      }
    }
  return 0;
  }


// create standard geometry (maeterial array) with solid boundary cells around each 
// subdomain and read a 'land mask' if filename!="None".
int ModelEvaluator::CreateGeometry(std::string filename)
  {
  HYMLS_PROF3(label_,"CreateGeometry");
  // the material array has (n+2)x(m+2)x(l+2) entries globally, and on each
  // domain as well

  int I0 = 0; int I1 = nx_+1; 
  int J0 = 0; int J1 = ny_+1;
  int K0 = 0; int K1 = nz_+1;
  

  // create an overlapping distributed map
  int i0 = domain_->FirstI()+1;   // 'grid-style' indexing is 1-based 
  int i1 = domain_->LastI()+1;
  int j0 = domain_->FirstJ()+1; 
  int j1 = domain_->LastJ()+1;
  int k0 = domain_->FirstK()+1;
  int k1 = domain_->LastK()+1;
  
  //add the boundary cells i=0,n+1 etc (this is independent of overlap)
  i0--; i1++; j0--; j1++; k0--; k1++;
  
  DEBUG("create material map with overlap...");  
  Teuchos::RCP<Epetra_Map> matmap_loc = domain_->CreateMap(i0,i1,j0,j1,k0,k1,
    I0,I1,J0,J1,K0,K1,1,*comm_);
  
  int nloc = i1-i0+1;
  int mloc = j1-j0+1;
  int lloc = k1-k0+1;

  material_ = Teuchos::rcp(new Epetra_IntVector(*matmap_loc));
  // now fill the overlapping mmaterial vector
  CHECK_ZERO(material_->PutValue((int)FLUID));
  for (int j=0;j<mloc;j++)
    for (int k=0;k<lloc;k++)
      {
      int lid0 = HYMLS::Tools::sub2ind(nloc,mloc,lloc,1,0,j,k,0);
      int lid1 = HYMLS::Tools::sub2ind(nloc,mloc,lloc,1,nloc-1,j,k,0);
      (*material_)[lid0]=(int)SOLID;
      (*material_)[lid1]=(int)SOLID;
      }
  for (int i=0;i<nloc;i++)
    for (int k=0;k<lloc;k++)
      {
      int lid0 = HYMLS::Tools::sub2ind(nloc,mloc,lloc,1,i,0,k,0);
      int lid1 = HYMLS::Tools::sub2ind(nloc,mloc,lloc,1,i,mloc-1,k,0);
      (*material_)[lid0]=(int)SOLID;
      (*material_)[lid1]=(int)SOLID;
      }
  for (int i=0;i<nloc;i++)
    for (int j=0;j<mloc;j++)
      {
      int lid0 = HYMLS::Tools::sub2ind(nloc,mloc,lloc,1,i,j,0,0);
      int lid1 = HYMLS::Tools::sub2ind(nloc,mloc,lloc,1,i,j,lloc-1,0);
      (*material_)[lid0]=(int)SOLID;
      (*material_)[lid1]=(int)SOLID;
      }
  if (filename!="None")
    {
    HYMLS::Tools::Error("reading material mask not implemented",
        __FILE__,__LINE__);
    }
  return 0;  
  }

Teuchos::RCP<Epetra_Vector> ModelEvaluator::Vertices1D(std::string sdir)
  {
  HYMLS_PROF3(label_,"Vertices1D");
  int dir = (sdir=="x")? 0: (sdir=="y")? 1: (sdir=="z")? 2: -1;
  if (dir==-1) HYMLS::Tools::Error("invalid input to Vertices1D function",
        __FILE__,__LINE__);

  // create a communicator in the coordinate direction 
  Teuchos::RCP<Epetra_Comm> comm1D = domain_->GetProcRow(dir);
  
  int i0, i1;
  
  // lower and upper index bounds (including overlap between subdomains)
  if (dir==0)
    {
    i0=domain_->FirstI();
    i1=domain_->LastI()+1;
    }
  else if (dir==1)
    {
    i0=domain_->FirstJ();
    i1=domain_->LastJ()+1;
    }
  else if (dir==2)
    {
    i0=domain_->FirstK();
    i1=domain_->LastK()+1;
    }

  int num_elts = i1-i0+1;
  int *my_elts = new int[num_elts];
  // TODO: handle periodic BC
  for (int i=0;i<num_elts;i++) my_elts[i]=i0+i;

  // create the map
  Teuchos::RCP<Epetra_Map> map1D = Teuchos::rcp(new Epetra_Map(-1,
            num_elts,my_elts,0,*comm1D));
            
  delete [] my_elts;
  
  DEBVAR(sdir);
  DEBVAR(*map1D);
              
  // now create the stretched coordinates
  Teuchos::ParameterList galeriList; 

  std::string lab="n"+sdir;
  int n = paramList_.get(lab,1);
  galeriList.set("nx",n+1);
  
  lab = sdir+"min";
  double xmin = paramList_.get(lab,0.0);

  lab = sdir+"max";
  double xmax = paramList_.get(lab,0.0);

  galeriList.set("lx",xmax-xmin);
  
  lab = "cell ratio ("+sdir+")";  
  galeriList.set("cell ratio (x)",paramList_.get(lab,1.0));
  
  Teuchos::RCP<Epetra_MultiVector> mv = GaleriExt::Utils::CreateCartesianCoordinates
        ("1D", map1D, galeriList);

  Teuchos::RCP<Epetra_Vector> coords
        = Teuchos::rcp(new Epetra_Vector(*map1D));
  for (int i=0; i<mv->MyLength();i++)
    {
    (*coords)[i] = (*mv)[0][i]+xmin;
    }
  return coords;
  }

int ModelEvaluator::write_vector(const Epetra_MultiVector& vec, std::string filename)
  {
  return HYMLS::MatrixUtils::mmwrite(filename,vec);
  }

// call fortran routine to read vector
int ModelEvaluator::read_vector(Epetra_Vector& vec, std::string filename)
  {
  return HYMLS::MatrixUtils::mmread(filename,vec);
  }

int ModelEvaluator::write_state(std::string asciifile, std::string binfile)
  {
  
  Teuchos::ParameterList stateList;
  stateList.set("State Vector",binfile);

    for (int i=0;i<p_names_->size();i++)
      {
      stateList.set((*p_names_)[i],(*p_values_)[i]);
      }
  
  if (comm_->MyPID()==0)
    {
    writeParameterListToXmlFile(stateList,asciifile);
    }
  return 0;
  }

int ModelEvaluator::read_state(std::string asciifile, std::string& binfile)
  {

  Teuchos::ParameterList stateList;
  Teuchos::updateParametersFromXmlFile(asciifile,Teuchos::ptr(&stateList));

  binfile=stateList.get("State Vector","undefined");

    for (int i=0;i<p_names_->size();i++)
      {
      (*p_values_)[i]=stateList.get((*p_names_)[i],(*p_values_)[i]);
      if (pVector_->isParameter((*p_names_)[i]))
        {
        pVector_->setValue((*p_names_)[i],(*p_values_)[i]);
        }
      else
        {
        pVector_->addParameter((*p_names_)[i],(*p_values_)[i]);
        }
      }
  
  last_backup_ = pVector_->getValue(cont_param_);
  last_output_ = pVector_->getValue(cont_param_);
  
  stateList.unused(std::cerr);
  
  return 0;
  }

///////////////////////////////////////////////////////////////////////

// enhanced interface: Interface (for NOX/LOCA)

LocaInterface::LocaInterface(Teuchos::ParameterList& plist, Teuchos::RCP<Epetra_Comm> comm,
    const Teuchos::RCP<LOCA::GlobalData>& globalData,
    Teuchos::RCP<Teuchos::ParameterList> lsParams)
      : ModelEvaluator(comm,plist),
        LOCA::Epetra::ModelEvaluatorInterface(globalData,rcp(this,false)),
        backup_filename_("Config.bin"), force_backup_(false),
        ownEigStrat_(false),linSys_(Teuchos::null),sharedParams_(Teuchos::null)
  {
  if (lsParams!=Teuchos::null)
    {
    // TODO: setup scaling
    // set parameters
    Teuchos::RCP<Teuchos::ParameterList> hymlsParams
        = Teuchos::rcp(&(lsParams->sublist("HYMLS")),false);
    
    // set the problem definition for the solver
     Teuchos::ParameterList& probList = hymlsParams->sublist("Problem");
    
        if (probList.isParameter("Equations")==false)
      {
      if (dof_==dim_+1)
        {
        probList.set("Equations","Stokes-C");
        }
      else if (dof_==dim_+2)
        {
        probList.set("Equations","Bous-C");
        }
      else
        {
        //HYMLS::Tools::Error("please manually set 'HYMLS'->'Problem'->'Equations' parameter,\n"
        //                    "unable to guess.",__FILE__,__LINE__);
        }
      }
    probList.set("Dimension",dim_);
    probList.set("nx",nx_);
    probList.set("ny",ny_);
    if (dim_>2)
      {
      probList.set("nz",nz_);
      }
    
    DEBUG("Create HYMLS solver");
    precPtr_ = Teuchos::rcp(new HYMLS::Preconditioner(jac_,hymlsParams));
    }
  else
    {
    precPtr_=Teuchos::null;
    }

  ownEigStrat_=paramList_.get("Adaptive Cayley",false);

  dumpLinearSystems_=paramList_.get("Dump all Linear Systems",false);
  numNewtonSolves_=0;
  numNewtonIters_=0;

  nullSpace_=Teuchos::rcp(new Epetra_MultiVector(*map_,1));
  (*nullSpace_)(0)->PutScalar(0.0);
  for (int i=0;i<nullSpace_->MyLength();i++)
    {
    int gid = map_->GID(i);
    if (MOD(gid,dof_)==dim_)
      {
      (*nullSpace_)[0][i]=1.0;
      }
    }

  }

LocaInterface::~LocaInterface()
  {
  }

///////////////////////////////////////////////////////////////////////

// get the mass matrix (can also be dene using evalModel, but
// we're too lazy for that
Teuchos::RCP<Epetra_CrsMatrix> LocaInterface::getMassMatrix() const
  {

  Teuchos::DataAccess c1 = Teuchos::Copy;
  Epetra_DataAccess c2 = static_cast<Epetra_DataAccess>(c1);
  Teuchos::RCP<Epetra_CrsMatrix> massMat = Teuchos::rcp(new 
          Epetra_CrsMatrix(c2, *map_,1,true));
  int gid;
  double val;
  for (int i=0; i<massMat->NumMyRows();i++)
    {
    gid = map_->GID(i);
    val = (*mass_)[i];
    CHECK_ZERO(massMat->InsertGlobalValues(gid, 1, &val,&gid))
    }
  CHECK_ZERO(massMat->FillComplete());
  return massMat;
  }



// compute preconditioner, which can then be retrieved by getPreconditioner()
bool LocaInterface::computePreconditioner(const Epetra_Vector& x,
                                         Epetra_Operator& Prec,
                                         Teuchos::ParameterList* p)
  {
  int result=0;
  DEBUG("enter LocaInterface::computePreconditioner");
  if (precPtr_ == Teuchos::null)
    {
    // no preconditioner parameters passed to constructor
    HYMLS::Tools::Error("No Preconditioner available!",__FILE__,__LINE__);
    }
  
  if (precPtr_->IsInitialized()==false)
    {
    CHECK_ZERO(result=precPtr_->Initialize());
    }

  if (result!=0)
    {
    HYMLS::Tools::Warning("Error code "+Teuchos::toString(result)+" returned when "+
                   " initializing the solver!",__FILE__,__LINE__);
    }
  
  
  if (result==0) 
    {
    Teuchos::RCP<HYMLS::Preconditioner> hymls = 
        Teuchos::rcp_dynamic_cast<HYMLS::Preconditioner>(precPtr_);
    result=precPtr_->Compute();
    if (result!=0)
      {
      HYMLS::Tools::Warning("Error code "+Teuchos::toString(result)+" returned when "+
                   " computing the solver!",__FILE__,__LINE__);
      }
    }
  
  DEBUG("leave LocaInterface::computePreconditioner");
  
  return (result==0);
  }

///////////////////////////////////////////////////////////////////////

// for XYZT output
void LocaInterface::dataForPrintSolution(const int conStep, const int timeStep,
                                  const int totalTimeSteps)
  {    
  }

////////////////
// Call user's own print routine for vector-parameter pair
void LocaInterface::printSolution(const Epetra_Vector& x,
                   double conParam) 
  {
  int dump_psiome=0;
  if (output_interval_>=0)
    {
    if (conParam-last_output_>output_interval_ || force_backup_)
      {
      std::string paramName=Teuchos::StrUtils::varSubstitute(cont_param_," ","_");
      paramName=Teuchos::StrUtils::varSubstitute(paramName,"-","_");
      std::string filename =  
        "solution_"+paramName+"_"+Teuchos::toString((int)conParam)+".mm";
      this->write_vector(x, filename);
      this->write_state("restart.xml",filename);
      last_output_=conParam;
      force_backup_=false;
      dump_psiome=1;
      }
    }

    
  if ((backup_interval_>=0)||force_backup_)
    {
      
    if ((conParam-last_backup_>backup_interval_)||force_backup_||(conParam==last_backup_))
      {
      // NOTE: we can restart from any file written as output
      }
    }
    
  // in every step, we compute the meridional and barotropic streamfunctions
  // and store their maximum in fort.7 (in the old THCM format)
  // At this point, 'grid' contains the complete solution in 3D array format
  this->Monitor(const_cast<Epetra_Vector&>(x),conParam,dump_psiome);

  // invalidate preconditioner after Time-/Continuation step
  }

// implementation of NOX::Abstract::PrePostOperator
// (functions that should be called before and after each nonlinear solve
// and nonlinear solver iteration, respectively)

// executed at the start of a call to iterate()
void LocaInterface::runPreIterate(const NOX::Solver::Generic& solver)
  {
 HYMLS_PROF2(label_,"runPreIterate");
  //HYMLS::Tools::Out("FVM: pre-process Newton step (does nothing right now)");
  if (dumpLinearSystems_)
    {
    std::string filename = 
      "input_"+Teuchos::toString(numNewtonSolves_+1)+"_"+Teuchos::toString(numNewtonIters_+1);
    const NOX::Abstract::Group& group = solver.getSolutionGroup();

    const NOX::Abstract::Vector& solNOX = group.getX();

    try 
     {
     const NOX::Epetra::Vector& solNOXEpetra =
       dynamic_cast<const NOX::Epetra::Vector&>(solNOX);
     const Epetra_Vector& sol = solNOXEpetra.getEpetraVector();
     this->write_vector(sol,filename+".mtx");
     } catch (...)
     {
//     HYMLS::Tools::Warning("could not cast to Epetra vectors, writing in ASCII",
//        __FILE__,__LINE__);
        
     std::ofstream ofs((filename+".txt").c_str());
      ofs<<std::setw(15) << std::setprecision(15);
     solNOX.print(ofs);
     }
    }
  }
    
// executed at the end of a call to iterate()
void LocaInterface::runPostIterate(const NOX::Solver::Generic& solver)
  {
  HYMLS_PROF2(label_,"runPostIterate");
  //HYMLS::Tools::Out("FVM: post-process Newton step (does nothing right now)");
  if (dumpLinearSystems_)
    {
    std::string filename1 = 
      "result_"+Teuchos::toString(numNewtonSolves_+1)+"_"+Teuchos::toString(numNewtonIters_+1);
    std::string filename2 = 
      "rhs_"+Teuchos::toString(numNewtonSolves_+1)+"_"+Teuchos::toString(numNewtonIters_+1);
    std::string filename3 = 
      "upd_"+Teuchos::toString(numNewtonSolves_+1)+"_"+Teuchos::toString(numNewtonIters_+1);
    std::string filename4 = 
      "jac_"+Teuchos::toString(numNewtonSolves_+1)+"_"+Teuchos::toString(numNewtonIters_+1);
      std::string filename5 = "massMatrix";
    const NOX::Abstract::Group& group = solver.getSolutionGroup();
    const NOX::Abstract::Group& old_group = solver.getPreviousSolutionGroup();
    const NOX::Abstract::Vector& solNOX = group.getX();
    const NOX::Abstract::Vector& rhsNOX = old_group.getF();
    const NOX::Abstract::Vector& updNOX = group.getNewton();

    try 
     {
     const NOX::Epetra::Vector& solNOXEpetra = 
        dynamic_cast<const NOX::Epetra::Vector&>(solNOX);
     const Epetra_Vector& sol = solNOXEpetra.getEpetraVector();
     
     const NOX::Epetra::Vector& rhsNOXEpetra = 
        dynamic_cast<const NOX::Epetra::Vector&>(rhsNOX);
     const Epetra_Vector& rhs = rhsNOXEpetra.getEpetraVector();
     
     const NOX::Epetra::Vector& updNOXEpetra = 
        dynamic_cast<const NOX::Epetra::Vector&>(updNOX);
     const Epetra_Vector& upd = updNOXEpetra.getEpetraVector();

    this->write_vector(sol,filename1+".mtx");
    this->write_vector(rhs,filename2+".mtx");
    this->write_vector(upd,filename3+".mtx");
    } catch (...)    
      {      
      //HYMLS::Tools::Warning("could not cast to Epetra vectors, writing in ASCII",
      //  __FILE__,__LINE__);
        
      std::ofstream ofs1((filename1+".txt").c_str());
      ofs1<<std::setw(15) << std::setprecision(15);
      solNOX.print(ofs1);
      ofs1.close();
      std::ofstream ofs2((filename2+".txt").c_str());
      ofs2<<std::setw(15) << std::setprecision(15);
      rhsNOX.print(ofs2);
      ofs2.close();
      std::ofstream ofs3((filename3+".txt").c_str());
      ofs3<<std::setw(15) << std::setprecision(15);
      updNOX.print(ofs3);
      ofs3.close();
      }
   
    HYMLS::MatrixUtils::Dump(*jac_,filename4+".txt",false);    
    HYMLS::MatrixUtils::Dump(*mass_,filename5+".txt",false);       
    }
  numNewtonIters_++;
  }
        
// executed at the start of a call to solve()
void LocaInterface::runPreSolve(const NOX::Solver::Generic& solver)
  {
  HYMLS_PROF2(label_,"runPreSolve");
  //HYMLS::Tools::Out("FVM: pre-process Newton solve (does nothing right now)");
  numNewtonIters_=0;
  }

// executed at the end of a call to solve()
void LocaInterface::runPostSolve(const NOX::Solver::Generic& solver)
  {
HYMLS_PROF2(label_,"runPostSolve");
  //HYMLS::Tools::Out("FVM: post-process Newton solve");
  numNewtonSolves_++;
  if (const_cast<NOX::Solver::Generic&>(solver).getStatus() == NOX::StatusTest::Unconverged)
    {
    if (linSys_!=Teuchos::null)
      {
      // destroy the preconditioner, assuming that it is too old
      HYMLS::Tools::Out("destroying preconditioner because Newton failed to converge");
      linSys_->destroyPreconditioner();
      }
    else
      {
      HYMLS::Tools::Warning("Newton did not converge, and I cannot destroy \n"
                            "the preconditioner because you did not give me\n"
                            "the pointer to the linear system object.",__FILE__,__LINE__);
      }
    }
  return;
  }

void LocaInterface::preProcessContinuationStep(
                             LOCA::Abstract::Iterator::StepStatus stepStatus,
                             LOCA::Epetra::Group& group)
  {
  //HYMLS::Tools::Out("FVM: pre-process continuation step (does nothing right now)");
  // do nothing
  }
// This function is called after LOCA/Anasazi computes eigenvalues
// at the end of a continuation step to store the eigenpairs in some
// way. We use this occasion to adjust the Cayley parameters for the
// next step if "Adaptive Cayley" is set in the "Model" sublist.
NOX::Abstract::Group::ReturnType
      LocaInterface::save(Teuchos::RCP< std::vector<double> >& evals_r,
           Teuchos::RCP< std::vector<double> >& evals_i,
           Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
           Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i)
  {
  HYMLS_PROF2(label_,"save eigs");
  if (comm_->MyPID()==0)
    {
    std::ofstream ofs("Eigenvalues.txt",ios::app);
    ofs << std::scientific << std::setw(8) << std::setprecision(8);
    
    ofs << "***************************************************"<<std::endl;
    ofs << "* Parameters                                      *"<<std::endl;
    ofs << "***************************************************"<<std::endl;
    for (int i=0;i<p_names_->size();i++)
      {
      ofs << (*p_names_)[i] << "\t" << (*p_values_)[i] << std::endl;
      }
    ofs << "***************************************************"<<std::endl;
    ofs << "* Eigenvalues                                     *"<<std::endl;
    ofs << "***************************************************"<<std::endl;
    for (int i=0;i<evals_r->size();i++)
      {
      ofs << (*evals_r)[i] << " + " << (*evals_i)[i]<<"i"<<std::endl;
      }
    ofs << "***************************************************"<<std::endl;
    
    ofs.close();    
    }
  // to avoid huge text files we overwrite the eigenvectors every time

      const NOX::Epetra::MultiVector& evrNOXEpetra = 
        dynamic_cast<const NOX::Epetra::MultiVector&>(*evecs_r);
     const Epetra_MultiVector& evr = evrNOXEpetra.getEpetraMultiVector();

  const NOX::Epetra::MultiVector& eviNOXEpetra = 
        dynamic_cast<const NOX::Epetra::MultiVector&>(*evecs_i);
     const Epetra_MultiVector& evi = eviNOXEpetra.getEpetraMultiVector();
  
  HYMLS::MatrixUtils::Dump(evr,"eigenVectors_re.txt");
  HYMLS::MatrixUtils::Dump(evi,"eigenVectors_im.txt");
  
  if (ownEigStrat_)
    {
    if (sharedParams_==Teuchos::null)
      {
      HYMLS::Tools::Error("cannot proceed, shared param list not set\n"
                          "You have to call setSharedData() in main program\n",
                          __FILE__,__LINE__);
      }

  // note: we inherit globalData from LOCA's ModelEvaluatorInterface
  LOCA::Parameter::SublistParser parsedParams(globalData);
  parsedParams.parseSublists(Teuchos::rcp(sharedParams_.get(),false));

  Teuchos::RCP<Teuchos::ParameterList> eigenParams = 
    parsedParams.getSublist("Eigensolver");
    string opType=eigenParams->get("Operator","Cayley");
    if (opType!="Cayley")
      {
      HYMLS::Tools::Error("'Adaptive Cayley' set with 'Operator' not 'Cayley'\n"
                          "but '"+opType+"'",
        __FILE__,__LINE__);
      }
    HYMLS::Tools::Out("***************************************************");
    HYMLS::Tools::Out("* Adaptive Cayley strategy                        *");
    HYMLS::Tools::Out("***************************************************");
    // the transformed operator is defined as (A-sigma B)\(A-mu B), where
    // sigma is the 'Cayley Pole' and mu the 'Cayley Zero'. Our strategy 
    // is to choose the centerline c between the two such that k eigenvalues
    // are between c and 0. Those eigenvalues are mapped to positions outside
    // the unit circle and will converge quickly in the next continuation step.
    double cayleyPole = eigenParams->get("Cayley Pole",0.0);
    double cayleyZero = eigenParams->get("Cayley Zero",-1.0);
    int numEigs = eigenParams->get("Num of Eigenvalues",10);
    HYMLS::Tools::out() << "current Cayley pole: "<<cayleyPole<<std::endl;
    HYMLS::Tools::out() << "current Cayley zero: "<<cayleyZero<<std::endl;
    
    HYMLS::Tools::Out("---------------------------------------------------");
    HYMLS::Tools::Out("| Current Eigenvalues                             |");
    HYMLS::Tools::Out("---------------------------------------------------");
    for (int i=0;i<evals_r->size();i++)
      {
      HYMLS::Tools::out() << (*evals_r)[i] << " + " << (*evals_i)[i]<<"i"<<std::endl;
      }
    HYMLS::Tools::Out("---------------------------------------------------");

    // select new shifts so that k eigenvalues are between c and 0.
    double c = 0.0;
    int count=0;
    double max_im=-1.0e99;
    double max_neg_re=-1.0e99;
    double min_re=1.0e99;
    for (int i=0;i<evals_r->size();i++)
      {
      if ((*evals_r)[i]<0.0)
        {
        max_neg_re=std::max(max_neg_re,(*evals_r)[i]);
        max_im=std::max(max_im,(*evals_i)[i]);
        }
      min_re = std::min(min_re,(*evals_r)[i]);
      }
    HYMLS::Tools::out() <<"smallest real part: "<<min_re<<std::endl;
    HYMLS::Tools::out() <<"largest neg. real part: "<<max_neg_re<<std::endl;
    HYMLS::Tools::out() <<"largest stable imaginary part: "<<max_im<<std::endl;
    
    cayleyPole = max_im;
    // We use the smallest computed eigenvalue as c,
    // where c is the centerline between sigma (pole) and mu (zero):
    //min_eig = (cayleyPole+cayleyZero)/2;
    cayleyZero=2.0*min_re-cayleyPole;
    
    if (numEigs<evals_r->size())
      {
      HYMLS::Tools::Out(" increase search range because not all eigs converged.");
      cayleyZero*=1.25;
      }
    
    HYMLS::Tools::out() <<"new Cayley pole: "<<cayleyPole<<std::endl;
    HYMLS::Tools::out() <<"new Cayley zero: "<<cayleyZero<<std::endl;
    HYMLS::Tools::Out("***************************************************");

    eigenParams->set("Cayley Pole",cayleyPole);
    eigenParams->set("Cayley Zero",cayleyZero);
    }

  return NOX::Abstract::Group::Ok;
  }

// call fortran routine to print vector
int LocaInterface::write_vector(const Epetra_MultiVector& vec, std::string filename)
  {
  return HYMLS::MatrixUtils::mmwrite(filename,vec);
  
  }

// call fortran routine to read vector
int LocaInterface::read_vector(Epetra_Vector& vec, std::string filename)
  {
   return HYMLS::MatrixUtils::mmread(filename,vec);
  }

int LocaInterface::write_state(std::string asciifile, std::string binfile)
  {
  
  Teuchos::ParameterList stateList;
  stateList.set("State Vector",binfile);

    for (int i=0;i<p_names_->size();i++)
      {
      stateList.set((*p_names_)[i],(*p_values_)[i]);
      }
  
  if (comm_->MyPID()==0)
    {
    writeParameterListToXmlFile(stateList,asciifile);
    }
  return 0;
  }

int LocaInterface::read_state(std::string asciifile, std::string& binfile)
  {

  Teuchos::ParameterList stateList;
  Teuchos::updateParametersFromXmlFile(asciifile,inoutArg(stateList));

  binfile=stateList.get("State Vector","undefined");

    for (int i=0;i<p_names_->size();i++)
      {
      (*p_values_)[i]=stateList.get((*p_names_)[i],(*p_values_)[i]);
      if (pVector_->isParameter((*p_names_)[i]))
        {
        pVector_->setValue((*p_names_)[i],(*p_values_)[i]);
        }
      else
        {
        pVector_->addParameter((*p_names_)[i],(*p_values_)[i]);
        }
      }
  
  last_backup_ = pVector_->getValue(cont_param_);
  last_output_ = pVector_->getValue(cont_param_);
  
  stateList.unused(std::cerr);
  
  return 0;
  }

  
}
///////////////////////////////////////////////////////////////////////

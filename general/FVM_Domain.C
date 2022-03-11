/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/

#include "FVM_Domain.H"
#include <fstream>
//#include <iostream>
#include "HYMLS_Tools.hpp"
#include "HYMLS_MatrixUtils.hpp"

#include "Epetra_DistObject.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif

namespace FVM {

      /* Constructor
          input: dimensions of the global box:
            - N: east-west
            - M: north-south
            - L: z-direction
      */
Domain::Domain(int N, int M, int L, int dof, 
                       double Xmin, double Xmax, 
                       double Ymin, double Ymax, 
                       double Zmin, double Zmax, 
                       Teuchos::RCP<Epetra_Comm> Comm)
  : m(M),n(N),l(L),xmin(Xmin),xmax(Xmax),ymin(Ymin),ymax(Ymax),zmin(Zmin),zmax(Zmax),
    dof_(dof),
    useDifferentSolveMap_(false),
    periodic(false),
    comm(Comm),
    label_("Domain")
  {
  HYMLS_PROF3(label_,"Constructor");
  //std::cerr << "n = "<<n<<std::endl<<std::flush;
  HYMLS_DEBUG(n);
  HYMLS_DEBUG(m);
  HYMLS_DEBVAR(l);
  HYMLS_DEBVAR(dof_);
  HYMLS_DEBVAR(xmin);
  HYMLS_DEBVAR(xmax);
  HYMLS_DEBVAR(ymin);
  HYMLS_DEBVAR(ymax);
  HYMLS_DEBVAR(zmin);
  HYMLS_DEBVAR(zmax);
  //TODO: the do we really need a globally collected 'colmap'?
  /*
  int dim = m*n*l*dof_;
  int *MyGlobalElements = new int[dim];
  for (int i=0;i<dim;i++) MyGlobalElements[i]=i;
  ColMap = Teuchos::rcp(new Epetra_Map(dim,dim,MyGlobalElements,0,*comm) );
  delete [] MyGlobalElements;
  */
  }
      
// Destructor
Domain::~Domain()
  {
  HYMLS_PROF3(label_,"~Domain");
  // destructor handled by Teuchos::rcp's
  }
      
                  
// decompose the domain for a 2D processor array.
// The xy-directions are split up.
int Domain::Decomp2D()
  {
  HYMLS_PROF2(label_,"Decomp2D");
  int nprocs = comm->NumProc();
  int pid = comm->MyPID();

////////////////////////////////////////////////////////////////////

  npL = 1;

// Factor the number of processors into two dimensions. (nprocs = npN*npM)

  int t1 = nprocs;
  int t2 = 1;
  npM=t1;
  npN=t2;

  double r;//remainder
  double r_min = 100;

  while (t1>0)
    {
    t2=(int)(nprocs/t1);
    r = std::abs(m/t1-n/t2);
    if (t1*t2==nprocs && r<=r_min)
      {
      r_min=r;
      npM=t1;
      npN=t2;
      }
//    (*info) << t1 << " " << t2 << " " << r << std::endl;
    t1--;
    }

  HYMLS::Tools::out()<<"+++ Domain decomposition +++"<<std::endl;
  HYMLS::Tools::out()<<"factoring, np = "<<nprocs<<std::endl;
  HYMLS::Tools::out()<<" n = "<< n<<std::endl;
  HYMLS::Tools::out()<<" m = "<< m<<std::endl;
  HYMLS::Tools::out()<<" npN = "<< npN<<std::endl;
  HYMLS::Tools::out()<<" npM = "<< npM << std::endl<<std::endl;

                        
  // find out where in the domain we are situated.
  // the subdomains are numbered in a row-major 'matrix' fashion
  //                  
  // P0 P1 P2  |      
  //           m,pidM 
  //           |      
  // P3 P4 P5  v      
  //  --n-->          
  //   pidN           
  //                  
  
  // note that this corresponds to the column-major ordering in Fortran,
  // i.e. i (n-direction) is the fastest index, and k (l-dir.) the slowest
        
  pidL = 0; // for reasons of generality
  pidN = pid%npN;
  pidM = (pid-pidN)/npN;

  // dimension of actual subdomain (without ghost-nodes)
        
  mloc0 = (int)(m/npM);
  nloc0 = (int)(n/npN);
  lloc0 = l;
  
  // offsets for local->global index conversion
  Loff0 = 0;
  Moff0 = pidM*(int)(m/npM);
  Noff0 = pidN*(int)(n/npN);

  // distribute remaining points among first few cpu's
  int remM=m%npM;
  int remN=n%npN;
  
  if (pidM<remM) mloc0++;
  if (pidN<remN) nloc0++;
  for (int i=0;i<std::min(remM,pidM);i++) Moff0++;
  for (int i=0;i<std::min(remN,pidN);i++) Noff0++;

  //subdomain dimensions/offsets including ghost-nodes
  // (will be added further down)
  mloc=mloc0;
  nloc=nloc0;
  lloc=lloc0;

  Loff = Loff0;
  Moff = Moff0;
  Noff = Noff0;

  // we need two layers of overlap (ghost-nodes) between subdomains:
  //  _________............+                                    
  // |      :  |           :                                    
  // | SD11 :  |  SD12     :                                    
  // |      :  |           :                                    
  // +---------+-----------+                                    
  // (for each variable, obviously)                             
  //							        
  // We need two because THCM assumes one layer of 'LAND' cells 
  // at the domain boundaries. These are discarded only when    
  // assembling global distributed vectors/matrices.            
  
  Teuchos::RCP<Epetra_Comm> xcomm = this->GetProcRow(0);
  xparallel = (xcomm->NumProc()>1); // if there is only one subdomain in the 
                                        // x-direction, periodicity is left to THCM

  if (pidM>0)      { mloc+=num_ghosts; Moff-=num_ghosts;}
  if (pidM<npM-1) { mloc+=num_ghosts;}
  if ((pidN>0)||(periodic&&xparallel))      { nloc+=num_ghosts; Noff-=num_ghosts;}
  if ((pidN<npN-1)||(periodic&&xparallel)) { nloc+=num_ghosts;}

// in the case of periodic boundary conditions the offsets may now be 
// negative or the local domain may exceed the global one. when       
// using nloc and noff, we therefore have to take mod(i,nglob)

  return CommonSetup();		
  }

// 3D domain decomposition (for Navier-Stokes)
// This implementation is more modern than the Decomp2D function
// because we have lots of functionality from HYMLS we can use now.
int Domain::Decomp3D()
  {
  HYMLS_PROF2(label_,"Decomp3D");
  int nprocs = comm->NumProc();
  int pid = comm->MyPID();
  HYMLS::Tools::SplitBox(n,m,l,nprocs,npN,npM,npL);
  HYMLS::Tools::ind2sub(npN,npM,npL,pid,pidN,pidM,pidL);

  Loff0 = pidL*(int)(l/npL);
  Moff0 = pidM*(int)(m/npM);
  Noff0 = pidN*(int)(n/npN);
       
 Noff = Noff0; 
 Moff = Moff0; 
 Loff = Loff0;

  mloc0 = (int)(m/npM);
  nloc0 = (int)(n/npN);
  lloc0 = (int)(l/npL);
 
 nloc = nloc0;
 mloc = mloc0;
 lloc = lloc0;

  if (pidL>0)      { lloc+=num_ghosts; Loff-=num_ghosts;}
  if (pidL<npL-1) { lloc+=num_ghosts;}
  if (pidM>0)      { mloc+=num_ghosts; Moff-=num_ghosts;}
  if (pidM<npM-1) { mloc+=num_ghosts;}
  if ((pidN>0)||(periodic&&xparallel))      { nloc+=num_ghosts; Noff-=num_ghosts;}
  if ((pidN<npN-1)||(periodic&&xparallel)) { nloc+=num_ghosts;}

  return CommonSetup();  
  }

  int Domain::CommonSetup()
    {
    HYMLS_PROF2(label_,"CommonSetup");

//HYMLS::Tools::out()<<"processor position: (N,M,L) = ("<<pidN<<","<<pidM<<","<<pidL<<")"<<std::endl;
//HYMLS::Tools::out()<<"subdomain offsets: "<<Noff0<<","<<Moff0<<","<<Loff0<<std::endl;
//HYMLS::Tools::out()<<"grid dimension on subdomain: "<<nloc0<<"x"<<mloc0<<"x"<<lloc0<<std::endl;
HYMLS_DEBUG("processor position: (N,M,L) = ("<<pidN<<","<<pidM<<","<<pidL<<")");
HYMLS_DEBUG("subdomain offsets: "<<Noff0<<","<<Moff0<<","<<Loff0);
HYMLS_DEBUG("grid dimension on subdomain: "<<nloc0<<"x"<<mloc0<<"x"<<lloc0);
HYMLS_DEBUG("+++ including ghost nodes: +++");
HYMLS_DEBUG("subdomain offsets: "<<Noff<<","<<Moff<<","<<Loff);
HYMLS_DEBUG("grid dimension on subdomain: "<<nloc<<"x"<<mloc<<"x"<<lloc);


  // create the maps:
       
  StandardMap = CreateStandardMap(dof_);
  AssemblyMap = CreateAssemblyMap(dof_);
  // no load-balancing object available, yet (has to be set by user)
  SolveMap = StandardMap;


#ifdef HYMLS_DEBUGGING
comm->Barrier();
HYMLS_DEBUG("create importers...");
#endif

  // finally make the Import/Export objects (transfer function
  // between the two maps)
  as2std = Teuchos::rcp(new Epetra_Import(*AssemblyMap,*StandardMap));
  std2sol = Teuchos::null;
/*
  HYMLS::Tools::out()<<"importer: "<<*as2std);
  Epetra_Vector test_as(*AssemblyMap);
  Epetra_Vector test_so(*SolveMap);
  
  for (int i=0;i<test_as.MyLength();i++)
    {
    test_as[i]=i;
    }
  this->Assembly2Solve(test_as,test_so);
  HYMLS::Tools::out()<<"assembly vector: "<<test_as);
  HYMLS::Tools::out()<<"solve vector: "<<test_so);
  this->Solve2Assembly(test_so,test_as);
  HYMLS::Tools::out()<<"new assembly vector: "<<test_as);
  

//  HYMLS_DEBUG("Importer: "<<std::endl)
//  HYMLS_DEBUG( (*Importer) );
*/

  // determine the physical bounds of the subdomain
  // (must be passed to THCM)
  
  // grid constants as computed in 'grid.f':
  double dx = (xmax-xmin)/n;
  double dy = (ymax-ymin)/m;
  double dz = (zmax-zmin)/l;
        
  xmin_loc = xmin + Noff*dx;
  xmax_loc = xmin + (Noff+nloc)*dx;
  ymin_loc = ymin + Moff*dy;
  ymax_loc = ymin + (Moff+mloc)*dy;
  zmin_loc = zmin + Loff*dz;
  zmax_loc = zmin + (Loff+lloc)*dz;

HYMLS_DEBVAR(xmin)
HYMLS_DEBVAR(xmax)
HYMLS_DEBVAR(ymin)
HYMLS_DEBVAR(ymax)
HYMLS_DEBVAR(xmin_loc)
HYMLS_DEBVAR(xmax_loc)
HYMLS_DEBVAR(ymin_loc)
HYMLS_DEBVAR(ymax_loc)
    return 0;
    }

  // find out wether a particular local index is on a ghost node
  bool Domain::IsGhost(int ind,int dof) const
    {
    int i,j,k,xx;
    int row = ind; 

    int nun = dof>0? dof:dof_;

    bool result = false;


//    HYMLS_DEBUG("GHOST CHECK: li (C++) = "<<ind);

    // find out where the point is
    HYMLS::Tools::ind2sub(nloc,mloc,lloc,nun,row,i,j,k,xx);
    
    // ghost nodes at periodic boundary only if more than one proc in x-direction
//    bool perio = periodic&&xparallel;
    bool perio=false; //TODO
    //HYMLS_DEBUG("GHOST CHECK: i,j,k,xx (C++) = "<<i<<", "<<j<<", "<<k<<", "<<xx);

    for (int ii=0;ii<num_ghosts;ii++)
      {
      result = result||(i==ii && ((pidN>0)||perio));
      result = result||(i==nloc-1-ii && ((pidN<npN-1)||perio));
      result = result||(j==ii && pidM>0);
      result = result||(j==mloc-1-ii && pidM<npM-1);
      result = result||(k==ii && pidL>0);
      result = result||(k==lloc-1-ii && pidL<npL-1);
      }

    //HYMLS_DEBVAR(result);
    return result;
    }

void Domain::SetSolveMap(Teuchos::RCP<Epetra_Map> const map)
  {
  useDifferentSolveMap_ = true;
  SolveMap = map;
  std2sol = Teuchos::rcp(new Epetra_Import(*StandardMap, *SolveMap));
  }

// public map creation function (can only create a limited range of maps)
Teuchos::RCP<Epetra_Map> Domain::CreateSolveMap(int nun, bool depth_av) const
  {
  HYMLS_PROF3(label_,"CreateSolveMap");
  Teuchos::RCP<Epetra_Map> M=Teuchos::null;
  if (!UseLoadBalancing()) // no load-balancing? use our own decomposition
    {
    M=CreateStandardMap(nun,depth_av);
    }
  else
    {
    HYMLS::Tools::Error("not implemented",__FILE__,__LINE__);
//    int l_=depth_av?1:l;
//    M = loadbal->createBalancedMap(l_,nun);
    }
  return M;
  }

Teuchos::RCP<Epetra_Map> Domain::CreateStandardMap(int nun, bool depth_av) const
  {
  HYMLS_PROF3(label_,"CreateStandardMap");
  Teuchos::RCP<Epetra_Map> M = Teuchos::null;
  if (depth_av)
    {
    M = CreateMap(Noff0, Moff0,0,nloc0, mloc0, 1, nun);
    }
  else
    {
    M = CreateMap(Noff0, Moff0,Loff0,nloc0, mloc0, lloc0, nun);
    }
  return M;
  }

// public map creation function (can only create a limited range of maps)
Teuchos::RCP<Epetra_Map> Domain::CreateAssemblyMap(int nun, bool depth_av) const
  {
  HYMLS_PROF3(label_,"CreateAssemblyMap");
  Teuchos::RCP<Epetra_Map> M=Teuchos::null;
  if (depth_av)
    {
    M = CreateMap(Noff, Moff,0,nloc, mloc, 1, nun);
    }
  else
    {
    M = CreateMap(Noff, Moff,Loff,nloc, mloc, lloc, nun);
    }
  return M;
  }
  
// this version is private and very general
Teuchos::RCP<Epetra_Map> Domain::CreateMap(int noff_, int moff_, int loff_,
                              int nloc_, int mloc_, int lloc_, int nun) const
  {
  HYMLS_PROF3(label_,"CreateMap");
  int i0=noff_,j0=moff_,k0=loff_;  
  int i1=i0+nloc_-1,j1=j0+mloc_-1,k1=k0+lloc_-1;
  return CreateMap(i0,i1,j0,j1,k0,k1,0,n-1,0,m-1,0,l-1,nun,*comm);
  }
 
Teuchos::RCP<Epetra_Map>
Domain::CreateMap(int i0, int i1, int j0, int j1, int k0, int k1,
  int I0, int I1, int J0, int J1, int K0, int K1,
  int dof, const Epetra_Comm& comm) const
  {
  HYMLS_PROF3(label_, "CreateMap (2)");
  Teuchos::RCP<Epetra_Map> result = Teuchos::null;

  HYMLS_DEBUG("CreateMap ");
  HYMLS_DEBUG("[" << i0 << ".." << i1 << "]");
  HYMLS_DEBUG("[" << j0 << ".." << j1 << "]");
  HYMLS_DEBUG("[" << k0 << ".." << k1 << "]");

  int n = std::max(i1 - i0 + 1, 0); int N = I1 - I0 + 1;
  int m = std::max(j1 - j0 + 1, 0); int M = J1 - J0 + 1;
  int l = std::max(k1 - k0 + 1, 0); int L = K1 - K0 + 1;

  HYMLS_DEBVAR(N);
  HYMLS_DEBVAR(M);
  HYMLS_DEBVAR(L);

  int NumMyElements = n*m*l*dof;
  int NumGlobalElements = -1; // note that there may be overlap
  int *MyGlobalElements = new int[NumMyElements];

  int pos = 0;
  for (int k = k0; k <= k1; k++)
    for (int j = j0; j <= j1; j++)
      for (int i = i0; i <= i1; i++)
        for (int var = 0; var < dof; var++)
          {
          MyGlobalElements[pos++] =
            HYMLS::Tools::sub2ind(N, M, L, dof, i, j, k, var);
          }

  result = Teuchos::rcp(new Epetra_Map(NumGlobalElements,
      NumMyElements, MyGlobalElements, 0, comm));
  delete [] MyGlobalElements;
  return result;
  }

// create communication groups
Teuchos::RCP<Epetra_Comm> Domain::GetProcRow(int dim)
  {
  HYMLS_PROF3(label_,"GetProcRow");
#ifndef HAVE_MPI
  return comm; // can't split anything in sequential mode
#else //{
  Teuchos::RCP<Epetra_MpiComm> mpi_comm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(comm);
  if (mpi_comm==Teuchos::null) HYMLS::Tools::Error("Bad Communicator encountered!",__FILE__,__LINE__);
  MPI_Comm old_comm = mpi_comm->GetMpiComm();
  HYMLS_DEBVAR(old_comm);
  MPI_Comm row_comm;
  MPI_Group old_group, row_group;

    int nproc_row, disp, offset;
                    
    if (dim==0)
      {
      nproc_row=npN;
      disp=1;
      offset = (pidL*npM+pidM)*npN;
      }
    else if (dim==1)
      {
      nproc_row=npM;
      disp=npN;
      offset = pidL*npN*npM+pidN;
      }
    else if (dim==2)
      {
      nproc_row=npL;
      disp=npN*npM;
      offset = pidM*npN+pidN;
      }
    else
      {
      HYMLS::Tools::out()<<"Bad dimension to be extracted from proc array: "<<dim<<std::endl;
      HYMLS::Tools::Error("Cannot split dimension",__FILE__,__LINE__);
      }
    int *row_ranks = new int[nproc_row];


//  Create vector of ranks for row processes: 
    HYMLS_DEBUG("My position: ("<<pidN<<", "<<pidM<<", "<<pidL<<")");
    HYMLS_DEBUG("extract communicator for dim="<<dim);
    HYMLS_DEBUG("Creating sub-communicator consisting of:")
    for (int k = 0; k < nproc_row; k++)
      {
      row_ranks[k] = offset+k*disp;
      HYMLS_DEBUG(k<<": "<<row_ranks[k]);
      }

  /* ----------------------- */
  /* create the global group */
  /* ----------------------- */
  int ierr=MPI_Comm_group(old_comm, &old_group);

  if (ierr!=0) 
    {
    HYMLS::Tools::Error("MPI_Comm_group returned ierr="+Teuchos::toString(ierr),
        __FILE__,__LINE__);
    }

  /* --------------------- */
  /* Create the row group  */
  /* --------------------- */
  ierr=MPI_Group_incl(old_group, nproc_row, row_ranks, &row_group);
  if (ierr!=0) HYMLS::Tools::Error("MPI call 'Group_incl' failed!",__FILE__,__LINE__);

  /* ---------------------------------------------------- */
  /* Create the new communicator; MPI assigns the context */
  /* ---------------------------------------------------- */
  
  ierr = MPI_Comm_create(old_comm, row_group, &row_comm);
  if (ierr!=0) HYMLS::Tools::Error("MPI call 'Comm_create' failed!",__FILE__,__LINE__);
  
  // finally create the Epetra comm object
  Teuchos::RCP<Epetra_Comm> new_comm = Teuchos::rcp(new Epetra_MpiComm(row_comm));
  
  delete [] row_ranks;
  return new_comm;
#endif //}
  }

/////////////////////////////////////////////////////////////////////////////////
// Data Migration functions                                                    //
/////////////////////////////////////////////////////////////////////////////////

      //
      int Domain::Assembly2Standard
          (const Epetra_DistObject& source, Epetra_DistObject& target) const
        {
//        HYMLS_DEBUG("Enter Assembly2Standard");
#ifdef TESTING
        if (!(source.Map().SameAs(*AssemblyMap)&&target.Map().SameAs(*StandardMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
          }
#endif
        CHECK_ZERO(target.Export(source,*as2std,Zero));
//        HYMLS_DEBUG("Leave Assembly2Standard");
        return 0;
        }


      //
      int Domain::Standard2Assembly
        (const Epetra_DistObject& source, Epetra_DistObject& target) const
        {
//        HYMLS_DEBUG("Enter Standard2Assembly");
#ifdef TESTING
        if (!(source.Map().SameAs(*StandardMap)&&target.Map().SameAs(*AssemblyMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
          }
#endif
        CHECK_ZERO(target.Import(source,*as2std,Insert));
//        HYMLS_DEBUG("Leave Standard2Assembly");
        return 0;
        }

      //
      int Domain::Standard2Solve
        (const Epetra_DistObject& source, Epetra_DistObject& target) const
        {
//        HYMLS_DEBUG("Enter Standard2Solve");
#ifdef TESTING
        if (!(source.Map().SameAs(*StandardMap)&&target.Map().SameAs(*SolveMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
          }
#endif
        if (!UseDifferentSolveMap())
          {
          HYMLS::Tools::Error("not implemented",__FILE__,__LINE__);
          // this copy doesn't work for abstract objects, I put it in
          // when the function was only available for vectors (TODO)
//          target = source;
          }
        else
          {
          CHECK_ZERO(target.Export(source,*std2sol,Insert));
          }
//        HYMLS_DEBUG("Leave Standard2Solve");
        return 0;
        }

      //
      int Domain::Solve2Standard
        (const Epetra_DistObject& source, Epetra_DistObject& target) const
        {
#ifdef TESTING
        if (!(source.Map().SameAs(*SolveMap)&&target.Map().SameAs(*StandardMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
          }
#endif
        if (!UseDifferentSolveMap())
          {
          HYMLS::Tools::Error("not implemented",__FILE__,__LINE__);
          // this copy doesn't work for abstract objects, I put it in
          // when the function was only available for vectors (TODO)
//          target = source;
          }
        else
          {
          CHECK_ZERO(target.Import(source,*std2sol,Insert));
          }
        return 0;
        }

      //
      int Domain::Solve2Assembly
         (const Epetra_DistObject& source, Epetra_DistObject& target) const
        {
#ifdef TESTING
        if (!(source.Map().SameAs(*SolveMap)&&target.Map().SameAs(*AssemblyMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
          }
#endif
        if (!UseDifferentSolveMap())
          {
          Standard2Assembly(source,target);
          }
        else
          {
          Epetra_Vector tmp(*StandardMap);
          Solve2Standard(source,tmp);
          Standard2Assembly(tmp,target);
          }
        return 0;
        }


      //
      int Domain::Assembly2Solve
        (const Epetra_DistObject& source, Epetra_DistObject& target) const
        {
#ifdef TESTING
        if (!(source.Map().SameAs(*AssemblyMap)&&target.Map().SameAs(*SolveMap)))
          {
          HYMLS::Tools::Error("Invalid Transfer Function called!",__FILE__,__LINE__);
          }
#endif 
        if (!UseDifferentSolveMap())
          {
          Assembly2Standard(source,target);
          }
        else
          {
          Epetra_Vector tmp(*StandardMap);
          Assembly2Standard(source,tmp);
          Standard2Solve(tmp,target);
          }
        return 0;
        }

}//namespace

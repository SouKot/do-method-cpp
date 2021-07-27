#include "GaleriExt_Utils.H"
#include "Teuchos_ParameterList.hpp"

#include "Galeri_Maps.h"
#include "Galeri_Utils.h"

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

namespace GaleriExt {

Teuchos::RCP<Epetra_MultiVector>
Utils::CreateCartesianCoordinates(const std::string CoordType,
                           Teuchos::RCP<const Epetra_BlockMap> vertexMap,
                           Teuchos::ParameterList& List)
  {
  int dim=(CoordType=="1D")? 1: (CoordType=="2D")? 2: (CoordType=="3D")? 3: -1;
  int nx,ny,nz;
  Teuchos::RCP<Epetra_MultiVector> vertexCoords;

  // number of vertices
  nx=List.get("nx",16);
  ny=List.get("ny",dim>1? nx:1);
  nz=List.get("nz",dim>2?nx:1);
  
  int n_[3];
  n_[0]=nx; n_[1]=ny; n_[2]=nz;

  double stretchx=List.get("cell ratio (x)",1.0);
  double stretchy=List.get("cell ratio (y)",1.0);
  double stretchz=List.get("cell ratio (z)",1.0);
  
  double lx = List.get("lx",1.0);
  double ly = List.get("ly",lx);
  double lz = List.get("lz",lx);
  
  double L[3];
  L[0]=lx; L[1]=ly; L[2]=lz;

  Teuchos::ParameterList galeriList;  

  galeriList.set("nx",nx);
  galeriList.set("ny",ny);
  galeriList.set("nz",nz);

  std::string mapType="Cartesian3D";
  
  galeriList.set("lx",1.0);
  galeriList.set("ly",1.0);
  galeriList.set("lz",1.0);
  
  // create the coordinates of the cell centres of a uniform structured grid
  vertexCoords = Teuchos::rcp(Galeri::CreateCartesianCoordinates
        (CoordType, vertexMap.get(), galeriList));

  // grid stretching - transform x->tanh(alpha*x) such that max(dx)/min(dx)~=stretch
  
  // determine alphax/y/z 
  double alpha[3];
  for (int i=0;i<dim;i++) alpha[i]=0.0; // no stretching
  if (stretchx!=1.0) alpha[0] = FindStretchParam(nx,stretchx);
  if (stretchy!=1.0) alpha[1] = FindStretchParam(ny,stretchy);
  if (stretchz!=1.0) alpha[2] = FindStretchParam(nz,stretchz);

  for (int j=0;j<dim;j++)
    {
    for (int i=0; i<vertexCoords->MyLength();i++)
      {
      (*vertexCoords)[j][i] = stretch((*vertexCoords)[j][i], alpha[j]);
      }
    }
  double minval[dim];
  double maxval[dim];
  double scale[dim];
  
  CHECK_ZERO(vertexCoords->MinValue(minval));
  CHECK_ZERO(vertexCoords->MaxValue(maxval));
  for (int i=0;i<dim;i++) scale[i]=L[i]/(maxval[i]-minval[i]);

  // transform to cell-centred coordinates with cell faces at 0 and 1
  for (int j=0;j<dim;j++)
    {
    for (int i=0; i<vertexCoords->MyLength();i++)
      {
      (*vertexCoords)[j][i] = ((*vertexCoords)[j][i] - minval[j])*scale[j];
      }
    }

  return vertexCoords;
  }

/////////////////////////////////////////////////////////////////////////////////////////////

Teuchos::RCP<Epetra_MultiVector> Utils::CreateCellCenterCoordinates
        (const Epetra_MultiVector& vertexCoords,
         Teuchos::ParameterList& List)
  {

  // number of grid cells
  int nx=List.get("nx",16);
  int ny=List.get("ny",nx);
  int nz=List.get("nz",nx);

if (nx*ny*nz != vertexCoords.MyLength())
  {
  Error("mismatch in parameter file",__FILE__,__LINE__);
  }

  Teuchos::ParameterList galeriList;  

  galeriList.set("nx",nx-1);
  galeriList.set("ny",ny-1);
  galeriList.set("nz",nz-1);
  
  int dim=vertexCoords.NumVectors();

  std::string mapType="Cartesian"+Teuchos::toString(dim)+"D";

  Teuchos::RCP<Epetra_BlockMap> cellMap = Teuchos::null;
  const Epetra_Comm& const_comm = vertexCoords.Comm();
  Epetra_Comm& comm = const_cast<Epetra_Comm&>(const_comm);
  try {
    cellMap= Teuchos::rcp(Galeri::CreateMap(mapType, comm, galeriList));
    } catch (Galeri::Exception G) {G.Print();}

  Teuchos::RCP<Epetra_MultiVector> cellCoords=Teuchos::null;        

  // create the coordinates of the cell centres
  cellCoords = Teuchos::rcp(new Epetra_MultiVector(*cellMap,dim));
 
  if (dim!=3) Error("not implemented!",__FILE__,__LINE__);
  
  int ix,iy,iz;
  int idx1, idx2,jdx1,jdx2,kdx1,kdx2;
  for (int i = 0; i<cellCoords->MyLength(); i++) 
    {
    int gid_i = cellMap->GID(i);
    ind2sub(nx,ny,nz, gid_i,ix,iy,iz);
    idx1 = sub2ind(nx+1,ny+1,nz+1,ix,iy,iz);
    idx2 = sub2ind(nx+1,ny+1,nz+1,ix+1,iy,iz);
    jdx1 = sub2ind(nx+1,ny+1,nz+1,ix,iy,iz);
    jdx2 = sub2ind(nx+1,ny+1,nz+1,ix,iy+1,iz);
    kdx1 = sub2ind(nx+1,ny+1,nz+1,ix,iy,iz);
    kdx2 = sub2ind(nx+1,ny+1,nz+1,ix,iy,iz+1);

    (*cellCoords)[0][i] = (vertexCoords[0][idx1]+vertexCoords[0][idx1])*0.5;
    (*cellCoords)[1][i] = (vertexCoords[1][jdx1]+vertexCoords[1][jdx1])*0.5;
    (*cellCoords)[2][i] = (vertexCoords[2][kdx1]+vertexCoords[2][kdx1])*0.5;
    }
  return cellCoords;
  } 
  
  
double Utils::FindStretchParam(int n, double ratio)
  {
  double tol = ratio*1.0e-3;
  double fval = 1.0;
  int itmax = 100;
  double left=1.0e-3;
  double right=1.0e+3;
  double alpha=1.0;
  double dx0=1.0/(n+1);
  double x0=0, x1=dx0, y0=0.5, y1=0.5+dx0;
  // simple bisection method to find alpha such that
  // the max(dx)/min(dx) ratio matches the user input
  // for the transformed mesh tanh(alpha*x)
  for (int it=0;it<itmax;it++)
    {
    if (fval<tol) break;
    alpha = (left+right)*0.5;
    fval= (stretch(y1,alpha)-stretch(y0,alpha))/(stretch(x1,alpha)-stretch(x0,alpha))-ratio;
    if (fval>0) right = alpha;
    else left = alpha;     
    }
  return alpha;
  }

}//namespace

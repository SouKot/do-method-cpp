#ifndef GALERIEXT_UTILS_H
#define GALERIEXT_UTILS_H

#include "Teuchos_RCP.hpp"
#include "HYMLS_Tools.hpp"
#include "HYMLS_Macros.hpp"
class Epetra_BlockMap;
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_LinearProblem;
using namespace std;
namespace Teuchos {
  class ParameterList;
}



namespace GaleriExt {

class Utils
{

public:

//! similar to the Galeri function, except it accepts
//! a parameter "cell ratio (x)" (and similar for y and z)
//! to generate a mesh stretched such that the ratio of the
//! most central cells to the near boundary cells is approximately
//! the value indicated in the list.
static Teuchos::RCP<Epetra_MultiVector>
CreateCartesianCoordinates(const std::string CoordType,
                           Teuchos::RCP<const Epetra_BlockMap> vertexMap,
                           Teuchos::ParameterList& List);


//! given a vector of vertex coordinates, return a vector of coordinates
//! for the cell centers. The list should be the same used to construct 
//! the vertex coordinates
static Teuchos::RCP<Epetra_MultiVector> CreateCellCenterCoordinates
        (const Epetra_MultiVector& vertexCoords,
         Teuchos::ParameterList& List);

static inline void Error(std::string msg, const char* file, int line)
  {
  HYMLS::Tools::Error(msg,file,line);
  }

private:

// given x between 0 and 1, return x stretched towards 0 and 1 by
// the function tanh(alpha*x).
static inline double stretch(double x,double alpha)
  {
  if (alpha==0.0) return x;
  double pi = 2.0*std::asin(1.0);
  return (std::tanh(alpha*(x-0.5)*pi*2.0)+1.0)/2.0;
  }

// determine approximate alpha such that stretch(x,alpha)
// gives the desired stretched mesh for uniform x.
static double FindStretchParam(int n, double ratio);

//! \name indexing functions borrowed from HYMLS
//!@{

  //! converts linear index to cartesian subscripts
  static inline void ind2sub(int nx, int ny, int nz, int dof, 
        int idx, int& i, int& j, int& k, int& var)
    {
#ifdef TESTING    
    if ((idx<0)||(idx>=(nx*ny*nz*dof)))
      {
      std::cerr << "dim=["<<nx<<","<<ny<<","<<nz<<"], dof="<<dof<<": ind="<<idx<<std::endl;
      Error("ind2sub: Index out of range!",__FILE__,__LINE__);
      }
#endif      
    int rem=idx;
    var=MOD(rem,dof);
    rem=(rem-var)/dof;
    i=MOD(rem,nx);
    rem=(rem-i)/nx;
    j=MOD(rem,ny);
    rem=(rem-j)/ny;
    k=MOD(rem,nz);
    }
  //! converts linear index to cartesian subscripts
  static inline void ind2sub(int nx, int ny, int nz, int idx, int& i, int& j, int& k)
    {
    int dummy;
    ind2sub(nx,ny,nz,1,idx,i,j,k,dummy);
    return;
    }

  //! converts cartesian subscripts to linear index
  static inline int sub2ind(int nx, int ny, int nz, int dof, int i, int j, int k, int var)
    {
#ifdef TESTING    
    if ((i<0)||(i>=nx))
      {
      Error("sub2ind: i-Index out of range!",__FILE__,__LINE__);
      }
    if ((j<0)||(j>=ny))
      {
      Error("sub2ind: j-Index out of range!",__FILE__,__LINE__);
      }
    if ((k<0)||(k>=nz))
      {
      Error("sub2ind: k-Index out of range!",__FILE__,__LINE__);
      }
    if ((k<0)||(var>=dof))
      {
      Error("sub2ind: var-Index out of range!",__FILE__,__LINE__);
      }
#endif      
    return ((k*ny+j)*nx+i)*dof+var;
    }

  //! converts cartesian subscripts to linear index
  static inline int sub2ind(int nx, int ny, int nz, int i, int j, int k)
    {
    return sub2ind(nx,ny,nz,1,i,j,k,0);
    }

//!@}
                                                                                
};


}//namespace
#endif

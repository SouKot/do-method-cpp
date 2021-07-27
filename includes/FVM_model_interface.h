#ifndef FVM_MODEL_INTERFACE_H
#define FVM_MODEL_INTERFACE_H

extern "C" {

// initialize model, allocate memory etc.        
// input: nx etc: local grid size nx/y/z,        
//        material: array of material properties 
//              (size nx*ny*nz, eg. FLUID/SOLID) 
//        x/y/z: coordinates of grid vertices    
//              (size nx/y/z+1)                  
//        npar:  number of model parameters,     
//        p_values: initial values of parameters.
void model_init(int* nx, int* ny, int* nz, double* x_ptr, double* y_ptr, double* z_ptr, int* nx_, int* ny_, int* nz_,
                double* xminloc,double* xmaxloc, double* yminloc, double* ymaxloc,double* xmin_,double* xmax_, double* ymin_, double* ymax_, 
                int* material, int* nrows, double* sol, int* ierr );

// free all memory allocated by the fortran code
void model_free();

// get number of parameters
void model_get_num_params(int* npar);

// get parameter name i (i 1-based), null-terminated so it can be used in C)
void model_get_param_name(int* i, int* maxlen, char *name, 
        double* default_value, int* ierr);

// set model parameters
void model_setparams(int* npar, double* p_values, int* ierr);

// compute rhs = F(x_in)
// input: length of the vector, rhs (output), x_in (input)
void model_rhs(int *nrows, double *rhs, double *x_in, int* ierr);
void model_rhs_s(int *nrows, int* nzmax, double *svalues, int* srows, int* scols, int* ierr);
// compute stochastic forcing
  void model_stoch_frc(int *nrows, double *rhs, int* col,  int* ierr);
// compute the jacobian jac = dF/dx(x_in)
// input: nzmax - amount of memory allocated for values and cols.
//        x_in - current solution
// output jacobian for all nx*ny*nz local rows in compressed
//      row storage.
//        ierr = 0 if success
void model_jac(int* nrows, int* nzmax, double *values, int* rows, int* cols,
        double* x_in, int* ierr);
void model_jaclin(int* nrows, int* nzmax, 
         int* ierr, double* linvalues, int* linrows, int* lincols);

// get diagonal of mass matrix as a 1D array.
// input: nrows (just for checking): length of allocated storage
//        x
// output: mass: diagonal values of M
void model_massmat(int* nrows, double* x, double* mass);

// get number of parameters
void model_get_num_stochvectors(int* nWvec);

void model_bil(double* u, double* v, double* uv );
// get an estimate of the amount of ints and doubles allocated by
// the fortran code
void model_memory_estimate(double* num_ints, double* num_doubles);

}

#endif

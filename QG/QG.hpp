#ifndef QG_CPP_H
#define QG_CPP_H

#include "QG_Data.hpp"

namespace QG
{

class QG
{
    int m_;
    int n_;
    int ndim_;
    int adim_;

    double hdim_;     // depth first layer (m)
    double f0dim_;    // Coriolis parameter (s_1)
    double beta0dim_; // gradient plan. vort. (ms)-1
    double Lxdim_;    // length basin (m)
    double Lydim_;    // width basin (m)
    double udim_;     // Sverdrup velocity (m/s)
    double gdim_;     // density 1st layer (kg/m3)
    double rhodim_;
    double taudim_;   // windstress (N/m2)
    double Ahdim_;    // lateral friction (m2/s)
    double bfdim_;    // bottom friction (1/s)

    double xmin_;
    double xmax_;
    double ymin_;
    double ymax_;

    double dx_;
    double dy_;

    Vector1D x_;
    Vector1D y_;

    Vector2D tx_;
    Vector2D ty_;

    Vector1D par_;

    Vector3D Llzz_; Vector3D Llzp_;
    Vector3D Nlzz_; Vector3D Nlzp_;
    Vector3D Nlpp_; Vector3D Nlpz_;
    Vector3D Tlzz_; Vector3D Tlzp_;
    Vector3D Llpz_; Vector3D Llpp_;
    Vector3D Alzz_; Vector3D Alzp_;
    Vector3D Alpz_; Vector3D Alpp_;

    Matrix A_;
    Matrix B_;

    Matrix lu_;

    Vector1D diag_;
    Vector1D scaling_;

public:
    QG() = delete;
    QG(int m, int n);
    QG(QG const &other) = default;

    virtual ~QG();

    QG &operator=(QG const &other) = delete;

    void jacob(double const *un, double sig);
  void jacobian(double const *un, double sig, int *beg, int *jco, double *co );
    void rhs(double const *un, double *b);
    void bilin(double const *un,double const *vn, double *b);
    void mass(double *M);
    void readfort(int irs, double *u);

    int compute_precon();
    int solve(double *x);

    double get_par(int par);
    void set_par(int par, double val);

    void writeA(char const *name, double const *un, double sig);
    void writeM(char const *name);

    void apply(double const *x, double *y);

protected:
    double curl(double x, double y);
    void compute_forcing();
    void compute_linear();

    void assembleA();
    void assembleB();
    void fillcolA();
    void fillcolB();

    void nlin_jac(Vector2D const &om, Vector2D const &ps);
    void nlin_rhs(Vector2D const &om, Vector2D const &ps);

    void lin();
    void nonlin(int type, Vector3D &atom, Vector2D const &om, Vector2D const &ps);
    void deriv(int type, Vector3D &atom);

    void timedep();

    void boundaries();
    void u_to_psi(double const *un, Vector2D &om, Vector2D &ps);
    void pack(int *beg, int *jco, double *co);

    int lusolve(double *x);
    double dot(double *x, double *y);
    int cgstab(double *b, double tol);

    void compute_scaling();
    void apply_scaling();
    void apply_scaling(double *x);
    void remove_scaling();

    void Asort();
};

}

#endif

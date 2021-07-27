#ifndef QG_DATA_H
#define QG_DATA_H

namespace QG
{

class Matrix
{
    int ndim_;
    int adim_;
public:
    int *beg;
    int *jco;
    double *co;
    int *di;

    Matrix() = delete;
    Matrix(int ndim, int adim);
    Matrix(Matrix const &other);
    Matrix &operator=(Matrix const &other) = delete;
    ~Matrix();

    void pack();
};

class Vector1D
{
    int n_;
    double *data_;
public:
    Vector1D() = delete;
    Vector1D(int n);
    Vector1D(Vector1D const &other);
    ~Vector1D();

    double &operator[](int i);
    double const &operator[](int i) const;
    double &operator()(int i);
    double const &operator()(int i) const;

    double *get();
    int size();
};

class Vector2D
{
    int m_;
    int n_;
    double *data_;
public:
    Vector2D() = delete;
    Vector2D(int m, int n);
    Vector2D(Vector2D const &other);
    ~Vector2D();

    double &operator[](int i);
    double const &operator[](int i) const;
    double &operator()(int i, int j);
    double const &operator()(int i, int j) const;

    double *get();
    int size();
};

class Vector3D
{
    int l_;
    int m_;
    int n_;
    double *data_;
public:
    Vector3D() = delete;
    Vector3D(int l, int m, int n);
    Vector3D(Vector3D const &other);
    ~Vector3D();

    double &operator[](int i);
    double const &operator[](int i) const;
    double &operator()(int i, int j, int k);
    double const &operator()(int i, int j, int k) const;

    double *get();
    int size();
};

}

#endif

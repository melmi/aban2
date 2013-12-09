/**
 * pressure.h
 *
 * by Mohammad Elmi
 * Created on 23 Nov 2013
 */

#ifndef _CPS_H_
#define _CPS_H_

#include "domain.h"
#include <vector>
#include <algorithm>
#include <Eigen32/SparseCore>
#include <Eigen32/IterativeLinearSolvers>

namespace aban2
{

class CPS
{
public:
    CPS(domain *_d);

    ~CPS();

private:
    typedef Eigen::Triplet<double> Triplet;
    typedef std::vector<Triplet> TripletVector;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, long> SparseMatrix;
    typedef Eigen::BiCGSTAB<SparseMatrix> PressureSolver;
    typedef PressureSolver GradientSolver;

    domain *d;
    double h2inv;
    TripletVector *raw_coeffs;
    SparseMatrix *pmatrix, *gmatrix;
    PressureSolver *psolver;
    GradientSolver *gsolver;

    // functions to make pressure coeffs matrix
    void make_pressure_matrix();

    void add_pressure_row(mesh_row *row);

    void apply_pressure_row_bc(size_t ix0, size_t ix1, bctype bct);

    // functions to make rhs of pressure equation
    double *get_pressure_rhs();

    void apply_pressure_rhs_bc(double *rhs);

    void apply_single_pressure_rhs_bc(mesh_row *row, double *rhs, bcside side);

    // functions to make pressure gradient coeffs matrix
    void make_grad_matrix();

    void add_grad_row(mesh_row *row);

    void apply_grad_row_bc(size_t ix0 , size_t ix1, bctype bct);

public:
    void solve_pressure();

    void solve_grads();

    void update_u();
};

}

#endif
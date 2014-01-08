/**
 * projection.h
 *
 * by Mohammad Elmi
 * Created on 5 Oct 2013
 */

#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include <vector>
#include <algorithm>
#include <Eigen32/SparseCore>
#include <Eigen32/IterativeLinearSolvers>

#include "domain.h"

namespace aban2
{

class projection
{
public:
    projection(domain *_d);
    ~projection();

private:
    typedef Eigen::Triplet<double> triplet_t;
    typedef std::vector<triplet_t> triplet_vector;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, long> matrix_t;
    typedef Eigen::BiCGSTAB<matrix_t> solver_t;

    domain *d;
    double h2inv;
    matrix_t *pmatrix;
    solver_t *psolver;

    void make_matrix();
    void add_row(mesh_row *row, triplet_vector *coeffs);
    void apply_row_bc(size_t ix0, size_t ix1,bcdesc desc, triplet_vector *coeffs);

    double *get_rhs();
    void apply_rhs_bc(double *rhs);
    void apply_single_rhs_bc(mesh_row *row, double *rhs, bcside side);

public:
    void solve_p();
    void update_u();
    void update_uf();
};

}

#endif
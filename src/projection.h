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
#include <seldon-5.2/Seldon.hxx>
#include <seldon-5.2/SeldonSolver.hxx>

#include "domain.h"

#define SELDON_DEBUG_LEVEL_0

namespace aban2
{

class projection
{
public:
    projection(domain *_d);
    ~projection();

private:
    domain *d;
    double h2inv;
    Seldon::Matrix<double, Seldon::General, Seldon::ArrayRowSparse> pmatrix;

    void make_matrix();
    void add_row(mesh_row *row);
    void apply_row_bc(size_t ix0, size_t ix1, bcdesc desc);

    double *get_rhs();
    void apply_rhs_bc(double *rhs);

public:
    void solve_p();
    void update_u();
    void update_uf();
};

}

#endif
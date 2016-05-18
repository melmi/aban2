/**
 * projection.h
 *
 * by Mohammad Elmi
 * Created on 5 Oct 2013
 */

#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include "common.h"
#include <vector>
#include <algorithm>
#include "lpw.h"

#include "domain.h"

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

    lpw::matrix_t *create_matrix();
    void add_row(lpw::matrix_t *pmatrix, mesh::row *row);
    void apply_row_bc(lpw::matrix_t *pmatrix, size_t no0, size_t no1, bcdesc desc, double rho_b);

    double *get_rhs();
    void apply_rhs_bc(double *rhs);

public:
    void solve_p();
    void update_u();
    void update_uf();
};

}

#endif

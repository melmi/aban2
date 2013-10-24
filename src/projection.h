/**
 * projection.h
 *
 * by Mohammad Elmi
 * Created on 5 Oct 2013
 */

#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include <iostream>
#include <algorithm>
#include <eigen3/Eigen/IterativeLinearSolvers>

#include "domain.h"
#include "pressure.h"

using namespace std;

namespace aban2
{

class projection
{
    static double *get_div_ustar(domain *d)
    {
        double *result = (double *)d->create_var(1);
        gradient::divergance(d, d->ustar, result, &flow_boundary::velbc);
        double rho_dt = d->rho / d->dt;
        for (int i = 0; i < d->n; ++i) result[i] *= rho_dt;
        return result;
    }

    static void apply_dirichlet_p_bc(domain *d, double *rhs)
    {
        bcond bc;
        double coeff = 1. / d->delta / d->delta / 4.;

        for (size_t dir = 0; dir < NDIRS; ++dir)
            for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            {
                mesh_row *row = d->rows[dir] + irow;

                bc = d->boundaries[row->start_code].pbc;
                if (bc.type == bctype::dirichlet)
                {
                    rhs[d->cellno(row, 0)] += -6.*bc.val * coeff;
                    rhs[d->cellno(row, 1)] += -2.*bc.val * coeff;
                }
                bc = d->boundaries[row->end_code].pbc;
                if (bc.type == bctype::dirichlet)
                {
                    rhs[d->cellno(row, row->n - 1)] += -6.*bc.val * coeff;
                    rhs[d->cellno(row, row->n - 2)] += -2.*bc.val * coeff;
                }
            }
    }

    static double *get_pressure_rhs(domain *d)
    {
        double *rhs = get_div_ustar(d);
        apply_dirichlet_p_bc(d, rhs);
        return rhs;
    }

public:
    typedef Eigen::BiCGSTAB<pressure::sparse_matrix> psolver;

    static void solve_p(domain *d, psolver *solver)
    {
        double *rhs = get_pressure_rhs(d);

        Eigen::VectorXd b(d->n);
        for (int i = 0; i < d->n; ++i) b[i] = rhs[i];
        delete[] rhs;

        Eigen::VectorXd x(d->n);
        x = solver->solve(b);

        for (int i = 0; i < d->n; ++i) d->p[i] = x[i];
    }

    static void update_u(domain *d)
    {
        double **grad_p = (double **)d->create_var(2);
        gradient::of_scalar(d, d->p, grad_p, &flow_boundary::pbc);

        for (int i = 0; i < d->n; ++i)
            for (int dir = 0; dir < NDIRS; ++dir)
                d->u[dir][i] = d->ustar[dir][i] - grad_p[dir][i] * d->dt;

        d->delete_var(2, grad_p);
    }
};

}

#endif
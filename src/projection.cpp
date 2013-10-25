/**
 * projection.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "projection.h"

#include "gradient.h"

namespace aban2
{

double *projection::get_div_ustar(domain *d)
{
    double *result = (double *)d->create_var(1);
    gradient::divergance(d, d->ustar, result, &flow_boundary::velbc);
    double rho_dt = d->rho / d->dt;
    for (int i = 0; i < d->n; ++i) result[i] *= rho_dt;
    return result;
}

void projection::apply_dirichlet_p_bc(domain *d, double *rhs)
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

double *projection::get_pressure_rhs(domain *d)
{
    double *rhs = get_div_ustar(d);
    apply_dirichlet_p_bc(d, rhs);
    return rhs;
}

void projection::solve_p(domain *d, psolver *solver)
{
    double *rhs = get_pressure_rhs(d);

    Eigen::VectorXd b(d->n);
    for (int i = 0; i < d->n; ++i) b[i] = rhs[i];
    delete[] rhs;

    Eigen::VectorXd x(d->n);
    x = solver->solve(b);

    for (int i = 0; i < d->n; ++i) d->p[i] = x[i];
}

void projection::update_u(domain *d)
{
    double **grad_p = (double **)d->create_var(2);
    gradient::of_scalar(d, d->p, grad_p, &flow_boundary::pbc);

    for (int i = 0; i < d->n; ++i)
        for (int dir = 0; dir < NDIRS; ++dir)
            d->u[dir][i] = d->ustar[dir][i] - grad_p[dir][i] * d->dt;

    d->delete_var(2, grad_p);
}

}

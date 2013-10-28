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
    gradient::divergance(d, d->ustar, result, &bcondition::u);
    double rho_dt = d->rho / d->dt;
    for (int i = 0; i < d->n; ++i) result[i] *= rho_dt;
    return result;
}

void apply_single_p_bc(domain *d, mesh_row *row, double *rhs, bcside side, double coeff)
{
    size_t idx0, idx1;
    bcondition *bc;
    double dx;

    if (side == bcside::start)
    {
        bc = d->boundaries[row->start_code];
        idx0 = d->cellno(row, 0);
        idx1 = d->cellno(row, 1);
        dx = -d->delta / 2.0;
    }
    else
    {
        bc = d->boundaries[row->end_code];
        idx0 = d->cellno(row, row->n - 1);
        idx1 = d->cellno(row, row->n - 2);
        dx = +d->delta / 2.0;
    }

    if (bc->ptype == bctype::dirichlet)
    {
        double p = bc->p(d, row, bcside::start, row->dir);
        rhs[idx0] -= 6.0 * p * coeff;
        rhs[idx1] -= 2.0 * p * coeff;
    }
    else
    {
        double rhog = d->rho * d->g.components[row->dir] * dx;
        rhs[idx0] -= 6.0 * rhog * coeff;
        rhs[idx1] -= 2.0 * rhog * coeff;
    }
}

void projection::apply_p_bc(domain *d, double *rhs)
{
    double coeff = 1.0 / d->delta / d->delta / 4.0;

    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;

            apply_single_p_bc(d, row, rhs, bcside::start, coeff);
            apply_single_p_bc(d, row, rhs, bcside::end, coeff);
        }
}

double *projection::get_pressure_rhs(domain *d)
{
    double *rhs = get_div_ustar(d);
    std::copy_n(rhs, d->n, d->u[0]);
    apply_p_bc(d, rhs);
    std::copy_n(rhs, d->n, d->u[1]);
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
    gradient::of_scalar(d, d->p, grad_p, &bcondition::p);

    for (int i = 0; i < d->n; ++i)
        for (int dir = 0; dir < NDIRS; ++dir)
            d->u[dir][i] = d->ustar[dir][i] - grad_p[dir][i] * d->dt;

    d->delete_var(2, grad_p);
}

}

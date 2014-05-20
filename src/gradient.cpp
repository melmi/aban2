/**
 * gradient.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "gradient.h"

namespace aban2
{

void gradient::add_1d_row(domain *d, mesh_row *row, double *phi, double *grad_dir, flowbc::member bc)
{
    double *phi_row = d->extract_scalars(row, phi);
    double *grad_row = d->extract_scalars(row, grad_dir);
    double twodx = d->delta * 2;
    size_t n = row->n;

    for (size_t i = 1; i < n - 1; ++i)
        grad_row[i] += (phi_row[i + 1] - phi_row[i - 1]) / twodx;

    auto startbc = flowbc::fval_start(d, row, phi, bc);
    auto endbc   = flowbc::fval_end  (d, row, phi, bc);

    grad_row[0] += (0.5 * (phi_row[0] + phi_row[1]) - startbc) / d->delta;
    grad_row[n - 1] += (endbc - 0.5 * (phi_row[n - 1] + phi_row[n - 2])) / d->delta;

    d->insert_scalars(row, grad_dir, grad_row);
    delete[] phi_row;
    delete[] grad_row;
}

double *gradient::of_scalar_dir_oriented(domain *d, double *phi, flowbc::member bc, size_t dir)
{
    double *result = (double *)d->create_var(1);
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            add_1d_row(d, d->rows[dir] + irow, phi, result, bc);
    return result;
}


double **gradient::of_scalar(domain *d, double *phi, flowbc::member bc)
{
    double **result = new double*[3];
    for (size_t dir = 0; dir < NDIRS; ++dir)
        result[dir] = of_scalar_dir_oriented(d, phi, bc, dir);
    return result;
}

double ** *gradient::of_vec(domain *d, double **phi, flowbc::member bc[3])
{
    double ***result = new double **[3];
    for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
        result[icmpnt] = of_scalar(d, phi[icmpnt], bc[icmpnt]);
    return result;
}

double *gradient::divergance(domain *d, double **phi, flowbc::member bc[3])
{
    double *result = (double *)d->create_var(1);
    for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
        for (size_t irow = 0; irow < d->nrows[icmpnt]; ++irow)
            add_1d_row(d, d->rows[icmpnt] + irow, phi[icmpnt], result, bc[icmpnt]);
    return result;
}

}


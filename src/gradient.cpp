/**
 * gradient.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "gradient.h"

namespace aban2
{

void gradient::add_1d_row(domain *d, mesh_row *row, double *phi, double *grad, bcondition::func bcfunc, size_t cmpnt)
{
    double dx = d->delta;
    size_t n = row->n;

    for (size_t i = 1; i < n - 1; ++i)
        grad[i] += (phi[i + 1] - phi[i - 1]) / (2.0 * dx);

    auto startbc = d->boundaries[row->start_code];
    auto endbc = d->boundaries[row->end_code];

    grad[0] +=
        (
            0.5 * (phi[0] + phi[1]) -
            startbc->row_face_val(bcfunc, phi, row, bcside::start, cmpnt)
        ) / dx;
    grad[n - 1] +=
        (
            endbc  ->row_face_val(bcfunc, phi, row, bcside::end  , cmpnt) -
            0.5 * (phi[n - 1] + phi[n - 2])
        ) / dx;
}

double *gradient::get_1d_row(domain *d, mesh_row *row, double *phi, bcondition::func bcfunc, size_t cmpnt)
{
    double *grad = new double[row->n];
    std::fill_n(grad, row->n, 0.0);
    add_1d_row(d, row, phi, grad, bcfunc, cmpnt);
    return grad;
}

void gradient::of_scalar(domain *d, double *phi, double **grad, bcondition::func bcfunc)
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            double *vals = d->extract_scalars(row, phi);
            double *g = get_1d_row(d, row, vals, bcfunc, dir);
            d->insert_scalars(row, grad[dir], g);
            delete[] g;
            delete[] vals;
        }
}

void gradient::of_vec(domain *d, double **phi, double ***grad, bcondition::func bcfunc)
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
            {
                double *cmpnt = d->extract_scalars(row, phi[icmpnt]);
                double *g = get_1d_row(d, row, cmpnt, bcfunc, icmpnt);
                d->insert_scalars(row, grad[icmpnt][dir], g);
                delete[] g;
                delete[] cmpnt;
            }
        }
}

void gradient::divergance(domain *d, double **phi, double *divergance, bcondition::func bcfunc)
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            double *vals = d->extract_scalars(row, phi[dir]);
            double *g = d->extract_scalars(row, divergance);
            add_1d_row(d, row, vals, g, bcfunc, dir);
            d->insert_scalars(row, divergance, g);
            delete[] vals;
        }
}

}


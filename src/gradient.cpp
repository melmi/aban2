/**
 * gradient.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "gradient.h"

namespace aban2
{

void gradient::add_1d_row(domain_t *d, row_t *row, double *phi, double *grad_dir, flowbc::bc_val_getter bc)
{
    double *phi_row = d->extract_scalars(row, phi);
    double *grad_row = d->extract_scalars(row, grad_dir);
    double twodx = d->delta * 2;
    size_t n = row->n;

    for (size_t i = 1; i < n - 1; ++i)
        grad_row[i] += (phi_row[i + 1] - phi_row[i - 1]) / twodx;

    auto startbc = bc(d, row, phi, bcside::start);
    auto endbc = bc(d, row, phi, bcside::end);

    grad_row[0] += (0.5 * (phi_row[0] + phi_row[1]) - startbc) / d->delta;
    grad_row[n - 1] += (endbc - 0.5 * (phi_row[n - 1] + phi_row[n - 2])) / d->delta;

    d->insert_scalars(row, grad_dir, grad_row);
    delete[] phi_row;
    delete[] grad_row;
}

double *gradient::of_scalar_dir(domain_t *d, double *phi, flowbc::bc_val_getter bc, size_t dir)
{
    double *result = (double *)d->create_var(1);
    for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        add_1d_row(d, d->rows[dir] + irow, phi, result, bc);
    return result;
}

double **gradient::of_scalar(domain_t *d, double *phi, flowbc::bc_val_getter bc)
{
    double **result = new double *[3];
    result[2] = nullptr;
    for (size_t dir = 0; dir < NDIRS; ++dir)
        result[dir] = of_scalar_dir(d, phi, bc, dir);
    return result;
}

double ***gradient::of_vec(domain_t *d, double **phi, flowbc::bc_val_getter bc[3])
{
    double ***result = new double **[3];
    result[2] = nullptr;
    for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
        result[icmpnt] = of_scalar(d, phi[icmpnt], bc[icmpnt]);
    return result;
}

double *avg_uf(domain_t *d, size_t dir)
{
    double *result = new double[d->n];

    for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
    {
        auto row = d->rows[dir] + irow;
        auto uf_row = d->extract_scalars(row, d->uf[dir]);
        auto uf_row_avg = new double[row->n];

        size_t n = row->n;
        for (size_t i = 1; i < n - 1; ++i)
            uf_row_avg[i] = (uf_row[i - 1] + uf_row[i]) / 2.0;
        uf_row_avg[0] = (flowbc::bc_u_getter[dir](d, row, d->u[dir], bcside::start) + uf_row[0]) / 2.0;
        uf_row_avg[n - 1] = (flowbc::bc_u_getter[dir](d, row, d->u[dir], bcside::end) + uf_row[n - 2]) / 2.0;

        d->insert_scalars(row, result, uf_row_avg);
        delete[] uf_row;
        delete[] uf_row_avg;
    }

    return result;
}

double *gradient::of_uf_dir(domain_t *d, size_t dir)
{
    double *result = new double[d->n];

    for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
    {
        auto row = d->rows[dir] + irow;
        auto uf_row = d->extract_scalars(row, d->uf[dir]);
        auto grad_row = new double[row->n];

        size_t n = row->n;
        for (size_t i = 1; i < n - 1; ++i)
            grad_row[i] = (uf_row[i] - uf_row[i - 1]) / d->delta;
        grad_row[0] = (uf_row[0] - flowbc::bc_u_getter[dir](d, row, d->u[dir], bcside::start)) / d->delta;
        grad_row[n - 1] = (flowbc::bc_u_getter[dir](d, row, d->u[dir], bcside::end) - uf_row[n - 2]) / d->delta;

        d->insert_scalars(row, result, grad_row);
        delete[] uf_row;
        delete[] grad_row;
    }

    return result;
}

double ***gradient::of_uf(domain_t *d)
{
    double *avg[3]{
        avg_uf(d, 0),
        avg_uf(d, 1),
        avg_uf(d, 2)};

    double ***result = new double **[3];
    result[2] = nullptr;
    for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
    {
        result[icmpnt] = of_scalar(d, avg[icmpnt], flowbc::bc_u_getter[icmpnt]);
        delete[] result[icmpnt][icmpnt];
        result[icmpnt][icmpnt] = of_uf_dir(d, icmpnt);
    }

    delete[] avg[0];
    delete[] avg[1];
    delete[] avg[2];

    return result;
}

double *gradient::divergance(domain_t *d)
{
    double *result = (double *)d->create_var(1);
    for (size_t dir = 0; dir < NDIRS; ++dir)
    {
        auto g_dir = of_uf_dir(d, dir);
        for (size_t i = 0; i < d->n; ++i)
            result[i] += g_dir[i];
        delete[] g_dir;
    }
    return result;
}

double *gradient::divergance_of(domain_t *d, double **phi, flowbc::bc_val_getter bc[3])
{
    double *result = (double *)d->create_var(1);
    for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
        for (size_t irow = 0; irow < d->nrows[icmpnt]; ++irow)
            add_1d_row(d, d->rows[icmpnt] + irow, phi[icmpnt], result, bc[icmpnt]);
    return result;
}
}

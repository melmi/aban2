/**
 * advection.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "advection.h"

#include <algorithm>

#include "domain.h"
#include "gradient.h"

namespace aban2
{

void advection::advect(mesh_row *row, double *phi, double *grad_dir, flowbc::member bc)
{
    double flux;
    double dx = d->delta, dt = d->dt;
    size_t n = row->n;

    double *u = d->extract_scalars(row, d->uf[row->dir]);
    double *phi_row = d->extract_scalars(row, phi);
    double *grad_row = d->extract_scalars(row, grad_dir);
    double *mass = new double[n];

    for (size_t i = 0; i < n; ++i) mass[i] = phi_row[i] * dx;

    for (size_t face = 0; face < n - 1; ++face)
    {
        size_t i = face, ii = face + 1;
        double uf = u[i];
        if (uf > 0)
            flux = (phi_row[i ] + (+dx / 2. - uf * dt / 2.) * grad_row[i ]) * uf * dt;
        else
            flux = (phi_row[ii] + (-dx / 2. - uf * dt / 2.) * grad_row[ii]) * uf * dt;
        mass[i ] -= flux;
        mass[ii] += flux;
    }

    double uf_start = flowbc::fval_start(d, row, d->u[row->dir], flowbc::umembers[row->dir]);
    double uf_end   = flowbc::fval_end  (d, row, d->u[row->dir], flowbc::umembers[row->dir]);

    mass[0    ] += flowbc::fval_start(d, row, phi, bc) * uf_start * dt;
    mass[n - 1] -= flowbc::fval_end  (d, row, phi, bc) * uf_end   * dt;

    for (size_t i = 0; i < n; ++i) phi_row[i] = mass[i] / dx;

    d->insert_scalars(row, phi, phi_row);

    delete[] phi_row;
    delete[] grad_row;
    delete[] mass;
    delete[] u;
}

void advection::advect_ustar()
{
    for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
        for (size_t dir = 0; dir < NDIRS; ++dir)
        {
            double *grad_dir = gradient::of_scalar_dir_oriented(
                                   d,
                                   d->ustar[icmpnt],
                                   flowbc::umembers[icmpnt],
                                   dir);
            for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            {
                mesh_row *row = d->rows[dir] + irow;
                advect(row, d->ustar[icmpnt], grad_dir, flowbc::umembers[icmpnt]);
            }
            delete[] grad_dir;
        }
}

}

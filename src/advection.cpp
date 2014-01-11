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

void advection::advect(mesh_row *row, double *var, bcondition::func bcfunc, size_t cmpnt)
{
    double flux;
    double dx = d->delta, dt = d->dt;
    size_t n = row->n;

    double *u = d->extract_scalars(row, d->uf[row->dir]);
    double *phi = d->extract_scalars(row, var);
    double *grad = gradient::get_1d_row(d, row, phi, bcfunc, cmpnt);
    double *mass = new double[n];

    for (size_t i = 0; i < n; ++i) mass[i] = phi[i] * dx;

    for (size_t face = 0; face < n - 1; ++face)
    {
        size_t i = face, ii = face + 1;
        double uf = u[i];
        if (uf > 0)
            flux = (phi[i ] + (+dx / 2. - uf * dt / 2.) * grad[i ]) * uf * dt;
        else
            flux = (phi[ii] + (-dx / 2. - uf * dt / 2.) * grad[ii]) * uf * dt;
        mass[i ] -= flux;
        mass[ii] += flux;
    }

    auto startbc = d->boundaries[row->start_code];
    auto endbc   = d->boundaries[row->end_code  ];

    double uf_start = startbc->face_val(&bcondition::u, d->u[row->dir], row, bcside::start, row->dir);
    double uf_end   = endbc  ->face_val(&bcondition::u, d->u[row->dir], row, bcside::end  , row->dir);

    mass[0    ] += startbc->face_val(bcfunc, var, row, bcside::start, cmpnt) * uf_start * dt;
    mass[n - 1] -= endbc  ->face_val(bcfunc, var, row, bcside::end  , cmpnt) * uf_end   * dt;

    for (size_t i = 0; i < n; ++i) phi[i] = mass[i] / dx;

    d->insert_scalars(row, var, phi);

    delete[] phi;
    delete[] grad;
    delete[] mass;
    delete[] u;
}

void advection::advect_ustar()
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
                advect(row, d->ustar[icmpnt], &bcondition::u, icmpnt);
        }
}

}

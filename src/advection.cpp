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

void advection::advect(domain *d, mesh_row *row, double *phi, double *u, bcondition::func bcfunc, size_t cmpnt)
{
    double flux;

    double dx = d->delta, dt = d->dt;
    size_t n = row->n;

    double *grad = gradient::get_1d_row(d, row, phi, bcfunc, cmpnt);
    double *mass = new double[n];

    for (size_t i = 0; i < n; ++i) mass[i] = phi[i] * dx;

    for (size_t face = 0; face < n - 1; ++face)
    {
        size_t i = face, ii = face + 1;
        double uf = 0.5 * (u[i] + u[ii]);
        if (uf > 0)
            flux = (phi[i] + (dx / 2. - uf * dt / 2.) * grad[i]) * uf * dt;
        else
            flux = (phi[ii] + (-dx / 2. - uf * dt / 2.) * grad[ii]) * uf * dt;
        mass[i] -= flux;
        mass[ii] += flux;
    }

    auto startbc = d->boundaries[row->start_code];
    auto endbc   = d->boundaries[row->end_code  ];

    mass[0    ] += (startbc->*bcfunc)(d, row, bcside::start, cmpnt) * u[0    ] * dt;
    mass[n - 1] -= (startbc->*bcfunc)(d, row, bcside::end  , cmpnt) * u[n - 1] * dt;

    for (size_t i = 0; i < n; ++i) phi[i] = mass[i] / dx;

    delete[] grad;
    delete[] mass;
}

void advection::advect_ustar(domain *d)
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            double *u = d->extract_scalars(row, d->u[dir]);
            for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
            {
                double *cmpnt = d->extract_scalars(row, d->ustar[icmpnt]);
                advect(d, row, cmpnt, u, &bcondition::u, dir);
                d->insert_scalars(row, d->ustar[icmpnt], cmpnt);
                delete[] cmpnt;
            }
            delete[] u;
        }
}

}

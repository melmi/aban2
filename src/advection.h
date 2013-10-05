/**
 * advection.h
 *
 * by Mohammad Elmi
 * Created on 4 Sep 2013
 */

#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include "domain.h"
#include "gradient.h"
#include <algorithm>

namespace aban2
{

class advection
{
public:
    static void advect(size_t n, double *phi, double *u, double dt, double dx, bcond startbc, bcond endbc)
    {
        double flux;

        double *grad = gradient::get_1d_row(n, phi, dx, startbc, endbc);
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

        if (startbc.type == bctype::dirichlet)
            mass[0] += startbc.val * u[0] * dt;
        else
            mass[0] += phi[0] * u[0] * dt;

        if (endbc.type == bctype::dirichlet)
            mass[n - 1] -= endbc.val * u[0] * dt;
        else
            mass[n - 1] -= phi[n - 1] * u[n - 1] * dt;

        for (size_t i = 0; i < n; ++i) phi[i] = mass[i] / dx;

        delete[] grad;
        delete[] mass;
    }

    static void advect_ustar(domain *d)
    {
        for (size_t dir = 0; dir < NDIRS; ++dir)
            for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            {
                mesh_row *row = d->rows[dir] + irow;
                double *u = d->extract_scalars(row, d->u[dir]);
                for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
                {
                    double *cmpnt = d->extract_scalars(row, d->ustar[icmpnt]);
                    bcond *start_bc = d->boundaries[row->start_bc].velbc + icmpnt;
                    bcond *end_bc = d->boundaries[row->end_bc].velbc + icmpnt;
                    advect(row->n, cmpnt, u, d->dt, d->delta, *start_bc, *end_bc);
                    d->insert_scalars(row, d->ustar[icmpnt], cmpnt);
                    delete[] cmpnt;
                }
                delete[] u;
            }
    }
};

}

#endif
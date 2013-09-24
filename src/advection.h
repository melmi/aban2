/**
 * advection.h
 *
 * by Mohammad Elmi
 * Created on 4 Sep 2013
 */

#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include "domain.h"
#include <algorithm>

namespace aban2
{

class advection
{
public:
    static double *get_grad(size_t n, double *phi, double dx)
    {
        double *grad = new double[n];
        for (size_t i = 1; i < n - 1; ++i)
            grad[i] = (phi[i + 1] - phi[i - 1]) / (2.*dx);

        grad[0] = (phi[1] - phi[0]) / dx;
        grad[n - 1] = (phi[n - 1] - phi[n - 2]) / dx;

        return grad;
    }

    static void advect(size_t n, double *phi, double *u, double dt, double dx, bcond startbc, bcond endbc)
    {
        double flux;
        double *grad = get_grad(n, phi, dx);
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

        if (startbc.type == bctype::neumann)
            mass[0] += startbc.val * dt;
        else
            mass[0] += u[0] * dt * phi[0];

        if (endbc.type == bctype::neumann)
            mass[n - 1] += endbc.val * dt;
        else
            mass[n - 1] += u[n - 1] * dt * phi[n - 1];

        for (size_t i = 0; i < n; ++i) phi[i] = mass[i] / dx;

        delete[] grad;
        delete[] mass;
    }

    static void advect_ustar(domain *d)
    {

        for (size_t ieq = 0; ieq < 3; ++ieq)
            for (size_t irow = 0; irow < d->nrows[ieq]; ++irow)
            {
                mesh_row *row = d->rows[ieq] + irow;
                double *u = d->extract_scalars(row, d->u[ieq]);
                for (size_t icmpnt = 0; icmpnt < 3; ++icmpnt)
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
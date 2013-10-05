/**
 * gradient.h
 *
 * by Mohammad Elmi
 * Created on 5 Oct 2013
 */

#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include "domain.h"
#include <algorithm>

namespace aban2
{

class gradient
{
public:
    static void add_1d_row(size_t n, double *phi, double *grad, double dx, bcond startbc, bcond endbc)
    {
        for (size_t i = 1; i < n - 1; ++i)
            grad[i] += (phi[i + 1] - phi[i - 1]) / (2.*dx);

        double phi0 = startbc.type == bctype::dirichlet ? startbc.val : phi[0];
        double phin = endbc.type == bctype::dirichlet ? endbc.val : phi[n - 1];

        grad[0] += (0.5 * (phi[0] + phi[1]) - phi0) / dx;
        grad[n - 1] += (phin - 0.5 * (phi[n - 1] + phi[n - 2])) / dx;
    }

    static double *get_1d_row(size_t n, double *phi, double dx, bcond startbc, bcond endbc)
    {
        double *grad = new double[n];
        std::fill_n(grad, n, 0.0);
        add_1d_row(n, phi, grad, dx, startbc, endbc);
        return grad;
    }

    static void of_scalar(domain *d, double *phi, double **grad, bcond flow_boundary::*bc)
    {
        for (size_t dir = 0; dir < NDIRS; ++dir)
            for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            {
                mesh_row *row = d->rows[dir] + irow;
                double *vals = d->extract_scalars(row, phi);
                bcond start_bc = d->boundaries[row->start_bc].*bc;
                bcond end_bc = d->boundaries[row->end_bc].*bc;
                double *g = get_1d_row(row->n, vals, d->delta, start_bc, end_bc);
                d->insert_scalars(row, grad[dir], g);
                delete[] g;
                delete[] vals;
            }
    }

    static void of_vec(domain *d, double **phi, double ***grad, bcond * flow_boundary::*bc)
    {
        for (size_t dir = 0; dir < NDIRS; ++dir)
            for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            {
                mesh_row *row = d->rows[dir] + irow;
                for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
                {
                    double *cmpnt = d->extract_scalars(row, phi[icmpnt]);
                    bcond *start_bc = d->boundaries[row->start_bc].*bc + icmpnt;
                    bcond *end_bc = d->boundaries[row->end_bc].*bc + icmpnt;
                    double *g = get_1d_row(row->n, cmpnt, d->delta, *start_bc, *end_bc);
                    d->insert_scalars(row, grad[icmpnt][dir], g);
                    delete[] g;
                    delete[] cmpnt;
                }
            }
    }

    static void divergance(domain *d, double **phi, double *divergance, bcond (flow_boundary::*bc)[3])
    {
        for (size_t dir = 0; dir < NDIRS; ++dir)
            for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            {
                mesh_row *row = d->rows[dir] + irow;
                double *vals = d->extract_scalars(row, phi[dir]);
                bcond *start_bc = d->boundaries[row->start_bc].*bc + dir;
                bcond *end_bc = d->boundaries[row->end_bc].*bc + dir;
                double *g = d->extract_scalars(row, divergance);
                add_1d_row(row->n, vals, g, d->delta, *start_bc, *end_bc);
                d->insert_scalars(row, divergance, g);
                delete[] vals;
            }
    }
};

}

#endif

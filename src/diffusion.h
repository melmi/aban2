/**
 * diffusion.h
 *
 * by Mohammad Elmi
 * Created on 5 Sep 2013
 */

#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include "domain.h"
#include <algorithm>

namespace aban2
{

class diffusion
{
public:
    static void solve_tridiagonal_in_place_destructive(double *x, const size_t N, const double *a, const double *b, double *c)
    {
        /* unsigned integer of same size as pointer */
        size_t in;

        /*
         solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
         note that contents of input vector c will be modified, making this a one-time-use function
         x[] - initially contains the input vector v, and returns the solution x. indexed from [0, ..., N - 1]
         N - number of equations
         a[] - subdiagonal (means it is the diagonal below the main diagonal) -- indexed from [1, ..., N - 1]
         b[] - the main diagonal, indexed from [0, ..., N - 1]
         c[] - superdiagonal (means it is the diagonal above the main diagonal) -- indexed from [0, ..., N - 2]
         */

        c[0] = c[0] / b[0];
        x[0] = x[0] / b[0];

        /* loop from 1 to N - 1 inclusive */
        for (in = 1; in < N; in++)
        {
            double m = 1.0 / (b[in] - a[in] * c[in - 1]);
            c[in] = c[in] * m;
            x[in] = (x[in] - a[in] * x[in - 1]) * m;
        }

        /* loop from N - 2 to 0 inclusive, safely testing loop end condition */
        for (in = N - 1; in-- > 0; )
            x[in] = x[in] - c[in] * x[in + 1];
    }

    static void diffuse(size_t n, double *phi, double d, double dt, double dx, bcond startbc, bcond endbc)
    {
        double coeff = d * dt / dx / dx;

        double *aa = new double[n], *bb = new double[n], *cc = new double[n];

        for (int i = 0; i < n; ++i)
        {
            aa[i] = cc[i] = -coeff;
            bb[i] = 1. + 2.*coeff;
        }

        // TODO: set boundary conditions here
        bb[0] = bb[n - 1] = 1. + coeff;
        // aa[0] = cc[n - 1] = 0;

        solve_tridiagonal_in_place_destructive(phi, n, aa, bb, cc);

        delete aa;
        delete bb;
        delete cc;
    }

    static void diffuse_ustar(domain d)
    {

        for (size_t ieq = 0; ieq < 3; ++ieq)
            for (size_t irow = 0; irow < d.nrows[ieq]; ++irow)
            {
                mesh_row *row = d.rows[ieq] + irow;
                for (size_t icmpnt = 0; icmpnt < 3; ++icmpnt)
                {
                    double *cmpnt = d.extract_scalars(*row, d.ustar[icmpnt]);
                    bcond *start_bc = d.boundaries[row->start_bc].velbc + icmpnt;
                    bcond *end_bc = d.boundaries[row->end_bc].velbc + icmpnt;
                    diffuse(row->n, cmpnt, d.mu, d.dt, d.delta, *start_bc, *end_bc);
                    d.insert_scalars(*row, d.ustar[icmpnt], cmpnt);
                    delete[] cmpnt;
                }
            }
    }
};

}

#endif
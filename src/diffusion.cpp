/**
 * diffusion.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "diffusion.h"

namespace aban2
{
void diffusion::solve_tridiagonal_in_place_destructive(double *x, const size_t N, const double *a, const double *b, double *c)
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

void diffusion::diffuse(size_t n, double *phi, double d, double dt, double dx, bcond startbc, bcond endbc)
{
    double coeff = d * dt / dx / dx;

    double *aa = new double[n], *bb = new double[n], *cc = new double[n];

    for (size_t i = 0; i < n; ++i)
    {
        aa[i] = cc[i] = -coeff;
        bb[i] = 1.0 + 2.0 * coeff;
    }

    if (startbc.type == bctype::dirichlet)
    {
        bb[0] = 1.0 + 3.0 * coeff;
        phi[0] += 2.0 * coeff *  startbc.val;
    }
    else
    {
        bb[0] = 1.0 + coeff;
    }

    if (endbc.type == bctype::dirichlet)
    {
        bb[n - 1] = 1.0 + 3.0 * coeff;
        phi[n - 1] +=  2.0 * coeff * endbc.val;
    }
    else
    {
        bb[n - 1] = 1.0 + coeff;
    }
    solve_tridiagonal_in_place_destructive(phi, n, aa, bb, cc);

    delete aa;
    delete bb;
    delete cc;
}

void diffusion::diffuse_ustar(domain *d)
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
            {
                double *cmpnt = d->extract_scalars(row, d->ustar[icmpnt]);
                bcond *start_code = d->boundaries[row->start_code].velbc + icmpnt;
                bcond *end_code = d->boundaries[row->end_code].velbc + icmpnt;
                diffuse(row->n, cmpnt, d->mu, d->dt, d->delta, *start_code, *end_code);
                d->insert_scalars(row, d->ustar[icmpnt], cmpnt);
                delete[] cmpnt;
            }
        }
}

}


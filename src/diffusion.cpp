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

void diffusion::diffuse(mesh_row *row, double *phi, double *D, flowbc::member mem)
{
    double *phi_row = d->extract_scalars(row, phi);
    double *d_row = d->extract_scalars(row, D);
    double dx = d->delta, dt = d->dt;
    size_t n = row->n;
    double dt_dx2 = dt / dx / dx;

    double *aa = new double[n], *bb = new double[n], *cc = new double[n];

    for (size_t i = 0; i < n; ++i)
    {
        aa[i] = cc[i] = -d_row[i] * dt_dx2;
        bb[i] = 1.0 + 2.0 * d_row[i] * dt_dx2;
    }

    auto startbc = (d->boundaries[row->start_code]->*mem)->desc(d->cellno(row, 0         ), row->dir);
    auto endbc   = (d->boundaries[row->end_code  ]->*mem)->desc(d->cellno(row, row->n - 1), row->dir);

    bb[0] = 1.0 - (2.0 * startbc.sw - 3.0) * d_row[0] * dt_dx2;
    phi_row[0] += 2.0 * d_row[0] * dt_dx2 * startbc.cte;

    bb[n - 1] = 1.0 - (2.0 * endbc.sw - 3.0) * d_row[n - 1] * dt_dx2;
    phi_row[n - 1] += 2.0 * d_row[n - 1] * dt_dx2 * endbc.cte;

    solve_tridiagonal_in_place_destructive(phi_row, n, aa, bb, cc);

    d->insert_scalars(row, phi, phi_row);
    delete[] phi_row;
    delete[] d_row;

    delete aa;
    delete bb;
    delete cc;
}

void diffusion::diffuse_ustar()
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
                diffuse(d->rows[dir] + irow, d->ustar[icmpnt], d->nu, flowbc::umembers[icmpnt]);
}

}


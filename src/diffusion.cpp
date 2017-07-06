/**
 * diffusion.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "diffusion.h"

namespace aban2
{

// https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
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
    for (in = N - 1; in-- > 0;)
        x[in] = x[in] - c[in] * x[in + 1];
}

void diffusion::diffuse(row_t *row, double *phi, double *D, flowbc::member mem)
{
    double *phi_row = d->extract_scalars(row, phi);
    double *d_row = d->extract_scalars(row, D);
    double dx = d->delta, dt = d->dt;
    size_t n = row->n;
    double dt_dx2 = dt / dx / dx;

    double *aa = new double[n], *bb = new double[n], *cc = new double[n];
    // a: subdiagonal, b: the main diagonal, c: superdiagonal

    for (size_t i = 1; i < n - 1; ++i)
    {
        double d_e = (d_row[i + 1] + d_row[i]) / 2;
        double d_w = (d_row[i] + d_row[i - 1]) / 2;

        aa[i] = -d_w * dt_dx2;
        bb[i] = 1.0 + (d_w + d_e) * dt_dx2;
        cc[i] = -d_e * dt_dx2;
    }

    // start bc
    {
        double d_e = (d_row[0] + d_row[1]) / 2;
        double d_w = d_row[0];

        auto startbc = (d->boundaries[(int)row->start_code]->*mem)->desc(d->cellno(row, 0), row->dir);
        cc[0] = -d_e * dt_dx2;
        bb[0] = 1.0 - ((d_e + d_w) * (startbc.sw - 1) - d_w) * dt_dx2;
        phi_row[0] += 2.0 * d_w * dt_dx2 * startbc.cte;
    }

    // end bc
    {
        double d_e = d_row[n - 1];
        double d_w = (d_row[n - 2] + d_row[n - 1]) / 2;

        auto endbc = (d->boundaries[(int)row->end_code]->*mem)->desc(d->cellno(row, row->n - 1), row->dir);
        aa[n - 1] = -d_w * dt_dx2;
        bb[n - 1] = 1.0 - ((d_e + d_w) * (endbc.sw - 1) - d_e) * dt_dx2;
        phi_row[n - 1] += 2.0 * d_e * dt_dx2 * endbc.cte;
    }

    solve_tridiagonal_in_place_destructive(phi_row, n, aa, bb, cc);

    d->insert_scalars(row, phi, phi_row);
    delete[] phi_row;
    delete[] d_row;

    delete[] aa;
    delete[] bb;
    delete[] cc;
}

void diffusion::diffuse_ustar()
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            for (size_t icmpnt = 0; icmpnt < NDIRS; ++icmpnt)
                diffuse(d->rows[dir] + irow, d->ustar[icmpnt], d->nu, flowbc::umembers[icmpnt]);
}
}

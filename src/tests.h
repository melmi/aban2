/**
 * tests.h
 *
 * by Mohammad Elmi
 * Created on 24 Sep 2013
 */

#ifndef _TESTS_H_
#define _TESTS_H_

#include <iostream>
#include <cmath>
#include <sstream>

#include "domain.h"

using namespace std;

namespace aban2
{

string rowtostr(mesh_row *r)
{
    std::stringstream s;
    s << "n: (" << r->n << ")  ";
    s << "codes: (" << r->start_code << ", " << r->end_code << ")  ";
    s << "principal dir: (" << r->dir << ")  ";
    s << "start: (" << r->start[0] << ", " << r->start[1] << ", " << r->start[2] << ")  ";
    s << "end: (" << r->end[0] << ", " << r->end[1] << ", " << r->end[2] << ")  ";
    return s.str();
}

void print_cell_nos(domain *d)
{
    for (int i = d->ndir[0]-1; i > 0; --i)
    {
        for (int j = 1; j < d->ndir[1]-1; ++j)
        {
            size_t x = d->idx(i, j, 0);
            cout.width(4);
            if (d->codes[x] == mesh::INSIDE)
                cout << d->cellnos[x];
            // else
            //     cout << d->codes[x];
        }

        cout << endl;
    }
}

void print_rows(domain *d)
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
    {
        cout << "======= dir: " << dir << endl;
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            cout << rowtostr(row) << endl;

            // auto a = d->get_row_idxs(row);
            for (int i = 0; i < row->n; ++i)
                cout << d->cellno(row, i) << " ";
            // cout << a[i] << " ";
            cout << endl;
            // delete[] a;
        }
        cout << endl;
    }
}


void adv_test1(domain *d)
{
    double x0 = 51, y0 = 51;
    for (int i = 0; i < d->ndir[0]; ++i)
        for (int j = 0; j < d->ndir[1]; ++j)
            for (int k = 0; k < d->ndir[2]; ++k)
                if (d->exists(i, j, k))
                {
                    size_t ix = d->cellnos[d->idx(i, j, k)];
                    // cout<<ix<<endl<<flush;
                    d->u[0][ix] = d->u[1][ix] = d->u[2][ix] = 1;
                    double x = i * d->delta + d->delta / 2.0;
                    double y = j * d->delta + d->delta / 2.0;
                    double z = k * d->delta + d->delta / 2.0;
                    double phi = exp(-( (x - x0) * (x - x0) + (y - y0) * (y - y0)) / 100.0);
                    d->ustar[0][ix] = d->ustar[1][ix] = d->ustar[2][ix] = phi;
                }
}

void diff_test1(domain *d)
{
    double x0 = 21, y0 = 21;
    for (int i = 0; i < d->ndir[0]; ++i)
        for (int j = 0; j < d->ndir[1]; ++j)
            for (int k = 0; k < d->ndir[2]; ++k)
                if (d->exists(i, j, k))
                {
                    size_t ix = d->cellnos[d->idx(i, j, k)];
                    // cout<<ix<<endl<<flush;
                    double x = i * d->delta + d->delta / 2.0;
                    double y = j * d->delta + d->delta / 2.0;
                    double z = k * d->delta + d->delta / 2.0;
                    double phi =  exp(-( (x - x0) * (x - x0) + (y - y0) * (y - y0)) / 100.0);
                    d->ustar[0][ix] = d->ustar[1][ix] = d->ustar[2][ix] = phi;
                }
}

void vortex(domain *d, double w, double x0, double y0)
{
    for (int i = 0; i < d->ndir[0]; ++i)
        for (int j = 0; j < d->ndir[1]; ++j)
            for (int k = 0; k < d->ndir[2]; ++k)
                if (d->exists(i, j, k))
                {
                    size_t ix = d->cellnos[d->idx(i, j, k)];

                    double x = i * d->delta + d->delta / 2.0;
                    double y = j * d->delta + d->delta / 2.0;
                    double z = k * d->delta + d->delta / 2.0;

                    d->u[0][ix] = -w * (y - y0);// + 2.0 * x;
                    d->u[1][ix] = w * (x - x0);
                    d->u[2][ix] = 0;
                }
}

}

#endif
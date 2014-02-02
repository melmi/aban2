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
#include "ireconst.h"

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
    for (int i = d->ndir[0] - 1; i > 0; --i)
    {
        for (int j = 1; j < d->ndir[1] - 1; ++j)
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
                    d->q[0][ix] = d->q[1][ix] = d->q[2][ix] = 1;
                    double x = i * d->delta + d->delta / 2.0;
                    double y = j * d->delta + d->delta / 2.0;
                    double z = k * d->delta + d->delta / 2.0;
                    double phi = exp(-( (x - x0) * (x - x0) + (y - y0) * (y - y0)) / 100.0);
                    d->qstar[0][ix] = d->qstar[1][ix] = d->qstar[2][ix] = phi;
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
                    d->qstar[0][ix] = d->qstar[1][ix] = d->qstar[2][ix] = phi;
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

                    d->q[0][ix] = -w * (y - y0);// + 2.0 * x;
                    d->q[1][ix] = w * (x - x0);
                    d->q[2][ix] = 0;
                }
}

void zalesak_disk_2d(domain *d)
{
    vector x0 {50, 75, 0};
    double r = 15;

    double h_2 = d->delta / 2, r2 = r * r, v = d->delta * d->delta * d->delta;
    vector half[] = {{h_2, h_2, 0}, {h_2, -h_2, 0}, { -h_2, h_2, 0}, { -h_2, -h_2, 0}};
    vector h {d->delta, d->delta, d->delta};

    // creating circle
    for (size_t j = 0; j < d->ndir[1]; ++j)
        for (size_t i = 0; i < d->ndir[0]; ++i)
        {
            size_t ix = d->cellnos[d->idx(i, j, 0)];
            vector x {d->delta * i, d->delta * j, 0};
            double l2[] = {(x + half[0] - x0).l2(), (x + half[1] - x0).l2(), (x + half[2] - x0).l2(), (x + half[3] - x0).l2()};
            if (*std::min_element(begin(l2), end(l2)) > r2)
            {
                d->vof[ix] = 0;
                continue;
            }
            if (*std::max_element(begin(l2), end(l2)) < r2)
            {
                d->vof[ix] = 1;
                continue;
            }

            double min_dist = std::sqrt(*std::min_element(begin(l2), end(l2)));
            vector n = x - x0;
            n.normalize();
            double alpha = r - min_dist;

            ireconst reconst = ireconst::from_alpha(h, n, alpha);
            d->vof[ix] = reconst.volume / v;
        }

    //removing slot
    double eps = 1e-4;
    vector c {50, 72.5, 0}, l {5 + eps, 25 + eps, 0};
    for (size_t j = 0; j < d->ndir[1]; ++j)
        for (size_t i = 0; i < d->ndir[0]; ++i)
        {
            size_t ix = d->cellnos[d->idx(i, j, 0)];
            vector x {d->delta * i, d->delta * j, 0};
            if (
                (x.x > c.x - l.x / 2.0) &&
                (x.x < c.x + l.x / 2.0) &&
                (x.y > c.y - l.y / 2.0) &&
                (x.y < c.y + l.y / 2.0)) d->vof[ix] = 0;
        }

    //setting velocities
    double pi=std::atan2(0, -1);
    double pi_314 = pi / 314.0;
    for (size_t j = 0; j < d->ndir[1]; ++j)
        for (size_t i = 0; i < d->ndir[0]; ++i)
        {
            size_t ix = d->cellnos[d->idx(i, j, 0)];
            vector x {d->delta *i + h_2, d->delta *j + h_2, 0}; // interface u is a half cell staggered

            d->uf[0][ix] = pi_314 * (50. - x.y);
            d->uf[1][ix] = pi_314 * (x.x - 50.);
        }
}

}

#endif
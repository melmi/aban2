/**
 * tests.h
 *
 * by Mohammad Elmi
 * Created on 24 Sep 2013
 */

#ifndef _TESTS_H_
#define _TESTS_H_

#include "common.h"
#include <iostream>
#include <cmath>
#include <sstream>

#include "domain.h"
#include "volreconst.h"
#include "vof.h"

using namespace std;

namespace aban2
{

//-------------------------- utility functions

template <typename T>
string to_string (T number)
{
    ostringstream ss;
    ss << number;
    return ss.str();
}

void check_continuety(domain *d)
{
    double *div = new double[d->n];
    std::fill_n(div, d->n, 0);

    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 1; irow < d->nrows[dir] - 1; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            for (size_t i = 1; i < row->n - 1; ++i)
            {
                double *u = d->extract_scalars(row, d->uf[row->dir]);

                size_t no = d->cellno(row, i);
                div[no] = u[i] - u[i - 1];

                delete[] u;
            }
        }

    std::for_each(div, div + d->n, [](double & x)
    {
        x = std::abs(x);
    });
    cout << "divergance: " << std::accumulate(div, div + d->n, 0) << std::endl;
}

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

template <class I>
std::size_t min_element_index ( I first, I last )
{
    I lowest = first;
    std::size_t index = 0;
    std::size_t i = 0;
    if (first == last) return index;
    while (++first != last)
    {
        ++i;
        if (*first < *lowest)
        {
            lowest = first;
            index = i;
        }
    }
    return index;
}

//-------------------------- tests

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

void circle(domain *d, vector x0, double r);
void vortex(domain *d, vector x0, double omega);

void vof_reconst_accuracy_test()
{
    string files[]
    {
        "mesh/cavity10x10.json",
        "mesh/cavity20x20.json",
        "mesh/cavity40x40.json",
        "mesh/cavity80x80.json",
        "mesh/cavity160x160.json"
    };
    // bool first = true;
    for (auto f : files)
    {
        cout << "========" << endl;
        cout << "input file: " << f << endl;
        cout << "Reading mesh" << endl;
        domain *d = domain::create_from_file(f);
        cout << "initializing" << endl;
        circle(d, {0.5, 0.5, 0}, 0.30);
        vof vofc(d);
        cout << "calculating normals" << endl;
        //vofc.calculate_normals();
        cout << "writing" << endl;
        // if (first)
        //     d->write_vtk("out/my-method.vtk");
        // first = false;
        cout << "comparing" << endl;
        cout << vof_err::compare(vofc) << endl;
        //delete d;
    }

    cout << "done" << endl;
}

void zalesak_disk(domain *);

void zalesak_disk_rotation_test()
{
    domain *d = domain::create_from_file("mesh/cavity100x100.json");
    vof _vof(d);
    vortex(d, {50, 50, 0}, 0.01);
    zalesak_disk(d);

    d->write_vtk("out/zalesak0.vtk");
    for (int i = 0; i < 700; ++i)
    {
        std::cout << "step " << i + 1 << std::endl << std::flush;
        _vof.advect();

        stringstream ss;
        ss << "out/zalesak" << i + 1 << ".vtk";
        d->write_vtk(ss.str());
    }


    delete d;
}

//-------------------------- vof shapes

void circle(domain *d, vector x0, double r)
{
    double h_2 = d->delta / 2, r2 = r * r, v = d->delta * d->delta * d->delta;
    vector half[] = {{h_2, h_2, 0}, {h_2, -h_2, 0}, { -h_2, h_2, 0}, { -h_2, -h_2, 0}};
    vector h {d->delta, d->delta, d->delta};

    size_t no;
    // creating circle
    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            if (d->exists(i, j, 0, no))
            {
                vector x {d->delta *i + h_2, d->delta *j + h_2, 0};
                vector ps[] = {x - x0 + half[0], x - x0 + half[1], x - x0 + half[2], x - x0 + half[3]};
                double l2[] = {ps[0].l2(), ps[1].l2(), ps[2].l2(), ps[3].l2()};

                vector n = x - x0;
                n.normalize();
                n.to_data(d->u, no);

                if (*std::min_element(begin(l2), end(l2)) > r2)
                {
                    d->vof[no] = 0;
                    continue;
                }
                if (*std::max_element(begin(l2), end(l2)) < r2)
                {
                    d->vof[no] = 1;
                    continue;
                }

                size_t idx = min_element_index(l2, l2 + 4);
                double alpha = r - ps[idx] * n;
                volreconst *reconst = volreconst::from_alpha(h, n, alpha);
                d->vof[no] = reconst->volume / v;
                delete reconst;
            }
}

void rectangel(domain *d, vector c, vector l, double value)
{
    double h_2 = d->delta / 2.0;
    size_t no;

    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            if (d->exists(i, j, 0, no))
            {
                vector x {d->delta *i + h_2, d->delta *j + h_2, 0};
                if (
                    (x.x > c.x - l.x / 2.0) &&
                    (x.x < c.x + l.x / 2.0) &&
                    (x.y > c.y - l.y / 2.0) &&
                    (x.y < c.y + l.y / 2.0))
                    d->vof[no] = value;
            }
}

void square(domain *d, vector c, double a)
{
    rectangel(d, c, {a, a, a}, 1);
}

void zalesak_disk(domain *d)
{
    circle(d, {50, 75, 0}, 15);
    rectangel(d, {50, 72.5, 0}, {5 , 25 , 0}, 0);
}

//-------------------------- velocity fields

void vortex(domain *d, vector x0, double omega)
{
    size_t no;
    double h_2 = d->delta / 2.0;
    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            if (d->exists(i, j, 0, no))
            {
                vector x {d->delta *(i + 1), d->delta *(j + 1), 0};

                d->uf[0][no] = omega * (x0.y - x.y);
                d->uf[1][no] = omega * (x.x - x0.x);
            }
}

}

#endif
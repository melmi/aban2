/**
 * tests.h
 *
 * by Mohammad Elmi
 * Created on 24 Sep 2013
 */

#ifndef _TESTS_H_
#define _TESTS_H_

#include "common.h"
#include "debug.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>

#include "domain.h"
#include "vof.h"
#include "volreconst.h"

using namespace std;

namespace aban2
{

//-------------------------- utility functions

void check_continuety(domain_t *d)
{
    double *div = new double[d->n];
    std::fill_n(div, d->n, 0);

    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 1; irow < d->nrows[dir] - 1; ++irow)
        {
            row_t *row = d->rows[dir] + irow;
            for (size_t i = 1; i < row->n - 1; ++i)
            {
                double *u = d->extract_scalars(row, d->uf[row->dir]);

                size_t no = d->cellno(row, i);
                div[no] = u[i] - u[i - 1];

                delete[] u;
            }
        }

    std::for_each(div, div + d->n, [](double &x) {
        x = std::abs(x);
    });
    cout << "divergance: " << std::accumulate(div, div + d->n, 0) << std::endl;
}

string rowtostr(row_t *r)
{
    std::stringstream s;
    s << "n: (" << r->n << ")  ";
    s << "codes: (" << r->start_code << ", " << r->end_code << ")  ";
    s << "principal dir: (" << r->dir << ")  ";
    s << "start: (" << r->start[0] << ", " << r->start[1] << ", " << r->start[2] << ")  ";
    s << "end: (" << r->end[0] << ", " << r->end[1] << ", " << r->end[2] << ")  ";
    return s.str();
}

void print_cell_nos(domain_t *d)
{
    for (size_t i = d->ndir[0] - 1; i > 0; --i)
    {
        for (size_t j = 1; j < d->ndir[1] - 1; ++j)
        {
            size_t x = d->idx(i, j, 0);
            cout.width(4);
            // if (d->codes[x] == mesh::INSIDE)
            cout << d->cellnos[x];
            // else
            //     cout << d->codes[x];
        }

        cout << endl;
    }
}

void print_rows(domain_t *d)
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
    {
        cout << "======= dir: " << dir << endl;
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            row_t *row = d->rows[dir] + irow;
            cout << rowtostr(row) << endl;

            // auto a = d->get_row_idxs(row);
            for (size_t i = 0; i < row->n; ++i)
                cout << d->cellno(row, i) << " ";
            // cout << a[i] << " ";
            cout << endl;
            // delete[] a;
        }
        cout << endl;
    }
}

template <class I>
std::size_t min_element_index(I first, I last)
{
    I lowest = first;
    std::size_t index = 0;
    std::size_t i = 0;
    if (first == last)
        return index;
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

void adv_test1(domain_t *d)
{
    double x0 = 51, y0 = 51;
    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            for (size_t k = 0; k < d->ndir[2]; ++k)
                if (d->is_inside(i, j, k))
                {
                    size_t ix = d->cellnos[d->idx(i, j, k)];
                    // cout<<ix<<endl<<flush;
                    d->u[0][ix] = d->u[1][ix] = d->u[2][ix] = 1;
                    double x = i * d->delta + d->delta / 2.0;
                    double y = j * d->delta + d->delta / 2.0;
                    /* double z = k * d->delta + d->delta / 2.0; */
                    double phi = exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / 100.0);
                    d->ustar[0][ix] = d->ustar[1][ix] = d->ustar[2][ix] = phi;
                }
}

void diff_test1(domain_t *d)
{
    double x0 = 21, y0 = 21;
    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            for (size_t k = 0; k < d->ndir[2]; ++k)
                if (d->is_inside(i, j, k))
                {
                    size_t ix = d->cellnos[d->idx(i, j, k)];
                    // cout<<ix<<endl<<flush;
                    double x = i * d->delta + d->delta / 2.0;
                    double y = j * d->delta + d->delta / 2.0;
                    /* double z = k * d->delta + d->delta / 2.0; */
                    double phi = exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / 100.0);
                    d->ustar[0][ix] = d->ustar[1][ix] = d->ustar[2][ix] = phi;
                }
}

void square(domain_t *d, vector c, double a);
void circle(domain_t *d, vector x0, double r);
void vortex(domain_t *d, vector x0, double omega);

void vof_reconst_accuracy_test()
{
    string files[]{
        "mesh/cavity10x10.json",
        "mesh/cavity20x20.json",
        "mesh/cavity40x40.json",
        "mesh/cavity80x80.json",
        "mesh/cavity160x160.json"};

    for (auto f : files)
    {
        cout << "========" << endl;
        cout << "input file: " << f << endl;
        cout << "Reading mesh" << endl;
        domain_t *d = domain_t::create_from_file(f);
        cout << "initializing" << endl;
        circle(d, {0.5, 0.5, 0}, 0.30);
        vof_t *vofc = new vof_t(d);
        cout << "calculating normals" << endl;
        vofc->calculate_normals();
        cout << "writing" << endl;
        d->write_vtk("out/" + f.substr(5) + ".vtk");
        cout << "comparing" << endl;
        cout << vof_err_t::compare(*vofc) << endl;
        delete vofc;
        delete d;
    }

    cout << "done" << endl;
}

void zalesak_disk(domain_t *d, vector center = {0.5, 0.75, 0});
void fill_func(domain_t *d, double *phi, std::function<double(vector)> f);

double error_in_rectangel(domain_t *d, vector c, vector l, double *phi1, double *phi2)
{
    double h_2 = d->delta / 2.0;
    size_t no;
    double result = 0;

    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            if (d->is_inside(i, j, 0, no))
            {
                vector x{d->delta * i - h_2, d->delta * j - h_2, 0};
                if (
                    (x.x > c.x - l.x / 2.0) &&
                    (x.x < c.x + l.x / 2.0) &&
                    (x.y > c.y - l.y / 2.0) &&
                    (x.y < c.y + l.y / 2.0))
                {
                    double e = std::abs(phi1[no] - phi2[no]);
                    result += e;
                }
            }
    return result;
}

void consistent_vof_test()
{
    string files[]{
        "mesh/cavity20x20.json",
        "mesh/cavity40x40.json",
        "mesh/cavity80x80.json",
        "mesh/cavity160x160.json"};

    for (auto f : files)
    {
        cout << "========" << endl;
        cout << "input file: " << f << endl;
        domain_t *d = domain_t::create_from_file(f);
        vof_t *vof = new vof_t(d);

        double dt = 1, dx = d->delta, CFL = 0.5;
        d->dt = dt;
        double u = CFL * dx / dt;
        double omega = u / 0.5; //0.25;
        int nSteps = 2 * 3.1415 / (omega * dt);

        vortex(d, {0.5, 0.5, 0}, omega);
        vector x0{0.5, 0.7, 0};
        zalesak_disk(d, x0);
        // circle(d, {0.5, 0.5, 0}, 0.15);
        double sigma = 0.15;
        double two_sigma2 = 2 * sigma * sigma;
        fill_func(d, d->u[0], [x0, two_sigma2](vector x) {
            return std::exp(-((x - x0) * (x - x0)) / two_sigma2);
        });
        fill_n(d->u[1], d->n, 0);
        for (size_t i = 0; i < d->n; ++i)
            d->rho[i] = d->rho_bar(d->vof[i]);

        double *dd = new double[d->n];
        double *vv = new double[d->n];
        copy_n(d->u[0], d->n, dd);
        copy_n(d->vof, d->n, vv);

        vof->calculate_normals();
        d->write_vtk("out/zalesak0.vtk");
        for (int i = 0; i < nSteps; ++i)
        {
            // std::cout << "step " << i + 1 << std::endl << std::flush;
            std::copy_n(d->u[0], d->n, d->ustar[0]);
            vof->advect();

            vof->calculate_normals();
            std::copy_n(d->ustar[0], d->n, d->u[0]);
            // d->write_vtk("out/zalesak" + to_string(i + 1) + ".vtk");
        }
        d->write_vtk("out/zalesak1.vtk");

        std::cout << "err dd: " << error_in_rectangel(d, x0, {0.4, 0.4, 0}, dd, d->u[0]) << std::endl;
        std::cout << "err vv: " << error_in_rectangel(d, x0, {0.4, 0.4, 0}, vv, d->vof) << std::endl;

        delete vof;
        delete d;
        delete dd;
    }
}

void zalesak_disk_rotation_test()
{
    domain_t *d = domain_t::create_from_file("mesh/cavity160x160.json");
    // for (size_t i = 0; i < d->n; ++i) d->p[i] = i;

    vof_t *_vof = new vof_t(d);
    // fill_func(d, d->uf[1], [](vector x) {return -0.004;});
    vortex(d, {0.5, 0.5, 0}, 0.01);
    // fill_func(d, d->ustar[0], [](vector x) {return x* vector{1, 1, 0};});
    // square(d, {0.5, 0.75, 0}, 0.15);
    zalesak_disk(d);
    fill_func(d, d->u[0], [](vector x) {
        vector x0{0.5, 0.75, 0};
        double sigma2 = 0.2 * 0.2;
        return std::exp(-((x - x0) * (x - x0)) / sigma2);
    });
    for (size_t i = 0; i < d->n; ++i)
        d->rho[i] = d->rho_bar(d->vof[i]);

    _vof->calculate_normals();
    d->write_vtk("out/zalesak0.vtk");
    for (int i = 0; i < 700; ++i)
    {
        std::cout << "step " << i + 1 << std::endl
                  << std::flush;
        std::copy_n(d->u[0], d->n, d->ustar[0]);
        _vof->advect();

        _vof->calculate_normals();
        std::copy_n(d->ustar[0], d->n, d->u[0]);
        d->write_vtk("out/zalesak" + to_string(i + 1) + ".vtk");
    }

    delete _vof;
    delete d;
}

//-------------------------- vof shapes

/*
       y
       ^
       |
       +---+---+---+---+
       |   |   |   |   |
       +---+---+---+---+
       |   |   |   |   |
       +---o---+---+---+-> x
       |   |   |   |   |
       +---+---+---+---+

*/

void circle(domain_t *d, vector x0, double r)
{
    double h_2 = d->delta / 2, r2 = r * r, v = d->delta * d->delta * d->delta;
    vector half[] = {{h_2, h_2, 0}, {h_2, -h_2, 0}, {-h_2, h_2, 0}, {-h_2, -h_2, 0}};
    vector h{d->delta, d->delta, d->delta};

    size_t no;
    // creating circle
    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            if (d->is_inside(i, j, 0, no))
            {
                vector x{d->delta * i - h_2, d->delta * j - h_2, 0};
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
                d->vof[no] = reconst->get_volume() / v;
                delete reconst;
            }
}

void rectangel(domain_t *d, vector c, vector l, double value)
{
    double h_2 = d->delta / 2.0;
    size_t no;

    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            if (d->is_inside(i, j, 0, no))
            {
                vector x{d->delta * i - h_2, d->delta * j - h_2, 0};
                if (
                    (x.x > c.x - l.x / 2.0) &&
                    (x.x < c.x + l.x / 2.0) &&
                    (x.y > c.y - l.y / 2.0) &&
                    (x.y < c.y + l.y / 2.0))
                    d->vof[no] = value;
            }
}

void square(domain_t *d, vector c, double a)
{
    rectangel(d, c, {a, a, a}, 1);
}

void zalesak_disk(domain_t *d, vector center)
{
    circle(d, center, 0.15);
    rectangel(d, center - vector{0, 0.1, 0}, {0.05, 0.3, 0}, 0);
}

//-------------------------- velocity fields

void fill_func(domain_t *d, double *phi, std::function<double(vector)> f)
{
    size_t no;
    /* double h_2 = d->delta / 2.0; */
    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            if (d->is_inside(i, j, 0, no))
            {
                vector x{d->delta * i, d->delta * j, 0};
                phi[no] = f(x);
            }
}

void vortex(domain_t *d, vector x0, double omega)
{
    fill_func(d, d->uf[0], [x0, omega](vector x) {
        return omega * (x0.y - x.y);
    });
    fill_func(d, d->uf[1], [x0, omega](vector x) {
        return omega * (x.x - x0.x);
    });
}
}

#endif

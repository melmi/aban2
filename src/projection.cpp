/**
 * projection.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include <iostream>
#include "projection.h"
#include "gradient.h"

namespace aban2
{

projection::projection(domain *_d): d(_d)
{
    h2inv = 1.0 / (d->delta * d->delta);
}

projection::~projection()
{
}

lpw::matrix_t *projection::create_matrix()
{
    lpw::matrix_t *pmatrix = new lpw::matrix_t(d->n);

    for (size_t idir = 0; idir < NDIRS; ++idir)
        for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
        {
            mesh::row *row = d->rows[idir] + irow;
            add_row(pmatrix, row);
        }

    return pmatrix;
}

void projection::add_row(lpw::matrix_t *pmatrix, mesh::row *row)
{
    size_t n = row->n;
    size_t *row_cellnos = d->get_row_cellnos(row);

    for (size_t i = 1; i < n - 1; ++i)
    {
        size_t no_w = row_cellnos[i - 1], no_p = row_cellnos[i], no_e = row_cellnos[i + 1];
        double rho_w = d->rho[no_w], rho_p = d->rho[no_p], rho_e = d->rho[no_e];
        double _rho_f_w = 2.0 / (rho_w + rho_p), _rho_f_e = 2.0 / (rho_e + rho_p);

        (*pmatrix)(no_p, no_w) += +h2inv * _rho_f_w;
        (*pmatrix)(no_p, no_p) += -h2inv * (_rho_f_e + _rho_f_w);
        (*pmatrix)(no_p, no_e) += +h2inv * _rho_f_e;
    }

    auto desc_start = d->boundaries[(int)row->start_code]->p->desc(d->cellno(row, 0         ), row->dir);
    auto desc_end   = d->boundaries[(int)row->end_code  ]->p->desc(d->cellno(row, row->n - 1), row->dir);
    double rho_start = d->rho_bar(flowbc::bc_vof_getter(d, row, d->vof, bcside::start));
    double rho_end   = d->rho_bar(flowbc::bc_vof_getter(d, row, d->vof, bcside::end  ));
    apply_row_bc(pmatrix, row_cellnos[0    ], row_cellnos[1    ], desc_start, rho_start);
    apply_row_bc(pmatrix, row_cellnos[n - 1], row_cellnos[n - 2], desc_end  , rho_end  );

    delete[] row_cellnos;
}

void projection::apply_row_bc(lpw::matrix_t *pmatrix, size_t no0 , size_t no1, bcdesc desc, double rho_b)
{
    double rho_0 = d->rho[no0],
           rho_1 = d->rho[no1],
           _rho_b = 1.0 / rho_b,
           _rho_f = 2.0 / (rho_0 + rho_1);

    double c0 = 2.0 * desc.sw - 2.0;
    double c1 = -1.0;

    (*pmatrix)(no0, no0) += (c0 * _rho_b + c1 * _rho_f) * h2inv;
    (*pmatrix)(no0, no1) += _rho_f * h2inv;
}

double *projection::get_rhs()
{
    auto rhs = gradient::divergance_of(d, d->ustar, flowbc::bc_u_getter);
    for (size_t i = 0; i < d->n; ++i) rhs[i] /= d->dt;
    apply_rhs_bc(rhs);
    return rhs;
}

void projection::apply_rhs_bc(double *rhs)
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh::row *row = d->rows[dir] + irow;

            size_t cellno_start = d->cellno(row, 0);
            size_t cellno_end   = d->cellno(row, row->n - 1);

            auto bcdesc_start = d->boundaries[(int)row->start_code]->p->desc(cellno_start, row->dir);
            auto bcdesc_end   = d->boundaries[(int)row->end_code  ]->p->desc(cellno_end  , row->dir);

            double _rho_start = 1.0 / d->rho_bar(flowbc::bc_vof_getter(d, row, d->vof, bcside::start));
            double _rho_end   = 1.0 / d->rho_bar(flowbc::bc_vof_getter(d, row, d->vof, bcside::end  ));

            rhs[cellno_start] -= 2.0 * bcdesc_start.cte * _rho_start * h2inv;
            rhs[cellno_end  ] -= 2.0 * bcdesc_end  .cte * _rho_end   * h2inv;

            // rhs[cellno_start] += bcdesc_start.grad * _rho_start / d->delta;
            // rhs[cellno_end  ] -= bcdesc_end  .grad * _rho_end   / d->delta;
        }
}

void projection::solve_p()
{
    auto pmatrix = create_matrix();
    auto rhs = get_rhs();
    auto niter = cg(pmatrix, d->p, rhs, d->n, lpw::precond_t::jacobi, 1.2);
    std::cout << "                   #iterations: " << niter << std::endl;
    delete[] rhs;
    delete pmatrix;
}

double p_f(double p_P, double rho_P, double p_N, double rho_N)
{
    return (p_P * rho_N + p_N * rho_P) / (rho_P + rho_N);
}

double **gradient_of_p(domain *d)
{
    double **result = new double*[3];
    result[2] = nullptr;
    for (size_t dir = 0; dir < NDIRS; ++dir)
    {
        result[dir] = (double *)d->create_var(1);
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            auto row = d->rows[dir] + irow;
            double *p_row = d->extract_scalars(row, d->p);
            double *rho_row = d->extract_scalars(row, d->rho);
            double *grad_row = d->extract_scalars(row, result[dir]);
            size_t n = row->n;

            for (size_t i = 1; i < n - 1; ++i)
                grad_row[i] += (p_f(p_row[i], rho_row[i], p_row[i + 1], rho_row[i + 1]) -
                                p_f(p_row[i], rho_row[i], p_row[i - 1], rho_row[i - 1])) / d->delta;

            auto startbc = flowbc::bc_p_getter(d, row, d->p, bcside::start);
            auto endbc   = flowbc::bc_p_getter(d, row, d->p, bcside::end);

            grad_row[0] += (p_f(p_row[0], rho_row[0], p_row[1], rho_row[1]) - startbc) / d->delta;
            grad_row[n - 1] += (endbc - p_f(p_row[n - 1], rho_row[n - 1], p_row[n - 2], rho_row[n - 2])) / d->delta;

            d->insert_scalars(row, result[dir], grad_row);
            delete[] p_row;
            delete[] rho_row;
            delete[] grad_row;
        }
    }
    return result;
}

void projection::update_u()
{
    double **gradp = gradient_of_p(d);

    for (size_t i = 0; i < d->n; ++i)
        for (int dir = 0; dir < NDIRS; ++dir)
            d->u[dir][i] = d->ustar[dir][i] - gradp[dir][i] * d->dt / d->rho[i];

    d->delete_var(2, gradp);
}

void projection::update_uf()
{
    for (size_t idir = 0; idir < NDIRS; ++idir)
        for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
        {
            mesh::row *row = d->rows[idir] + irow;

            double *ustar = d->extract_scalars(row, d->ustar[idir]);
            double *p = d->extract_scalars(row, d->p);
            double *rho = d->extract_scalars(row, d->rho);
            double *uf = new double[row->n];

            for (size_t i = 0; i < row->n - 1; ++i)
                uf[i] = (ustar[i] + ustar[i + 1]) / 2.0
                        - (p[i + 1] - p[i]) / d->delta * d->dt / ((rho[i] + rho[i + 1]) / 2.0);

            d->insert_scalars(row, d->uf[idir], uf);

            delete[] ustar;
            delete[] p;
            delete[] uf;
            delete[] rho;
        }
}

}

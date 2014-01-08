/**
 * projection.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "projection.h"
#include "gradient.h"

namespace aban2
{

projection::projection(domain *_d): d(_d)
{
    h2inv = 1.0 / (d->delta * d->delta);

    pmatrix = new matrix_t(d->n, d->n);

    make_matrix();

    psolver = new solver_t(*pmatrix);
}

projection::~projection()
{
    delete pmatrix;
    delete psolver;
}

void projection::make_matrix()
{
    triplet_vector *coeffs = new triplet_vector();

    for (size_t idir = 0; idir < NDIRS; ++idir)
        for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
        {
            mesh_row *row = d->rows[idir] + irow;
            add_row(row, coeffs);
        }

    pmatrix->setFromTriplets(coeffs->begin(), coeffs->end());

    delete coeffs;
}

void projection::add_row(mesh_row *row, triplet_vector *coeffs)
{
    size_t n = row->n;
    size_t *row_idxs = d->get_row_idxs(row);
    auto start_bc_desc = d->boundaries[row->start_code]->p(row, bcside::start, 0);
    auto end_bc_desc   = d->boundaries[row->end_code  ]->p(row, bcside::end  , 0);

    for (size_t i = 1; i < n - 1; ++i)
    {
        size_t ix = row_idxs[i];
        coeffs->push_back({ix, row_idxs[i - 1], +1.0 * h2inv});
        coeffs->push_back({ix, ix             , -2.0 * h2inv});
        coeffs->push_back({ix, row_idxs[i + 1], +1.0 * h2inv});
    }

    apply_row_bc(row_idxs[0    ], row_idxs[1    ], start_bc_desc, coeffs);
    apply_row_bc(row_idxs[n - 1], row_idxs[n - 2], end_bc_desc  , coeffs);

    delete[] row_idxs;
}

void projection::apply_row_bc(size_t ix0 , size_t ix1, bcdesc desc, triplet_vector *coeffs)
{
    coeffs->push_back({ix0, ix0, (2.0 * desc.coeff - 3.0)*h2inv});
    coeffs->push_back({ix0, ix1, +1.0 * h2inv});
}

double *projection::get_rhs()
{
    double *rhs = (double *)d->create_var(1);
    gradient::divergance(d, d->ustar, rhs, &bcondition::u);
    double rho_dt = d->rho / d->dt;
    for (int i = 0; i < d->n; ++i) rhs[i] *= rho_dt;
    apply_rhs_bc(rhs);
    return rhs;
}

void projection::apply_rhs_bc(double *rhs)
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;

            apply_single_rhs_bc(row, rhs, bcside::start);
            apply_single_rhs_bc(row, rhs, bcside::end);
        }
}

void projection::apply_single_rhs_bc(mesh_row *row, double *rhs, bcside side)
{
    bcondition *bc;
    if (side == bcside::start)
        bc = d->boundaries[row->start_code];
    else
        bc = d->boundaries[row->end_code];

    auto desc = bc->p(row, bcside::start, row->dir);

    rhs[desc.cellno] -= 2.0 * desc.val * h2inv;
}

void projection::solve_p()
{
    double *rhs = get_rhs();
    Eigen::VectorXd b(d->n), x(d->n);
    for (int i = 0; i < d->n; ++i) b[i] = rhs[i];
    delete[] rhs;
    x = psolver->solve(b);
    for (int i = 0; i < d->n; ++i) d->p[i] = x[i];
}

void projection::update_u()
{
    double **gradp = (double **)d->create_var(2);
    gradient::of_scalar(d, d->p, gradp, &bcondition::p);

    for (int i = 0; i < d->n; ++i)
        for (int dir = 0; dir < NDIRS; ++dir)
            d->u[dir][i] = d->ustar[dir][i] - gradp[dir][i] / d->rho * d->dt;

    d->delete_var(2, gradp);
}

void projection::update_uf()
{
    for (size_t idir = 0; idir < NDIRS; ++idir)
        for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
        {
            mesh_row *row = d->rows[idir] + irow;
            size_t n = row->n;

            double *ustar = d-> extract_scalars(row, d->ustar[idir]);
            double *p = d-> extract_scalars(row, d->p);
            double *uf = new double[n];

            for (size_t i = 0; i < n - 1; ++i)
                uf[i] = (ustar[i] + ustar[i + 1]) / 2.0 + (p[i + 1] - p[i]) / d->delta;

            d->insert_scalars(row, d->uf[idir], uf);

            delete[] ustar;
            delete[] p;
            delete[] uf;
        }
}

}

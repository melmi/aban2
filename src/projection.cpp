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

    raw_coeffs = new triplet_vector();
    make_matrix();
    delete raw_coeffs;

    psolver = new solver_t(*pmatrix);
}

projection::~projection()
{
    delete pmatrix;
    delete psolver;
}

/*
 * functions to make pressure coeffs matrix
 */

void projection::make_matrix()
{
    for (size_t idir = 0; idir < NDIRS; ++idir)
        for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
        {
            mesh_row *row = d->rows[idir] + irow;
            add_row(row);
        }

    pmatrix->setFromTriplets(raw_coeffs->begin(), raw_coeffs->end());
}

void projection::add_row(mesh_row *row)
{
    size_t n = row->n;
    size_t *row_idxs = d->get_row_idxs(row);
    auto start_bc_type = d->boundaries[row->start_code]->ptype;
    auto end_bc_type = d->boundaries[row->end_code]->ptype;

    for (size_t i = 1; i < n - 1; ++i)
    {
        size_t ix = row_idxs[i];
        raw_coeffs->push_back(triplet_t(ix, row_idxs[i - 1], +1.0 * h2inv));
        raw_coeffs->push_back(triplet_t(ix, ix             , -2.0 * h2inv));
        raw_coeffs->push_back(triplet_t(ix, row_idxs[i + 1], +1.0 * h2inv));
    }

    apply_row_bc(row_idxs[0    ], row_idxs[1    ], start_bc_type);
    apply_row_bc(row_idxs[n - 1], row_idxs[n - 2], end_bc_type  );

    delete[] row_idxs;
}

void projection::apply_row_bc(size_t ix0 , size_t ix1, bctype bct)
{
    if (bct == bctype::neumann)
        raw_coeffs->push_back(triplet_t(ix0, ix0, -1.0 * h2inv));
    else
        raw_coeffs->push_back(triplet_t(ix0, ix0, -3.0 * h2inv));

    raw_coeffs->push_back(triplet_t(ix0, ix1, +1.0 * h2inv));
}

/*
 * functions to make rhs of pressure equation
 */

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
    size_t idx;
    bcondition *bc;
    double dx;

    if (side == bcside::start)
    {
        bc = d->boundaries[row->start_code];
        idx = d->cellno(row, 0);
        dx = -d->delta / 2.0;
    }
    else
    {
        bc = d->boundaries[row->end_code];
        idx = d->cellno(row, row->n - 1);
        dx = +d->delta / 2.0;
    }

    if (bc->ptype == bctype::dirichlet)
    {
        double p = bc->p(d, row, bcside::start, row->dir);
        rhs[idx] -= 2.0 * p * h2inv;
    }
    else
    {
        double rhog = d->rho * d->g.components[row->dir] * dx;
        rhs[idx] -= 2.0 * rhog * h2inv;
    }
}

/*
 * public functions
 */

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
    // for (int i = 0; i < d->n; ++i)
    //     for (int dir = 0; dir < NDIRS; ++dir)
    //         d->u[dir][i] = d->ustar[dir][i] + (- d->gradp[dir][i] / d->rho + d->g.components[dir]) * d->dt;
}

}

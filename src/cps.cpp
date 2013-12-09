/**
 * cps.cpp
 *
 * by Mohammad Elmi
 * Created on 23 Nov 2013
 */

#include "cps.h"
#include "gradient.h"

namespace aban2
{

CPS::CPS(domain *_d): d(_d)
{
    h2inv = 1.0 / (d->delta * d->delta);

    pmatrix = new SparseMatrix(d->n, d->n);
    gmatrix = new SparseMatrix(d->n, d->n);

    raw_coeffs = new TripletVector();
    make_pressure_matrix();
    raw_coeffs->clear();
    make_grad_matrix();
    delete raw_coeffs;

    psolver = new PressureSolver(*pmatrix);
    gsolver = new GradientSolver(*gmatrix);
}

CPS::~CPS()
{
    delete pmatrix;
    delete psolver;

    delete pmatrix;
    delete gmatrix;
}

/*
 * functions to make pressure coeffs matrix
 */

void CPS::make_pressure_matrix()
{
    for (size_t idir = 0; idir < NDIRS; ++idir)
        for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
        {
            mesh_row *row = d->rows[idir] + irow;
            add_pressure_row(row);
        }

    pmatrix->setFromTriplets(raw_coeffs->begin(), raw_coeffs->end());
}

void CPS::add_pressure_row(mesh_row *row)
{
    size_t n = row->n;
    size_t *row_idxs = d->get_row_idxs(row);
    auto start_bc_type = d->boundaries[row->start_code]->ptype;
    auto end_bc_type = d->boundaries[row->end_code]->ptype;

    for (size_t i = 1; i < n - 1; ++i)
    {
        size_t ix = row_idxs[i];
        raw_coeffs->push_back(Triplet(ix, row_idxs[i - 1], +1.0 * h2inv));
        raw_coeffs->push_back(Triplet(ix, ix             , -2.0 * h2inv));
        raw_coeffs->push_back(Triplet(ix, row_idxs[i + 1], +1.0 * h2inv));
    }

    apply_pressure_row_bc(row_idxs[0    ], row_idxs[1    ], start_bc_type);
    apply_pressure_row_bc(row_idxs[n - 1], row_idxs[n - 2], end_bc_type  );

    delete[] row_idxs;
}

void CPS::apply_pressure_row_bc(size_t ix0 , size_t ix1, bctype bct)
{
    if (bct == bctype::neumann)
        raw_coeffs->push_back(Triplet(ix0, ix0, -1.0 * h2inv));
    else
        raw_coeffs->push_back(Triplet(ix0, ix0, -3.0 * h2inv));

    raw_coeffs->push_back(Triplet(ix0, ix1, +1.0 * h2inv));
}

/*
 * functions to make rhs of pressure equation
 */

double *CPS::get_pressure_rhs()
{
    double *rhs = (double *)d->create_var(1);
    gradient::divergance(d, d->ustar, rhs, &bcondition::u);
    double rho_dt = d->rho / d->dt;
    for (int i = 0; i < d->n; ++i) rhs[i] *= rho_dt;
    apply_pressure_rhs_bc(rhs);
    return rhs;
}

void CPS::apply_pressure_rhs_bc(double *rhs)
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;

            apply_single_pressure_rhs_bc(row, rhs, bcside::start);
            apply_single_pressure_rhs_bc(row, rhs, bcside::end);
        }
}

void CPS::apply_single_pressure_rhs_bc(mesh_row *row, double *rhs, bcside side)
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
 * functions to make gradient coeff matrix
 */

void CPS::make_grad_matrix()
{
    for (size_t idir = 0; idir < NDIRS; ++idir)
        for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
        {
            mesh_row *row = d->rows[idir] + irow;
            add_grad_row(row);
        }

    gmatrix->setFromTriplets(raw_coeffs->begin(), raw_coeffs->end());
}

void CPS::add_grad_row(mesh_row *row)
{
    size_t n = row->n;
    size_t *row_idxs = d->get_row_idxs(row);
    auto start_bc_type = d->boundaries[row->start_code]->ptype;
    auto end_bc_type = d->boundaries[row->end_code]->ptype;
    double coeff = 1.0 / (4.0 * NDIRS);

    for (size_t i = 1; i < n - 1; ++i)
    {
        size_t ix = row_idxs[i];
        raw_coeffs->push_back(Triplet(ix, row_idxs[i - 1], 1.0 * coeff));
        raw_coeffs->push_back(Triplet(ix, ix             , 2.0 * coeff));
        raw_coeffs->push_back(Triplet(ix, row_idxs[i + 1], 1.0 * coeff));
    }

    apply_grad_row_bc(row_idxs[0]    , row_idxs[1]    , start_bc_type);
    apply_grad_row_bc(row_idxs[n - 1], row_idxs[n - 2], end_bc_type  );

    delete[] row_idxs;
}

void CPS::apply_grad_row_bc(size_t ix0 , size_t ix1, bctype bct)
{
    /*
     * if (bct == bctype::neumann)  grad p = rho *g else grad p = 0
     * so grad p has a known value for both conditions
     */

    double coeff = 1.0 / (4.0 * NDIRS);

    raw_coeffs->push_back(Triplet(ix0, ix0, +1.0 * coeff));
    raw_coeffs->push_back(Triplet(ix0, ix1, +1.0 * coeff));
}

/*
 * public functions
 */

void CPS::solve_pressure()
{
    double *rhs = get_pressure_rhs();
    Eigen::VectorXd b(d->n), x(d->n);
    for (int i = 0; i < d->n; ++i) b[i] = rhs[i];
    delete[] rhs;
    x = psolver->solve(b);
    for (int i = 0; i < d->n; ++i) d->p[i] = x[i];
}

void CPS::solve_grads()
{
    // double **rhs = (double **)d->create_var(2);
    // gradient::of_scalar(d, d->p, rhs, &bcondition::p);

    // Eigen::VectorXd b(d->n), x(d->n);
    // for (int dir = 0; dir < NDIRS; ++dir)
    // {
    //     double *rhs_dir = rhs[dir];
    //     for (int i = 0; i < d->n; ++i) b[i] = rhs_dir[i];
    //     x = gsolver->solve(b);
    //     double *gradp_dir = d->gradp[dir];
    //     for (int i = 0; i < d->n; ++i) gradp_dir[i] = x[i];
    // }

    // d->delete_var(2, rhs);
}

void CPS::update_u()
{
    // for (int i = 0; i < d->n; ++i)
    //     for (int dir = 0; dir < NDIRS; ++dir)
    //         d->u[dir][i] = d->ustar[dir][i] + (- d->gradp[dir][i] / d->rho + d->g.components[dir]) * d->dt;
}

}
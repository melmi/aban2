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

    pmatrix.Reallocate(d->n, d->n);

    make_matrix();
}

projection::~projection()
{
}

void projection::make_matrix()
{
    for (size_t idir = 0; idir < NDIRS; ++idir)
        for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
        {
            mesh_row *row = d->rows[idir] + irow;
            add_row(row);
        }
}

void projection::add_row(mesh_row *row)
{
    size_t n = row->n;
    size_t *row_idxs = d->get_row_idxs(row);
    auto start_bc_desc = d->boundaries[row->start_code]->p(row, bcside::start, 0);
    auto end_bc_desc   = d->boundaries[row->end_code  ]->p(row, bcside::end  , 0);

    for (size_t i = 1; i < n - 1; ++i)
    {
        size_t ix = row_idxs[i];
        pmatrix.AddInteraction(ix, row_idxs[i - 1], +1.0 * h2inv);
        pmatrix.AddInteraction(ix, ix             , -2.0 * h2inv);
        pmatrix.AddInteraction(ix, row_idxs[i + 1], +1.0 * h2inv);
    }

    apply_row_bc(row_idxs[0    ], row_idxs[1    ], start_bc_desc);
    apply_row_bc(row_idxs[n - 1], row_idxs[n - 2], end_bc_desc  );

    delete[] row_idxs;
}

void projection::apply_row_bc(size_t ix0 , size_t ix1, bcdesc desc)
{
    pmatrix.AddInteraction(ix0, ix0, (2.0 * desc.sw - 3.0)*h2inv);
    pmatrix.AddInteraction(ix0, ix1, +1.0 * h2inv);
}

double *projection::get_rhs()
{
    double *rhs = (double *)d->create_var(1);
    gradient::divergance(d, d->qstar, rhs, &bcondition::q);
    for (int i = 0; i < d->n; ++i) rhs[i] /= d->dt;
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
    bcondition *bc = side == bcside::start ? d->boundaries[row->start_code] : d->boundaries[row->end_code];
    auto desc = bc->p(row, side, row->dir);
    rhs[desc.cellno] -= 2.0 * desc.cte * h2inv;
}

void projection::solve_p()
{
    double *rhs = get_rhs();
    Seldon::Vector<double> b(d->n), x(d->n);
    x.SetData(d->n, d->p);
    b.SetData(d->n, rhs);

    Seldon::Preconditioner_Base precond;

    Seldon::Iteration<double> iter(1000, 1e-3);
    iter.HideMessages();
    iter.SetInitGuess(false);

    int error = Seldon::BiCgStab(pmatrix, x, b, precond, iter);

    std::cout << "                   #iterations: " << iter.GetNumberIteration() << std::endl;
    std::cout << "                   successful:  " << (error == 0) << std::endl;

    x.Nullify();
    b.Nullify();

    delete[] rhs;
}

void projection::update_q()
{
    double **gradp = (double **)d->create_var(2);
    gradient::of_scalar(d, d->p, gradp, &bcondition::p);

    for (int i = 0; i < d->n; ++i)
        for (int dir = 0; dir < NDIRS; ++dir)
            d->q[dir][i] = d->qstar[dir][i] - gradp[dir][i] * d->dt;

    d->delete_var(2, gradp);
}

void projection::update_uf()
{
    for (size_t idir = 0; idir < NDIRS; ++idir)
        for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
        {
            mesh_row *row = d->rows[idir] + irow;

            double *qstar = d->extract_scalars(row, d->qstar[idir]);
            double *p = d->extract_scalars(row, d->p);
            double *uf = new double[row->n];

            for (size_t i = 0; i < row->n - 1; ++i)
                uf[i] = (qstar[i] + qstar[i + 1]) / d->_rho / 2.0 - (p[i + 1] - p[i]) / d->delta / d->_rho;

            d->insert_scalars(row, d->uf[idir], uf);

            delete[] qstar;
            delete[] p;
            delete[] uf;
        }
}

}

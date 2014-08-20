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
    auto desc_start = d->boundaries[row->start_code]->p->desc(d->cellno(row, 0         ), row->dir);
    auto desc_end   = d->boundaries[row->end_code  ]->p->desc(d->cellno(row, row->n - 1), row->dir);

    for (size_t i = 1; i < n - 1; ++i)
    {
        size_t ix = row_idxs[i];
        pmatrix.AddInteraction(ix, row_idxs[i - 1], +1.0 * h2inv);
        pmatrix.AddInteraction(ix, ix             , -2.0 * h2inv);
        pmatrix.AddInteraction(ix, row_idxs[i + 1], +1.0 * h2inv);
    }

    apply_row_bc(row_idxs[0    ], row_idxs[1    ], desc_start);
    apply_row_bc(row_idxs[n - 1], row_idxs[n - 2], desc_end  );

    delete[] row_idxs;
}

void projection::apply_row_bc(size_t ix0 , size_t ix1, bcdesc desc)
{
    pmatrix.AddInteraction(ix0, ix0, (2.0 * desc.sw - 3.0)*h2inv);
    pmatrix.AddInteraction(ix0, ix1, +1.0 * h2inv);
}

double *projection::get_rhs()
{
    auto rhs = gradient::divergance(d, d->ustar, flowbc::bc_u_getter);
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

            size_t cellno_start = d->cellno(row, 0);
            size_t cellno_end   = d->cellno(row, row->n - 1);

            auto bcdesc_start = d->boundaries[row->start_code]->p->desc(cellno_start, row->dir);
            auto bcdesc_end   = d->boundaries[row->end_code  ]->p->desc(cellno_end  , row->dir);

            rhs[cellno_start] -= 2.0 * bcdesc_start.cte * h2inv;
            rhs[cellno_end  ] -= 2.0 * bcdesc_end  .cte * h2inv;
        }
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

void projection::update_u()
{
    double **gradp = gradient::of_scalar(d, d->p, flowbc::bc_p_getter);

    for (int i = 0; i < d->n; ++i)
        for (int dir = 0; dir < NDIRS; ++dir)
            d->u[dir][i] = d->ustar[dir][i] - gradp[dir][i] * d->dt / d->rho[i];

    d->delete_var(2, gradp);
}

void projection::update_uf()
{
    for (size_t idir = 0; idir < NDIRS; ++idir)
        for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
        {
            mesh_row *row = d->rows[idir] + irow;

            double *ustar = d->extract_scalars(row, d->ustar[idir]);
            double *p = d->extract_scalars(row, d->p);
            double *uf = new double[row->n];

            for (size_t i = 0; i < row->n - 1; ++i)
                uf[i] = (ustar[i] + ustar[i + 1]) / 2.0 - (p[i + 1] - p[i]) / d->delta;

            d->insert_scalars(row, d->uf[idir], uf);

            delete[] ustar;
            delete[] p;
            delete[] uf;
        }
}

}

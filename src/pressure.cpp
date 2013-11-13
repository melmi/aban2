/**
 * pressure.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "pressure.h"

namespace aban2
{

void pressure::apply_pressure_row_bc(size_t ix0 , size_t ix1 , size_t ix2 , size_t ix3 , triplet_vector *v, bctype bct, double coeff)
{
    //0
    if (bct == bctype::neumann)
        v->push_back(triplet(ix0, ix0, -2.0 * coeff));
    else
        v->push_back(triplet(ix0, ix0, -8.0 * coeff));
    v->push_back(triplet(ix0, ix1, +1.0 * coeff));
    v->push_back(triplet(ix0, ix2, +1.0 * coeff));

    //1
    if (bct == bctype::neumann)
        v->push_back(triplet(ix1, ix0, +1.0 * coeff));
    else
        v->push_back(triplet(ix1, ix0, -1.0 * coeff));
    v->push_back(triplet(ix1, ix1, -2.0 * coeff));
    v->push_back(triplet(ix1, ix3, +1.0 * coeff));
}

void pressure::add_pressure_row(mesh_row *row, size_t *row_idxs, triplet_vector *v, bctype start_bc_type, bctype end_bc_type, double coeff)
{
    size_t n = row->n;

    for (size_t i = 2; i < n - 2; ++i)
    {
        size_t ix = row_idxs[i];
        v->push_back(triplet(ix, row_idxs[i - 2], 1.0 * coeff));
        v->push_back(triplet(ix, ix, -2.0 * coeff));
        v->push_back(triplet(ix, row_idxs[i + 2], 1.0 * coeff));
    }

    apply_pressure_row_bc(row_idxs[0], row_idxs[1], row_idxs[2], row_idxs[3], v, start_bc_type, coeff);
    apply_pressure_row_bc(row_idxs[n - 1], row_idxs[n - 2], row_idxs[n - 3], row_idxs[n - 4], v, end_bc_type, coeff);
}

pressure::sparse_matrix *pressure::make_pressure_matrix(domain *d)
{
    double coeff = 1.0 / (4.0 * d->delta * d->delta);
    triplet_vector v(d->n * 5);

    for (size_t idir = 0; idir < NDIRS; ++idir)
        for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
        {
            mesh_row *row = d->rows[idir] + irow;
            size_t *row_idxs = d->get_row_idxs(row);
            auto start_bc_type = d->boundaries[row->start_code]->ptype;
            auto end_bc_type = d->boundaries[row->end_code]->ptype;
            add_pressure_row(row, row_idxs, &v, start_bc_type, end_bc_type, coeff);
            delete[] row_idxs;
        }

    sparse_matrix *m = new sparse_matrix(d->n, d->n);
    m->setFromTriplets(v.begin(), v.end());

    return m;
}

}

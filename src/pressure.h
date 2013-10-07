/**
 * pressure.h
 *
 * by Mohammad Elmi
 * Created on 3 Sep 2013
 */

#ifndef _PRESSURE_H_
#define _PRESSURE_H_

#include "domain.h"
#include <vector>
#include <eigen3/Eigen/SparseCore>

namespace aban2
{

class pressure
{
    typedef Eigen::Triplet<double> triplet;
    typedef std::vector<triplet> triplet_vector;

    static void apply_pressure_row_bc(size_t ix0 , size_t ix1 , size_t ix2 , size_t ix3 , triplet_vector *v, bcond bc, double coeff)
    {
        //0
        if (bc.type == bctype::neumann)
            v->push_back(triplet(ix0, ix0, -2.0 * coeff));
        else
            v->push_back(triplet(ix0, ix0, -4.0 * coeff));
        v->push_back(triplet(ix0, ix1, +1.0 * coeff));
        v->push_back(triplet(ix0, ix2, +1.0 * coeff));

        //1
        if (bc.type == bctype::neumann)
            v->push_back(triplet(ix1, ix0, +1.0 * coeff));
        else
            v->push_back(triplet(ix1, ix0, -1.0 * coeff));
        v->push_back(triplet(ix1, ix1, -2.0 * coeff));
        v->push_back(triplet(ix1, ix3, +1.0 * coeff));
    }

    static void add_pressure_row(size_t n, size_t *row_idxs, triplet_vector *v, bcond startbc, bcond endbc, double coeff)
    {
        for (size_t i = 2; i < n - 2; ++i)
        {
            size_t ix = row_idxs[i];
            v->push_back(triplet(ix, row_idxs[i - 2], 1.0 * coeff));
            v->push_back(triplet(ix, ix, -2.0 * coeff));
            v->push_back(triplet(ix, row_idxs[i + 2], 1.0 * coeff));
        }

        apply_pressure_row_bc(row_idxs[0], row_idxs[1], row_idxs[2], row_idxs[3], v, startbc, coeff);
        apply_pressure_row_bc(row_idxs[n - 1], row_idxs[n - 2], row_idxs[n - 3], row_idxs[n - 4], v, endbc, coeff);
    }

public:
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, long> sparse_matrix;

    sparse_matrix *make_pressure_matrix(domain *d)
    {
        double coeff = 1.0 / (4.0 * d->delta * d->delta);
        triplet_vector v(d->n * 5);

        for (size_t idir = 0; idir < 3; ++idir)
            for (size_t irow = 0; irow < d->nrows[idir]; ++irow)
            {
                mesh_row *row = d->rows[idir] + irow;
                size_t *row_idxs = d->get_row_idxs(row);
                bcond start_bc = d->boundaries[row->start_bc].pbc;
                bcond end_bc = d->boundaries[row->end_bc].pbc;
                add_pressure_row(row->n, row_idxs, &v, start_bc, end_bc, coeff);
                delete[] row_idxs;
            }

        sparse_matrix *m = new sparse_matrix(d->n, d->n);
        m->setFromTriplets(v.begin(), v.end());
        return m;
    }
};

}

#endif
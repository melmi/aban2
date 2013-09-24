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

    static void add_pressure_row(size_t n, size_t *row_idxs, double dx, triplet_vector &v)
    {
        for (size_t i = 2; i < n - 2; ++i)
        {
            size_t ix = row_idxs[i];
            v.push_back(triplet(ix, row_idxs[i - 2], 1.0));
            v.push_back(triplet(ix, ix, 2.0));
            v.push_back(triplet(ix, row_idxs[i + 2], 1.0));
        }
    }

public:
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, long> sparse_matrix;

    sparse_matrix *make_pressure_matrix(domain d)
    {
        triplet_vector v(d.n * 5);

        for (size_t idir = 0; idir < 3; ++idir)
            for (size_t irow = 0; irow < d.nrows[idir]; ++irow)
            {
                mesh_row *row = d.rows[idir] + irow;
                size_t *row_idxs = d.get_row_idxs(*row);
                add_pressure_row(row->n, row_idxs, d.delta, v);
                delete[] row_idxs;
            }

        sparse_matrix * m = new sparse_matrix(d.n, d.n);
        m->setFromTriplets(v.begin(), v.end());
        return m;
    }
};

}

#endif
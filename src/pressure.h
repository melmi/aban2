/**
 * pressure.h
 *
 * by Mohammad Elmi
 * Created on 3 Sep 2013
 */

#ifndef _PRESSURE_H_
#define _PRESSURE_H_

#include <vector>
#include <eigen3/Eigen/SparseCore>

#include "domain.h"
#include "bcondition.h"

namespace aban2
{

class pressure
{
    typedef Eigen::Triplet<double> triplet;
    typedef std::vector<triplet> triplet_vector;

    static void apply_pressure_row_bc(
        size_t ix0, size_t ix1, size_t ix2, size_t ix3,
        triplet_vector *v, bctype bct, double coeff);

    static void add_pressure_row(
        mesh_row *row,
        size_t *row_idxs, triplet_vector *v,
        bctype start_bc_type, bctype end_bc_type, 
        double coeff);

public:
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, long> sparse_matrix;

    static sparse_matrix *make_pressure_matrix(domain *d);
};

}

#endif
/**
 * projection.h
 *
 * by Mohammad Elmi
 * Created on 5 Oct 2013
 */

#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include <algorithm>
#include <eigen3/Eigen/IterativeLinearSolvers>

#include "domain.h"
#include "pressure.h"

namespace aban2
{

class projection
{
    static double *get_div_ustar(domain *d);

    static void apply_p_bc(domain *d, double *rhs);

    static double *get_pressure_rhs(domain *d);

public:
    typedef Eigen::BiCGSTAB<pressure::sparse_matrix> psolver;

    static void solve_p(domain *d, psolver *solver);

    static void update_u(domain *d);

};

}

#endif
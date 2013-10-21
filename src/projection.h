/**
 * projection.h
 *
 * by Mohammad Elmi
 * Created on 5 Oct 2013
 */

#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include <algorithm>

#include "domain.h"
#include "pressure.h"

namespace aban2
{

class projection
{
    static double *get_div_ustar(domain *d)
    {
        double *result = (double **)d->create_var(1);
        gradient::divergance(d, d->ustar, result, flow_boundary::velbc);
        double rho_dt = d->rho / d->dt;
        for (int i = 0; i < d->n; ++i) result[i] *= rho_dt;
        apply_dirichlet_p_bc_on_laplac_p_rhs(d, result);
        return result;
    }

    static void apply_dirichlet_p_bc_on_laplac_p_rhs(domain *d, double *rhs)
    {
        bcond bc;

        for (size_t dir = 0; dir < NDIRS; ++dir)
            for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            {
                mesh_row *row = d->rows[dir] + irow;

                bc = d->boundaries[row->startbc].pbc;
                if (bc.type == bctype::dirichlet)
                    rhd[row->start[row->ii]] += -bc.value / d->dx / d->dx;

                bc = d->boundaries[row->endbc].pbc;
                if (bc.type == bctype::dirichlet)
                    rhd[row->start[row->ii]] += -bc.value / d->dx / d->dx;
            }
    }

    static double *get_pressure_rhs(domain *d)
    {
        double *prhs = get_div_ustar(d);
        apply_dirichlet_p_bc_on_laplac_p_rhs(d, rhs);
        return rhs;
    }

public:

    void solve_p(domain *d, pressure::sparse_matrix *m)
    {
        double *rhs = get_pressure_rhs(d);
        //m->
        delete[] rhs;
    }
};

}

#endif
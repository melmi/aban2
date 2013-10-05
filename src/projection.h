/**
 * projection.h
 *
 * by Mohammad Elmi
 * Created on 5 Oct 2013
 */

#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include "domain.h"
#include <algorithm>

namespace aban2
{

class projection
{
    static double *get_div_ustar(domain *d)
    {
        double *result = (double **)d->create_var(1);
        gradient::divergance(d, d->ustar, result, flow_boundary::velbc);
        return result;
    }

    static void apply_dirichlet_p_bc_on_div_ustar(domain* d, double* div_ustar)
    {

    }

    static double* get_pressure_rhs(domain* d)
    {
        
    }

public:
};

}

#endif
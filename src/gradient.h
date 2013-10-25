/**
 * gradient.h
 *
 * by Mohammad Elmi
 * Created on 5 Oct 2013
 */

#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include "domain.h"
#include <algorithm>

namespace aban2
{

class gradient
{
public:
    static void add_1d_row(size_t n, double *phi, double *grad, double dx, bcond startbc, bcond endbc);

    static double *get_1d_row(size_t n, double *phi, double dx, bcond startbc, bcond endbc);

    static void of_scalar(domain *d, double *phi, double **grad, bcond flow_boundary::*bc);

    static void of_vec(domain *d, double **phi, double ***grad, bcond * flow_boundary::*bc);

    static void divergance(domain *d, double **phi, double *divergance, bcond (flow_boundary::*bc)[3]);
};

}

#endif

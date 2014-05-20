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
#include "bcondition.h"

namespace aban2
{

class gradient
{
public:
    static void add_1d_row(domain *d, mesh_row *row, double *phi, double *grad_dir, flowbc::member bc);

    static double *of_scalar_dir_oriented(domain *d, double *phi, flowbc::member bc, size_t dir);

    static double **of_scalar(domain *d, double *phi, flowbc::member bc);

    static double ** *of_vec(domain *d, double **phi, flowbc::member bc[3]);

    static double *divergance(domain *d, double **phi, flowbc::member bc[3]);
};

}

#endif

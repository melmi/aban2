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
    static void add_1d_row(domain *d, mesh_row *row, double *phi, double *grad, bcondition::func bcfunc, size_t cmpnt);

    static double *get_1d_row(domain *d, mesh_row *row, double *phi, bcondition::func bcfunc, size_t cmpnt);

    static void of_scalar(domain *d, double *phi, double **grad, bcondition::func bcfunc);

    static void of_vec(domain *d, double **phi, double ***grad, bcondition::func bcfunc);

    static void divergance(domain *d, double **phi, double *divergance, bcondition::func bcfunc);
};

}

#endif

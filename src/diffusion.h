/**
 * diffusion.h
 *
 * by Mohammad Elmi
 * Created on 5 Sep 2013
 */

#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include <algorithm>

#include "domain.h"
#include "bcondition.h"
 
namespace aban2
{

class diffusion
{
public:
    static void solve_tridiagonal_in_place_destructive(double *x, const size_t N, const double *a, const double *b, double *c);

    static void diffuse(domain *d, mesh_row *row, double *phi, double D, bcondition::type bctype, bcondition::func bcfunc, size_t cmpnt);

    static void diffuse_ustar(domain *d);
};

}

#endif
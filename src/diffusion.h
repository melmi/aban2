/**
 * diffusion.h
 *
 * by Mohammad Elmi
 * Created on 5 Sep 2013
 */

#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include "domain.h"
#include <algorithm>

namespace aban2
{

class diffusion
{
public:
    static void solve_tridiagonal_in_place_destructive(double *x, const size_t N, const double *a, const double *b, double *c);
    static void diffuse(size_t n, double *phi, double d, double dt, double dx, bcond startbc, bcond endbc);
    static void diffuse_ustar(domain *d);
};

}

#endif
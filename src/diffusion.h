/**
 * diffusion.h
 *
 * by Mohammad Elmi
 * Created on 5 Sep 2013
 */

#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include "common.h"
#include <algorithm>

#include "domain.h"
#include "bcondition.h"
 
namespace aban2
{

class diffusion
{
private:
    domain* d;

    static void solve_tridiagonal_in_place_destructive(double *x, const size_t N, const double *a, const double *b, double *c);

public:
    diffusion(domain *_d): d(_d) {}
    void diffuse(mesh::row *row, double *phi, double* D, flowbc::member mem);
    void diffuse_ustar();
};

}

#endif

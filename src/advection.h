/**
 * advection.h
 *
 * by Mohammad Elmi
 * Created on 4 Sep 2013
 */

#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include "domain.h"
#include "gradient.h"
#include <algorithm>

namespace aban2
{

class advection
{
public:
    static void advect(size_t n, double *phi, double *u, double dt, double dx, bcond startbc, bcond endbc);
    static void advect_ustar(domain *d);
};

}

#endif
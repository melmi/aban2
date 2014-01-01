/**
 * advection.h
 *
 * by Mohammad Elmi
 * Created on 4 Sep 2013
 */

#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include "bcondition.h"

namespace aban2
{

class domain;
class mesh_row;

class advection
{
public:
    static void advect(domain *d, mesh_row *row, double *phi, double *u, bcondition::func bcfunc, size_t cmpnt);

    static void advect_ustar(domain *d);
};

}

#endif
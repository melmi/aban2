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
    domain *d;

public:
    advection(domain *_d): d(_d) {}

    void advect(mesh_row *row, double *phi, bcondition::func bcfunc, size_t cmpnt);

    void advect_qstar();
};

}

#endif
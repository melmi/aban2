/**
 * advection.h
 *
 * by Mohammad Elmi
 * Created on 4 Sep 2013
 */

#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include "common.h"
#include "bcondition.h"
#include "domain.h"

namespace aban2
{

class advection
{
    domain_t *d;

public:
    advection(domain_t *_d): d(_d) {}

    void advect(row_t *row, double *phi, double *grad_dir, flowbc::bc_val_getter bc);

    void advect_ustar();
};

}

#endif

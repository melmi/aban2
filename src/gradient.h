/**
 * gradient.h
 *
 * by Mohammad Elmi
 * Created on 5 Oct 2013
 */

#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include "common.h"
#include "domain.h"
#include <algorithm>
#include "bcondition.h"

namespace aban2
{

class gradient
{
    static void add_1d_row(domain_t *d, row_t *row, double *phi, double *grad_dir, flowbc::bc_val_getter bc);
public:
    static double *of_scalar_dir(domain_t *d, double *phi, flowbc::bc_val_getter bc, size_t dir);
    static double **of_scalar(domain_t *d, double *phi, flowbc::bc_val_getter bc);
    static double ** *of_vec(domain_t *d, double **phi, flowbc::bc_val_getter bc[3]);
    static double ** *of_uf(domain_t *d);
    static double *of_uf_dir(domain_t *d, size_t dir);
    static double *divergance(domain_t *d);
    static double *divergance_of(domain_t *d, double **phi, flowbc::bc_val_getter bc[3]);
};

}

#endif

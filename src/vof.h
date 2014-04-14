/**
 * vof.h
 *
 * by Mohammad Elmi
 * Created on 10 Nov 2013
 */

#ifndef _VOF_H_
#define _VOF_H_

#include "domain.h"
#include "diffusion.h"

namespace aban2
{

class voset;

class vof
{
    domain *d;
    voset *vos;
    double *vof_mass;
    double h3;
    bool xfirst = false;
    double get_vof_flux(mesh_row *row, size_t i, double udt);
    void reconstruct_full_and_empty_cells();
    void advect_row(mesh_row *row);
public:
    vof(domain *_d);
    ~vof();

    void advect();
};

}

#endif
/**
 * solver.h
 *
 * by Mohammad Elmi
 * Created on 3 Sep 2013
 */

#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <algorithm>
#include <sstream>
 
#include "domain.h"
#include "advection.h"
#include "diffusion.h"
#include "pressure.h"
#include "projection.h"
#include "bodyforce.h"

namespace aban2
{

class solver
{
public:
    domain *d;
    pressure::sparse_matrix *pmatrix;
    projection::psolver *psolver;
    string out_path;

    void write_step(size_t step)
    {
        stringstream s;
        s << out_path << step << ".vtk";
        d->write_vtk(s.str());
    }

    solver(domain *_d, string _out_path): d(_d), out_path(_out_path)
    {
        pmatrix = pressure::make_pressure_matrix(d);
        psolver = new projection::psolver(*pmatrix);
    }

    void step()
    {
        write_step(0);

        for (int i = 0; i < NDIRS; ++i)
            std::copy_n(d->u[i], d->n, d->ustar[i]);

        write_step(1);
        advection::advect_ustar(d);
        write_step(2);
        diffusion::diffuse_ustar(d);
        write_step(3);
        bodyforce::add_g(d);
        write_step(4);
        projection::solve_p(d, psolver);
        write_step(5);
        projection::update_u(d);
        write_step(6);
    }
};

}

#endif
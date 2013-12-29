/**
 * solver.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "solver.h"

#include <iostream>
#include <Eigen32/Eigenvalues>

namespace aban2
{

void solver::write_step(size_t step)
{
    std::stringstream s;
    s << out_path << step << ".vtk";
    d->write_vtk(s.str());
}

solver::solver(domain *_d, std::string _out_path): d(_d), out_path(_out_path)
{
    projector = new projection(d);
}

solver::~solver()
{
    delete projector;
}

void solver::step()
{
    for (int i = 0; i < NDIRS; ++i)
        std::copy_n(d->u[i], d->n, d->ustar[i]);

    advection::advect_ustar(d);
    diffusion::diffuse_ustar(d);
    projector->solve_p();
    projector->update_u();

    d->t += d->dt;
}

void solver::run(double tend)
{
    write_step(0);
    size_t nsteps = (tend - d->t) / d->dt;

    for (int it = 0; it < nsteps; ++it)
    {
        std::cout << "running step " << it + 1
                  // << " [ iterations: " << psolver->iterations() << "  error: " << psolver->error() << " ]"
                  << std::endl << std::flush;
        step();

        if ((it + 1) % d->step_write == 0)
            write_step((it + 1) / d->step_write);
    }
}

}

/**
 * solver.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "solver.h"

#include <iostream>

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
    advector = new advection(d);
    diffusor = new diffusion(d);
}

solver::~solver()
{
    delete projector;
    delete advector;
    delete diffusor;
}

void solver::step()
{
    for (int i = 0; i < NDIRS; ++i)
        std::copy_n(d->q[i], d->n, d->qstar[i]);

    std::cout << "      advection " << std::endl;
    advector->advect_qstar();
    std::cout << "      diffusion " << std::endl;
    diffusor->diffuse_qstar();
    std::cout << "      source terms " << std::endl;
    apply_source_terms();
    std::cout << "      pressure " << std::endl;
    projector->solve_p();
    std::cout << "      updating " << std::endl;
    projector->update_q();
    projector->update_uf();

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

void solver::apply_source_terms()
{
    for (size_t i = 0; i < d->n; ++i)
    {
        d->qstar[0][i] += d->g.components[0] * d->_rho * d->dt;
        d->qstar[1][i] += d->g.components[1] * d->_rho * d->dt;
        d->qstar[2][i] += d->g.components[2] * d->_rho * d->dt;
    }
}

}

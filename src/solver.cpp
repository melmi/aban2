/**
 * solver.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "solver.h"
#include "gradient.h"

#include <iostream>
#include <cmath>
#include <iomanip>

namespace aban2
{

void solver::write_step(size_t step)
{
    std::stringstream s;
    s << out_path
      << std::setfill('0') << std::setw(6) << step
      << ".vtk";
    d->write_vtk(s.str());
}

solver::solver(domain *_d, std::string _out_path): d(_d), out_path(_out_path)
{
    projector = new projection(d);
    _vof = new vof(d);
    diffusor = new diffusion(d);
}

solver::~solver()
{
    delete projector;
    delete _vof;
    delete diffusor;
}

void solver::step()
{
    for (int i = 0; i < NDIRS; ++i)
        std::copy_n(d->u[i], d->n, d->ustar[i]);

    std::cout << "      advection " << std::endl;
    _vof->advect();
    for (size_t i = 0; i < d->n; ++i)
        d->nu[i] = d->nu_bar(d->vof[i]);
    std::cout << "      diffusion " << std::endl;
    diffusor->diffuse_ustar();
    std::cout << "      source terms " << std::endl;
    apply_source_terms();
    std::cout << "      pressure " << std::endl;
    projector->solve_p();
    std::cout << "      updating " << std::endl;
    projector->update_u();
    std::cout << "      calculating uf " << std::endl;
    projector->update_uf();
    std::cout //<< std::scientific
            << "      divergance: " << divergance()
            << std::endl;

    d->t += d->dt;
}

double solver::divergance()
{
    auto grad_uf = gradient::of_uf(d);

    double result = 0;
    for (size_t i = 0; i < d->n; ++i)
    #ifdef THREE_D
        result += std::abs(grad_uf[0][0][i] + grad_uf[1][1][i] + grad_uf[2][2][i]);
    #else
        result += std::abs(grad_uf[0][0][i] + grad_uf[1][1][i]);
    #endif

    domain::delete_var(3, grad_uf);

    return result;
}

void solver::run(size_t nsteps)
{
    for (size_t i = 0; i < d->n; ++i) d->rho[i] = d->rho_bar(d->vof[i]);

    write_step(0);

    for (size_t it = 0; it < nsteps; ++it)
    {
        std::cout << "running step " << it + 1
                  << std::endl << std::flush;
        step();

        if ((it + 1) % d->write_interval == 0)
            write_step((it + 1) / d->write_interval);
    }
}

void solver::apply_source_terms()
{
    for (size_t i = 0; i < d->n; ++i)
    {
        d->ustar[0][i] += d->g.cmpnt[0] * d->dt;
        d->ustar[1][i] += d->g.cmpnt[1] * d->dt;
        d->ustar[2][i] += d->g.cmpnt[2] * d->dt;
    }
}

}

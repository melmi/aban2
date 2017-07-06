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
    domain->write_vtk(s.str());
}

solver::solver(domain_t *_d, std::string _out_path): domain(_d), out_path(_out_path)
{
    projector = new projection(domain);
    vof = new vof_t(domain);
    diffusor = new diffusion(domain);
}

solver::~solver()
{
    delete projector;
    delete vof;
    delete diffusor;
}

void solver::step()
{
    for (int i = 0; i < NDIRS; ++i)
        std::copy_n(domain->u[i], domain->n, domain->ustar[i]);

    std::cout << "      advection " << std::endl;
    vof->advect();
    std::cout << "      diffusion " << std::endl;
    for (size_t i = 0; i < domain->n; ++i)
        domain->nu[i] = domain->nu_bar(domain->vof[i], domain->rho[i]);
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

    domain->t += domain->dt;
}

double solver::divergance()
{
    auto grad_uf = gradient::of_uf(domain);

    double result = 0;
    for (size_t i = 0; i < domain->n; ++i)
    #ifdef THREE_D
        result += std::abs(grad_uf[0][0][i] + grad_uf[1][1][i] + grad_uf[2][2][i]);
    #else
        result += std::abs(grad_uf[0][0][i] + grad_uf[1][1][i]);
    #endif

    domain_t::delete_var(3, grad_uf);

    return result;
}

void solver::run(size_t nsteps)
{
    for (size_t i = 0; i < domain->n; ++i) domain->rho[i] = domain->rho_bar(domain->vof[i]);

    write_step(0);

    for (size_t it = 0; it < nsteps; ++it)
    {
        std::cout << "running step " << it + 1
                  << std::endl << std::flush;
        step();

        if ((it + 1) % domain->write_interval == 0)
            write_step((it + 1) / domain->write_interval);
    }
}

void solver::apply_source_terms()
{
    for (size_t i = 0; i < domain->n; ++i)
    {
        domain->ustar[0][i] += domain->g.cmpnt[0] * domain->dt;
        domain->ustar[1][i] += domain->g.cmpnt[1] * domain->dt;
        domain->ustar[2][i] += domain->g.cmpnt[2] * domain->dt;
    }
}

}

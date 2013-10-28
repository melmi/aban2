/**
 * solver.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "solver.h"

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
    pmatrix = pressure::make_pressure_matrix(d);
    psolver = new projection::psolver(*pmatrix);
}

void solver::step()
{
    write_step(0);

    for (int i = 0; i < NDIRS; ++i)
        std::copy_n(d->u[i], d->n, d->ustar[i]);

    write_step(1);
    advection::advect_ustar(d);
    write_step(2);
    diffusion::diffuse_ustar(d);
    write_step(3);
    projection::solve_p(d, psolver);
    write_step(4);
    projection::update_u(d);
    write_step(5);
}

}

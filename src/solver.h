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

namespace aban2
{

class solver
{
public:
    domain *d;
    pressure::sparse_matrix *pmatrix;
    projection::psolver *psolver;
    std::string out_path;

    void write_step(size_t step);

    solver(domain *_d, std::string _out_path);

    void step();

    void run(double tend);

};

}

#endif
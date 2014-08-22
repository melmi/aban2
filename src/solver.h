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

#include "common.h"
#include "domain.h"
#include "vof.h"
#include "diffusion.h"
#include "projection.h"

namespace aban2
{

class solver
{
private:
    domain *d;
    projection *projector;
    vof *_vof;
    diffusion *diffusor;
    std::string out_path;

    void apply_source_terms();
    double divergance();
public:
    solver(domain *_d, std::string _out_path);
    ~solver();

    void write_step(size_t step);
    void step();
    void run(double tend);
};

}

#endif
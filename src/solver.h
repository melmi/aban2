/**
 * solver.h
 *
 * by Mohammad Elmi
 * Created on 3 Sep 2013
 */

#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <jsoncpp/json/json.h>

#include "domain.h"
#include "advection.h"
#include "diffusion.h"
#include "pressure.h"
#include "projection.h"
 
namespace aban2
{

class solver: public domain
{
public:
    // for (size_t i = 0; i < 3; ++i)
    //     std::copy_n(d.u[i], d.n, d.ustar[i]);
};

}

#endif
#include <algorithm>
#include <iostream>
#include <sstream>
#include <jsoncpp/json/json.h>
#include "mesh.h"
#include "domain.h"
#include "advection.h"
#include "diffusion.h"
#include "projection.h"
#include "tests.h"
#include "solver.h"
#include "vof.h"
#include "vector.h"

using namespace std;
using namespace aban2;

int main(int argc, char const *argv[])
{
    domain *d = domain::create_from_file("mesh/cavity100x100.json");

    std::fill_n(d->vof, d->n, 1.0);
    std::fill_n(d->rho, d->n, d->rho1);
    std::fill_n(d->nu, d->n, d->nu1);
    std::cout << "before run" << std::endl << std::flush;
    solver s(d, "out/cav");
    s.run(10000);

    delete d;
    
    return 0;
}
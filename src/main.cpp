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
    // for (int i = 0; i < d->n; ++i)d->p[i] = i;
    // zalesak_disk_2d(d);
    square_2d(d);
    vof vofc(d);
    for (int i = 0; i < 650; ++i)
    {
        std::cout << "time step: " << i << std::endl;
        d->write_vtk("out/vof" + to_string(i) + ".vtk");
        vofc.advect();
    }

    delete d;
    return 0;
}
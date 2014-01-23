#include <algorithm>
#include <iostream>
#include <jsoncpp/json/json.h>
#include "mesh.h"
#include "domain.h"
#include "advection.h"
#include "diffusion.h"
#include "projection.h"
#include "tests.h"
#include "solver.h"

using namespace std;
using namespace aban2;

int main(int argc, char const *argv[])
{
    cout << "Reading mesh" << endl;
    domain *d = domain::create_from_file("mesh/cavity100x100.json");
    zalesak_disk_2d(d);
    d->write_vtk("out/test.vtk");
    // solver s(d, "out/out");
    // s.step();
    // s.run(10);

    return 0;
}
#include <iostream>
#include <jsoncpp/json/json.h>
#include "mesh.h"
#include "domain.h"
#include "advection.h"
#include "diffusion.h"
#include "pressure.h"
#include "tests.h"

using namespace std;
using namespace aban2;

int main(int argc, char const *argv[])
{
    domain *d = domain::create_from_file("mesh/slope.json");
    cout << NDIRS << endl;

    cout << d->nrows[0] << "    " << d->rows[0][0].start[0] << endl
         << d->nrows[1] << "    " << d->rows[1][0].start[0] << endl;

     adv_test1(d);

     diffusion::diffuse_ustar(d);

     d->write_vtk("/home/mohammad/Desktop/x.vtk");

    return 0;
}
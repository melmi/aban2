#include <iostream>
#include <jsoncpp/json/json.h>
#include "mesh.h"
#include "domain.h"
#include "advection.h"
#include "diffusion.h"

using namespace std;
using namespace aban2;

int main(int argc, char const *argv[])
{
    domain *d = domain::create_from_file("mesh/cavity10x10.json");

    cout << d->nrows[0] << "    " << d->rows[0][0].start[0] << endl
         << d->nrows[1] << "    " << d->rows[1][0].start[1] << endl;

     d->write_vtk("/home/mohammad/Desktop/x.vtk");

    return 0;
}
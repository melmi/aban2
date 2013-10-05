#include <iostream>
#include <sstream>
#include <jsoncpp/json/json.h>
#include "mesh.h"
#include "domain.h"
#include "advection.h"
#include "diffusion.h"
#include "pressure.h"
#include "tests.h"

using namespace std;
using namespace aban2;

string get_fname(string path, size_t step)
{
    stringstream s;
    s << path << step << ".vtk";
    return s.str();
}

int main(int argc, char const *argv[])
{
    cout << "Reading mesh" << endl;
    domain *d = domain::create_from_file("mesh/cavity100x100.json");

    cout << "Initializing domain" << endl;
    // diff_test1(d);
    vortex(d, 1, 51, 51);
    //double *dd = (double *)d->create_var(1);
    gradient::divergance(d, d->u, d->p, &flow_boundary::velbc);
    auto r = std::minmax_element(d->p, d->p + d->n);
    cout << "min: " << *r.first << endl << "max: " << *r.second << endl;
    // d->delete_var(1, dd);
    d->write_vtk(get_fname( "out/out", 0));

    // d->write_vtk(get_fname( "out/out", 0));
    // for (int it = 0; it < 100; ++it)
    // {
    //     cout << "step " << it << endl;
    //     diffusion::diffuse_ustar(d);
    //     d->write_vtk(get_fname( "out/out", it + 1));
    // }

    return 0;
}
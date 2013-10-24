#include <iostream>
#include <jsoncpp/json/json.h>
#include "mesh.h"
#include "domain.h"
#include "advection.h"
#include "diffusion.h"
#include "pressure.h"
#include "projection.h"
#include "tests.h"
#include "solver.h"

using namespace std;
using namespace aban2;

int main(int argc, char const *argv[])
{
    cout << "Reading mesh" << endl;
    domain *d = domain::create_from_file("mesh/cavity4x4.json");
    cout << "mesh has " << d->n << " inner nodes" << endl;
    // print_cell_nos(d);
    // print_rows(d);
    // print_cell_nos(d);

    // cout << "Initializing domain" << endl;
    // // diff_test1(d);
    // vortex(d, 1, 51, 51);
    // //double *dd = (double *)d->create_var(1);
    // gradient::divergance(d, d->u, d->p, &flow_boundary::velbc);
    // auto r = std::minmax_element(d->p, d->p + d->n);
    // cout << "min: " << *r.first << endl << "max: " << *r.second << endl;
    // // d->delete_var(1, dd);
    // d->write_vtk(get_fname( "out/out", 0));

    // // d->write_vtk(get_fname( "out/out", 0));
    // // for (int it = 0; it < 100; ++it)
    // // {
    // //     cout << "step " << it << endl;
    // //     diffusion::diffuse_ustar(d);
    // //     d->write_vtk(get_fname( "out/out", it + 1));
    // // }

    solver s(d, "out/out");
    s.step();


    // print_rows(d);
    return 0;
}
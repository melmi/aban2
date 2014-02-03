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
#include "fsm.h"
#include "vector.h"

using namespace std;
using namespace aban2;

template <typename T>
string to_string (T number)
{
    ostringstream ss;
    ss << number;
    return ss.str();
}

int main(int argc, char const *argv[])
{
    cout << "Reading mesh" << endl;
    domain *d = domain::create_from_file("mesh/cavity100x100.json");
    square_2d(d);
    vof vofc(d);
    d->write_vtk("out/test0.vtk");
    for (int i = 0; i < 628 / d->dt; ++i)
    {
        cout << "timestep " << i + 1 << endl;
        vofc.advect();
        d->write_vtk("out/test" + to_string(i + 1) + ".vtk");
    }
    cout << "done" << endl;

    // solver s(d, "out/out");
    // s.step();
    // s.run(10);

    // cout << endl;
    // aban2::vector n {1, 1, 0};
    // n.normalize();
    // ireconst aa = ireconst::from_alpha({0, 1, 1}, n, 1.414);
    // cout << aa.alpha_max << endl;
    // cout << aa.volume << endl;

    // cout << endl;
    // ireconst a = ireconst::from_volume({3, 1, 1}, n, 1.5);
    // cout << "alpha: " << a.alpha << endl;
    // cout << "n: " << n.x << ", " << n.y << ", " << n.z << endl;
    // cout << endl;
    // cout << a.get_flux(0, 1, n) << endl;
    // cout << a.get_flux(0, 2, n) << endl;
    // cout << a.get_flux(0, 3, n) << endl;
    // cout << endl;
    // cout << a.get_flux(0, -1, n) << endl;
    // cout << a.get_flux(0, -2, n) << endl;
    // cout << a.get_flux(0, -3, n) << endl;

    // ireconst a;
    // a = ireconst::from_volume({1, 1, 1}, n, 0.00); cout << a.alpha << endl;
    // a = ireconst::from_volume({1, 1, 1}, n, 0.25); cout << a.alpha << endl;
    // a = ireconst::from_volume({1, 1, 1}, n, 0.50); cout << a.alpha << endl;
    // a = ireconst::from_volume({1, 1, 1}, n, 1.00); cout << a.alpha << endl;
    // cout << endl;
    // a = ireconst::from_alpha({1, 1, 1}, n, 1.414 * 0.00); cout << a.volume << endl;
    // a = ireconst::from_alpha({1, 1, 1}, n, 1.414 * 0.25); cout << a.volume << endl;
    // a = ireconst::from_alpha({1, 1, 1}, n, 1.414 * 0.50); cout << a.volume << endl;
    // a = ireconst::from_alpha({1, 1, 1}, n, 1.414 * 1.00); cout << a.volume << endl;

    return 0;
}
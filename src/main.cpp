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
    domain *d = domain::create_from_file("mesh/cavity100x100.json");

    diff_test1(d);

    d->write_vtk(get_fname( "out/out", 0));
    for (int it = 0; it < 100; ++it)
    {
        cout << "step " << it << endl;
        diffusion::diffuse_ustar(d);
        d->write_vtk(get_fname( "out/out", it + 1));
    }

    return 0;
}
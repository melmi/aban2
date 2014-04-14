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

void check_continuety(domain *d)
{
    double *div = new double[d->n];
    std::fill_n(div, d->n, 0);

    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 1; irow < d->nrows[dir] - 1; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            for (size_t i = 1; i < row->n - 1; ++i)
            {
                double *u = d->extract_scalars(row, d->uf[row->dir]);

                size_t no = d->cellno(row, i);
                div[no] = u[i] - u[i - 1];

                delete[] u;
            }
        }

    std::for_each(div, div + d->n, [](double & x)
    {
        x = std::abs(x);
    });
    cout << "divergance: " << std::accumulate(div, div + d->n, 0) << std::endl;
}

int main(int argc, char const *argv[])
{
    cout << "Reading mesh" << endl;
    domain *d = domain::create_from_file("mesh/cavity100x100.json");
    square_2d(d);
    check_continuety(d);
    vof vofc(d);
    d->write_vtk("out/test0.vtk");
    for (int i = 0; i < 628 / d->dt; ++i)
    {
        cout << "timestep " << i + 1 << endl;
        vofc.advect();
        d->write_vtk("out/test" + to_string(i + 1) + ".vtk");
    }
    // d->write_vtk("out/test1.vtk");
    cout << "done" << endl;

    // solver s(d, "out/out");
    // s.step();
    // s.run(10);

    return 0;
}
#include "advection.h"
#include "diffusion.h"
#include "domain.h"
#include "projection.h"
#include "solver.h"
#include "tests.h"
#include "vector.h"
#include "vof.h"
#include <algorithm>
#include <iostream>
#include <jsoncpp/json/json.h>
#include <sstream>

using namespace std;
using namespace aban2;

int main(int argc, char const *argv[])
{
    domain_t *d = domain_t::create_from_file("mesh/solitary.json");
    solitary(d, 1, 10, 70, 10);
    for (size_t i = 0; i < d->n; ++i)
    {
        d->rho[i] = d->rho_bar(d->vof[i]);
        d->nu[i] = d->nu_bar(d->vof[i], d->rho[i]);
    }

    std::cout << "before run" << std::endl << std::flush;
    // d->write_vtk("out/solitary0.vtk");
    solver *s = new solver(d, "out/sol");
    s->run(3000);
    delete s;
    delete d;

    // dambreak
    // domain_t *d = domain_t::create_from_file("mesh/dambreak.json");
    // std::fill_n(d->vof, d->n, 0.0);
    // // std::iota(d->no, d->no + d->n, 0);
    // rectangel_by_bounds(d, {0, 0, 0}, {0.05715, 0.05715, 0}, 1);
    // for (size_t i = 0; i < d->n; ++i)
    // {
    //     d->rho[i] = d->rho_bar(d->vof[i]);
    //     d->nu[i] = d->nu_bar(d->vof[i], d->rho[i]);
    // }

    // std::cout << "before run" << std::endl << std::flush;
    // solver *s = new solver(d, "out/dmbrk");
    // s->run(3000);

    // delete s;
    // delete d;

    // consistent_vof_test();
    // zalesak_disk_rotation_test();

    // if (false)
    // {
    //     auto v = aban2::volreconst::from_alpha({2, 2, 2}, {1, 0, 0}, 0.6);
    //     std::cout << "{" << std::endl
    //               << "\tc: " << v->c << "," << std::endl
    //               << "\tm: " << v->m << "," << std::endl
    //               << "\tvolume: { expected: " << 2.4 << ", calculated: " << v->volume << " }" << std::endl
    //               << "}" << std::endl;
    // }

    return 0;
}

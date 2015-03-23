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
    domain *d = domain::create_from_file("mesh/dambreak.json");
    std::fill_n(d->vof, d->n, 0.0);
    // std::iota(d->no, d->no + d->n, 0);
    rectangel(d, {0.1, 0.2, 0}, {0.2, 0.4, 0}, 1);
    // rectangel(d, {0.5, 0.5, 0}, {5, 0.6, 0}, 1);
    for (size_t i = 0; i < d->n; ++i)
    {
        d->rho[i] = d->rho_bar(d->vof[i]);
        d->nu[i] = d->nu_bar(d->vof[i]);
    }

    std::cout << "before run" << std::endl << std::flush;
    solver *s = new solver(d, "out/dmbrk");
    s->run(2);

    delete s;
    delete d;

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

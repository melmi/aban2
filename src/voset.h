/**
 * voset.h
 *
 * by Mohammad Elmi
 * Created on 1 Feb 2014
 */

#ifndef _VOSET_H
#define _VOSET_H

#include "domain.h"

namespace aban2
{

class fsm;
class ireconst;

enum class fullness
{
    full, empty, half
};

class voset
{
    domain *d;
    double *old_ls;
    double h3;
    vector celldims;

    fsm *redistancer;

    // methods to detect interface cells
    void detect_on_interface_cells();
    bool is_on_interface(size_t i, size_t j, size_t k, size_t no);
    fullness get_fullness(size_t no);

    //methods for gauss-sidel iter
    double err();
    double calculate_normals(double *phi);
    void calculate_interface_distances();

public:
    bool *on_interface;
    fullness *fullnesses;
    ireconst *reconsts;

    voset(domain *_d);
    ~voset();

    void make_ls();
};

} // aban2

#endif // _VOSET_H
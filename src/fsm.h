/**
 * fsm.h
 *
 * by Mohammad Elmi
 * Created on 31 Jan 2014
 */

#include "domain.h"

#ifndef _FSM_H
#define _FSM_H

namespace aban2
{

// based on "A FAST SWEEPING METHOD FOR EIKONAL EQUATIONS" by "HONGKAI ZHAO"

class fsm
{
    domain *d;
    const double max_dist_coeff = 3.0;
    double max_dist, h2, inf;

    // methods to detect interface cells
    enum class fullness
    {
        full, empty, half
    };
    bool *on_interface;
    fullness *fullnesses;
    void detect_on_interface_cells();
    bool is_on_interface(size_t i, size_t j, size_t k, size_t no);
    fullness get_fullness(size_t no);

    // methods to redistance
    void init_ls();
    void set_new_phi(size_t i, size_t j, size_t k, size_t no);
    void set_for_bounds(size_t dir, bool asc, size_t &start, int &limit, int &step);
    double calculate_phis(bool iasc, bool jasc, bool kasc);

public:
    fsm(domain *_d);
    ~fsm();

    void redist();
};

} // aban2


#endif // _FSM_H
/**
 * fsm.h
 *
 * by Mohammad Elmi
 * Created on 31 Jan 2014
 */

#include "domain.h"
#include "voset.h"

#ifndef _FSM_H
#define _FSM_H

namespace aban2
{

// based on "A FAST SWEEPING METHOD FOR EIKONAL EQUATIONS" by "HONGKAI ZHAO"

class fsm
{
    domain *d;
    double h2, inf;

    bool *on_interface;
    fullness *fullnesses;

    // methods to redistance
    void init_ls();
    void set_new_phi(size_t i, size_t j, size_t k, size_t no);
    void set_for_bounds(size_t dir, bool asc, size_t &start, int &limit, int &step);
    double calculate_phis(bool iasc, bool jasc, bool kasc);

public:
    fsm(domain *_d, bool *_on_interface, fullness *_fullnesses);
    ~fsm();

    void redist();
};

} // aban2


#endif // _FSM_H
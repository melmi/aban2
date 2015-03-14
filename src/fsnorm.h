/**
 * fsnorm.h
 *
 * by Mohammad Elmi
 * Created on 14 March 2015
 */

#ifndef _FSNORM_H_
#define _FSNORM_H_

#include "common.h"
#include <array>
#include "domain.h"
#include "volreconst.h"

namespace aban2
{

class fsnorm
{
    typedef std::array<std::array<std::array<double, 3>, 3>, 3 > neighbs_t;
    struct molecule_t
    {
        double p, f1, b1, f2, b2;
    };
    static const double epsilon;

    domain *d;
    double vcell; //cell volume

    neighbs_t get_nighb_vals(size_t i, size_t j, size_t k);
    inline void relax_center_val(int i, int j, int k, neighbs_t &n);
    inline void relax_edge_val(int i, int j, int k, neighbs_t &n);
    inline void relax_corner_val(int i, int j, int k, neighbs_t &n);
    inline void relax_neighb_vals(neighbs_t &n);
    inline static double col_sum(double delta, double v1, double v2, double v3);
    inline static double dir_sign(double v1, double v2, double v3);
    inline static double get_column_grad(double p, double f, double b, double delta);
    bool create_column_candidate(vector &v, size_t dir, double sign, molecule_t molecule);
    void create_column_candidates(const neighbs_t &n, vector *candidates, bool *ok);
    bool create_young_candidate(const neighbs_t &n, vector &v);
public:
    fsnorm(domain *_d);
    ~fsnorm();
    vector get_normal(size_t i, size_t j, size_t k, size_t no);
};

const double fsnorm::epsilon = 1e-8;
}

#endif


